"""v140: bimodal + U-Net ensemble on PROTEAS.

The v139 LOPO experiment showed:
  - Hand-crafted bimodal kernel (sigma=7): 60.07% overall, 22.57% outgrowth
  - Learned U-Net: 10.95% overall (-49 pp), 38.79% outgrowth (+16.22 pp)
  - The two are COMPLEMENTARY.

v140 builds the natural ensemble:
  heat = max(bimodal_at_sigma_7, U_Net_logits_sigmoid)
and tests whether the ensemble beats both individually on
both overall AND outgrowth.

Reuses the v139 LOPO setup (stride=4, 11 test patients, 36 follow-ups).

Reads the v139 per-patient CSV directly to extract per-fu coverage
values; this avoids re-training and re-saving raw model outputs.

NOTE: The v139 CSV has the per-patient bimodal and learned-U-Net
coverages already. To compute the ensemble we need the ACTUAL voxel
predictions, not just the coverage values. So this script re-runs
the v139 pipeline, but stores the ensemble predictions inline.

For brevity, instead of re-training, we approximate the ensemble
result analytically:
  ensemble = max(bimodal, learned)
For overall coverage:
  ensemble_overall_coverage = coverage(future, max(bim_heat, lrn_heat) >= 0.5)
                             >= max(bimodal_overall, learned_overall)
For outgrowth coverage similarly.

To get a tight estimate we re-run U-Net on the v139 LOPO test patients,
but with checkpoint reuse from v139 IF available. Since we don't have
saved checkpoints, we re-train per-fold (5 min/fold, 11 folds = ~14 min).

Outputs:
    source_data/v140_bimodal_unet_ensemble.json
    source_data/v140_bimodal_unet_ensemble_per_patient.csv
"""
from __future__ import annotations

import csv
import gc
import json
import re
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.ndimage import gaussian_filter, zoom

ROOT = Path(__file__).parent.parent
DATA_ZIP = Path(r"C:\Users\kamru\Downloads\Datasets\PKG - PROTEAS-brain-mets-zenodo-17253793.zip")
OUT_JSON = ROOT / "source_data" / "v140_bimodal_unet_ensemble.json"
OUT_CSV = ROOT / "source_data" / "v140_bimodal_unet_ensemble_per_patient.csv"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

SIGMA_BROAD = 7.0
EPOCHS = 30
LR = 1e-3
LOPO_STRIDE = 4
TARGET_SHAPE = (32, 64, 64)
N_BOOT = 5000
RNG = np.random.default_rng(14001)


def parse_rx(value):
    if pd.isna(value): return None
    m = re.search(r"(\d+(?:\.\d+)?)", str(value))
    return float(m.group(1)) if m else None


def parse_fractions(value):
    if pd.isna(value): return None
    try: return int(float(value))
    except (ValueError, TypeError): return None


def load_clinical_table(outer, tmpdir):
    name = "PROTEAS-Clinical_and_demographic_data.xlsx"
    out = tmpdir / name
    out.write_bytes(outer.read(name))
    df = pd.read_excel(out, sheet_name="PROTEAS")
    df["pid"] = df["Patient ID (Zenodo)"].astype(str).str.strip()
    return df


def load_nii_from_inner(inner, name, tmpdir):
    out = tmpdir / Path(name).name
    out.write_bytes(inner.read(name))
    img = nib.load(str(out))
    return np.asanyarray(img.dataobj).astype(np.float32)


def find_patient_entries(inner, pid):
    names = inner.namelist()
    prefix = f"{pid}/"
    seg_dirs = [f"{prefix}tumor segmentation/", f"{prefix}tumor_segmentation/"]
    baseline = next(
        (f"{seg_dir}{pid}_tumor_mask_baseline.nii.gz" for seg_dir in seg_dirs
         if f"{seg_dir}{pid}_tumor_mask_baseline.nii.gz" in names),
        f"{prefix}tumor segmentation/{pid}_tumor_mask_baseline.nii.gz",
    )
    followups = sorted([
        n for n in names
        if any(n.startswith(f"{seg_dir}{pid}_tumor_mask_fu") for seg_dir in seg_dirs)
        and n.endswith(".nii.gz")
    ])
    return {"baseline": baseline, "followups": followups}


def heat_constant(mask, sigma):
    if mask.sum() == 0: return np.zeros_like(mask, dtype=np.float32)
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0: h = h / h.max()
    return h.astype(np.float32)


def heat_bimodal(mask, sigma_broad):
    persistence = mask.astype(np.float32)
    h_broad = heat_constant(mask, sigma_broad)
    return np.maximum(persistence, h_broad)


def resize_to_target(arr, target_shape):
    factors = [t/s for t, s in zip(target_shape, arr.shape)]
    if arr.dtype == bool or np.array_equal(arr, arr.astype(bool).astype(arr.dtype)):
        return zoom(arr.astype(np.float32), factors, order=0).astype(np.float32)
    return zoom(arr, factors, order=1).astype(np.float32)


def overall_coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    fut = future_mask.astype(bool); base = baseline_mask.astype(bool)
    out = fut & (~base)
    if out.sum() == 0: return float("nan")
    return float((out & region_mask.astype(bool)).sum() / out.sum())


def vectorised_cluster_bootstrap(values, alpha=0.05):
    values = np.asarray(values, dtype=float)
    valid = ~np.isnan(values)
    if valid.sum() == 0:
        return float("nan"), float("nan"), float("nan")
    v = values[valid]
    n = len(v)
    sample_idx = RNG.integers(0, n, size=(N_BOOT, n))
    boot_means = v[sample_idx].mean(axis=1)
    lo = float(np.percentile(boot_means, 100 * alpha / 2))
    hi = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return float(boot_means.mean()), lo, hi


class UNet3D(nn.Module):
    def __init__(self, in_ch=2, base=24):
        super().__init__()
        self.enc1 = self._block(in_ch, base)
        self.enc2 = self._block(base, base * 2)
        self.enc3 = self._block(base * 2, base * 4)
        self.dec2 = self._block(base * 4 + base * 2, base * 2)
        self.dec1 = self._block(base * 2 + base, base)
        self.out = nn.Conv3d(base, 1, 1)
        self.pool = nn.MaxPool3d(2)
        self.up = nn.Upsample(scale_factor=2, mode="trilinear", align_corners=False)

    def _block(self, in_ch, out_ch):
        return nn.Sequential(
            nn.Conv3d(in_ch, out_ch, 3, padding=1),
            nn.GroupNorm(8, out_ch), nn.GELU(),
            nn.Conv3d(out_ch, out_ch, 3, padding=1),
            nn.GroupNorm(8, out_ch), nn.GELU(),
        )

    def forward(self, x):
        e1 = self.enc1(x)
        e2 = self.enc2(self.pool(e1))
        e3 = self.enc3(self.pool(e2))
        d2 = self.dec2(torch.cat([self.up(e3), e2], dim=1))
        d1 = self.dec1(torch.cat([self.up(d2), e1], dim=1))
        return self.out(d1)


def focal_dice_loss(logits, target, alpha=0.95, gamma=2.0, smooth=1e-5):
    p = torch.sigmoid(logits)
    p_t = p * target + (1 - p) * (1 - target)
    alpha_t = alpha * target + (1 - alpha) * (1 - target)
    focal = -alpha_t * (1 - p_t) ** gamma * torch.log(p_t.clamp(1e-7, 1 - 1e-7))
    return focal.mean() + (1 - (2 * (p * target).sum() + smooth) /
                              (p.sum() + target.sum() + smooth))


def collect_patient_data(outer, entry, work):
    pid = Path(entry.filename).stem
    nested_path = work / entry.filename
    with outer.open(entry) as src, open(nested_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    patient_tmp = work / f"{pid}_files"
    patient_tmp.mkdir(exist_ok=True)
    out = []
    try:
        with zipfile.ZipFile(nested_path) as inner:
            entries = find_patient_entries(inner, pid)
            if entries["baseline"] not in inner.namelist(): return []
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            if base_mask.sum() == 0: return []
            heat = heat_bimodal(base_mask, SIGMA_BROAD)
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any(): continue
                outgrowth = fu_mask & (~base_mask)
                m_r = resize_to_target(base_mask.astype(np.float32), TARGET_SHAPE) > 0.5
                fu_r = resize_to_target(fu_mask.astype(np.float32), TARGET_SHAPE) > 0.5
                out_r = resize_to_target(outgrowth.astype(np.float32), TARGET_SHAPE) > 0.5
                heat_r = resize_to_target(heat, TARGET_SHAPE)
                out.append({"pid": pid, "fu_name": Path(fu_name).stem,
                            "mask_r": m_r, "fu_r": fu_r,
                            "outgrowth_r": out_r, "heat_bimodal_r": heat_r})
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out


def main():
    print("=" * 78, flush=True)
    print("v140 BIMODAL + U-NET ENSEMBLE ON PROTEAS", flush=True)
    print(f"  device: {DEVICE}; epochs={EPOCHS}; LOPO stride={LOPO_STRIDE}", flush=True)
    print(f"  ensemble: heat = max(bimodal_sigma_7, U_Net_sigmoid)", flush=True)
    print("=" * 78, flush=True)

    print("\nLoading PROTEAS data...", flush=True)
    all_data = []
    with tempfile.TemporaryDirectory(prefix="proteas_v140_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} patient zips", flush=True)
            t0 = time.time()
            for i, entry in enumerate(entries, 1):
                all_data.extend(collect_patient_data(outer, entry, work))
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} ({time.time()-t0:.0f}s)", flush=True)
        print(f"\nLoaded {len(all_data)} follow-ups across {len(set(d['pid'] for d in all_data))} patients", flush=True)

        unique_pids = sorted(set(d["pid"] for d in all_data))
        test_pids = unique_pids[::LOPO_STRIDE]
        print(f"\nLOPO test patients (stride={LOPO_STRIDE}): {len(test_pids)} of {len(unique_pids)}", flush=True)

        results = []
        for fold_idx, test_pid in enumerate(test_pids):
            t0 = time.time()
            train_data = [d for d in all_data if d["pid"] != test_pid]
            test_data = [d for d in all_data if d["pid"] == test_pid]
            if not test_data: continue

            X_train = np.stack([np.stack([d["mask_r"], d["heat_bimodal_r"]], axis=0)
                                  for d in train_data]).astype(np.float32)
            Y_train = np.stack([d["outgrowth_r"] for d in train_data]).astype(np.float32)[:, None]
            Xt = torch.from_numpy(X_train).to(DEVICE)
            Yt = torch.from_numpy(Y_train).to(DEVICE)

            model = UNet3D(in_ch=2, base=24).to(DEVICE)
            opt = torch.optim.AdamW(model.parameters(), lr=LR)
            n_train = len(train_data)
            for ep in range(EPOCHS):
                model.train()
                perm = np.random.permutation(n_train)
                for idx in perm:
                    xb = Xt[idx:idx+1]; yb = Yt[idx:idx+1]
                    logits = model(xb)
                    loss = focal_dice_loss(logits, yb)
                    opt.zero_grad(); loss.backward(); opt.step()

            model.eval()
            with torch.no_grad():
                for d in test_data:
                    x = np.stack([d["mask_r"], d["heat_bimodal_r"]], axis=0).astype(np.float32)
                    xt = torch.from_numpy(x[None]).to(DEVICE)
                    logits = model(xt)
                    pred = torch.sigmoid(logits).cpu().numpy()[0, 0]
                    fu = d["fu_r"]; base = d["mask_r"]; bim = d["heat_bimodal_r"]
                    learned_region = pred >= 0.5
                    bimodal_region = bim >= 0.5
                    ensemble_heat = np.maximum(bim, pred)
                    ensemble_region = ensemble_heat >= 0.5
                    results.append({
                        "test_pid": test_pid, "fu_name": d["fu_name"],
                        "learned_overall": overall_coverage(fu, learned_region),
                        "learned_outgrowth": outgrowth_coverage(fu, base, learned_region),
                        "bimodal_overall": overall_coverage(fu, bimodal_region),
                        "bimodal_outgrowth": outgrowth_coverage(fu, base, bimodal_region),
                        "ensemble_overall": overall_coverage(fu, ensemble_region),
                        "ensemble_outgrowth": outgrowth_coverage(fu, base, ensemble_region),
                    })
            del model, Xt, Yt; torch.cuda.empty_cache(); gc.collect()
            print(f"  fold {fold_idx+1}/{len(test_pids)} pid={test_pid}: "
                  f"{time.time()-t0:.0f}s ({len(test_data)} fu)", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if results:
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            w.writeheader(); w.writerows(results)
        print(f"\nWrote per-patient CSV: {OUT_CSV}", flush=True)

    print(f"\n=== AGGREGATE on {len(results)} test follow-ups ===", flush=True)
    out = {"version": "v140",
           "experiment": "Bimodal + U-Net ensemble on PROTEAS",
           "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_test_followups": len(results), "ci_method": "vectorised cluster-bootstrap (per-fu)",
           "n_bootstrap": N_BOOT, "methods": {}}
    for method in ["bimodal", "learned", "ensemble"]:
        ovr = np.array([r[f"{method}_overall"] for r in results], dtype=float)
        out_v = np.array([r[f"{method}_outgrowth"] for r in results], dtype=float)
        m_o, lo_o, hi_o = vectorised_cluster_bootstrap(ovr)
        m_x, lo_x, hi_x = vectorised_cluster_bootstrap(out_v)
        out["methods"][method] = {
            "overall_mean_pct": round(m_o * 100, 2),
            "overall_ci95_pct": [round(lo_o * 100, 2), round(hi_o * 100, 2)],
            "outgrowth_mean_pct": round(m_x * 100, 2),
            "outgrowth_ci95_pct": [round(lo_x * 100, 2), round(hi_x * 100, 2)],
        }
        print(f"  {method:10s}: overall {m_o*100:5.2f}% [{lo_o*100:.2f}, {hi_o*100:.2f}]  "
              f"outgrowth {m_x*100:5.2f}% [{lo_x*100:.2f}, {hi_x*100:.2f}]", flush=True)

    print(f"\n=== Paired delta CIs ===", flush=True)
    out["paired_deltas"] = {}
    bim_o = np.array([r["bimodal_overall"] for r in results], dtype=float)
    bim_x = np.array([r["bimodal_outgrowth"] for r in results], dtype=float)
    lrn_o = np.array([r["learned_overall"] for r in results], dtype=float)
    lrn_x = np.array([r["learned_outgrowth"] for r in results], dtype=float)
    ens_o = np.array([r["ensemble_overall"] for r in results], dtype=float)
    ens_x = np.array([r["ensemble_outgrowth"] for r in results], dtype=float)

    for label, d_o, d_x in [
        ("ensemble - bimodal", ens_o - bim_o, ens_x - bim_x),
        ("ensemble - learned", ens_o - lrn_o, ens_x - lrn_x),
    ]:
        m_o, lo_o, hi_o = vectorised_cluster_bootstrap(d_o)
        m_x, lo_x, hi_x = vectorised_cluster_bootstrap(d_x)
        sig_o = "**SIG**" if (lo_o > 0 or hi_o < 0) else ""
        sig_x = "**SIG**" if (lo_x > 0 or hi_x < 0) else ""
        out["paired_deltas"][label] = {
            "overall_pp": [round(m_o * 100, 2), round(lo_o * 100, 2), round(hi_o * 100, 2),
                           bool(lo_o > 0 or hi_o < 0)],
            "outgrowth_pp": [round(m_x * 100, 2), round(lo_x * 100, 2), round(hi_x * 100, 2),
                              bool(lo_x > 0 or hi_x < 0)],
        }
        print(f"  {label:25s}: overall {m_o*100:+6.2f} pp [{lo_o*100:+.2f}, {hi_o*100:+.2f}] "
              f"{sig_o:7s} outgrowth {m_x*100:+6.2f} pp [{lo_x*100:+.2f}, {hi_x*100:+.2f}] {sig_x}",
              flush=True)

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
