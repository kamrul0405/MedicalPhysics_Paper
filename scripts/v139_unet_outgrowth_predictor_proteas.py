"""v139: GPU U-Net trained on PROTEAS LOPO for outgrowth prediction.

Tests whether a learned 3D U-Net beats the v133 hand-crafted bimodal
kernel (max(persistence, sigma=7)) at predicting outgrowth voxels
on PROTEAS-brain-mets.

If learned matches or beats the bimodal kernel, that validates the
hand-crafted inductive bias and adds a deep-learning extension. If
learned fails (likely on small cohort), that's evidence FOR
hand-crafted physics-grounded priors over deep learning on small
cohorts.

Architecture:
  - 3D U-Net with downsample 32 -> 64 -> 128 channels
  - Input: 2 channels (baseline mask + bimodal heat map at sigma=7)
  - Output: per-voxel probability of being a future-lesion outgrowth voxel
  - Loss: focal BCE (alpha=0.95, gamma=2) + Dice loss
  - LOPO: leave-one-patient-out, but with limited folds (every 5th patient)
    to keep total training time tractable

This compares LEARNED outgrowth prediction to HAND-CRAFTED bimodal
prediction on the SAME held-out patients.

Outputs:
    source_data/v139_unet_outgrowth_predictor.json
    source_data/v139_unet_outgrowth_per_patient.csv
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
OUT_JSON = ROOT / "source_data" / "v139_unet_outgrowth_predictor.json"
OUT_CSV = ROOT / "source_data" / "v139_unet_outgrowth_per_patient.csv"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

SIGMA_BROAD = 7.0
EPOCHS = 30
BATCH_SIZE = 1
LR = 1e-3
HEAT_THRESHOLDS = [0.5]
LOPO_STRIDE = 4   # train on n-1 patients, test on every 4th patient
TARGET_SHAPE = (32, 64, 64)


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
    """Trilinear resize for floats; nearest for binary."""
    factors = [t/s for t, s in zip(target_shape, arr.shape)]
    if arr.dtype == bool or np.array_equal(arr, arr.astype(bool).astype(arr.dtype)):
        return zoom(arr.astype(np.float32), factors, order=0).astype(np.float32)
    return zoom(arr, factors, order=1).astype(np.float32)


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    fut = future_mask.astype(bool); base = baseline_mask.astype(bool)
    out = fut & (~base)
    if out.sum() == 0: return float("nan")
    return float((out & region_mask.astype(bool)).sum() / out.sum())


def overall_coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


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
    """Focal BCE + Dice. Target is 0/1 outgrowth mask."""
    p = torch.sigmoid(logits)
    p_t = p * target + (1 - p) * (1 - target)
    alpha_t = alpha * target + (1 - alpha) * (1 - target)
    focal = -alpha_t * (1 - p_t) ** gamma * torch.log(p_t.clamp(1e-7, 1 - 1e-7))
    focal_loss = focal.mean()
    intersection = (p * target).sum()
    dice_loss = 1 - (2 * intersection + smooth) / (p.sum() + target.sum() + smooth)
    return focal_loss + dice_loss


def collect_patient_data(outer, entry, work):
    """Returns list of (mask_resized, fu_resized, outgrowth_resized,
    bimodal_resized, baseline_orig_shape) per follow-up."""
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
                # Resize to fixed shape
                m_r = resize_to_target(base_mask.astype(np.float32), TARGET_SHAPE) > 0.5
                fu_r = resize_to_target(fu_mask.astype(np.float32), TARGET_SHAPE) > 0.5
                out_r = resize_to_target(outgrowth.astype(np.float32), TARGET_SHAPE) > 0.5
                heat_r = resize_to_target(heat, TARGET_SHAPE)
                out.append({
                    "pid": pid, "fu_name": Path(fu_name).stem,
                    "mask_r": m_r, "fu_r": fu_r,
                    "outgrowth_r": out_r, "heat_bimodal_r": heat_r,
                })
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out


def main():
    print("=" * 78, flush=True)
    print("v139 GPU U-NET OUTGROWTH PREDICTOR ON PROTEAS", flush=True)
    print(f"  device: {DEVICE}; epochs={EPOCHS}; LOPO stride={LOPO_STRIDE}", flush=True)
    print(f"  target shape: {TARGET_SHAPE}; sigma_broad bimodal feat: {SIGMA_BROAD}", flush=True)
    print("=" * 78, flush=True)

    print("\nLoading PROTEAS data + computing bimodal heat features...", flush=True)
    all_data = []
    with tempfile.TemporaryDirectory(prefix="proteas_v139_") as td:
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

        # LOPO: leave-one-patient-out, every LOPO_STRIDE-th patient
        unique_pids = sorted(set(d["pid"] for d in all_data))
        test_pids = unique_pids[::LOPO_STRIDE]
        print(f"\nLOPO test patients (stride={LOPO_STRIDE}): {len(test_pids)} of {len(unique_pids)}", flush=True)

        results = []
        for fold_idx, test_pid in enumerate(test_pids):
            t0 = time.time()
            train_data = [d for d in all_data if d["pid"] != test_pid]
            test_data = [d for d in all_data if d["pid"] == test_pid]
            if not test_data:
                continue

            # Build tensors
            X_train = np.stack([
                np.stack([d["mask_r"], d["heat_bimodal_r"]], axis=0)
                for d in train_data
            ]).astype(np.float32)
            Y_train = np.stack([d["outgrowth_r"] for d in train_data]).astype(np.float32)[:, None]

            Xt = torch.from_numpy(X_train).to(DEVICE)
            Yt = torch.from_numpy(Y_train).to(DEVICE)

            # Train
            model = UNet3D(in_ch=2, base=24).to(DEVICE)
            opt = torch.optim.AdamW(model.parameters(), lr=LR)

            n_train = len(train_data)
            for ep in range(EPOCHS):
                model.train()
                perm = np.random.permutation(n_train)
                ep_loss = 0.0
                for idx in perm:
                    xb = Xt[idx:idx+1]
                    yb = Yt[idx:idx+1]
                    logits = model(xb)
                    loss = focal_dice_loss(logits, yb)
                    opt.zero_grad(); loss.backward(); opt.step()
                    ep_loss += float(loss.item())

            # Evaluate on test patient(s)
            model.eval()
            with torch.no_grad():
                for d in test_data:
                    x = np.stack([d["mask_r"], d["heat_bimodal_r"]], axis=0).astype(np.float32)
                    xt = torch.from_numpy(x[None]).to(DEVICE)
                    logits = model(xt)
                    pred = torch.sigmoid(logits).cpu().numpy()[0, 0]
                    fu_mask = d["fu_r"]
                    base_mask = d["mask_r"]
                    # Threshold pred at 0.5
                    pred_region = pred >= 0.5
                    bimodal_region = d["heat_bimodal_r"] >= 0.5
                    learned_overall = overall_coverage(fu_mask, pred_region)
                    learned_outgrowth = outgrowth_coverage(fu_mask, base_mask, pred_region)
                    bimodal_overall = overall_coverage(fu_mask, bimodal_region)
                    bimodal_outgrowth = outgrowth_coverage(fu_mask, base_mask, bimodal_region)
                    results.append({
                        "fold_idx": fold_idx,
                        "test_pid": test_pid,
                        "fu_name": d["fu_name"],
                        "learned_overall": learned_overall,
                        "learned_outgrowth": learned_outgrowth,
                        "bimodal_overall": bimodal_overall,
                        "bimodal_outgrowth": bimodal_outgrowth,
                        "delta_overall_learned_minus_bimodal": (learned_overall - bimodal_overall) if not np.isnan(learned_overall) and not np.isnan(bimodal_overall) else float("nan"),
                        "delta_outgrowth_learned_minus_bimodal": (learned_outgrowth - bimodal_outgrowth) if not np.isnan(learned_outgrowth) and not np.isnan(bimodal_outgrowth) else float("nan"),
                    })
            del model, Xt, Yt
            torch.cuda.empty_cache()
            gc.collect()
            print(f"  fold {fold_idx+1}/{len(test_pids)} pid={test_pid}: "
                  f"{time.time()-t0:.0f}s ({len(test_data)} fu)", flush=True)

    # Save per-patient CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if results:
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            w.writeheader(); w.writerows(results)
        print(f"\nWrote per-patient CSV: {OUT_CSV}", flush=True)

    # Aggregate
    n = len(results)
    if n == 0:
        print("No results", flush=True)
        return
    learned_overall = np.array([r["learned_overall"] for r in results])
    learned_outgrowth = np.array([r["learned_outgrowth"] for r in results])
    bimodal_overall = np.array([r["bimodal_overall"] for r in results])
    bimodal_outgrowth = np.array([r["bimodal_outgrowth"] for r in results])
    d_overall = np.array([r["delta_overall_learned_minus_bimodal"] for r in results])
    d_outgrowth = np.array([r["delta_outgrowth_learned_minus_bimodal"] for r in results])

    print(f"\n=== AGGREGATE on {n} test follow-ups ===", flush=True)
    print(f"  Bimodal kernel (sigma=7) at heat>=0.5:", flush=True)
    print(f"    overall:   {np.nanmean(bimodal_overall)*100:.2f}%", flush=True)
    print(f"    outgrowth: {np.nanmean(bimodal_outgrowth)*100:.2f}%", flush=True)
    print(f"  Learned U-Net at thresh>=0.5:", flush=True)
    print(f"    overall:   {np.nanmean(learned_overall)*100:.2f}%", flush=True)
    print(f"    outgrowth: {np.nanmean(learned_outgrowth)*100:.2f}%", flush=True)
    print(f"  Paired delta (learned - bimodal):", flush=True)
    print(f"    overall:   {np.nanmean(d_overall)*100:+.2f} pp", flush=True)
    print(f"    outgrowth: {np.nanmean(d_outgrowth)*100:+.2f} pp", flush=True)

    out = {
        "version": "v139",
        "experiment": "GPU U-Net learned outgrowth predictor vs hand-crafted bimodal kernel on PROTEAS",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "device": str(DEVICE),
        "epochs": EPOCHS,
        "lopo_stride": LOPO_STRIDE,
        "target_shape": list(TARGET_SHAPE),
        "n_test_followups": n,
        "n_test_patients": int(len(set(r["test_pid"] for r in results))),
        "aggregate": {
            "learned_overall_mean_pct": float(np.nanmean(learned_overall) * 100),
            "learned_outgrowth_mean_pct": float(np.nanmean(learned_outgrowth) * 100),
            "bimodal_overall_mean_pct": float(np.nanmean(bimodal_overall) * 100),
            "bimodal_outgrowth_mean_pct": float(np.nanmean(bimodal_outgrowth) * 100),
            "delta_overall_mean_pp": float(np.nanmean(d_overall) * 100),
            "delta_outgrowth_mean_pp": float(np.nanmean(d_outgrowth) * 100),
            "delta_overall_sd_pp": float(np.nanstd(d_overall, ddof=1) * 100),
            "delta_outgrowth_sd_pp": float(np.nanstd(d_outgrowth, ddof=1) * 100),
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
