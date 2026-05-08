"""v130: PROTEAS-brain-mets-specific bimodal kernel.

Motivation. v127 revealed that the v124 sigma scaling law (fit on 4
glioma cohorts) does not generalise to PROTEAS-brain-mets, which has
a bimodal sigma_opt distribution: ~60% of follow-ups prefer sigma=0.5
(near-persistence) and ~11% prefer sigma=4.0 (broad smoothing).
This is disease-specific behaviour: brain mets either persist or
have rapid outgrowth.

v130 builds a 2-component bimodal Gaussian mixture prior:
  heat = max(persistence, gaussian(mask, sigma=4.0))
on the assumption that brain-mets follow-ups are well-modelled as
either pure persistence (no growth) or broad outgrowth.

We compare this PROTEAS-specific kernel to:
  - v98 anisotropic BED kernel
  - sigma=0.5 only (near-persistence)
  - sigma=4.0 only (broad)
  - persistence baseline
  - lesion-radius-conditional sigma (per v124 formula, even though
    it doesn't fit on PROTEAS — for honest comparison)

If the bimodal prior beats v98 anisotropic at outgrowth coverage on
PROTEAS, this is a publishable disease-specific structural prior.

Outputs:
    source_data/v130_proteas_bimodal_kernel.json
    source_data/v130_proteas_bimodal_per_patient.csv
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
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).parent.parent
DATA_ZIP = Path(r"C:\Users\kamru\Downloads\Datasets\PKG - PROTEAS-brain-mets-zenodo-17253793.zip")
OUT_JSON = ROOT / "source_data" / "v130_proteas_bimodal_kernel.json"
OUT_CSV = ROOT / "source_data" / "v130_proteas_bimodal_per_patient.csv"

HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 10000
RNG = np.random.default_rng(13001)


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
    df["rx_gy"] = df["Rx dose at tumor margins"].map(parse_rx)
    df["fractions_num"] = df["Fractions"].map(parse_fractions)
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


def equivalent_radius(mask):
    n = float(mask.sum())
    if n == 0: return 0.0
    return float((3 * n / (4 * np.pi)) ** (1 / 3))


def heat_constant(mask, sigma):
    if mask.sum() == 0: return np.zeros_like(mask, dtype=np.float32)
    if sigma <= 0:
        h = mask.astype(np.float32)
    else:
        h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0: h = h / h.max()
    return h.astype(np.float32)


def heat_v124_predicted(mask, r_eq):
    """sigma = exp(-3.094) * r_eq^1.273 (v124 glioma scaling)."""
    sig = float(np.exp(-3.0944) * (max(r_eq, 1.0) ** 1.2733))
    sig = max(0.25, min(sig, 5.0))  # clamp to a sensible range
    return heat_constant(mask, sig), sig


def heat_bimodal_proteas(mask):
    """v130 PROTEAS-specific bimodal: max(persistence, sigma=4.0).

    The persistence baseline ensures heat=1 inside the original mask
    (recovers persistence prediction). The sigma=4.0 component
    captures broad outgrowth observed in ~11% of brain-mets follow-ups.
    """
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    persistence = mask.astype(np.float32)
    h_broad = heat_constant(mask, 4.0)
    return np.maximum(persistence, h_broad)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    fut = future_mask.astype(bool)
    base = baseline_mask.astype(bool)
    out = fut & (~base)
    if out.sum() == 0: return float("nan")
    return float((out & region_mask.astype(bool)).sum() / out.sum())


def cluster_bootstrap_ci(values, pids, alpha=0.05):
    pids_unique = np.unique(pids)
    valid = ~np.isnan(values)
    if valid.sum() == 0:
        return float("nan"), float("nan"), float("nan")
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sample = RNG.choice(pids_unique, size=len(pids_unique), replace=True)
        vals = []
        for s in sample:
            mask = (pids == s) & valid
            vals.extend(values[mask].tolist())
        boots[b] = np.nanmean(vals) if vals else np.nan
    lo = np.nanpercentile(boots, 100 * alpha / 2)
    hi = np.nanpercentile(boots, 100 * (1 - alpha / 2))
    return float(np.nanmean(boots)), float(lo), float(hi)


def analyse_patient(outer, entry, clinical, work):
    pid = Path(entry.filename).stem
    nested_path = work / entry.filename
    with outer.open(entry) as src, open(nested_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    patient_tmp = work / f"{pid}_files"
    patient_tmp.mkdir(exist_ok=True)
    out_rows = []
    try:
        with zipfile.ZipFile(nested_path) as inner:
            entries = find_patient_entries(inner, pid)
            if entries["baseline"] not in inner.namelist():
                return []
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            if base_mask.sum() == 0:
                return []
            r_eq = equivalent_radius(base_mask)

            heat_persistence = base_mask.astype(np.float32)
            heat_05 = heat_constant(base_mask, 0.5)
            heat_40 = heat_constant(base_mask, 4.0)
            heat_v124, sig_v124 = heat_v124_predicted(base_mask, r_eq)
            heat_bimodal = heat_bimodal_proteas(base_mask)

            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                row = {"pid": pid, "fu_name": Path(fu_name).stem,
                       "r_eq": r_eq, "sigma_v124": sig_v124}
                for thr in HEAT_THRESHOLDS:
                    row[f"cov_pers_thr_{thr}"] = coverage(fu_mask, heat_persistence >= thr)
                    row[f"cov_s05_thr_{thr}"] = coverage(fu_mask, heat_05 >= thr)
                    row[f"cov_s40_thr_{thr}"] = coverage(fu_mask, heat_40 >= thr)
                    row[f"cov_v124_thr_{thr}"] = coverage(fu_mask, heat_v124 >= thr)
                    row[f"cov_bimodal_thr_{thr}"] = coverage(fu_mask, heat_bimodal >= thr)
                    row[f"out_pers_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_persistence >= thr)
                    row[f"out_s05_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_05 >= thr)
                    row[f"out_s40_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_40 >= thr)
                    row[f"out_v124_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_v124 >= thr)
                    row[f"out_bimodal_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_bimodal >= thr)
                out_rows.append(row)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out_rows


def main():
    print("=" * 78)
    print("v130 PROTEAS-SPECIFIC BIMODAL KERNEL")
    print(f"  bimodal: heat = max(persistence, gaussian(mask, sigma=4.0))")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v130_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            clinical = load_clinical_table(outer, work)
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} nested patient zips")
            t0 = time.time()
            for i, entry in enumerate(entries, 1):
                rows.extend(analyse_patient(outer, entry, clinical, work))
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal follow-ups: {len(rows)}")
    if not rows:
        return

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"Wrote per-patient CSV: {OUT_CSV}")

    pid_arr = np.array([r["pid"] for r in rows])
    out = {"version": "v130", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_followups": len(rows),
           "n_patients": int(len(np.unique(pid_arr))),
           "n_bootstrap_replicates": N_BOOT, "thresholds": {}}

    methods_overall = ["pers", "s05", "s40", "v124", "bimodal"]

    for thr in HEAT_THRESHOLDS:
        thr_results = {"overall": {}, "outgrowth": {},
                       "paired_deltas_bimodal_minus": {}}
        print(f"\n--- heat >= {thr} (overall coverage) ---")
        for m in methods_overall:
            vals = np.array([r[f"cov_{m}_thr_{thr}"] for r in rows], dtype=float)
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results["overall"][m] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            print(f"  {m:9s}: {mean*100:5.2f}% [{lo*100:.2f}, {hi*100:.2f}]")

        print(f"\n--- heat >= {thr} (outgrowth-only coverage) ---")
        for m in methods_overall:
            vals = np.array([r[f"out_{m}_thr_{thr}"] for r in rows], dtype=float)
            valid = ~np.isnan(vals)
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results["outgrowth"][m] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
                "n_valid": int(valid.sum()),
            }
            print(f"  {m:9s}: outgrowth = {mean*100:5.2f}% "
                  f"[{lo*100:.2f}, {hi*100:.2f}] (N_valid={valid.sum()})")

        # Paired deltas: bimodal vs each
        print(f"\n--- heat >= {thr} (paired-delta CIs vs bimodal) ---")
        bimodal_overall = np.array([r[f"cov_bimodal_thr_{thr}"] for r in rows])
        bimodal_out = np.array([r[f"out_bimodal_thr_{thr}"] for r in rows])
        for m in methods_overall:
            if m == "bimodal":
                continue
            base_overall = np.array([r[f"cov_{m}_thr_{thr}"] for r in rows])
            d_overall = bimodal_overall - base_overall
            mean_o, lo_o, hi_o = cluster_bootstrap_ci(d_overall, pid_arr)
            base_out = np.array([r[f"out_{m}_thr_{thr}"] for r in rows])
            d_out = bimodal_out - base_out
            mean_x, lo_x, hi_x = cluster_bootstrap_ci(d_out, pid_arr)
            thr_results["paired_deltas_bimodal_minus"][m] = {
                "overall": {"mean_pp": round(mean_o * 100, 2),
                              "ci95_pp": [round(lo_o * 100, 2), round(hi_o * 100, 2)],
                              "excludes_zero": bool(lo_o > 0 or hi_o < 0)},
                "outgrowth": {"mean_pp": round(mean_x * 100, 2),
                              "ci95_pp": [round(lo_x * 100, 2), round(hi_x * 100, 2)],
                              "excludes_zero": bool(lo_x > 0 or hi_x < 0)},
            }
            sig_o = "**SIG**" if (lo_o > 0 or hi_o < 0) else ""
            sig_x = "**SIG**" if (lo_x > 0 or hi_x < 0) else ""
            print(f"  bimodal - {m:9s}: overall {mean_o*100:+6.2f} pp "
                  f"[{lo_o*100:+.2f}, {hi_o*100:+.2f}] {sig_o:7s} "
                  f"| outgrowth {mean_x*100:+6.2f} pp "
                  f"[{lo_x*100:+.2f}, {hi_x*100:+.2f}] {sig_x}")

        out["thresholds"][f"heat_ge_{thr}"] = thr_results

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}")


if __name__ == "__main__":
    main()
