"""v95: α/β sensitivity sweep for the BED-aware spatially-varying heat-kernel.

Tests whether the BED-aware coverage advantage (§3.9) is robust to the choice
of α/β ratio. Sweeps tumour α/β ∈ {8, 10, 12} Gy and re-evaluates per-voxel
BED-aware kernel coverage on PROTEAS-brain-mets.

Outputs:
  source_data/v95_alpha_beta_sensitivity.json
"""
from __future__ import annotations

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
OUT_JSON = ROOT / "source_data" / "v95_alpha_beta_sensitivity.json"

SIGMA_LO = 1.5
SIGMA_HI = 4.0
SIGMA_BASE = 2.5
ALPHA_BETA_GRID = [8.0, 10.0, 12.0]
HEAT_THRESHOLDS = [0.5, 0.8]


def parse_rx(value):
    if pd.isna(value):
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", str(value))
    return float(match.group(1)) if match else None


def load_clinical_table(outer, tmpdir):
    name = "PROTEAS-Clinical_and_demographic_data.xlsx"
    out = tmpdir / name
    out.write_bytes(outer.read(name))
    df = pd.read_excel(out, sheet_name="PROTEAS")
    df["pid"] = df["Patient ID (Zenodo)"].astype(str).str.strip()
    df["rx_gy"] = df["Rx dose at tumor margins"].map(parse_rx)
    df["fractions_num"] = pd.to_numeric(df["Fractions"], errors="coerce")
    return df


def load_nii_from_inner(inner, name, tmpdir):
    out = tmpdir / Path(name).name
    out.write_bytes(inner.read(name))
    img = nib.load(str(out))
    return np.asanyarray(img.dataobj).astype(np.float32)


def find_patient_entries(inner, pid):
    names = inner.namelist()
    prefix = f"{pid}/"
    dose = f"{prefix}{pid}_RTP.nii.gz"
    seg_dirs = [f"{prefix}tumor segmentation/", f"{prefix}tumor_segmentation/"]
    baseline = next(
        (f"{seg_dir}{pid}_tumor_mask_baseline.nii.gz" for seg_dir in seg_dirs
         if f"{seg_dir}{pid}_tumor_mask_baseline.nii.gz" in names),
        f"{prefix}tumor segmentation/{pid}_tumor_mask_baseline.nii.gz",
    )
    followups = sorted([
        n for n in names
        if any(n.startswith(f"{seg_dir}{pid}_tumor_mask_fu") for seg_dir in seg_dirs) and n.endswith(".nii.gz")
    ])
    return {"dose": dose, "baseline": baseline, "followups": followups}


def constant_heat(mask, sigma=SIGMA_BASE):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def bed_aware_heat(mask, dose_total, n_fx, ab):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    h_lo = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_LO)
    h_hi = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_HI)
    if h_lo.max() > 0:
        h_lo = h_lo / h_lo.max()
    if h_hi.max() > 0:
        h_hi = h_hi / h_hi.max()
    dose_per_fx = dose_total / max(n_fx, 1)
    bed = n_fx * dose_per_fx * (1 + dose_per_fx / ab)
    bed_max = bed.max() if bed.max() > 0 else 1.0
    bed_norm = np.clip(bed / bed_max, 0.0, 1.0)
    w = bed_norm
    heat = w * h_lo + (1 - w) * h_hi
    if heat.max() > 0:
        heat = heat / heat.max()
    return heat.astype(np.float32)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0:
        return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def analyse_patient(outer, entry, clinical, work, ab):
    pid = Path(entry.filename).stem
    row = clinical.loc[clinical["pid"] == pid]
    if row.empty:
        return []
    rx = parse_rx(row.iloc[0]["Rx dose at tumor margins"])
    n_fx = float(row.iloc[0]["fractions_num"]) if pd.notna(row.iloc[0]["fractions_num"]) else float("nan")
    if rx is None or not (np.isfinite(rx) and np.isfinite(n_fx) and n_fx > 0):
        return []
    n_fx = int(n_fx)
    nested_path = work / entry.filename
    with outer.open(entry) as src, open(nested_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)
    patient_tmp = work / f"{pid}_files"
    patient_tmp.mkdir(exist_ok=True)
    out_rows = []
    try:
        with zipfile.ZipFile(nested_path) as inner:
            entries = find_patient_entries(inner, pid)
            if entries["dose"] not in inner.namelist() or entries["baseline"] not in inner.namelist():
                return []
            dose = load_nii_from_inner(inner, entries["dose"], patient_tmp)
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            if dose.shape != base_mask.shape:
                return []
            heat_const = constant_heat(base_mask)
            heat_bed = bed_aware_heat(base_mask, dose, n_fx, ab)
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                for thr in HEAT_THRESHOLDS:
                    cov_const = coverage(fu_mask, heat_const >= thr)
                    cov_bed = coverage(fu_mask, heat_bed >= thr)
                    cov_dose95 = coverage(fu_mask, dose >= 0.95 * rx)
                    out_rows.append({
                        "pid": pid,
                        "alpha_beta": ab,
                        "threshold": thr,
                        "coverage_const_sigma": cov_const,
                        "coverage_bed_aware": cov_bed,
                        "coverage_dose_95pct_rx": cov_dose95,
                    })
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try:
            nested_path.unlink()
        except OSError:
            pass
    return out_rows


def main():
    import gc
    print("=" * 78)
    print("v95 ALPHA/BETA SENSITIVITY SWEEP (BED-aware kernel)")
    print(f"  Sweeping tumour alpha/beta in {ALPHA_BETA_GRID} Gy")
    print("=" * 78)

    out_all = {"version": "v95", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "alpha_beta_grid_gy": ALPHA_BETA_GRID, "results": {}}

    with tempfile.TemporaryDirectory(prefix="proteas_v95ab_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            clinical = load_clinical_table(outer, work)
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} nested patient zips")
            for ab in ALPHA_BETA_GRID:
                print(f"\n  --- alpha/beta = {ab} Gy ---")
                t0 = time.time()
                rows = []
                for i, entry in enumerate(entries, 1):
                    pr = analyse_patient(outer, entry, clinical, work, ab)
                    rows.extend(pr)
                    del pr
                    gc.collect()
                    if i % 10 == 0 or i == len(entries):
                        print(f"    {i}/{len(entries)} processed  ({time.time()-t0:.0f}s)", flush=True)
                # Aggregate per threshold
                per_thr = {}
                for thr in HEAT_THRESHOLDS:
                    sub = [r for r in rows if r["threshold"] == thr]
                    cov_const = np.array([r["coverage_const_sigma"] for r in sub if np.isfinite(r["coverage_const_sigma"])])
                    cov_bed = np.array([r["coverage_bed_aware"] for r in sub if np.isfinite(r["coverage_bed_aware"])])
                    per_thr[f"heat_ge_{thr}"] = {
                        "n_rows": len(sub),
                        "constant_sigma_mean_pct": round(float(cov_const.mean() * 100), 2),
                        "bed_aware_mean_pct": round(float(cov_bed.mean() * 100), 2),
                        "delta_bed_minus_const_pp": round(float((cov_bed.mean() - cov_const.mean()) * 100), 2),
                    }
                out_all["results"][f"alpha_beta_{ab}"] = per_thr
                print(f"    -> heat>=0.50: bed-aware {per_thr['heat_ge_0.5']['bed_aware_mean_pct']:.2f}% vs const {per_thr['heat_ge_0.5']['constant_sigma_mean_pct']:.2f}% (Delta={per_thr['heat_ge_0.5']['delta_bed_minus_const_pp']:+.2f} pp)")
                print(f"    -> heat>=0.80: bed-aware {per_thr['heat_ge_0.8']['bed_aware_mean_pct']:.2f}% vs const {per_thr['heat_ge_0.8']['constant_sigma_mean_pct']:.2f}% (Delta={per_thr['heat_ge_0.8']['delta_bed_minus_const_pp']:+.2f} pp)")
                # Save intermediate state
                OUT_JSON.write_text(json.dumps(out_all, indent=2))
                del rows
                gc.collect()

    print(f"\nFinal saved {OUT_JSON}")


if __name__ == "__main__":
    main()
