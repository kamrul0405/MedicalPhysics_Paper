"""v109: Heat-equation evolution-time σ sweep on PROTEAS-brain-mets.

Sweeps the heat-kernel scale parameter σ ∈ {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}
voxels — corresponding to heat-equation evolution times t = σ²/2 ∈
{0.5, 1.13, 2.0, 3.13, 4.5, 6.13, 8.0} voxel-time units — and re-evaluates
future-lesion coverage at heat ≥ 0.50 and heat ≥ 0.80 thresholds.

Empirically maps the optimal σ across the PROTEAS cohort and validates the
scale-space framework (Lindeberg 1994; Witkin 1983) for choosing σ.

Outputs:
    source_data/v109_heat_equation_sigma_sweep.json
"""
from __future__ import annotations

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
OUT_JSON = ROOT / "source_data" / "v109_heat_equation_sigma_sweep.json"

SIGMA_GRID = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]  # voxels (= 1 mm at 1 mm isotropic)
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


def heat_kernel(mask, sigma):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0:
        return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def analyse_patient(outer, entry, clinical, work):
    pid = Path(entry.filename).stem
    row = clinical.loc[clinical["pid"] == pid]
    if row.empty:
        return []
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
            heat_maps = {sig: heat_kernel(base_mask, sig) for sig in SIGMA_GRID}
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                for sig in SIGMA_GRID:
                    for thr in HEAT_THRESHOLDS:
                        out_rows.append({
                            "pid": pid,
                            "sigma_voxels": sig,
                            "evolution_time_t": sig ** 2 / 2,
                            "threshold": thr,
                            "coverage": coverage(fu_mask, heat_maps[sig] >= thr),
                        })
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try:
            nested_path.unlink()
        except OSError:
            pass
    return out_rows


def main():
    print("=" * 78)
    print("v109 HEAT-EQUATION EVOLUTION-TIME SIGMA SWEEP ON PROTEAS")
    print(f"  sigma grid: {SIGMA_GRID} voxels")
    print(f"  evolution-time t = sigma^2 / 2: {[round(s**2/2, 2) for s in SIGMA_GRID]} voxel-time")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v109_") as td:
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
                pr = analyse_patient(outer, entry, clinical, work)
                rows.extend(pr)
                del pr
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} processed  ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal rows: {len(rows)}")

    out = {"version": "v109", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "sigma_grid": SIGMA_GRID,
           "evolution_times": [s ** 2 / 2 for s in SIGMA_GRID],
           "results": []}
    for sig in SIGMA_GRID:
        per_thr = {}
        for thr in HEAT_THRESHOLDS:
            sub = [r for r in rows if r["sigma_voxels"] == sig and r["threshold"] == thr]
            cov = np.array([r["coverage"] for r in sub if np.isfinite(r["coverage"])])
            per_thr[f"heat_ge_{thr}"] = {
                "n_rows": len(sub),
                "mean_pct": round(float(cov.mean() * 100), 2),
            }
        out["results"].append({
            "sigma_voxels": sig,
            "evolution_time_t": sig ** 2 / 2,
            **per_thr,
        })

    # Find optimal sigma per threshold
    print("\nResults summary:")
    for r in out["results"]:
        print(f"  sigma={r['sigma_voxels']} (t={r['evolution_time_t']}):  "
              f"heat>=0.50={r['heat_ge_0.5']['mean_pct']}%  "
              f"heat>=0.80={r['heat_ge_0.8']['mean_pct']}%")
    best_05 = max(out["results"], key=lambda r: r["heat_ge_0.5"]["mean_pct"])
    best_08 = max(out["results"], key=lambda r: r["heat_ge_0.8"]["mean_pct"])
    print(f"\n  Optimal sigma at heat>=0.50: {best_05['sigma_voxels']} ({best_05['heat_ge_0.5']['mean_pct']}%)")
    print(f"  Optimal sigma at heat>=0.80: {best_08['sigma_voxels']} ({best_08['heat_ge_0.8']['mean_pct']}%)")

    out["optimal_sigma_at_heat_05"] = best_05["sigma_voxels"]
    out["optimal_sigma_at_heat_08"] = best_08["sigma_voxels"]
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
