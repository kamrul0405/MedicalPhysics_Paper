"""v101: Anisotropic BED kernel σ_par/σ_perp sensitivity sweep.

Robustness check for the v98 anisotropic BED-aware kernel breakthrough
(+19.3 pp coverage gain at heat ≥ 0.80). Tests four (σ_par, σ_perp)
combinations:
    (1.0, 3.5)
    (1.0, 4.0)
    (2.0, 4.0)
    (2.0, 4.5)
to verify the +12.31 pp gain over isotropic BED is not a cherry-picked
parameter setting.

Outputs:
    source_data/v101_anisotropic_sensitivity_sweep.json
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
OUT_JSON = ROOT / "source_data" / "v101_anisotropic_sensitivity_sweep.json"

ALPHA_BETA = 10.0
HEAT_THRESHOLDS = [0.5, 0.8]

SIGMA_PAIRS = [
    (1.0, 3.5),
    (1.0, 4.0),
    (1.5, 4.0),  # v98 baseline
    (2.0, 4.0),
    (2.0, 4.5),
]


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


def anisotropic_bed_heat(mask, dose, n_fx, sigma_par, sigma_perp, ab=ALPHA_BETA):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    gx, gy, gz = np.gradient(dose.astype(np.float32))
    g_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2) + 1e-6
    w_x = np.abs(gx) / g_mag
    w_y = np.abs(gy) / g_mag
    w_z = np.abs(gz) / g_mag
    m = mask.astype(np.float32)
    h_par_x = gaussian_filter(m, sigma=(sigma_par, 0, 0))
    h_perp_x = gaussian_filter(m, sigma=(sigma_perp, 0, 0))
    h_par_y = gaussian_filter(m, sigma=(0, sigma_par, 0))
    h_perp_y = gaussian_filter(m, sigma=(0, sigma_perp, 0))
    h_par_z = gaussian_filter(m, sigma=(0, 0, sigma_par))
    h_perp_z = gaussian_filter(m, sigma=(0, 0, sigma_perp))
    h_x = w_x * h_par_x + (1 - w_x) * h_perp_x
    h_y = w_y * h_par_y + (1 - w_y) * h_perp_y
    h_z = w_z * h_par_z + (1 - w_z) * h_perp_z
    heat = np.cbrt(np.maximum(h_x, 1e-6) * np.maximum(h_y, 1e-6) * np.maximum(h_z, 1e-6))
    dose_per_fx = dose / max(n_fx, 1)
    bed = n_fx * dose_per_fx * (1 + dose_per_fx / ab)
    bed_max = bed.max() if bed.max() > 0 else 1.0
    bed_norm = np.clip(bed / bed_max, 0.0, 1.0)
    heat = heat * (1 + 0.5 * bed_norm)
    if heat.max() > 0:
        heat = heat / heat.max()
    return heat.astype(np.float32)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0:
        return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def analyse_patient(outer, entry, clinical, work, sigma_pairs):
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
            heat_maps = {}
            for sp, spp in sigma_pairs:
                heat_maps[(sp, spp)] = anisotropic_bed_heat(base_mask, dose, n_fx, sp, spp)
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                for sp, spp in sigma_pairs:
                    for thr in HEAT_THRESHOLDS:
                        out_rows.append({
                            "pid": pid,
                            "sigma_par": sp,
                            "sigma_perp": spp,
                            "threshold": thr,
                            "coverage": coverage(fu_mask, heat_maps[(sp, spp)] >= thr),
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
    print("v101 ANISOTROPIC BED-aware sigma_par/sigma_perp SENSITIVITY SWEEP")
    print(f"  testing {len(SIGMA_PAIRS)} parameter combinations")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v101_") as td:
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
                pr = analyse_patient(outer, entry, clinical, work, SIGMA_PAIRS)
                rows.extend(pr)
                del pr
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} processed  ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal rows: {len(rows)}")

    out = {"version": "v101", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "sigma_pairs": [{"sigma_par": sp, "sigma_perp": spp} for sp, spp in SIGMA_PAIRS],
           "results": []}
    for sp, spp in SIGMA_PAIRS:
        per_thr = {}
        for thr in HEAT_THRESHOLDS:
            sub = [r for r in rows if r["sigma_par"] == sp and r["sigma_perp"] == spp and r["threshold"] == thr]
            cov = np.array([r["coverage"] for r in sub if np.isfinite(r["coverage"])])
            per_thr[f"heat_ge_{thr}"] = {
                "n_rows": len(sub),
                "mean_pct": round(float(cov.mean() * 100), 2),
            }
        out["results"].append({"sigma_par": sp, "sigma_perp": spp, **per_thr})
        print(f"\n  sigma_par={sp}, sigma_perp={spp}:")
        print(f"    heat>=0.50: {per_thr['heat_ge_0.5']['mean_pct']}%")
        print(f"    heat>=0.80: {per_thr['heat_ge_0.8']['mean_pct']}%")

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
