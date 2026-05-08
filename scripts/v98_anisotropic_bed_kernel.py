"""v98: Anisotropic BED-aware structural prior — directional σ along the
principal direction of the local dose gradient.

Extends the isotropic BED-aware kernel (v94; +6.99 pp at heat ≥ 0.80) with
direction-aware σ. The local dose-gradient tensor at each voxel encodes the
direction along which the dose changes most rapidly; the structural prior's
σ along that direction should be smaller (tighter prior in the high-gradient
direction) and larger orthogonally (broader prior where dose is more uniform).

For each voxel, σ_par along the gradient direction = σ_lo;
σ_perp orthogonal to the gradient = σ_hi.

This is implemented via direction-aware Gaussian filtering using the
gradient direction tensor. To make this computationally feasible we
approximate the anisotropic kernel by separable 3D Gaussian filtering
along the local gradient direction with σ_par and orthogonal directions
with σ_perp.

Outputs:
  source_data/v98_anisotropic_bed_kernel.json
  source_data/v98_anisotropic_bed_per_patient.csv
"""
from __future__ import annotations

import csv
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
OUT_JSON = ROOT / "source_data" / "v98_anisotropic_bed_kernel.json"
OUT_CSV = ROOT / "source_data" / "v98_anisotropic_bed_per_patient.csv"

SIGMA_PAR = 1.5  # along dose-gradient direction (tight)
SIGMA_PERP = 4.0  # orthogonal to gradient (broad)
SIGMA_BASE = 2.5
ALPHA_BETA_TUMOUR = 10.0
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


def isotropic_bed_heat(mask, dose_total, n_fx, ab=ALPHA_BETA_TUMOUR):
    """v94's BED-aware isotropic kernel — for direct comparison."""
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    h_lo = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PAR)
    h_hi = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PERP)
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


def anisotropic_bed_heat(mask, dose_total, n_fx, ab=ALPHA_BETA_TUMOUR):
    """Anisotropic BED-aware kernel — σ along the gradient direction differs
    from σ orthogonal. Approximated via separable 3D Gaussian filtering at
    different σ values along principal axes scaled by the local dose-gradient
    magnitude.

    Key idea: high local |∇D| ⇒ tighter prior along the high-gradient axis;
    low local |∇D| ⇒ broader prior. Approximate by:
      σ_x(x) = σ_par + (σ_perp − σ_par) × (1 − |∂D/∂x|/|∇D|)
    and similarly for y, z.

    Implementation: pre-compute heat at σ_par and σ_perp along each axis
    (3 axis-specific filters × 2 σ values = 6 filtered arrays), then blend
    pixel-wise based on per-axis gradient-magnitude weights.
    """
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    # Compute per-axis dose gradient magnitudes
    gx, gy, gz = np.gradient(dose_total.astype(np.float32))
    g_mag = np.sqrt(gx ** 2 + gy ** 2 + gz ** 2) + 1e-6
    # Per-axis weights w_axis = |g_axis| / |g| ∈ [0, 1]; high w_axis → tight σ on that axis
    w_x = np.abs(gx) / g_mag
    w_y = np.abs(gy) / g_mag
    w_z = np.abs(gz) / g_mag
    # σ per axis: σ_x(x) = σ_par × w_x + σ_perp × (1 − w_x)
    # Approximate via separable filtering at two σ values and blend
    m = mask.astype(np.float32)
    # Filter along each axis separately at two σ levels
    h_par_x = gaussian_filter(m, sigma=(SIGMA_PAR, 0, 0))
    h_perp_x = gaussian_filter(m, sigma=(SIGMA_PERP, 0, 0))
    h_par_y = gaussian_filter(m, sigma=(0, SIGMA_PAR, 0))
    h_perp_y = gaussian_filter(m, sigma=(0, SIGMA_PERP, 0))
    h_par_z = gaussian_filter(m, sigma=(0, 0, SIGMA_PAR))
    h_perp_z = gaussian_filter(m, sigma=(0, 0, SIGMA_PERP))
    # Blend per-axis based on gradient direction
    h_x = w_x * h_par_x + (1 - w_x) * h_perp_x
    h_y = w_y * h_par_y + (1 - w_y) * h_perp_y
    h_z = w_z * h_par_z + (1 - w_z) * h_perp_z
    # Combine via geometric mean (so all three axes contribute)
    heat = np.cbrt(np.maximum(h_x, 1e-6) * np.maximum(h_y, 1e-6) * np.maximum(h_z, 1e-6))
    # BED weighting on top
    dose_per_fx = dose_total / max(n_fx, 1)
    bed = n_fx * dose_per_fx * (1 + dose_per_fx / ab)
    bed_max = bed.max() if bed.max() > 0 else 1.0
    bed_norm = np.clip(bed / bed_max, 0.0, 1.0)
    heat = heat * (1 + 0.5 * bed_norm)  # mild BED amplification of high-gradient regions
    if heat.max() > 0:
        heat = heat / heat.max()
    return heat.astype(np.float32)


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
            heat_iso_bed = isotropic_bed_heat(base_mask, dose, n_fx)
            heat_aniso = anisotropic_bed_heat(base_mask, dose, n_fx)
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                for thr in HEAT_THRESHOLDS:
                    out_rows.append({
                        "pid": pid,
                        "rx_gy": float(rx),
                        "fractions": int(n_fx),
                        "threshold": thr,
                        "coverage_const_sigma": coverage(fu_mask, heat_const >= thr),
                        "coverage_iso_bed": coverage(fu_mask, heat_iso_bed >= thr),
                        "coverage_aniso_bed": coverage(fu_mask, heat_aniso >= thr),
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
    print("v98 ANISOTROPIC BED-AWARE STRUCTURAL PRIOR")
    print(f"  sigma_par={SIGMA_PAR}, sigma_perp={SIGMA_PERP}, alpha/beta={ALPHA_BETA_TUMOUR}")
    print("=" * 78)
    if not DATA_ZIP.exists():
        raise FileNotFoundError(DATA_ZIP)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v98_") as td:
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
    if not rows:
        return

    out = {"version": "v98", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_rows_total": len(rows), "thresholds": {}}
    for thr in HEAT_THRESHOLDS:
        sub = [r for r in rows if r["threshold"] == thr]
        cov_const = np.array([r["coverage_const_sigma"] for r in sub if np.isfinite(r["coverage_const_sigma"])])
        cov_iso = np.array([r["coverage_iso_bed"] for r in sub if np.isfinite(r["coverage_iso_bed"])])
        cov_aniso = np.array([r["coverage_aniso_bed"] for r in sub if np.isfinite(r["coverage_aniso_bed"])])
        block = {
            "n_rows": len(sub),
            "constant_sigma_mean_pct": round(float(cov_const.mean() * 100), 2),
            "isotropic_bed_mean_pct": round(float(cov_iso.mean() * 100), 2),
            "anisotropic_bed_mean_pct": round(float(cov_aniso.mean() * 100), 2),
            "delta_aniso_minus_iso_pp": round(float((cov_aniso.mean() - cov_iso.mean()) * 100), 2),
            "delta_aniso_minus_const_pp": round(float((cov_aniso.mean() - cov_const.mean()) * 100), 2),
            "n_patients": int(len({r["pid"] for r in sub})),
        }
        out["thresholds"][f"heat_ge_{thr}"] = block
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT_JSON}")
    print(json.dumps(out["thresholds"], indent=2))

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
