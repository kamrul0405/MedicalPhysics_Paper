"""v116: Fair comparison of v98 anisotropic BED-aware kernel against
the LESION-PERSISTENCE BASELINE on PROTEAS-brain-mets.

v115 revealed that sub-voxel sigma (sigma=0.25) wins at heat>=0.80
on all 4 cache_3d cohorts because the kernel collapses to ~the
binary mask itself, recovering lesion persistence. This raises a
fair-comparison concern for the v98 +12.33 pp anisotropic gain:
does the anisotropic kernel actually outperform the persistence
baseline, or is its gain partly explained by sigma rescaling?

This script computes future-lesion coverage on PROTEAS at:
  - sigma = 2.5 (legacy default)
  - sigma = 1.0 (v109 optimum from sigma>=1 grid)
  - sigma = 0.5 (sub-voxel)
  - sigma = 0.25 (sub-voxel; near-mask)
  - Pure lesion persistence (heat = mask itself)
  - v98 anisotropic BED-aware kernel

Then computes 10,000-replicate cluster-bootstrap CIs on the paired
delta of (anisotropic - persistence) and (anisotropic - sigma=1.0).
If the anisotropic CIs exclude zero against persistence, the v98
finding is genuinely bulletproof. If not, we report this honestly
and reframe the contribution.

Outputs:
    source_data/v116_anisotropic_vs_persistence.json
    source_data/v116_anisotropic_vs_persistence_per_patient.csv
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
OUT_JSON = ROOT / "source_data" / "v116_anisotropic_vs_persistence.json"
OUT_CSV = ROOT / "source_data" / "v116_anisotropic_vs_persistence_per_patient.csv"

SIGMA_GRID = [0.25, 0.5, 1.0, 2.5]
HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 10000
RNG = np.random.default_rng(11601)

# Anisotropic BED-aware kernel parameters (matching v98)
SIGMA_PAR = 1.5
SIGMA_PERP = 4.0
ALPHA_BETA_TUMOUR = 10.0  # Gy
BED_REF = 60.0  # Gy reference for normalisation


def parse_rx(value):
    if pd.isna(value):
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", str(value))
    return float(match.group(1)) if match else None


def parse_fractions(value):
    if pd.isna(value):
        return None
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


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
    dose = f"{prefix}{pid}_RTP.nii.gz"
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
    return {"dose": dose, "baseline": baseline, "followups": followups}


def heat_kernel_constant(mask, sigma):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    if sigma <= 0:
        h = mask.astype(np.float32)
    else:
        h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def heat_kernel_persistence(mask):
    """The 'persistence baseline' is just the binary mask itself,
    treated as a 0/1 heat map. heat >= 0.5 == heat >= 0.8 == mask."""
    return mask.astype(np.float32)


def heat_kernel_anisotropic_bed(mask, dose, rx_dose, fractions):
    """Anisotropic BED-aware kernel matching v98."""
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    # BED computation
    if fractions and fractions > 0:
        d = dose / fractions  # dose per fraction
        bed = dose * (1 + d / ALPHA_BETA_TUMOUR)
    else:
        bed = dose * (1 + dose / ALPHA_BETA_TUMOUR)  # SRS single fraction
    bed_norm = np.clip(bed / BED_REF, 0.0, 2.0)

    # Per-axis gradient magnitudes
    gx = np.gradient(bed, axis=0)
    gy = np.gradient(bed, axis=1)
    gz = np.gradient(bed, axis=2)
    gmag = np.sqrt(gx**2 + gy**2 + gz**2) + 1e-9
    wx = np.abs(gx) / gmag
    wy = np.abs(gy) / gmag
    wz = np.abs(gz) / gmag

    # Per-axis filtered components
    h_par = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PAR)
    h_perp = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PERP)

    # Geometric blend along each axis
    eps = 1e-9
    log_par = np.log(np.maximum(h_par, eps))
    log_perp = np.log(np.maximum(h_perp, eps))

    # Combine: per-axis weights determine local par/perp blend
    log_h = (wx + wy + wz) / 3.0 * log_par + (1.0 - (wx + wy + wz) / 3.0) * log_perp
    h = np.exp(log_h)

    # Mild BED amplification
    h = h * (1.0 + 0.1 * bed_norm)

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
    rx_dose = row.iloc[0]["rx_gy"]
    fractions = row.iloc[0]["fractions_num"]
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
            try:
                dose = load_nii_from_inner(inner, entries["dose"], patient_tmp)
            except Exception:
                dose = None

            # Pre-compute heat maps for the comparison set
            heat_maps = {}
            for sig in SIGMA_GRID:
                heat_maps[f"sigma_{sig}"] = heat_kernel_constant(base_mask, sig)
            heat_maps["persistence"] = heat_kernel_persistence(base_mask)
            if dose is not None and dose.shape == base_mask.shape:
                heat_maps["aniso_bed"] = heat_kernel_anisotropic_bed(
                    base_mask, dose, rx_dose, fractions)
            else:
                heat_maps["aniso_bed"] = None

            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                row_data = {"pid": pid, "fu_name": Path(fu_name).stem}
                for thr in HEAT_THRESHOLDS:
                    for sig in SIGMA_GRID:
                        row_data[f"cov_sigma_{sig}_thr_{thr}"] = coverage(
                            fu_mask, heat_maps[f"sigma_{sig}"] >= thr)
                    row_data[f"cov_persistence_thr_{thr}"] = coverage(
                        fu_mask, heat_maps["persistence"] >= thr)
                    if heat_maps["aniso_bed"] is not None:
                        row_data[f"cov_aniso_bed_thr_{thr}"] = coverage(
                            fu_mask, heat_maps["aniso_bed"] >= thr)
                    else:
                        row_data[f"cov_aniso_bed_thr_{thr}"] = float("nan")
                out_rows.append(row_data)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try:
            nested_path.unlink()
        except OSError:
            pass
    return out_rows


def cluster_bootstrap_ci(values, pids, alpha=0.05):
    pids_unique = np.unique(pids)
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sample = RNG.choice(pids_unique, size=len(pids_unique), replace=True)
        vals = []
        for s in sample:
            mask = pids == s
            vals.extend(values[mask].tolist())
        boots[b] = np.nanmean(vals) if vals else np.nan
    lo = np.nanpercentile(boots, 100 * alpha / 2)
    hi = np.nanpercentile(boots, 100 * (1 - alpha / 2))
    return float(np.nanmean(boots)), float(lo), float(hi)


def main():
    print("=" * 78)
    print("v116 ANISOTROPIC BED vs LESION-PERSISTENCE BASELINE on PROTEAS")
    print(f"  comparison set: sigma in {SIGMA_GRID}, persistence, anisotropic BED")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v116_") as td:
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
                    print(f"  {i}/{len(entries)} processed  ({time.time()-t0:.0f}s)",
                          flush=True)

    print(f"\nTotal follow-up rows: {len(rows)}")
    if not rows:
        return

    # Save per-patient CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote per-patient CSV: {OUT_CSV}")

    # Build cluster-bootstrap CIs for each method at each threshold
    pid_arr = np.array([r["pid"] for r in rows])
    out = {"version": "v116", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_rows": len(rows), "n_patients": int(len(np.unique(pid_arr))),
           "n_bootstrap_replicates": N_BOOT, "alpha": 0.05,
           "thresholds": {}}

    for thr in HEAT_THRESHOLDS:
        thr_results = {"methods": {}, "paired_deltas": {}}
        method_cols = (
            [(f"sigma_{s}", f"cov_sigma_{s}_thr_{thr}") for s in SIGMA_GRID]
            + [("persistence", f"cov_persistence_thr_{thr}"),
               ("anisotropic_bed", f"cov_aniso_bed_thr_{thr}")]
        )
        # Per-method coverage
        for label, col in method_cols:
            vals = np.array([r[col] for r in rows], dtype=float)
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results["methods"][label] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            print(f"  thr {thr}  {label:18s}: {mean*100:5.2f}% "
                  f"[{lo*100:.2f}, {hi*100:.2f}]")

        # Paired-delta CIs: anisotropic vs each baseline
        aniso_col = f"cov_aniso_bed_thr_{thr}"
        aniso_vals = np.array([r[aniso_col] for r in rows], dtype=float)
        for label, col in method_cols:
            if label == "anisotropic_bed":
                continue
            base_vals = np.array([r[col] for r in rows], dtype=float)
            diff = aniso_vals - base_vals
            mean, lo, hi = cluster_bootstrap_ci(diff, pid_arr)
            thr_results["paired_deltas"][f"aniso_minus_{label}"] = {
                "mean_pp": round(mean * 100, 2),
                "ci95_pp": [round(lo * 100, 2), round(hi * 100, 2)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            sig = "**" if (lo > 0 or hi < 0) else "  "
            print(f"  thr {thr}  delta vs {label:14s}: {mean*100:+6.2f} pp "
                  f"[{lo*100:+.2f}, {hi*100:+.2f}] {sig}")
        out["thresholds"][f"heat_ge_{thr}"] = thr_results
        print()

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
