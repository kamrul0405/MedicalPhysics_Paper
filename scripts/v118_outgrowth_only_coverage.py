"""v118: Outgrowth-only coverage on PROTEAS-brain-mets.

Motivation. v117 revealed that the persistence baseline trivially
achieves 51.95% coverage at heat>=0.80 because ~52% of future-lesion
voxels are already in the baseline mask. The persistence prediction
is uninformative for the clinical question that actually matters
in radiation oncology: WHERE WILL NEW LESION APPEAR? — i.e., the
outgrowth voxels (future-lesion voxels OUTSIDE the baseline mask).

This script redefines the primary endpoint as outgrowth-only
coverage: the fraction of future-lesion-OUTGROWTH voxels covered
by each prior at heat>=0.50 and 0.80. The persistence baseline is
zero-by-construction on this metric (heat = baseline mask, so
heat AND NOT baseline = empty).

If anisotropic BED achieves significant outgrowth coverage where
persistence achieves 0%, this is a clean publishable result that
supports Proposal A as a high-impact-journal headline finding.

Outputs:
    source_data/v118_outgrowth_only_coverage.json
    source_data/v118_outgrowth_only_per_patient.csv
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
OUT_JSON = ROOT / "source_data" / "v118_outgrowth_only_coverage.json"
OUT_CSV = ROOT / "source_data" / "v118_outgrowth_only_per_patient.csv"

SIGMA_GRID = [0.5, 1.0, 2.5]
HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 10000
RNG = np.random.default_rng(11801)

SIGMA_PAR = 1.5
SIGMA_PERP = 4.0
ALPHA_BETA_TUMOUR = 10.0
BED_REF = 60.0


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
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def heat_kernel_iso_bed(mask, dose, fractions):
    """Isotropic BED-aware kernel (matching v94)."""
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    if fractions and fractions > 0:
        d = dose / fractions
        bed = dose * (1 + d / ALPHA_BETA_TUMOUR)
    else:
        bed = dose * (1 + dose / ALPHA_BETA_TUMOUR)
    bed_norm = np.clip(bed / BED_REF, 0.0, 2.0)
    sig_local = 1.5 + 1.0 * (1.0 - np.clip(bed_norm, 0, 1))
    sig_mean = float(np.mean(sig_local[mask]))
    h = gaussian_filter(mask.astype(np.float32), sigma=sig_mean)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def heat_kernel_aniso_bed(mask, dose, fractions):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    if fractions and fractions > 0:
        d = dose / fractions
        bed = dose * (1 + d / ALPHA_BETA_TUMOUR)
    else:
        bed = dose * (1 + dose / ALPHA_BETA_TUMOUR)
    bed_norm = np.clip(bed / BED_REF, 0.0, 2.0)

    gx = np.gradient(bed, axis=0)
    gy = np.gradient(bed, axis=1)
    gz = np.gradient(bed, axis=2)
    gmag = np.sqrt(gx**2 + gy**2 + gz**2) + 1e-9
    wx = np.abs(gx) / gmag
    wy = np.abs(gy) / gmag
    wz = np.abs(gz) / gmag

    h_par = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PAR)
    h_perp = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PERP)
    eps = 1e-9
    log_par = np.log(np.maximum(h_par, eps))
    log_perp = np.log(np.maximum(h_perp, eps))
    blend = (wx + wy + wz) / 3.0
    log_h = blend * log_par + (1.0 - blend) * log_perp
    h = np.exp(log_h)
    h = h * (1.0 + 0.1 * bed_norm)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    """Fraction of OUTGROWTH voxels (future and not baseline) in region_mask."""
    fut = future_mask.astype(bool)
    base = baseline_mask.astype(bool)
    outgrowth = fut & (~base)
    if outgrowth.sum() == 0:
        return float("nan")
    return float((outgrowth & region_mask.astype(bool)).sum() / outgrowth.sum())


def analyse_patient(outer, entry, clinical, work):
    pid = Path(entry.filename).stem
    row = clinical.loc[clinical["pid"] == pid]
    if row.empty:
        return []
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
            heat_maps = {}
            for sig in SIGMA_GRID:
                heat_maps[f"sigma_{sig}"] = heat_kernel_constant(base_mask, sig)
            heat_maps["persistence"] = base_mask.astype(np.float32)
            if dose is not None and dose.shape == base_mask.shape:
                heat_maps["iso_bed"] = heat_kernel_iso_bed(base_mask, dose, fractions)
                heat_maps["aniso_bed"] = heat_kernel_aniso_bed(base_mask, dose, fractions)
            else:
                heat_maps["iso_bed"] = None
                heat_maps["aniso_bed"] = None
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                outgrowth_mask = fu_mask & (~base_mask)
                outgrowth_count = int(outgrowth_mask.sum())
                row_data = {
                    "pid": pid, "fu_name": Path(fu_name).stem,
                    "outgrowth_voxels": outgrowth_count,
                }
                for thr in HEAT_THRESHOLDS:
                    for sig in SIGMA_GRID:
                        row_data[f"out_sigma_{sig}_thr_{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, heat_maps[f"sigma_{sig}"] >= thr)
                    row_data[f"out_persistence_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, heat_maps["persistence"] >= thr)
                    if heat_maps["iso_bed"] is not None:
                        row_data[f"out_iso_bed_thr_{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, heat_maps["iso_bed"] >= thr)
                    else:
                        row_data[f"out_iso_bed_thr_{thr}"] = float("nan")
                    if heat_maps["aniso_bed"] is not None:
                        row_data[f"out_aniso_bed_thr_{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, heat_maps["aniso_bed"] >= thr)
                    else:
                        row_data[f"out_aniso_bed_thr_{thr}"] = float("nan")
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
    print("v118 OUTGROWTH-ONLY COVERAGE on PROTEAS")
    print(f"  metric: outgrowth = future_mask AND NOT baseline_mask")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v118_") as td:
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

    print(f"\nTotal follow-up rows with outgrowth voxels: "
          f"{sum(1 for r in rows if r['outgrowth_voxels'] > 0)}")
    print(f"Total follow-up rows (incl. zero-outgrowth): {len(rows)}")
    if not rows:
        return

    # Filter to follow-ups with at least one outgrowth voxel
    rows_with_outgrowth = [r for r in rows if r["outgrowth_voxels"] > 0]
    print(f"Analysing only follow-ups with outgrowth_voxels > 0 (N = "
          f"{len(rows_with_outgrowth)})")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote per-patient CSV: {OUT_CSV}")

    pid_arr = np.array([r["pid"] for r in rows_with_outgrowth])
    out = {"version": "v118", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_total_followups": len(rows),
           "n_followups_with_outgrowth": len(rows_with_outgrowth),
           "n_patients": int(len(np.unique(pid_arr))),
           "n_bootstrap_replicates": N_BOOT, "alpha": 0.05,
           "thresholds": {}}

    print(f"\n{len(np.unique(pid_arr))} unique patients with at least one "
          f"outgrowth follow-up\n")

    for thr in HEAT_THRESHOLDS:
        thr_results = {"methods": {}, "paired_deltas_vs_aniso": {}}

        method_keys = [(f"sigma_{s}", f"out_sigma_{s}_thr_{thr}") for s in SIGMA_GRID] + [
            ("persistence", f"out_persistence_thr_{thr}"),
            ("iso_bed", f"out_iso_bed_thr_{thr}"),
            ("aniso_bed", f"out_aniso_bed_thr_{thr}"),
        ]

        print(f"--- heat >= {thr} ---")
        for label, col in method_keys:
            vals = np.array([r[col] for r in rows_with_outgrowth], dtype=float)
            if np.all(np.isnan(vals)):
                continue
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results["methods"][label] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            print(f"  {label:18s}: outgrowth coverage = {mean*100:5.2f}% "
                  f"[{lo*100:.2f}, {hi*100:.2f}]")
        print()

        # Paired delta vs aniso BED
        aniso_col = f"out_aniso_bed_thr_{thr}"
        aniso_vals = np.array([r[aniso_col] for r in rows_with_outgrowth], dtype=float)
        for label, col in method_keys:
            if label == "aniso_bed":
                continue
            base_vals = np.array([r[col] for r in rows_with_outgrowth], dtype=float)
            diff = aniso_vals - base_vals
            mean, lo, hi = cluster_bootstrap_ci(diff, pid_arr)
            thr_results["paired_deltas_vs_aniso"][f"aniso_minus_{label}"] = {
                "mean_pp": round(mean * 100, 2),
                "ci95_pp": [round(lo * 100, 2), round(hi * 100, 2)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            sig = "**SIG**" if (lo > 0 or hi < 0) else ""
            print(f"  delta aniso - {label:14s}: {mean*100:+6.2f} pp "
                  f"[{lo*100:+.2f}, {hi*100:+.2f}] {sig}")
        print()

        out["thresholds"][f"heat_ge_{thr}"] = thr_results

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
