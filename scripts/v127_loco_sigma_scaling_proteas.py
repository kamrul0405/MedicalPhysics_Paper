"""v127: LOCO validation of the v124 sigma scaling law on PROTEAS-brain-mets.

The v124 mixed-effects model was fit on UCSF, MU, RHUH, LUMIERE
(N = 505 patient observations). PROTEAS-brain-mets was NOT in the
training set; v109 only saw it at sigma >= 1.0, so it has no
sub-voxel optimum.

This script:
  (a) Computes per-patient sigma_opt on PROTEAS-brain-mets at the
      v124 grid {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}
      at heat >= 0.50 (and 0.80 as sanity check).
  (b) Predicts each patient's sigma_opt from the v124 formula:
      sigma_pred = exp(-3.0944) * r_eq^1.2733
  (c) Compares predicted vs actual: R^2, RMSE, paired correlation.
  (d) Tests coverage of v124's 95% CI for the held-out cohort.

If the v124 formula generalises to PROTEAS within reasonable error,
this is the bulletproof LOCO validation that makes Proposal H
high-impact-publishable.

Outputs:
    source_data/v127_loco_sigma_scaling_proteas.json
    source_data/v127_loco_sigma_scaling_per_patient.csv
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
from scipy.stats import pearsonr

ROOT = Path(__file__).parent.parent
DATA_ZIP = Path(r"C:\Users\kamru\Downloads\Datasets\PKG - PROTEAS-brain-mets-zenodo-17253793.zip")
OUT_JSON = ROOT / "source_data" / "v127_loco_sigma_scaling_proteas.json"
OUT_CSV = ROOT / "source_data" / "v127_loco_sigma_scaling_per_patient.csv"

# v124 fitted formula at heat >= 0.50:
V124_INTERCEPT = -3.0944
V124_SLOPE = 1.2733
V124_INTERCEPT_CI = [-3.5, -2.7]   # rough; use exact from json
V124_SLOPE_CI = [1.158, 1.389]

SIGMA_GRID = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
HEAT_THRESHOLDS = [0.5, 0.8]


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


def heat_kernel(mask, sigma):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0:
        h = h / h.max()
    return h.astype(np.float32)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def equivalent_radius(mask):
    n = float(mask.sum())
    if n == 0:
        return 0.0
    return float((3 * n / (4 * np.pi)) ** (1 / 3))


def predict_sigma_v124(r_eq, intercept=V124_INTERCEPT, slope=V124_SLOPE):
    return float(np.exp(intercept + slope * np.log(max(r_eq, 1.0))))


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
            sigma_pred = predict_sigma_v124(r_eq)

            heat_maps = {sig: heat_kernel(base_mask, sig) for sig in SIGMA_GRID}

            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any():
                    continue
                row_data = {"pid": pid, "fu_name": Path(fu_name).stem,
                            "r_eq": r_eq, "log_r_eq": float(np.log(r_eq)),
                            "sigma_pred_v124": sigma_pred}
                # Per-fu sigma_opt at each threshold
                for thr in HEAT_THRESHOLDS:
                    best_sig = None
                    best_cov = -1.0
                    for sig in SIGMA_GRID:
                        cov = coverage(fu_mask, heat_maps[sig] >= thr)
                        if np.isnan(cov):
                            continue
                        if cov > best_cov:
                            best_cov = cov
                            best_sig = sig
                    row_data[f"sigma_opt_actual_thr_{thr}"] = best_sig
                    row_data[f"best_coverage_thr_{thr}"] = best_cov
                    if best_sig is not None and best_sig > 0:
                        row_data[f"log_sigma_opt_thr_{thr}"] = float(np.log(best_sig))
                    else:
                        row_data[f"log_sigma_opt_thr_{thr}"] = None
                out_rows.append(row_data)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out_rows


def main():
    print("=" * 78)
    print("v127 LOCO VALIDATION OF v124 SIGMA SCALING LAW ON PROTEAS-BRAIN-METS")
    print(f"  v124 formula: sigma_pred = exp({V124_INTERCEPT}) * r_eq^{V124_SLOPE}")
    print(f"  v124 slope 95% CI: {V124_SLOPE_CI}")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v127_") as td:
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

    print(f"\nTotal follow-ups with data: {len(rows)}")
    if not rows:
        return

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote per-follow-up CSV: {OUT_CSV}")

    # === Validation analysis at heat >= 0.50 ===
    out = {"version": "v127", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "v124_formula": {"intercept": V124_INTERCEPT, "slope": V124_SLOPE,
                              "slope_95ci": V124_SLOPE_CI},
           "n_followups": len(rows), "thresholds": {}}

    for thr in HEAT_THRESHOLDS:
        print(f"\n=== heat >= {thr} ===")
        # Filter rows with valid sigma_opt
        valid = [r for r in rows
                 if r[f"sigma_opt_actual_thr_{thr}"] is not None
                 and r[f"sigma_opt_actual_thr_{thr}"] > 0]
        if len(valid) < 10:
            print(f"  too few valid rows ({len(valid)}); skipping")
            continue

        sigma_actual = np.array([r[f"sigma_opt_actual_thr_{thr}"] for r in valid])
        sigma_pred = np.array([r["sigma_pred_v124"] for r in valid])
        r_eq = np.array([r["r_eq"] for r in valid])
        log_sigma_actual = np.log(sigma_actual)
        log_sigma_pred = np.log(sigma_pred)

        rmse_log = float(np.sqrt(np.mean((log_sigma_actual - log_sigma_pred) ** 2)))
        mae_log = float(np.mean(np.abs(log_sigma_actual - log_sigma_pred)))
        # R^2 of predicted vs actual on log scale
        ss_res = np.sum((log_sigma_actual - log_sigma_pred) ** 2)
        ss_tot = np.sum((log_sigma_actual - log_sigma_actual.mean()) ** 2)
        r2 = float(1 - ss_res / max(ss_tot, 1e-9))
        # Pearson correlation between log_sigma_actual and log_r_eq (within PROTEAS)
        r_pearson_log_r, p_pearson_log_r = pearsonr(np.log(r_eq), log_sigma_actual)
        # Pearson correlation between actual and predicted
        r_pearson_pred, p_pearson_pred = pearsonr(log_sigma_pred, log_sigma_actual)

        # Refit a within-PROTEAS slope and compare
        x = np.log(r_eq)
        y = log_sigma_actual
        n = len(x)
        x_mean = x.mean(); y_mean = y.mean()
        slope_proteas = float(np.sum((x - x_mean) * (y - y_mean)) /
                                max(np.sum((x - x_mean) ** 2), 1e-9))
        intercept_proteas = float(y_mean - slope_proteas * x_mean)
        # SE on slope
        residuals = y - (intercept_proteas + slope_proteas * x)
        sse = float(np.sum(residuals ** 2))
        se_slope = float(np.sqrt(sse / (n - 2) / max(np.sum((x - x_mean) ** 2), 1e-9)))
        slope_proteas_ci = [slope_proteas - 1.96 * se_slope,
                              slope_proteas + 1.96 * se_slope]

        # Test if v124 slope (1.273) is within PROTEAS slope's 95% CI
        v124_slope_in_proteas_ci = (slope_proteas_ci[0] <= V124_SLOPE
                                      <= slope_proteas_ci[1])

        # Sigma_opt distribution counts
        counts = {s: int((sigma_actual == s).sum()) for s in SIGMA_GRID}
        unique_count = sum(1 for v in counts.values() if v > 0)

        thr_results = {
            "n_followups": int(len(valid)),
            "rmse_log_sigma": rmse_log,
            "mae_log_sigma": mae_log,
            "r_squared": r2,
            "pearson_r_log_actual_vs_log_r_eq": float(r_pearson_log_r),
            "pearson_p_log_actual_vs_log_r_eq": float(p_pearson_log_r),
            "pearson_r_pred_vs_actual": float(r_pearson_pred),
            "proteas_within_cohort_slope": slope_proteas,
            "proteas_within_cohort_slope_se": se_slope,
            "proteas_within_cohort_slope_ci": slope_proteas_ci,
            "proteas_within_cohort_intercept": intercept_proteas,
            "v124_slope_in_proteas_ci": v124_slope_in_proteas_ci,
            "sigma_opt_distribution": counts,
            "n_unique_sigma_opt": unique_count,
            "median_r_eq": float(np.median(r_eq)),
        }
        out["thresholds"][f"heat_ge_{thr}"] = thr_results

        print(f"  N follow-ups with valid optimum: {len(valid)}")
        print(f"  Median r_eq: {np.median(r_eq):.2f} voxels")
        print(f"  Sigma_opt distribution: {counts}")
        print(f"  Unique sigma_opt values: {unique_count}")
        print(f"  RMSE(log sigma): {rmse_log:.3f}")
        print(f"  MAE(log sigma): {mae_log:.3f}")
        print(f"  R^2 (predicted vs actual): {r2:+.3f}")
        print(f"  Pearson r (log r_eq vs log sigma_actual): "
              f"{r_pearson_log_r:+.3f} (p = {p_pearson_log_r:.4f})")
        print(f"  Pearson r (predicted vs actual): "
              f"{r_pearson_pred:+.3f} (p = {p_pearson_pred:.4f})")
        print(f"  PROTEAS within-cohort slope: {slope_proteas:+.3f} +/- {se_slope:.3f}")
        print(f"  PROTEAS within-cohort slope 95% CI: "
              f"[{slope_proteas_ci[0]:+.3f}, {slope_proteas_ci[1]:+.3f}]")
        print(f"  v124 slope ({V124_SLOPE}) within PROTEAS CI? "
              f"{v124_slope_in_proteas_ci}")

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}")


if __name__ == "__main__":
    main()
