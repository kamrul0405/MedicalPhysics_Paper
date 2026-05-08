"""v122: Ensemble prior heat = max(persistence, aniso_bed) on PROTEAS.

Motivation. v117 showed that the persistence baseline dominates at
heat>=0.80 (51.95% vs aniso 49.44%) while aniso BED dominates at
heat>=0.50 (52.84% vs persistence 51.87%). A natural clinically
deployable prior is the union: heat = max(persistence, aniso_bed).
This guarantees:
  - Recovery of all persistence voxels (heat=1.0 inside baseline mask).
  - Outgrowth-aware extension via the aniso BED kernel.

If the ensemble beats aniso BED alone, this is a publishable
clinically-actionable structural prior.

Outputs:
    source_data/v122_ensemble_prior.json
    source_data/v122_ensemble_prior_per_patient.csv
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
OUT_JSON = ROOT / "source_data" / "v122_ensemble_prior.json"
OUT_CSV = ROOT / "source_data" / "v122_ensemble_prior_per_patient.csv"

HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 10000
RNG = np.random.default_rng(12201)
SIGMA_PAR = 1.5
SIGMA_PERP = 4.0
ALPHA_BETA_TUMOUR = 10.0
BED_REF = 60.0


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


def heat_kernel_aniso_bed(mask, dose, fractions):
    if mask.sum() == 0:
        return np.zeros_like(mask, dtype=np.float32)
    if fractions and fractions > 0:
        d = dose / fractions
        bed = dose * (1 + d / ALPHA_BETA_TUMOUR)
    else:
        bed = dose * (1 + dose / ALPHA_BETA_TUMOUR)
    bed_norm = np.clip(bed / BED_REF, 0.0, 2.0)
    gx, gy, gz = np.gradient(bed, axis=0), np.gradient(bed, axis=1), np.gradient(bed, axis=2)
    gmag = np.sqrt(gx**2 + gy**2 + gz**2) + 1e-9
    blend = (np.abs(gx) + np.abs(gy) + np.abs(gz)) / (3.0 * gmag)
    h_par = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PAR)
    h_perp = gaussian_filter(mask.astype(np.float32), sigma=SIGMA_PERP)
    log_par = np.log(np.maximum(h_par, 1e-9))
    log_perp = np.log(np.maximum(h_perp, 1e-9))
    log_h = blend * log_par + (1.0 - blend) * log_perp
    h = np.exp(log_h) * (1.0 + 0.1 * bed_norm)
    if h.max() > 0: h = h / h.max()
    return h.astype(np.float32)


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


def analyse_patient(outer, entry, clinical, work):
    pid = Path(entry.filename).stem
    row = clinical.loc[clinical["pid"] == pid]
    if row.empty: return []
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
            if entries["baseline"] not in inner.namelist(): return []
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            try:
                dose = load_nii_from_inner(inner, entries["dose"], patient_tmp)
            except Exception:
                dose = None

            persistence = base_mask.astype(np.float32)  # heat = mask itself
            if dose is not None and dose.shape == base_mask.shape:
                aniso = heat_kernel_aniso_bed(base_mask, dose, fractions)
            else:
                aniso = None
            if aniso is not None:
                ensemble = np.maximum(persistence, aniso).astype(np.float32)
            else:
                ensemble = None

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
                    row_data[f"cov_pers_thr_{thr}"] = coverage(fu_mask, persistence >= thr)
                    if aniso is not None:
                        row_data[f"cov_aniso_thr_{thr}"] = coverage(fu_mask, aniso >= thr)
                    else:
                        row_data[f"cov_aniso_thr_{thr}"] = float("nan")
                    if ensemble is not None:
                        row_data[f"cov_ens_thr_{thr}"] = coverage(fu_mask, ensemble >= thr)
                    else:
                        row_data[f"cov_ens_thr_{thr}"] = float("nan")
                    # Outgrowth coverage
                    row_data[f"out_pers_thr_{thr}"] = outgrowth_coverage(
                        fu_mask, base_mask, persistence >= thr)
                    if aniso is not None:
                        row_data[f"out_aniso_thr_{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, aniso >= thr)
                    else:
                        row_data[f"out_aniso_thr_{thr}"] = float("nan")
                    if ensemble is not None:
                        row_data[f"out_ens_thr_{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, ensemble >= thr)
                    else:
                        row_data[f"out_ens_thr_{thr}"] = float("nan")
                out_rows.append(row_data)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
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
    print("v122 ENSEMBLE PRIOR max(persistence, aniso_bed) on PROTEAS")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)
    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v122_") as td:
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
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"Wrote {OUT_CSV}")

    pid_arr = np.array([r["pid"] for r in rows])
    out = {"version": "v122", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_followups": len(rows),
           "n_patients": int(len(np.unique(pid_arr))),
           "n_bootstrap_replicates": N_BOOT, "alpha": 0.05,
           "thresholds": {}}

    for thr in HEAT_THRESHOLDS:
        thr_results = {}
        print(f"\n--- heat >= {thr} (overall future-lesion coverage) ---")
        cols = [("persistence", f"cov_pers_thr_{thr}"),
                ("aniso_bed", f"cov_aniso_thr_{thr}"),
                ("ensemble", f"cov_ens_thr_{thr}")]
        for label, col in cols:
            vals = np.array([r[col] for r in rows], dtype=float)
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results[f"overall_{label}"] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            print(f"  overall {label:12s}: {mean*100:5.2f}% [{lo*100:.2f}, {hi*100:.2f}]")

        # Paired delta ensemble vs aniso, ensemble vs persistence
        ens_vals = np.array([r[f"cov_ens_thr_{thr}"] for r in rows], dtype=float)
        for label, col in cols[:-1]:
            base_vals = np.array([r[col] for r in rows], dtype=float)
            diff = ens_vals - base_vals
            mean, lo, hi = cluster_bootstrap_ci(diff, pid_arr)
            thr_results[f"overall_delta_ensemble_minus_{label}"] = {
                "mean_pp": round(mean * 100, 2),
                "ci95_pp": [round(lo * 100, 2), round(hi * 100, 2)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            sig = "**SIG**" if (lo > 0 or hi < 0) else ""
            print(f"  delta ensemble - {label:12s}: {mean*100:+6.2f} pp "
                  f"[{lo*100:+.2f}, {hi*100:+.2f}] {sig}")

        # Outgrowth-only
        print(f"\n--- heat >= {thr} (outgrowth-only coverage) ---")
        out_cols = [("persistence", f"out_pers_thr_{thr}"),
                    ("aniso_bed", f"out_aniso_thr_{thr}"),
                    ("ensemble", f"out_ens_thr_{thr}")]
        for label, col in out_cols:
            vals = np.array([r[col] for r in rows], dtype=float)
            valid = ~np.isnan(vals)
            if valid.sum() == 0: continue
            mean, lo, hi = cluster_bootstrap_ci(vals[valid], pid_arr[valid])
            thr_results[f"outgrowth_{label}"] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
                "n_with_outgrowth": int(valid.sum()),
            }
            print(f"  outgrowth {label:12s}: {mean*100:5.2f}% [{lo*100:.2f}, {hi*100:.2f}]"
                  f" (N={valid.sum()})")

        out["thresholds"][f"heat_ge_{thr}"] = thr_results

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
