"""v142: Time-stratified bimodal kernel coverage on PROTEAS.

Stratify PROTEAS follow-ups by time-since-baseline (early/mid/late
based on follow-up index 1, 2, 3+) and test whether the bimodal
kernel's advantage over persistence is stable, increases, or decays
with follow-up time.

The PROTEAS dataset names follow-ups as fu1, fu2, fu3, ... in
chronological order. We treat the FU index as a proxy for time
(since exact dates are not always available).

Stratification:
  - Early: fu1 (first follow-up after baseline)
  - Mid: fu2 (second follow-up)
  - Late: fu3+ (third or later)

For each stratum, compute bimodal vs persistence overall + outgrowth
coverage with cluster-bootstrap CIs.

This is the temporal-deployment-window analysis required for clinical
journals.

Outputs:
    source_data/v142_time_stratified_bimodal.json
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
OUT_JSON = ROOT / "source_data" / "v142_time_stratified_bimodal.json"

SIGMA_BROAD = 7.0
HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 5000
RNG = np.random.default_rng(14201)


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


def overall_coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    fut = future_mask.astype(bool); base = baseline_mask.astype(bool)
    out = fut & (~base)
    if out.sum() == 0: return float("nan")
    return float((out & region_mask.astype(bool)).sum() / out.sum())


def vectorised_cluster_bootstrap(values, alpha=0.05):
    values = np.asarray(values, dtype=float)
    valid = ~np.isnan(values)
    if valid.sum() == 0:
        return float("nan"), float("nan"), float("nan")
    v = values[valid]
    n = len(v)
    sample_idx = RNG.integers(0, n, size=(N_BOOT, n))
    boot_means = v[sample_idx].mean(axis=1)
    lo = float(np.percentile(boot_means, 100 * alpha / 2))
    hi = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    return float(boot_means.mean()), lo, hi


def parse_fu_index(fu_name: str) -> int:
    """Extract the integer follow-up index from names like
    'P01_tumor_mask_fu1', 'P03_tumor_mask_fu12', etc."""
    m = re.search(r"_fu(\d+)", fu_name)
    return int(m.group(1)) if m else 0


def analyse_patient(outer, entry, work):
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
            if entries["baseline"] not in inner.namelist(): return []
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            if base_mask.sum() == 0: return []
            persistence = base_mask.astype(np.float32)
            heat = heat_bimodal(base_mask, SIGMA_BROAD)
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any(): continue
                fu_index = parse_fu_index(Path(fu_name).stem)
                row = {"pid": pid, "fu_name": Path(fu_name).stem,
                       "fu_index": fu_index}
                for thr in HEAT_THRESHOLDS:
                    row[f"cov_pers_thr_{thr}"] = overall_coverage(fu_mask, persistence >= thr)
                    row[f"cov_bim_thr_{thr}"] = overall_coverage(fu_mask, heat >= thr)
                    row[f"out_pers_thr_{thr}"] = outgrowth_coverage(fu_mask, base_mask, persistence >= thr)
                    row[f"out_bim_thr_{thr}"] = outgrowth_coverage(fu_mask, base_mask, heat >= thr)
                out_rows.append(row)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out_rows


def main():
    print("=" * 78, flush=True)
    print("v142 TIME-STRATIFIED BIMODAL KERNEL COVERAGE on PROTEAS", flush=True)
    print(f"  sigma_broad = {SIGMA_BROAD}; N_BOOT = {N_BOOT}", flush=True)
    print("=" * 78, flush=True)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v142_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} patient zips", flush=True)
            t0 = time.time()
            for i, entry in enumerate(entries, 1):
                rows.extend(analyse_patient(outer, entry, work))
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal follow-ups: {len(rows)}", flush=True)
    if not rows: return

    # Stratify by fu_index
    strata = {
        "early_fu1": [r for r in rows if r["fu_index"] == 1],
        "mid_fu2": [r for r in rows if r["fu_index"] == 2],
        "late_fu3plus": [r for r in rows if r["fu_index"] >= 3],
    }
    out = {"version": "v142",
           "experiment": "Time-stratified bimodal kernel coverage on PROTEAS",
           "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_total_followups": len(rows),
           "n_bootstrap": N_BOOT,
           "strata": {}}

    for stratum_name, sub_rows in strata.items():
        n = len(sub_rows)
        n_pids = len(set(r["pid"] for r in sub_rows))
        print(f"\n--- {stratum_name}: N = {n} follow-ups across {n_pids} patients ---",
              flush=True)
        if n == 0:
            continue
        stratum_results = {"n_followups": n, "n_patients": n_pids, "thresholds": {}}
        for thr in HEAT_THRESHOLDS:
            print(f"  threshold heat >= {thr}:", flush=True)
            cov_pers = np.array([r[f"cov_pers_thr_{thr}"] for r in sub_rows], dtype=float)
            cov_bim = np.array([r[f"cov_bim_thr_{thr}"] for r in sub_rows], dtype=float)
            out_pers = np.array([r[f"out_pers_thr_{thr}"] for r in sub_rows], dtype=float)
            out_bim = np.array([r[f"out_bim_thr_{thr}"] for r in sub_rows], dtype=float)

            mp_o, lp_o, hp_o = vectorised_cluster_bootstrap(cov_pers)
            mb_o, lb_o, hb_o = vectorised_cluster_bootstrap(cov_bim)
            mp_x, lp_x, hp_x = vectorised_cluster_bootstrap(out_pers)
            mb_x, lb_x, hb_x = vectorised_cluster_bootstrap(out_bim)

            d_o = cov_bim - cov_pers
            d_x = out_bim - out_pers
            md_o, ld_o, hd_o = vectorised_cluster_bootstrap(d_o)
            md_x, ld_x, hd_x = vectorised_cluster_bootstrap(d_x)

            sig_o = "**SIG**" if (ld_o > 0 or hd_o < 0) else ""
            sig_x = "**SIG**" if (ld_x > 0 or hd_x < 0) else ""

            print(f"    persistence overall {mp_o*100:5.2f}%  bimodal overall {mb_o*100:5.2f}%  "
                  f"delta {md_o*100:+5.2f} pp [{ld_o*100:+.2f}, {hd_o*100:+.2f}] {sig_o}",
                  flush=True)
            print(f"    persistence outgrwh {mp_x*100:5.2f}%  bimodal outgrwh {mb_x*100:5.2f}%  "
                  f"delta {md_x*100:+5.2f} pp [{ld_x*100:+.2f}, {hd_x*100:+.2f}] {sig_x}",
                  flush=True)

            stratum_results["thresholds"][f"heat_ge_{thr}"] = {
                "persistence_overall_mean_pct": round(mp_o * 100, 2),
                "bimodal_overall_mean_pct": round(mb_o * 100, 2),
                "persistence_outgrowth_mean_pct": round(mp_x * 100, 2),
                "bimodal_outgrowth_mean_pct": round(mb_x * 100, 2),
                "delta_overall_pp": round(md_o * 100, 2),
                "delta_overall_ci95_pp": [round(ld_o * 100, 2), round(hd_o * 100, 2)],
                "delta_overall_excludes_zero": bool(ld_o > 0 or hd_o < 0),
                "delta_outgrowth_pp": round(md_x * 100, 2),
                "delta_outgrowth_ci95_pp": [round(ld_x * 100, 2), round(hd_x * 100, 2)],
                "delta_outgrowth_excludes_zero": bool(ld_x > 0 or hd_x < 0),
            }
        out["strata"][stratum_name] = stratum_results

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
