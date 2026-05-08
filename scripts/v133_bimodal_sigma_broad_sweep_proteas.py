"""v133: Sigma_broad sweep for the v130 bimodal kernel on PROTEAS.

The v130 kernel used heat = max(persistence, gaussian(mask, sigma=4))
with sigma_broad = 4.0 chosen by inspection of the v127 sigma_opt
distribution. v133 sweeps sigma_broad in {1, 2, 3, 4, 5, 6, 7} to
identify the data-driven optimum and provide cluster-bootstrap CIs.

Outputs:
    source_data/v133_bimodal_sigma_broad_sweep.json
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
OUT_JSON = ROOT / "source_data" / "v133_bimodal_sigma_broad_sweep.json"

SIGMA_BROAD_GRID = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
HEAT_THRESHOLDS = [0.5, 0.8]
N_BOOT = 10000
RNG = np.random.default_rng(13301)


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


def heat_constant(mask, sigma):
    if mask.sum() == 0: return np.zeros_like(mask, dtype=np.float32)
    h = gaussian_filter(mask.astype(np.float32), sigma=sigma)
    if h.max() > 0: h = h / h.max()
    return h.astype(np.float32)


def heat_bimodal(mask, sigma_broad):
    persistence = mask.astype(np.float32)
    h_broad = heat_constant(mask, sigma_broad)
    return np.maximum(persistence, h_broad)


def coverage(future_mask, region_mask):
    fm = future_mask.astype(bool)
    if fm.sum() == 0: return float("nan")
    return float((fm & region_mask.astype(bool)).sum() / fm.sum())


def outgrowth_coverage(future_mask, baseline_mask, region_mask):
    fut = future_mask.astype(bool); base = baseline_mask.astype(bool)
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
            if entries["baseline"] not in inner.namelist(): return []
            baseline = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = baseline > 0
            if base_mask.sum() == 0: return []
            heat_maps = {sb: heat_bimodal(base_mask, sb) for sb in SIGMA_BROAD_GRID}
            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any(): continue
                row = {"pid": pid, "fu_name": Path(fu_name).stem}
                for thr in HEAT_THRESHOLDS:
                    for sb in SIGMA_BROAD_GRID:
                        row[f"cov_sb{sb}_thr{thr}"] = coverage(fu_mask, heat_maps[sb] >= thr)
                        row[f"out_sb{sb}_thr{thr}"] = outgrowth_coverage(
                            fu_mask, base_mask, heat_maps[sb] >= thr)
                out_rows.append(row)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out_rows


def main():
    print("=" * 78)
    print("v133 BIMODAL SIGMA_BROAD SWEEP on PROTEAS")
    print(f"  sigma_broad grid: {SIGMA_BROAD_GRID}")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v133_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            clinical = load_clinical_table(outer, work)
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} nested patient zips", flush=True)
            t0 = time.time()
            for i, entry in enumerate(entries, 1):
                rows.extend(analyse_patient(outer, entry, clinical, work))
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal follow-ups: {len(rows)}", flush=True)
    if not rows: return

    pid_arr = np.array([r["pid"] for r in rows])
    out = {"version": "v133", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_followups": len(rows),
           "sigma_broad_grid": SIGMA_BROAD_GRID,
           "n_bootstrap_replicates": N_BOOT,
           "thresholds": {}}

    for thr in HEAT_THRESHOLDS:
        thr_results = {"overall": {}, "outgrowth": {}}
        print(f"\n--- heat >= {thr} ---", flush=True)
        for sb in SIGMA_BROAD_GRID:
            vals = np.array([r[f"cov_sb{sb}_thr{thr}"] for r in rows], dtype=float)
            mean, lo, hi = cluster_bootstrap_ci(vals, pid_arr)
            thr_results["overall"][f"sigma_broad_{sb}"] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            vals_out = np.array([r[f"out_sb{sb}_thr{thr}"] for r in rows], dtype=float)
            valid = ~np.isnan(vals_out)
            mean_o, lo_o, hi_o = cluster_bootstrap_ci(vals_out, pid_arr)
            thr_results["outgrowth"][f"sigma_broad_{sb}"] = {
                "mean_pct": round(mean_o * 100, 2),
                "ci95_pct": [round(lo_o * 100, 2), round(hi_o * 100, 2)],
                "n_with_outgrowth": int(valid.sum()),
            }
            print(f"  sigma_broad={sb}: overall {mean*100:5.2f}% "
                  f"[{lo*100:.2f}, {hi*100:.2f}]  | outgrowth {mean_o*100:5.2f}% "
                  f"[{lo_o*100:.2f}, {hi_o*100:.2f}]")

        # Find optimum on overall
        best_overall = max(thr_results["overall"].items(),
                            key=lambda kv: kv[1]["mean_pct"])
        best_outgrowth = max(thr_results["outgrowth"].items(),
                              key=lambda kv: kv[1]["mean_pct"])
        thr_results["best_sigma_broad_overall"] = best_overall[0]
        thr_results["best_sigma_broad_outgrowth"] = best_outgrowth[0]
        print(f"\n  -> best sigma_broad on overall: {best_overall[0]} "
              f"({best_overall[1]['mean_pct']:.2f}%)")
        print(f"  -> best sigma_broad on outgrowth: {best_outgrowth[0]} "
              f"({best_outgrowth[1]['mean_pct']:.2f}%)")

        out["thresholds"][f"heat_ge_{thr}"] = thr_results

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}")


if __name__ == "__main__":
    main()
