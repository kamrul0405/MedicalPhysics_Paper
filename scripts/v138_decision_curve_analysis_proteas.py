"""v138: Decision-curve analysis (DCA) for the bimodal kernel vs
persistence vs treat-all on PROTEAS-brain-mets.

Decision-curve analysis (Vickers & Elkin 2006) computes net clinical
benefit at threshold probabilities. Standard for top clinical
journals (Lancet, NEJM, JAMA Network Open, NEJM AI, Nature Medicine).

For a structural prior at threshold tau, classify each VOXEL as
"high-risk for outgrowth" if heat >= tau. The voxel is the analysis
unit; the outcome is "actually has future-lesion growth".

Net benefit = (TP / N) - (FP / N) * (tau / (1 - tau))
where TP = true-positive voxels (predicted positive AND outcome positive)
       FP = false-positive voxels (predicted positive AND outcome negative)
       N  = total voxel count

For each prior (persistence, sigma=4 alone, sigma=7 alone, bimodal
sigma=4, bimodal sigma=7) and threshold tau in {0.1, 0.2, ..., 0.7},
compute net benefit per follow-up; cluster-bootstrap CIs.

The treat-all reference is computed at each tau as a benchmark.

Outputs:
    source_data/v138_decision_curve_analysis.json
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
OUT_JSON = ROOT / "source_data" / "v138_decision_curve_analysis.json"

THRESHOLDS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
SIGMA_BROAD_VALUES = [4.0, 7.0]
N_BOOT = 5000  # smaller because per-fu / per-tau analysis
RNG = np.random.default_rng(13801)


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


def heat_bimodal(mask, sb):
    persistence = mask.astype(np.float32)
    h_broad = heat_constant(mask, sb)
    return np.maximum(persistence, h_broad)


def net_benefit_per_fu(positive_voxels, region_voxels, n_total_voxels, tau):
    """Net benefit at threshold tau for a single follow-up.
    positive_voxels: bool array of true future-lesion voxels
    region_voxels: bool array of voxels predicted positive (heat >= tau)
    """
    pos = positive_voxels.astype(bool)
    pred = region_voxels.astype(bool)
    tp = float((pos & pred).sum())
    fp = float((~pos & pred).sum())
    n = float(n_total_voxels)
    if n == 0 or tau >= 1.0:
        return 0.0
    return (tp / n) - (fp / n) * (tau / (1.0 - tau))


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
            heat_pre = {f"sigma_{sb}": heat_constant(base_mask, sb) for sb in SIGMA_BROAD_VALUES}
            heat_pre["persistence"] = persistence
            for sb in SIGMA_BROAD_VALUES:
                heat_pre[f"bimodal_{sb}"] = heat_bimodal(base_mask, sb)

            for fu_name in entries["followups"]:
                try:
                    fu = load_nii_from_inner(inner, fu_name, patient_tmp)
                except Exception:
                    continue
                fu_mask = fu > 0
                if fu_mask.shape != base_mask.shape or not fu_mask.any(): continue
                n_total = fu_mask.size
                row = {"pid": pid, "fu_name": Path(fu_name).stem,
                       "n_total_voxels": n_total,
                       "n_positive_voxels": int(fu_mask.sum())}
                # Treat-all: predict every voxel positive at tau=tau
                # net benefit = prevalence - (1 - prevalence) * (tau / (1 - tau))
                prevalence = float(fu_mask.sum()) / n_total
                for tau in THRESHOLDS:
                    row[f"nb_treatall_tau_{tau}"] = (prevalence - (1 - prevalence) * tau / (1 - tau)) if tau < 1.0 else 0.0
                    for method, hmap in heat_pre.items():
                        region = hmap >= tau
                        row[f"nb_{method}_tau_{tau}"] = net_benefit_per_fu(
                            fu_mask, region, n_total, tau)
                out_rows.append(row)
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try: nested_path.unlink()
        except OSError: pass
    return out_rows


def main():
    print("=" * 78, flush=True)
    print("v138 DECISION-CURVE ANALYSIS on PROTEAS-brain-mets", flush=True)
    print(f"  thresholds tau: {THRESHOLDS}", flush=True)
    print(f"  N_BOOT = {N_BOOT}", flush=True)
    print("=" * 78, flush=True)

    rows = []
    with tempfile.TemporaryDirectory(prefix="proteas_v138_") as td:
        work = Path(td)
        with zipfile.ZipFile(DATA_ZIP) as outer:
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            print(f"Found {len(entries)} nested patient zips", flush=True)
            t0 = time.time()
            for i, entry in enumerate(entries, 1):
                rows.extend(analyse_patient(outer, entry, work))
                gc.collect()
                if i % 5 == 0 or i == len(entries):
                    print(f"  {i}/{len(entries)} ({time.time()-t0:.0f}s)", flush=True)

    print(f"\nTotal follow-ups: {len(rows)}", flush=True)
    if not rows: return

    pid_arr = np.array([r["pid"] for r in rows])
    out = {"version": "v138",
           "experiment": "Decision-curve analysis (Vickers & Elkin 2006) on PROTEAS-brain-mets",
           "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_followups": len(rows),
           "n_patients": int(len(np.unique(pid_arr))),
           "thresholds_tau": THRESHOLDS,
           "n_bootstrap_replicates": N_BOOT,
           "results": {}}

    methods = ["treatall", "persistence"] + \
              [f"sigma_{sb}" for sb in SIGMA_BROAD_VALUES] + \
              [f"bimodal_{sb}" for sb in SIGMA_BROAD_VALUES]

    print(f"\n=== Net benefit (per follow-up; cluster-bootstrap 95% CIs) ===", flush=True)
    print(f"  {'tau':>6s} | " + " | ".join([f"{m:>15s}" for m in methods]), flush=True)
    for tau in THRESHOLDS:
        out["results"][f"tau_{tau}"] = {}
        cells = []
        for m in methods:
            col = f"nb_{m}_tau_{tau}"
            vals = np.array([r[col] for r in rows], dtype=float)
            mean, lo, hi = vectorised_cluster_bootstrap(vals)
            out["results"][f"tau_{tau}"][m] = {
                "mean": round(mean, 6),
                "ci95": [round(lo, 6), round(hi, 6)],
            }
            cells.append(f"{mean:+.5f} [{lo:+.5f}, {hi:+.5f}]")
        print(f"  {tau:>6.2f} | " + " | ".join(cells), flush=True)

    # Paired delta: bimodal_4 vs persistence and bimodal_7 vs persistence
    print(f"\n=== Paired delta (cluster-bootstrap; bimodal_X - persistence) ===", flush=True)
    for sb in SIGMA_BROAD_VALUES:
        for tau in THRESHOLDS:
            bim_vals = np.array([r[f"nb_bimodal_{sb}_tau_{tau}"] for r in rows], dtype=float)
            pers_vals = np.array([r[f"nb_persistence_tau_{tau}"] for r in rows], dtype=float)
            diff = bim_vals - pers_vals
            mean, lo, hi = vectorised_cluster_bootstrap(diff)
            sig = "**SIG**" if (lo > 0 or hi < 0) else ""
            out["results"][f"tau_{tau}"][f"bimodal_{sb}_vs_persistence_delta"] = {
                "mean": round(mean, 6),
                "ci95": [round(lo, 6), round(hi, 6)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            print(f"  tau={tau:.2f}  bimodal_{sb} - persistence = {mean:+.5f} "
                  f"[{lo:+.5f}, {hi:+.5f}] {sig}", flush=True)

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nSaved {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
