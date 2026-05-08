"""v117: Paired cluster-bootstrap test of the v98 anisotropic BED-aware
kernel vs the LESION-PERSISTENCE BASELINE on PROTEAS-brain-mets,
matching follow-ups by (pid, fu order index).

Uses the original v98_anisotropic_bed_per_patient.csv (the actual
v98 anisotropic coverage values) joined with the v116
per-patient CSV (which contains the persistence baseline coverage).

This tests whether the v98 +12.33 pp anisotropic gain over constant
sigma=2.5 survives against the most aggressive structural baseline:
the lesion mask itself ('persistence baseline').

Outputs:
    source_data/v117_aniso_vs_persistence_paired.json
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
V98_CSV = ROOT / "source_data" / "v98_anisotropic_bed_per_patient.csv"
V116_CSV = ROOT / "source_data" / "v116_anisotropic_vs_persistence_per_patient.csv"
OUT_JSON = ROOT / "source_data" / "v117_aniso_vs_persistence_paired.json"

N_BOOT = 10000
RNG = np.random.default_rng(11701)


def fnum(s, default=np.nan):
    s = (s or "").strip()
    if not s:
        return default
    try:
        return float(s)
    except ValueError:
        return default


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
    print("v117 PAIRED v98 ANISOTROPIC vs PERSISTENCE BASELINE on PROTEAS")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    v98 = list(csv.DictReader(open(V98_CSV)))
    v116 = list(csv.DictReader(open(V116_CSV)))

    print(f"v98 CSV: {len(v98)} rows ({len(v98)//2} follow-ups x 2 thresholds)")
    print(f"v116 CSV: {len(v116)} rows ({len(v116)} follow-ups)")

    # Index v98 by (pid, threshold) -> {const, iso, aniso}
    v98_by_pid_thr = {}
    for r in v98:
        key = (r["pid"], float(r["threshold"]))
        v98_by_pid_thr.setdefault(key, []).append({
            "const": fnum(r["coverage_const_sigma"]),
            "iso": fnum(r["coverage_iso_bed"]),
            "aniso": fnum(r["coverage_aniso_bed"]),
        })

    # Index v116 by (pid, fu_order_within_pid) -> {sigma_*, persistence}
    v116_by_pid = {}
    for r in v116:
        v116_by_pid.setdefault(r["pid"], []).append(r)

    # v98 stores coverage per patient per follow-up as separate rows
    # within v98_by_pid_thr[(pid, thr)]. The follow-ups are listed in
    # patient order. Same in v116. Match by index within patient.
    matched = []  # (pid, threshold, persistence, aniso, const, iso)
    for pid in v116_by_pid.keys():
        v116_rows = v116_by_pid[pid]
        for thr in [0.5, 0.8]:
            v98_rows = v98_by_pid_thr.get((pid, thr), [])
            n = min(len(v116_rows), len(v98_rows))
            for i in range(n):
                matched.append({
                    "pid": pid,
                    "threshold": thr,
                    "persistence": fnum(v116_rows[i][f"cov_persistence_thr_{thr}"]),
                    "sigma_1.0": fnum(v116_rows[i][f"cov_sigma_1.0_thr_{thr}"]),
                    "sigma_2.5": fnum(v116_rows[i][f"cov_sigma_2.5_thr_{thr}"]),
                    "aniso_bed": v98_rows[i]["aniso"],
                    "const_sigma": v98_rows[i]["const"],
                    "iso_bed": v98_rows[i]["iso"],
                })
    print(f"\nMatched {len(matched)} (pid, follow-up, threshold) tuples\n")

    out = {"version": "v117", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_matched": len(matched), "n_bootstrap_replicates": N_BOOT,
           "alpha": 0.05, "thresholds": {}}

    pid_arr_full = np.array([m["pid"] for m in matched])
    print(f"Unique patients in matched set: {len(np.unique(pid_arr_full))}\n")

    for thr in [0.5, 0.8]:
        sub = [m for m in matched if m["threshold"] == thr]
        sub_pids = np.array([s["pid"] for s in sub])
        print(f"--- heat >= {thr} (N = {len(sub)}; "
              f"{len(np.unique(sub_pids))} patients) ---")

        thr_results = {"methods": {}, "paired_deltas_vs_aniso": {},
                       "paired_deltas_vs_persistence": {}}

        # Per-method point estimates with CIs
        method_keys = ["persistence", "sigma_1.0", "sigma_2.5",
                       "aniso_bed", "const_sigma", "iso_bed"]
        for key in method_keys:
            vals = np.array([s[key] for s in sub], dtype=float)
            mean, lo, hi = cluster_bootstrap_ci(vals, sub_pids)
            thr_results["methods"][key] = {
                "mean_pct": round(mean * 100, 2),
                "ci95_pct": [round(lo * 100, 2), round(hi * 100, 2)],
            }
            print(f"  {key:18s}: {mean*100:5.2f}% [{lo*100:5.2f}, {hi*100:5.2f}]")
        print()

        # Paired delta vs anisotropic BED (the v98 result; "is it better than each baseline?")
        aniso_vals = np.array([s["aniso_bed"] for s in sub], dtype=float)
        for key in ["persistence", "sigma_1.0", "sigma_2.5",
                    "const_sigma", "iso_bed"]:
            base_vals = np.array([s[key] for s in sub], dtype=float)
            diff = aniso_vals - base_vals
            mean, lo, hi = cluster_bootstrap_ci(diff, sub_pids)
            thr_results["paired_deltas_vs_aniso"][f"aniso_minus_{key}"] = {
                "mean_pp": round(mean * 100, 2),
                "ci95_pp": [round(lo * 100, 2), round(hi * 100, 2)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            sig = "**SIG**" if (lo > 0 or hi < 0) else ""
            print(f"  delta aniso - {key:14s}: {mean*100:+6.2f} pp "
                  f"[{lo*100:+.2f}, {hi*100:+.2f}] {sig}")
        print()

        # Paired delta vs persistence (i.e. "does each method beat the persistence baseline?")
        pers_vals = np.array([s["persistence"] for s in sub], dtype=float)
        for key in ["sigma_1.0", "sigma_2.5", "aniso_bed", "iso_bed", "const_sigma"]:
            method_vals = np.array([s[key] for s in sub], dtype=float)
            diff = method_vals - pers_vals
            mean, lo, hi = cluster_bootstrap_ci(diff, sub_pids)
            thr_results["paired_deltas_vs_persistence"][f"{key}_minus_persistence"] = {
                "mean_pp": round(mean * 100, 2),
                "ci95_pp": [round(lo * 100, 2), round(hi * 100, 2)],
                "excludes_zero": bool(lo > 0 or hi < 0),
            }
            sig = "**SIG**" if (lo > 0 or hi < 0) else ""
            print(f"  delta {key:14s} - persistence: {mean*100:+6.2f} pp "
                  f"[{lo*100:+.2f}, {hi*100:+.2f}] {sig}")
        print()

        out["thresholds"][f"heat_ge_{thr}"] = thr_results

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
