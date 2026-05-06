"""v86: Fractionation-stratified primary-endpoint analysis on PROTEAS-brain-mets.

Stratifies the primary coverage endpoints by fractionation schedule:
  - Single-fraction (fx=1, Rx 18 or 20 Gy)
  - Multi-fraction hypofractionated (fx=7, Rx 35 Gy = 5 Gy x 7)

Computes per-stratum:
  - Mean future-lesion coverage by dose >=95% Rx, dose >=100% Rx, heat >=0.50,
    heat >=0.80 (with cluster-bootstrap 95% CI)
  - Failure rate (proportion of follow-up rows with <50% coverage)
  - BED per schedule under linear-quadratic with alpha/beta = 10 Gy (tumour)
    and alpha/beta = 2 Gy (late-responding brain)

Outputs:
  source_data/v86_fractionation_strata.json
"""
import csv
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
CSV_PATH = ROOT / "source_data" / "v77_proteas_rtdose_patient_metrics.csv"
OUT = ROOT / "source_data" / "v86_fractionation_strata.json"

N_BOOT = 10000
RNG = np.random.default_rng(8601)


def cluster_bootstrap_ci(values, pids, alpha=0.05):
    """95% CI via cluster bootstrap over patients."""
    pids_unique = np.unique(pids)
    n_cl = len(pids_unique)
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sample = RNG.choice(pids_unique, size=n_cl, replace=True)
        vals = []
        for s in sample:
            mask = pids == s
            vals.extend(values[mask].tolist())
        boots[b] = np.mean(vals) if vals else np.nan
    lo = np.nanpercentile(boots, 100 * alpha / 2)
    hi = np.nanpercentile(boots, 100 * (1 - alpha / 2))
    return float(np.mean(boots[~np.isnan(boots)])), float(lo), float(hi)


def bed(rx_gy, n_fx, ab):
    d = rx_gy / n_fx
    return n_fx * d * (1 + d / ab)


def main():
    rows = list(csv.DictReader(open(CSV_PATH)))

    def fnum(s, default=np.nan):
        s = (s or "").strip()
        if not s:
            return default
        try:
            return float(s)
        except ValueError:
            return default

    rows = [r for r in rows if (r.get("future_inside_dose95_pct") or "").strip() != ""]
    pid_arr = np.array([r["pid"] for r in rows])
    fx_arr = np.array([fnum(r["fractions"]) for r in rows])
    rx_arr = np.array([fnum(r["rx_gy"]) for r in rows])
    cov_dose95 = np.array([fnum(r["future_inside_dose95_pct"]) for r in rows])
    cov_dose100 = np.array([fnum(r["future_inside_dose100_pct"]) for r in rows])
    cov_heat80 = np.array([fnum(r["future_inside_heat80_pct"]) for r in rows])
    cov_heat50 = np.array([fnum(r["future_inside_heat50_pct"]) for r in rows])

    out = {"version": "v86", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "strata": {}}

    for stratum_name, mask_fn in [
        ("all", lambda: np.ones_like(fx_arr, dtype=bool)),
        ("single_fraction", lambda: fx_arr == 1.0),
        ("seven_fraction_hypofx", lambda: fx_arr == 7.0),
    ]:
        m = mask_fn()
        if m.sum() == 0:
            continue
        pids_str = pid_arr[m]
        n_pts = len(np.unique(pids_str))
        n_rows = m.sum()
        block = {
            "n_rows": int(n_rows),
            "n_patients": int(n_pts),
            "rx_distribution": {f"{rx:.0f}_gy": int((rx_arr[m] == rx).sum())
                                for rx in np.unique(rx_arr[m])},
        }

        for label, vals in [
            ("dose_ge_95pct_rx", cov_dose95[m]),
            ("dose_ge_100pct_rx", cov_dose100[m]),
            ("heat_ge_0p80", cov_heat80[m]),
            ("heat_ge_0p50", cov_heat50[m]),
        ]:
            mean, lo, hi = cluster_bootstrap_ci(vals, pids_str)
            failure_rate = float((vals < 50.0).mean())
            block[label] = {
                "mean_pct": round(mean, 2),
                "ci95": [round(lo, 2), round(hi, 2)],
                "failure_rate": round(failure_rate, 3),
            }

        out["strata"][stratum_name] = block

    # BED summary
    bed_summary = {}
    for fx, rx in [(1, 18), (1, 20), (7, 35)]:
        bed10 = bed(rx, fx, 10.0)
        bed2 = bed(rx, fx, 2.0)
        n_match = int(((fx_arr == fx) & (rx_arr == rx)).sum())
        bed_summary[f"{rx}gy_in_{fx}_fx"] = {
            "n_rows": n_match,
            "BED_alpha_beta_10gy_tumour": round(bed10, 1),
            "BED_alpha_beta_2gy_late_brain": round(bed2, 1),
        }
    out["BED_summary"] = bed_summary

    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT}")
    print(json.dumps(out["strata"], indent=2))


if __name__ == "__main__":
    main()
