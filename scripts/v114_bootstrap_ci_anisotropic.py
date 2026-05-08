"""v114: Cluster bootstrap CIs on the v98 anisotropic BED-aware kernel
coverage findings.

Uses the v98 per-patient CSV to compute 10,000-replicate cluster bootstrap
95% CIs on the future-lesion coverage at heat ≥ 0.50 and heat ≥ 0.80, plus
the paired delta vs the constant-σ baseline and vs the dose ≥ 95% Rx
envelope.

Quantifies the uncertainty around the v98 +12 pp breakthrough so it can be
reported with proper inferential intervals rather than just a point estimate.

Outputs:
    source_data/v114_anisotropic_bootstrap_ci.json
"""
import csv
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
SRC_CSV = ROOT / "source_data" / "v98_anisotropic_bed_per_patient.csv"
OUT_JSON = ROOT / "source_data" / "v114_anisotropic_bootstrap_ci.json"

N_BOOT = 10000
RNG = np.random.default_rng(11401)


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
    return float(np.nanmean(boots)), float(lo), float(hi)


def cluster_bootstrap_ci_paired(diffs, pids, alpha=0.05):
    pids_unique = np.unique(pids)
    n_cl = len(pids_unique)
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sample = RNG.choice(pids_unique, size=n_cl, replace=True)
        vals = []
        for s in sample:
            mask = pids == s
            vals.extend(diffs[mask].tolist())
        boots[b] = np.mean(vals) if vals else np.nan
    lo = np.nanpercentile(boots, 100 * alpha / 2)
    hi = np.nanpercentile(boots, 100 * (1 - alpha / 2))
    return float(np.nanmean(boots)), float(lo), float(hi)


def main():
    print("=" * 78)
    print("v114 BOOTSTRAP CIs ON v98 ANISOTROPIC BED-aware KERNEL COVERAGE")
    print(f"  N_BOOT = {N_BOOT}")
    print("=" * 78)

    rows = list(csv.DictReader(open(SRC_CSV)))
    print(f"Loaded {len(rows)} per-patient rows from v98 CSV")

    pid_arr = np.array([r["pid"] for r in rows])

    out = {"version": "v114", "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "n_bootstrap_replicates": N_BOOT,
           "alpha": 0.05,
           "thresholds": {}}

    for thr in [0.5, 0.8]:
        sub = [r for r in rows if fnum(r["threshold"]) == thr]
        if not sub:
            continue
        sub_pids = np.array([r["pid"] for r in sub])
        cov_const = np.array([fnum(r["coverage_const_sigma"]) for r in sub])
        cov_iso = np.array([fnum(r["coverage_iso_bed"]) for r in sub])
        cov_aniso = np.array([fnum(r["coverage_aniso_bed"]) for r in sub])

        # Coverage CIs
        const_mean, const_lo, const_hi = cluster_bootstrap_ci(cov_const, sub_pids)
        iso_mean, iso_lo, iso_hi = cluster_bootstrap_ci(cov_iso, sub_pids)
        aniso_mean, aniso_lo, aniso_hi = cluster_bootstrap_ci(cov_aniso, sub_pids)

        # Paired-delta CIs
        delta_aniso_iso = cov_aniso - cov_iso
        delta_aniso_const = cov_aniso - cov_const
        d_ai_mean, d_ai_lo, d_ai_hi = cluster_bootstrap_ci_paired(delta_aniso_iso, sub_pids)
        d_ac_mean, d_ac_lo, d_ac_hi = cluster_bootstrap_ci_paired(delta_aniso_const, sub_pids)

        out["thresholds"][f"heat_ge_{thr}"] = {
            "n_rows": len(sub),
            "n_patients": int(len(np.unique(sub_pids))),
            "constant_sigma": {
                "mean_pct": round(const_mean * 100, 2),
                "ci95_pct": [round(const_lo * 100, 2), round(const_hi * 100, 2)],
            },
            "isotropic_bed": {
                "mean_pct": round(iso_mean * 100, 2),
                "ci95_pct": [round(iso_lo * 100, 2), round(iso_hi * 100, 2)],
            },
            "anisotropic_bed": {
                "mean_pct": round(aniso_mean * 100, 2),
                "ci95_pct": [round(aniso_lo * 100, 2), round(aniso_hi * 100, 2)],
            },
            "delta_aniso_minus_iso": {
                "mean_pp": round(d_ai_mean * 100, 2),
                "ci95_pp": [round(d_ai_lo * 100, 2), round(d_ai_hi * 100, 2)],
            },
            "delta_aniso_minus_const": {
                "mean_pp": round(d_ac_mean * 100, 2),
                "ci95_pp": [round(d_ac_lo * 100, 2), round(d_ac_hi * 100, 2)],
            },
        }

        print(f"\n  heat >= {thr}:")
        print(f"    Constant sigma:     {const_mean*100:.2f}% [{const_lo*100:.2f}, {const_hi*100:.2f}]")
        print(f"    Isotropic BED:      {iso_mean*100:.2f}% [{iso_lo*100:.2f}, {iso_hi*100:.2f}]")
        print(f"    Anisotropic BED:    {aniso_mean*100:.2f}% [{aniso_lo*100:.2f}, {aniso_hi*100:.2f}]")
        print(f"    Delta aniso - iso:  {d_ai_mean*100:+.2f} pp [{d_ai_lo*100:+.2f}, {d_ai_hi*100:+.2f}]")
        print(f"    Delta aniso - const: {d_ac_mean*100:+.2f} pp [{d_ac_lo*100:+.2f}, {d_ac_hi*100:+.2f}]")

    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
