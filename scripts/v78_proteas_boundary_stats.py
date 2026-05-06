"""v78_proteas_boundary_stats.py

Cluster-bootstrap boundary statistics for the PROTEAS RT-dose audit.

This is a lightweight CPU follow-up to v77. It converts the empirical RT-dose
audit into reviewer-readable boundary statements for the NBE manuscript:
which regions cover future lesions, how often they fail, and whether heat-risk
regions should be interpreted as clinical benefit or as a conditional-use
stress test.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "05_results"
CSV = RESULTS / "v77_proteas_rtdose_patient_metrics.csv"
OUT = RESULTS / "v78_proteas_boundary_stats.json"
SEED = 7802
BOOT = 10000


def _cluster_sums(vals: pd.DataFrame, column: str) -> tuple[np.ndarray, np.ndarray]:
    grouped = vals.groupby("pid")[column].agg(["sum", "count"]).dropna()
    return grouped["sum"].to_numpy(dtype=float), grouped["count"].to_numpy(dtype=float)


def mean_ci_cluster(df: pd.DataFrame, column: str) -> dict:
    rng = np.random.default_rng(SEED)
    vals = df[["pid", column]].dropna()
    sums, counts = _cluster_sums(vals, column)
    if len(sums) == 0:
        return {"mean": None, "ci95_lo": None, "ci95_hi": None, "n_rows": 0, "n_patients": 0}
    pids = vals["pid"].unique()
    draws = []
    for _ in range(BOOT):
        sampled = rng.integers(0, len(sums), size=len(sums))
        draws.append(float(sums[sampled].sum() / counts[sampled].sum()))
    return {
        "mean": float(vals[column].mean()),
        "ci95_lo": float(np.percentile(draws, 2.5)),
        "ci95_hi": float(np.percentile(draws, 97.5)),
        "n_rows": int(len(vals)),
        "n_patients": int(len(pids)),
    }


def paired_delta_cluster(df: pd.DataFrame, a: str, b: str) -> dict:
    rng = np.random.default_rng(SEED)
    vals = df[["pid", a, b]].dropna().copy()
    vals["delta"] = vals[a] - vals[b]
    sums, counts = _cluster_sums(vals, "delta")
    if len(sums) == 0:
        return {"mean_delta": None, "ci95_lo": None, "ci95_hi": None, "p_delta_gt_0": None}
    pids = vals["pid"].unique()
    draws = []
    for _ in range(BOOT):
        sampled = rng.integers(0, len(sums), size=len(sums))
        draws.append(float(sums[sampled].sum() / counts[sampled].sum()))
    try:
        wilcoxon_p = float(stats.wilcoxon(vals[a], vals[b], zero_method="wilcox").pvalue)
    except Exception:
        wilcoxon_p = None
    return {
        "comparison": f"{a}_minus_{b}",
        "mean_delta": float(vals["delta"].mean()),
        "ci95_lo": float(np.percentile(draws, 2.5)),
        "ci95_hi": float(np.percentile(draws, 97.5)),
        "p_delta_gt_0": float(np.mean(np.array(draws) > 0)),
        "wilcoxon_p": wilcoxon_p,
        "n_rows": int(len(vals)),
        "n_patients": int(len(pids)),
    }


def threshold_failure(df: pd.DataFrame, column: str, thresholds: list[float]) -> dict:
    vals = df[["pid", column]].dropna()
    out = {}
    for threshold in thresholds:
        out[f"pct_rows_below_{threshold:g}"] = float(100 * (vals[column] < threshold).mean())
    return out


def main() -> None:
    start = time.time()
    if not CSV.exists():
        raise FileNotFoundError(f"Missing v77 PROTEAS audit CSV: {CSV}")
    df = pd.read_csv(CSV)
    follow = df[df["followup"].notna()].copy()

    columns = [
        "future_inside_dose95_pct",
        "future_inside_dose100_pct",
        "future_inside_heat80_pct",
        "future_inside_heat50_pct",
        "future_dmean_gy",
        "gtv_v95_pct",
        "heat80_v95_pct",
        "heat50_v95_pct",
    ]
    clustered_means = {c: mean_ci_cluster(follow if c.startswith("future_") else df.drop_duplicates("pid"), c) for c in columns}
    paired = {
        "dose95_vs_heat80": paired_delta_cluster(follow, "future_inside_dose95_pct", "future_inside_heat80_pct"),
        "dose95_vs_heat50": paired_delta_cluster(follow, "future_inside_dose95_pct", "future_inside_heat50_pct"),
        "heat50_vs_heat80": paired_delta_cluster(follow, "future_inside_heat50_pct", "future_inside_heat80_pct"),
        "dose100_vs_heat80": paired_delta_cluster(follow, "future_inside_dose100_pct", "future_inside_heat80_pct"),
    }
    failures = {
        c: threshold_failure(follow, c, [25, 50, 75])
        for c in [
            "future_inside_dose95_pct",
            "future_inside_dose100_pct",
            "future_inside_heat80_pct",
            "future_inside_heat50_pct",
        ]
    }

    patient_level = df.drop_duplicates("pid").copy()
    try:
        corr_heat_dose = stats.spearmanr(patient_level["heat80_v95_pct"], patient_level["gtv_v95_pct"], nan_policy="omit")
        heat_dose_corr = {
            "spearman_r": float(corr_heat_dose.statistic),
            "p": float(corr_heat_dose.pvalue),
            "n": int(patient_level[["heat80_v95_pct", "gtv_v95_pct"]].dropna().shape[0]),
        }
    except Exception:
        heat_dose_corr = {"spearman_r": None, "p": None, "n": 0}

    output = {
        "version": "v78_proteas_boundary_stats",
        "status": "PROTEAS_RT_DOSE_BOUNDARY_STATISTICS",
        "date": time.strftime("%Y-%m-%d"),
        "source_csv": str(CSV),
        "n_followup_rows": int(len(follow)),
        "n_patients": int(follow["pid"].nunique()),
        "cluster_bootstrap_replicates": BOOT,
        "clustered_means": clustered_means,
        "paired_deltas": paired,
        "coverage_failure_rates": failures,
        "heat80_v95_vs_gtv_v95_spearman": heat_dose_corr,
        "editorial_interpretation": (
            "PROTEAS provides real RT-dose maps, but the coverage statistics do not justify "
            "clinical-benefit or deployment-readiness claims. They support a conditional-use "
            "boundary: high-risk heat regions are anatomically tied to baseline GTV, while "
            "future-lesion coverage is incomplete and disease-context transfer is limited."
        ),
        "runtime_s": round(time.time() - start, 2),
    }
    OUT.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({
        "n_followup_rows": output["n_followup_rows"],
        "n_patients": output["n_patients"],
        "dose95_minus_heat80": paired["dose95_vs_heat80"],
        "heat50_minus_heat80": paired["heat50_vs_heat80"],
        "runtime_s": output["runtime_s"],
    }, indent=2))


if __name__ == "__main__":
    main()
