"""v77_proteas_rtdose_audit.py

Empirical PROTEAS RT-dose audit for the NBE paper.

The previous manuscript had a proxy DVH section. This script reads the local
PROTEAS archive, extracts each patient's NIfTI dose map (RTP), prescription
dose from the clinical workbook, and tumour masks, then computes dose-volume
coverage for heat-kernel high-risk regions and baseline/follow-up lesions.

Output:
  05_results/v77_proteas_rtdose_audit.json
  05_results/v77_proteas_rtdose_patient_metrics.csv
"""

from __future__ import annotations

import json
import argparse
import re
import shutil
import tempfile
import time
import zipfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation, gaussian_filter


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(r"C:\Users\kamru\Downloads\Datasets")
RESULTS = ROOT / "05_results"
PROTEAS_ZIP = DATA / "PKG - PROTEAS-brain-mets-zenodo-17253793.zip"
OUT_JSON = RESULTS / "v77_proteas_rtdose_audit.json"
OUT_CSV = RESULTS / "v77_proteas_rtdose_patient_metrics.csv"
SIGMA = 2.5
HEAT_THRESHOLD = 0.80
BOOT = 4000
SEED = 7702


def parse_rx(value) -> float | None:
    if pd.isna(value):
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", str(value))
    return float(match.group(1)) if match else None


def load_clinical_table(outer: zipfile.ZipFile, tmpdir: Path) -> pd.DataFrame:
    name = "PROTEAS-Clinical_and_demographic_data.xlsx"
    out = tmpdir / name
    out.write_bytes(outer.read(name))
    df = pd.read_excel(out, sheet_name="PROTEAS")
    df["pid"] = df["Patient ID (Zenodo)"].astype(str).str.strip()
    df["rx_gy"] = df["Rx dose at tumor margins"].map(parse_rx)
    df["fractions_num"] = pd.to_numeric(df["Fractions"], errors="coerce")
    return df


def load_nii_from_inner(inner: zipfile.ZipFile, name: str, tmpdir: Path) -> tuple[np.ndarray, tuple[float, float, float]]:
    out = tmpdir / Path(name).name
    out.write_bytes(inner.read(name))
    img = nib.load(str(out))
    arr = np.asanyarray(img.dataobj).astype(np.float32)
    zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
    return arr, zooms


def heat_map(mask: np.ndarray) -> np.ndarray:
    h = gaussian_filter(mask.astype(np.float32), sigma=SIGMA)
    mx = float(h.max())
    return (h / mx).astype(np.float32) if mx > 0 else h.astype(np.float32)


def dvh_metrics(dose: np.ndarray, structure: np.ndarray, rx: float) -> dict[str, float]:
    vals = dose[structure > 0].astype(np.float64)
    if len(vals) == 0:
        return {
            "n_vox": 0, "volume_cm3": 0.0, "d98_gy": float("nan"),
            "d95_gy": float("nan"), "dmean_gy": float("nan"),
            "v95_pct": float("nan"), "v100_pct": float("nan"),
            "v107_pct": float("nan"),
        }
    return {
        "n_vox": int(len(vals)),
        "volume_cm3": float(len(vals) * 0.001),  # PROTEAS NIfTI headers are 1 mm^3
        "d98_gy": float(np.percentile(vals, 2)),
        "d95_gy": float(np.percentile(vals, 5)),
        "dmean_gy": float(np.mean(vals)),
        "v95_pct": float(np.mean(vals >= 0.95 * rx) * 100),
        "v100_pct": float(np.mean(vals >= rx) * 100),
        "v107_pct": float(np.mean(vals >= 1.07 * rx) * 100),
    }


def bootstrap_ci(values: list[float]) -> dict[str, float]:
    vals = np.array([v for v in values if np.isfinite(v)], dtype=float)
    if len(vals) == 0:
        return {"mean": float("nan"), "ci95_lo": float("nan"), "ci95_hi": float("nan"), "n": 0}
    rng = np.random.default_rng(SEED)
    means = rng.choice(vals, size=(BOOT, len(vals)), replace=True).mean(axis=1)
    return {
        "mean": float(vals.mean()),
        "ci95_lo": float(np.percentile(means, 2.5)),
        "ci95_hi": float(np.percentile(means, 97.5)),
        "n": int(len(vals)),
    }


def find_patient_entries(inner: zipfile.ZipFile, pid: str) -> dict[str, object]:
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
        if any(n.startswith(f"{seg_dir}{pid}_tumor_mask_fu") for seg_dir in seg_dirs) and n.endswith(".nii.gz")
    ])
    return {"dose": dose, "baseline": baseline, "followups": followups}


def analyse_patient(outer: zipfile.ZipFile, entry: zipfile.ZipInfo, clinical: pd.DataFrame, work: Path) -> list[dict[str, object]]:
    pid = Path(entry.filename).stem
    row = clinical.loc[clinical["pid"] == pid]
    if row.empty:
        return []
    rx = parse_rx(row.iloc[0]["Rx dose at tumor margins"])
    if rx is None or not np.isfinite(rx):
        return []

    nested_path = work / entry.filename
    with outer.open(entry) as src, open(nested_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)

    patient_tmp = work / f"{pid}_files"
    patient_tmp.mkdir(exist_ok=True)
    metrics = []
    try:
        with zipfile.ZipFile(nested_path) as inner:
            entries = find_patient_entries(inner, pid)
            if entries["dose"] not in inner.namelist() or entries["baseline"] not in inner.namelist():
                return []
            dose, zooms = load_nii_from_inner(inner, entries["dose"], patient_tmp)
            base, _ = load_nii_from_inner(inner, entries["baseline"], patient_tmp)
            base_mask = base > 0
            if dose.shape != base_mask.shape:
                return []
            heat = heat_map(base_mask)
            high = heat >= HEAT_THRESHOLD
            high_05 = heat >= 0.50
            gtv_dvh = dvh_metrics(dose, base_mask, rx)
            high_dvh = dvh_metrics(dose, high, rx)
            high05_dvh = dvh_metrics(dose, high_05, rx)
            margin5 = binary_dilation(base_mask, iterations=5)
            common = {
                "pid": pid,
                "rx_gy": float(rx),
                "fractions": float(row.iloc[0]["fractions_num"]) if pd.notna(row.iloc[0]["fractions_num"]) else float("nan"),
                "voxel_mm3": float(np.prod(zooms)),
                "dose_max_gy": float(np.max(dose)),
                "baseline_gtv_cm3": float(base_mask.sum() * np.prod(zooms) / 1000.0),
                "heat80_cm3": float(high.sum() * np.prod(zooms) / 1000.0),
                "heat50_cm3": float(high_05.sum() * np.prod(zooms) / 1000.0),
                "heat80_inside_gtv_pct": float(np.mean(base_mask[high]) * 100) if high.any() else float("nan"),
                "heat80_inside_gtv_plus5mm_pct": float(np.mean(margin5[high]) * 100) if high.any() else float("nan"),
                "gtv_d95_gy": gtv_dvh["d95_gy"],
                "gtv_dmean_gy": gtv_dvh["dmean_gy"],
                "gtv_v95_pct": gtv_dvh["v95_pct"],
                "gtv_v100_pct": gtv_dvh["v100_pct"],
                "heat80_d95_gy": high_dvh["d95_gy"],
                "heat80_dmean_gy": high_dvh["dmean_gy"],
                "heat80_v95_pct": high_dvh["v95_pct"],
                "heat80_v100_pct": high_dvh["v100_pct"],
                "heat50_d95_gy": high05_dvh["d95_gy"],
                "heat50_dmean_gy": high05_dvh["dmean_gy"],
                "heat50_v95_pct": high05_dvh["v95_pct"],
                "heat50_v100_pct": high05_dvh["v100_pct"],
            }
            for fu_name in entries["followups"]:
                fu, _ = load_nii_from_inner(inner, fu_name, patient_tmp)
                fu_mask = fu > 0
                if not fu_mask.any():
                    continue
                fu_id = re.search(r"_mask_(fu\d+)", fu_name)
                metrics.append({
                    **common,
                    "followup": fu_id.group(1) if fu_id else Path(fu_name).stem,
                    "future_lesion_cm3": float(fu_mask.sum() * np.prod(zooms) / 1000.0),
                    "future_inside_dose95_pct": float(np.mean(dose[fu_mask] >= 0.95 * rx) * 100),
                    "future_inside_dose100_pct": float(np.mean(dose[fu_mask] >= rx) * 100),
                    "future_inside_heat80_pct": float(np.mean(high[fu_mask]) * 100),
                    "future_inside_heat50_pct": float(np.mean(high_05[fu_mask]) * 100),
                    "future_dmean_gy": float(np.mean(dose[fu_mask])),
                })
            if not metrics:
                metrics.append({**common, "followup": "none"})
    finally:
        shutil.rmtree(patient_tmp, ignore_errors=True)
        try:
            nested_path.unlink()
        except OSError:
            pass
    return metrics


def summarise(rows: list[dict[str, object]]) -> dict[str, object]:
    patient_rows = {}
    for row in rows:
        patient_rows.setdefault(row["pid"], row)
    patient_level = list(patient_rows.values())
    keys_patient = [
        "baseline_gtv_cm3", "heat80_cm3", "heat80_inside_gtv_pct",
        "heat80_inside_gtv_plus5mm_pct", "gtv_v95_pct", "gtv_v100_pct",
        "heat80_d95_gy", "heat80_v95_pct", "heat80_v100_pct",
        "heat50_v95_pct", "heat50_v100_pct",
    ]
    keys_pair = [
        "future_inside_dose95_pct", "future_inside_dose100_pct",
        "future_inside_heat80_pct", "future_inside_heat50_pct",
        "future_dmean_gy", "future_lesion_cm3",
    ]
    return {
        "patient_level": {k: bootstrap_ci([float(r[k]) for r in patient_level if k in r]) for k in keys_patient},
        "followup_pair_level": {k: bootstrap_ci([float(r[k]) for r in rows if k in r]) for k in keys_pair},
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PROTEAS RT-dose NIfTI audit")
    parser.add_argument("--limit", type=int, default=0, help="Optional first-N nested patient zips for smoke tests")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.time()
    if not PROTEAS_ZIP.exists():
        raise FileNotFoundError(PROTEAS_ZIP)
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="proteas_v77_") as td:
        work = Path(td)
        with zipfile.ZipFile(PROTEAS_ZIP) as outer:
            clinical = load_clinical_table(outer, work)
            entries = sorted(
                [e for e in outer.infolist() if re.fullmatch(r"P\d+[ab]?\.zip", e.filename)],
                key=lambda e: e.filename,
            )
            if args.limit and args.limit > 0:
                entries = entries[:args.limit]
            for i, entry in enumerate(entries, 1):
                patient_rows = analyse_patient(outer, entry, clinical, work)
                rows.extend(patient_rows)
                print(f"Analysed {i}/{len(entries)} {entry.filename}: {len(patient_rows)} follow-up rows")

    if not rows:
        raise RuntimeError("No PROTEAS rows analysed")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    summary = summarise(rows)
    out = {
        "version": "v77_proteas_rtdose_audit",
        "status": "EMPIRICAL_RTDose_NIFTI_AUDIT",
        "date": time.strftime("%Y-%m-%d"),
        "dataset": str(PROTEAS_ZIP),
        "n_patients": int(df["pid"].nunique()),
        "n_followup_pairs": int(len(df.dropna(subset=["followup"]))),
        "heat_threshold_primary": HEAT_THRESHOLD,
        "bootstrap_replicates": BOOT,
        "summary": summary,
        "interpretation": (
            "These are retrospective dose-map compatibility metrics from PROTEAS RT plan NIfTI maps. "
            "They replace prior proxy-DVH language but do not constitute prospective clinical utility, "
            "reader-study evidence, or glioma-specific radiotherapy validation."
        ),
        "runtime_s": round(time.time() - start, 2),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "n_patients": out["n_patients"],
        "n_followup_pairs": out["n_followup_pairs"],
        "heat80_v95_pct_mean": out["summary"]["patient_level"]["heat80_v95_pct"]["mean"],
        "future_inside_dose95_pct_mean": out["summary"]["followup_pair_level"]["future_inside_dose95_pct"]["mean"],
        "future_inside_heat80_pct_mean": out["summary"]["followup_pair_level"]["future_inside_heat80_pct"]["mean"],
        "runtime_s": out["runtime_s"],
    }, indent=2))


if __name__ == "__main__":
    main()
