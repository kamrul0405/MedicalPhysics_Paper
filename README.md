# Conditional-Use Boundary Stress Testing of Longitudinal Brain-Tumour Risk Maps

**Target journal:** Nature Biomedical Engineering
**Manuscript version:** v8.2 (submission-ready, 2026-05-06)
**Status:** Submission-ready

This repository contains the source manuscript, figures, source-data files and reproducibility scripts for the Nature Biomedical Engineering submission.

## Repository structure

```
Nature_BME_paper/
├── manuscript/
│   ├── Final_Nature_BME.md       # Markdown source (v8.2; 168-word abstract; ~3,237 main words)
│   ├── Final_Nature_BME.pdf      # Submission PDF (43 KB)
│   └── CoverLetter_NBE.md        # Cover letter for Nature Biomedical Engineering editor
├── figures/
│   ├── main/                     # Main Figure 1 (.png + .tif at 300 DPI)
│   └── extended_data/            # 11 Extended Data figures (.png at 300 DPI)
├── source_data/                  # Versioned JSON/CSV underlying every numerical claim
└── scripts/                      # Python scripts that produced source_data files
```

## Headline results

- **Real patient-specific RTDOSE audit** on PROTEAS-brain-mets cohort (43 patients, 122 follow-up pairs) with patient-clustered bootstrap CIs.
- **Future-lesion coverage**: dose ≥95% Rx covers 37.8% of future lesion volume; heat ≥0.80 covers 30.1%; heat ≥0.50 covers 47.3%; heat ≥0.80 fails (<50% coverage) in 71.1% of follow-up rows.
- **Heat regions track plan geometry**: Spearman ρ=0.979 between heat ≥0.80 V95 and GTV V95.
- **Threshold-sensitivity sweep** (heat 0.5–0.9 × dose 0.5–1.0) confirms findings are not threshold artefacts.
- **V@80 review burden**: heat 0.388 cm³ vs static prior 0.710 cm³ (45% reduction; CIs non-overlapping).
- **Two-axis pre-deployment audit**: pi*=0.43 endpoint-composition screen + Yale-derived label-free acquisition-shift screen (AUROC 0.847).
- **Honest negatives** retained: survival C-index 0.565, NHS NPV −£1,479/year at 100 cases/week.

## Engineering contribution

The conditional-use-boundary deployment framework is the engineering deliverable. It defines, before any clinical activation, the **seven decision axes** any AI deployment must answer (calibration, anatomical alignment, dosimetric envelope, decision burden, conditional-use boundary, temporal degradation, fairness/survival/decision-support gates). The heat-kernel risk map is a deliberately simple, fully-auditable instantiation; the framework can be applied to any AI tool.

## Reproducibility

Every numerical claim in the manuscript maps to a versioned `source_data/*.json` or `*.csv` file.

### Key source-data files

| File | Description |
|---|---|
| `v77_proteas_rtdose_audit.json` | Real PROTEAS RTDOSE audit (Result 3; Table 2) |
| `v77_proteas_rtdose_patient_metrics.csv` | Per-patient RTDOSE coverage metrics |
| `v78_proteas_boundary_stats.json` | Cluster-bootstrap CIs on coverage statistics |
| `v78_nbe_boundary_source_data.csv` | Figure source data |
| `v81_proteas_threshold_sensitivity.json` | Threshold sensitivity sweep (Result 3) |
| `v60_yale_expansion.json` | Yale label-free acquisition-shift audit (Result 11) |
| `v69_bme_upgrade.json` | DECIDE-AI/TRIPOD-AI/ESTRO-AAPM compliance, RANO 2.0, fairness, temporal |
| `master_neurooncology_dataset_index.csv` | 8-cohort master index (Methods) |

## Hardware / software

Analyses run on Intel Core i7 CPU + NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM), CUDA 12.8. Python 3.11.9; nibabel 5.4.2; NumPy 2.4.4; SciPy 1.17.1; statsmodels for meta-regression.

## Companion manuscript

A companion submission is in preparation for *Nature Machine Intelligence* (`Nature_MI_paper` repository) establishing the theoretical mixture-weighted Brier framework that underlies the conditional-use boundary.

## Citation

[To be added at acceptance.]

## License

Manuscript and figures: CC BY 4.0 (proposed at acceptance).
Code: MIT.

## Contact

Authors blinded for review.
