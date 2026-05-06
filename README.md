# RTO_paper

**Manuscript:** *Future-lesion coverage by AI heat-kernel risk maps relative to prescription-dose envelopes in brain-metastasis stereotactic radiotherapy: a 43-patient multi-axis recurrence-pattern analysis*

**Target journal:** Radiotherapy and Oncology (the Green Journal; Elsevier; ISSN 0167-8140)
**Submission version:** v85 (2026-05-06)
**Status:** Submission-ready (formatted per Elsevier RT&O Guide for Authors)

---

## What this repository contains

Public companion to a clinical-radiotherapy submission to *Radiotherapy and Oncology* (the ESTRO official journal). Manuscript, figures, source-data files (JSON + CSV), and reproducibility scripts for every numerical claim in the paper.

```
RTO_paper/
├── manuscript/
│   ├── Manuscript_v85_for_RTandO.md/.pdf     <- PRIMARY SUBMISSION
│   ├── Manuscript_for_RTandO.md/.pdf         (v8.2 formatted; superseded by v85)
│   ├── Manuscript_v83_for_NBE.md/.pdf        (OCUMS companion / NBE alt)
│   ├── Final_Nature_BME.md/.pdf              (v8.2 master)
│   └── CoverLetter_NBE.md                    (cover letter; rename to CoverLetter_RTandO at submission)
├── figures/
│   ├── main/                                 (Main Fig 1: PNG + 300 DPI TIFF)
│   └── extended_data/                        (11 Extended Data figures at 300 DPI)
├── source_data/                              (Versioned JSON/CSV underlying every numerical claim)
└── scripts/                                  (Python scripts producing source_data files)
```

## Headline empirical findings

- **Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients** (122 follow-up paired evaluations) from PROTEAS-brain-mets (Zenodo open-access).
- **Future-lesion coverage by dose ≥95% Rx:** 37.8% (cluster-bootstrap 95% CI 27.8–48.2). Fails <50% in 60.3% of follow-up rows.
- **Heat ≥0.50:** 47.3% (CI 37.4–57.2). **Heat ≥0.80:** 30.1% (CI 22.1–38.5). Fails <50% coverage in 71.1% of follow-up rows for heat ≥0.80.
- **Both regions cover 30–47% of future-lesion volume on average**, consistent with the published 30–40% biological out-of-envelope recurrence baseline for brain-metastasis SRS (Lehrer 2022; Bohoudi 2020; Wiggenraad 2011).
- **Heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%)** with Spearman ρ=0.979.
- **Threshold sensitivity sweep** (5 heat × 6 dose thresholds) confirms boundary stability.
- **Conformal three-regime classifier** achieves 1.00 empirical coverage at α=0.05 nominal across N=7 cohort calibration set.

## Clinical context

Brain-metastasis stereotactic radiotherapy (SRS) is a well-established treatment with single-fraction (18–24 Gy) or hypofractionated (5×6 Gy = 30 Gy) delivery to GTV+0–2 mm margin. The 30–40% future-lesion volume falling outside the planned dose envelope is the well-characterised clinical reality of SRS recurrence patterns, driven by micro-metastatic seeding, dose conformity-vs-coverage trade-offs, and recurrence at the GTV margin and in distant satellite locations.

## Reproducibility

| File | Description | Maps to |
|---|---|---|
| `v77_proteas_rtdose_audit.json` | Real PROTEAS RTDOSE coverage audit | §3.1 / Table 1 |
| `v77_proteas_rtdose_patient_metrics.csv` | Per-patient RTDOSE coverage metrics | Methods |
| `v78_proteas_boundary_stats.json` | Cluster-bootstrap CIs on coverage | §3.1 |
| `v78_nbe_boundary_source_data.csv` | Figure 1 source data | Fig 1 |
| `v81_proteas_threshold_sensitivity.json` | Threshold sensitivity sweep | §3.3 |
| `v60_yale_expansion.json` | Yale label-free acquisition-shift audit | §3.4 |
| `v84_E3_conformal_coverage.json` | Conformal coverage empirical validation | §3.4 |
| `master_neurooncology_dataset_index.csv` | 8-cohort master index | Methods §2.1 |

## Hardware / software

Intel Core i7 CPU + NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM), CUDA 12.8. Python 3.11.9; nibabel 5.4.2; NumPy; SciPy; statsmodels.

## Companion repository

`MedIA_Paper` — companion methodology submission to *Medical Image Analysis*.

## License

Manuscript and figures: CC BY 4.0 (proposed at acceptance). Code: MIT.

## Contact

[Authors blinded for review.]
