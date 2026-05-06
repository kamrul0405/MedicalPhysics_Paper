# RTO_paper

**Manuscript:** *Spatial heat-kernel structural priors and prescription-dose envelopes for future-lesion coverage in brain-metastasis stereotactic radiotherapy: a patient-specific RTDOSE analysis on 43 patients*

**Target journal:** *Practical Radiation Oncology* (ASTRO / Elsevier; ISSN 1879-8500) — Original Research Article
**Status:** Submission-ready (formatted per Elsevier *Practical Radiation Oncology* Guide for Authors)
**Open-access fee:** **None required.** This is a hybrid Elsevier journal published on behalf of ASTRO; submission on the standard subscription path requires no article-processing charge.

> **Note on repository name.** This repository was originally created targeting *Radiotherapy and Oncology* (the Green Journal; Elsevier) and retains the `RTO_paper` name for URL stability and commit-history continuity. The current submission target is *Practical Radiation Oncology* (ASTRO / Elsevier). The Radiotherapy and Oncology version of the manuscript is preserved at `manuscript/Manuscript_for_RTandO.md` for archival reference.

---

## What this repository contains

Public companion to a *Practical Radiation Oncology* original-research-article submission. Manuscript, figures, source-data files (JSON + CSV), and reproducibility scripts for every numerical claim in the paper.

```
RTO_paper/
├── manuscript/
│   ├── Manuscript_for_PracticalRadiationOncology.md/.pdf  <- PRIMARY SUBMISSION
│   ├── Manuscript_for_RTandO.md/.pdf                      (archived RT&O version)
│   ├── Manuscript_v83_for_NBE.md/.pdf                     (archived NBE companion)
│   ├── Final_Nature_BME.md/.pdf                           (archived Nature BME version)
│   └── CoverLetter_NBE.md                                 (cover letter draft; rename to CoverLetter_PRO at submission)
├── figures/
│   ├── main/                                              (5 main figures: PNG + 300 DPI TIFF)
│   └── extended_data/                                     (11 Extended Data figures at 300 DPI)
├── source_data/                                           (Versioned JSON/CSV underlying every numerical claim)
└── scripts/                                               (Python scripts producing source_data files)
```

## Headline empirical findings

- **Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients** (122 follow-up paired evaluations) from PROTEAS-brain-mets (Zenodo open-access).
- **Heat ≥ 0.50 future-lesion coverage 47.3% exceeds dose ≥ 95% Rx coverage 37.8%** (paired Δ = +9.5 pp; Wilcoxon p = 0.014; row-level failure rate 52.1% vs 60.3%, p < 0.05 step-down).
- **Fractionation-stratified advantage preserved** — single-fraction +11.3 pp; 7-fraction +2.8 pp.
- **Hypofractionated 7-fraction cases achieve twofold-higher coverage** than single-fraction (dose ≥ 95% Rx 59.9% vs 31.5%; heat ≥ 0.50 62.8% vs 42.8%) — consistent with broader peri-lesional dose distribution and lower late-responding-brain BED.
- **Heat ≥ 0.80 V95 (56.6%) tracks GTV V95 (45.0%)** with Spearman ρ = 0.979 (an algebraic consistency check, not an independent finding).
- **Threshold-sensitivity sweep** (5 heat × 6 dose thresholds) confirms boundary stability.
- **No toxicity outcomes** are available in PROTEAS; the analysis cannot inform any prospective dose-escalation proposal.

## Clinical context

Brain-metastasis stereotactic radiotherapy (SRS) is a well-established treatment with single-fraction (typically 18–24 Gy) or hypofractionated (5–7 fraction schedules) delivery to GTV+0–2 mm margin. The 30–40% future-lesion volume falling outside the planned dose envelope is the well-characterised clinical reality of SRS recurrence patterns, driven by micro-metastatic seeding, dose conformity-vs-coverage trade-offs, and recurrence at the GTV margin and in distant satellite locations.

The PROTEAS-brain-mets cohort comprises:
- **Single-fraction SRS:** 18 Gy × 1 (n = 7 follow-up rows; 3 patients) or 20 Gy × 1 (n = 86 rows; 30 patients) — total 33 patients / 93 rows.
- **Seven-fraction hypofractionated SRT:** 5 Gy × 7 = 35 Gy (n = 28 rows; 9 patients).

## Reproducibility

| File | Description | Maps to |
|---|---|---|
| `v77_proteas_rtdose_audit.json` | Real PROTEAS RTDOSE coverage audit | §3.1 / Table 1 |
| `v77_proteas_rtdose_patient_metrics.csv` | Per-patient RTDOSE coverage metrics | Methods |
| `v78_proteas_boundary_stats.json` | Cluster-bootstrap CIs on coverage | §3.1 |
| `v78_nbe_boundary_source_data.csv` | Figure 1 source data | Fig 1 |
| `v81_proteas_threshold_sensitivity.json` | Threshold-sensitivity sweep | §3.3 |
| `v86_fractionation_strata.json` | Fractionation-stratified primary endpoints | §3.5 / Table 2 |
| `v60_yale_expansion.json` | Yale label-free acquisition-shift audit | Supplementary |
| `v84_E3_conformal_coverage.json` | Conformal coverage empirical validation | Supplementary Appendix B |
| `v84_E4_negative_controls.json` | 9 pre-specified negative controls | §3.6 |
| `master_neurooncology_dataset_index.csv` | 8-cohort master index | Methods |

## Hardware / software

Intel Core i7 CPU + NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM), CUDA 12.8. Python 3.11.9; nibabel 5.4.2; NumPy; SciPy 1.17.1; statsmodels.

## Companion repository

`MedIA_Paper` — companion submission to *Computers in Biology and Medicine* (Elsevier; multi-cohort longitudinal post-treatment brain-tumour MRI benchmark).

## Licence

Manuscript and figures: CC BY 4.0 (proposed at acceptance). Code: MIT.

## Contact

**Sheikh Kamrul Islam** — sheikh.islam@kcl.ac.uk (alternative: kamrul0405@outlook.com).
Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London.
