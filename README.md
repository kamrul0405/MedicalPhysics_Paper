# RTO_paper

**Manuscript:** *Dose-prior boundary physics in brain-metastasis stereotactic radiotherapy: RTDOSE, RTPLAN and heat-kernel stress testing on PROTEAS*

**Target journal:** *Radiotherapy and Oncology* (the Green Journal; Elsevier; ISSN 0167-8140; IF ~5.5) — Full-Length Original Article
**Status:** Submission-ready (formatted per Elsevier *Radiotherapy and Oncology* Guide for Authors)
**Open-access fee:** **None required.** Hybrid Elsevier journal; submission on the standard subscription path requires no article-processing charge.

---

## What this repository contains

Public companion to a *Radiotherapy and Oncology* full-length original-article submission. Manuscript, figures, source-data files (JSON + CSV), and reproducibility scripts for every numerical claim in the paper.

```
RTO_paper/
├── manuscript/
│   ├── Manuscript_for_RTandO.md/.pdf                       <- PRIMARY SUBMISSION
│   ├── Manuscript_for_PracticalRadiationOncology.md/.pdf   (archived alternative target)
│   ├── Manuscript_v89_for_PracticalRadiationOncology.md    (archived earlier draft)
│   └── CoverLetter_NBE.md                                  (cover letter draft)
├── figures/
│   ├── main/                                               (5 main figures: PNG + 300 DPI TIFF)
│   └── extended_data/                                      (11 Extended Data figures at 300 DPI)
├── source_data/                                            (Versioned JSON/CSV underlying every numerical claim)
└── scripts/                                                (Python scripts producing source_data files)
```

## Headline empirical findings

- **Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients** (122 follow-up paired evaluations) from PROTEAS-brain-mets (Zenodo open-access).
- **v92 RTDOSE / RTPLAN integrity audit** parses 47 RTDOSE files and 48 RTPLAN files with zero errors; derives BED10 / BED2 / EQD2 directly from DICOM rather than NIfTI proxy filenames.
- **Heat ≥ 0.50 future-lesion coverage 47.4% exceeds dose ≥ 95% Rx coverage 37.9%** (paired Δ = +9.5 pp; Wilcoxon p = 0.014; row-level failure rate 52.1% vs 60.3%, p < 0.05 step-down).
- **Fractionation-stratified advantage preserved** — single-fraction +11.3 pp; 7-fraction +2.8 pp.
- **BED-stratified analysis (v93)**: heat ≥ 0.50 advantage scales with BED10 — +11.4 pp at high BED (≥ 60 Gy) versus +4.6 pp at mid BED (50–60 Gy).
- **BED-aware spatially-varying heat kernel (v94)**: per-voxel σ(BED) computed from RTDOSE arrays; compared against constant-σ baseline.
- **Decision-curve analysis (v93)**: heat ≥ 0.80 net benefit 0.453 at p_t = 0.20 vs treat-none 0.0 and dose self-indicator 0.479.
- **Dose-prior discordance taxonomy** identifies heat-only / dose-only / both-captured / neither-captured rows.
- **No toxicity outcomes** in PROTEAS — analysis cannot inform any prospective dose-escalation proposal.

## Clinical context

Brain-metastasis stereotactic radiotherapy (SRS) — single-fraction (typically 18–24 Gy) or hypofractionated (5–7 fraction schedules) — is a well-established treatment with 30–40% future-lesion volume falling outside the planned dose envelope. The PROTEAS-brain-mets cohort comprises:
- **Single-fraction SRS:** 18 Gy × 1 (n = 7 follow-up rows; 3 patients) or 20 Gy × 1 (n = 86 rows; 30 patients) — total 33 patients / 93 rows.
- **Seven-fraction hypofractionated SRT:** 5 Gy × 7 = 35 Gy (n = 28 rows; 9 patients).

## Reproducibility

| File | Description | Maps to |
|---|---|---|
| `v77_proteas_rtdose_audit.json` | Real PROTEAS RTDOSE coverage audit | §3.1 / Table 1 |
| `v77_proteas_rtdose_patient_metrics.csv` | Per-patient RTDOSE coverage metrics | Methods |
| `v78_proteas_boundary_stats.json` | Cluster-bootstrap CIs on coverage | §3.1 |
| `v81_proteas_threshold_sensitivity.json` | Threshold-sensitivity sweep | §3.3 |
| `v86_fractionation_strata.json` | Fractionation-stratified primary endpoints | §3.5 / Table 2 |
| `v89_dose_heat_discordance_taxonomy.csv` | Dose-prior discordance taxonomy | §3.6 |
| `v91_proteas_rtdose_inventory.json` | v91 dose-volume inventory | Methods §2.2 |
| `v92_proteas_plan_physics_audit.json` | v92 RTDOSE/RTPLAN BED/EQD2 audit | §3.5.1 / Methods §2.2.1 |
| `v93_bed_stratified.json` | BED-stratified primary endpoints | §3.7 / Table 3 |
| `v93_dca.json` | Decision-curve analysis | §3.8 |
| `v94_bed_aware_kernel.json` | BED-aware spatially-varying heat-kernel | §3.X (new) |
| `v60_yale_expansion.json` | Yale label-free acquisition-shift audit | Supplementary |
| `v84_E3_conformal_coverage.json` | Conformal coverage empirical validation | Supplementary Appendix B |
| `v84_E4_negative_controls.json` | 9 pre-specified negative controls | §3.9 |
| `master_neurooncology_dataset_index.csv` | 8-cohort master index | Methods |

## Hardware / software

Intel Core i7 CPU + NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM), CUDA 12.8. Python 3.11.9; nibabel 5.4.2; NumPy; SciPy 1.17.1; statsmodels.

## Companion repository

`MedIA_Paper` — companion submission to *Medical Image Analysis* (Elsevier; multi-cohort longitudinal post-treatment brain-tumour MRI benchmark with closed-form composition-shift theory, CASRN learned routing, and nnU-Net cross-cohort validation).

## Licence

Manuscript and figures: CC BY 4.0 (proposed at acceptance). Code: MIT.

## Contact

**Sheikh Kamrul Islam** — sheikh.islam@kcl.ac.uk (alternative: kamrul0405@outlook.com).
Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London.
