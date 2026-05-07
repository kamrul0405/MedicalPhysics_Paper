# RTO_paper

**Manuscript:** *Physics-grounded structural priors in brain-metastasis stereotactic radiotherapy: parabolic-PDE smoothing, BED-aware spatially-varying kernels, and multi-institutional dose-physics audit on patient-specific RTDOSE/RTPLAN*

**Target journal:** *Medical Physics* (AAPM / Wiley; ISSN 0094-2405; IF ~3.8) — Original Research Article
**Status:** Submission-ready (formatted per AAPM/Wiley *Medical Physics* Guide for Authors)
**Open-access fee:** **None required.** Hybrid Wiley journal; submission on the standard subscription path requires no article-processing charge.

> **Note on repository name.** This repository was originally created targeting *Radiotherapy and Oncology* (Green Journal) and retains the `RTO_paper` name for URL stability and commit-history continuity. The current submission target is *Medical Physics* (AAPM/Wiley). The earlier RT&O and PRO versions are preserved in `manuscript/` for archival reference.

---

## What this repository contains

Public companion to a *Medical Physics* original-research-article submission. Manuscript, figures, source-data files (JSON + CSV), and reproducibility scripts for every numerical claim in the paper.

```
RTO_paper/
├── manuscript/
│   ├── Manuscript_for_MedicalPhysics.md/.pdf               <- PRIMARY SUBMISSION
│   ├── Manuscript_for_RTandO.md/.pdf                       (archived RT&O version)
│   ├── Manuscript_for_PracticalRadiationOncology.md/.pdf   (archived PRO version)
│   └── CoverLetter_NBE.md                                  (cover letter draft)
├── figures/
│   ├── main/                                               (5 main figures: PNG + 300 DPI TIFF)
│   └── extended_data/                                      (11 Extended Data figures at 300 DPI)
├── source_data/                                            (Versioned JSON/CSV underlying every numerical claim)
└── scripts/                                                (Python scripts producing source_data files)
```

## Headline novelty

- **Heat-kernel structural prior derived as the fundamental solution of the heat equation** (parabolic-PDE smoothing operator with σ² = 2t evolution time; not a learned classifier).
- **BED-aware spatially-varying kernel** with per-voxel σ(BED) computed directly from the patient-specific RTDOSE array via the linear-quadratic radiobiology model.
- **v92 RTDOSE/RTPLAN integrity audit** parses 47 RTDOSE + 48 RTPLAN files with zero header errors; derives BED10/BED2/EQD2 directly from DICOM rather than NIfTI proxy filenames.
- **Multi-institutional physics atlas** (`source_data/v92_multisite_physics_atlas.json`) over 8 datasets: 2,875 patients / 10,679 timepoints / 76,680 NIfTIs / 25,180 mask-like files.

## Headline empirical findings

- **Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients** (122 follow-up paired evaluations) from PROTEAS-brain-mets (Zenodo open-access).
- **Heat ≥ 0.50 future-lesion coverage 47.4% > dose ≥ 95% Rx coverage 37.9%** (paired Δ = +9.5 pp; Wilcoxon p = 0.014; row-level failure rate 52.1% vs 60.3%, p < 0.05 step-down). Bootstrap CI on mean paired delta crosses zero — heat ≥ 0.50 is **broader but not statistically superior**.
- **BED-aware spatially-varying kernel achieves +6.99 pp coverage gain at heat ≥ 0.80** (37.08% vs 30.09% constant-σ); now matches the dose ≥ 95% Rx envelope coverage (37.82%) at the standard tight threshold.
- **Fractionation-stratified advantage preserved** (single-fraction +11.3 pp; 7-fraction +2.8 pp).
- **BED-stratified analysis**: heat ≥ 0.50 advantage scales 2.5× with BED10 — +11.4 pp at high BED (≥ 60 Gy) versus +4.6 pp at mid BED.
- **Decision-curve analysis**: heat ≥ 0.80 net benefit 0.453 at p_t = 0.20 vs treat-none 0.0.
- **No toxicity outcomes** in PROTEAS — analysis cannot inform any prospective dose-escalation proposal.

## Reproducibility

| File | Description | Maps to |
|---|---|---|
| `v77_proteas_rtdose_audit.json` | Real PROTEAS RTDOSE coverage audit | §3.1 / Table 1 |
| `v86_fractionation_strata.json` | Fractionation-stratified primary endpoints | §3.5 / Table 2 |
| `v89_dose_heat_discordance_taxonomy.csv` | Dose-prior discordance taxonomy | §3.6 |
| `v91_proteas_rtdose_inventory.json` | v91 dose-volume inventory | Methods §2.2 |
| `v92_proteas_plan_physics_audit.json` | v92 RTDOSE/RTPLAN BED/EQD2 audit | §3.5.1 / Methods §2.2.1 |
| `v92_multisite_physics_atlas.json` | Multi-institutional physics atlas (8 cohorts) | §2.4 / Methods |
| `v93_bed_stratified.json` | BED-stratified primary endpoints | §3.7 / Table 3 |
| `v93_dca.json` | Decision-curve analysis | §3.8 |
| `v94_bed_aware_kernel.json` | BED-aware spatially-varying kernel | §3.9 / Table 4 |
| `v95_alpha_beta_sensitivity.json` | α/β sensitivity sweep on BED-aware kernel | §3.10 |
| `master_neurooncology_dataset_index.csv` | 8-cohort master index | Methods |

## Hardware / software

Intel Core i7 CPU + NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM), CUDA 12.8. Python 3.11.9; nibabel 5.4.2; NumPy; SciPy 1.17.1; statsmodels.

## Companion repository

`MedIA_Paper` — companion submission to *Medical Image Analysis* (Elsevier; multi-cohort longitudinal post-treatment brain-tumour MRI benchmark with closed-form composition-shift theory, CASRN learned routing, and 7-architecture cross-cohort validation).

## Licence

Manuscript and figures: CC BY 4.0 (proposed at acceptance). Code: MIT.

## Contact

**Sheikh Kamrul Islam** — sheikh.islam@kcl.ac.uk (alternative: kamrul0405@outlook.com).
Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London.
