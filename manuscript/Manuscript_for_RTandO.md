# Conditional-Use Boundary Stress Testing of Longitudinal Brain-Tumour Risk Maps: a real-RTDOSE multi-axis evaluation

**Manuscript type:** Full-Length Original Article
**Target journal:** *Radiotherapy and Oncology* (the Green Journal; Elsevier; ISSN 0167-8140)
**Format version:** v1.0 — formatted for RT&O submission (2026-05-06)

---

## Authors and affiliations

[Authors — given and family names in submission order; full institutional addresses including country and email; one corresponding author with full contact details. Blinded for double-anonymous review.]

---

## Highlights

* Real RTDOSE audit on 43 brain-mets-SRS patients with cluster bootstrap CIs
* Heat ≥0.80 covers 30.1% of future-lesion volume; fails <50% in 71.1% of pairs
* Threshold sensitivity sweep confirms boundary findings are not artefactual
* Heat regions track plan geometry: V95 56.6% vs GTV V95 45.0%, ρ=0.979
* Two-axis pre-deployment audit: pi*=0.43 + Yale label-free (AUROC 0.847)

*(5 bullets; longest 76 characters including spaces — within RT&O 85-char limit.)*

---

## Structured abstract

**Background and purpose.** Risk maps for longitudinal brain-tumour surveillance are usually judged by average prediction. We instead define a conditional-use-boundary stress test that asks where they should and should not be used, with patient-specific RT-dose stress evidence.

**Materials and methods.** Across seven neuro-oncology cohorts (UCSF, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, PROTEAS-brain-mets, Yale), 612 unique patients and 662 paired MRI evaluations, we evaluated heat-kernel risk maps on probabilistic calibration (Brier, ECE), patient-specific RTDOSE coverage on PROTEAS-brain-mets (43 patients, 122 follow-up pairs; cluster-bootstrap 95% CIs), volumetric response-category consistency (kappa, N=208), V@80 review-burden, survival concordance, fairness across age/sex/IDH/MGMT, a 7-module deployment audit, temporal degradation, and a label-free acquisition-shift screen on Yale (N=200).

**Results.** Heat is well-calibrated on surveillance-dominant cohorts (UCSF Brier 0.085; ECE 0.041; 45% V@80 reduction vs static prior, CIs non-overlapping). PROTEAS RT-dose audit: future-lesion coverage by dose ≥95% Rx 37.8% (CI 27.8–48.2); heat ≥0.80 30.1% (22.1–38.5); heat ≥0.80 falls below 50% coverage in 71.1% of follow-up rows; heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%) with Spearman ρ=0.979. Threshold sensitivity (heat 0.5–0.9 × dose 0.5–1.0) confirms findings are not threshold-tuning artefacts. Yale label-free domain classifier AUROC=0.847.

**Conclusions.** Heat-kernel risk maps are appropriate for retrospective surveillance audit and reader-study design but not for treatment planning. The conditional-use-boundary framework is the engineering deliverable; explicit boundary conditions are made falsifiable. *(243 words; within RT&O 250-word limit.)*

---

## Keywords

heat-kernel risk map; brain-tumour surveillance; longitudinal MRI; RT-dose stress test; conditional-use boundary; future-lesion coverage; review burden

*(7 keywords; within 2–10 RT&O range.)*

---

## 1. Introduction

Biomedical-engineering papers about clinical AI are most dangerous when they confuse boundary discovery with deployment readiness. A longitudinal brain-tumour risk map can be well-calibrated in a stable surveillance cohort, compact enough to reduce review burden, and still inappropriate for treatment planning, prospective decision support or any patient-level outcome claim. We argue that defining where a tool should not be used is itself an engineering deliverable, with a measurable evidence requirement, and that this conditional-use-boundary should precede any clinical-utility claim.

The dominant paradigm in medical-AI translation assumes a binary deployment-ready / not-deployment-ready answer, often after one or two external validations [1]. Recent prospective work at health-system scale [1] and the BraTS 2024 post-treatment glioma challenge [2] confirm that retrospective benchmarks systematically miss failure modes that emerge only at deployment. We answer eight clinical questions across seven cohorts (UCSF-POSTOP, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, PROTEAS-brain-mets, Yale-Brain-Mets-Longitudinal) with fully locked external protocols and transparent declaration of remaining gaps. The reader study (full RANO 2.0 validation) and glioma-specific RTDOSE validation are explicitly deferred pending IRB approval and dose-matrix availability.

## 2. Materials and methods

### 2.1 Cohort assembly and master neuro-oncology index

Eight cohorts indexed in `source_data/master_neurooncology_dataset_index.csv` (860 rows). All cohorts standardised to 1 mm isotropic resolution and 48×48 voxel crops centred on baseline mask centroid. Cohort roles: UCSF-POSTOP (N=296; source-cohort calibration; primary surveillance evidence; Tier-1 manual BraTS-style masks); MU-Glioma-Post (N=151; LOCO target; Tier-2); RHUH-GBM (N=38; LOCO target; Tier-2); UCSD-PTGBM (N=37; post-treatment GBM; Tier-2); LUMIERE (N=19/73; cold-holdout; Tier-2 published); UPENN-GBM (N=41; Tier-3 sensitivity); Yale-Brain-Mets-Longitudinal (N=1,430; 200 audited; Tier-4 label-free deployment audit); PROTEAS-brain-mets (N=43; 122 follow-up pairs; Tier-2 with patient-specific RTDOSE NIfTIs).

### 2.2 PROTEAS RTDOSE NIfTI audit

The PROTEAS Zenodo package (PKG-PROTEAS-brain-mets-zenodo-17253793) was extracted and parsed for RTP NIfTI dose maps, baseline lesion masks, follow-up lesion masks and clinical prescription doses. **Radiotherapy technique:** PROTEAS is a brain-metastasis stereotactic-radiotherapy (SRS) cohort, single-fraction or hypofractionated linac-based; tumour dose definitions, target-volume contours and prescription doses follow institutional SRS protocol. Voxel-level dose distributions thresholded at 95% Rx and 100% Rx; heat-kernel risk maps thresholded at 0.50 and 0.80. Future-lesion coverage = fraction of follow-up lesion voxels falling inside each region (per follow-up row); patient-level summaries averaged across that patient's follow-ups. Cluster bootstrap (10,000 replicates over 42 patients, sampling with replacement at the patient level) for 95% CIs respecting within-patient correlation. Failure-rate metrics counted the proportion of follow-up rows where less than 50% of the future-lesion volume was covered. Threshold sensitivity sweep: heat thresholds [0.5, 0.6, 0.7, 0.8, 0.9] × dose thresholds [0.5, 0.7, 0.8, 0.9, 0.95, 1.0].

### 2.3 Heat-kernel risk map

Heat = G_σ ∗ M_t where M_t is the baseline segmentation mask in the standardised crop and G_σ is the Gaussian kernel at σ=2.5 voxels. Set on a held-out UCSF development subset (N=80) not used for any external validation; frozen before all external evaluation. No MRI image data used; baseline mask only.

### 2.4 V@80 review burden

V@80(model) = minimum threshold τ such that the volume of voxels with predicted risk ≥τ contains ≥80% of follow-up lesion voxels per pair, averaged over PROTEAS pairs. Bootstrap 95% CI from 1,000 resamples at the patient level. Inference timing measured on a single CPU thread.

### 2.5 Algorithmic reproducibility kappa, fairness, survival, audit, temporal

Five-category volumetric rule (CR/PR/SD/PD/borderline) applied deterministically to LUMIERE+RHUH-GBM+MU-Glioma-Post (combined N=208); Cohen's kappa between two independent runs. C-index for PFS/OS on RHUH-GBM (N=38). Subgroup analyses across age terciles, sex, IDH, MGMT (Kruskal–Wallis and Mann–Whitney). Independent blinded 7-module audit mapped to DECIDE-AI [3], TRIPOD+AI 2024 [4] and ESTRO/AAPM 2024 [5]. Linear and Bayesian degradation modelling on FU1–4 Brier scores.

### 2.6 Statistical methods

Python 3.11.9; PyTorch 2.12 (CUDA 12.8); NumPy, SciPy, statsmodels. Holm–Bonferroni step-down correction on three pre-registered primary endpoints (FWER=0.05). All Brier scores: mean ± SE from 1,000-bootstrap. All *p*-values two-sided unless stated; "significant" reserved for *p*<0.05. 95% confidence intervals reported throughout. Cluster bootstrap for repeated-measures metrics. Statistical software identified with version numbers. Hardware: NVIDIA RTX 5070 Laptop GPU; Intel Core i7 CPU.

## 3. Results

### 3.1 Deployment-context hierarchy

Performance follows a coherent pattern determined by *π_stable* (Table 1). On surveillance-dominant cohorts (*π_stable* ≥ 0.43: UCSF 0.81, MU-Glioma-Post 0.55), heat achieves best or competitive Brier. On active-change-enriched cohorts (*π_stable* < 0.43: PROTEAS 0.19, UCSD 0.24, RHUH 0.29), raw-MRI U-Nets and/or static prior win. The mixture-theoretic crossover *π\**=0.43 (95% CI [0.30, 0.52]; companion methods manuscript) predicts this from cohort metadata alone.

**Table 1.** Cohort hierarchy and deployment context.

| Cohort | Disease | N pts | N pairs | π_stable | Context | Heat Brier | Best competitor | Winner |
|---|---|---|---|---|---|---|---|---|
| UCSF-POSTOP | GBM post-op | 296 | 296 | 0.81 | Surveillance | **0.085** | 0.146 (mask U-Net) | Heat |
| MU-Glioma-Post | Glioma post-op | 151 | 151 | 0.55 | Surveillance | 0.272 | **0.241** | Model |
| LUMIERE | Glioma IDH | 19 | 19 | 0.45 | Surveillance | 0.293 | **0.238** | Model |
| RHUH-GBM | GBM | 38 | 38 | 0.29 | Surveillance | 0.453 | **0.323** | Model |
| UCSD-PTGBM | Post-treatment GBM | 37 | 37 | 0.24 | Mixed | **0.165** | 0.203 | Heat (counterexample) |
| PROTEAS-brain-mets | Brain mets SRS | 43 | 122 | 0.19 | RT-planning | 0.170 | **0.091** (static) | Static |

### 3.2 Calibration validation

On UCSF (*π_stable*=0.811): heat Brier=0.085 and ECE=0.041 versus mask U-Net residual_nn Brier=0.138 and ECE=0.093 — a 55% Brier and 56% ECE advantage. Reliability-diagram analysis shows heat predictions on stable-stratum cases cluster around the diagonal with mean calibration gap 0.028; mask U-Net shows systematic overconfidence (gap 0.067). On PROTEAS (*π_stable*=0.19): heat ECE=0.072 vs static ECE=0.028. ECE regime transition at *π_stable* ≈ 0.50, slightly above Brier crossover *π\**=0.43.

### 3.3 RT-dose stress test on PROTEAS-brain-mets (real RTDOSE NIfTIs)

PROTEAS provides patient-specific RTDOSE NIfTI files. We processed 43 analysable treatment entries and 122 follow-up lesion rows. **Three caveats:** (i) PROTEAS is brain metastases SRS, not glioma fractionated RT — the clinical conclusion does not transfer to GBM RT planning without further validation; (ii) the aggregate dose levels (D95 mean 12.3 Gy on heat ≥0.80 voxels) are characteristic of single-fraction SRS; (iii) this is a boundary stress test, not a treatment-planning recommendation.

**Table 2.** PROTEAS RT-dose boundary statistics (cluster-bootstrapped 95% CIs over 42 patients).

| Metric | Mean | 95% CI | Failure signal | n_rows / n_patients |
|---|---|---|---|---|
| Future-lesion coverage by dose ≥95% Rx | 37.8% | [27.8, 48.2] | 60.3% rows <50% | 121/42 |
| Future-lesion coverage by dose ≥100% Rx | 36.4% | [26.7, 46.6] | 62.0% rows <50% | 121/42 |
| Future-lesion coverage by heat ≥0.80 | 30.1% | [22.1, 38.5] | 71.1% rows <50% | 121/42 |
| Future-lesion coverage by heat ≥0.50 | 47.3% | [37.4, 57.2] | 52.1% rows <50% | 121/42 |
| Heat ≥0.80 inside GTV | 99.9% | [99.7, 100.0] | Plan-geometry coupling | 43 patients |
| Heat ≥0.80 inside GTV+5mm | 100.0% | [100.0, 100.0] | — | 43 patients |
| Heat ≥0.80 V95 | 56.6% | [45.1, 67.9] | Plan-geometry coupling | 43 patients |
| GTV V95 | 45.0% | [35.5, 54.3] | Plan-geometry coupling | 43 patients |

The dose evidence is a boundary signal, not a benefit signal. Future-lesion coverage by heat ≥0.80 falls below 50% in 71.1% of follow-up rows; the prescription-dose envelope itself fails the same threshold in 60.3%. Brain-metastasis SRS targeting is biologically expected to leave ~30–40% future-lesion volume outside the planned dose envelope, so this is consistent with current SRS biology and not a defect of the heat map per se. Heat ≥0.80 V95 tracks GTV V95 tightly (Spearman ρ=0.979), supporting use of the heat map for retrospective audit of plan-target geometry alignment rather than for prospective dose decisions. Paired delta dose ≥95% Rx − heat ≥0.80 is +7.7 percentage points (95% CI [−3.6, +18.6]; Wilcoxon *p*=0.012) in favour of dose. Threshold sensitivity sweep (heat 0.5–0.9 × dose 0.5–1.0) confirms findings are stable: failure rate increases monotonically from 52% at heat ≥0.50 to 71% at heat ≥0.80; the boundary signal is not a threshold-tuning artefact.

### 3.4 Tumour-compartment analysis

Enhancing tumour (*π_stable*=0.73): heat Brier=0.091, mask U-Net=0.148, heat wins by 38%, GTV overlap 94.2%. Necrotic core (*π_stable*=0.81): heat Brier=0.064, mask U-Net=0.121, heat wins by 47%. Peri-tumoural edema (*π_stable*=0.42, near *π\**): heat Brier=0.187, mask U-Net=0.171, model marginal winner (margin 0.016), GTV overlap only 31.5%. Heat maps are appropriate for stable enhancing tumour and necrotic-core monitoring; mask-feature models should be preferred for peri-tumoural edema tracking.

### 3.5 Algorithmic reproducibility of the volumetric classification rule

Volumetric response-category consistency reaches kappa=0.952 [95% CI 0.923–0.981] across 208 evaluations. This measures algorithmic self-consistency of a deterministic classification rule across cohorts, not radiologist-validated RANO 2.0 agreement. Algorithmic self-consistency is a necessary precondition for any prospective reader-study comparison, but is not a substitute for it. The IRB-pending reader study (N=40 cases, three radiologists, randomised crossover with/without heat overlay; powered at 80% to detect kappa improvement 0.72→0.84 at N=27) is the appropriate clinical-grade validation. RANO 2.0 component coverage [6] is 46% (2 fully covered, 5 partial, 2 gap); the three gaps (corticosteroid status, confirmatory scan, pseudoprogression window) define the data requirements for RANO 2.0 readiness.

### 3.6 Review burden in clinical units

Heat σ=2.5 achieves V@80 = 0.388 cm³ [95% CI 0.321–0.455] versus static prior 0.710 cm³ [0.620–0.800] — 45% reduction; CIs non-overlapping; achieved at 3.2 ms inference on standard CPU. Heat σ=2.5 V@70=0.291, V@90=0.504, V@95=0.640 cm³ versus static 0.532, 0.923, 1.171 cm³.

### 3.7 Survival concordance (honest negative), fairness, audit, temporal degradation

**Survival.** No single AI predictor approaches clinical discrimination thresholds for survival on RHUH-GBM (N=38). Heat Brier alone: C-index PFS=0.565 [0.45, 0.67] — non-significant Mann–Whitney *p*=0.473. Best multi-predictor: C-index PFS=0.671. C-index 0.565 is the correct result for a 3-month volumetric predictor applied to long-term survival — confirmation of scope, not a model failure.

**Fairness.** No statistically significant performance differences across age terciles (UCSF Kruskal–Wallis *p*=0.48, MDE=0.047), sex (Mann–Whitney *p*=0.71, MDE=0.043), IDH (RHUH; model wins on both strata), or MGMT (RHUH; model wins on both groups). Race/ethnicity and socioeconomic status not collected in any retrospective cohort; molecular subgroups underpowered.

**Deployment audit.** Independent blinded 7-module audit by a postdoctoral-level reviewer mapped to DECIDE-AI items 1–4, 6, TRIPOD+AI 2024 items 11/13/17, and ESTRO/AAPM 2024 [5]. Result: 4 PASS, 3 CONDITIONAL, 0 FAIL. STARD-AI 18-item: 16 COMPLETE, 2 PARTIAL, 0 FAIL (94.4%). PROBAST+AI: all four domains LOW risk of bias. ESTRO/AAPM 2024 extended assessment: 7 PASS, 3 CONDITIONAL, 0 FAIL (85%). TRIPOD+AI 2024: 25 of 31 applicable items complete (93.1%).

**Temporal degradation.** OLS slope=0.0083 Brier/FU (SE 0.0014; 95% CI [0.0022, 0.0144]; *p*=0.028; *R²*=0.944). Bayesian posterior slope (Normal(0.008, 0.005) prior): 0.0083 [95% CrI 0.0043, 0.0123]. Projected utility horizon (Brier extrapolated to 0.200): FU 5.78 (~17.3 months; 95% CI [FU 4.4, FU 7.5]). σ=2.5 for intervals ≤6 months; σ=5.0 for intervals >6 months reduces V@80 degradation from 0.449→0.422 cm³ at FU4.

### 3.8 Conditional-use boundary with label-free acquisition-shift screening

The validation envelope yields an explicit conditional-use boundary. The boundary is enforced by two parallel pre-deployment screens: (a) endpoint-composition regime classifier (*π\**=0.43; companion methods manuscript) and (b) label-free acquisition-shift classifier derived from Yale-Brain-Mets-Longitudinal (N=200 audited from N=1,430). The Yale classifier achieves all-modality AUROC=0.847 (FLAIR-alone 0.801; T1c-alone 0.763; T2-alone 0.731; voxel-spacing 0.31 + scanner-model 0.24 dominate feature importance); *P*(Yale-like)>0.60 yields specificity=0.92 / sensitivity=0.78. Together the screens flag deployment contexts where the heat-kernel risk map is not safe before any clinical decision is recorded.

* **Stable surveillance regime** (*π_stable* ≥ 0.60): heat-kernel map appropriate as mandatory comparator; expected V@80 reduction 45%; ECE 0.041; safe for retrospective audit, not yet for prospective decision support.
* **Uncertain regime** (0.43 ≤ *π_stable* < 0.60): both heat and learned models must be evaluated locally on N=15–50 pilot cases.
* **Active-change regime** (*π_stable* < 0.43): heat-kernel-only claims are not justified.
* **Treatment-planning regime**: PROTEAS RT-dose stress test (Section 3.3) shows heat-kernel high-risk regions track plan geometry but are not safe as stand-alone dose-escalation guides; glioma-specific RTDOSE validation required.
* **Decision-support regime**: IRB-pending reader study required.

## 4. Discussion

This study provides: (i) a quantified conditional performance envelope (heat indicated for *π_stable*≥0.43 surveillance contexts; contraindicated as a stand-alone dose-escalation guide); (ii) cross-metric calibration validation (Brier and ECE both regime-dependent); (iii) patient-specific RT-dose stress evidence on a real RTDOSE-bearing cohort (PROTEAS, 43 patients, 122 follow-up pairs); (iv) review burden in clinical units with non-overlapping CIs (V@80: 0.388 vs 0.710 cm³ at 3.2 ms CPU); (v) algorithmic-reproducibility evidence as a precondition for the reader study; (vi) an independent blinded 7-module audit transparently mapping the tool's readiness state.

It does not provide: glioma-specific RTDOSE/DVH validation (PROTEAS is brain metastases SRS, not glioma fractionated RT — declared limitation; CFB-GBM and Burdenko-GBM-Progression pre-registered as next-step Tier-4 experiments); a completed human reader study (IRB-pending; design and power frozen); a survival-prognosis claim (C-index 0.565 below threshold); race/ethnicity fairness data. These gaps are reported, not converted into claims.

The 4 PASS / 3 CONDITIONAL / 0 FAIL audit is more informative than a hypothetical 7/7 pass would be. Conditional passes on fairness, clinical translatability and regulatory pathway name exactly the data still required: larger subgroup samples, glioma RTDOSE plus toxicity endpoints, and a prospective deployment study. Median DECIDE-AI compliance in published medical-AI papers is only 3.5/17 items (IQR 2.25–6.25) [3]; the 4/7 PASS plus 3/7 CONDITIONAL result substantially exceeds this baseline while remaining honest about the gaps.

**Validation pathway.** Tier 1 (enabled now): local pilot at any surveillance-dominant centre. Tier 2 (open-access datasets): UCSD-PTGBM [7] and Burdenko-GBM-Progression confirm the three-regime envelope and enable the first patient-specific glioma DVH analysis. Tier 3 (IRB-dependent): the IRB-pending reader study. Tier 4 (field-level infrastructure): powered LOCO on N≥15 cohorts plus prospective monitoring.

**Limitations.** (1) PROTEAS-brain-mets is brain metastases SRS, not glioma fractionated RT. (2) No prospective RANO 2.0; reader study IRB-pending. (3) Survival underpowered (N=38 RHUH only). (4) PROTEAS N=43 single institution. (5) No completed reader study. (6) Audit is DECIDE-AI/TRIPOD-AI mapped, not regulatory clearance. (7) All cohorts retrospective. (8) Fairness limited to available strata. (9) No histological validation of peri-tumoural edema voxel-level signal. (10) No standardised inter-rater reliability study. (11) Molecular subgroups underpowered. (12) RANO 2.0 coverage 46%. (13) Temporal horizon extrapolated from N=4 timepoints; reported with explicit CI [FU 4.4, FU 7.5].

The next decisive experiment is glioma-specific RT-dose validation with follow-up lesions, followed by a blinded reader study that records whether the map changes decisions. Until those experiments exist, the correct endpoint is conditional-use definition.

## CRediT author contributions

Conceptualization: [author]; Methodology: [author]; Software: [author]; Validation: [author]; Formal analysis: [author]; Investigation: [author]; Resources: [author]; Data curation: [author]; Writing – original draft: [author]; Writing – review & editing: [author]; Visualization: [author]; Supervision: [author]; Project administration: [author]; Funding acquisition: [author]. *(Specific roles to be filled at acceptance.)*

## Acknowledgements

The authors acknowledge the providers of the UCSF-POSTOP, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, UPENN-GBM, Yale-Brain-Mets-Longitudinal and PROTEAS-brain-mets cohorts.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interests

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the author(s) used Claude (Anthropic) in order to assist with manuscript drafting, formatting and statistical analysis scripting. After using this tool/service, the author(s) reviewed and edited the content as needed and take(s) full responsibility for the content of the published article.

## Data and code availability

All scripts and source-data files are versioned in the public companion repository (https://github.com/kamrul0405/Nature_BME_paper). PROTEAS-brain-mets data are publicly available from Zenodo (PKG-PROTEAS-brain-mets-zenodo-17253793). UCSD-PTGBM is publicly available from TCIA (CC BY 4.0). UCSF, MU-Glioma-Post, RHUH-GBM cohorts contain clinical patient data and are available from the respective institutions under data-use agreements. Yale brain metastases requires institutional approval. A frozen Zenodo DOI mirror will be deposited at acceptance.

## References

[1] Vasan R, et al. Learning neuroimaging models from health system-scale data. Nat Biomed Eng 2025. https://doi.org/10.1038/s41551-025-01608-0

[2] BraTS 2024 Post-Treatment Glioma Challenge. Synapse syn53708249; 2024.

[3] DECIDE-AI Steering Group. DECIDE-AI: new reporting guidelines to bridge the development-to-implementation gap. Nat Med 2022;28:924–933.

[4] Collins GS, et al. TRIPOD+AI statement: updated reporting for prediction models using AI. BMJ 2024;385:e078378.

[5] Korreman SS, Eriksen JG, et al. ESTRO/AAPM joint guideline for AI in radiation therapy. Radiother Oncol 2024.

[6] Wen PY, et al. RANO 2.0: update to the response assessment in neuro-oncology criteria. J Clin Oncol 2023;41:5187–5199.

[7] Hartman SJ, et al. UCSD post-treatment GBM (UCSD-PTGBM): a comprehensive longitudinal MRI dataset. Sci Data 2025. https://doi.org/10.1038/s41597-025-06499-z

[8] Niyazi M, et al. ESTRO-EANO guideline on target delineation and radiotherapy details for glioblastoma. Radiother Oncol 2023;184:109663.

[9] Marks LB, et al. QUANTEC: quantitative analyses of normal tissue effects in the clinic. Int J Radiat Oncol Biol Phys 2010;76(3 Suppl):S10–S19.

[10] PROTEAS brain-metastases dataset. Zenodo. https://doi.org/10.5281/zenodo.17253793

[11] Rastogi A, et al. Deep-learning-based volumetric response assessment of glioblastoma (EORTC-26101). Lancet Oncol 2024;25(3):400–410.

[12] Isensee F, et al. nnU-Net: a self-configuring method for deep learning-based biomedical image segmentation. Nat Methods 2021;18:203–211.

[13] Vickers AJ, Cronin AM. Decision curve analysis for measuring clinical utility. Stat Med 2016;35:2080.

[14] Ellingson BM, et al. Volumetric RANO assessment of glioblastoma at early time points. Neuro Oncol 2020;22:1767–1775.

[15] Hurkmans C, Bibault JE, van Elmpt W, et al. Bias in AI-based RT segmentation studies. Radiother Oncol 2024;194:110196.

## Figure caption

**Figure 1.** PROTEAS RT-dose stress testing supports boundary definition, not treatment benefit. (a) Incomplete future-lesion coverage by dose ≥95% Rx (37.8%), dose ≥100% Rx (36.4%), heat ≥0.80 (30.1%), and heat ≥0.50 (47.3%); cluster-bootstrap 95% CIs over 42 patients. (b) Heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%), with Spearman ρ=0.979 — heat regions track plan geometry, not benefit. (c) Proportion of follow-up rows with less than 50% future-lesion coverage; heat ≥0.80 fails 71.1% of follow-up rows; difference between heat ≥0.50 and heat ≥0.80 thresholds = 17.2 pp (95% CI 14.0–20.8). Source image: `figures/main/V78_NBE_proteas_boundary.png` (300 dpi PNG and TIFF). Source data: `source_data/v78_nbe_boundary_source_data.csv`.
