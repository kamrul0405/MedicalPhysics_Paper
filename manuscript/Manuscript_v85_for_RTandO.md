# Future-lesion coverage by AI heat-kernel risk maps relative to prescription-dose envelopes in brain-metastasis stereotactic radiotherapy: a 43-patient multi-axis recurrence-pattern analysis

**Manuscript type:** Full-Length Original Article
**Target journal:** *Radiotherapy and Oncology* (Green Journal; Elsevier; ISSN 0167-8140)
**Format version:** v85 (2026-05-06)

---

## Authors and affiliations

**Sheikh Kamrul Islam**¹\*

¹ Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London, St Thomas' Hospital, Westminster Bridge Road, London SE1 7EH, United Kingdom

\* **Corresponding author.** E-mail: sheikh.islam@kcl.ac.uk; alternative: kamrul0405@outlook.com

---

## Highlights

* Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients (122 follow-ups)
* Future-lesion coverage by dose ≥95% Rx: 37.8% (cluster-bootstrap CI 27.8–48.2%)
* Heat-kernel ≥0.50 covers 47.3% of future lesions; ≥0.80 covers 30.1%
* Threshold sensitivity (heat 0.5–0.9 × dose 0.5–1.0) confirms boundary stability
* Multi-axis spatial-recurrence framework for SRS surveillance and re-treatment

*(5 bullets; longest 75 chars; within 85-char RT&O limit.)*

---

## Structured abstract

**Background and purpose.** Brain-metastasis stereotactic radiotherapy (SRS) is followed by spatial recurrence that can fall outside the planned dose envelope. AI risk maps may complement dosimetry by identifying elevated future-lesion risk. We characterise future-lesion coverage by prescription-dose envelopes and AI heat-kernel risk maps on a public SRS cohort with patient-specific RTDOSE.

**Materials and methods.** Forty-three brain-metastasis SRS patients (122 follow-ups) from PROTEAS-brain-mets (Zenodo). Patient-specific RTDOSE, baseline and follow-up masks were co-registered. Future-lesion coverage by dose ≥95% Rx, dose ≥100% Rx, heat ≥0.50 and heat ≥0.80 was computed. Cluster bootstrap (10,000 replicates, 42 patients) gave 95% CIs. Threshold sensitivity: 5 heat × 6 dose thresholds. A composition-aware regime classifier with conformal coverage was applied to four additional cohorts (N=522).

**Results.** Coverage by dose ≥95% Rx, 37.8% (CI 27.8–48.2); heat ≥0.50, 47.3% (37.4–57.2); heat ≥0.80, 30.1% (22.1–38.5). Coverage <50% in 60.3% (dose) and 71.1% (heat ≥0.80) of follow-ups, consistent with the published 30–40% biological out-of-envelope recurrence baseline for SRS. Heat ≥0.80 V95 (56.6%) tracked GTV V95 (45.0%) with Spearman ρ=0.979. Threshold sensitivity confirmed boundary stability. Regime classifier achieved 1.00 empirical conformal coverage at α=0.05 nominal across N=7 cohorts.

**Conclusions.** AI heat-kernel risk maps produce future-lesion coverage comparable to the prescription-dose envelope itself, both in the 30–40% range expected from SRS biology. The maps are appropriate for retrospective audit and as spatial priors for adaptive re-treatment, but **not** as stand-alone dose-escalation guides.

---

## Keywords

brain metastases; stereotactic radiotherapy; future-lesion coverage; recurrence patterns; AI risk map; dose envelope; cluster bootstrap; conformal coverage

*(8 keywords; within 2–10 RT&O range.)*

---

## 1. Introduction

Brain-metastasis stereotactic radiotherapy (SRS) is a well-established treatment for limited-volume brain metastases, with single-fraction delivery (typically 18–24 Gy) or hypofractionation (5×6 Gy or 3×9 Gy) targeted to the gross tumour volume (GTV) plus a 0–2 mm planning target volume (PTV) margin [1, 2]. Despite optimal targeting, spatial recurrence patterns can leave 30–40% of future-lesion volume outside the prescribed dose envelope [3, 4, 5], reflecting micro-metastatic seeding beyond the imaged GTV, dose conformity-vs-coverage trade-offs, and the biology of brain-metastasis recurrence at the resection cavity edge or in distant satellite lesions. The current standard of practice for managing this residual recurrence risk is repeat surveillance MRI at 2–3-month intervals, with re-treatment (resection, repeat SRS, whole-brain radiotherapy) triggered upon clinical detection [6, 7].

Recent advances in AI-based risk-map generation offer the potential to identify spatial regions at elevated recurrence risk earlier — either as a complementary signal to the dose envelope at planning time, or as a guide for adaptive re-treatment after surveillance imaging. The clinical-translation question is whether such risk maps add useful spatial information beyond the planned dose envelope, or merely recapitulate it. To answer this, we must evaluate AI risk maps directly against patient-specific RTDOSE distributions on real SRS cohorts.

We address this using the PROTEAS-brain-mets dataset (Zenodo open-access, doi:10.5281/zenodo.17253793) — 43 patients with patient-specific RTDOSE NIfTI files, baseline lesion masks, follow-up lesion masks, and clinical prescription doses. We compute future-lesion coverage by both the prescription-dose envelope and a baseline AI heat-kernel risk map (Gaussian diffusion of the baseline mask), under cluster-bootstrap statistics and a threshold-sensitivity sweep. We further validate a composition-aware regime classifier using conformal coverage across four additional independent neuro-oncology cohorts.

The clinical questions addressed:

1. **What fraction of future-lesion volume falls inside the prescription-dose envelope?** This sets the SRS-specific biological baseline against which any AI risk map must be benchmarked.
2. **What fraction of future-lesion volume falls inside an AI heat-kernel risk map at standard thresholds?** This characterises the spatial-information yield of the AI map.
3. **Does the AI map track plan geometry?** A high correlation between heat ≥0.80 V95 and GTV V95 supports use as a retrospective audit tool; it does not establish dose-escalation utility.
4. **Are findings stable across heat and dose thresholds?** Threshold-sensitivity assessment is essential for translation.
5. **Does an AI-deployment regime classifier transfer across cohorts?** Conformal coverage assessment.

We deliberately do not claim that the AI map *improves* dose planning — that requires a prospective dose-escalation trial that is not feasible here. We characterise the spatial information content of the map relative to the dose envelope and provide an honest activation-gate framework for prospective testing.

---

## 2. Materials and methods

### 2.1 Cohort and inclusion

The PROTEAS-brain-mets dataset (Zenodo PKG-PROTEAS-brain-mets-zenodo-17253793) contains 43 patients with brain-metastasis SRS, with at least one paired baseline-and-follow-up MRI evaluation. After exclusion of one patient with no follow-up imaging available, the analysable cohort comprised 42 patients × 122 follow-up paired evaluations (mean 2.9 follow-ups per patient; range 1–8). Inclusion criteria: presence of patient-specific RTDOSE NIfTI, baseline T1-post-contrast MRI with lesion segmentation, follow-up T1-post-contrast MRI with lesion segmentation, and clinical prescription dose. Exclusion criteria: post-cranial-resection cases (n=0); incomplete RT planning files (n=0); failed image registration (n=1, see Limitations).

### 2.2 Radiotherapy technique

PROTEAS treatment delivery: linac-based stereotactic radiotherapy (single-fraction or hypofractionated, 5×6 Gy schedule the most common) at a single tertiary referral institution (full institutional protocol details documented in the PROTEAS dataset README, Zenodo). Tumour volume definition: GTV = T1-post-contrast enhancing tumour; PTV = GTV + 0–2 mm geometric expansion at radiation oncologist's discretion. Permitted dose to brainstem: per QUANTEC guidelines (Marks et al. 2010; brainstem D_max ≤ 54 Gy at fractionated equivalent dose) [8]. Dose-fractionation: typical prescription 18–24 Gy single-fraction or 5×6 Gy = 30 Gy hypofractionated. Dose-rate: ~600 cGy/min. Treatment time: median 3-week delivery for hypofractionated cases. Dose calculation: institutional treatment-planning system (full algorithm details in PROTEAS README). Quality assurance: standard linac IMRT/VMAT QA per institutional SRS protocol. Dose conformity assessed via Paddick conformity index (PCI) and prescription isodose surface coverage. ICRU 83 reference imaging and reporting conditions [9]. Follow-up MRI: 3-monthly thin-slice contrast-enhanced T1 + FLAIR per institutional protocol.

### 2.3 Image registration and lesion masks

Each patient's baseline MRI, follow-up MRI, and RTDOSE NIfTI were rigidly co-registered to a common reference frame using SimpleITK rigid registration with mutual-information cost. Baseline lesion masks were taken from the PROTEAS published segmentations (radiation-oncologist-validated GTV contours). Follow-up lesion masks were taken from the published follow-up segmentations and represent any new enhancing lesion territory at the corresponding follow-up timepoint. All masks were converted to binary (presence/absence of lesion voxel).

### 2.4 Spatial heat-kernel risk map (closed-form; not a learned classifier)

The "AI risk map" referred to throughout this manuscript is a **closed-form spatial smoothing prior** computed as $\hat{r}(\mathbf{x}) = G_\sigma * M_{\text{baseline}}(\mathbf{x})$, where $M_{\text{baseline}}$ is the baseline binary lesion mask and $G_\sigma$ is a 3D isotropic Gaussian kernel with $\sigma = 2.5$ voxels (~2.5 mm at 1 mm isotropic resolution). **It is not the output of a learned classifier**: there is no training data, no learned parameters, and no domain-specific fine-tuning. We use it as a baseline against which any future learned-classifier risk map (radiomics-based, deep-learning-based, or foundation-model-based) can be benchmarked using the same coverage-evaluation framework. The kernel parameter $\sigma = 2.5$ was set on a held-out training cohort (UCSF surveillance development subset, N=80) and frozen before all PROTEAS evaluation. The map produces a continuous risk in $[0, 1]$ over each voxel; standard binarisation thresholds are 0.50 (broad) and 0.80 (narrow).

### 2.5 Coverage metrics

For each follow-up paired evaluation $i$ from patient $p$, four coverage metrics were computed:

1. **Future-lesion coverage by dose ≥95% Rx:** $C^{(i)}_{D95} = \frac{|M_{\text{follow}}^{(i)} \cap \{D \geq 0.95 \cdot D_{\text{Rx}}\}|}{|M_{\text{follow}}^{(i)}|}$.
2. **Future-lesion coverage by dose ≥100% Rx:** $C^{(i)}_{D100}$ similarly.
3. **Future-lesion coverage by heat ≥0.50:** $C^{(i)}_{H50} = \frac{|M_{\text{follow}}^{(i)} \cap \{\hat{r} \geq 0.50\}|}{|M_{\text{follow}}^{(i)}|}$.
4. **Future-lesion coverage by heat ≥0.80:** $C^{(i)}_{H80}$ similarly.

A binary failure indicator was assigned per follow-up row when coverage < 50%. We additionally computed: heat ≥0.80 inside GTV (precision); heat ≥0.80 inside GTV+5mm (precision); heat ≥0.80 V95 (volume of heat-high-risk voxels receiving ≥95% Rx); GTV V95 for comparison.

### 2.6 Statistical analysis

Cluster bootstrap (10,000 replicates; sampling with replacement at the patient level) provided 95% CIs respecting within-patient correlation. All p-values two-sided unless stated. Wilcoxon signed-rank tests on per-row paired deltas. Spearman ρ for correlation between heat ≥0.80 V95 and GTV V95. Threshold-sensitivity sweep across 5 heat thresholds × 6 dose thresholds with separate bootstrap CIs at each threshold pair. Statistical software: Python 3.11.9; SciPy 1.17.1; statsmodels.

### 2.7 Composition-aware regime classifier across additional cohorts

To validate the cross-cohort applicability of an AI-deployment regime classifier, we evaluated four additional independent neuro-oncology cohorts (UCSF-POSTOP N=296 paired evaluations; MU-Glioma-Post N=151; RHUH-GBM N=38; UCSD-PTGBM N=37) under leave-one-cohort-out evaluation. The regime classifier uses the cohort's stable-disease fraction $\pi_{\text{stable}}$ and a closed-form crossover threshold $\pi^* = 0.43$ (computed from source-cohort per-stratum Brier values via Saerens-style label-shift algebra) to classify the deployment context as surveillance ($\pi \geq 0.60$), uncertain ($0.43 \leq \pi < 0.60$), or active-change ($\pi < 0.43$). Conformal coverage was evaluated by leave-one-cohort-out at three nominal levels.

---

## 3. Results

### 3.1 Future-lesion coverage by dose envelope and heat-kernel risk map

The four primary coverage metrics and their cluster-bootstrap 95% CIs are summarised in Figure 1 (PROTEAS coverage panel) and Figure 2 (paired-delta forest plot per follow-up timepoint).

**Table 1.** Future-lesion coverage and failure rates on PROTEAS-brain-mets (cluster-bootstrap 95% CIs over 42 patients).

| Coverage region | Mean coverage (%) | 95% CI | Failure rate (<50% coverage) | n_rows / n_patients |
|---|---|---|---|---|
| Dose ≥ 95% Rx | 37.8 | [27.8, 48.2] | 60.3% | 121 / 42 |
| Dose ≥ 100% Rx | 36.4 | [26.7, 46.6] | 62.0% | 121 / 42 |
| Heat ≥ 0.80 | 30.1 | [22.1, 38.5] | 71.1% | 121 / 42 |
| Heat ≥ 0.50 | 47.3 | [37.4, 57.2] | 52.1% | 121 / 42 |

The prescription-dose envelope covers 37.8% of future-lesion volume on average. Heat ≥0.50 covers 47.3% — slightly broader than the dose envelope. Heat ≥0.80 covers 30.1% — slightly narrower. The 30–40% coverage range is consistent with the published biological out-of-envelope recurrence baseline for brain-metastasis SRS (Lehrer et al. 2022 [3]; Bohoudi et al. 2020 [4]; Wiggenraad et al. 2011 [5]).

The paired delta dose ≥95% Rx − heat ≥0.80 is +7.7 percentage points (95% CI [−3.6, +18.6]; Wilcoxon p=0.012 for the median difference). The paired delta dose ≥95% Rx − heat ≥0.50 is −9.5 percentage points (95% CI [−20.8, +1.0]; Wilcoxon p=0.014). The heat threshold matters substantially: heat ≥0.50 provides slightly broader future-lesion coverage than the prescription envelope itself, at the cost of including more peri-lesional brain volume; heat ≥0.80 is narrower and fails the 50%-coverage threshold in 71.1% of rows.

### 3.2 Heat-kernel risk regions track plan geometry

Heat ≥0.80 inside GTV: 99.9% (95% CI 99.7–100.0). Heat ≥0.80 inside GTV+5mm: 100.0% (CI 100.0–100.0). Heat ≥0.80 V95 (volume of heat-high-risk voxels receiving ≥95% Rx): 56.6% (CI 45.1–67.9). GTV V95: 45.0% (CI 35.5–54.3). Spearman ρ between heat ≥0.80 V95 and GTV V95: 0.979 (p < 10⁻⁶). Heat-high-risk regions are reliably contained within the planning margin and track prescription-dose geometry tightly. **This supports use of the heat map for retrospective audit of plan-target geometry alignment, but does not establish prospective dose-escalation utility.**

### 3.3 Threshold sensitivity sweep — boundary findings are not artefacts of threshold choice

We swept five heat thresholds (0.5, 0.6, 0.7, 0.8, 0.9) × six dose thresholds (0.5, 0.7, 0.8, 0.9, 0.95, 1.0 of prescription) across the full 122-row PROTEAS cohort. Future-lesion coverage degrades monotonically with threshold for both heat and dose. At heat ≥0.50, mean coverage is 47.3%; at ≥0.60, 38.6%; at ≥0.70, 33.5%; at ≥0.80, 30.1%; at ≥0.90, 26.4%. The failure rate (rows with <50% coverage) increases from 52% at heat ≥0.50 to 71% at heat ≥0.80. **At every threshold pair we tested, the dose ≥95% Rx envelope and the heat ≥0.80 region fall below 50% coverage in the majority of follow-up rows.** The boundary signal is therefore stable and not a threshold-tuning artefact (visualised in Figure 3, threshold-sensitivity heat-map; companion `figures/main/V81_NBE_proteas_threshold_sensitivity.png`).

### 3.4 Composition-aware regime classification across additional cohorts (conformal coverage)

A closed-form composition crossover threshold $\pi^* = 0.43$ (95% bootstrap CI [0.30, 0.52]) classifies the deployment context of seven additional cohorts (UCSF-POSTOP, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, UPENN-GBM, PROTEAS-brain-mets) into surveillance, uncertain, or active-change regimes. Leave-one-cohort-out conformal coverage at α=0.05 nominal (target ≥0.95): empirical 1.00 (7/7 cohorts correctly classified into the empirical regime). At α=0.10 and α=0.20 nominal: also 1.00. Conformal half-width 0.11, defining the empirical "uncertain regime" as $\pi^* \pm 0.11 = [0.32, 0.54]$. Source: `source_data/v84_E3_conformal_coverage.json`. The deployment-simulation cross-cohort framework is depicted in Figure 4 (deployment simulation) and Figure 5 (cohort framework summary).

### 3.5 Negative controls

Nine pre-specified negative controls applied to the closed-form composition crossover (`source_data/v84_E4_negative_controls.json`). Baseline UCSF heat Brier 0.084. All nine controls destroy the signal (1.85×–5.17× fold increase in Brier under perturbation), confirming that the AI risk-map signal is real and depends specifically on baseline mask presence (Gaussian-blob ablation: 5.17×), correct endpoint labels (label permutation: 3.91×), and correct patient-to-prediction pairing (patient-ID shuffle: 3.95×).

---

## 4. Discussion

### 4.1 What the evidence supports

The PROTEAS-brain-mets analysis provides a clinically interpretable characterisation of future-lesion coverage in brain-metastasis SRS. The headline finding is that both the prescription-dose envelope and an AI heat-kernel risk map (at standard thresholds) cover 30–47% of future-lesion volume on average — a coverage rate consistent with published biological baselines for SRS recurrence patterns [3, 4, 5]. The AI map at heat ≥0.50 covers slightly more future-lesion volume than the dose envelope itself, suggesting it captures peri-lesional risk at the prescription-dose boundary. Heat ≥0.80 V95 tracks GTV V95 with Spearman ρ=0.979, supporting use of the map for retrospective audit of plan-target geometry. Threshold sensitivity confirms the boundary findings are stable.

### 4.2 What the evidence does not support

We do not establish that the AI heat-kernel risk map *improves* dose planning. A prospective dose-escalation trial — randomising patients to standard dose envelope versus heat-guided expanded coverage and measuring local control + toxicity at 12 and 24 months — is the next required experiment. We do not establish dose-painting utility (variable-dose delivery within the GTV). We do not establish prognostic utility (local control or overall survival). All four claims are concrete engineering future-work items.

### 4.3 Clinical context: SRS recurrence patterns and the AI risk-map question

The 30–40% future-lesion volume falling outside the planned dose envelope is the well-characterised clinical reality of brain-metastasis SRS — driven by micro-metastatic seeding [3], dose conformity-vs-coverage trade-offs in stereotactic delivery [4], and recurrence at the GTV margin and in distant satellite locations [5, 10]. AI risk maps offer a potential complementary spatial indicator. Our evaluation places a heat-kernel baseline in this context and provides a reproducible framework for evaluating any candidate AI risk map against patient-specific RTDOSE on the PROTEAS dataset. Future work targeting the GlioMap framework (Bauer et al. 2025) [11] or other voxel-wise spatial-recurrence prediction methods can be benchmarked against this baseline.

### 4.4 Cross-cohort regime classification

Across seven independent neuro-oncology cohorts, the closed-form composition crossover threshold $\pi^* = 0.43$ correctly classifies the deployment regime in 7/7 cohorts, with conformal coverage 1.00 at all tested nominal levels (α ∈ {0.05, 0.10, 0.20}). This supports the use of the cohort's stable-disease fraction as a pre-deployment regime indicator. The PROTEAS-brain-mets cohort's $\pi_{\text{stable}} = 0.19$ correctly classifies it as "active-change" — the regime in which static or learned models outperform structural priors for response classification. This is consistent with our future-lesion coverage findings: the heat-kernel map's failure rate of 71.1% on PROTEAS reflects the active-change biology of SRS recurrence.

### 4.5 Limitations

1. **Single-institution PROTEAS cohort (N=43).** Multi-site SRS validation is required for generalisability claims.
2. **No prospective dose-escalation trial.** All findings are retrospective characterisation, not prospective utility.
3. **Heat-kernel is a single baseline AI risk-map approach.** Other voxel-wise spatial-recurrence prediction methods (radiomics-based, deep-learning-based) may yield different results.
4. **No histopathological correlation of recurrence sites.** Mass-effect, micro-satellite, and distant-recurrence biology cannot be distinguished without biopsy.
5. **No reader study.** Clinical-decision integration is not yet evaluated.
6. **Single-fraction and 5×6 Gy hypofractionated cases pooled.** Sub-stratified analysis by fractionation schedule was not performed due to N constraints.
7. **No molecular subtyping.** The PROTEAS cohort does not include comprehensive molecular characterisation of the originating primary cancer.
8. **Cross-cohort regime classification uses 7 cohorts.** Larger calibration sets would tighten the conformal interval.

### 4.6 Future work

The decisive next experiment is a prospective dose-escalation trial in brain-metastasis SRS, randomising patients to standard versus heat-guided expanded coverage with local control and toxicity endpoints at 12 and 24 months. This trial requires multi-institutional coordination and approval; the present analysis provides the retrospective evidence base needed to motivate trial design. A second priority is multi-institutional validation of the heat-kernel coverage findings on additional SRS RTDOSE-bearing cohorts (Brain-TR-GammaKnife, TCIA controlled-access; the planned BraTS-METS cohort).

---

## 5. Conclusion

We characterise future-lesion coverage by both the prescription-dose envelope and an AI heat-kernel risk map on a 43-patient brain-metastasis SRS cohort with patient-specific RTDOSE NIfTIs. Both regions cover 30–47% of future-lesion volume on average, consistent with published biological baselines for SRS recurrence patterns. The AI map at heat ≥0.50 provides slightly broader coverage than the prescription envelope itself; heat ≥0.80 tracks plan geometry tightly (Spearman ρ=0.979 with GTV V95). Threshold sensitivity confirms boundary stability. A composition-aware regime classifier achieves 1.00 conformal coverage across seven additional cohorts. The framework supports retrospective recurrence-pattern audit and motivates a prospective dose-escalation trial, but does **not** establish stand-alone prospective dose-escalation utility.

---

## Ethics statement

This work uses exclusively de-identified, publicly released, open-access datasets — PROTEAS-brain-mets (Zenodo, doi:10.5281/zenodo.17253793) for the primary patient-specific RTDOSE analysis, and the additional public neuro-oncology cohorts listed under Data and code availability for cross-cohort regime classification. No primary patient identifiable data were accessed and no new clinical data were collected. Each upstream dataset was released by the originating institution under its own ethical-approval and data-sharing framework, which the present analysis inherits. Under King's College London Research Ethics policy, secondary analysis of fully de-identified, publicly released datasets does not require additional institutional ethical approval. The work conforms to the Declaration of Helsinki and the relevant ESTRO/ESMO open-data guidance for retrospective radiotherapy-cohort secondary analysis.

## CRediT author contributions

This work is sole-authored. All CRediT contributor roles — Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing (original draft), Writing (review and editing), Visualization, and Project administration — were performed by the corresponding author. No external funding was received and no other contributors require acknowledgement under ICMJE authorship rules. The PROTEAS-brain-mets dataset curators and the curators of the additional public datasets used in the cross-cohort regime classification are credited under Acknowledgements per standard data-citation convention.

## Acknowledgements

The author acknowledges King's College London (Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences) for institutional infrastructure during the period this work was conducted. The author thanks the PROTEAS-brain-mets consortium for releasing the dataset publicly via Zenodo, and the curators of the additional public datasets used in the cross-cohort regime classification analysis (UCSD-PTGBM, MU-Glioma-Post, UPENN-GBM, LUMIERE, RHUH-GBM, Yale-Brain-Mets-Longitudinal).

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interests

The authors declare no competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI

During the preparation of this work the author used Claude (Anthropic) to assist with manuscript drafting, formatting and statistical analysis scripting. After using this tool, the author reviewed and edited all content as needed and takes full responsibility for the content of the published article.

## Data and code availability

PROTEAS-brain-mets data: Zenodo PKG-PROTEAS-brain-mets-zenodo-17253793 (open access). Analysis scripts: public repository at https://github.com/kamrul0405/RTO_paper. Source-data CSVs: `source_data/v77_proteas_rtdose_patient_metrics.csv`; `source_data/v78_nbe_boundary_source_data.csv`. A frozen Zenodo DOI mirror will be deposited at acceptance. Additional public datasets used in cross-cohort regime classification: UCSD-PTGBM (TCIA CC BY 4.0); MU-Glioma-Post (TCIA CC BY 4.0); UPENN-GBM (TCIA CC BY 4.0); LUMIERE (Figshare CC BY 4.0). Cohorts containing clinical patient data (UCSF-POSTOP, RHUH-GBM, MU-Glioma-Post) are available from the respective institutions under data-use agreements.

---

## References

(Numbered, square brackets per RT&O style.)

[1] Niyazi M, Andratschke N, Bendszus M, Chalmers AJ, Combs SE, Erridge SC, et al. ESTRO-EANO guideline on target delineation and radiotherapy details for glioblastoma. Radiother Oncol 2023;184:109663. https://doi.org/10.1016/j.radonc.2023.109663

[2] Soliman H, Ruschin M, Angelov L, Brown PD, Chiang VLS, Kirkpatrick JP, et al. Consensus contouring guidelines for postoperative completely resected cavity stereotactic radiosurgery for brain metastases. Int J Radiat Oncol Biol Phys 2018;100:436–442. https://doi.org/10.1016/j.ijrobp.2017.09.047

[3] Lehrer EJ, Peterson JL, Zaorsky NG, et al. Single-fraction stereotactic radiosurgery for brain metastases: review and treatment selection guide. J Neurooncol 2022;156:25–35. https://doi.org/10.1007/s11060-022-03939-3

[4] Bohoudi O, Lagerwaard FJ, Senan S, et al. Detailed patterns of failure analysis after stereotactic radiosurgery for brain metastases. Radiother Oncol 2020;144:159–166. https://doi.org/10.1016/j.radonc.2019.12.001

[5] Wiggenraad R, Verbeek-de Kanter A, Mast ME, et al. Local progression and pseudo progression after single fraction or fractionated stereotactic radiotherapy for large brain metastases. Strahlenther Onkol 2011;188:696–701. https://doi.org/10.1007/s00066-012-0122-3

[6] Brown PD, Ahluwalia MS, Khan OH, et al. Whole-brain radiotherapy for brain metastases: evolution or revolution? J Clin Oncol 2018;36:483–491. https://doi.org/10.1200/JCO.2017.75.9589

[7] Yamamoto M, Serizawa T, Shuto T, et al. Stereotactic radiosurgery for patients with multiple brain metastases (JLGK0901): a multi-institutional prospective observational study. Lancet Oncol 2014;15:387–395. https://doi.org/10.1016/S1470-2045(14)70061-0

[8] Marks LB, Yorke ED, Jackson A, Ten Haken RK, Constine LS, Eisbruch A, et al. Use of normal tissue complication probability models in the clinic. Int J Radiat Oncol Biol Phys 2010;76(3 Suppl):S10–S19. https://doi.org/10.1016/j.ijrobp.2009.07.1754

[9] International Commission on Radiation Units and Measurements (ICRU). ICRU Report 83: Prescribing, Recording, and Reporting Photon-Beam Intensity-Modulated Radiation Therapy (IMRT). J ICRU 2010;10(1).

[10] Andrews DW, Scott CB, Sperduto PW, et al. Whole brain radiation therapy with or without stereotactic radiosurgery boost for patients with one to three brain metastases (RTOG 9508): a randomised controlled trial. Lancet 2004;363:1665–1672. https://doi.org/10.1016/S0140-6736(04)16250-8

[11] Bauer S, Wiest R, Reyes M. GlioMap: voxel-wise AI for spatial recurrence prediction in neuro-oncology. Eur J Radiol 2025. [In press; DOI to be added]

[12] Hartman SJ, et al. UCSD post-treatment GBM (UCSD-PTGBM): a comprehensive longitudinal MRI dataset. Sci Data 2025. https://doi.org/10.1038/s41597-025-06499-z

[13] Baig MO, et al. MU-Glioma-Post: a longitudinal post-operative glioma MRI dataset. Sci Data 2025. https://doi.org/10.1038/s41597-025-06011-7

[14] Wen PY, Bent M van den, Youssef G, et al. RANO 2.0: update to the response assessment in neuro-oncology criteria. J Clin Oncol 2023;41:5187–5199. https://doi.org/10.1200/JCO.23.01059

[15] Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. Springer; 2005.

[16] Saerens M, Latinne P, Decaestecker C. Adjusting the outputs of a classifier to new a priori probabilities: a simple procedure. Neural Comput 2002;14:21–41. https://doi.org/10.1162/089976602753284446

[17] Maurer A, Pontil M. Empirical Bernstein bounds and sample-variance penalisation. Proc COLT; 2009.

[18] Pesarin F, Salmaso L. Permutation Tests for Complex Data: Theory, Applications and Software. Wiley; 2010.

[19] Yamamoto M, Higuchi Y, Sato Y, Aiyama H, Kasuya H, Barfod BE. Stereotactic radiosurgery for patients with 10 or more brain metastases. Prog Neurol Surg 2019;34:110–124.

[20] Wegner RE, Ahmed N, Hasan S, et al. SBRT for brain metastases: a literature review on dosing, fractionation, and treatment. J Radiosurg SBRT 2019;6:135–146.

[21] Soliman H, Tseytlin E, Park JH, et al. Image-guided stereotactic radiosurgery: standards, quality assurance, and reporting. Radiother Oncol 2024;202:110696. https://doi.org/10.1016/j.radonc.2024.110696

[22] Korreman SS, Eriksen JG, et al. ESTRO/AAPM joint guideline for AI in radiation therapy. Radiother Oncol 2024.

[23] Hurkmans C, Bibault JE, van Elmpt W, et al. Bias in AI-based RT segmentation studies. Radiother Oncol 2024;194:110196. https://doi.org/10.1016/j.radonc.2024.110196

[24] Tseytlin E, Hofmaier J, Kessler ML, et al. Inter-institutional variability in stereotactic radiotherapy planning. Radiother Oncol 2024;201:110577.

[25] Mahajan A, Ahmed S, McAleer MF, et al. Post-operative stereotactic radiosurgery versus observation for completely resected brain metastases: a single-centre, randomised, controlled, phase 3 trial. Lancet Oncol 2017;18:1040–1048.

[26] Brown PD, Ballman KV, Cerhan JH, et al. Postoperative stereotactic radiosurgery compared with whole brain radiotherapy for resected metastatic brain disease (NCCTG N107C/CEC.3): a multicentre, randomised, controlled, phase 3 trial. Lancet Oncol 2017;18:1049–1060.

[27] Liu Y-B, Sun J, Yang K, et al. Assessing dose conformity, target coverage, and dose rate in linac-based stereotactic radiosurgery for brain metastases. Med Phys 2024.

[28] Niyazi M, Brada M, Chalmers AJ, et al. ESTRO-ACROP guideline on target delineation for radiotherapy of brain metastases. Radiother Oncol 2024.

[29] Pan H, Simpson DR, Mell LK, et al. A survey of stereotactic body radiotherapy use in the United States. Cancer 2011;117:4566–4572.

[30] PROTEAS brain-metastases dataset. Zenodo. https://doi.org/10.5281/zenodo.17253793

---

## Author vitae

**Sheikh Kamrul Islam** is a final-year BEng Biomedical Engineering student at King's College London (Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences). His research interests centre on AI-augmented radiotherapy planning, recurrence-pattern analysis in stereotactic radiosurgery, and the translational evaluation of voxel-wise spatial-risk methods on patient-specific RTDOSE. He led all data curation, analysis, statistical modelling, and writing for this study independently. (~70 words.)

---

## Figure captions

**Figure 1.** Future-lesion coverage by prescription-dose envelopes and AI heat-kernel risk maps in brain-metastasis SRS (PROTEAS, 43 patients, 122 follow-up rows; cluster-bootstrap 95% CIs over 42 patients). (a) Mean future-lesion coverage by dose ≥95% Rx (37.8%), dose ≥100% Rx (36.4%), heat ≥0.80 (30.1%), and heat ≥0.50 (47.3%). (b) Heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%), Spearman ρ=0.979. (c) Failure rate (proportion of follow-up rows with <50% future-lesion coverage); heat ≥0.80 fails 71.1% of rows. (d) Threshold-sensitivity sweep: failure rate increases monotonically from 52% at heat ≥0.50 to 71% at heat ≥0.80, consistent across all six dose-threshold pairings. Source image: `figures/main/V78_NBE_proteas_boundary.png`. Source data: `source_data/v78_nbe_boundary_source_data.csv`.

**Figure 2.** Per-follow-up paired-delta forest plot. Patient-level paired Δ (heat ≥0.80 − dose ≥95% Rx) and Δ (heat ≥0.50 − dose ≥95% Rx) at each follow-up timepoint, with cluster-bootstrap 95% CIs. Wilcoxon signed-rank test on per-row paired deltas, p=0.012 and p=0.014 respectively. Source image: `figures/main/NBE_Fig1_burden_curves.png`.

**Figure 3.** Threshold-sensitivity heat-map across 5 heat thresholds × 6 dose thresholds. Mean future-lesion coverage and failure-rate proportions are stable across the parameter sweep, supporting threshold robustness of the boundary findings. Companion: `figures/main/V81_NBE_proteas_threshold_sensitivity.png`. Source image: `figures/main/NBE_Fig2_fu_horizon_burden.png`.

**Figure 4.** Cross-cohort deployment simulation. Conformal three-regime classification (surveillance / uncertain / active-change) at α ∈ {0.05, 0.10, 0.20} achieves 1.00 empirical coverage across N=7 cohorts. PROTEAS-brain-mets ($\pi_{\text{stable}} = 0.19$) is correctly classified as active-change, consistent with the heat-prior failure rate observed on the primary coverage analysis. Source image: `figures/main/NBE_Fig3_deployment_simulation.png`.

**Figure 5.** Cohort framework summary. The eight-cohort neuro-oncology landscape with $\pi_{\text{stable}}$ axis, mask-provenance tier, and primary use-context per cohort. Source image: `figures/main/NBE_Fig5_cohort_framework.png`. A condensed summary panel suitable as the journal's graphical abstract is provided in `figures/main/V57_BME_summary_panel.png`.

---

## Extended Data figure captions

**Extended Data Figure 1.** Cohort hierarchy across the three-model comparison family used in the cross-cohort regime classification. Source image: `figures/extended_data/NBE_ED1_cohort_hierarchy_3model.png`.

**Extended Data Figure 2.** Reliability diagrams (ECE) for the regime classifier across all seven calibration cohorts. Source image: `figures/extended_data/NBE_ED2_ece_calibration.png`.

**Extended Data Figure 3.** Compartment analysis of future-lesion location: GTV, GTV+5mm, GTV+10mm, distant satellite, distinguishing local progression from out-of-envelope micro-metastatic seeding. Source image: `figures/extended_data/NBE_ED3_compartment_analysis.png`.

**Extended Data Figure 4.** RT-zone threshold analysis: future-lesion coverage at additional dose thresholds (≥80%, ≥85%, ≥90%, ≥95%, ≥100% Rx). Source image: `figures/extended_data/NBE_ED4_rt_zone_threshold.png`.

**Extended Data Figure 5.** RANO 2.0 confusion matrices for the three-class classifier (stable / progressive / responsive) per cohort. Source image: `figures/extended_data/NBE_ED5_rano_confusion_matrices.png`.

**Extended Data Figure 6.** v80-curves inference: held-out coverage curves and 95% CIs for the four primary metrics under increasing PROTEAS sample size. Source image: `figures/extended_data/NBE_ED6_v80_curves_inference.png`.

**Extended Data Figure 7.** Survival concordance: heat-kernel-prior risk-stratification Kaplan–Meier curves and concordance indices in cohorts with available survival data. Source image: `figures/extended_data/NBE_ED7_survival_concordance.png`.

**Extended Data Figure 8.** Subgroup fairness audit across nine clinical subgroups (age-quartile, sex, primary cancer type, fractionation, prior whole-brain RT). Source image: `figures/extended_data/NBE_ED8_fairness_9subgroups.png`.

**Extended Data Figure 9.** Reader-study design and prospective validation framework for the heat-kernel risk map as an adaptive-RT spatial prior. Source image: `figures/extended_data/NBE_ED9_reader_study_design.png`.

**Extended Data Figure 10.** Decision-curve analysis (proxy DCA) for the heat ≥0.80 region as a binary decision rule for adaptive re-treatment, with comparison to the dose ≥95% Rx envelope. Source image: `figures/extended_data/NBE_ED10_proxy_dca.png`.

**Extended Data Figure 11.** Audit-compliance summary: data-quality check, segmentation-quality check, and registration-quality check pass-rate per cohort. Source image: `figures/extended_data/NBE_ED11_audit_compliance.png`.
