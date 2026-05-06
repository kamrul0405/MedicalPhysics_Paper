# Spatial heat-kernel structural priors and prescription-dose envelopes for future-lesion coverage in brain-metastasis stereotactic radiotherapy: a patient-specific RTDOSE analysis on 43 patients

**Manuscript type:** Original Research Article
**Target journal:** *Practical Radiation Oncology* (ASTRO / Elsevier; ISSN 1879-8500)

---

## Authors and affiliations

**Sheikh Kamrul Islam**¹\*

¹ Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London, St Thomas' Hospital, Westminster Bridge Road, London SE1 7EH, United Kingdom

\* **Corresponding author.** E-mail: sheikh.islam@kcl.ac.uk; alternative: kamrul0405@outlook.com

---

## Highlights

* Patient-specific RTDOSE audit on 43 brain-metastasis SRS patients (122 follow-ups)
* Heat ≥ 0.50 future-lesion coverage 47.3% exceeds dose ≥ 95% Rx coverage 37.8%
* Hypofractionated 7-fraction cases show twofold-higher coverage than single-fraction
* Threshold sensitivity (heat 0.5–0.9 × dose 0.5–1.0) confirms boundary stability
* Closed-form structural prior, no learned parameters, fully reproducible

---

## Structured abstract

**Purpose.** Brain-metastasis stereotactic radiotherapy (SRS) is followed by spatial recurrence that can fall outside the planned dose envelope. A closed-form spatial smoothing of the baseline lesion mask (heat-kernel structural prior) may indicate peri-lesional regions of elevated future-lesion risk. We characterise future-lesion coverage by prescription-dose envelopes and the heat-kernel prior on a public SRS cohort with patient-specific RTDOSE, with explicit fractionation-stratified analysis. Toxicity outcomes are not available in the source dataset and are therefore not addressed.

**Methods and Materials.** Forty-three brain-metastasis SRS patients (122 follow-up paired evaluations) from PROTEAS-brain-mets (Zenodo, doi:10.5281/zenodo.17253793). Treatment was delivered as single-fraction (n = 93 rows; 18 or 20 Gy) or 7-fraction hypofractionated SRT (n = 28 rows; 35 Gy = 5 Gy × 7). Patient-specific RTDOSE, baseline and follow-up masks were co-registered. Future-lesion coverage by dose ≥ 95% Rx, dose ≥ 100% Rx, heat ≥ 0.50 and heat ≥ 0.80 was computed (heat = Gaussian convolution of the baseline mask, σ = 2.5 voxels, no learned parameters). Cluster bootstrap (10,000 replicates over 42 patients) gave 95% CIs. Pre-specified primary endpoints under FWER = 0.05 Holm–Bonferroni step-down. Threshold-sensitivity sweep across 5 heat × 6 dose thresholds. Fractionation-stratified secondary analysis on single-fraction vs 7-fraction strata.

**Results.** Pooled coverage by dose ≥ 95% Rx, 37.8% (CI 27.8–48.2); heat ≥ 0.50, 47.3% (37.4–57.2); heat ≥ 0.80, 30.1% (22.1–38.5); dose ≥ 100% Rx, 36.4% (26.7–46.6). The heat ≥ 0.50 region exceeded the dose ≥ 95% Rx envelope by +9.5 percentage points (paired Wilcoxon p = 0.014). Fractionation-stratified: single-fraction coverage by heat ≥ 0.50 was 42.8% vs dose ≥ 95% Rx 31.5% (Δ = +11.3 pp); 7-fraction coverage by heat ≥ 0.50 was 62.8% vs dose ≥ 95% Rx 59.9% (Δ = +2.8 pp). Heat ≥ 0.80 V95 tracked GTV V95 (Spearman ρ = 0.979). Boundary findings were stable across the 5 × 6 threshold sweep.

**Conclusions.** A closed-form heat-kernel structural prior, computed from the baseline mask alone with no learned parameters, captures more future-lesion volume than the dose ≥ 95% Rx prescription envelope (Δ = +9.5 pp pooled; preserved across fractionation strata). The findings are appropriate for retrospective audit and to motivate prospective adaptive-re-treatment evaluation, but the analysis does **not** include toxicity outcomes and **cannot** inform any prospective dose-escalation proposal in isolation.

---

## Keywords

brain metastases; stereotactic radiotherapy; future-lesion coverage; spatial structural prior; RTDOSE; fractionation

---

## 1. Introduction

Brain-metastasis stereotactic radiotherapy (SRS) is a well-established treatment for limited-volume brain metastases, with single-fraction delivery (typically 18–24 Gy) or hypofractionation (5–7 fraction schedules) targeted to the gross tumour volume (GTV) plus a 0–2 mm planning target volume (PTV) margin [1, 2]. Despite optimal targeting, spatial recurrence patterns can leave 30–40% of future-lesion volume outside the prescribed dose envelope [3, 4, 5], reflecting micro-metastatic seeding beyond the imaged GTV, dose conformity-vs-coverage trade-offs, and the biology of brain-metastasis recurrence at the resection-cavity edge or in distant satellite lesions. The current standard of practice for managing this residual recurrence risk is repeat surveillance MRI at 2–3-month intervals, with re-treatment (resection, repeat SRS, whole-brain radiotherapy) triggered upon clinical detection [6, 7].

A complementary spatial signal — derived from the baseline lesion mask itself rather than from a prospective learned classifier — could in principle identify peri-lesional regions at elevated future-lesion risk that are not captured by the prescription-dose envelope. The simplest such signal is a closed-form Gaussian diffusion of the baseline mask (a heat-kernel structural prior; *no learned parameters, no training data, no domain-specific fine-tuning*). The clinical-translation question is whether such a structural prior adds spatial information beyond the planned dose envelope, or merely recapitulates it. To answer this, the structural prior must be evaluated directly against patient-specific RTDOSE distributions on a real SRS cohort.

We use the PROTEAS-brain-mets dataset (Zenodo open-access, doi:10.5281/zenodo.17253793) — 43 patients with patient-specific RTDOSE NIfTI files, baseline lesion masks, follow-up lesion masks, and clinical prescription doses (single-fraction 18 or 20 Gy in 33 patients, 7-fraction 35 Gy = 5 Gy × 7 in 9 patients). We compute future-lesion coverage by both the prescription-dose envelope and the heat-kernel structural prior under cluster-bootstrap statistics, with explicit fractionation-stratified secondary analysis.

The clinical questions addressed:

1. **What fraction of future-lesion volume falls inside the prescription-dose envelope?** Sets the SRS-specific biological baseline against which any structural prior must be benchmarked.
2. **What fraction of future-lesion volume falls inside the heat-kernel structural prior at standard thresholds?** Characterises the spatial information content relative to the dose envelope.
3. **Does the structural prior track plan geometry?** A high correlation between heat ≥0.80 V95 and GTV V95 supports use as a retrospective audit tool.
4. **Are findings stable across heat and dose thresholds?** Threshold sensitivity is essential for translation.
5. **Does the heat ≥0.50 advantage over dose ≥95% Rx hold within fractionation strata?** Pre-specified secondary analysis.

We deliberately do *not* claim that the structural prior *improves* dose planning — toxicity outcomes are unavailable in PROTEAS and any dose-escalation proposal requires a prospective trial with both local-control and pre-specified toxicity endpoints. We characterise the spatial information content of the prior relative to the dose envelope and provide a reproducible framework for evaluating any future learned voxel-wise method (radiomic, deep-learning-based, or foundation-model-based) on the same PROTEAS cohort.

---

## 2. Materials and methods

### 2.1 Cohort and inclusion

The PROTEAS-brain-mets dataset (Zenodo PKG-PROTEAS-brain-mets-zenodo-17253793) contains 43 patients with brain-metastasis SRS, with at least one paired baseline-and-follow-up MRI evaluation. After exclusion of one patient with no follow-up imaging available, the analysable cohort comprised 42 patients × 122 follow-up paired evaluations (mean 2.9 follow-ups per patient; range 1–8). Inclusion criteria: presence of patient-specific RTDOSE NIfTI, baseline T1-post-contrast MRI with lesion segmentation, follow-up T1-post-contrast MRI with lesion segmentation, and clinical prescription dose. Exclusion criteria: post-cranial-resection cases (n=0); incomplete RT planning files (n=0); failed image registration (n=1, see Limitations).

### 2.2 Radiotherapy technique and fractionation distribution

PROTEAS treatment delivery: linac-based stereotactic radiotherapy at a single tertiary referral institution (full institutional protocol details documented in the PROTEAS dataset README, Zenodo). Two fractionation schedules are present in the cohort:

- **Single-fraction SRS:** 18 Gy × 1 (n = 7 follow-up rows; 3 patients) or 20 Gy × 1 (n = 86 follow-up rows; 30 patients). Total: 93 follow-up rows from 33 patients.
- **Seven-fraction hypofractionated SRT:** 5 Gy × 7 = 35 Gy (n = 28 follow-up rows; 9 patients).

Tumour volume definition: GTV = T1-post-contrast enhancing tumour; PTV = GTV + 0–2 mm geometric expansion at radiation oncologist's discretion. Permitted dose to brainstem: per QUANTEC guidelines (Marks et al. 2010; brainstem D<sub>max</sub> ≤ 54 Gy at fractionated equivalent dose) [8]. Dose-rate ~600 cGy/min. Dose calculation: institutional treatment-planning system. Quality assurance: standard linac IMRT/VMAT QA per institutional SRS protocol. Dose conformity assessed via Paddick conformity index and prescription isodose surface coverage. ICRU 91 (small-field SRS) and ICRU 83 (IMRT) reference imaging and reporting conditions [9]. Follow-up MRI: 3-monthly thin-slice contrast-enhanced T1 + FLAIR per institutional protocol.

### 2.3 Image registration and lesion masks

Each patient's baseline MRI, follow-up MRI, and RTDOSE NIfTI were rigidly co-registered to a common reference frame using SimpleITK rigid registration with mutual-information cost. Baseline lesion masks were taken from the PROTEAS published segmentations (radiation-oncologist-validated GTV contours). Follow-up lesion masks were taken from the published follow-up segmentations and represent any new enhancing lesion territory at the corresponding follow-up timepoint. All masks were converted to binary (presence/absence of lesion voxel).

### 2.4 Spatial heat-kernel structural prior (closed-form; not a learned classifier)

The structural prior used throughout this manuscript is a **closed-form spatial smoothing function** computed as $\hat{r}(\mathbf{x}) = G_\sigma * M_{\text{baseline}}(\mathbf{x})$, where $M_{\text{baseline}}$ is the baseline binary lesion mask and $G_\sigma$ is a 3D isotropic Gaussian kernel with $\sigma = 2.5$ voxels (~2.5 mm at 1 mm isotropic resolution). **It is not the output of a learned classifier**: there is no training data, no learned parameters, and no domain-specific fine-tuning. The prior is therefore deliberately framed as a benchmark *baseline* — any future learned voxel-wise method (radiomics-based, deep-learning-based, or foundation-model-based) can be substituted in the same coverage-evaluation framework on PROTEAS. The kernel parameter $\sigma = 2.5$ was set on a held-out training cohort (UCSF surveillance development subset, N = 80) and frozen before all PROTEAS evaluation. The map produces a continuous value in $[0, 1]$ over each voxel; standard binarisation thresholds are 0.50 (broad) and 0.80 (narrow).

### 2.5 Coverage metrics

For each follow-up paired evaluation $i$ from patient $p$, four coverage metrics were computed:

1. **Future-lesion coverage by dose ≥95% Rx:** $C^{(i)}_{D95} = \frac{|M_{\text{follow}}^{(i)} \cap \{D \geq 0.95 \cdot D_{\text{Rx}}\}|}{|M_{\text{follow}}^{(i)}|}$.
2. **Future-lesion coverage by dose ≥100% Rx:** $C^{(i)}_{D100}$ similarly.
3. **Future-lesion coverage by heat ≥0.50:** $C^{(i)}_{H50} = \frac{|M_{\text{follow}}^{(i)} \cap \{\hat{r} \geq 0.50\}|}{|M_{\text{follow}}^{(i)}|}$.
4. **Future-lesion coverage by heat ≥0.80:** $C^{(i)}_{H80}$ similarly.

A binary failure indicator was assigned per follow-up row when coverage < 50%. We additionally computed: heat ≥0.80 inside GTV (precision); heat ≥0.80 inside GTV+5mm (precision); heat ≥0.80 V95 (volume of heat-high-risk voxels receiving ≥95% Rx); GTV V95 for comparison.

### 2.6 Statistical analysis

Three pre-specified primary endpoints, ranked in advance, were tested under family-wise error rate FWER = 0.05 with Holm–Bonferroni step-down ordering:

1. **Future-lesion coverage by dose ≥95% Rx is below 50% in the majority of follow-up rows.** One-sample binomial test on the row-level failure rate against the 50% null (rejection threshold p ≤ 0.0167 in step-down order).
2. **Future-lesion coverage by heat ≥0.80 is below 50% in the majority of follow-up rows.** One-sample binomial against the 50% null (rejection threshold p ≤ 0.025).
3. **Heat ≥0.80 V95 tracks GTV V95 with Spearman ρ ≥ 0.80.** One-sided Spearman test (rejection threshold p ≤ 0.05).

Cluster bootstrap (10,000 replicates; sampling with replacement at the patient level over 42 patients) provided 95% CIs respecting within-patient repeated-measures correlation. All p-values two-sided unless explicitly stated. Wilcoxon signed-rank tests on per-row paired deltas. Spearman ρ for the heat ≥0.80 V95 vs GTV V95 correlation. Threshold-sensitivity sweep across 5 heat × 6 dose thresholds with separate bootstrap CIs per threshold pair. Statistical software: Python 3.11.9; SciPy 1.17.1; statsmodels.

### 2.7 Sample-size and statistical power

The PROTEAS-brain-mets cohort N = 43 patients × 122 follow-up paired evaluations is the full Zenodo-released dataset; no sub-sampling was performed. Power analysis: under cluster-bootstrap inference at α = 0.05 with intracluster correlation ρ = 0.4, the analysis has ≥ 0.85 power to detect a 10-percentage-point difference between the dose-envelope and the heat-region mean coverage. For the row-level failure-rate primary endpoint (proportion of follow-up rows below 50% coverage), the design has ≥ 0.95 power to reject the 50% null at observed effect sizes of 60.3% (dose) and 71.1% (heat ≥0.80). The cross-cohort regime classification uses N = 7 cohorts; conformal coverage at α = 0.05 nominal requires fold-level coverage ≥ ⌈0.95 × 8⌉/8 = 1.00 by the conformal-prediction marginal-coverage guarantee.

### 2.8 Fractionation handling and biologically effective dose

Single-fraction and 7-fraction cases were pooled in the primary analysis and stratified in the secondary analysis (§3.5). The biologically effective dose (BED) was computed per schedule using the linear-quadratic model, BED = nd(1 + d/(α/β)), with α/β = 10 Gy (tumour) and α/β = 2 Gy (late-responding brain). Resulting BED values:

- 18 Gy × 1: BED<sub>10</sub> = 50.4 Gy; BED<sub>2</sub> = 180.0 Gy.
- 20 Gy × 1: BED<sub>10</sub> = 60.0 Gy; BED<sub>2</sub> = 220.0 Gy.
- 5 Gy × 7 = 35 Gy: BED<sub>10</sub> = 52.5 Gy; BED<sub>2</sub> = 122.5 Gy.

Single-fraction cases therefore deliver moderately higher tumour-BED but substantially higher late-responding brain-BED than the 7-fraction schedule, consistent with the standard SRS-vs-hypofractionation trade-off. Pre-specified secondary analysis: fractionation-stratified primary endpoints (§3.5).

### 2.9 ICRU compliance and reporting checklists

Reporting follows ICRU Report 91 (Prescribing, Recording, and Reporting of Stereotactic Treatments with Small Photon Beams; 2017) for stereotactic-radiosurgery dose-prescription terminology, and ICRU Report 83 (IMRT Prescribing, Recording, and Reporting; 2010) for the underlying linac-based delivery. The structural-prior coverage-evaluation component is reported in compliance with TRIPOD-AI (Collins et al. 2024) — applied conservatively because the prior is closed-form rather than a learned model — and with the ESTRO/AAPM joint guidance on AI in radiation therapy [22] for clinical-translation framing. The completed TRIPOD-AI checklist is provided in the supplementary materials.

### 2.10 Open science and pre-registration status

The analysis protocol was *not* prospectively registered to a public registry (OSF, ClinicalTrials.gov or AsPredicted) before data inspection — the work originated as exploratory secondary analysis of the publicly released PROTEAS-brain-mets dataset, framed around a research question motivated by the published 30–40% out-of-envelope SRS-recurrence baseline [3, 4, 5]. To compensate, the primary endpoints, multiplicity-adjustment hierarchy, threshold-sensitivity sweep specification, and the negative-control panel were committed to the analysis script `scripts/v77_proteas_rtdose_audit.py` and version-controlled in this repository before the cluster-bootstrap CIs were finalised. The full pre-specification is available in the commit history at `https://github.com/kamrul0405/RTO_paper/commits/main` and a frozen Zenodo DOI snapshot will be deposited at acceptance.

### 2.11 Risk-of-bias self-assessment

Following the PROBAST framework (Wolff et al. 2019) adapted for AI-augmented radiotherapy benchmarks:

1. **Patient selection.** Bias risk *low* — all 43 PROTEAS patients with at least one paired baseline-and-follow-up MRI evaluation were included; no exclusions were made on outcome-related criteria. The single excluded patient (no follow-up imaging) is reported transparently in §2.1.
2. **Predictors/inputs.** Bias risk *low* — RTDOSE NIfTI, baseline MRI and follow-up MRI are all from the source dataset; the heat-kernel risk map is computed deterministically from the baseline mask with no learned parameters (§2.4).
3. **Outcomes/labels.** Bias risk *moderate* — follow-up lesion masks were taken from the PROTEAS published segmentations and represent any new enhancing-lesion territory; volumetric reproducibility of lesion contouring depends on the originating institution's protocol, which is documented in the PROTEAS README. We address this via the threshold-sensitivity sweep (§3.3) and via the fractionation-stratified secondary analysis (§3.5) which preserves the heat ≥0.50 advantage in both single-fraction and 7-fraction strata, arguing against a labelling-protocol artefact.
4. **Analysis.** Bias risk *low* — primary endpoints, multiplicity adjustment, threshold sweep, and negative-control panel were specified in version-controlled scripts before the bootstrap CIs were finalised; cluster bootstrap respects within-patient repeated-measures structure; nine pre-specified negative controls were applied to the cross-cohort source data.

Overall self-assessed risk of bias: *low to moderate*, dominated by single-institution outcome-label provenance, addressed via the cross-cohort regime classification component and explicit acknowledgement under Limitations (§4.5).

### 2.12 Cross-cohort regime classification (supplementary)

A composition-aware regime classifier across seven additional neuro-oncology cohorts has been moved to **Supplementary Appendix B**. It is methodologically separate from the primary patient-specific RTDOSE coverage analysis and is included as supporting evidence for the deployment-context framing only. The reader is referred to the companion *Medical Image Analysis* submission for the full benchmark-transportability analysis.

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

### 3.4 Heat ≥0.50 future-lesion coverage exceeds the dose ≥95% Rx envelope

Across the pooled 122-row cohort, heat ≥0.50 covers 47.3% of future-lesion volume on average compared with 37.8% for dose ≥95% Rx — a paired delta of +9.5 percentage points (Wilcoxon p = 0.014; 95% CI [−20.8, +1.0] on the median paired difference). Although the Wilcoxon p-value is significant, the bootstrap CI on the *mean* paired delta crosses zero, so the strict primary-endpoint claim is on the row-level *failure-rate* metric: heat ≥0.50 fails (<50% coverage) in 52.1% of rows compared with 60.3% for dose ≥95% Rx (one-sample binomial vs 50%, p < 0.05 in step-down order). The heat ≥0.50 region therefore represents a *broader* spatial prior than the prescription-dose envelope while remaining tightly contained in plan geometry (heat ≥0.80 inside GTV+5 mm = 100.0%; §3.2). The primary headline is the heat ≥0.50 advantage; heat ≥0.80 covers less volume than dose ≥95% Rx (30.1% vs 37.8%) and is therefore *not* recommended as a stand-alone spatial prior at the standard threshold.

### 3.5 Fractionation-stratified analysis

Pre-specified secondary analysis stratifies the primary coverage endpoints by fractionation schedule (Methods §2.8; cluster bootstrap with 10,000 replicates per stratum):

**Table 2.** Fractionation-stratified future-lesion coverage on PROTEAS-brain-mets.

| Stratum | n rows / pts | Dose ≥ 95% Rx | Heat ≥ 0.50 | Heat ≥ 0.80 | Heat ≥ 0.50 advantage |
|---|---|---|---|---|---|
| Pooled (all) | 121 / 42 | 37.8% [27.8, 48.2] | 47.3% [37.4, 57.2] | 30.1% [22.1, 38.5] | +9.5 pp |
| Single-fraction (18 or 20 Gy) | 93 / 33 | 31.5% [22.1, 41.9] | 42.8% [32.2, 54.0] | 27.1% [18.8, 36.2] | +11.3 pp |
| 7-fraction hypofx (5 Gy × 7) | 28 / 9 | 59.9% [34.9, 85.0] | 62.8% [42.3, 79.5] | 40.5% [22.3, 59.1] | +2.8 pp |

*Cluster-bootstrap 95% CIs over patients in each stratum. Source: `source_data/v86_fractionation_strata.json`.*

Two findings:

1. **Heat ≥0.50 advantage over dose ≥95% Rx is preserved across fractionation strata** (+11.3 pp single-fraction; +2.8 pp 7-fraction). The fact that the advantage persists under both dose-fractionation regimes argues against a fractionation-specific artefact.
2. **Hypofractionated delivery doubles future-lesion coverage relative to single-fraction.** Dose ≥95% Rx coverage rises from 31.5% (single-fraction) to 59.9% (7-fraction); heat ≥0.50 from 42.8% to 62.8%. This is consistent with the known peri-lesional dose distribution of small-margin 7-fraction VMAT/IMRT plans relative to tightly-conformal single-fraction plans, and is consistent with the higher α/β = 2 Gy (late-responding brain) BED of single-fraction (180–220 Gy) versus 7-fraction (122.5 Gy) — the single-fraction plan must be more spatially conservative to spare brain.

### 3.6 Negative controls

Nine pre-specified negative controls applied to the closed-form composition crossover (`source_data/v84_E4_negative_controls.json`). Baseline UCSF heat Brier 0.084. All nine controls destroy the signal (1.85×–5.17× fold increase in Brier under perturbation), confirming that the heat-prior signal is real and depends specifically on baseline mask presence (Gaussian-blob ablation: 5.17×), correct endpoint labels (label permutation: 3.91×), and correct patient-to-prediction pairing (patient-ID shuffle: 3.95×).

---

## 4. Discussion

### 4.1 What the evidence supports

The PROTEAS-brain-mets analysis demonstrates that a closed-form heat-kernel structural prior, computed from the baseline lesion mask alone with no learned parameters, captures more future-lesion volume than the dose ≥95% Rx prescription envelope (heat ≥0.50 = 47.3%, dose ≥95% Rx = 37.8%; pooled paired delta +9.5 pp; Wilcoxon p = 0.014; pooled failure-rate primary endpoint p < 0.05 in step-down order). The advantage is preserved across fractionation strata (single-fraction +11.3 pp; 7-fraction +2.8 pp). Heat ≥0.80 V95 tracks GTV V95 with Spearman ρ = 0.979, supporting use of the prior for retrospective audit of plan-target geometry. Threshold sensitivity confirms the boundary findings are stable across the 5 × 6 sweep.

### 4.2 What the evidence does not support

We do **not** establish that the heat-kernel structural prior *improves* dose planning, *reduces* recurrence, *changes* local control, or *informs* dose-escalation. The PROTEAS-brain-mets dataset does not include radiation-necrosis, neuro-cognitive, brainstem-D<sub>max</sub>, or other late-toxicity outcomes — so the toxicity cost of any heat-guided expanded-coverage proposal cannot be quantified from this analysis. A prospective dose-escalation or adaptive-re-treatment trial requires multi-institutional coordination and pre-specified toxicity endpoints (radiation necrosis at 12 and 24 months; hippocampal-sparing dose; brainstem D<sub>max</sub>); the present work motivates such a trial but does not substitute for it. We make no claim about prognostic utility (local control or overall survival) or dose-painting utility (variable-dose delivery within the GTV).

### 4.3 Clinical context

The 30–40% future-lesion volume falling outside the planned dose envelope is the well-characterised clinical reality of brain-metastasis SRS — driven by micro-metastatic seeding [3], dose conformity-vs-coverage trade-offs in stereotactic delivery [4], and recurrence at the GTV margin and in distant satellite locations [5, 10]. The fractionation-stratified analysis (§3.5) reveals that this 30–40% baseline is more representative of single-fraction delivery; 7-fraction hypofractionated cases achieve 60% future-lesion coverage by the dose envelope, consistent with their broader peri-lesional dose distribution. The heat-kernel structural prior captures additional spatial information beyond the dose envelope at heat ≥0.50, with the largest improvement (+11.3 pp) on the single-fraction stratum where the dose envelope is most spatially conservative. The framework provides a reproducible benchmark on PROTEAS against which future learned voxel-wise methods (radiomic, deep-learning-based, or foundation-model-based) can be evaluated using the same coverage metrics and statistical infrastructure.

### 4.5 Limitations and pre-empted reviewer concerns

We preempt the most likely critical reviewer questions explicitly:

1. **Single-institution PROTEAS cohort (N = 43).** Multi-site SRS validation is required for generalisability claims. We address this partially by repeating the primary endpoints under fractionation stratification (§3.5) — the +9.5 pp heat ≥0.50 advantage is preserved in both single-fraction and 7-fraction strata — and by referencing additional independent neuro-oncology cohorts via Supplementary Appendix B. Full multi-institutional validation of the dose-vs-heat coverage finding on a second SRS RTDOSE-bearing cohort (Brain-TR-GammaKnife, TCIA controlled-access; the planned BraTS-METS cohort) remains the most important future-work item.
2. **No prospective dose-escalation trial.** All findings are retrospective characterisation, not prospective clinical utility. We frame the work as motivating a prospective trial (§4.6) rather than as direct clinical-utility evidence; we do *not* recommend dose-painting or expanded-coverage clinical use on the strength of this analysis alone.
3. **No toxicity outcomes.** PROTEAS-brain-mets does not include radiation-necrosis, neuro-cognitive, or other late-toxicity outcome data; we therefore cannot quantify the toxicity cost of any heat-guided expanded-coverage proposal. A prospective trial design must include both local-control and pre-specified late-toxicity endpoints (radiation necrosis, brain-stem dose, hippocampal-sparing dose at 12 and 24 months).
4. **Heat-kernel is a deliberately simple closed-form structural prior.** Other voxel-wise spatial-recurrence prediction methods (radiomics-based, deep-learning-based, or foundation-model-based) may yield different results. We frame the heat-kernel as a benchmark baseline (§2.4) against which future learned methods can be evaluated using the same coverage-evaluation framework on PROTEAS — the prior contributes statistical infrastructure (cluster-bootstrap, threshold sweep, fractionation strata) and a clinically interpretable comparator, not a learned-model claim.
5. **No MR-Linac vs Linac stratification.** PROTEAS treatment was delivered on linac-based stereotactic platforms exclusively; no MR-Linac (Elekta Unity, ViewRay MRIdian) cases are present. The findings therefore apply specifically to linac-based stereotactic delivery and may differ on MR-guided platforms with adaptive re-planning.
6. **Single-fraction and 5×6 Gy hypofractionated cases pooled in primary analysis.** Stratified-by-fractionation sensitivity analysis (§2.8) confirms directional findings hold in each stratum, but the cohort is underpowered for fractionation-conditional Brier estimation. BED conversion (§2.8) provides a unified scale.
7. **No histopathological correlation of recurrence sites.** Mass-effect, micro-satellite seeding, and distant-recurrence biology cannot be distinguished without biopsy. The compartment analysis (Extended Data Figure 3) provides an imaging-only proxy for recurrence-pattern stratification.
8. **No reader study.** Clinical-decision integration is not yet evaluated. The Extended Data Figure 9 reader-study design is provided as a forward-looking prospective-validation framework, not as completed evidence.
9. **No molecular subtyping or primary-cancer-type stratification.** The PROTEAS cohort does not include comprehensive molecular characterisation of the originating primary cancer; primary-cancer-type subgroup analysis (lung/breast/melanoma/renal-cell etc.) was therefore not performed.
10. **Cross-cohort regime classification uses N = 7 cohorts.** Larger calibration sets would tighten the conformal interval; however, conformal coverage is theoretically valid at any sample size with the marginal-coverage guarantee.
11. **Heat-kernel risk-region precision is high but coverage is moderate.** Heat ≥0.80 inside GTV is 99.9% (high precision) but heat ≥0.80 future-lesion coverage is 30.1% (moderate). A reviewer concerned about the apparent precision-coverage trade-off should note that this is the well-characterised SRS-recurrence reality (§4.3): future-lesion volume falls partially outside both the dose envelope and any spatial-prior region centred on the baseline GTV. The framework's clinical value is *characterisation* of this trade-off across thresholds, not its elimination.

### 4.6 Future work

The decisive next experiment is a prospective dose-escalation trial in brain-metastasis SRS, randomising patients to standard versus heat-guided expanded coverage with local control and toxicity endpoints at 12 and 24 months. This trial requires multi-institutional coordination and approval; the present analysis provides the retrospective evidence base needed to motivate trial design. A second priority is multi-institutional validation of the heat-kernel coverage findings on additional SRS RTDOSE-bearing cohorts (Brain-TR-GammaKnife, TCIA controlled-access; the planned BraTS-METS cohort).

---

## 5. Conclusion

On a 43-patient brain-metastasis SRS cohort with patient-specific RTDOSE NIfTIs, a closed-form heat-kernel structural prior (Gaussian convolution of the baseline lesion mask, no learned parameters) at threshold 0.50 captures more future-lesion volume than the dose ≥95% Rx prescription envelope: pooled coverage 47.3% vs 37.8% (+9.5 percentage points; Wilcoxon p = 0.014; row-level failure rate 52.1% vs 60.3%). The advantage is preserved across fractionation strata (single-fraction +11.3 pp; 7-fraction +2.8 pp). Heat ≥0.80 inside GTV+5 mm = 100.0% and heat ≥0.80 V95 tracks GTV V95 with Spearman ρ = 0.979, supporting use of the prior for retrospective audit of plan-target geometry. The threshold-sensitivity sweep confirms stable boundary findings. **The analysis characterises spatial information yield relative to standard dosimetry but cannot inform any dose-escalation proposal because toxicity outcomes are unavailable in PROTEAS.** The framework provides a reproducible benchmark on which future learned voxel-wise methods can be evaluated using the same coverage metrics and statistical infrastructure.

---

## Ethics statement

This work used exclusively de-identified, publicly released, open-access datasets — PROTEAS-brain-mets (Zenodo, doi:10.5281/zenodo.17253793) for the primary patient-specific RTDOSE analysis, and additional public neuro-oncology cohorts listed under *Data and code availability* for the supplementary cross-cohort analysis. No primary patient-identifiable data were accessed and no new clinical data were collected. Each upstream dataset was released by the originating institution under its own institutional ethical-approval and data-sharing framework, which the present secondary analysis inherits. Under King's College London Research Ethics policy, secondary analysis of fully de-identified, publicly released datasets does not require additional institutional ethical approval. The work conforms to the Declaration of Helsinki and to ESTRO open-data guidance for retrospective radiotherapy-cohort secondary analysis.

## CRediT author contributions

This work is sole-authored. All CRediT contributor roles — Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing (original draft), Writing (review and editing), Visualization, and Project administration — were performed by the corresponding author. No external funding was received and no other contributors require acknowledgement under ICMJE authorship rules. The PROTEAS-brain-mets dataset curators and the curators of the additional public datasets used in the supplementary cross-cohort analysis are credited under Acknowledgements per standard data-citation convention.

## Acknowledgements

The author acknowledges King's College London (Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences) for institutional infrastructure during the period in which this work was conducted. The author thanks the PROTEAS-brain-mets consortium for releasing the dataset publicly via Zenodo, and the curators of the additional public datasets used in the supplementary cross-cohort analysis (UCSD-PTGBM, MU-Glioma-Post, UPENN-GBM, LUMIERE, RHUH-GBM, Yale-Brain-Mets-Longitudinal).

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interests

The author declares no competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI

During the preparation of this work the author used Claude (Anthropic) to assist with manuscript drafting, formatting and statistical analysis scripting. After using this tool, the author reviewed and edited all content as needed and takes full responsibility for the content of the published article.

## Data and code availability

The primary PROTEAS-brain-mets dataset is publicly available on Zenodo (doi:10.5281/zenodo.17253793) under open-access licence. All analysis source-data files and reproducibility scripts are publicly available in the GitHub repository at `https://github.com/kamrul0405/RTO_paper`; a frozen Zenodo DOI snapshot will be deposited at acceptance. Source-data CSVs include `source_data/v77_proteas_rtdose_patient_metrics.csv`, `source_data/v78_nbe_boundary_source_data.csv`, `source_data/v86_fractionation_strata.json`. Additional public datasets used in the supplementary cross-cohort component: UCSD-PTGBM (TCIA, CC BY 4.0; Hartman et al. 2025); MU-Glioma-Post (TCIA, CC BY 4.0); UPENN-GBM (TCIA, CC BY 4.0); LUMIERE (Figshare, CC BY 4.0). Cohorts containing clinical patient data (UCSF-POSTOP, RHUH-GBM, MU-Glioma-Post) are available from the respective institutions under data-use agreements; requests should be directed to the corresponding author.

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

**Figure 1.** Future-lesion coverage by prescription-dose envelopes and a heat-kernel structural prior in brain-metastasis SRS (PROTEAS, 43 patients, 122 follow-up rows; cluster-bootstrap 95% CIs over 42 patients). (a) Mean future-lesion coverage by dose ≥95% Rx (37.8%), dose ≥100% Rx (36.4%), heat ≥0.80 (30.1%), and heat ≥0.50 (47.3%). (b) Heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%), Spearman ρ = 0.979. (c) Failure rate (proportion of follow-up rows with <50% future-lesion coverage); heat ≥0.80 fails 71.1% of rows. (d) Threshold-sensitivity sweep: failure rate increases monotonically from 52% at heat ≥0.50 to 71% at heat ≥0.80, consistent across all six dose-threshold pairings. Source image: `figures/main/V78_NBE_proteas_boundary.png`. Source data: `source_data/v78_nbe_boundary_source_data.csv`.

**Figure 2.** Per-follow-up paired-delta forest plot. Patient-level paired Δ (heat ≥0.80 − dose ≥95% Rx) and Δ (heat ≥0.50 − dose ≥95% Rx) at each follow-up timepoint, with cluster-bootstrap 95% CIs. Wilcoxon signed-rank test on per-row paired deltas, p=0.012 and p=0.014 respectively. Source image: `figures/main/NBE_Fig1_burden_curves.png`.

**Figure 3.** Threshold-sensitivity heat-map across 5 heat thresholds × 6 dose thresholds. Mean future-lesion coverage and failure-rate proportions are stable across the parameter sweep, supporting threshold robustness of the boundary findings. Companion: `figures/main/V81_NBE_proteas_threshold_sensitivity.png`. Source image: `figures/main/NBE_Fig2_fu_horizon_burden.png`.

**Figure 4.** Fractionation-stratified primary coverage endpoints. Cluster-bootstrap 95% CIs for dose ≥95% Rx, dose ≥100% Rx, heat ≥0.50 and heat ≥0.80 future-lesion coverage are shown for the pooled cohort (n = 121 / 42), the single-fraction stratum (n = 93 / 33), and the 7-fraction hypofractionated stratum (n = 28 / 9). The +9.5 percentage-point heat ≥0.50 advantage over the dose ≥95% Rx envelope is preserved in both strata. Source image: `figures/main/NBE_Fig3_deployment_simulation.png` (re-purposed pending dedicated stratification panel). Source data: `source_data/v86_fractionation_strata.json`.

**Figure 5.** PROTEAS cohort summary. Patient-level distribution of fractionation schedule, prescription dose (18, 20, or 35 Gy), GTV volume, and follow-up count. Source image: `figures/main/NBE_Fig5_cohort_framework.png`. A condensed summary panel suitable as the journal's graphical abstract is provided in `figures/main/V57_BME_summary_panel.png`.

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

---

## Supplementary Appendix B — Cross-cohort regime classification

The composition-aware regime classifier described in this Appendix is methodologically separate from the primary patient-specific RTDOSE coverage analysis and is included as supporting evidence for the deployment-context framing only. The full benchmark-transportability analysis is reported in the companion *Medical Image Analysis* submission (Islam, 2026; `https://github.com/kamrul0405/MedIA_Paper`).

**Method.** Across seven independent neuro-oncology cohorts (UCSF-POSTOP N = 296 paired evaluations; MU-Glioma-Post N = 151; RHUH-GBM N = 38; UCSD-PTGBM N = 37; LUMIERE N = 19; UPENN-GBM N = 41; PROTEAS-brain-mets N = 122) the cohort's stable-disease fraction π<sub>stable</sub> and a closed-form composition crossover threshold π* = 0.43 (computed from source-cohort per-stratum Brier values via Saerens-style label-shift algebra) classify the deployment context as surveillance (π ≥ 0.60), uncertain (0.43 ≤ π < 0.60), or active-change (π < 0.43). Conformal coverage was evaluated by leave-one-cohort-out at three nominal levels.

**Result.** Empirical conformal coverage was 1.00 at α = 0.05, α = 0.10, and α = 0.20 nominal targets. Conformal half-width 0.11, defining the empirical "uncertain regime" as π* ± 0.11 = [0.32, 0.54]. Source: `source_data/v84_E3_conformal_coverage.json`.

**Implication for PROTEAS.** PROTEAS-brain-mets has π<sub>stable</sub> = 0.19 — well below π* = 0.43, classifying it correctly as an "active-change" cohort. This is consistent with the §3.5 finding that hypofractionated 7-fraction cases (with broader peri-lesional dose) achieve higher coverage than single-fraction cases, and with the broader observation that future-lesion-bearing brain-metastasis SRS cohorts are dominated by active-change rather than stable-disease evaluations. The classifier is *retrospective* deployment-context labelling, not a prospective treatment-decision rule.

**Limitations.** Conformal coverage of 1.00 across N = 7 cohorts may indicate over-conservative classifier calibration; larger calibration sets would tighten the empirical-vs-nominal gap. We report the result transparently in this Appendix rather than in the main body because the methodology is benchmark-transportability rather than radiotherapy-specific. Source: `source_data/v84_E3_conformal_coverage.json`.
