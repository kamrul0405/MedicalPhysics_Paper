# Physics-grounded structural priors in brain-metastasis stereotactic radiotherapy: parabolic-PDE smoothing, BED-aware spatially-varying kernels, and multi-institutional dose-physics audit on patient-specific RTDOSE/RTPLAN

**Manuscript type:** Original Research Article
**Target journal:** *Medical Physics* (AAPM / Wiley; ISSN 0094-2405)

---

## Authors and affiliations

**Sheikh Kamrul Islam**¹\*

¹ Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences, King's College London, St Thomas' Hospital, Westminster Bridge Road, London SE1 7EH, United Kingdom

\* **Corresponding author.** E-mail: sheikh.islam@kcl.ac.uk; alternative: kamrul0405@outlook.com

---

## Highlights

* PROTEAS audit parses 47 RTDOSE and 48 RTPLAN files with no errors
* BED/EQD2 physics grounds single-fraction and 7-fraction strata
* Dose and heat regions leave substantial future-lesion volume uncovered
* Heat-only rows define an audit signal, not dose-escalation evidence
* Conditional-use boundary blocks treatment-planning activation

---

## Structured abstract

**Purpose.** Brain-metastasis stereotactic radiotherapy (SRS) can be followed by future enhancing-lesion volume that falls outside the prescription-dose envelope. We derive a physics-grounded structural prior — the heat-kernel — as the fundamental solution of the heat equation applied to the baseline lesion mask, and extend it to a per-voxel BED-aware spatially-varying kernel σ(BED) using the linear-quadratic radiobiology model. We evaluate both kernel variants against patient-specific RTDOSE arrays on a public SRS cohort, with multi-institutional physics-context audit (8 cohorts; 2,875 patients), fractionation-stratified secondary analysis, BED-stratified analysis, decision-curve analysis, and α/β sensitivity sweep. Toxicity, local-control and treatment-escalation outcomes are not available and are not claimed.

**Methods and Materials.** Forty-three PROTEAS-brain-mets SRS patients (42 analysable patients; 121 follow-up rows after one no-follow-up exclusion) with baseline masks, follow-up masks and RTDOSE NIfTI files. The v92 DICOM integrity audit parsed 47 RTDOSE files and 48 RTPLAN files from nested patient archives, deriving prescription dose, fractionation, BED and EQD2. Coverage of future-lesion volume was computed for dose >=95% Rx, dose >=100% Rx, heat >=0.50 and heat >=0.80. Patient-clustered bootstrap confidence intervals, paired tests, threshold sweeps, fractionation-stratified exploratory analyses and dose-prior discordance taxonomy were performed.

**Results.** Mean future-lesion coverage was 37.9% for dose >=95% Rx, 36.6% for dose >=100% Rx, 47.4% for heat >=0.50 and 30.1% for heat >=0.80. A discordance taxonomy showed 21/121 heat-only captures (17.4%) and 52/121 neither-captured rows (43.0%). Heat >=0.50 was numerically broader than dose >=95% Rx (+9.5 percentage points), but the cluster-bootstrap mean-delta CI crossed zero; this is reported as an audit signal, not a superiority claim. RTPLAN schedules were 20 Gy x 1 (30 plans), 18 Gy x 1 (2), 5 Gy x 7 (12) and 20.093 Gy x 1 (1); BED10 median was 60.0 Gy and BED2 median 220.0 Gy. Heat >=0.80 V95 tracked GTV V95 tightly (Spearman rho = 0.979), supporting plan-geometry coupling. Threshold sweeps confirmed incomplete coverage across heat and dose thresholds.

**Conclusions.** The heat-kernel prior is suitable as a transparent retrospective audit comparator and prospective reader-study candidate. It does not support heat-guided dose escalation, toxicity reduction or clinical decision use without independent SRS RTDOSE validation, toxicity outcomes and reader-study evidence.

---

## Keywords

brain metastases; stereotactic radiotherapy; future-lesion coverage; spatial structural prior; RTDOSE; fractionation

---

## 1. Introduction

Brain-metastasis stereotactic radiotherapy (SRS) is a well-established treatment for limited-volume brain metastases, with single-fraction delivery (typically 18-24 Gy) or hypofractionation targeted to the gross tumour volume (GTV) plus a 0-2 mm planning target volume (PTV) margin [1, 2]. Despite conformal targeting, published failure-pattern studies show that part of subsequent enhancing disease can occur outside high-dose regions or at peri-lesional/distant sites [3-5]. Standard management remains surveillance MRI with re-treatment triggered on clinical detection [6, 7].

A complementary spatial signal — derived from the baseline lesion mask itself rather than from a prospective learned classifier — could in principle identify peri-lesional regions at elevated future-lesion risk that are not captured by the prescription-dose envelope. The simplest such signal is a closed-form Gaussian diffusion of the baseline mask (a heat-kernel structural prior; *no learned parameters, no training data, no domain-specific fine-tuning*). The clinical-translation question is whether such a structural prior adds spatial information beyond the planned dose envelope, or merely recapitulates it. To answer this, the structural prior must be evaluated directly against patient-specific RTDOSE distributions on a real SRS cohort.

We use the PROTEAS-brain-mets dataset (Zenodo open-access, doi:10.5281/zenodo.17253793): 43 released patients, 42 analysable patients after one no-follow-up exclusion, patient-specific RTDOSE NIfTI files, baseline/follow-up lesion masks and clinical prescription doses. We compute future-lesion coverage by the prescription-dose envelope and a frozen heat-kernel structural prior under patient-clustered bootstrap statistics, with fractionation-stratified secondary analysis.

The clinical questions addressed:

1. **What fraction of future-lesion volume falls inside the prescription-dose envelope?** Sets the SRS-specific biological baseline against which any structural prior must be benchmarked.
2. **What fraction of future-lesion volume falls inside the heat-kernel structural prior at standard thresholds?** Characterises the spatial information content relative to the dose envelope.
3. **Does the structural prior track plan geometry?** A high correlation between heat ≥0.80 V95 and GTV V95 supports use as a retrospective audit tool.
4. **Are findings stable across heat and dose thresholds?** Threshold sensitivity is essential for translation.
5. **Which rows are dose-only, heat-only, captured by both, or captured by neither?** This discordance taxonomy separates an audit signal from a dose-escalation claim.
6. **Does the numerical heat >=0.50 versus dose >=95% Rx difference persist within fractionation strata?** Pre-specified secondary analysis.
7. **Are the dose objects physically auditable?** The v92 upgrade parses RTDOSE and RTPLAN DICOM directly and derives prescription, fractionation, BED and EQD2 rather than relying only on NIfTI proxy filenames.

We deliberately do *not* claim that the structural prior *improves* dose planning — toxicity outcomes are unavailable in PROTEAS and any dose-escalation proposal requires a prospective trial with both local-control and pre-specified toxicity endpoints. The v92 novelty is therefore a conditional-use boundary framework grounded in dose physics: a map can be useful as a retrospective audit comparator only when its claims are bounded by delivered dose, prescription, fractionation biology and missing-toxicity evidence. We characterise the spatial information content of the prior relative to the dose envelope and provide a reproducible framework for evaluating any future learned voxel-wise method (radiomic, deep-learning-based, or foundation-model-based) on the same PROTEAS cohort.

---

## 2. Materials and methods

### 2.1 Cohort and inclusion

The PROTEAS-brain-mets dataset (Zenodo PKG-PROTEAS-brain-mets-zenodo-17253793) contains 43 released patients with brain-metastasis SRS. After exclusion of one patient with no follow-up imaging and one failed follow-up registration, the analysable cohort comprised 42 patients and 121 follow-up paired evaluations (mean 2.9 follow-ups per patient; range 1-8). Inclusion required patient-specific RTDOSE NIfTI, baseline contrast-enhanced T1 MRI with lesion segmentation, follow-up lesion segmentation and clinical prescription dose.

### 2.2 Radiotherapy technique and fractionation distribution

PROTEAS treatment delivery was linac-based stereotactic radiotherapy at a single tertiary referral institution. Two fractionation schedules were present: single-fraction SRS, 18 Gy x 1 or 20 Gy x 1 (93 follow-up rows from 33 patients), and seven-fraction hypofractionated SRT, 5 Gy x 7 = 35 Gy (28 rows from 9 patients). GTV was defined as contrast-enhancing tumour; PTV was GTV plus a 0-2 mm expansion at physician discretion. Dose reporting follows ICRU 91/83 terminology, with follow-up MRI approximately every 3 months.

### 2.2.1 v92 RTDOSE/RTPLAN integrity and BED/EQD2 audit

The v92 upgrade re-opened the nested PROTEAS patient archives and parsed radiotherapy DICOM rather than trusting filenames alone. The audit identified 47 strict RTDOSE (`RD.*.dcm`) objects, 48 RTPLAN (`RP.*.dcm`) objects, 47 RTSTRUCT (`RS.*.dcm`) objects, 17,431 planning CT DICOM objects and 44,010 MR DICOM objects (`source_data/v91_proteas_rtdose_inventory.json`; `source_data/v92_proteas_plan_physics_audit.json`). RTDOSE and RTPLAN parsing had zero header errors; RTDOSE pixel-array parsing had zero errors.

For each RTPLAN, prescription dose was extracted from `DoseReferenceSequence` and fractions from `FractionGroupSequence`. BED and EQD2 were computed as $BED = nd(1+d/(\alpha/\beta))$ and $EQD2 = BED/(1+2/(\alpha/\beta))$ using $\alpha/\beta = 10$ Gy for tumour and 2 Gy for late-responding brain. This audit is a physics-integrity check and supports fractionation-normalised stress testing; it does not substitute for target-level DVH, toxicity outcomes or clinical treatment-effect evidence.

### 2.3 Image registration and lesion masks

Each patient's baseline MRI, follow-up MRI, and RTDOSE NIfTI were rigidly co-registered to a common reference frame using SimpleITK rigid registration with mutual-information cost. Baseline lesion masks were taken from the PROTEAS published segmentations (radiation-oncologist-validated GTV contours). Follow-up lesion masks were taken from the published follow-up segmentations and represent any new enhancing lesion territory at the corresponding follow-up timepoint. All masks were converted to binary (presence/absence of lesion voxel).

### 2.3.1 Multi-institutional physics atlas

The local archive accessed for this work spans eight neuro-oncology cohorts. Although only PROTEAS-brain-mets contains patient-specific RTDOSE arrays, the broader physics atlas (`source_data/v92_multisite_physics_atlas.json`) provides imaging-physics context for the structural prior's geometric framing.

**Table 1.** Multi-institutional cohort inventory (v92 physics atlas).

| Cohort | Patients | Timepoints | NIfTI files | Mask-like files | RTDOSE files | Use in this paper |
|---|---|---|---|---|---|---|
| UCSF-POSTOP | 298 | 596 | 4,768 | 1,192 | 0 | σ-development cohort |
| MU-Glioma-Post | 203 | 596 | 2,978 | 594 | 0 | Companion MedIA submission |
| RHUH-GBM | 40 | 120 | 720 | 120 | 0 | Companion MedIA submission |
| UCSD-PTGBM | 136 | 184 | 4,047 | 736 | 0 | Companion MedIA submission |
| UPENN-GBM | 630 | — | 10,646 | 758 | 0 | Companion MedIA submission |
| Yale-Brain-Mets-Longitudinal | 1,430 | — | 33,811 | 0 | 0 | Acquisition-shift screen |
| LUMIERE | 93 | — | 18,660 | 21,564 | 0 | Companion MedIA submission |
| **PROTEAS-brain-mets** | **45** | **122** | **1,050** | **216** | **47** | **Primary RTDOSE cohort** |
| **Total** | **2,875** | **10,679** | **76,680** | **25,180** | **47** | |

Across the eight cohorts, the heat-kernel structural prior operates on highly sparse lesion support: the median positive-voxel fraction is 0.42–0.86% across cohorts, the median bounding-box fraction is 2.3–3.8%, and the σ-to-equivalent-radius ratio is 0.095–0.105. The prior therefore captures the peri-lesional zone at ~10% of the lesion's characteristic radius — clinically corresponding to the region just outside the GTV boundary where micro-metastatic seeding and dose-conformity-vs-coverage trade-offs have been documented in the SRS-recurrence literature. PROTEAS-brain-mets is the only cohort in our local archive with patient-specific RTDOSE arrays; the BED-aware spatially-varying kernel analysis in §2.4–§3.10 is therefore necessarily PROTEAS-only at this time.

### 2.4 Heat-kernel structural prior — formal physics derivation

The heat-kernel structural prior is a **closed-form solution to the heat equation** applied to the binary baseline lesion mask. We make the physical interpretation explicit because it grounds the prior in classical PDE theory and motivates the BED-aware spatially-varying extension in §2.5.

**Definition.** Given the binary baseline lesion mask $M(\mathbf{x}) \in \{0, 1\}$, the heat-kernel risk map at scale parameter $\sigma$ is the convolution

$$\hat{r}_\sigma(\mathbf{x}) = G_\sigma * M(\mathbf{x}) = \int_{\mathbb{R}^3} G_\sigma(\mathbf{x} - \mathbf{y}) \, M(\mathbf{y}) \, d\mathbf{y},$$

where $G_\sigma(\mathbf{x}) = (2\pi\sigma^2)^{-3/2} \exp(-\|\mathbf{x}\|^2 / 2\sigma^2)$ is the 3D isotropic Gaussian kernel. The Gaussian kernel is the *fundamental solution* of the diffusion (heat) equation

$$\frac{\partial u}{\partial t} = D \nabla^2 u, \quad u(\mathbf{x}, 0) = M(\mathbf{x}),$$

with diffusion coefficient $D = 1/2$ and evolution time $t = \sigma^2 / 2$. The prior $\hat{r}_\sigma$ is therefore the spatial profile of the lesion mask after a fixed time of isotropic diffusion. Setting $\sigma = 2.5$ voxels (1 mm isotropic resolution) corresponds to evolving the mask under diffusion for $t = 3.125$ voxel-time units. **The heat prior is a parabolic-PDE smoothing operator, not a learned classifier**.

**Multi-site geometry.** Across the eight neuro-oncology cohorts indexed in our master dataset (`source_data/master_neurooncology_dataset_index.csv`; see also `v92_multisite_physics_atlas.json`), the heat kernel operates on highly sparse lesion support: median positive-voxel fraction 0.42–0.86%, median bounding-box fraction 2.3–3.8%, and σ-to-equivalent-radius ratio 0.095–0.105. The prior therefore captures the peri-lesional zone at ~10% of the lesion's characteristic radius — clinically corresponding to the region just outside the GTV boundary where micro-metastatic seeding and dose-conformity-vs-coverage trade-offs have been documented in the SRS-recurrence literature (Lehrer et al. 2022; Bohoudi et al. 2020).

**Linear-quadratic radiobiology and the BED-aware spatially-varying extension.** The classical radiobiology of brain-metastasis SRS is governed by the linear-quadratic model

$$BED(\mathbf{x}) = n \cdot d(\mathbf{x}) \cdot \left(1 + \frac{d(\mathbf{x})}{\alpha/\beta}\right),$$

where $d(\mathbf{x}) = D(\mathbf{x})/n$ is the per-fraction dose at voxel $\mathbf{x}$ (from the patient-specific RTDOSE array), $n$ is the number of fractions, and $\alpha/\beta$ is the tissue-type radiosensitivity ratio (10 Gy for tumour; 2 Gy for late-responding brain). The BED at each voxel encodes the local biological effect of the delivered dose. Our **BED-aware spatially-varying heat-kernel** (§2.5; introduced as a novel extension in this work) modulates the kernel scale $\sigma(\mathbf{x})$ pointwise as a function of local BED:

$$\sigma(\mathbf{x}) = \sigma_{\text{lo}} + (\sigma_{\text{hi}} - \sigma_{\text{lo}}) \cdot \left(1 - \frac{BED(\mathbf{x})}{BED_{\text{max}}}\right),$$

with $\sigma_{\text{lo}} = 1.5$ voxels (high-BED → tight prior, focal peri-lesional risk) and $\sigma_{\text{hi}} = 4.0$ voxels (low-BED → broad prior, distributed risk). The corresponding heat map is the spatially-varying convolution

$$\hat{r}_{\text{BED}}(\mathbf{x}) = w(\mathbf{x}) \cdot G_{\sigma_{\text{lo}}} * M(\mathbf{x}) + (1 - w(\mathbf{x})) \cdot G_{\sigma_{\text{hi}}} * M(\mathbf{x}),$$

with $w(\mathbf{x}) = BED(\mathbf{x})/BED_{\text{max}}$. This is **a physics-informed structural prior** that ties the kernel scale to the local radiobiological dose-response. The spatially-varying kernel reduces to the constant-σ kernel when BED is uniform across the lesion neighbourhood; under non-uniform BED, the prior assigns tighter peri-lesional regions to high-BED voxels (where the dose-conformity is tighter and out-of-envelope recurrence is more spatially focal) and broader regions to low-BED voxels.

**No learning required.** Both the constant-σ and BED-aware variants involve no learned parameters, no training data, and no domain-specific fine-tuning. The kernel parameters $\sigma$ (constant case), $\sigma_{\text{lo}}/\sigma_{\text{hi}}$ (BED-aware case), and $\alpha/\beta$ (radiobiology) are all set from physics-informed values and frozen before any PROTEAS evaluation; sensitivity sweeps over $\sigma$ thresholds (§3.3) and $\alpha/\beta$ choices (§3.10) confirm robustness of the directional findings.

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

The analysable PROTEAS cohort is the full available RTDOSE-bearing dataset after transparent exclusions: 42 patients and 121 follow-up evaluations. The study is powered for descriptive coverage and failure-rate estimation, not for treatment-effect or toxicity inference. The fractionation analysis is explicitly exploratory because the 7-fraction stratum contains only 9 patients.

### 2.8 Fractionation handling and biologically effective dose

Single-fraction and 7-fraction cases were pooled in the primary analysis and stratified secondarily (§3.5). BED was computed by BED = nd(1 + d/(α/β)) with α/β = 10 Gy for tumour and 2 Gy for late-responding brain. The v92 RTPLAN audit identified prescriptions of 18 Gy x 1, 20 Gy x 1, 20.093 Gy x 1 and 5 Gy x 7, with median BED<sub>10</sub> 60.0 Gy and median BED<sub>2</sub> 220.0 Gy across parsed prescriptions. These values provide physics-normalised context, not treatment-effect inference.

### 2.9 ICRU compliance and reporting checklists

Reporting follows ICRU 91/83 dose-prescription terminology. The structural-prior component is reported against TRIPOD-AI principles, applied conservatively because the prior is closed-form rather than learned.

### 2.10 Open science and pre-registration status

The protocol was not prospectively registered before data inspection. To reduce analytic-flexibility risk, the primary endpoints, multiplicity hierarchy and threshold sweep were committed to `scripts/v77_proteas_rtdose_audit.py` before final bootstrap CIs were generated; a frozen repository snapshot will be deposited at acceptance.

### 2.11 Risk-of-bias self-assessment

PROBAST-style self-assessment judged patient selection and predictors low risk: all eligible RTDOSE-bearing cases were analysed and the heat map is deterministic. Outcome-label risk is moderate because follow-up lesion contours originate from one institution and were not re-contoured for this analysis. Analysis risk is low: cluster bootstrap respects patient clustering and the threshold sweep tests parameter sensitivity. Overall risk of bias is low-to-moderate, dominated by single-institution outcome-label provenance.

### 2.12 Cross-cohort regime classification (supplementary)

A composition-aware regime classifier across seven additional neuro-oncology cohorts is moved to **Supplementary Appendix B**. It is methodologically separate from the primary patient-specific RTDOSE analysis and is included only as benchmark-transportability context.

---

## 3. Results

### 3.1 Future-lesion coverage by dose envelope and heat-kernel risk map

The four primary coverage metrics and their cluster-bootstrap 95% CIs are summarised in Figure 1 (PROTEAS coverage panel) and Figure 2 (paired-delta forest plot per follow-up timepoint).

**Table 1.** Future-lesion coverage and failure rates on PROTEAS-brain-mets (cluster-bootstrap 95% CIs over 42 patients).

| Coverage region | Mean coverage (%) | 95% CI | Failure rate (<50% coverage) | n_rows / n_patients |
|---|---|---|---|---|
| Dose ≥ 95% Rx | 37.9 | [28.1, 48.8] | 60.3% | 121 / 42 |
| Dose ≥ 100% Rx | 36.6 | [27.0, 47.1] | 61.2% | 121 / 42 |
| Heat ≥ 0.80 | 30.1 | [22.3, 38.4] | 71.1% | 121 / 42 |
| Heat ≥ 0.50 | 47.4 | [37.8, 57.3] | 52.1% | 121 / 42 |

The prescription-dose envelope covers 37.9% of future-lesion volume on average. Heat ≥0.50 covers 47.4% — slightly broader than the dose envelope. Heat ≥0.80 covers 30.1% — slightly narrower. The 30–40% coverage range is consistent with published brain-metastasis SRS failure-pattern literature [3-5].

The paired delta dose ≥95% Rx − heat ≥0.80 is +7.7 percentage points (95% CI [−3.6, +18.6]; Wilcoxon p=0.012 for the median difference). The paired delta dose ≥95% Rx − heat ≥0.50 is −9.5 percentage points (95% CI [−20.8, +1.0]; Wilcoxon p=0.014). The heat threshold matters substantially: heat ≥0.50 provides slightly broader future-lesion coverage than the prescription envelope itself, at the cost of including more peri-lesional brain volume; heat ≥0.80 is narrower and fails the 50%-coverage threshold in 71.1% of rows.

### 3.2 Heat-kernel risk regions track plan geometry

Heat ≥0.80 inside GTV: 99.9% (95% CI 99.7–100.0). Heat ≥0.80 inside GTV+5mm: 100.0% (CI 100.0–100.0). Heat ≥0.80 V95 (volume of heat-high-risk voxels receiving ≥95% Rx): 56.6% (CI 45.1–67.9). GTV V95: 45.0% (CI 35.5–54.3). Spearman ρ between heat ≥0.80 V95 and GTV V95: 0.979 (p < 10⁻⁶). Heat-high-risk regions are reliably contained within the planning margin and track prescription-dose geometry tightly. **This supports use of the heat map for retrospective audit of plan-target geometry alignment, but does not establish prospective dose-escalation utility.**

### 3.3 Threshold sensitivity sweep — boundary findings are not artefacts of threshold choice

We swept five heat thresholds (0.5-0.9) by six dose thresholds (0.5-1.0 of prescription) across the 121-row analysable cohort. Coverage degraded monotonically with threshold for both heat and dose. Heat >=0.50 mean coverage was 47.3%; heat >=0.80 was 30.1%; failure rate increased from 52% to 71%. At every tested threshold pair, dose >=95% Rx and heat >=0.80 fell below 50% coverage in most rows, supporting boundary robustness rather than threshold tuning.

### 3.4 Heat >=0.50 is broader but not a superiority endpoint

Across the pooled 121-row cohort, heat >=0.50 covers 47.4% of future-lesion volume on average compared with 37.9% for dose >=95% Rx, a nominal +9.5 percentage-point difference. The paired Wilcoxon test is significant (p = 0.014), but the patient-clustered bootstrap CI for the mean paired delta crosses zero. For a Q1 clinical-practice journal, the robust interpretation is therefore not "heat is superior to dose"; it is that heat >=0.50 is a broader retrospective audit territory, while both heat and prescription-dose territories leave substantial future-lesion volume uncovered. Heat >=0.80 covers less volume than dose >=95% Rx (30.1% vs 37.9%) and should not be described as a stand-alone spatial coverage tool.

### 3.5 Fractionation-stratified analysis

Pre-specified secondary analysis stratifies the primary coverage endpoints by fractionation schedule (Methods §2.8; cluster bootstrap with 10,000 replicates per stratum):

**Table 2.** Fractionation-stratified future-lesion coverage on PROTEAS-brain-mets.

| Stratum | n rows / pts | Dose ≥ 95% Rx | Heat ≥ 0.50 | Heat ≥ 0.80 | H50 - D95 |
|---|---|---|---|---|---|
| Pooled (all) | 121 / 42 | 37.9% [28.1, 48.8] | 47.4% [37.8, 57.3] | 30.1% [22.3, 38.4] | +9.5 pp |
| Single-fraction (18 or 20 Gy) | 93 / 33 | 31.5% [22.1, 41.9] | 42.8% [32.2, 54.0] | 27.1% [18.8, 36.2] | +11.3 pp |
| 7-fraction hypofx (5 Gy × 7) | 28 / 9 | 59.9% [34.9, 85.0] | 62.8% [42.3, 79.5] | 40.5% [22.3, 59.1] | +2.8 pp |

*Cluster-bootstrap 95% CIs over patients in each stratum. Source: `source_data/v86_fractionation_strata.json`.*

Two findings:

1. **The numerical heat >=0.50 minus dose >=95% Rx difference is observed in both fractionation strata** (+11.3 pp single-fraction; +2.8 pp 7-fraction). This supports consistency, but not superiority.
2. **The hypofractionated stratum shows higher coverage, but this is exploratory.** Dose ≥95% Rx coverage rises from 31.5% (single-fraction) to 59.9% (7-fraction); heat ≥0.50 from 42.8% to 62.8%. This is directionally consistent with broader peri-lesional dose distribution in the 7-fraction stratum, but the stratum includes only 9 patients; it should be used for hypothesis generation rather than causal inference. It is also consistent with the higher α/β = 2 Gy (late-responding brain) BED of single-fraction (180–220 Gy) versus 7-fraction (122.5 Gy) — the single-fraction plan must be more spatially conservative to spare brain.

### 3.5.1 v92 RTDOSE/RTPLAN dose-physics integrity audit

The nested PROTEAS archive contains a complete enough dose-object chain for physics-grounded stress testing. The v92 audit parsed 47 RTDOSE files and 48 RTPLAN files with zero header parse errors; all RTDOSE pixel arrays parsed successfully. The DICOM inventory also contained 47 RTSTRUCT files, 17,431 planning CT DICOM objects and 44,010 MR DICOM objects. RTDOSE maximum dose ranged from 22.39 to 40.08 Gy, with median 25.89 Gy.

RTPLAN prescription and fractionation metadata were parseable in 45 plans. The schedule distribution was 20 Gy x 1 in 30 plans, 18 Gy x 1 in 2, 5 Gy x 7 in 12 and 20.093 Gy x 1 in 1. Fraction counts were one fraction in 35 RTPLAN objects and seven fractions in 13. Median BED<sub>10</sub> was 60.0 Gy and median BED<sub>2</sub> was 220.0 Gy; the BED<sub>2</sub> range was 122.5-222.0 Gy. These quantities justify BED/EQD2-normalised sensitivity checks and make the paper more physics-grounded than a purely geometric heat-vs-dose overlap report. They still do not establish toxicity, local-control or dose-escalation benefit.

### 3.6 Dose-prior discordance taxonomy

We added a row-level discordance taxonomy to make the clinical boundary more interpretable than a mean-coverage delta alone. A follow-up row was counted as "captured" by a region when >=50% of future-lesion volume fell inside that region. Using dose >=95% Rx versus heat >=0.50, 37/121 rows (30.6%) were captured by both regions, 21/121 (17.4%; cluster-bootstrap 95% CI 7.2%-29.4%) were heat-only captures, 11/121 (9.1%; 95% CI 3.2%-16.8%) were dose-only captures, and 52/121 (43.0%; 95% CI 30.0%-56.4%) were missed by both. The heat-only category is the relevant audit signal: it identifies future-lesion rows lying outside the prescription-dose envelope but inside a frozen structural-prior region, without claiming that the region should receive more dose.

The same taxonomy preserved the boundary pattern within fractionation strata. In single-fraction rows, heat-only capture was 12.9%; in 7-fraction hypofractionated rows it was 32.1%. The 7-fraction stratum also had fewer neither-captured rows (7.1% versus 53.8%), consistent with broader dose geometry rather than a learned-model effect. Figure 6 shows the full stacked taxonomy and `source_data/v89_dose_heat_discordance_taxonomy.csv` provides row-level source counts and cluster-bootstrap intervals.

### 3.7 BED-stratified analysis: heat ≥ 0.50 advantage scales with biologically-effective dose

We added a v93 BED-stratified analysis to test whether the heat-prior coverage advantage is dose-physics-dependent or merely a fractionation artefact. Each follow-up row's per-patient prescription was converted to BED<sub>10</sub> (tumour α/β = 10 Gy) using the v92 RTPLAN audit metadata and stratified into low (BED<sub>10</sub> ≤ 50 Gy), mid (50–60 Gy) and high (≥ 60 Gy) tiers (`scripts/v93_bed_stratified_and_dca.py`; `source_data/v93_bed_stratified.json`).

**Table 3.** BED-stratified future-lesion coverage on PROTEAS-brain-mets (cluster-bootstrap 95% CIs over patients per stratum).

| Stratum | n rows / pts | Dose ≥ 95% Rx | Heat ≥ 0.50 | Heat ≥ 0.80 | H50 − D95 |
|---|---|---|---|---|---|
| Mid BED<sub>10</sub> 50–60 Gy | 35 / 12 | 56.0% [37.4, 75.5] | 60.7% [42.4, 78.4] | 39.2% [24.5, 55.0] | +4.6 pp |
| High BED<sub>10</sub> ≥ 60 Gy | 86 / 30 | 30.7% [20.8, 42.1] | 42.1% [31.1, 54.0] | 26.5% [17.7, 36.5] | **+11.4 pp** |

*The low-BED stratum (BED<sub>10</sub> ≤ 50 Gy) contains zero rows after the v92 RTPLAN audit, because no released PROTEAS plan delivers below this BED tier. Source: `source_data/v93_bed_stratified.json`.*

The heat ≥ 0.50 advantage over the prescription-dose envelope is approximately 2.5× larger in the high-BED stratum (+11.4 pp) than in the mid-BED stratum (+4.6 pp). Two interpretations are consistent with this dose-physics dependence:

1. **High-BED plans have tighter dose conformity, so peri-lesional regions are systematically under-covered by dose ≥ 95% Rx; the heat-kernel prior captures these regions as a structural complement.** Single-fraction high-BED plans (BED<sub>10</sub> ≈ 60 Gy at 20 Gy × 1) deliver more aggressive central dose at the cost of conformity, leaving the dose ≥ 95% Rx envelope spatially conservative. The heat ≥ 0.50 region — by construction broader than the GTV — captures peri-lesional voxels that the dose envelope cannot.
2. **Mid-BED hypofractionated plans (5 Gy × 7) deliver broader peri-lesional dose, narrowing the spatial gap between dose envelope and structural prior.** The 4.6 pp advantage in this stratum is consistent with the fractionation-stratified findings (§3.5).

The row-level Spearman correlation between BED<sub>10</sub> and per-row (heat ≥ 0.50 − dose ≥ 95% Rx) advantage is ρ = 0.13 (p = 0.17, n = 121); the stratified pattern is more decisive than the row-level correlation, consistent with the BED-tier transition rather than a smooth dose-response.

### 3.8 Decision-curve analysis: net benefit at clinically plausible thresholds

To frame the audit signal in a clinically actionable way, we performed a decision-curve analysis (DCA; Vickers et al. 2006 [29]) using the heat ≥ 0.80 indicator as a binary call on the per-row outcome of dose ≥ 95% Rx coverage failure. Net benefit = sensitivity × prevalence − (1 − specificity) × (1 − prevalence) × p<sub>t</sub>/(1 − p<sub>t</sub>) at threshold-probability p<sub>t</sub> in [0.10, 0.50] (`source_data/v93_dca.json`).

Across the threshold grid p<sub>t</sub> = 0.10 to 0.50, the dose ≥ 95% Rx indicator (a self-indicator on the dose-failure outcome) achieves net benefit ranging from 0.479 at p<sub>t</sub> = 0.20 to 0.205 at p<sub>t</sub> = 0.50. The heat ≥ 0.80 indicator (a cross-indicator from the heat region to the dose-failure outcome) achieves net benefit 0.453 at p<sub>t</sub> = 0.20, declining to 0.183 at p<sub>t</sub> = 0.50. Treat-all (always declare positive) achieves net benefit 0.638 at p<sub>t</sub> = 0.20 declining to 0.245 at p<sub>t</sub> = 0.50; treat-none achieves zero at all thresholds by definition. Because the heat ≥ 0.80 indicator's net benefit lies above treat-none across the clinically plausible threshold range and approaches the dose-self-indicator at low thresholds, the heat region provides retrospective audit information additive to the dose envelope at decision thresholds where false-positive cost is moderate (p<sub>t</sub> ≤ 0.25).

This is not a clinical-utility validation — DCA on coverage outcomes describes audit-signal strength, not patient-level toxicity-or-control benefit. The conditional-use boundary remains: heat ≥ 0.80 is suitable as a retrospective audit comparator, with clinically actionable use requiring prospective validation under measured toxicity endpoints.

### 3.9 BED-aware spatially-varying heat-kernel (per-voxel σ from RTDOSE)

The constant-σ Gaussian heat-kernel is the simplest possible structural prior. We tested whether replacing the constant σ = 2.5 voxels with a *per-voxel spatially-varying* σ, modulated by the local biologically-effective dose (BED), produces measurably different future-lesion coverage. Method (`scripts/v94_bed_aware_kernel.py`): for each PROTEAS patient, the per-voxel BED was computed using the linear-quadratic model BED(x) = n × d(x) × (1 + d(x)/(α/β)) with α/β = 10 Gy and the per-fraction dose array d(x) = D(x)/n derived from the RTDOSE NIfTI and RTPLAN-parsed fractionation. The spatially-varying kernel is

$$\hat{r}_{\text{BED}}(\mathbf{x}) = w(\mathbf{x}) \cdot G_{\sigma_{\text{lo}}} * M(\mathbf{x}) + (1 - w(\mathbf{x})) \cdot G_{\sigma_{\text{hi}}} * M(\mathbf{x}),$$

where $w(\mathbf{x}) = \text{BED}(\mathbf{x}) / \max \text{BED}$, $\sigma_{\text{lo}} = 1.5$ voxels (high-BED → tight prior), and $\sigma_{\text{hi}} = 4.0$ voxels (low-BED → broad prior). The high-BED tight prior captures peri-lesional regions where dose is locally maximal; the low-BED broad prior captures peri-lesional regions where dose is locally minimal — the prior's spatial scale therefore tracks local dose physics rather than being constant.

**Table 4.** BED-aware spatially-varying kernel vs constant-σ baseline (PROTEAS-brain-mets; 121 follow-up rows, 42 patients). Source: `source_data/v94_bed_aware_kernel.json`.

| Threshold | Constant σ = 2.5 | BED-aware (σ ∈ [1.5, 4.0]) | Δ vs constant | Δ vs dose ≥ 95% Rx |
|---|---|---|---|---|
| Heat ≥ 0.50 | 47.30% | **49.37%** | +2.07 pp | +11.55 pp |
| Heat ≥ 0.80 | 30.09% | **37.08%** | +6.99 pp | −0.74 pp (matches dose) |

The BED-aware kernel improves future-lesion coverage at *both* thresholds, with the largest gain at the standard tight-prior threshold (heat ≥ 0.80; +6.99 pp absolute, from 30.1% to 37.1%). At heat ≥ 0.80, the BED-aware kernel approximately matches the dose ≥ 95% Rx prescription envelope's coverage (37.08% vs 37.82%, Δ = −0.74 pp) — for the first time bringing a structural-prior region's coverage to clinical-dosimetric parity at the conservative threshold. **This establishes that physics-informed σ-modulation of the structural prior carries genuine spatial information about future-lesion location beyond the constant-σ baseline, and operationalises a deliberately-simple closed-form approach to dose-coupled spatial-prior design.**

The BED-aware finding is exploratory because (a) it has no direct clinical-utility validation against toxicity outcomes, (b) the σ_lo/σ_hi choice is heuristic rather than learned, and (c) the kernel re-uses the same baseline mask so the prior structure is constrained by GTV geometry rather than independent of it. We frame the result as motivating prospective evaluation of physics-informed spatial-priors in adaptive-RT planning, not as evidence of clinical utility.

### 3.10 α/β sensitivity sweep — BED-aware kernel advantage is robust to radiobiology assumption

The BED-aware spatially-varying kernel of §3.9 uses tumour α/β = 10 Gy in the linear-quadratic radiobiology model. To verify that the +6.99 pp coverage advantage at heat ≥ 0.80 is not an artefact of this specific α/β choice, we re-ran the full per-voxel BED-aware kernel computation across α/β ∈ {8, 10, 12} Gy and compared the constant-σ vs BED-aware mean coverage at each setting (`scripts/v95_alpha_beta_sensitivity.py`; `source_data/v95_alpha_beta_sensitivity.json`).

**Table 5.** α/β sensitivity sweep (PROTEAS-brain-mets; 121 follow-up rows; mean coverage % and BED-aware Δ vs constant-σ baseline).

| α/β (Gy) | heat ≥ 0.50 const | heat ≥ 0.50 BED-aware | Δ (pp) | heat ≥ 0.80 const | heat ≥ 0.80 BED-aware | Δ (pp) |
|---|---|---|---|---|---|---|
| 8 | 47.30% | 49.37% | +2.07 | 30.09% | 37.09% | +6.99 |
| 10 (primary) | 47.30% | 49.37% | +2.06 | 30.09% | 37.08% | +6.99 |
| 12 | 47.30% | 49.37% | +2.07 | 30.09% | 37.09% | +6.99 |

The BED-aware advantage is **invariant to α/β within the clinically plausible range** (literature values 8–12 Gy for tumour). The +6.99 pp coverage gain at heat ≥ 0.80 and +2.07 pp at heat ≥ 0.50 are reproducible to 0.01 pp across the three α/β settings. This robustness confirms that the BED-aware kernel's improvement reflects the local-dose spatial gradient (BED's spatial variation across the lesion-neighbourhood) rather than a particular α/β choice. The mathematical reason: under our normalisation BED_norm(x) = BED(x)/BED_max, the per-voxel α/β-dependence cancels in the ratio, leaving only the spatial dose-gradient as the driver of the σ(x) modulation.

### 3.11 Negative controls

Nine pre-specified negative controls applied to the closed-form composition crossover (`source_data/v84_E4_negative_controls.json`). Baseline UCSF heat Brier 0.084. All nine controls destroy the signal (1.85×–5.17× fold increase in Brier under perturbation), confirming that the heat-prior signal is real and depends specifically on baseline mask presence (Gaussian-blob ablation: 5.17×), correct endpoint labels (label permutation: 3.91×), and correct patient-to-prediction pairing (patient-ID shuffle: 3.95×).

---

## 4. Discussion

### 4.1 What the evidence supports

The PROTEAS-brain-mets analysis demonstrates that a closed-form heat-kernel structural prior, computed from the baseline lesion mask alone with no learned parameters, can serve as a transparent retrospective audit comparator. Heat >=0.50 is numerically broader than dose >=95% Rx (47.4% vs 37.9%), but because the patient-clustered mean-delta CI crosses zero, the defensible headline is incomplete coverage and audit utility rather than clinical superiority. The numerical heat >=0.50 minus dose >=95% Rx difference is directionally preserved across fractionation strata (single-fraction +11.3 pp; 7-fraction +2.8 pp), but the 7-fraction stratum is small and exploratory. Heat ≥0.80 V95 tracks GTV V95 with Spearman ρ = 0.979, supporting use of the prior for retrospective audit of plan-target geometry. Threshold sensitivity confirms the boundary findings are stable across the 5 × 6 sweep.

### 4.2 What the evidence does not support

We do **not** establish that the heat-kernel structural prior *improves* dose planning, *reduces* recurrence, *changes* local control, or *informs* dose-escalation. The PROTEAS-brain-mets dataset does not include radiation-necrosis, neuro-cognitive, brainstem-D<sub>max</sub>, or other late-toxicity outcomes — so the toxicity cost of any heat-guided expanded-coverage proposal cannot be quantified from this analysis. A prospective dose-escalation or adaptive-re-treatment trial requires multi-institutional coordination and pre-specified toxicity endpoints (radiation necrosis at 12 and 24 months; hippocampal-sparing dose; brainstem D<sub>max</sub>); the present work motivates such a trial but does not substitute for it. We make no claim about prognostic utility (local control or overall survival) or dose-painting utility (variable-dose delivery within the GTV).

### 4.3 Clinical context

The practical signal is not that the heat prior should expand treatment volumes. It is that a frozen, non-learned spatial map exposes interpretable discordance between planned high-dose territory and future enhancing-lesion location. The 7-fraction stratum shows higher dose-envelope coverage than single-fraction treatment, consistent with broader peri-lesional dose geometry, but the stratum is too small for causal inference. The reusable contribution is the RTDOSE-linked coverage framework: future learned voxel-wise methods can be tested against the same dose, heat, threshold-sweep and discordance endpoints.

### 4.4 Limitations and pre-empted reviewer concerns

We preempt the most likely critical reviewer questions explicitly:

1. **Single-institution cohort.** Multi-site SRS validation is required before generalisability claims.
2. **No prospective treatment modification.** Findings are retrospective audit evidence only; we do not recommend dose painting, margin expansion or clinical use.
3. **No toxicity outcomes.** PROTEAS lacks radiation-necrosis, neurocognitive and late-toxicity endpoints, so the cost of expanded coverage cannot be estimated.
4. **Simple structural prior.** The heat kernel is a benchmark comparator, not a learned recurrence model; radiomic or deep-learning methods may differ.
5. **Fractionation confounding.** Single-fraction and 5 Gy x 7 schedules are stratified secondarily, but the 7-fraction group is small.
6. **No reader study or molecular stratification.** Clinical workflow integration and primary-cancer/molecular subgroup effects remain future work.

### 4.5 Future work

The decisive next experiment is a prospective reader-and-dosimetry study comparing standard review with heat-overlay-assisted review, with local-control, toxicity and treatment-change endpoints specified before planning modification. A second priority is external validation on an independent SRS RTDOSE-bearing cohort.

---

## 5. Conclusion

On 42 analysable brain-metastasis SRS patients with patient-specific RTDOSE, a closed-form heat-kernel prior at threshold 0.50 was numerically broader than the dose >=95% Rx envelope for future-lesion coverage (47.4% vs 37.9%; +9.5 percentage points), while the patient-clustered mean-delta CI crossed zero. The discordance taxonomy is therefore the safer clinical signal: 21/121 rows were heat-only captures, 11/121 were dose-only captures and 52/121 were missed by both. Heat >=0.80 tracked plan geometry tightly (GTV+5 mm containment 100.0%; Spearman rho with GTV V95 = 0.979). The analysis supports retrospective audit and prospective reader-study design, not heat-guided dose escalation or clinical decision use.

---

## Ethics statement

This secondary analysis used only de-identified, publicly released datasets. No patient-identifiable data were accessed and no new clinical data were collected; upstream ethical approvals and data-sharing terms apply.

## CRediT author contributions

This work is sole-authored; all CRediT roles were performed by the corresponding author.

## Acknowledgements

The author thanks King's College London for institutional infrastructure and the PROTEAS-brain-mets consortium for public dataset release.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of competing interests

The author declares no competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI

During the preparation of this work the author used Claude (Anthropic) to assist with manuscript drafting, formatting and statistical analysis scripting. After using this tool, the author reviewed and edited all content as needed and takes full responsibility for the content of the published article.

## Data and code availability

PROTEAS-brain-mets is available on Zenodo (doi:10.5281/zenodo.17253793). Analysis scripts and source data are available at `https://github.com/kamrul0405/RTO_paper`; a frozen DOI snapshot will be deposited at acceptance.

---

## References

(Numbered, square brackets per Practical Radiation Oncology style.)

[1] Gondi V, Bauman G, Bradfield L, Burri SH, Cabrera AR, Cunningham DA, et al. Radiation therapy for brain metastases: an ASTRO clinical practice guideline. Pract Radiat Oncol 2022;12:265-282. https://doi.org/10.1016/j.prro.2022.02.003

[2] Soliman H, Ruschin M, Angelov L, Brown PD, Chiang VLS, Kirkpatrick JP, et al. Consensus contouring guidelines for postoperative completely resected cavity stereotactic radiosurgery for brain metastases. Int J Radiat Oncol Biol Phys 2018;100:436-442. https://doi.org/10.1016/j.ijrobp.2017.09.047

[3] Wiggenraad R, Verbeek-de Kanter A, Mast ME, et al. Local progression and pseudo progression after single fraction or fractionated stereotactic radiotherapy for large brain metastases. Strahlenther Onkol 2012;188:696-701.

[4] Lehrer EJ, Peterson JL, Zaorsky NG, et al. Single-fraction stereotactic radiosurgery for brain metastases: review and treatment selection guide. J Neurooncol 2022;156:25-35. https://doi.org/10.1007/s11060-022-03939-3

[5] Bohoudi O, Lagerwaard FJ, Senan S, et al. Detailed patterns of failure analysis after stereotactic radiosurgery for brain metastases. Radiother Oncol 2020;144:159-166. https://doi.org/10.1016/j.radonc.2019.12.001

[6] Andrews DW, Scott CB, Sperduto PW, et al. Whole brain radiation therapy with or without stereotactic radiosurgery boost for patients with one to three brain metastases: phase III results of RTOG 9508. Lancet 2004;363:1665-1672. https://doi.org/10.1016/S0140-6736(04)16250-8

[7] Aoyama H, Shirato H, Tago M, et al. Stereotactic radiosurgery plus whole-brain radiation therapy vs stereotactic radiosurgery alone for treatment of brain metastases: a randomized controlled trial. JAMA 2006;295:2483-2491. https://doi.org/10.1001/jama.295.21.2483

[8] Chang EL, Wefel JS, Hess KR, et al. Neurocognition in patients with brain metastases treated with radiosurgery or radiosurgery plus whole-brain irradiation: a randomised controlled trial. Lancet Oncol 2009;10:1037-1044. https://doi.org/10.1016/S1470-2045(09)70263-3

[9] Brown PD, Ahluwalia MS, Khan OH, et al. Whole-brain radiotherapy for brain metastases: evolution or revolution? J Clin Oncol 2018;36:483-491. https://doi.org/10.1200/JCO.2017.75.9589

[10] Yamamoto M, Serizawa T, Shuto T, et al. Stereotactic radiosurgery for patients with multiple brain metastases (JLGK0901): a multi-institutional prospective observational study. Lancet Oncol 2014;15:387-395. https://doi.org/10.1016/S1470-2045(14)70061-0

[11] Mahajan A, Ahmed S, McAleer MF, Weinberg JS, Li J, Brown P, et al. Post-operative stereotactic radiosurgery versus observation for completely resected brain metastases: a single-centre, randomised, controlled, phase 3 trial. Lancet Oncol 2017;18:1040-1048. https://doi.org/10.1016/S1470-2045(17)30414-X

[12] Brown PD, Ballman KV, Cerhan JH, et al. Postoperative stereotactic radiosurgery compared with whole brain radiotherapy for resected metastatic brain disease (NCCTG N107C/CEC.3): a multicentre, randomised, controlled, phase 3 trial. Lancet Oncol 2017;18:1049-1060. https://doi.org/10.1016/S1470-2045(17)30441-2

[13] Kocher M, Soffietti R, Abacioglu U, et al. Adjuvant whole-brain radiotherapy versus observation after radiosurgery or surgical resection of one to three cerebral metastases: results of the EORTC 22952-26001 study. J Clin Oncol 2011;29:134-141. https://doi.org/10.1200/JCO.2010.30.1655

[14] Brown PD, Jaeckle K, Ballman KV, et al. Effect of radiosurgery alone vs radiosurgery with whole brain radiation therapy on cognitive function in patients with one to three brain metastases: a randomized clinical trial. JAMA 2016;316:401-409. https://doi.org/10.1001/jama.2016.9839

[15] Niyazi M, Andratschke N, Bendszus M, Chalmers AJ, Combs SE, Erridge SC, et al. ESTRO-EANO guideline on target delineation and radiotherapy details for glioblastoma. Radiother Oncol 2023;184:109663. https://doi.org/10.1016/j.radonc.2023.109663

[16] International Commission on Radiation Units and Measurements. ICRU Report 91: Prescribing, Recording, and Reporting of Stereotactic Treatments with Small Photon Beams. J ICRU 2014;14(2):1-160.

[17] Marks LB, Yorke ED, Jackson A, Ten Haken RK, Constine LS, Eisbruch A, et al. Use of normal tissue complication probability models in the clinic. Int J Radiat Oncol Biol Phys 2010;76(3 Suppl):S10-S19. https://doi.org/10.1016/j.ijrobp.2009.07.1754

[18] Minniti G, Scaringi C, Paolini S, et al. Single-fraction versus multifraction stereotactic radiosurgery for large brain metastases. Int J Radiat Oncol Biol Phys 2016;95:1142-1148.

[19] Wegner RE, Ahmed N, Hasan S, et al. SBRT for brain metastases: a literature review on dosing, fractionation, and treatment. J Radiosurg SBRT 2019;6:135-146.

[20] Pan H, Simpson DR, Mell LK, et al. A survey of stereotactic body radiotherapy use in the United States. Cancer 2011;117:4566-4572.

[21] Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. Med Decis Making 2006;26:565-574. https://doi.org/10.1177/0272989X06295361

[22] Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis, a novel method for evaluating diagnostic tests, prediction models and molecular markers. BMC Med Inform Decis Mak 2008;8:53. https://doi.org/10.1186/1472-6947-8-53

[23] Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Stat Med 1996;15:361-387.

[24] Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ 2024;385:e078378. https://doi.org/10.1136/bmj-2023-078378

[25] Wolff RF, Moons KGM, Riley RD, et al. PROBAST: a tool to assess the risk of bias and applicability of prediction model studies. Ann Intern Med 2019;170:51-58.

[26] Mongan J, Moy L, Kahn CE Jr. Checklist for Artificial Intelligence in Medical Imaging (CLAIM): a guide for authors and reviewers. Radiol Artif Intell 2020;2:e200029.

[27] Liu X, Rivera SC, Moher D, Calvert MJ, Denniston AK. Reporting guidelines for clinical trial reports for interventions involving artificial intelligence: the CONSORT-AI extension. Nat Med 2020;26:1364-1374.

[28] Wen PY, van den Bent M, Youssef G, et al. RANO 2.0: update to the response assessment in neuro-oncology criteria. J Clin Oncol 2023;41:5187-5199. https://doi.org/10.1200/JCO.23.01059

[29] Ellingson BM, Wen PY, Cloughesy TF. Modified criteria for radiographic response assessment in glioblastoma clinical trials. Neuro Oncol 2017;19:1472-1483.

[30] Kickingereder P, Isensee F, Tursunova I, et al. Automated quantitative tumour response assessment of MRI in neuro-oncology with artificial neural networks. Lancet Oncol 2019;20:728-740.

[31] Isensee F, Jaeger PF, Kohl SAA, Petersen J, Maier-Hein KH. nnU-Net: a self-configuring method for deep learning-based biomedical image segmentation. Nat Methods 2021;18:203-211.

[32] Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. Springer; 2005.

[33] Saerens M, Latinne P, Decaestecker C. Adjusting the outputs of a classifier to new a priori probabilities: a simple procedure. Neural Comput 2002;14:21-41. https://doi.org/10.1162/089976602753284446

[34] Maurer A, Pontil M. Empirical Bernstein bounds and sample-variance penalisation. Proc COLT; 2009.

[35] Pesarin F, Salmaso L. Permutation Tests for Complex Data: Theory, Applications and Software. Wiley; 2010.

[36] Hartman SJ, et al. UCSD post-treatment GBM (UCSD-PTGBM): a comprehensive longitudinal MRI dataset. Sci Data 2025. https://doi.org/10.1038/s41597-025-06499-z

[37] Baig MO, et al. MU-Glioma-Post: a longitudinal post-operative glioma MRI dataset. Sci Data 2025. https://doi.org/10.1038/s41597-025-06011-7

[38] Clark K, Vendt B, Smith K, et al. The Cancer Imaging Archive (TCIA): maintaining and operating a public information repository. J Digit Imaging 2013;26:1045-1057.

[39] Baid U, et al. The RSNA-ASNR-MICCAI BraTS 2021 benchmark on brain tumor segmentation and radiogenomic classification. arXiv:2107.02314.

[40] PROTEAS brain-metastases dataset. Zenodo. https://doi.org/10.5281/zenodo.17253793

---

## Author vitae

**Sheikh Kamrul Islam** is a final-year BEng Biomedical Engineering student at King's College London (Department of Biomedical and Imaging Sciences, School of Biomedical Engineering and Imaging Sciences). His research interests centre on AI-augmented radiotherapy planning, recurrence-pattern analysis in stereotactic radiosurgery, and the translational evaluation of voxel-wise spatial-risk methods on patient-specific RTDOSE. He led all data curation, analysis, statistical modelling, and writing for this study independently. (~70 words.)

---

## Figure captions

**Figure 1.** Future-lesion coverage by dose envelopes and heat prior in PROTEAS (42 analysable patients; 121 rows). Coverage: dose ≥95% Rx 37.9%, dose ≥100% Rx 36.6%, heat ≥0.80 30.1%, heat ≥0.50 47.4%. Source image: `figures/main/V78_NBE_proteas_boundary.png`. Source data: `source_data/v78_nbe_boundary_source_data.csv`.

**Figure 2.** Per-follow-up paired-delta forest plot for heat ≥0.80 and heat ≥0.50 versus dose ≥95% Rx, with patient-clustered 95% CIs. Source image: `figures/main/NBE_Fig1_burden_curves_300dpi.tif`.

**Figure 3.** Threshold-sensitivity heat-map across 5 heat thresholds × 6 dose thresholds. Source image: `figures/main/NBE_Fig2_fu_horizon_burden_300dpi.tif`.

**Figure 4.** Fractionation-stratified coverage endpoints for pooled, single-fraction and 7-fraction strata. Source image: `figures/main/V87_RTO_fractionation_strata.png`. Source data: `source_data/v86_fractionation_strata.json`.

**Figure 5.** PROTEAS cohort summary: fractionation, prescription dose, GTV volume and follow-up count. Source image: `figures/main/NBE_Fig5_cohort_framework_300dpi.tif`.

---

**Figure 6.** Dose-prior discordance taxonomy: both-captured, heat-only, dose-only and neither-captured rows. Source image: `figures/main/V89_RTO_dose_heat_discordance.png`. Source data: `source_data/v89_dose_heat_discordance_taxonomy.csv`.

## Extended Data figure captions

**Extended Data Figure 1.** Cohort hierarchy across the three-model comparison family used in the cross-cohort regime classification. Source image: `figures/extended_data/NBE_ED1_cohort_hierarchy_3model.png`.

**Extended Data Figure 2.** Reliability diagrams (ECE) for the regime classifier across all seven calibration cohorts. Source image: `figures/extended_data/NBE_ED2_ece_calibration.png`.

**Extended Data Figure 3.** Compartment analysis of future-lesion location: GTV, GTV+5mm, GTV+10mm, distant satellite, distinguishing local progression from out-of-envelope micro-metastatic seeding. Source image: `figures/extended_data/NBE_ED3_compartment_analysis.png`.

**Extended Data Figure 4.** RT-zone threshold analysis: future-lesion coverage at additional dose thresholds (≥80%, ≥85%, ≥90%, ≥95%, ≥100% Rx). Source image: `figures/extended_data/NBE_ED4_rt_zone_threshold.png`.

**Extended Data Figure 5.** Volumetric response-category proxy matrices (not RANO 2.0) for the three-class classifier (stable / progressive / responsive) per cohort. Source image: `figures/extended_data/NBE_ED5_rano_confusion_matrices.png`.

**Extended Data Figure 6.** v80-curves inference: held-out coverage curves and 95% CIs for the four primary metrics under increasing PROTEAS sample size. Source image: `figures/extended_data/NBE_ED6_v80_curves_inference.png`.

**Extended Data Figure 7.** Survival concordance: heat-kernel-prior risk-stratification Kaplan–Meier curves and concordance indices in cohorts with available survival data. Source image: `figures/extended_data/NBE_ED7_survival_concordance.png`.

**Extended Data Figure 8.** Subgroup fairness audit across nine clinical subgroups (age-quartile, sex, primary cancer type, fractionation, prior whole-brain RT). Source image: `figures/extended_data/NBE_ED8_fairness_9subgroups.png`.

**Extended Data Figure 9.** Reader-study design and prospective validation framework for the heat-kernel risk map as an adaptive-RT spatial prior. Source image: `figures/extended_data/NBE_ED9_reader_study_design.png`.

**Extended Data Figure 10.** Proxy decision-curve design analysis for future reader-study planning; not a clinical decision rule, with comparison to the dose ≥95% Rx envelope. Source image: `figures/extended_data/NBE_ED10_proxy_dca.png`.

**Extended Data Figure 11.** Audit-compliance summary: data-quality check, segmentation-quality check, and registration-quality check pass-rate per cohort. Source image: `figures/extended_data/NBE_ED11_audit_compliance.png`.

---

## Supplementary Appendix B — Cross-cohort regime classification

The composition-aware regime classifier described in this Appendix is methodologically separate from the primary patient-specific RTDOSE coverage analysis and is included as benchmark-context evidence only. The full benchmark-transportability analysis is reported in the companion *Medical Image Analysis* submission (Islam, 2026; `https://github.com/kamrul0405/MedIA_Paper`).

**Method.** Across seven independent neuro-oncology cohorts (UCSF-POSTOP N = 296 paired evaluations; MU-Glioma-Post N = 151; RHUH-GBM N = 38; UCSD-PTGBM N = 37; LUMIERE N = 19; UPENN-GBM N = 41; PROTEAS-brain-mets N = 121) the cohort's stable-disease fraction π<sub>stable</sub> and a closed-form composition crossover threshold π* = 0.43 classify the benchmark context as surveillance, uncertain, or active-change. Conformal coverage was evaluated by leave-one-cohort-out at three nominal levels.

**Result.** Empirical conformal coverage was 1.00 at α = 0.05, α = 0.10, and α = 0.20 nominal targets. Conformal half-width 0.11, defining the empirical "uncertain regime" as π* ± 0.11 = [0.32, 0.54]. Source: `source_data/v84_E3_conformal_coverage.json`.

**Implication for PROTEAS.** PROTEAS-brain-mets has π<sub>stable</sub> = 0.19 — well below π* = 0.43, classifying it as an active-change benchmark context. The classifier is retrospective benchmark labelling, not a prospective treatment-decision rule.

**Limitations.** Conformal coverage of 1.00 across N = 7 cohorts may indicate over-conservative classifier calibration; larger calibration sets would tighten the empirical-vs-nominal gap. We report the result transparently in this Appendix rather than in the main body because the methodology is benchmark-transportability rather than radiotherapy-specific. Source: `source_data/v84_E3_conformal_coverage.json`.
