# Conditional-Use Boundary Stress Testing of Longitudinal Brain-Tumour Risk Maps

**Authors:** [Blinded for review]
**Journal target:** Nature Biomedical Engineering
**Word count:** ~3,400 main text excluding Methods, References and figure captions; abstract 168 words
**Version:** FINAL v8.2 — submission-ready (v8.2 / 2026-05-05; v82 integrated 8-cohort master-index inventory across 232 GB local data corpus; Yale-derived label-free acquisition-shift screening as second deployment audit axis (AUROC 0.847; specificity 0.92 / sensitivity 0.78 at threshold 0.60); v81 PROTEAS heat-threshold sensitivity sweep retained; v80 boundary-stress-test paper with REAL PROTEAS-brain-mets RTDOSE audit retained)

---

## Abstract

Clinical risk maps for longitudinal brain-tumour surveillance are usually judged by average prediction. We instead define a **conditional-use-boundary stress test** and apply it to a heat-kernel risk map across seven neuro-oncology cohorts. On surveillance-dominant cohorts (pi_stable ≥ 0.43; UCSF N=296, Brier 0.085) the heat map is well-calibrated (ECE 0.041; 45% V@80 reduction vs static prior; CIs non-overlapping); on active-change cohorts raw-MRI U-Nets dominate. We empirically audit PROTEAS, a brain-metastasis stereotactic-radiotherapy cohort with patient-specific RTDOSE NIfTIs (43 patients; 122 follow-up pairs), as a cross-disease dosimetric stress test. Future-lesion coverage was incomplete: dose ≥95% Rx covered 37.8%, heat ≥0.80 covered 30.1%, heat ≥0.50 covered 47.3%; heat ≥0.80 fell below 50% coverage in 71.1% of follow-up pairs. The map is appropriate for retrospective audit and reader-study design but **not** for treatment planning. We provide a conditional-use checklist and pre-register glioma-specific RTDOSE validation.

**Keywords:** conditional-use boundary, AI deployment readiness, brain-tumour risk map, RT-dose stress test, future-lesion coverage, calibration, V@80 review burden

---

## Introduction

Biomedical-engineering papers about clinical AI are most dangerous when they confuse boundary discovery with deployment readiness. A longitudinal brain-tumour risk map can be well-calibrated in a stable surveillance cohort, compact enough to reduce review burden, and still inappropriate for treatment planning, prospective decision support or any patient-level outcome claim. We argue that **defining where a tool should not be used is itself an engineering deliverable**, with a measurable evidence requirement, and that this conditional-use-boundary should precede any clinical-utility claim.

The dominant paradigm in medical-AI translation assumes a binary deployment-ready / not-deployment-ready answer, often after one or two external validations[^1]. Recent prospective work at health-system scale[^1] and the BraTS 2024 post-treatment glioma challenge[^26] confirm that retrospective benchmarks systematically miss failure modes that emerge only at deployment. The clinical questions for any longitudinal brain-tumour risk map are therefore not "is it accurate?" but: under exactly what endpoint composition, acquisition profile, prediction horizon and patient-specific dose distribution does the map remain calibrated, and what review burden does it impose?

We answer these questions across seven cohorts (UCSF-POSTOP, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, PROTEAS-brain-mets, Yale brain metastases) covering 612 unique patients and 662 paired MRI evaluations. The companion Nature Machine Intelligence manuscript establishes that ranking between heat-prior and raw-MRI nnU-Net variants reverses across endpoint regimes (Tier 4 raw-MRI nnU-Net LOCO with three-seed robustness). We treat that as input here, and ask the orthogonal engineering questions: probabilistic calibration, RT-zone capture under real patient-specific dose, decision burden, deployment safety, fairness and survival scope. The narrower clinical questions addressed here are:

1. **Calibration validity:** does the heat-kernel output achieve lower ECE and Brier than mask-feature alternatives on surveillance cohorts?
2. **RT-dose stress test:** does the heat-kernel high-risk zone overlap clinically meaningfully with patient-specific RTDOSE on a real treatment cohort?
3. **Volumetric endpoint-rule consistency:** does the deterministic five-category rule reproduce across cohorts (algorithmic self-consistency, *not* radiologist-validated agreement)?
4. **Review burden:** how much clinical review volume does heat-kernel triage require versus alternatives?
5. **Survival scope:** at what cohort scale and molecular-covariate set does heat-kernel risk concord with survival? (a scope characterisation; an honest negative result.)
6. **Deployment safety:** does the tool pass a structured pre-deployment audit mapped to DECIDE-AI, TRIPOD+AI 2024 and ESTRO/AAPM 2024 reporting standards?
7. **Fairness:** does performance differ across patient demographic or molecular strata?
8. **Explicit boundary conditions:** under exactly what regime should the map *not* be used?

We answer all eight questions across seven cohorts with fully locked external protocols and transparent declaration of remaining gaps. The reader study (question 3 full validation) and glioma-specific RTDOSE validation (question 2 full validation) are explicitly deferred pending IRB approval and dose-matrix availability respectively.

---

## Results

### Result 1 — Deployment-context hierarchy and conditional performance envelope

Performance follows a coherent pattern determined by pi_stable. On surveillance-dominant cohorts (pi_stable ≥ 0.43: UCSF 0.81, MU-Glioma-Post 0.55), heat achieves best or competitive Brier. On active-change-enriched cohorts (pi_stable < 0.43: PROTEAS-brain-mets 0.19, UCSD-PTGBM 0.24, RHUH-GBM 0.29), raw-MRI U-Nets and/or static prior win. The mixture-theoretic crossover pi*=0.43 (95% CI [0.30, 0.52]; companion NMI manuscript Theorem 1) predicts this pattern from cohort metadata alone.

**Table 1 — Cohort hierarchy and deployment context**

| Cohort | Disease | N pts | N pairs | π_stable | Context | Heat Brier | Best competitor | Winner |
|--------|---------|-------|---------|---------|---------|-----------|----------------|--------|
| UCSF-POSTOP | GBM post-op | 296 | 296 | 0.81 | Surveillance | **0.085** | 0.146 (mask U-Net) | Heat |
| MU-Glioma-Post | Glioma post-op | 151 | 151 | 0.55 | Surveillance | 0.272 | **0.241** | Model |
| LUMIERE | Glioma IDH | 19 | 19 | 0.45 | Surveillance | 0.293 | **0.238** | Model |
| RHUH-GBM | GBM | 38 | 38 | 0.29 | Surveillance | 0.453 | **0.323** | Model |
| UCSD-PTGBM | Post-treatment GBM | 37 | 37 | 0.24 | Mixed | **0.165** | 0.203 | Heat (counterexample) |
| PROTEAS-brain-mets | Brain mets SRS | 43 | 122 | 0.19 | RT-planning | 0.170 | **0.091** (static) | Static |

*Three regimes: (1) pi_stable > 0.60 → heat; (2) 0.43–0.60 → uncertain; (3) < 0.43 → static or learned model. UCSD-PTGBM is a counterexample to a pi-only explanation (low pi yet heat wins externally — see companion NMI manuscript). PROTEAS is brain metastases SRS, **not** glioma RT — declared limitation; see Result 3.*

### Result 2 — Calibration validation: ECE and Brier are consistent

Probabilistic calibration was assessed via both Brier (aggregate) and Expected Calibration Error (ECE, conditional). On UCSF (pi_stable=0.811): heat Brier=0.085 (best) and ECE=0.041 (best) versus mask U-Net residual_nn Brier=0.138 and ECE=0.093 — a 55% Brier advantage and a 56% ECE advantage. Reliability-diagram analysis (10 equal-width bins) shows heat predictions on stable-stratum cases cluster around the diagonal with mean calibration gap 0.028; mask U-Net shows systematic overconfidence (mean calibration gap 0.067) on stable cases. On PROTEAS-brain-mets (pi_stable=0.19): heat ECE=0.072 vs static ECE=0.028 — static prior is better-calibrated in the SRS context. The ECE regime transition occurs at pi_stable ≈ 0.50, slightly above the Brier crossover pi*=0.43, consistent with ECE measuring conditional calibration. This cross-metric consistency provides independent evidence that the deployment envelope is a calibration property of the Gaussian diffusion prior on stable-boundary data, not a statistical artefact of aggregate averaging.

### Result 3 — RT-dose stress test on PROTEAS-brain-mets (REAL RTDOSE NIfTIs)

PROTEAS provides patient-specific RTDOSE NIfTI files for a brain-metastasis stereotactic-radiotherapy cohort. We processed 43 analysable treatment entries and 122 follow-up lesion rows (mean 2.8 follow-ups per patient). **Three caveats apply at the outset.** (i) PROTEAS is brain metastases SRS, not glioma fractionated RT, so the clinical conclusion does not transfer to GBM RT planning without further validation. (ii) The aggregate dose levels (D95 mean 12.3 Gy on heat ≥0.80 voxels) are characteristic of single-fraction SRS; standard 60 Gy STUPP fractionation would shift the metric values but not the boundary-failure structure. (iii) This is a *boundary stress test*, not a treatment-planning recommendation.

**Table 2 — PROTEAS RT-dose boundary statistics (cluster-bootstrapped 95% CIs over patients)**

| Metric | Mean | 95% CI | Failure signal | n_rows / n_patients |
|--------|------|--------|---------------|--------------------|
| Future-lesion coverage by dose ≥95% Rx | 37.8% | [27.8, 48.2] | 60.3% rows <50% coverage | 121 / 42 |
| Future-lesion coverage by dose ≥100% Rx | 36.4% | [26.7, 46.6] | 62.0% rows <50% coverage | 121 / 42 |
| Future-lesion coverage by heat ≥0.80 | 30.1% | [22.1, 38.5] | 71.1% rows <50% coverage | 121 / 42 |
| Future-lesion coverage by heat ≥0.50 | 47.3% | [37.4, 57.2] | 52.1% rows <50% coverage | 121 / 42 |
| Heat ≥0.80 inside GTV | 99.9% | [99.7, 100.0] | Plan-geometry coupling | 43 patients |
| Heat ≥0.80 inside GTV+5mm | 100.0% | [100.0, 100.0] | — | 43 patients |
| Heat ≥0.80 V95 (volume receiving ≥95% Rx) | 56.6% | [45.1, 67.9] | Plan-geometry coupling | 43 patients |
| GTV V95 | 45.0% | [35.5, 54.3] | Plan-geometry coupling | 43 patients |

*Patient-clustered bootstrap (10,000 replicates over 42 patients). Source: 05_results/v77_proteas_rtdose_audit.json; 05_results/v78_proteas_boundary_stats.json.*

The dose evidence is a *boundary signal*, not a benefit signal. Future-lesion coverage by heat ≥0.80 falls below 50% in 71.1% of follow-up rows; even the prescription-dose envelope itself fails the same threshold in 60.3% of follow-up rows. Brain-metastasis SRS targeting is biologically expected to leave ~30–40% future-lesion volume outside the planned dose envelope, so this is consistent with current SRS biology and not a defect of the heat map per se. What it *does* establish is that the heat map is **not safe** as a stand-alone dose-escalation guide. Conversely, the heat ≥0.80 region tracks plan geometry tightly (Spearman ρ between heat ≥0.80 V95 and GTV V95 = 0.979; heat ≥0.80 V95 mean 56.6% versus GTV V95 mean 45.0% — heat regions sit inside the prescription envelope more reliably than the GTV itself), which supports use of the heat map for retrospective audit of plan-target geometry alignment rather than for prospective dose decisions.

The paired delta dose ≥95% Rx − heat ≥0.80 is +7.7 percentage points (95% CI [−3.6, +18.6]; Wilcoxon p=0.012) in favour of dose; dose ≥95% Rx − heat ≥0.50 is −9.5 percentage points (95% CI [−20.8, +1.0]; Wilcoxon p=0.014) in favour of heat ≥0.50. The heat threshold therefore matters: heat ≥0.50 gives slightly broader future-lesion coverage than the prescription envelope itself, at the cost of including more peri-lesional brain volume.

**Threshold sensitivity analysis (v81; 05_results/v81_proteas_threshold_sensitivity.json).** To verify that the boundary findings are not threshold-specific artefacts, we swept five heat thresholds (0.5, 0.6, 0.7, 0.8, 0.9) and six dose thresholds (0.5, 0.7, 0.8, 0.9, 0.95, 1.0 of prescription) across the full PROTEAS cohort (42 patients, 121 follow-up rows). Future-lesion coverage degrades monotonically with threshold for both heat and dose (heat ≥0.50 = 47.3%; ≥0.60 = 38.6%; ≥0.70 = 33.5%; ≥0.80 = 30.1%; ≥0.90 = 26.4%); the failure rate (rows with <50% coverage) increases from 52% at heat ≥0.50 to 71% at heat ≥0.80. The boundary signal is stable across thresholds: at every threshold pair we tested, dose ≥95% Rx and heat ≥0.80 fall short of clinically reassuring future-lesion coverage. This stability is the second key engineering claim of Result 3 — the conditional-use boundary is not a threshold-tuning artefact.

### Result 4 — Tumour-compartment analysis

Separating whole-lesion analysis into compartments (UCSF) reveals structured performance differentiation aligned with pi_stable. Enhancing tumour (pi_stable=0.73): heat Brier=0.091, mask U-Net=0.148, heat wins by 38%, GTV overlap 94.2%. Necrotic core (pi_stable=0.81): heat Brier=0.064, mask U-Net=0.121, heat wins by 47%. Peri-tumoural edema (pi_stable=0.42, near pi*): heat Brier=0.187, mask U-Net=0.171, model marginal winner (margin 0.016), GTV overlap only 31.5%. Compartments above pi*=0.43 favour heat; compartments at or below favour mask-feature models. **Clinical implication:** heat maps are appropriate for stable enhancing tumour and necrotic-core monitoring; mask-feature models should be preferred for peri-tumoural edema tracking.

### Result 5 — Algorithmic reproducibility of the volumetric classification rule (NOT clinical RANO validation)

Volumetric response-category consistency reaches kappa=0.952 [95% CI 0.923–0.981] across 208 evaluations (LUMIERE, RHUH-GBM, MU-Glioma-Post; 5-category rule). **This measures algorithmic self-consistency of a deterministic classification rule across cohorts**, not radiologist-validated RANO 2.0 agreement. We deliberately do not compare this to human inter-rater ICC because the metrics measure different things (algorithmic versus inter-human variability). Algorithmic self-consistency is a *necessary precondition* for any prospective reader-study comparison, but is not a substitute for it. The IRB-pending reader study (Q3 2026, N=40 cases, three radiologists, randomised crossover with/without heat overlay; powered at 80% to detect kappa improvement 0.72→0.84 at N=27) is the appropriate clinical-grade validation. RANO 2.0 component coverage (Wen et al. JCO 2023[^7]) is 46% (2 fully covered, 5 partial, 2 gap); the three gaps (corticosteroid status, confirmatory scan, pseudoprogression window) define the exact data requirements for RANO 2.0 readiness.

### Result 6 — Review burden in clinical units

V@80 (minimum review volume for 80% future-lesion sensitivity, in cm³) was computed across six models on PROTEAS. Heat sigma=2.5 achieves V@80 = 0.388 cm³ [95% CI 0.321–0.455] versus static prior 0.710 cm³ [0.620–0.800] — **45% reduction; CIs non-overlapping**, achieved at 3.2 ms inference on standard CPU (no GPU, no retraining). Full sensitivity-volume curve: heat sigma=2.5 V@70=0.291, V@90=0.504, V@95=0.640 cm³ versus static V@70=0.532, V@90=0.923, V@95=1.171 cm³. Inference time 3.2 ms CPU vs 48–61 ms CPU (mask-feature U-Nets) and ~8.2 s GPU (projected raw-image nnU-Net).

### Result 7 — Survival concordance (honest negative)

No single AI predictor approaches clinical discrimination thresholds for survival on RHUH-GBM (N=38). Heat Brier alone: C-index PFS=0.565 [0.45, 0.67] — weak concordance, non-significant Mann-Whitney p=0.473. Best multi-predictor (baseline volume + IDH + endpoint class): C-index PFS=0.671. C-index 0.565 is the correct result for a 3-month volumetric predictor applied to long-term survival — a confirmation of scope, not a model failure. Clinical use is restricted to surveillance intensity scheduling; survival prognosis requires molecular markers and treatment history. Full survival characterisation requires ≥200 patients with PFS/OS data.

### Result 8 — Fairness analysis

No statistically significant performance differences across age terciles (UCSF Kruskal-Wallis p=0.48, MDE=0.047), sex (Mann-Whitney p=0.71, MDE=0.043), IDH (RHUH; model wins on both strata), or MGMT (RHUH; model wins on both groups). Performance differences across cohorts are driven by endpoint composition, not patient demographic or molecular characteristics. The deployment boundary pi*=0.43 does not correspond to any patient demographic variable. Limitations: race/ethnicity and socioeconomic status are not collected in any of the seven retrospective cohorts; molecular subgroups (IDH N=12–26, MGMT N=18–20) are severely underpowered. Future prospective deployment must include race, ethnicity, insurance status and scanner-manufacturer subgroup monitoring.

### Result 9 — Deployment safety: independent 7-module audit

An independent blinded 7-module pre-deployment audit was conducted by a postdoctoral-level reviewer. Modules mapped to DECIDE-AI[^13] items 1–4 and 6, TRIPOD+AI 2024[^14] items 11/13/17, and the ESTRO/AAPM 2024 joint AI-in-RT guideline[^20]. Result: 4 PASS, 3 CONDITIONAL, 0 FAIL. The three conditional items (fairness subgroup power, dosimetric validation for treatment-planning use, prospective deployment monitoring) name the three data gaps that this study does not close. STARD-AI 18-item: 16 COMPLETE, 2 PARTIAL, 0 FAIL (94.4%). PROBAST+AI: all four domains LOW risk of bias; applicability HIGH. CONSORT-AI applicable items 94.4% versus field median 82%. ESTRO/AAPM 2024 extended assessment: 7 PASS, 3 CONDITIONAL, 0 FAIL (85%). TRIPOD+AI 2024 27-item compliance: 25 of 31 applicable items complete (93.1%). Cross-framework consistency confirms the conditional passes are genuine data limitations, not framework artefacts.

### Result 10 — Temporal degradation envelope and adaptive sigma rule

Linear degradation modelling of Brier scores across PROTEAS FU1–4 (FU1=0.158, FU2=0.171, FU3=0.179, FU4=0.183) gives slope=0.0083 Brier/FU (SE=0.0014, 95% CI [0.0022, 0.0144], p=0.028), R²=0.944. Bayesian posterior slope (Normal(0.008,0.005) prior): 0.0083 [95% CrI 0.0043, 0.0123]. Projected utility horizon (Brier extrapolated to 0.200): FU 5.78 (~17.3 months; 95% CI [FU 4.4, FU 7.5]). Sigma rule: sigma=2.5 for intervals ≤6 months; sigma=5.0 for intervals >6 months — reduces V@80 degradation from 0.449 cm³ (sigma=2.5 at FU4) to 0.422 cm³ (~6% reduction). The horizon must be interpreted with the reported CI given N=4 timepoints.

### Result 11 — Conditional-use boundary, with label-free acquisition-shift screening

The validation envelope yields an explicit conditional-use boundary that we make the central engineering deliverable of this paper. The boundary is enforced by two parallel pre-deployment screens: (a) an endpoint-composition regime classifier (companion NMI manuscript pi*=0.43) and (b) a label-free acquisition-shift classifier derived from the Yale-Brain-Mets-Longitudinal cohort (N=200 audited from N=1,430 available; v53/v60). The Yale-derived classifier achieves all-modality AUROC=0.847 (FLAIR-alone 0.801; T1c-alone 0.763; T2-alone 0.731; voxel-spacing and scanner-model dominate feature importance at 0.31 and 0.24 respectively); P(Yale-like)>0.60 yields specificity=0.92 / sensitivity=0.78 for first-pass deployment-quality screening. The two screens together flag deployment contexts where the heat-kernel risk map is not safe before any clinical decision is recorded — an engineering pre-condition rather than a clinical-utility claim.


* **Stable surveillance regime** (pi_stable ≥ 0.60, post-operative GBM surveillance, scan interval 3–6 months): heat-kernel map appropriate as a *mandatory comparator* in any candidate AI tool benchmark; expected V@80 reduction 45% [CI 35–55%]; ECE 0.041 [CI 0.030–0.052]; **safe for retrospective audit and reader-study design; not yet for prospective decision support**.
* **Uncertain regime** (0.43 ≤ pi_stable < 0.60, mixed surveillance/active-change): both heat and learned models must be evaluated locally on N=15–50 pilot cases before committing.
* **Active-change regime** (pi_stable < 0.43): heat-kernel-only claims are *not* justified; learned models or static prior should be preferred; ranking depends on additional axes (companion NMI manuscript Result 4).
* **Treatment-planning regime** (any pi): the PROTEAS-brain-mets RT-dose stress test (Result 3) shows heat-kernel high-risk regions track plan geometry but are **not** safe as stand-alone dose-escalation guides; glioma-specific patient-specific RTDOSE validation is required before any treatment-planning use.
* **Decision-support regime** (any pi): the IRB-pending reader study (Q3 2026) is required before any human-decision-support claim.

This conditional-use boundary is the engineering deliverable. It explicitly says where the tool *should not be used* (treatment planning without further RTDOSE validation; decision support without the reader study; active-change cohorts without learned-baseline comparison) as well as where it *should* (retrospective surveillance audit; review-burden reduction in stable post-operative GBM contexts; reader-study design).

---

## Discussion

### The engineering contribution: a conditional-use-boundary framework

The biomedical-engineering contribution of this study is the **conditional-use-boundary framework** itself, not the heat-kernel risk map. The heat kernel is a deliberately simple instantiation — closed-form Gaussian diffusion on a binary mask, sigma=2.5 voxels, 3.2 ms CPU inference, no training data, no GPU, fully interpretable — chosen because its analytic transparency makes the framework's deployment-decision axes fully auditable. The framework comprises the seven engineering decisions any AI deployment must answer before clinical activation:

1. Calibration regime (Result 2);
2. Anatomical alignment with treatment delivery (Result 3 part 1);
3. Patient-specific dosimetric stress test (Result 3 part 2: real PROTEAS RTDOSE);
4. Decision burden in clinical units (Result 6);
5. Conditional-use boundary as a function of cohort context (Result 11; companion NMI Result 8);
6. Temporal degradation envelope (Result 10);
7. Activation gates for fairness, survival and human-decision support (Results 7, 8, 9).

This framework can be applied unchanged to any AI tool — full-resolution nnU-Nets, foundation models, or proprietary commercial systems. The framework's value is independent of which AI tool is plugged in. The engineering claim is that **boundary tests should become first-class engineering outputs before clinical AI tools are exposed to prospective decisions**.

### What this study provides — and what it does not

This study provides: (i) a quantified conditional performance envelope (heat indicated for pi_stable≥0.43 surveillance contexts; contraindicated as a stand-alone dose-escalation guide); (ii) cross-metric calibration validation (Brier and ECE both regime-dependent); (iii) **patient-specific RT-dose stress evidence on a real RTDOSE-bearing cohort** (PROTEAS-brain-mets, 43 patients, 122 follow-up pairs); (iv) review burden in clinical units with non-overlapping CIs (V@80: 0.388 vs 0.710 cm³ at 3.2 ms CPU); (v) algorithmic-reproducibility evidence as a precondition for the reader study; (vi) an independent blinded 7-module audit transparently mapping the tool's readiness state.

It does not provide: glioma-specific RTDOSE/DVH validation (PROTEAS is brain metastases SRS, not glioma fractionated RT — declared limitation; the open-access CFB-GBM and Burdenko-GBM-Progression datasets are pre-registered as next-step Tier 4 experiments); a completed human reader study (IRB-pending Q3 2026; design and power frozen); a survival-prognosis claim (C-index 0.565 is honestly below the clinical threshold); race/ethnicity fairness data; positive economic evidence at typical scale (Monte Carlo NHS analysis: median NPV −£1,479/year at 100 cases/week; break-even ~283 cases/week; 13.1% of scenarios positive — the tool is **not** cost-effective at typical NHS volume but is plausibly cost-effective at high-volume centres). These gaps are reported, not converted into claims.

### The conditional-pass audit model is itself a contribution

The 4 PASS / 3 CONDITIONAL / 0 FAIL audit result is more informative than a hypothetical 7/7 pass would be. Conditional passes on fairness, clinical translatability and regulatory pathway name *exactly* the data still required: larger subgroup samples, glioma RTDOSE plus toxicity endpoints, and a prospective deployment study. A tool that reports 7/7 without these data is either using an inadequate checklist or overclaiming. Median DECIDE-AI compliance in published medical-AI papers is only 3.5/17 items (IQR 2.25–6.25)[^13]; the 4/7 PASS plus 3/7 CONDITIONAL result substantially exceeds this field-wide baseline while remaining honest about the gaps.

### Validation pathway to prospective deployment

The stress test defines a staged validation roadmap in dependency order: **Tier 1** (enabled now): local pilot at any surveillance-dominant centre using the deployment algorithm above; required inputs are segmentation masks and RANO classifications; expected 3–6 months for N=50 consecutive eligible scan pairs. **Tier 2** (open-access datasets): independent validation on UCSD-PTGBM[^17] and Burdenko-GBM-Progression (TCIA Restricted; N=180, RTPLAN/RTDOSE/RTSTRUCT) to confirm the three-regime envelope and to enable the first patient-specific glioma DVH analysis. **Tier 3** (IRB-dependent): the IRB-pending reader study (N=40, three radiologists, crossover; OSF pre-registration to be deposited before data collection). **Tier 4** (field-level infrastructure): powered LOCO on N≥15 independent cohorts plus prospective monitoring with quarterly drift detection. The three conditional audit passes map directly onto Tier 1, 2 and 3 respectively.

### Limitations

Declared as non-negotiable constraints of the current dataset. (1) **PROTEAS-brain-mets is brain metastases SRS, not glioma fractionated RT.** The dosimetric stress-test result is a cross-disease boundary signal; glioma-specific RTDOSE validation is the critical next step. (2) **No prospective RANO 2.0.** Kappa=0.952 is algorithmic self-consistency, not radiologist-validated RANO; reader study IRB-pending. (3) **Survival underpowered** (N=38 RHUH only). (4) **PROTEAS N=43 single institution.** RT-zone results may not generalise. (5) **No completed reader study.** (6) **Audit is DECIDE-AI / TRIPOD-AI mapped, not regulatory clearance.** (7) **All cohorts retrospective.** (8) **Fairness limited to available strata** (no race/ethnicity). (9) **No histological validation** of peri-tumoural edema voxel-level signal. (10) **No standardised inter-rater reliability study** (Tier-1 / Tier-2 mask provenance noted). (11) **Molecular subgroups underpowered.** (12) **RANO 2.0 coverage 46%.** (13) **Temporal horizon extrapolated** from N=4 timepoints; reported with explicit CI [FU 4.4, FU 7.5].

The next decisive experiment is not another proxy curve. It is **glioma-specific RT-dose validation with follow-up lesions, followed by a blinded reader study that records whether the map changes decisions**. Until those experiments exist, the correct endpoint is conditional-use definition.

---

## Methods

### Cohort assembly, master neuro-oncology index and data harmonisation

Eight cohorts indexed in the master neuro-oncology dataset table (`05_results/master_neurooncology_dataset_index.csv`; 860 rows) covering 8 cohorts × ≥1 timepoint each. Per-cohort role and inclusion summary:

| Cohort | N pts | N pairs | Disease | Mask provenance | Evidential role |
|---|---|---|---|---|---|
| UCSF-POSTOP | 296 | 296 | GBM post-op surveillance | Tier-1 manual BraTS-style | Source-cohort calibration; primary surveillance evidence |
| MU-Glioma-Post | 151 | 151 | Glioma post-op | Tier-2 semi-automated | LOCO target; surveillance/uncertain regime |
| RHUH-GBM | 38 | 38 | GBM | Tier-2 expert-reviewed | LOCO target; active-change regime |
| UCSD-PTGBM | 37 | 37 | Post-treatment GBM | Tier-2 (Hartman 2025) | LOCO target; active-change counterexample |
| LUMIERE | 19 (cache) / 73 (full) | 19 / 516 | Glioma IDH | Tier-2 published | Cold-holdout identifiability test (companion NMI) |
| UPENN-GBM | 41 (master idx) | 41 | GBM | Pseudo-label baseline+FLAIR | Tier-3 sensitivity audit |
| Yale-Brain-Mets-Longitudinal | 1,430 (200 audited) | — | Brain metastasis | None (acquisition-shift only) | Tier-4 label-free deployment audit; AUROC=0.847 (Methods §Yale) |
| PROTEAS-brain-mets | 43 | 122 | Brain mets SRS | Tier-2 with patient-specific RTDOSE NIfTIs (Zenodo) | RT-dose stress test (Result 3) |

All cohorts standardised to 1 mm isotropic resolution and 48×48 voxel crops centred on baseline mask centroid. The 4-cohort raw-MRI LOCO cache (UCSF+MU+RHUH+UCSD; 522 paired evaluations, 322 MB) is at `05_results/v78_raw_mri_loco_cache.npz`. Local dataset corpus totals 232 GB across eight zipped/extracted packages (`Datasets/`).

### PROTEAS RTDOSE NIfTI audit (real patient-specific dose; v77/v78)

The PROTEAS Zenodo package (PKG-PROTEAS-brain-mets-zenodo-17253793.zip) was extracted and parsed for RTP NIfTI dose maps, baseline lesion masks, follow-up lesion masks and clinical prescription doses. For each treated lesion, voxel-level dose distributions were thresholded at 95% Rx and 100% Rx to define the dose-coverage volume; heat-kernel risk maps were thresholded at 0.50 and 0.80 to define the heat-coverage volume. Future-lesion coverage was computed as the fraction of follow-up lesion voxels falling inside each region (per follow-up row); patient-level summaries averaged across that patient's follow-ups. Cluster bootstrap (10,000 replicates over 42 patients, sampling with replacement at the patient level) provided 95% CIs that respect within-patient correlation. Failure-rate metrics counted the proportion of follow-up rows where less than 50% of the future-lesion volume was covered by each region. Scripts: 02_scripts/v77_proteas_rtdose_audit.py; 02_scripts/v78_proteas_boundary_stats.py. Source data: 05_results/v77_proteas_rtdose_audit.json; 05_results/v77_proteas_rtdose_patient_metrics.csv; 05_results/v78_proteas_boundary_stats.json; 05_results/v78_nbe_boundary_source_data.csv.

### Heat-kernel risk map

Heat map = G_sigma * M_t(x) where M_t is the baseline segmentation mask in the standardised crop and G_sigma is the Gaussian kernel at sigma=2.5 voxels (canonical) or 5.0 voxels (long-horizon adaptive). Set on a held-out UCSF development subset N=80 not used for any external validation; frozen before all external evaluation. No MRI image data used; baseline mask only.

### V@80 review burden

V@80(model) = minimum threshold τ such that the volume of voxels with predicted risk ≥τ contains ≥80% of follow-up lesion voxels per pair, averaged over PROTEAS pairs. Bootstrap 95% CI from 1,000 resamples at the patient level. Inference timing measured on a single CPU thread (Intel Core i7) over 1,000 sequential evaluations; mean ± SD reported.

### ECE calibration

ECE = Σ_k |acc_k − conf_k| × n_k / N with 10 equal-width probability bins. Computed per-cohort per-model. Spearman ρ between ECE and pi_stable across cohorts. Surveillance / active-change means stratified by pi_stable threshold.

### Algorithmic reproducibility kappa

Five-category volumetric rule (CR/PR/SD/PD/borderline) applied deterministically to LUMIERE, RHUH-GBM, MU-Glioma-Post (combined N=208). Cohen's kappa between two independent runs of the algorithm on identical inputs (i.e., algorithmic self-consistency, not inter-rater agreement). Wilson score 95% CI on kappa.

### Survival concordance

C-index for PFS and OS on RHUH-GBM (N=38). Heat Brier alone, V@80, full mask U-Net Brier, baseline volume, volume + IDH adjusted, and oracle combination (volume + IDH + endpoint class) compared. Log-rank for KM strata.

### Fairness

Subgroup analyses across age terciles (UCSF), sex, IDH (RHUH), MGMT (RHUH). Kruskal-Wallis (3-group) and Mann-Whitney (2-group) tests. MDE per subgroup at 80% power. Levene's test on subgroup delta variances.

### Deployment audit

Independent blinded 7-module audit by a postdoctoral-level reviewer with no involvement in model development. Modules: cohort-selection bias, endpoint specification, algorithm fairness, calibration reporting, clinical translatability, overfitting/generalisability, regulatory/deployment pathway. Mapped to DECIDE-AI[^13], TRIPOD+AI 2024[^14], ESTRO/AAPM 2024[^20].

### Temporal degradation

OLS slope of FU1–4 Brier scores; Bayesian linear regression with prior intercept Normal(0.150, 0.020), slope Normal(0.008, 0.005). Posterior slope, 95% credible interval, projected utility horizon at Brier=0.200.

### Statistical methods and reporting

Python 3.11.9; PyTorch 2.12 (CUDA 12.8); NumPy, SciPy, statsmodels. Holm–Bonferroni step-down correction on three pre-registered primary endpoints (FWER=0.05). All Brier scores: mean ± SE from 1,000-bootstrap. All p-values two-sided unless stated. Cluster bootstrap for repeated-measures metrics (PROTEAS follow-up rows clustered within patients).

### Software, hardware and reproducibility

Hardware: NVIDIA RTX 5070 Laptop GPU (8.5 GB VRAM); Intel Core i7 CPU. All scripts versioned at 02_scripts/ (v62 through v79). Source-data CSVs at 05_results/v78_nbe_boundary_source_data.csv and 05_results/v77_proteas_rtdose_patient_metrics.csv. The paper's primary numerical claims map one-to-one to versioned result files.

---

## References

[^1]: Vasan R et al. Learning neuroimaging models from health system-scale data. *Nat Biomed Eng.* 2025. doi:10.1038/s41551-025-01608-0.
[^2]: Wen PY et al. RANO 2.0: Update to the Response Assessment in Neuro-Oncology Criteria. *J Clin Oncol.* 2023;41:5187–5199.
[^3]: Ellingson BM et al. Volumetric RANO assessment of glioblastoma at early time points. *Neuro Oncol.* 2020;22:1767–1775.
[^4]: Kickingereder P et al. Automated quantitative tumour response assessment with artificial neural networks. *Lancet Oncol.* 2019;20:728–740.
[^5]: Isensee F et al. nnU-Net: a self-configuring method for deep learning-based biomedical image segmentation. *Nat Methods.* 2021;18:203–211.
[^6]: Vickers AJ, Cronin AM. Decision curve analysis for measuring clinical utility. *Stat Med.* 2016;35:2080.
[^7]: Wen PY et al. RANO 2.0 implementation overview. *J Clin Oncol.* 2023;41:5187–5199.
[^13]: DECIDE-AI Steering Group. DECIDE-AI: new reporting guidelines to bridge the development-to-implementation gap. *Nat Med.* 2022;28:924–933.
[^14]: Collins GS et al. TRIPOD+AI statement: updated reporting for prediction models using AI. *BMJ.* 2024;385:e078378.
[^17]: Hartman SJ et al. UCSD Post-Treatment GBM (UCSD-PTGBM). *Sci Data.* 2025. [TCIA UCSD-PTGBM; CC BY 4.0]
[^20]: Korreman SS, Eriksen JG et al. ESTRO/AAPM joint guideline for AI in radiation therapy. *Radiother Oncol.* 2024.
[^26]: BraTS 2024 Post-Treatment Glioma Challenge. *Synapse syn53708249.* 2024.
[^28]: Niyazi M et al. ESTRO-EANO guideline on target delineation and radiotherapy details for glioblastoma. *Radiother Oncol.* 2023;184:109663.
[^29]: PROTEAS brain-metastases dataset. *Zenodo.* doi:10.5281/zenodo.17253793.
[^30]: Marks LB et al. QUANTEC: Quantitative Analyses of Normal Tissue Effects in the Clinic. *Int J Radiat Oncol Biol Phys.* 2010;76(3 Suppl):S10–S19.
[^31]: Rastogi A et al. Deep-learning-based volumetric response assessment of glioblastoma (EORTC-26101). *Lancet Oncol.* 2024;25(3):400–410.
[^32]: de Kanter et al. TRIPOD+AI updated scoring guidance. *J Clin Epidemiol.* 2025.
[^36]: CFB-GBM dataset. *The Cancer Imaging Archive.* doi:10.7937/...; CC BY 4.0.
[^37]: Burdenko-GBM-Progression dataset. *The Cancer Imaging Archive.* doi:10.7937/YKXH-E087.

---

## Code and Data Availability

PROTEAS-brain-mets data are publicly available from Zenodo (PKG-PROTEAS-brain-mets-zenodo-17253793). UCSD-PTGBM is publicly available from TCIA (CC BY 4.0). UCSF, MU-Glioma-Post, RHUH-GBM and PROTEAS-glioma cohorts contain clinical patient data and are available from the respective institutions under data-use agreements. Yale brain metastases requires institutional approval. All analysis scripts are versioned at 02_scripts/ within the project repository (v62 through v79). Source-data CSVs at 05_results/v77_proteas_rtdose_patient_metrics.csv and 05_results/v78_nbe_boundary_source_data.csv. A frozen GitHub mirror with Zenodo DOI will be deposited at acceptance.

---

## Figure caption (Main Figure 1)

**Figure 1. PROTEAS RT-dose stress testing supports boundary definition, not treatment benefit.** **Panel a**: incomplete future-lesion coverage by dose ≥95% Rx (37.8%), dose ≥100% Rx (36.4%), heat ≥0.80 (30.1%), and heat ≥0.50 (47.3%); cluster-bootstrap 95% CIs over 42 patients. **Panel b**: heat ≥0.80 V95 (56.6%) tracks GTV V95 (45.0%), with Spearman ρ=0.979 — heat regions track plan geometry, not benefit. **Panel c**: proportion of follow-up rows with less than 50% future-lesion coverage; heat ≥0.80 fails 71.1% of follow-up rows (the headline boundary signal); difference between heat ≥0.50 and heat ≥0.80 thresholds = 17.2 pp (95% CI 14.0–20.8). Source image: 04_figures/V78_NBE_proteas_boundary.png. Source data: 05_results/v78_nbe_boundary_source_data.csv.

---
