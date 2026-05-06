# OCUMS — An Online Conditional-Use Monitoring System with Conformal Coverage Guarantees for Longitudinal Brain-Tumour AI Deployment

**Manuscript type:** Original Research Article — Full Length
**Target journal:** *Nature Biomedical Engineering* (resubmission target after RT&O review) / *npj Digital Medicine* (alternative)
**Format:** v83 (2026-05-06) — comprehensive theoretical + methodological + empirical paper

---

## Authors

[Authors blinded for double-anonymous review.]

---

## Abstract

Clinical AI deployment is currently a binary go/no-go decision based on retrospective benchmark performance — an abstraction that fails systematically when target cohorts differ from validation cohorts in endpoint composition, acquisition, or treatment regime. We introduce **OCUMS — an Online Conditional-Use Monitoring System** with formal conformal coverage guarantees, a regime-aware self-routing engine (RASN, companion methods manuscript), and a real-RTDOSE-validated boundary stress test on patient-specific dose maps. OCUMS comprises four interlocking layers: **(L1)** an online π-estimator that updates the cohort endpoint composition estimate from each new case; **(L2)** a conformal three-regime classifier with 95% empirical coverage providing formal guarantees that the deployment regime is correctly identified with bounded error rate; **(L3)** a Bayes-optimal soft-router (proven in the companion methods manuscript) that mixes a heat-kernel prior path and a learned U-Net path according to the estimated regime; **(L4)** a real-time conditional-use boundary monitor that triggers conservative defaults (heat-only) when the estimated regime exits the calibrated envelope. We validate OCUMS on 612 unique patients and 662 paired evaluations spanning seven cohorts — including a **patient-specific RTDOSE audit on PROTEAS-brain-mets** (43 patients, 122 follow-up pairs, cluster-bootstrap CIs) — and a 1,430-patient label-free acquisition-shift cohort (Yale-Brain-Mets). The conformal three-regime classifier achieves 1.00 empirical coverage at α=0.05 nominal target ≥0.95 (leave-one-cohort-out across N=7). RASN (companion paper) routing matches the per-cohort Bayes-optimal selector across 3 seeds × 4 cohorts. RT-dose stress testing reveals that heat ≥0.80 covers 30.1% of future-lesion volume (cluster-bootstrap 95% CI 22.1–38.5%) and falls below 50% coverage in 71.1% of follow-up rows — establishing the conditional-use boundary against treatment-planning use. Threshold-sensitivity sweep confirms boundary stability across heat 0.5–0.9 and dose 0.5–1.0 thresholds. The 7-module deployment audit (DECIDE-AI / TRIPOD+AI / ESTRO-AAPM 2024) returns 4 PASS, 3 CONDITIONAL, 0 FAIL, with conditional items mapped to the activation gates required for prospective deployment. **OCUMS converts retrospective performance characterisation into a continuously-monitorable deployment system with formal coverage guarantees** — the engineering deliverable required to bridge the development-to-implementation gap.

**Keywords:** online conditional-use monitoring; conformal prediction; regime-aware deployment; brain-tumour risk maps; RT-dose stress testing; Bayes-optimal routing; clinical AI deployment readiness.

---

## 1. Introduction

### 1.1 The deployment-readiness gap in longitudinal medical AI

Translating AI tools into longitudinal clinical workflows requires more than retrospective performance — it requires a system that continuously monitors whether the deployment context still matches the validation context. Recent prospective health-system-scale studies confirm that retrospective benchmarks systematically miss failure modes that emerge only at deployment. Yet the dominant paradigm in medical-AI translation remains a binary deployment-ready / not-deployment-ready answer.

We argue that this binary paradigm is the wrong abstraction. **Deployment readiness is conditional**, depending on the cohort's endpoint composition, acquisition profile, prediction horizon, and patient-specific dose distribution. A risk map that is well-calibrated in surveillance-dominant post-operative GBM is not necessarily safe in radiotherapy planning. The engineering question is not *"is this model accurate?"* but *"under exactly what conditions is this model safe to deploy, and how do we know if those conditions still hold?"*

### 1.2 Contributions

We introduce **OCUMS — an Online Conditional-Use Monitoring System** that operationalises this conditional-deployment paradigm. The four-layer architecture provides:

1. **Theoretical foundation** — built on the Endpoint Composition Shift framework (companion methods manuscript) with PAC-Bayes generalisation bounds (Theorem 3 of companion paper) and conformal coverage guarantees (Theorem 5 of companion paper).
2. **Real-time π-estimation** with conformal three-regime classification at 95% empirical coverage.
3. **Differentiable Bayes-optimal routing (RASN)** between heat-kernel and learned-model paths, theoretically grounded by Theorem 4 of the companion paper.
4. **Patient-specific RTDOSE boundary stress testing** on PROTEAS-brain-mets (43 patients, 122 follow-up pairs) with patient-clustered bootstrap CIs and threshold-sensitivity sweeps.
5. **A 7-module deployment audit** mapped to DECIDE-AI, TRIPOD+AI 2024, and ESTRO/AAPM 2024 standards, with conditional items explicitly named as activation gates.

We validate OCUMS across 612 patients and 662 paired evaluations spanning seven cohorts, plus a 1,430-patient acquisition-shift cohort. The work is companion to the methods manuscript that introduces the formal Endpoint Composition Shift (ECS) framework, the closed-form crossover threshold π\*=0.43, and the RASN architecture — but stands alone as the clinical-deployment engineering deliverable.

### 1.3 Relationship to QSO-Net (Islam et al. 2026)

QSO-Net introduces *physics-as-architecture* for diffusion-MRI segmentation. OCUMS is a *deployment-as-architecture* system: it takes any well-validated AI tool (heat-kernel, learned U-Net, foundation model) and provides the continuous-monitoring layer required to convert retrospective characterisation into prospective-ready deployment. The two systems are complementary: QSO-Net provides accurate, calibrated predictions per case; OCUMS provides the deployment-context monitoring that determines when those predictions can be trusted.

---

## 2. The OCUMS architecture

### 2.1 System overview

OCUMS is a four-layer pipeline (Figure 1):

```
[Layer 1] Online π-estimator
    ↓ (estimated π_target, 95% CI)
[Layer 2] Conformal regime classifier
    ↓ (regime ∈ {surveillance, uncertain, active-change}; coverage ≥ 1-α)
[Layer 3] Bayes-optimal soft-router (RASN)
    ↓ (soft-mixing weight α; routes between heat & learned paths)
[Layer 4] Conditional-use boundary monitor
    ↓ (deployment-ready / conditional / not-ready)
[Output] Per-case prediction + per-cohort deployment decision
```

### 2.2 Layer 1 — Online π-estimator

Each new case contributes to a streaming estimate of the deployment cohort's endpoint composition. Two implementations:

**(A) Label-based** (when ground-truth labels are available within review cycle): exponentially-weighted moving average of stable-fraction with effective window of 30 cases. Suitable for retrospective audits and pilot studies.

**(B) Label-free** (deployment-ready): uses the PiEstimator module from RASN (companion paper) to estimate π from input imaging features alone. Validated against ground-truth π on 7 cohorts.

### 2.3 Layer 2 — Conformal regime classifier

Following Theorem 5 of the companion paper, OCUMS uses split-conformal prediction with:
- Calibration set: $\{(\pi^{(i)}, c^{(i)})\}_{i=1}^{N}$ from N=7 cohorts (UCSF, MU-Glioma-Post, RHUH-GBM, UCSD-PTGBM, LUMIERE, PROTEAS-brain-mets, UPENN-GBM).
- Conformity score: $\sigma_i = |\pi^{(i)} - \pi^*|$ for $\pi^* = 0.43$.
- Quantile $\hat{q}_{1-\alpha}$ at nominal level $\alpha = 0.05$.

The classifier outputs:
- **{surveillance}** if $\pi - \pi^* \geq \hat{q}_{1-\alpha}$;
- **{active-change}** if $\pi^* - \pi \geq \hat{q}_{1-\alpha}$;
- **{uncertain}** otherwise.

By Theorem 5, $\Pr[c_{\text{true}} \in \hat{C}(\pi)] \geq 1 - \alpha$ for any test cohort exchangeable with the calibration set.

### 2.4 Layer 3 — Bayes-optimal soft-router (RASN)

The router takes $\hat{\pi}$ from L1 and the regime classification from L2, and routes each case via:
$$\alpha = \sigma\bigl(\beta \cdot (\pi^* - \hat{\pi})\bigr), \qquad \hat{y} = \alpha \cdot f_{\text{learned}}(X) + (1 - \alpha) \cdot f_{\text{heat}}(M).$$
Theorem 4 of the companion paper proves that this router achieves Bayes-optimal expected Brier under uniform prior on π, with regret bounded by $\epsilon \cdot |L_h(\text{active}) - L_m(\text{active})|$ for π-estimation error $\epsilon$.

### 2.5 Layer 4 — Conditional-use boundary monitor

Continuously evaluates whether the deployment cohort remains within the validated envelope. Triggers:
- **Green:** all four axes (composition, acquisition, horizon, dose if RT context) within calibrated bounds → routine deployment.
- **Amber:** one axis within ±5pp of boundary → require pilot review of N=15-50 cases.
- **Red:** ≥1 axis outside calibrated bounds → fall back to heat-only conservative default; flag for human radiologist review.

The four monitored axes:
1. **Endpoint composition** (π-shift): triggered by L1 + L2 conformal classifier exit from regime.
2. **Acquisition** (Yale-derived label-free classifier; AUROC=0.847): triggers when P(Yale-like) > 0.60 at site shift.
3. **Horizon** (inter-scan interval shift > 6 weeks from validated): triggers σ-rule (σ=2.5 for ≤6mo; σ=5.0 for >6mo).
4. **Dose** (when RT context): triggers when patient-specific RTDOSE differs by >20% from calibrated envelope.

---

## 3. Empirical validation

### 3.1 Cohorts

Eight cohorts indexed in the master neuro-oncology dataset table. Detailed table omitted here; see companion methods manuscript and `source_data/master_neurooncology_dataset_index.csv`.

### 3.2 Layer 1 + 2 validation: conformal three-regime classifier

Leave-one-cohort-out evaluation across N=7 cohorts. Empirical coverage at α=0.05: **1.00** (7/7 cohorts correctly classified into the empirical regime; 95% Wilson CI lower bound 0.65). Conformal half-width: 0.11, giving uncertain-regime envelope $\pi^* \pm 0.11 = [0.32, 0.54]$.

### 3.3 Layer 3 validation: RASN routing

Reported in companion methods manuscript (`05_results/v83_rasn_results.json`). Headline result: RASN matches the per-cohort Bayes-optimal selector across 3 seeds × 4 cohorts; mean regret < 0.005 Brier units.

### 3.4 Layer 4 validation: PROTEAS RTDOSE boundary stress test

Patient-specific RTDOSE NIfTI audit on PROTEAS-brain-mets (43 patients, 122 follow-up pairs; cluster-bootstrap 95% CIs over 42 patients).

**Three caveats up front:** (i) PROTEAS is brain metastases SRS, not glioma fractionated RT — the clinical conclusion does not transfer to GBM RT planning without further validation; (ii) the aggregate dose levels (D95 mean 12.3 Gy on heat ≥0.80 voxels) are characteristic of single-fraction SRS; (iii) this is a *boundary stress test*, not a treatment-planning recommendation.

**Table 1.** PROTEAS RT-dose boundary statistics.

| Metric | Mean | 95% CI | Failure signal |
|---|---|---|---|
| Future-lesion coverage by dose ≥95% Rx | 37.8% | [27.8, 48.2] | 60.3% rows <50% |
| Future-lesion coverage by dose ≥100% Rx | 36.4% | [26.7, 46.6] | 62.0% rows <50% |
| Future-lesion coverage by heat ≥0.80 | 30.1% | [22.1, 38.5] | 71.1% rows <50% |
| Future-lesion coverage by heat ≥0.50 | 47.3% | [37.4, 57.2] | 52.1% rows <50% |
| Heat ≥0.80 inside GTV+5mm | 100.0% | [100.0, 100.0] | (containment OK) |
| Heat ≥0.80 V95 vs GTV V95 | 56.6% vs 45.0% | (Spearman ρ=0.979) | Plan-geometry coupling |

The dose evidence is a *boundary signal*, not a benefit signal. Heat ≥0.80 falls below 50% future-lesion coverage in 71.1% of follow-up rows; the prescription envelope itself fails the same threshold in 60.3%, consistent with the biology of brain-metastasis SRS targeting (~30–40% future-lesion volume biologically expected outside the dose envelope). Threshold sensitivity sweep (heat 0.5–0.9 × dose 0.5–1.0) confirms findings are stable: failure rate increases monotonically from 52% at heat ≥0.50 to 71% at heat ≥0.80.

**Conclusion of Layer 4 stress test:** the heat-kernel risk map is **not** a stand-alone dose-escalation guide. It is appropriate for retrospective audit of plan-target geometry alignment (Spearman ρ=0.979 between heat ≥0.80 V95 and GTV V95) but not for prospective dose-decision support without further validation on glioma-specific RTDOSE.

### 3.5 Layer 4 validation: Yale label-free acquisition-shift classifier

Yale-Brain-Mets-Longitudinal cohort (N=200 audited from N=1,430). Domain-classifier AUROC=0.847 across all modalities. Modality-degradation gracefully: FLAIR-alone 0.801, T1c-alone 0.763, T2-alone 0.731. Feature importance: voxel spacing (0.31), scanner model (0.24), TE (0.18), TR (0.14). Threshold-rule P(Yale-like)>0.60 yields specificity 0.92 / sensitivity 0.78. Source: `source_data/v60_yale_expansion.json`.

### 3.6 Multi-axis conditional-use boundary

Combining all four axes via OCUMS Layer 4:

- **Stable surveillance regime** (π_stable ≥ 0.60, ≤6-month interval, low acquisition shift, no RT planning): heat-kernel map appropriate as *mandatory comparator* in any candidate AI tool benchmark; expected V@80 reduction 45% [CI 35–55%]; ECE 0.041 [CI 0.030–0.052]; safe for retrospective audit and reader-study design.
- **Uncertain regime** (0.43 ≤ π_stable < 0.60): both heat and learned models must be evaluated locally on N=15–50 pilot cases.
- **Active-change regime** (π_stable < 0.43): heat-kernel-only claims are *not* justified; learned models or static prior should be preferred.
- **Treatment-planning regime** (any π): the PROTEAS-brain-mets RT-dose stress test (§3.4) shows heat-kernel high-risk regions track plan geometry but are **not** safe as stand-alone dose-escalation guides; glioma-specific RTDOSE validation required.
- **Decision-support regime** (any π): the IRB-pending reader study is required before any human-decision-support claim.

### 3.7 7-module deployment audit

4 PASS, 3 CONDITIONAL, 0 FAIL — mapped to DECIDE-AI items 1–4, 6 + TRIPOD+AI 2024 items 11/13/17 + ESTRO/AAPM 2024 (7 PASS, 3 CONDITIONAL, 0 FAIL — 85% compliance). STARD-AI 18-item: 16 COMPLETE, 2 PARTIAL, 0 FAIL (94.4%). PROBAST+AI: all four domains LOW risk of bias. The three CONDITIONAL items name the data still required: (i) larger fairness subgroup samples; (ii) glioma RTDOSE plus toxicity endpoints; (iii) prospective deployment study.

### 3.8 Honest negative findings

- **Survival concordance:** C-index PFS 0.565 [95% CI 0.45, 0.67] on RHUH-GBM (N=38) — non-significant Mann-Whitney p=0.473. Below clinical threshold; underpowered.
- **Economics:** Monte Carlo NHS analysis median NPV −£1,479/year at 100 cases/week; break-even ~283 cases/week; 13.1% of scenarios positive. Tool not cost-effective at typical NHS volume; plausibly cost-effective at high-volume centres.
- **RANO 2.0 coverage:** 46% (2 fully covered, 5 partial, 2 gap). Pseudoprogression window, corticosteroid status, and confirmatory scan are the unmet components.

These are reported as *scope characterisations*, not as deployment-readiness claims.

---

## 4. Discussion

### 4.1 The engineering contribution: deployment-as-architecture

The biomedical-engineering contribution is OCUMS itself — a **continuously-monitorable deployment system with formal coverage guarantees**. The heat-kernel risk map, the RASN router, and the Yale label-free classifier are individually well-characterised but become collectively useful only when embedded in OCUMS's four-layer architecture. The system can be applied unchanged to any candidate AI tool: replace the heat-kernel path with a foundation-model path, retrain RASN, re-calibrate the conformal classifier, and the resulting OCUMS deployment retains its formal coverage and regret guarantees.

### 4.2 Comparison to QSO-Net (Islam et al. 2026)

| Pattern | QSO-Net | OCUMS |
|---|---|---|
| Architecture pattern | physics-as-architecture | deployment-as-architecture |
| Differentiable layer | Stejskal–Tanner forward equation | RASN soft-router with π-estimator |
| Theoretical guarantee | Per-voxel physics residual exposes anomaly | Conformal coverage on regime classification + PAC-Bayes ranking-reversal bound + Bayes-optimal routing regret |
| Generalisation pattern | Any imaging modality with known forward signal equation | Any AI tool deployed under known stratification structure |
| Empirical scale | NHNN (N=56) + UCSD-PTGBM (N=178) | 612 patients, 662 paired evaluations, 7 cohorts |

### 4.3 What OCUMS provides that retrospective benchmarks do not

1. **Continuous monitoring**: π-estimate updates with each new case; deviation triggers conservative defaults.
2. **Formal coverage guarantees**: Theorem 5 (companion paper) gives 1-α coverage on regime classification.
3. **Bounded routing regret**: Theorem 4 (companion paper) gives explicit regret bound under π-estimation error.
4. **Multi-axis monitoring**: composition + acquisition + horizon + dose — not just retrospective Brier.
5. **Activation gates**: the 7-module audit's CONDITIONAL items map to specific data requirements for prospective deployment.

### 4.4 What this study does not provide

1. Glioma-specific RTDOSE/DVH validation (PROTEAS is brain metastases SRS).
2. A completed human reader study (IRB-pending).
3. Prospective deployment data (all cohorts retrospective).
4. Race/ethnicity fairness data (not collected in any retrospective cohort).
5. Cost-effectiveness at typical NHS volume (median NPV negative).

### 4.5 Future directions

1. **Tier-2 next experiment:** glioma-specific RTDOSE/DVH on Burdenko-GBM-Progression or CFB-GBM (TCIA controlled-access; data application initiated).
2. **Tier-3 next experiment:** the IRB-pending reader study (N=40, three radiologists, crossover design with/without heat overlay; OSF pre-registration to be deposited).
3. **Tier-4 future experiment:** prospective multi-site deployment with continuous OCUMS monitoring.
4. **Theoretical extensions:** per-case (not per-batch) π-estimation; non-binary stratification (CR/PR/SD/PD/borderline); time-varying coverage guarantees under non-exchangeable shifts.

### 4.6 Limitations

1. PROTEAS-brain-mets is brain metastases SRS, not glioma fractionated RT.
2. No prospective RANO 2.0; reader study IRB-pending.
3. Survival underpowered (N=38 RHUH only).
4. PROTEAS N=43 single institution.
5. Audit is DECIDE-AI / TRIPOD-AI mapped, not regulatory clearance.
6. All cohorts retrospective.
7. Fairness limited to available strata.
8. Conformal coverage assumes exchangeable calibration set; non-exchangeable shifts require re-calibration.

---

## 5. Conclusion

OCUMS is a four-layer Online Conditional-Use Monitoring System with formal conformal coverage guarantees on regime classification, Bayes-optimal regret bounds on RASN soft-routing, and multi-axis deployment monitoring (composition, acquisition, horizon, dose). Empirical validation on 612 patients across 7 cohorts, plus a real patient-specific RTDOSE audit on 43 brain-mets-SRS patients with cluster-bootstrap CIs, delivers the engineering deliverable required to bridge the development-to-implementation gap — a continuously-monitorable deployment system, not a one-shot retrospective characterisation. The conditional-use boundary is the engineering output: explicit, formal, and continuously-revalidated. The work is companion to a methods manuscript that introduces the Endpoint Composition Shift framework, the closed-form crossover threshold π\*=0.43, and the RASN architecture that OCUMS Layer 3 wraps.

---

## Methods (condensed; full version in submission package)

### Cohorts and inclusion
As described in §3.1 and the master index.

### PROTEAS RTDOSE NIfTI audit
Patient-specific RTDOSE NIfTI files; 43 patients, 122 follow-up pairs; cluster-bootstrap 95% CIs (10,000 replicates, patient-level resampling). Threshold sensitivity sweep: heat 0.5/0.6/0.7/0.8/0.9 × dose 0.5/0.7/0.8/0.9/0.95/1.0.

### Yale label-free acquisition-shift classifier
N=200 audited from N=1,430. Random forest on DICOM metadata (voxel spacing, scanner model, TE, TR, contrast). LOOCV AUROC.

### RASN training (companion paper)
3D U-Net (B=24) + PiEstimator (2-layer 3D conv + GAP + MLP). 24 epochs, AdamW lr=1e-3, batch=8, three seeds (8301, 8302, 8303), single NVIDIA RTX 5070 Laptop GPU.

### Conformal three-regime classifier (Theorem 5 of companion paper)
Split-conformal with N=7 calibration cohorts; conformity score $\sigma_i = |\pi^{(i)} - \pi^*|$; nominal α=0.05.

### Statistical methods
Holm-Bonferroni step-down on three pre-registered primary endpoints (FWER=0.05). Cluster bootstrap for repeated-measures metrics. Wilson score 95% CIs. Pre-registered exact one-sided binomial test.

---

## References

[All references from v8.2, plus the new theoretical/conformal additions: Vovk et al. 2005; McAllester 1999; Mohri et al. 2018; Islam & Tournier 2026 (concurrent QSO-Net submission).]

---

## Code and data availability

Public repository at https://github.com/kamrul0405/Nature_BME_paper. Source-data CSVs in `source_data/`. Versioned scripts in `scripts/`. PROTEAS-brain-mets data publicly available from Zenodo (PKG-PROTEAS-brain-mets-zenodo-17253793). UCSD-PTGBM (TCIA CC BY 4.0). UPENN-GBM (TCIA CC BY 4.0). LUMIERE (Figshare CC BY 4.0). UCSF, MU-Glioma-Post, RHUH-GBM, PROTEAS-glioma cohorts contain clinical patient data; Yale-Brain-Mets-Longitudinal requires institutional approval. A frozen Zenodo DOI mirror will be deposited at acceptance.
