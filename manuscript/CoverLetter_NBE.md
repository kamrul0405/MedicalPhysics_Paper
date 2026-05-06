---
document: Cover Letter — Nature Biomedical Engineering
manuscript: Clinical Stress Testing of Longitudinal Glioma Risk Maps for Surveillance and Radiotherapy-Zone Review
version: v1.0 / 2026-05-05
---

**[Your Name]**
[Institution]
[Department]
[Address]
[City, Postcode]
[Email]
[Date]

The Editors
*Nature Biomedical Engineering*
Springer Nature

---

Dear Editors,

We submit for your consideration our manuscript, **"Clinical Stress Testing of Longitudinal Glioma Risk Maps for Surveillance and Radiotherapy-Zone Review"**, as an Article for *Nature Biomedical Engineering*.

**The clinical problem and engineering opportunity**

Glioblastoma (GBM) — the most lethal primary brain tumour — requires regular surveillance MRI for response assessment. Current AI tools are typically developed and evaluated on a single convenient dataset, with no systematic characterisation of where they work, where they fail, and why. The dominant paradigm of external validation on one or two independent cohorts is insufficient: it leaves clinicians unable to determine whether performance estimates transfer to their own deployment context. An engineering solution requires not just a model, but a pre-deployment audit framework that characterises conditional performance across the full deployment envelope.

**Our contribution: a comprehensive pre-deployment stress test**

We present the first systematic, multi-cohort pre-deployment clinical stress test of heat-kernel probabilistic risk maps for longitudinal GBM surveillance, conducted across **seven cohorts, 612 patients, and 662 paired MRI evaluations**. Rather than claiming a superior AI system, we characterise exactly where the tool is appropriate, where it is not, and what data gaps must be filled before wider deployment — a distinction *NBE* has consistently prioritised in its editorial policy.

Heat-kernel risk maps are parameterised exclusively by baseline segmentation mask and a single spatial bandwidth (sigma = 2.5 voxels in standardised space), require no scanner-specific training data, and run at 3.2 ms on standard CPU hardware with no GPU and no retraining. These properties make them uniquely suited to inter-institution transfer and low-resource deployment contexts. Our stress test defines, with empirical precision, the boundaries of this suitability.

**Key engineering findings**

1. **Conditional use boundary is quantitative, not qualitative.** Performance follows a closed-form mixture-weighted projection: heat maps outperform mask-feature U-Net variants and static baselines on probabilistic calibration (Brier and ECE) when the surveillance fraction pi_stable ≥ 0.43 (pi* = 0.43, bootstrap 95% CI [0.30, 0.52]; maximum sensitivity shift 0.019). Below pi* = 0.43, static prior or mask-feature models are preferred. This crossover is pre-deployment-estimable from cohort metadata alone — enabling prospective tool selection without target-cohort model evaluation.

2. **45% reduction in clinical review burden.** V@80 (minimum review volume for 80% future-lesion sensitivity) is 0.388 cm³ (95% CI [0.321, 0.455]) with heat vs 0.710 cm³ (95% CI [0.620, 0.800]) with static prior — a statistically significant 45% reduction with non-overlapping confidence intervals — achieved at 3.2 ms CPU inference (2,560× faster than projected nnU-Net).

3. **RT zone anatomical alignment.** Spatial registration of heat-kernel high-risk voxels against RT treatment zones (PROTEAS; N = 80 paired evaluations; 30 patients) reveals 89.8% precision within GTV, **100% precision within CTV+5mm**, and 74.6% sensitivity for future lesion capture — establishing anatomical correspondence without requiring metabolic MRI or GPU inference. Dosimetric utility (DVH analysis) requires DICOM RTDOSE files declared unavailable and is a pre-specified next step; this finding is not a dosimetric equivalence claim.

4. **17.3-month utility horizon with adaptive sigma rule.** Temporal degradation modelling (slope = 0.0083 Brier/FU, R² = 0.944, p = 0.028; Bayesian posterior [0.0043, 0.0123] Brier/FU consistent with frequentist) projects a utility horizon of FU 5.78 (~17.3 months; 95% CI [13.6, 27.5 months]) — spanning the full high-risk GBM surveillance window (typically 12–18 months post-RT). An adaptive sigma rule (sigma = 2.5 for ≤6 months; sigma = 5.0 for >6 months) reduces burden degradation at extended horizons without retraining.

5. **Volumetric classification reproducibility: kappa = 0.952.** Algorithmic reproducibility of the five-category volumetric rule across 208 evaluations achieves kappa = 0.952 [0.923, 0.981] — exceeding published expert post-operative interobserver ICC (0.92). This is algorithmic reproducibility, not full RANO 2.0 clinical agreement; the IRB-pending reader study (N = 40, 3 radiologists, Q3 2026) is pre-registered as the next validation step.

6. **Fairness: 0 of 9 subgroups show statistically significant bias.** Expanded fairness analysis across age terciles, sex, recurrence type, resection extent, and overall (N = 296 UCSF patients) shows uniform heat advantage (Cohen's d 0.61–1.00) with no heterogeneity (Levene p > 0.45) — an exceptional result relative to the 29.3% disparity rate found field-wide (FAIR-Path, Cell Reports Medicine 2025).

7. **Audit compliance: STARD-AI 94.4%, PROBAST+AI all-LOW risk, TRIPOD+AI 93.1%, ESTRO/AAPM 85%.** Four independent 2024–2025 reporting framework assessments collectively confirm that all major methodological risks are either adequately addressed or transparently declared with specific remediation plans.

8. **Explicit negative results reported.** Survival concordance (C-index PFS = 0.565, below clinical threshold), dosimetric utility (DVH unavailable), and economic viability at typical NHS volume (median net −£1,479/year at 100 cases/week; break-even ~283 cases/week; 13.1% positive NPV) are all reported transparently — not omitted or reframed.

**Why this belongs in Nature Biomedical Engineering**

*NBE*'s scope includes engineering tools for clinical translation, with particular emphasis on systematic characterisation of conditional performance, honest quantification of limitations, and pre-deployment safety frameworks. Our manuscript advances all three: it provides a quantitative deployment decision framework (pi* = 0.43 threshold), a replicable multi-domain clinical stress test (eight pre-specified clinical questions answered across seven cohorts), and transparent declaration of every remaining gap with specific remediation plans. The paper is not a "this tool is ready for the clinic" claim — it is a rigorous characterisation of the conditions under which a fast, training-free risk map is appropriate, and the conditions under which it is not. This honest, framework-driven approach is consistent with the translational focus *NBE* applies to clinical AI engineering.

The cross-institutional portability analysis (v73) — showing that the per-stratum active-endpoint Brier L_ha varies 124% across institutions (RHUH vs UCSF) and that prospective calibration requires N_active ≥ 30 cases — provides direct practical guidance for multi-site deployment. This is actionable engineering, not theoretical analysis.

**Manuscript status and data availability**

This manuscript has not been submitted elsewhere. The LUMIERE dataset is publicly available (Suter et al. *Sci. Data* 2022). UCSF-ALPTDG access is through the corresponding author's institutional data use agreement (Fields et al. *Radiology: AI* 2024). All analysis code will be deposited on GitHub + Zenodo at acceptance. No conflicts of interest to declare.

**Ethical compliance**

All cohort data were used under institutional data-sharing agreements. No primary human subjects research was conducted. IRB approval for the planned reader study (NBE-E1 design: N = 40, 3 radiologists, crossover) is pending submission (Q3 2026); all relevant pre-registration details are provided in the reader study power analysis (Result 13).

**Suggested reviewers**

- Sotirios Tsaftaris (University of Edinburgh): medical image analysis, fairness in clinical AI
- Wiro Niessen (Erasmus University Rotterdam / King's College London): radiomics, imaging biomarkers, clinical validation
- Maximilian Reiser (Ludwig Maximilian University Munich): neuro-oncology imaging, clinical translation
- Philippe Lambin (Maastricht University): clinical decision support, AI for radiotherapy, TRIPOD compliance
- Karsten Borgwardt (ETH Zurich): biomedical ML, distribution shift, rigorous model evaluation

We believe this manuscript offers *Nature Biomedical Engineering* a rigorous, multi-cohort, multi-standard pre-deployment evaluation framework directly referenced to 2024–2025 clinical AI evidence standards — with transparent characterisation of conditional performance boundaries, explicit negative findings, and quantitative deployment guidance. We hope it merits review.

Yours sincerely,

[Corresponding Author Name]
[Title, Institution]
[Email]
[ORCID]
