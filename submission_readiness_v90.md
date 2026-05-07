# v90 Practical Radiation Oncology submission readiness

Generated: 2026-05-07

## Manuscript hardening

- Corrected cohort accounting to 42 analysable patients and 121 follow-up evaluations.
- Brought the manuscript under the PRO Original Report limit: 3,850 words for abstract + body + figure captions.
- Reduced highlights to 5 bullets, all <=85 characters.
- Reframed heat >=0.50 as a retrospective audit territory, not a superiority or dose-escalation result.
- Replaced deployment language with retrospective audit and prospective reader-study language.
- Repaired early failure-pattern citations so the first cited SRS failure-pattern claims resolve to references [3-5].

## Additional checks run

- CPU rerun of fractionation-stratified analysis via `scripts/v86_fractionation_strata.py`.
- Updated rounded manuscript values from the rerun: D95 37.9%, D100 36.6%, H50 47.4%, H80 30.1%.
- Figure/source-data audit across all cited main and extended-data figures.

## Outputs

- Markdown: `manuscript/Manuscript_v90_for_PracticalRadiationOncology.md`
- PDF: `manuscript/Manuscript_v90_for_PracticalRadiationOncology.pdf`
- Audit JSON: `submission_readiness_v90.json`

## Audit result

- PDF rebuilt successfully: 21 pages.
- PRO word-count limit: pass, 3,850 / 4,000.
- Missing source images: 0.
- Missing source data files: 0.
- Highlight count/length: pass.
- Display-item count: 6 figures + 2 tables = 8, within PRO Original Report limit.
- Main residual limitation: single-institution retrospective RTDOSE evidence with no toxicity/local-control endpoint; manuscript now states this explicitly.
