# Analyses added after initial results

This summary follows Supplementary Table S23 of the current manuscript. Dates are specification dates reported there, not dates assigned during this public-code preparation. An analysis specified before its **own** calculation after an earlier result was known is not an initial preregistration.

| Analysis | Date reported in Table S23 | Reason and scope |
| --- | --- | --- |
| D4 arm-specific support diagnostics | 2026-09-24 | Separated sparse RRT records from within-arm weight concentration after the original pooled-ESS result. |
| D4 allocation-only simulation | 2026-09-24 | Examined treatment frequency, sample size, and allocation selection; did not test the full stopping rule. |
| D3 parameter recovery and L1 deployment | 2026-09-25 | Recovered parameters after the earlier reconstructed-rule analysis; external semantic equivalence remains unestablished. |
| D4 alternative-partition comparison | 2026-09-25 | Tested the renal-input mechanism hypothesis; the prediction was not met. |
| D4 no-SOFA sensitivity | 2026-09-25 | Addressed timing of a baseline severity covariate; residual imbalance persisted. |
| Loop-diuretic timing check | 2026-09-25 | Post hoc feasibility check stopped before propensity or treatment-effect modeling. |
| D3 hospital-held-out mortality recalibration | 2026-09-26 | Added after unrecalibrated external results; source models and label rule remained fixed. |
| D3 hospital-prevalence funnel analysis | 2026-10-04 | Added after the unadjusted 69/199 hospital statistic and disclosure review; tested variation against binomial sampling limits without changing labels or mortality models. |
| Corrected Raw-model completion | 2026-10-08 | Authorized post hoc reconstruction with historical prediction agreement, unchanged L2 Raw recipe, corrected-temperature evaluation and current patient-bootstrap paired contrasts; current Figure 4C uses the completed Raw delivery. R 4.6.0 model runtime is separate from R 4.5.3 graphics. |
| External temperature-source correction | 2026-10-08 | Post hoc nursing C/F priority and periodic fallback, with unchanged fixed source coordinates/centroids. Base + K2 external predictions, bootstrap, recalibration and funnel sources are updated. At this earlier stage Raw-EN/Raw-RF were retained historical references; the separately listed same-date Raw completion replaces their current external sources. |
| External missing-pattern transfer | 2026-10-08 | Frozen before its own calculation but after prior results; technical mask-induced reassignment, not external accuracy or outcome-model evaluation. |

The fuller definitions, pre-analysis specification records where available, and resulting interpretation are in Supplementary Table S23. The post hoc loop-diuretic feasibility check stopped at its prespecified timing screen; no support diagnostics or treatment-effect models were fitted.
