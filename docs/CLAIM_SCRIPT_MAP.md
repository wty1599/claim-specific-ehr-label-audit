# Script-to-claim map

This map identifies why each retained script family is present. Shared utilities and quality-control scripts support the listed analysis; they are not independent positive findings. Numerical outputs require credentialed inputs that are not in this repository.

| Manuscript claim or artifact | Retained code | Scope |
| --- | --- | --- |
| Fixed K=2 construction and source-cohort definition | `sql/`, `R/partition/`, `R/common/` | Source extraction, preprocessing, MICE completed data, and fixed Lloyd labels. |
| Composition analysis | `R/step0/` | Cohort profile, renal and acid-base loading, identity completion, ablation, KDIGO and time-origin checks. |
| D1: assignment separation is not proof of discrete classes | `R/domain1/empirical/`, `R/domain1/sensitivity/` | Copula reference, shape/separation, split-allocation, Gap/SigClust-like and specified sensitivities. Files `05` and `06` document the historical comparison that motivated the final rule. |
| D2: incremental information depends on the comparator | `R/domain2/` | Locked-label, imputation-specific, fold-wise and landmark comparisons; files `05` and `06` assemble their reported canonical outputs. |
| D3: corrected technical assignment and limited local updating | `R/domain3/source_correction/02_assign_corrected_external.R` through `05_recalibrate_corrected_external.R` | Post hoc temperature correction, Base + K2 prediction update and hospital-held-out mortality recalibration; fixed MIMIC-IV labels/recipe/centroids. This is the earlier Base + K2-only stage; final current Raw evaluation and paired contrasts come from the subsequent completion. |
| D3: current corrected Raw outcome-model transport | `R/domain3/raw_temperature_completion/scripts/01_restore_and_evaluate_raw.R`, `03_finalize_historical_comparison.R`, `02_verify_completion.R` | Verified source-fit reconstruction with the archived L2 recipe, corrected Raw predictions for both external outcomes, patient-bootstrap AUC intervals and current paired comparisons. Historical CSV-readback diagnostics are not the validated main comparison. |
| D3: source reconciliation | `sql/source_audit/`, `R/domain3/source_correction/01_snapshot_consistency.R` | Read-only temperature/GCS checks and restricted export; existing GCS and modified SOFA remain unchanged in the temperature correction. |
| D3: hospital prevalence beyond binomial sampling | `R/domain3/source_correction/06_funnel_corrected_external.R`, `07_prepare_public_funnel_corrected.R` | Post hoc corrected-input funnel and disclosure-screened plotting sources; not a causal hospital effect. |
| D3: mechanical assignment under external missingness patterns | `R/domain3/source_correction/11_transfer_eicu_missing_patterns.R` | Frozen post hoc mask-transfer check, not external accuracy, semantic validation, or attribution of observed prevalence differences. |
| D3: retained historical sources | `R/domain3/generator/`, `deployed_rule/`, `parameter_recovery/`, `recalibration/`, `historical_recipe/` | Earlier L1 preparation and L2 comparators. Historical external results are not automatically temperature-corrected; source_correction is the current correction path. The two historical Raw fitting scripts do not persist their fitted outcome models; archived predictions are not saved fits. |
| D4: original RRT comparison support | `sql/domain4/`, `R/domain4/landmark/` | Hour-24 risk set, event ledger, allocation scores, prespecified criteria and independent QC. No RRT effect model. |
| D4: later arm-specific and mechanism checks | `R/domain4/arm_support/`, `R/domain4/posthoc/`, `R/domain4/allocation_simulation/` | Distinguishes scarcity, weight concentration, balance and timing; allocation-only simulations do not evaluate the complete original stopping rule. The loop-diuretic code is a timing feasibility check, not HTE. |
| Known-truth controls and supplementary sensitivities | `R/simulation/` | D1-D4 control behavior and reported SA analyses. The simulated patient records themselves are excluded. |
| Current article and supplementary figures | `figures/main/scripts/`, `figures/D3/scripts/`, `figures/supplement/` | Final renderers and required shared styling/preparation functions; see `FIGURE_CROSSWALK.md`. No plotting datasets or rendered images are included. |
| Release-source verification | `tests/` | Syntax and file-boundary checks only. |

The retained code includes analysis utilities and historical comparators explicitly reported in the supplement. It does not imply a one-command replay or that all historical versions are published.

Date: 2026-10-08. This is a single sepsis-partition illustration, not multi-partition or cross-disease validation. No analysis was executed during code packaging; original D2 fold assignments remain unavailable.
