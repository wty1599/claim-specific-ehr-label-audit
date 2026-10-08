# Public Source-Path Crosswalk

Date: 2026-10-08. Source labels identify authorized code/specification roles,
not patient datasets. No local workstation path or private repository name
is disclosed. Exact source/destination hashes are held in the separate local
review record; this document is an execution/source map, not an audit log.

## Correction Sources

Each source path below is relative to its authorized source family.

| Source family / relative source | Public destination | Adaptation |
| --- | --- | --- |
| Correction / `scripts/01_snapshot_consistency.R` | `R/domain3/source_correction/01_snapshot_consistency.R` | Unchanged code; explicit restricted input arguments. |
| Correction / `scripts/02_assign_corrected_external.R` | `R/domain3/source_correction/02_assign_corrected_external.R` | Public path/dependency adapter only. |
| Correction / `scripts/03_predict_corrected_external.R` | `R/domain3/source_correction/03_predict_corrected_external.R` | Public paths; only Base + K2 predictions updated. |
| Correction / `scripts/04_bootstrap_corrected_external.R` | `R/domain3/source_correction/04_bootstrap_corrected_external.R` | Public paths; archived raw-model comparators retained. |
| Correction / `scripts/05_recalibrate_corrected_external.R` | `R/domain3/source_correction/05_recalibrate_corrected_external.R` | Public paths; fixed hospital split and numerical procedures preserved. |
| Correction / `scripts/06_funnel_corrected_external.R` | `R/domain3/source_correction/06_funnel_corrected_external.R` | Environment-based defaults, external-output guard and public specification hash. |
| Correction / `scripts/07_prepare_public_funnel_corrected.R` | `R/domain3/source_correction/07_prepare_public_funnel_corrected.R` | Same public path/specification adapter; pooling logic preserved. |
| Correction / `scripts/08_prepare_figure4_corrected.R` | `R/domain3/source_correction/08_prepare_figure4_corrected.R` | External source directory; existing figure tables are prerequisites. |
| Correction / `scripts/11_transfer_eicu_missing_patterns.R` | `R/domain3/source_correction/11_transfer_eicu_missing_patterns.R` | Public paths and shared numerical utility. |
| Correction / `sql/01_eicu_temperature_gcs_readonly.sql` | `sql/source_audit/01_eicu_temperature_gcs_readonly.sql` | Read-only source query. |
| Correction / `sql/02_mimic_gcs_sofa_readonly.sql` | `sql/source_audit/02_mimic_gcs_sofa_readonly.sql` | Read-only source query. |
| Correction / `sql/03_eicu_feature_reconciliation_readonly.sql` | `sql/source_audit/03_eicu_feature_reconciliation_readonly.sql` | Read-only source query. |
| Correction / `sql/04_eicu_gcs_reconciliation_readonly.sql` | `sql/source_audit/04_eicu_gcs_reconciliation_readonly.sql` | Read-only source query. |
| Correction / `sql/05_export_eicu_temperature_private.sql` | `sql/source_audit/05_export_eicu_temperature_private.sql` | `temperature_export_path` psql variable replaces the local output location. |
| Current figure / `Figure4_project/02_draw_figure4.R` | `figures/D3/scripts/Figure4.R` | Authoritative final source is byte-identical to the earlier source; the public copy adapts paths/consumed-input records, not plotting calculations. |
| Correction / `D3_EXTERNAL_TEMPERATURE_CORRECTION_SPEC.md` | `R/domain3/source_correction/spec/D3_EXTERNAL_TEMPERATURE_CORRECTION_SPEC.md` | Original dated scientific specification. |
| Correction / `D3_MISSING_PATTERN_TRANSFER_SPEC_20261008.md` | `R/domain3/source_correction/spec/D3_MISSING_PATTERN_TRANSFER_SPEC_20261008.md` | Original dated scientific specification. |
| Funnel / `FUNNEL_POSTHOC_SPEC.md` | `R/domain3/source_correction/spec/FUNNEL_POSTHOC_SPEC.md` | Original statistical plan with explicit corrected-input note and removed private administrative paths. |

Earlier public tracked sources retain their public relative locations except
the three superseded funnel/histogram preparation files, which are replaced
by steps 06-08. No outputs or Git metadata are copied.

## Historical Raw Sources

| Relative source identifier | Retained public source | Persistence scope |
| --- | --- | --- |
| `04_analyze_d3_generator_transport_v1.R` | `R/domain3/generator/04_analyze_d3_generator_transport_v1.R` | Fits are assigned to in-memory `fit_store`; no fitted outcome-model archive is written. Historical prediction outputs are restricted prerequisites. |
| `19_external_headtohead_mice_labels_v2.R` | `R/domain3/historical_recipe/01_external_surrogate_headtohead_current.R` | Fits are local to `run()`; its return and final RDS exclude fitted outcome models. Preprocessing/partition provenance is not a fitted Raw model. |

The historical scripts explain why original Raw fitted objects were unavailable. Their archived predictions are the reproduction gate, not current corrected Raw results.

## Raw Completion Sources

| Source family / relative source | Public destination | Adaptation |
| --- | --- | --- |
| Raw completion / `scripts/01_restore_and_evaluate_raw.R` | `R/domain3/raw_temperature_completion/scripts/01_restore_and_evaluate_raw.R` | Public path/dependency mapping, public assignment-function name, and distributed addendum hash; numerical bodies and RNG sequence preserved. |
| Raw completion / `scripts/02_verify_completion.R` | Same public family / `scripts/02_verify_completion.R` | External run root and public code locator; deterministic checks preserved. |
| Raw completion / `scripts/03_finalize_historical_comparison.R` | Same public family / `scripts/03_finalize_historical_comparison.R` | Required by 02; external root only, no refit or bootstrap. |
| Raw completion / `run.ps1` | Same public family / `run.ps1` | External roots, configured executable, ordered completion/finalization/verification/staging. |
| Raw scientific plan / `RAW_MODEL_COMPLETION_SPEC.md` | Same public family / `spec/RAW_MODEL_COMPLETION_SPEC.md` | Original scientific bytes retained. |
| Raw completion / `CALIBRATION_DEFINITION_ADDENDUM_20261008.md` | Same public family / `spec/CALIBRATION_DEFINITION_ADDENDUM_20261008.md` | Public specification locator replaces workspace-relative locator; definitions unchanged. |
| Public path adapter | Same public family / `00_raw_paths.R` | Distributed historical generator/configuration/utilities, external scientific inputs and outputs. |
| Public final-source staging | Same public family / `scripts/04_stage_current_figure4.R` | Byte-copies nine unchanged prerequisites and the two completed panel-C files, without metric/plot calculations. |

Current performance: `raw_old_corrected_performance_validated.csv`; current paired comparison: `raw_vs_corrected_base_k2_paired_auc.csv`. These excluded outputs are under the external Raw-completion aggregate directory. Raw preprocessing is the unchanged archived L2 recipe, never the recovered L1 label-assignment scale.

## Public Function Mapping

| Required functions | Distributed definition |
| --- | --- |
| `apply_frozen_recipe`, `assign_nearest_centroid` | `R/domain3/generator/00_d3_generator_utils_v1.R` |
| `parse_age`, `parse_sex`, `fast_auc`, `calibration_slope_intercept` | Same shared generator utility. |
| `cluster_bootstrap_indices`, `percentile_interval` | Same shared generator utility. |
| Repository/data/output paths | `R/common/00_repo_paths.R`; correction adapter `R/domain3/source_correction/00_correction_paths.R`. |
| Main/current Figure 4 selection and actual consumed-input fingerprints | `R/common/02_figure_paths.R`; one explicit Figure 4 argument takes precedence over environment/default. |
| External-output boundary | `R/common/01_assert_external_output.R`. |
| D1 Lloyd utilities and strict decision rule | `R/domain1/empirical/00_domain1_revision_utils_v4_lloyd.R`, `01_audit_discreteness_v2.R`. |
| D1 definition-only engine loader/configuration | `R/domain1/empirical/dependencies/00_domain1_locked_engine_v3.R`, `00_domain1_revision_config.R`. |
| D1/S7 simulation source definitions | `R/simulation/core/00_main_S1_S5_formal.R`, `01_S7_discrete_outcome_null_formal.R`. |

## Restricted Artifact Locations

These are paths relative to the configured external roots, not included files.

| Root / relative location | Input or output role |
| --- | --- |
| `EHR_AUDIT_RESTRICTED_DATA_ROOT` / `eicu_external.csv`, `final_full.csv`, `baseline_covars.csv` | Authorized extracts. |
| `EHR_AUDIT_OUTPUT_ROOT` / `domain3_parameter_recovery/outputs/` | Recovered parameter object. |
| Same output root / `domain3_L1/private_local/` | Earlier fixed labels, fitted L1 models and predictions. |
| Same output root / `domain3_deployable_generator_v1_20260728/` | Historical comparator predictions and aggregate metrics. |
| Same output root / `domain3_recalibration/private_local/` | Archived patient-linked hospital split mapping. |
| Same output root / `domain3_source_correction/private_local/` | Corrected restricted temperature export, labels, predictions and hospital counts. |
| Same output root / `domain3_source_correction/aggregate_outputs/` | Corrected metrics, recalibration and funnel outputs. |
| Same output root / `domain3_source_correction/figure4/source_data/` | Earlier Base + K2-only figure stage and unchanged A/B/D prerequisites; optional `EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT`. |
| Same output root / `domain3_raw_temperature_completion/` | Verified reconstructed fits/predictions (restricted), current aggregates and two panel-C delivery files (excluded). |
| Same output root / `domain3_raw_temperature_completion/figure4/source_data/` | Complete current eleven-file Figure 4 source directory; optional `EHR_AUDIT_FIGURE4_SOURCE_ROOT`. |
| Same output root / `main_figures/source_data/` | Main-figure sources; optional `EHR_AUDIT_MAIN_FIGURE_SOURCE_ROOT`. |
| Same output root / `main_figures/figures/`, `main_figures/qa/`, `domain3_source_correction/figure4/figures/`, `domain3_source_correction/figure4/qa/` | Generated figures and actual consumed-input records, never repository members. |

Run the historical L1 preparation before correction steps 02-08, using retained
authorized inputs and their fixed fingerprints. The correction sequence is
02 assignment, 03 prediction, 04 bootstrap, 05 recalibration, 06 funnel,
07 screened funnel source, 08 earlier figure preparation. Then run Raw completion 01 (preflight/recover/evaluate), 03 deterministic historical finalization, 02 verification and 04 complete-source staging, followed by the unchanged Figure 4 renderer. Model computation uses R 4.6.0; graphics use R 4.5.3.
Step 01 is an archive snapshot check; step 11 is a separate sensitivity check.
None of these commands were executed during packaging.

`EHR_AUDIT_WORK_ROOT` supplies excluded historical work artifacts for older
families. English folder aliases such as `analysis_archive/simulations/`
replace earlier language-specific directory labels consistently; arranging
those artifacts is an authorized local staging step, not public data release.
Some historical registry/calibration inputs remain unavailable from this
code-only directory. Do not fabricate them or silently reconstruct D2 folds.
