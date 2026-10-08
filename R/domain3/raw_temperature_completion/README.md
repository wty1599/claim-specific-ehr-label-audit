# Corrected-Temperature Raw Model Completion

Date: 2026-10-08. Authorized post hoc completion after the temperature-source
correction and its Base + K2 results were known. The original fitted Raw objects
were not saved; the completion reconstructs source fits and gates evaluation on
patient-aligned reproduction of all four models for both outcomes within 1e-12.

## Sources And Runtime

The original four-model order, ten-fold loops and intervening external AUC
bootstrap calls preserve the random-number sequence before the second outcome.
Raw-EN and Raw-RF use Base covariates plus 33 features, the unchanged archived
L2 `object$recipe`, elastic-net `lambda.min`, and RF class `"1"` probability.
The recovered L1 label-assignment scale is not the Raw preprocessing recipe.
No MICE, reclustering, new tuning search or external model fitting is performed.

Model recovery/evaluation requires R 4.6.0, glmnet 5.0, ranger 0.18.0 and
data.table 1.18.4 as asserted by the source. Figure rendering uses the separate
R 4.5.3 graphics runtime. Runtime separation does not change plot calculations.

## Execution Contract

Set existing external `EHR_AUDIT_RESTRICTED_DATA_ROOT` and
`EHR_AUDIT_OUTPUT_ROOT`, plus `EHR_AUDIT_REPO_ROOT`. The data root contains
`final_full.csv`, `baseline_covars.csv` and `eicu_external.csv`.
Historical artifacts reside below the output root's
`domain3_deployable_generator_v1_20260728/`; corrected Base + K2 artifacts
and intermediate figure sources reside below `domain3_source_correction/`.
`00_raw_paths.R` maps distributed code dependencies and excluded artifacts.

The ordered stages in `run.ps1` are preflight, recover, evaluate, deterministic
historical finalization (03), separate completion verification (02), and current
Figure 4 source staging (04). Select the model executable through
`EHR_AUDIT_MODEL_RSCRIPT` or `-Rscript`; no local executable path is embedded.
This public packaging did not execute any of those stages.

Restricted fits, folds and predictions remain under the output root's
`domain3_raw_temperature_completion/private_local/`. No data, aggregates,
figure inputs, logs or provenance records are distributed.

## Current Sources

The authoritative current performance/comparison outputs are
`aggregate_outputs/raw_old_corrected_performance_validated.csv`,
`raw_old_corrected_comparison_validated.csv`, and
`raw_vs_corrected_base_k2_paired_auc.csv`. The first-pass historical CSV-readback
summaries are diagnostics: tiny RF rank-tie differences are resolved by 03 using
fully reproduced historical estimates and original intervals, without refitting
or repeating bootstrap. Current corrected estimates are not changed by 03.

The archived displayed calibration intercept is offset-only CITL; its slope
comes from a joint logistic regression. The separately named joint-fit
intercept does not replace the Figure 4 convention. See the original scientific
specification and the calibration-definition addendum in `spec/`.

Current AUC and paired intervals use 1,000 patient (`uniquepid`) bootstrap
replicates per statistic, conditional on fixed recovered source fits/predictions.
Current Base + K2 and held-out hospital recalibration are not refitted here.

Step 04 byte-copies the two completed panel-C files and the remaining nine
correction-stage Figure 4 prerequisites into the external current figure source
directory. The renderer reads that complete directory, not the two-file
`figure_source/` delivery alone. Figure 4A/B/D and plotting logic are unchanged.
