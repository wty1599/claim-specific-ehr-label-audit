options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(sub("^--file=", "", arg)), "../../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/raw_temperature_completion/00_raw_paths.R"))
root <- ROOT
agg <- file.path(root, "aggregate_outputs")
first_pass <- fread(file.path(agg, "raw_old_corrected_performance.csv"))
restored <- fread(file.path(agg, "reconstructed_old_performance.csv"))
manifest <- fread(file.path(root, "provenance/input_manifest.csv"))
archive <- fread(manifest[role == "archived_performance", path])
raw <- c("feat_pen", "feat_rf")
metric_fields <- c("external_auc", "external_auc_ci_low", "external_auc_ci_high",
                  "external_auc_bootstrap_mcse", "calibration_slope",
                  "calibration_intercept", "brier_score")
check <- merge(restored, archive, by = c("outcome", "model"), suffixes = c("_restored", "_archived"))
for (field in c("internal_cv_auc", metric_fields)) {
  stopifnot(max(abs(check[[paste0(field, "_restored")]] -
                   check[[paste0(field, "_archived")]])) <= 1e-12)
}
stopifnot(nrow(check) == 8L)
final <- copy(first_pass)
final[, auc_interval_origin := "current paired patient bootstrap"]
rounding_qc <- list()
for (i in which(final$input_version == "old")) {
  z <- restored[outcome == final$outcome[i] & model == final$model[i]]
  stopifnot(nrow(z) == 1L, z$model %in% raw)
  rounding_qc[[length(rounding_qc) + 1L]] <- data.table(
    outcome = z$outcome, model = z$model,
    archived_csv_rank_auc = final$external_auc[i],
    full_precision_restored_auc = z$external_auc,
    csv_minus_full_precision_auc = final$external_auc[i] - z$external_auc,
    full_precision_archive_reproduced = TRUE)
  for (field in metric_fields) final[[field]][i] <- z[[field]]
  final$bootstrap_valid[i] <- NA_integer_
  final$bootstrap_failed[i] <- NA_integer_
  final$bootstrap_seed[i] <- NA_integer_
  final$auc_interval_origin[i] <- "original continuous RNG stream; full-precision archive exactly reproduced"
}
old <- final[input_version == "old"]
new <- final[input_version == "corrected"]
comparison <- merge(old, new, by = c("outcome", "model"), suffixes = c("_old", "_corrected"))
for (field in c("external_auc", "brier_score", "calibration_slope", "calibration_intercept")) {
  comparison[[paste0(field, "_change")]] <-
    comparison[[paste0(field, "_corrected")]] - comparison[[paste0(field, "_old")]]
}
rep <- fread(file.path(agg, "patient_bootstrap_replicate_metrics.csv"))[input_version == "corrected"]
stopifnot(nrow(rep) == 8000L)
out <- list(
  raw_old_corrected_performance_validated.csv = final,
  raw_old_corrected_comparison_validated.csv = comparison,
  archived_csv_precision_auc_qc.csv = rbindlist(rounding_qc),
  current_patient_bootstrap_replicate_metrics.csv = rep
)
for (name in names(out)) {
  path <- file.path(agg, name)
  stopifnot(!file.exists(path))
  fwrite(out[[name]], path)
}
record <- list(recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  reason = "CSV rounding perturbs a few random-forest rank ties; use exactly recovered full-precision historical point estimates and original archived intervals",
  statistical_method_changed = FALSE, source_models_refitted = FALSE,
  bootstrap_repeated = FALSE, new_corrected_estimates_changed = FALSE,
  source_prediction_tolerance_unchanged = TRUE,
  primary_performance = "aggregate_outputs/raw_old_corrected_performance_validated.csv",
  primary_comparison = "aggregate_outputs/raw_old_corrected_comparison_validated.csv",
  first_pass_outputs_retained_as_diagnostics = TRUE)
path <- file.path(root, "provenance/historical_comparison_finalization.json")
stopifnot(!file.exists(path))
write_json(record, path, pretty = TRUE, auto_unbox = TRUE)
print(rbindlist(rounding_qc))
cat("DETERMINISTIC_HISTORICAL_COMPARISON_FINALIZED\n")
