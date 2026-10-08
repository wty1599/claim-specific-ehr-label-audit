options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(digest); library(jsonlite)})
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(sub("^--file=", "", arg)), "../../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/raw_temperature_completion/00_raw_paths.R"))
root <- ROOT
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
manifest <- fread(file.path(root, "provenance/input_manifest.csv"))
stopifnot(all(vapply(manifest$path, sha, character(1)) == manifest$sha256))
qc <- fread(file.path(root, "aggregate_outputs/old_prediction_reproduction.csv"))
stopifnot(nrow(qc) == 8L, all(qc$failures_above_tolerance == 0L),
          all(qc$max_absolute_error <= 1e-12))
perf <- fread(file.path(root, "aggregate_outputs/raw_old_corrected_performance_validated.csv"))
pair <- fread(file.path(root, "aggregate_outputs/raw_vs_corrected_base_k2_paired_auc.csv"))
rep <- fread(file.path(root, "aggregate_outputs/current_patient_bootstrap_replicate_metrics.csv"))
failed <- fread(file.path(root, "aggregate_outputs/bootstrap_failures.csv"))
stopifnot(nrow(perf) == 8L, nrow(pair) == 4L, nrow(rep) == 8000L,
          all(perf$bootstrap_requested == 1000L), all(pair$bootstrap_requested == 1000L))
for (name in c("external_auc", "brier_score", "calibration_slope", "calibration_intercept")) {
  stopifnot(all(is.finite(perf[[name]])))
}
old_path <- manifest[role == "archived_performance", path]
archive <- fread(old_path)
historical <- merge(perf[input_version == "old"], archive,
                    by = c("outcome", "model"), suffixes = c("_restored", "_archived"))
stopifnot(nrow(historical) == 4L)
for (name in c("external_auc", "brier_score", "calibration_slope", "calibration_intercept")) {
  stopifnot(max(abs(historical[[paste0(name, "_restored")]] -
                   historical[[paste0(name, "_archived")]])) <= 1e-12)
}
stats <- rep[, {
  v <- value[is.finite(value)]
  q <- quantile(v, c(.025, .975), names = FALSE, type = 7)
  .(requested = .N, valid = length(v), failed = sum(!is.finite(value)),
    ci_low = q[1], ci_high = q[2])
}, by = .(outcome, model, input_version, metric)]
stopifnot(nrow(stats) == 8L, all(stats$requested == 1000L))
for (i in which(perf$input_version == "corrected")) {
  z <- perf[i]
  q <- stats[outcome == z$outcome & model == z$model & input_version == z$input_version & metric == "auc"]
  stopifnot(nrow(q) == 1L, q$valid == z$bootstrap_valid, q$failed == z$bootstrap_failed,
            abs(q$ci_low - z$external_auc_ci_low) <= 1e-12,
            abs(q$ci_high - z$external_auc_ci_high) <= 1e-12)
}
for (i in seq_len(nrow(pair))) {
  z <- pair[i]
  q <- stats[outcome == z$outcome & model == z$model_a & metric == "paired_delta_auc"]
  stopifnot(nrow(q) == 1L, q$valid == z$bootstrap_valid, q$failed == z$bootstrap_failed,
            abs(q$ci_low - z$ci_low) <= 1e-12, abs(q$ci_high - z$ci_high) <= 1e-12)
}
stopifnot(sum(stats$failed) == nrow(failed))
figure <- fread(file.path(root, "figure_source/Figure4_panelC_L1_transport.csv"))
old_figure <- fread(manifest[role == "figure_transport", path])
stopifnot(nrow(figure) == 8L,
          identical(figure[!model %in% c("feat_pen", "feat_rf")],
                    old_figure[!model %in% c("feat_pen", "feat_rf")]),
          sha(file.path(root, "figure_source/Figure4_panelC_L1_calibration.csv")) ==
            manifest[role == "figure_calibration", sha256])
for (i in seq_len(nrow(perf[input_version == "corrected"]))) {
  z <- perf[input_version == "corrected"][i]
  f <- figure[outcome == z$outcome & model == z$model]
  for (name in c("external_auc", "external_auc_ci_low", "external_auc_ci_high",
                 "calibration_slope", "calibration_intercept", "brier_score")) {
    stopifnot(abs(f[[name]] - z[[name]]) <= 1e-12)
  }
}
sizes <- perf[, unique(.SD), .SDcols = c("outcome", "n_external_stays", "n_external_patients", "n_external_events")]
stopifnot(sizes[outcome == "mortality", n_external_stays] == 17284L,
          sizes[outcome == "mortality", n_external_events] == 2940L,
          sizes[outcome == "make", n_external_stays] == 17333L,
          sizes[outcome == "make", n_external_events] == 5714L)
code <- readLines(file.path(RAW_CODE_ROOT, "scripts/01_restore_and_evaluate_raw.R"), warn = FALSE)
stopifnot(!any(grepl("(kmeans|mice|KM)\\s*\\(", code)))
result <- list(status = "PASS", checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
               old_prediction_rows = sum(qc$n_stays), tolerance = 1e-12,
               input_hashes_unchanged = TRUE, figure_non_raw_rows_unchanged = TRUE,
               figure_calibration_bytes_unchanged = TRUE,
               bootstrap_statistics = nrow(stats), bootstrap_metric_rows = nrow(rep),
               bootstrap_failures = nrow(failed),
               historical_point_estimates_match_archive = TRUE,
               procedure = "separate deterministic source/summary check; not an independent model rerun")
output <- file.path(root, "provenance/completion_qc.json")
stopifnot(!file.exists(output))
write_json(result, output, pretty = TRUE, auto_unbox = TRUE)
print(result)
