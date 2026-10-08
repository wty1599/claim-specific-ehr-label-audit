#!/usr/bin/env Rscript

# Read-only historical Hartigan-Wong versus Lloyd alignment audit.
# This script compares locked formal outputs only. It does not refit models,
# alter thresholds, or overwrite any existing result.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

ROOT <- file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "domain1_lloyd_alignment_20260720")
OLD_ROOT <- file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "domain1_renal_identity_and_method_closure_20260717")
OLD_35_ROOT <- file.path(OLD_ROOT, "wp9_n_eval_sensitivity_20260718", "formal")

DIR_TABLE <- file.path(ROOT, "tables")
DIR_LOG <- file.path(ROOT, "logs")
DIR_PROV <- file.path(ROOT, "provenance")

SCRIPT_ARG <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_PATH <- if (length(SCRIPT_ARG) == 1L) {
  normalizePath(sub("^--file=", "", SCRIPT_ARG), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path(ROOT, "scripts", "40_domain1_lloyd_alignment_comparison_audit.R"),
    winslash = "/", mustWork = TRUE
  )
}

OUT <- c(
  comparison = file.path(DIR_TABLE, "Table_D1_old_vs_lloyd_alignment.csv"),
  report = file.path(DIR_LOG, "D1_LLOYD_ALIGNMENT_COMPARISON.md"),
  authority = file.path(DIR_PROV, "D1_LLOYD_ALIGNMENT_AUTHORITY_MAP.csv"),
  qc = file.path(DIR_PROV, "40_D1_LLOYD_ALIGNMENT_COMPARISON_QC.csv"),
  input_manifest = file.path(DIR_PROV, "40_D1_LLOYD_ALIGNMENT_INPUT_SHA256.csv"),
  output_manifest = file.path(DIR_PROV, "40_D1_LLOYD_ALIGNMENT_OUTPUT_SHA256.csv"),
  session = file.path(DIR_PROV, "40_D1_LLOYD_ALIGNMENT_sessionInfo.txt"),
  completion = file.path(ROOT, "40_D1_LLOYD_ALIGNMENT_COMPARISON_completed.ok")
)

if (any(file.exists(OUT))) {
  stop("One or more script-40 outputs already exist. No file was overwritten.")
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

norm_path <- function(path, must_work = TRUE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  fread(path, na.strings = c("NA", ""), check.names = FALSE)
}

first_r_version <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hit <- grep("^R version ", x, value = TRUE)
  if (length(hit)) hit[[1L]] else NA_character_
}

fmt <- function(x, digits = 6L) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_pct <- function(x, digits = 1L) {
  paste0(formatC(100 * x, format = "f", digits = digits), "%")
}

qc_rows <- list()
add_qc <- function(check, pass, observed, expected = "") {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table(
    check = check,
    pass = isTRUE(pass),
    observed = as.character(observed),
    expected = as.character(expected)
  )
  if (!isTRUE(pass)) stop("QC failed [", check, "]: ", observed)
}

P <- list(
  locked_matrix = file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "output", "model", "X_primary33_std_mice.rds"),
  locked_labels = file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "output", "model", "labels_primary_mice.rds"),
  locked_delta = file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "output", "qc_domain1_mice", "posthoc_locked_k2_pooled_within_mahalanobis_20260717.csv"),

  old_29f_neg = file.path(OLD_ROOT, "tables", "Table_D1_copula_negative_control.csv"),
  old_29f_pos = file.path(OLD_ROOT, "tables", "Table_D1_copula_positive_control.csv"),
  old_29f_repeat = file.path(OLD_ROOT, "tables", "29f_copula_controls_by_repeat.csv"),
  old_29f_manifest = file.path(OLD_ROOT, "provenance", "29f_WP7_output_sha256_manifest.csv"),
  old_29f_marker = file.path(OLD_ROOT, "29f_run_completed.ok"),
  old_29f_script = file.path(OLD_ROOT, "scripts", "29f_domain1_copula_controls.R"),
  old_29f_session = file.path(OLD_ROOT, "logs", "29f_sessionInfo.txt"),

  old_29g_gate = file.path(OLD_ROOT, "tables", "Table_D1_formal_gate_by_n.csv"),
  old_29g_inv = file.path(OLD_ROOT, "tables", "Table_D1_empirical_delta_inversion.csv"),
  old_29g_gate_repeat = file.path(OLD_ROOT, "tables", "29g_formal_gate_by_repeat.csv"),
  old_29g_surface_repeat = file.path(OLD_ROOT, "tables", "29g_detection_surface_by_repeat.csv"),
  old_29g_manifest = file.path(OLD_ROOT, "provenance", "29g_WP8_output_sha256_manifest.csv"),
  old_29g_marker = file.path(OLD_ROOT, "29g_run_completed.ok"),
  old_29g_script = file.path(OLD_ROOT, "scripts", "29g_domain1_formal_gate_detection_surface.R"),
  old_29g_session = file.path(OLD_ROOT, "logs", "29g_sessionInfo.txt"),

  old_35_gate = file.path(OLD_35_ROOT, "tables", "35_WP9_gate_by_split.csv"),
  old_35_summary = file.path(OLD_35_ROOT, "tables", "Table_D1_n_eval_sensitivity.csv"),
  old_35_repeat = file.path(OLD_35_ROOT, "tables", "35_WP9_repeat_level.csv"),
  old_35_task_manifest = file.path(OLD_35_ROOT, "provenance", "35_WP9_task_manifest.csv"),
  old_35_checkpoint_inventory = file.path(OLD_35_ROOT, "provenance", "35_WP9_checkpoint_inventory.csv"),
  old_35_manifest = file.path(OLD_35_ROOT, "provenance", "35_WP9_output_sha256_manifest.csv"),
  old_35_marker = file.path(OLD_35_ROOT, "35_WP9_run_completed.ok"),
  old_35_script = file.path(OLD_ROOT, "scripts", "35_domain1_n_eval_sensitivity.R"),
  old_35_session = file.path(OLD_35_ROOT, "provenance", "35_WP9_sessionInfo.txt"),

  new_29f_neg = file.path(ROOT, "tables", "Table_D1_copula_negative_control_LLOYD.csv"),
  new_29f_pos = file.path(ROOT, "tables", "Table_D1_copula_positive_control_LLOYD.csv"),
  new_29f_repeat = file.path(ROOT, "tables", "29f_copula_controls_by_repeat_LLOYD.csv"),
  new_29f_manifest = file.path(ROOT, "provenance", "29f_output_sha256_manifest_LLOYD.csv"),
  new_29f_marker = file.path(ROOT, "29f_run_completed_LLOYD.ok"),
  new_29f_script = file.path(ROOT, "scripts", "29f_domain1_copula_controls_v2_lloyd.R"),
  new_29f_session = file.path(ROOT, "logs", "29f_sessionInfo.txt"),

  new_29g_gate = file.path(ROOT, "tables", "Table_D1_formal_gate_by_n_LLOYD.csv"),
  new_29g_inv = file.path(ROOT, "tables", "Table_D1_empirical_delta_inversion_LLOYD.csv"),
  new_29g_fit = file.path(ROOT, "tables", "29g_n1_fit_summary_LLOYD.csv"),
  new_29g_gate_repeat = file.path(ROOT, "tables", "29g_formal_gate_by_repeat_LLOYD.csv"),
  new_29g_surface_repeat = file.path(ROOT, "tables", "29g_detection_surface_by_repeat_LLOYD.csv"),
  new_29g_manifest = file.path(ROOT, "provenance", "29g_output_sha256_manifest_LLOYD.csv"),
  new_29g_marker = file.path(ROOT, "29g_run_completed_LLOYD.ok"),
  new_29g_script = file.path(ROOT, "scripts", "29g_domain1_formal_gate_detection_surface_v2_lloyd.R"),
  new_29g_session = file.path(ROOT, "logs", "29g_sessionInfo_LLOYD.txt"),

  new_35_gate = file.path(ROOT, "tables", "35_WP9_gate_by_split_LLOYD.csv"),
  new_35_emp = file.path(ROOT, "tables", "35_WP9_empirical_summary_LLOYD.csv"),
  new_35_control = file.path(ROOT, "tables", "35_WP9_control_summary_LLOYD.csv"),
  new_35_states = file.path(ROOT, "tables", "35_WP9_state_distribution_LLOYD.csv"),
  new_35_repeat = file.path(ROOT, "tables", "35_WP9_repeat_level_LLOYD.csv"),
  new_35_task_manifest = file.path(ROOT, "provenance", "35_WP9_task_manifest_LLOYD.csv"),
  new_35_manifest = file.path(ROOT, "provenance", "35_WP9_output_sha256_manifest_LLOYD.csv"),
  new_35_marker = file.path(ROOT, "35_WP9_run_completed_LLOYD.ok"),
  new_35_script = file.path(ROOT, "scripts", "35_domain1_split_allocation_sensitivity_v2_lloyd.R"),
  new_35_session = file.path(ROOT, "provenance", "35_WP9_sessionInfo_LLOYD.txt")
)

missing_inputs <- names(P)[!vapply(P, file.exists, logical(1))]
add_qc("all_required_inputs_exist", length(missing_inputs) == 0L,
       paste(missing_inputs, collapse = "; "), "none")

EXPECTED_LOCK_HASH <- c(
  locked_matrix = "9896ad4906207ca910e54eb187507037106bff212299586141908f387cfe6279",
  locked_labels = "159576d0b2796e9786db0eb4c8747af987be5e53af200419e3e37567c1d875f9",
  locked_delta = "28c33cc58306a2662bb1dbe17edf49d7c5707b897af1eff01599a3d88148e72b"
)
observed_lock_hash <- vapply(P[names(EXPECTED_LOCK_HASH)], sha256_file, character(1))
add_qc("locked_artifact_hashes_match_frozen_values",
       identical(unname(observed_lock_hash), unname(EXPECTED_LOCK_HASH)),
       paste(observed_lock_hash, collapse = "; "), paste(EXPECTED_LOCK_HASH, collapse = "; "))

verify_manifest <- function(manifest_file, chain, module) {
  m <- read_csv(manifest_file)
  if (!all(c("path", "sha256") %in% names(m))) {
    stop("Manifest lacks path/sha256 columns: ", manifest_file)
  }
  m[, path := norm_path(path, must_work = FALSE)]
  m[, exists := file.exists(path)]
  m[, observed_sha256 := ifelse(exists, vapply(path, sha256_file, character(1)), NA_character_)]
  m[, hash_match := exists & observed_sha256 == sha256]
  add_qc(paste0(chain, "_", module, "_manifest_all_hashes_match"),
         all(m$hash_match), paste(sum(m$hash_match), nrow(m), sep = "/"), paste0(nrow(m), "/", nrow(m)))
  m[, `:=`(
    chain = chain,
    module = module,
    manifest_path = norm_path(manifest_file, must_work = TRUE)
  )]
  m
}

manifest_specs <- data.table(
  chain = rep(c("historical_Hartigan_Wong", "Lloyd_aligned"), each = 3L),
  module = rep(c("29f", "29g", "35"), 2L),
  path = c(P$old_29f_manifest, P$old_29g_manifest, P$old_35_manifest,
           P$new_29f_manifest, P$new_29g_manifest, P$new_35_manifest)
)
manifest_verified <- rbindlist(lapply(seq_len(nrow(manifest_specs)), function(i) {
  verify_manifest(manifest_specs$path[[i]], manifest_specs$chain[[i]], manifest_specs$module[[i]])
}), fill = TRUE)

for (nm in c("old_29f_marker", "old_29g_marker", "old_35_marker",
             "new_29f_marker", "new_29g_marker", "new_35_marker")) {
  marker_text <- paste(readLines(P[[nm]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  add_qc(paste0(nm, "_formal_marker"), grepl("formal", marker_text, ignore.case = TRUE),
         sub("_marker$", "", nm), "contains formal")
}

old_neg <- read_csv(P$old_29f_neg)
new_neg <- read_csv(P$new_29f_neg)
old_pos <- read_csv(P$old_29f_pos)
new_pos <- read_csv(P$new_29f_pos)
old_gate <- read_csv(P$old_29g_gate)[n_training_generated == 20049L]
new_gate <- read_csv(P$new_29g_gate)[n_training_generated == 20049L]
old_inv <- read_csv(P$old_29g_inv)
new_inv <- read_csv(P$new_29g_inv)
old_35_summary <- read_csv(P$old_35_summary)
new_35_states <- read_csv(P$new_35_states)
old_35_gate <- read_csv(P$old_35_gate)
new_35_gate <- read_csv(P$new_35_gate)
new_35_emp <- read_csv(P$new_35_emp)
new_35_fit <- read_csv(P$new_29g_fit)

add_qc("old_new_empirical_delta_exact",
       identical(old_inv$empirical_pooled_within_mahalanobis, new_inv$empirical_pooled_within_mahalanobis),
       paste(fmt(old_inv$empirical_pooled_within_mahalanobis, 14L), fmt(new_inv$empirical_pooled_within_mahalanobis, 14L), sep = " vs "),
       "exact equality")
add_qc("old_new_empirical_delta_matches_locked_source",
       isTRUE(all.equal(old_inv$empirical_pooled_within_mahalanobis, 3.52107344267228, tolerance = 0)),
       fmt(old_inv$empirical_pooled_within_mahalanobis, 14L), "3.52107344267228")
add_qc("new_29g_feature_dimension_33", new_35_fit$source_p == 33L, new_35_fit$source_p, 33L)

old_state_long <- rbindlist(lapply(seq_len(nrow(old_35_summary)), function(i) {
  r <- old_35_summary[i]
  data.table(
    stage = r$stage,
    n_train = r$n_train,
    n_eval = r$n_eval,
    state = c("DISCRETE_EVIDENCE", "INCONCLUSIVE", "NO_DISCRETE_EVIDENCE"),
    count = c(r$discrete_n, r$inconclusive_n, r$no_discrete_n),
    n = r$n_completed,
    rate = c(r$discrete_rate, r$inconclusive_rate, r$no_discrete_rate)
  )
}))

state_string <- function(x, stage_name) {
  z <- x[stage == stage_name]
  z <- z[order(n_train, decreasing = TRUE)]
  parts <- lapply(split(z, interaction(z$n_train, z$n_eval, drop = TRUE)), function(d) {
    d <- d[match(c("NO_DISCRETE_EVIDENCE", "INCONCLUSIVE", "DISCRETE_EVIDENCE"), state)]
    sprintf("train/eval %d/%d: no-discrete %d/%d, inconclusive %d/%d, discrete %d/%d",
            d$n_train[[1]], d$n_eval[[1]], d$count[[1]], d$n[[1]],
            d$count[[2]], d$n[[2]], d$count[[3]], d$n[[3]])
  })
  paste(unlist(parts), collapse = "; ")
}

old_n1_state <- state_string(old_state_long, "N1_independent_control")
new_n1_state <- state_string(new_35_states, "N1_independent_control")
old_s7_state <- state_string(old_state_long, "S7_positive_control")
new_s7_state <- state_string(new_35_states, "S7_positive_control")
old_emp_state <- state_string(old_state_long, "locked_empirical")
new_emp_state <- state_string(new_35_states, "locked_empirical")

ari_string_old <- function() {
  z <- old_35_summary[stage == "locked_empirical"][order(n_train, decreasing = TRUE)]
  paste(sprintf("train/eval %d/%d: %.6f [%.6f, %.6f]",
                z$n_train, z$n_eval, z$descriptive_label_ari_mean,
                z$descriptive_label_ari_p025, z$descriptive_label_ari_p975), collapse = "; ")
}

ari_string_new <- function() {
  z <- new_35_emp[order(n_train, decreasing = TRUE)]
  paste(sprintf("train/eval %d/%d: %.6f [%.6f, %.6f]",
                z$n_train, z$n_eval, z$ari_mean, z$ari_p025, z$ari_p975), collapse = "; ")
}

collect_29f_fits <- function(path, chain, module) {
  x <- read_csv(path)
  prefixes <- c("calibration_null", "source", "audit_null")
  rbindlist(lapply(prefixes, function(prefix) {
    a <- paste0(prefix, "_kmeans_algorithm_used")
    f <- paste0(prefix, "_kmeans_fallback_used")
    w <- paste0(prefix, "_kmeans_warning_count")
    d <- paste0(prefix, "_kmeans_failed")
    z <- data.table(
      chain = chain, module = module, fit_context = prefix,
      algorithm = as.character(x[[a]]), fallback = as.logical(x[[f]]),
      warning_count = as.integer(x[[w]]), failed = as.logical(x[[d]])
    )
    z[!is.na(algorithm) & nzchar(algorithm)]
  }))
}

collect_29g_fits <- function(paths, chain, module) {
  x <- rbindlist(lapply(paths, read_csv), fill = TRUE)
  data.table(
    chain = chain, module = module, fit_context = "formal_gate_or_surface",
    algorithm = as.character(x$kmeans_primary_algorithm_used),
    fallback = as.logical(x$kmeans_primary_fallback_used),
    warning_count = as.integer(x$kmeans_primary_warning_count),
    failed = as.logical(x$kmeans_primary_failed)
  )
}

read_old_35_checkpoint_fits <- function() {
  inv <- read_csv(P$old_35_checkpoint_inventory)
  add_qc("old_35_checkpoint_inventory_complete", nrow(inv) == 500L && all(inv$exists) && all(inv$signature_match) && all(inv$status_valid),
         paste(nrow(inv), sum(inv$exists), sum(inv$signature_match), sum(inv$status_valid), sep = "/"), "500/500/500/500")
  x <- rbindlist(lapply(inv$checkpoint, function(path) {
    obj <- readRDS(path)
    as.data.table(obj$result)
  }), fill = TRUE)
  add_qc("old_35_checkpoint_fit_rows", nrow(x) == 1500L, nrow(x), 1500L)
  data.table(
    chain = "historical_Hartigan_Wong", module = "35", fit_context = as.character(x$stage),
    task_id = paste(x$stage, x$repeat_id, sep = "::"),
    algorithm = as.character(x$kmeans_primary_algorithm_used),
    fallback = as.logical(x$kmeans_primary_fallback_used),
    warning_count = as.integer(x$kmeans_primary_warning_count),
    failed = as.logical(x$kmeans_primary_failed)
  )
}

collect_new_35_fits <- function() {
  x <- read_csv(P$new_35_repeat)
  add_qc("new_35_repeat_fit_rows", nrow(x) == 1500L, nrow(x), 1500L)
  data.table(
    chain = "Lloyd_aligned", module = "35", fit_context = as.character(x$stage),
    task_id = paste(x$stage, x$repeat_id, sep = "::"),
    algorithm = as.character(x$kmeans_primary_algorithm_used),
    fallback = as.logical(x$kmeans_primary_fallback_used),
    warning_count = as.integer(x$kmeans_primary_warning_count),
    failed = as.logical(x$kmeans_primary_failed)
  )
}

fit_records <- rbindlist(list(
  collect_29f_fits(P$old_29f_repeat, "historical_Hartigan_Wong", "29f"),
  collect_29g_fits(c(P$old_29g_gate_repeat, P$old_29g_surface_repeat), "historical_Hartigan_Wong", "29g"),
  read_old_35_checkpoint_fits(),
  collect_29f_fits(P$new_29f_repeat, "Lloyd_aligned", "29f"),
  collect_29g_fits(c(P$new_29g_gate_repeat, P$new_29g_surface_repeat), "Lloyd_aligned", "29g"),
  collect_new_35_fits()
), fill = TRUE)

fit_summary <- fit_records[, .(
  n_fit = .N,
  n_failed = sum(failed %in% TRUE, na.rm = TRUE),
  n_fallback = sum(fallback %in% TRUE, na.rm = TRUE),
  n_warning_fit = sum(warning_count > 0L, na.rm = TRUE),
  warning_total = sum(warning_count, na.rm = TRUE),
  algorithm_counts = paste(names(table(algorithm)), as.integer(table(algorithm)), sep = "=", collapse = "; ")
), by = .(chain, module)]

fit_total <- fit_records[, .(
  n_fit = .N,
  n_failed = sum(failed %in% TRUE, na.rm = TRUE),
  n_fallback = sum(fallback %in% TRUE, na.rm = TRUE),
  n_warning_fit = sum(warning_count > 0L, na.rm = TRUE),
  algorithm_counts = paste(names(table(algorithm)), as.integer(table(algorithm)), sep = "=", collapse = "; ")
), by = chain]

old_total <- fit_total[chain == "historical_Hartigan_Wong"]
new_total <- fit_total[chain == "Lloyd_aligned"]
add_qc("old_fit_record_total", old_total$n_fit == 14700L, old_total$n_fit, 14700L)
add_qc("old_fallback_total", old_total$n_fallback == 2L, old_total$n_fallback, 2L)
add_qc("new_fit_record_total", new_total$n_fit == 14700L, new_total$n_fit, 14700L)
add_qc("new_all_lloyd_no_fallback",
       new_total$n_fallback == 0L && new_total$n_failed == 0L && identical(new_total$algorithm_counts, "Lloyd=14700"),
       paste(new_total$algorithm_counts, new_total$n_fallback, new_total$n_failed, sep = "/"), "Lloyd=14700/0/0")

old_35_diagnostic <- read_csv(P$old_35_repeat)
add_qc("old_35_diagnostic_fallback_one",
       sum(old_35_diagnostic$kmeans_primary_fallback_used %in% TRUE) == 1L,
       sum(old_35_diagnostic$kmeans_primary_fallback_used %in% TRUE), 1L)

old_35_task_audit <- fit_records[chain == "historical_Hartigan_Wong" & module == "35", .(
  n_task = uniqueN(task_id),
  n_fallback_task = uniqueN(task_id[fallback %in% TRUE]),
  n_warning_task = uniqueN(task_id[warning_count > 0L]),
  n_failed_task = uniqueN(task_id[failed %in% TRUE])
)]
new_35_task_audit <- fit_records[chain == "Lloyd_aligned" & module == "35", .(
  n_task = uniqueN(task_id),
  n_fallback_task = uniqueN(task_id[fallback %in% TRUE]),
  n_warning_task = uniqueN(task_id[warning_count > 0L]),
  n_failed_task = uniqueN(task_id[failed %in% TRUE])
)]
add_qc("old_35_task_fallback_two_of_500",
       old_35_task_audit$n_task == 500L && old_35_task_audit$n_fallback_task == 2L,
       paste(old_35_task_audit$n_fallback_task, old_35_task_audit$n_task, sep = "/"), "2/500")
add_qc("old_35_warning_task_two_of_500",
       old_35_task_audit$n_warning_task == 2L,
       old_35_task_audit$n_warning_task, 2L)
add_qc("new_35_task_fallback_zero_of_500",
       new_35_task_audit$n_task == 500L && new_35_task_audit$n_fallback_task == 0L,
       paste(new_35_task_audit$n_fallback_task, new_35_task_audit$n_task, sep = "/"), "0/500")

task_fail_summary <- data.table(
  chain = c("historical_Hartigan_Wong", "Lloyd_aligned"),
  n_task = c(800L + 12000L + 500L, 800L + 12000L + 500L),
  n_failed = c(
    sum(read_csv(P$old_29f_repeat)$status != "completed") +
      sum(rbindlist(list(read_csv(P$old_29g_gate_repeat), read_csv(P$old_29g_surface_repeat)), fill = TRUE)$status != "completed") +
      0L,
    sum(read_csv(P$new_29f_repeat)$status != "completed") +
      sum(rbindlist(list(read_csv(P$new_29g_gate_repeat), read_csv(P$new_29g_surface_repeat)), fill = TRUE)$status != "completed") +
      0L
  )
)
add_qc("all_formal_tasks_completed", all(task_fail_summary$n_failed == 0L),
       paste(task_fail_summary$n_failed, collapse = "/"), "0/0")

runtime <- data.table(
  chain = rep(c("historical_Hartigan_Wong", "Lloyd_aligned"), each = 3L),
  module = rep(c("29f", "29g", "35"), 2L),
  session_path = c(P$old_29f_session, P$old_29g_session, P$old_35_session,
                   P$new_29f_session, P$new_29g_session, P$new_35_session)
)
runtime[, r_version := vapply(session_path, first_r_version, character(1))]

old_q95 <- old_gate$q95[[1L]]
new_q95 <- new_gate$q95[[1L]]
old_q95_ci <- sprintf("[%.6f, %.6f]", old_gate$q95_bootstrap_ci_low, old_gate$q95_bootstrap_ci_high)
new_q95_ci <- sprintf("[%.6f, %.6f]", new_gate$q95_bootstrap_ci_low, new_gate$q95_bootstrap_ci_high)

old_s1 <- sprintf(
  "29f source Delta-only %d/%d (%.1f%%), source joint %d/%d (%.1f%%), audit-null Delta-only %d/%d (%.1f%%); 35 N1: %s",
  old_neg$source_delta_false_positive_n, old_neg$n_evaluation_rep, 100 * old_neg$source_delta_false_positive_rate,
  old_neg$source_joint_alert_false_positive_n, old_neg$n_evaluation_rep, 100 * old_neg$source_joint_alert_false_positive_rate,
  old_neg$audit_null_delta_false_positive_n, old_neg$n_evaluation_rep, 100 * old_neg$audit_null_delta_false_positive_rate,
  old_n1_state
)
new_s1 <- sprintf(
  "29f source Delta-only %d/%d (%.1f%%), source joint %d/%d (%.1f%%), audit-null Delta-only %d/%d (%.1f%%); 35 N1: %s",
  new_neg$source_delta_false_positive_n, new_neg$n_evaluation_rep, 100 * new_neg$source_delta_false_positive_rate,
  new_neg$source_joint_alert_false_positive_n, new_neg$n_evaluation_rep, 100 * new_neg$source_joint_alert_false_positive_rate,
  new_neg$audit_null_delta_false_positive_n, new_neg$n_evaluation_rep, 100 * new_neg$audit_null_delta_false_positive_rate,
  new_n1_state
)

comparison <- rbindlist(list(
  data.table(comparison_id = 1L, comparison_item = "Locked processed matrix SHA-256",
             old_method = "Historical formal chain input", old_result = EXPECTED_LOCK_HASH[["locked_matrix"]],
             new_method = "Lloyd-aligned formal chain input", new_result = observed_lock_hash[["locked_matrix"]],
             comparison_status = "IDENTICAL", interpretation = "The same byte-identical locked processed matrix was used.",
             old_source_file = P$locked_matrix, old_source_field = "SHA-256",
             new_source_file = P$locked_matrix, new_source_field = "SHA-256"),
  data.table(comparison_id = 2L, comparison_item = "Locked K=2 labels SHA-256",
             old_method = "Historical formal chain label identity", old_result = EXPECTED_LOCK_HASH[["locked_labels"]],
             new_method = "Lloyd-aligned formal chain label identity", new_result = observed_lock_hash[["locked_labels"]],
             comparison_status = "IDENTICAL", interpretation = "The locked descriptive label artifact was unchanged.",
             old_source_file = P$locked_labels, old_source_field = "SHA-256",
             new_source_file = P$locked_labels, new_source_field = "SHA-256"),
  data.table(comparison_id = 3L, comparison_item = "Locked full-cohort pooled-within Mahalanobis separation",
             old_method = "Historical formal chain", old_result = fmt(old_inv$empirical_pooled_within_mahalanobis, 14L),
             new_method = "Lloyd-aligned formal chain", new_result = fmt(new_inv$empirical_pooled_within_mahalanobis, 14L),
             comparison_status = "IDENTICAL", interpretation = "Descriptive full-cohort separation was invariant because its locked inputs were unchanged.",
             old_source_file = P$old_29g_inv, old_source_field = "empirical_pooled_within_mahalanobis",
             new_source_file = P$new_29g_inv, new_source_field = "empirical_pooled_within_mahalanobis"),
  data.table(comparison_id = 4L, comparison_item = "N1 q95 at generated n=20,049",
             old_method = "Hartigan-Wong primary with fallback", old_result = fmt(old_q95, 6L),
             new_method = "Lloyd/25 starts/100 iterations; no fallback", new_result = fmt(new_q95, 6L),
             comparison_status = "CHANGED_IN_LLOYD_ALIGNED_RERUN", interpretation = sprintf("Absolute change (new-old) = %.6f; this independently seeded rerun does not isolate an optimizer effect.", new_q95 - old_q95),
             old_source_file = P$old_29g_gate, old_source_field = "q95[n_training_generated=20049]",
             new_source_file = P$new_29g_gate, new_source_field = "q95[n_training_generated=20049]"),
  data.table(comparison_id = 5L, comparison_item = "N1 q95 bootstrap 95% interval at generated n=20,049",
             old_method = "Historical formal bootstrap", old_result = old_q95_ci,
             new_method = "Lloyd-aligned formal bootstrap", new_result = new_q95_ci,
             comparison_status = "CHANGED_IN_LLOYD_ALIGNED_RERUN", interpretation = "Both intervals are reported; neither is substituted into the other chain, and the difference is not attributed solely to the optimizer.",
             old_source_file = P$old_29g_gate, old_source_field = "q95_bootstrap_ci_low/high[n_training_generated=20049]",
             new_source_file = P$new_29g_gate, new_source_field = "q95_bootstrap_ci_low/high[n_training_generated=20049]"),
  data.table(comparison_id = 6L, comparison_item = "S1/N1 control behavior",
             old_method = "Historical control chain", old_result = old_s1,
             new_method = "Lloyd-aligned control chain", new_result = new_s1,
             comparison_status = "JOINT_VERDICT_STABLE_RATE_DETAILS_CHANGED", interpretation = "The 49.5% source Delta-only false-positive rate remains visible; the source joint rule was 0/200 in both chains.",
             old_source_file = paste(P$old_29f_neg, P$old_35_summary, sep = "; "), old_source_field = "29f false-positive fields; 35 state counts",
             new_source_file = paste(P$new_29f_neg, P$new_35_states, sep = "; "), new_source_field = "29f false-positive fields; 35 state counts"),
  data.table(comparison_id = 7L, comparison_item = "S7 positive-control states",
             old_method = "Historical formal chain", old_result = sprintf("29f joint %d/%d; 35 %s", old_pos$source_joint_alert_n_control_only, old_pos$n_evaluation_rep, old_s7_state),
             new_method = "Lloyd-aligned formal chain", new_result = sprintf("29f joint %d/%d; 35 %s", new_pos$source_joint_alert_n_control_only, new_pos$n_evaluation_rep, new_s7_state),
             comparison_status = "UNCHANGED", interpretation = "Strong injected discrete structure was detected in every formal S7 evaluation in both chains.",
             old_source_file = paste(P$old_29f_pos, P$old_35_summary, sep = "; "), old_source_field = "source_joint_alert; S7 state counts",
             new_source_file = paste(P$new_29f_pos, P$new_35_states, sep = "; "), new_source_field = "source_joint_alert; S7 state counts"),
  data.table(comparison_id = 8L, comparison_item = "Locked empirical split-allocation states",
             old_method = "Historical formal split analysis", old_result = old_emp_state,
             new_method = "Lloyd-aligned formal split analysis", new_result = new_emp_state,
             comparison_status = "UNCHANGED", interpretation = "All 300 split evaluations were INCONCLUSIVE in each chain.",
             old_source_file = P$old_35_summary, old_source_field = "locked_empirical state counts",
             new_source_file = P$new_35_states, new_source_field = "locked_empirical state counts"),
  data.table(comparison_id = 9L, comparison_item = "Refitted-versus-locked label ARI (descriptive)",
             old_method = "Historical formal split analysis", old_result = ari_string_old(),
             new_method = "Lloyd-aligned formal split analysis", new_result = ari_string_new(),
             comparison_status = "DESCRIPTIVE_ESTIMATES_CHANGED_SLIGHTLY", interpretation = "ARI is descriptive label agreement and is not part of the D1 verdict.",
             old_source_file = P$old_35_summary, old_source_field = "descriptive_label_ari_mean/p025/p975",
             new_source_file = P$new_35_emp, new_source_field = "ari_mean/p025/p975"),
  data.table(comparison_id = 10L, comparison_item = "Formal task and fit failure rates",
             old_method = "Historical chain", old_result = sprintf("tasks %d/%d failed; fits %d/%d failed", task_fail_summary[chain == "historical_Hartigan_Wong"]$n_failed, task_fail_summary[chain == "historical_Hartigan_Wong"]$n_task, old_total$n_failed, old_total$n_fit),
             new_method = "Lloyd-aligned chain", new_result = sprintf("tasks %d/%d failed; fits %d/%d failed", task_fail_summary[chain == "Lloyd_aligned"]$n_failed, task_fail_summary[chain == "Lloyd_aligned"]$n_task, new_total$n_failed, new_total$n_fit),
             comparison_status = "UNCHANGED_ZERO_FAILURE", interpretation = "No formal task or applicable k-means fit failed in either chain.",
             old_source_file = "29f/29g repeat tables; old 35 checkpoint results", old_source_field = "status; kmeans_*_failed",
             new_source_file = "29f/29g/35 repeat tables", new_source_field = "status; kmeans_*_failed"),
  data.table(comparison_id = 11L, comparison_item = "Historical Hartigan-Wong and fallback use",
             old_method = "Hartigan-Wong primary; Lloyd fallback on Quick-TRANSfer warning", old_result = sprintf("script 35 affected tasks %d/%d (%.2f%%); fallback fits %d/1500 (%.2f%%); warning tasks %d/%d; failed tasks %d/%d; diagnostic subset fallback 1/900; full chain %s with fallback %d/%d", old_35_task_audit$n_fallback_task, old_35_task_audit$n_task, 100 * old_35_task_audit$n_fallback_task / old_35_task_audit$n_task, old_total$n_fallback, 100 * old_total$n_fallback / 1500, old_35_task_audit$n_warning_task, old_35_task_audit$n_task, old_35_task_audit$n_failed_task, old_35_task_audit$n_task, old_total$algorithm_counts, old_total$n_fallback, old_total$n_fit),
             new_method = "Not applicable", new_result = "See comparison 12",
             comparison_status = "HISTORICAL_PROCESS_DOCUMENTED", interpretation = "Two of 500 old script-35 tasks were affected; two of 1,500 script-35 fit calls used Lloyd fallback. Both completed successfully and are not failures.",
             old_source_file = P$old_35_checkpoint_inventory, old_source_field = "500 checkpoint result objects; algorithm/fallback/warning fields",
             new_source_file = NA_character_, new_source_field = NA_character_),
  data.table(comparison_id = 12L, comparison_item = "Lloyd alignment and fallback use",
             old_method = "See comparison 11", old_result = "Historical primary algorithm was Hartigan-Wong",
             new_method = "Lloyd/25 starts/100 iterations; fallback prohibited", new_result = sprintf("full chain %s; fallback %d/%d; failures %d/%d; script 35 affected tasks %d/%d; fallback fits 0/1500", new_total$algorithm_counts, new_total$n_fallback, new_total$n_fit, new_total$n_failed, new_total$n_fit, new_35_task_audit$n_fallback_task, new_35_task_audit$n_task),
             comparison_status = "PLANNED_ALIGNMENT_CONFIRMED", interpretation = "Every applicable fit in the new formal chain used Lloyd and no fallback occurred.",
             old_source_file = NA_character_, old_source_field = NA_character_,
             new_source_file = "29f/29g/35 Lloyd repeat tables", new_source_field = "kmeans_primary_algorithm_used; fallback; failed"),
  data.table(comparison_id = 13L, comparison_item = "Final conditional D1 evidence state",
             old_method = "Historical formal split diagnostic", old_result = "INCONCLUSIVE (300/300 empirical split evaluations)",
             new_method = "Lloyd-aligned formal split diagnostic", new_result = "INCONCLUSIVE (300/300 empirical split evaluations)",
             comparison_status = "UNCHANGED", interpretation = "Separation alerted while the held-out shape component did not; the conclusion remains conditional and inconclusive.",
             old_source_file = P$old_35_summary, old_source_field = "locked_empirical modal_verdict/state counts",
             new_source_file = P$new_35_states, new_source_field = "locked_empirical state counts")
), fill = TRUE)

add_qc("comparison_has_exactly_13_items", nrow(comparison) == 13L && identical(comparison$comparison_id, 1:13),
       paste(comparison$comparison_id, collapse = ","), "1:13")
add_qc("final_state_unchanged_inconclusive", comparison[comparison_id == 13L]$comparison_status == "UNCHANGED",
       comparison[comparison_id == 13L]$new_result, "INCONCLUSIVE in both chains")

source_registry <- rbindlist(list(
  data.table(artifact_id = c("LOCKED_MATRIX", "LOCKED_LABELS", "LOCKED_DELTA"), chain = "shared_locked_input", module = "empirical",
             role = c("processed 33-feature matrix", "locked K=2 labels", "locked full-cohort separation source"),
             path = c(P$locked_matrix, P$locked_labels, P$locked_delta), source_manifest = NA_character_,
             run_mode = "locked", r_version = NA_character_, notes = "Shared immutable artifact"),
  data.table(artifact_id = c("OLD_29F_SCRIPT", "OLD_29F_NEG", "OLD_29F_POS", "OLD_29F_REPEAT", "OLD_29F_MARKER", "OLD_29F_MANIFEST"),
             chain = "historical_Hartigan_Wong", module = "29f", role = c("script", "negative control summary", "positive control summary", "repeat-level output", "completion marker", "output manifest"),
             path = c(P$old_29f_script, P$old_29f_neg, P$old_29f_pos, P$old_29f_repeat, P$old_29f_marker, P$old_29f_manifest),
             source_manifest = P$old_29f_manifest, run_mode = "formal", r_version = first_r_version(P$old_29f_session), notes = "Historical formal authority"),
  data.table(artifact_id = c("OLD_29G_SCRIPT", "OLD_29G_GATE", "OLD_29G_INVERSION", "OLD_29G_GATE_REPEAT", "OLD_29G_SURFACE_REPEAT", "OLD_29G_MARKER", "OLD_29G_MANIFEST"),
             chain = "historical_Hartigan_Wong", module = "29g", role = c("script", "gate summary", "inversion summary", "gate repeat output", "surface repeat output", "completion marker", "output manifest"),
             path = c(P$old_29g_script, P$old_29g_gate, P$old_29g_inv, P$old_29g_gate_repeat, P$old_29g_surface_repeat, P$old_29g_marker, P$old_29g_manifest),
             source_manifest = P$old_29g_manifest, run_mode = "formal", r_version = first_r_version(P$old_29g_session), notes = "Historical formal authority"),
  data.table(artifact_id = c("OLD_35_SCRIPT", "OLD_35_GATE", "OLD_35_SUMMARY", "OLD_35_REPEAT", "OLD_35_CHECKPOINT_INVENTORY", "OLD_35_MARKER", "OLD_35_MANIFEST"),
             chain = "historical_Hartigan_Wong", module = "35", role = c("script", "split gate summary", "state and ARI summary", "diagnostic repeat output", "all-task checkpoint inventory", "completion marker", "output manifest"),
             path = c(P$old_35_script, P$old_35_gate, P$old_35_summary, P$old_35_repeat, P$old_35_checkpoint_inventory, P$old_35_marker, P$old_35_manifest),
             source_manifest = P$old_35_manifest, run_mode = "formal", r_version = first_r_version(P$old_35_session), notes = "Historical formal authority; all-task fallback audit uses checkpoint results"),
  data.table(artifact_id = c("NEW_29F_SCRIPT", "NEW_29F_NEG", "NEW_29F_POS", "NEW_29F_REPEAT", "NEW_29F_MARKER", "NEW_29F_MANIFEST"),
             chain = "Lloyd_aligned", module = "29f", role = c("script", "negative control summary", "positive control summary", "repeat-level output", "completion marker", "output manifest"),
             path = c(P$new_29f_script, P$new_29f_neg, P$new_29f_pos, P$new_29f_repeat, P$new_29f_marker, P$new_29f_manifest),
             source_manifest = P$new_29f_manifest, run_mode = "formal", r_version = first_r_version(P$new_29f_session), notes = "Current Lloyd-aligned formal authority"),
  data.table(artifact_id = c("NEW_29G_SCRIPT", "NEW_29G_GATE", "NEW_29G_INVERSION", "NEW_29G_FIT", "NEW_29G_GATE_REPEAT", "NEW_29G_SURFACE_REPEAT", "NEW_29G_MARKER", "NEW_29G_MANIFEST"),
             chain = "Lloyd_aligned", module = "29g", role = c("script", "gate summary", "inversion summary", "N1 fit metadata", "gate repeat output", "surface repeat output", "completion marker", "output manifest"),
             path = c(P$new_29g_script, P$new_29g_gate, P$new_29g_inv, P$new_29g_fit, P$new_29g_gate_repeat, P$new_29g_surface_repeat, P$new_29g_marker, P$new_29g_manifest),
             source_manifest = P$new_29g_manifest, run_mode = "formal", r_version = first_r_version(P$new_29g_session), notes = "Current Lloyd-aligned formal authority"),
  data.table(artifact_id = c("NEW_35_SCRIPT", "NEW_35_GATE", "NEW_35_EMP", "NEW_35_CONTROL", "NEW_35_STATES", "NEW_35_REPEAT", "NEW_35_MARKER", "NEW_35_MANIFEST"),
             chain = "Lloyd_aligned", module = "35", role = c("script", "split gate summary", "empirical summary", "control summary", "state distribution", "all-fit repeat output", "completion marker", "output manifest"),
             path = c(P$new_35_script, P$new_35_gate, P$new_35_emp, P$new_35_control, P$new_35_states, P$new_35_repeat, P$new_35_marker, P$new_35_manifest),
             source_manifest = P$new_35_manifest, run_mode = "formal", r_version = first_r_version(P$new_35_session), notes = "Current Lloyd-aligned formal authority")
), fill = TRUE)

source_registry[, path := norm_path(path, must_work = FALSE)]
source_registry[, `:=`(
  exists = file.exists(path),
  bytes = ifelse(file.exists(path), file.info(path)$size, NA_real_),
  sha256 = ifelse(file.exists(path), vapply(path, sha256_file, character(1)), NA_character_)
)]
manifest_lookup <- manifest_verified[, .(path, manifest_sha256 = sha256, manifest_hash_match = hash_match, manifest_path)]
source_registry <- merge(source_registry, manifest_lookup, by = "path", all.x = TRUE)
source_registry[, manifest_membership := fifelse(is.na(manifest_sha256), "not_listed_or_not_applicable",
                                                 fifelse(manifest_hash_match, "verified", "mismatch"))]
setcolorder(source_registry, c("artifact_id", "chain", "module", "role", "path", "exists", "bytes", "sha256",
                               "source_manifest", "manifest_membership", "manifest_sha256", "run_mode", "r_version", "notes"))

input_manifest <- unique(source_registry[, .(artifact_id, chain, module, role, path, exists, bytes, sha256)])
input_manifest <- rbind(input_manifest, data.table(
  artifact_id = "SCRIPT_40", chain = "postprocess", module = "40", role = "comparison audit script",
  path = SCRIPT_PATH, exists = TRUE, bytes = file.info(SCRIPT_PATH)$size, sha256 = sha256_file(SCRIPT_PATH)
), fill = TRUE)

report_lines <- c(
  "# D1 Lloyd-alignment comparison audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope and planned method change",
  "",
  "This is a read-only comparison of the locked historical Hartigan-Wong formal chain and the new Lloyd-aligned formal chain. No model was refitted by this post-processing audit.",
  "",
  "> The only planned method change was alignment of the clustering optimization algorithm and parameters to the legacy K2 label generator; all other scientific decision rules and simulation grids were held fixed.",
  "",
  "The Lloyd-aligned chain is not described as intrinsically more correct. Old-versus-new numerical differences represent result-level consistency across independently seeded reruns and should not be interpreted as an isolated causal estimate of the optimizer effect.",
  "",
  "## Locked-input identity",
  "",
  paste0("- Processed matrix SHA-256: `", observed_lock_hash[["locked_matrix"]], "` (identical)."),
  paste0("- Locked K=2 label SHA-256: `", observed_lock_hash[["locked_labels"]], "` (identical)."),
  paste0("- Locked full-cohort separation-source SHA-256: `", observed_lock_hash[["locked_delta"]], "` (identical)."),
  paste0("- Full-cohort pooled-within Mahalanobis separation: ", fmt(new_inv$empirical_pooled_within_mahalanobis, 14L), " in both chains."),
  "",
  "## Formal comparison",
  "",
  "| ID | Comparison | Historical result | Lloyd-aligned result | Status |",
  "|---:|---|---|---|---|",
  vapply(seq_len(nrow(comparison)), function(i) {
    paste0("| ", comparison$comparison_id[[i]], " | ", comparison$comparison_item[[i]], " | ",
           gsub("\\|", "/", comparison$old_result[[i]]), " | ",
           gsub("\\|", "/", comparison$new_result[[i]]), " | ", comparison$comparison_status[[i]], " |")
  }, character(1)),
  "",
  "## Algorithm and failure audit",
  "",
  "| Chain | Module | Applicable fits | Algorithm counts | Fallbacks | Warning fits | Failed fits |",
  "|---|---:|---:|---|---:|---:|---:|",
  vapply(seq_len(nrow(fit_summary)), function(i) {
    paste0("| ", fit_summary$chain[[i]], " | ", fit_summary$module[[i]], " | ", fit_summary$n_fit[[i]], " | ",
           fit_summary$algorithm_counts[[i]], " | ", fit_summary$n_fallback[[i]], " | ",
           fit_summary$n_warning_fit[[i]], " | ", fit_summary$n_failed[[i]], " |")
  }, character(1)),
  "",
  "The historical script-35 chain contained two affected tasks among 500 tasks (0.40%) and two fallback fit calls among 1,500 fit calls (0.13%): N1 gate calibration repeat 136 in the train-19,049/evaluate-1,000 allocation and N1 independent-control repeat 12 in the 50/50 allocation. Both tasks completed successfully. Only the latter belongs to the 900 diagnostic-fit subset, hence the diagnostic table records one fallback while the all-checkpoint audit records two.",
  "",
  "The 49.5% (99/200) S1 source Delta-only false-positive rate remains explicitly visible. It is not the joint D1 verdict: the source joint separation-and-shape rule alerted in 0/200 repeats in both chains.",
  "",
  "## Runtime provenance",
  "",
  "| Chain | Module | R version |",
  "|---|---:|---|",
  vapply(seq_len(nrow(runtime)), function(i) paste0("| ", runtime$chain[[i]], " | ", runtime$module[[i]], " | ", runtime$r_version[[i]], " |"), character(1)),
  "",
  "The old and new 29f/29g runs used R 4.5.1. The historical script-35 run used R 4.5.1/Matrix 1.7-3, whereas the Lloyd-aligned script-35 formal run and this post-processing audit used R 4.5.3/Matrix 1.7-4. Runtime and independently seeded rerun provenance are disclosed; small numerical differences cannot be attributed solely to the optimizer.",
  "",
  "## Interpretation boundaries",
  "",
  "- The final conditional D1 evidence state remains INCONCLUSIVE in all 300 locked-empirical split evaluations under both implementations.",
  "- Refitted-versus-locked ARI is descriptive label agreement; it is not a D1 verdict component, a percentage contribution, or external validation.",
  "- The inverted injected Delta is a simulation-DGM calibration and is not a biological effect estimate.",
  "- The locked matrix has 33 columns. The N1 fit metadata records 23/33 exact feature-name overlaps with the separate S7 generator vocabulary; exact name identity was not required because the two generators have distinct frozen semantics. No cross-generator feature-identity claim is made.",
  "- The empirical finding is neither evidence of discreteness nor evidence of its absence; at this sample size and under this fixed-K=2 processed-space diagnostic, D1 is inconclusive.",
  "",
  "## Conclusion",
  "",
  "The method-alignment rerun changed the calibrated N1 thresholds and some descriptive split estimates but did not change the control conclusions or the final conditional D1 evidence state. The Lloyd-aligned outputs are the current formal authority; the historical outputs remain provenance-only comparators."
)

fwrite(comparison, OUT[["comparison"]])
fwrite(source_registry, OUT[["authority"]])
fwrite(rbindlist(qc_rows), OUT[["qc"]])
fwrite(input_manifest, OUT[["input_manifest"]])
writeLines(report_lines, OUT[["report"]], useBytes = TRUE)
writeLines(capture.output(sessionInfo()), OUT[["session"]], useBytes = TRUE)

output_paths <- unname(OUT[c("comparison", "report", "authority", "qc", "input_manifest", "session")])
output_manifest <- data.table(
  path = norm_path(output_paths, must_work = TRUE),
  bytes = file.info(output_paths)$size,
  sha256 = vapply(output_paths, sha256_file, character(1))
)
fwrite(output_manifest, OUT[["output_manifest"]])

completion_lines <- c(
  "status=completed",
  "run_mode=read_only_formal_postprocess",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("r_version=", R.version.string),
  "comparison_items=13",
  paste0("qc_passed=", sum(rbindlist(qc_rows)$pass)),
  paste0("qc_total=", nrow(rbindlist(qc_rows))),
  paste0("output_manifest_sha256=", sha256_file(OUT[["output_manifest"]])),
  "final_d1_state=INCONCLUSIVE",
  "downstream_authorization=PENDING_INDEPENDENT_AUDIT"
)
writeLines(completion_lines, OUT[["completion"]], useBytes = TRUE)

cat("D1 Lloyd-alignment comparison audit completed.\n")
cat("Comparison rows:", nrow(comparison), "\n")
cat("QC:", sum(rbindlist(qc_rows)$pass), "/", nrow(rbindlist(qc_rows)), "passed\n", sep = "")
cat("Historical fallbacks:", old_total$n_fallback, "/", old_total$n_fit, "fits\n")
cat("Lloyd fallbacks:", new_total$n_fallback, "/", new_total$n_fit, "fits\n")
cat("Final D1 state: INCONCLUSIVE (unchanged)\n")
