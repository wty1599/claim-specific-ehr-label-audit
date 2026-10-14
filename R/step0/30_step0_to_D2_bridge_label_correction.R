#!/usr/bin/env Rscript

## WP10-D2 deterministic method-label correction.
## No estimate, confidence interval, cohort, label, clustering, or model is changed.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
AUDIT_ROOT <- file.path(
  PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"
)
SCRIPT_PATH <- file.path(
  AUDIT_ROOT, "scripts", "36b_domain1_to_domain2_bridge_method_label_correction.R"
)
RUN_ROOT <- file.path(
  AUDIT_ROOT, "wp10_domain2_domain4_bridge_20260718", "domain2"
)
TABLE_DIR <- file.path(RUN_ROOT, "tables")
LOG_DIR <- file.path(RUN_ROOT, "logs")
PROVENANCE_DIR <- file.path(RUN_ROOT, "provenance")

INPUT <- c(
  original_bridge = file.path(TABLE_DIR, "Table_D1_D2_bridge.csv"),
  original_completion = file.path(RUN_ROOT, "36_WP10_D2_run_completed.ok"),
  mice_pooled = file.path(PROJECT_ROOT, "output", "qc", "mice_b2_rubin_only.csv"),
  mice_source_script = file.path(PROJECT_ROOT, "r", "23b_mice_b2_rubin_only.R")
)
EXPECTED_SHA256 <- c(
  original_bridge = "334AC7A0ECFB6BDF5E2855FEBB41F3C314218D4F50D5183DCDAE8382785A416F",
  original_completion = "4633A3BB09C3672D34B97E8AB75FEC6D5C514A55C2C3A02B0280EB383CFBDDB3",
  mice_pooled = "EDCD4A7A9490BE05526F8FD15476CD3718479F0916DB0F2DFD3FA91103F4B488",
  mice_source_script = "4B33ECFB7BF9427F084CF3F07DD92A76B3088455C2A49988C385042C04262DA3"
)

OUTPUT_TABLE <- file.path(TABLE_DIR, "Table_D1_D2_bridge_method_corrected_v2.csv")
QC_FILE <- file.path(LOG_DIR, "36b_WP10_D2_method_correction_QC.csv")
LOG_FILE <- file.path(LOG_DIR, "36b_WP10_D2_method_correction_log.md")
INPUT_MANIFEST <- file.path(
  PROVENANCE_DIR, "36b_WP10_D2_method_correction_input_sha256_manifest.csv"
)
OUTPUT_MANIFEST <- file.path(
  PROVENANCE_DIR, "36b_WP10_D2_method_correction_output_sha256_manifest.csv"
)
COMPLETION <- file.path(RUN_ROOT, "36b_WP10_D2_method_correction_completed.ok")

for (p in c(SCRIPT_PATH, INPUT)) {
  if (!file.exists(p)) stop("Missing required WP10-D2 correction input: ", p, call. = FALSE)
}
if (any(file.exists(c(
  OUTPUT_TABLE, QC_FILE, LOG_FILE, INPUT_MANIFEST, OUTPUT_MANIFEST, COMPLETION
)))) {
  stop("Versioned WP10-D2 correction outputs already exist; refusing overwrite.", call. = FALSE)
}

sha256_file <- function(path) {
  toupper(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(x, tmp, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path, call. = FALSE)
}

atomic_write_lines <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic text write failed: ", path, call. = FALSE)
}

input_manifest <- rbindlist(lapply(names(INPUT), function(nm) {
  p <- INPUT[[nm]]
  data.table(
    input_name = nm,
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(p),
    expected_sha256 = EXPECTED_SHA256[[nm]],
    hash_match = sha256_file(p) == EXPECTED_SHA256[[nm]],
    bytes = file.info(p)$size
  )
}))
if (!all(input_manifest$hash_match)) {
  stop("At least one WP10-D2 correction input hash differs from the audited version.",
       call. = FALSE)
}
atomic_fwrite(input_manifest, INPUT_MANIFEST)

original <- fread(INPUT[["original_bridge"]])
mice <- fread(INPUT[["mice_pooled"]])
source_text <- paste(
  readLines(INPUT[["mice_source_script"]], warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)

required_bridge <- c(
  "bridge_section", "evidence_id", "analysis", "estimate", "ci_low", "ci_high",
  "component_auc", "reference_auc", "n", "folds_or_imputations"
)
required_mice <- c(
  "outcome", "m", "decisive_dAUC_feat_plus_label", "dAUC_lo", "dAUC_hi",
  "within_imp_var_W", "between_imp_var_B", "total_var_T"
)
if (length(setdiff(required_bridge, names(original)))) {
  stop("Original bridge lacks required columns.", call. = FALSE)
}
if (length(setdiff(required_mice, names(mice)))) {
  stop("Authoritative MICE table lacks variance fields.", call. = FALSE)
}

mice_rows <- grepl("^D2_mice_", original$evidence_id)
old_label <- "Raw-EN + K2 minus Raw-EN; paired DeLong within imputation and Rubin pooling"
new_label <- paste(
  "Raw-EN + K2 minus Raw-EN; 500 nonparametric bootstrap resamples within",
  "each imputation and Rubin-style variance pooling"
)
variance_label <- paste(
  "Within-imputation variance: nonparametric bootstrap (B=500);",
  "total variance: W + (1 + 1/m)B; interval: estimate +/- 1.96 SE"
)

if (sum(mice_rows) != 2L || !all(original$analysis[mice_rows] == old_label)) {
  stop("The expected two outdated MICE method labels were not found.", call. = FALSE)
}
if (!all(mice$m == 5L) || !grepl("BOOT_B\\s*<-\\s*500L", source_text, perl = TRUE) ||
    !grepl("Tt\\s*<-\\s*W\\s*\\+\\s*\\(1\\s*\\+\\s*1\\s*/\\s*MICE_M_USE\\)\\s*\\*\\s*Bv",
           source_text, perl = TRUE)) {
  stop("Could not verify B=500 and Rubin-style total-variance formula upstream.",
       call. = FALSE)
}

corrected <- copy(original)
corrected[mice_rows, analysis := new_label]
corrected[, variance_method := NA_character_]
corrected[mice_rows, variance_method := variance_label]
corrected[, correction_status := fifelse(
  mice_rows,
  "METHOD_LABEL_CORRECTED_VALUES_UNCHANGED",
  "UNCHANGED_FROM_ORIGINAL_BRIDGE"
)]

unchanged_columns <- setdiff(names(original), "analysis")
values_unchanged <- isTRUE(all.equal(
  original[, ..unchanged_columns],
  corrected[, ..unchanged_columns],
  check.attributes = TRUE
))

qc <- rbindlist(list(
  data.table(
    check = "locked_input_hashes_match", pass = all(input_manifest$hash_match),
    observed = sum(input_manifest$hash_match), expected = nrow(input_manifest)
  ),
  data.table(
    check = "exactly_two_mice_rows_corrected", pass = sum(mice_rows) == 2L,
    observed = sum(mice_rows), expected = 2L
  ),
  data.table(
    check = "upstream_bootstrap_and_pooling_verified",
    pass = all(mice$m == 5L) && all(is.finite(mice$within_imp_var_W)) &&
      all(is.finite(mice$between_imp_var_B)) && all(is.finite(mice$total_var_T)),
    observed = paste(unique(mice$m), nrow(mice), sep = "/"), expected = "5/2"
  ),
  data.table(
    check = "all_nonmethod_values_unchanged", pass = values_unchanged,
    observed = values_unchanged, expected = TRUE
  ),
  data.table(
    check = "outdated_DeLong_label_removed",
    pass = !any(grepl("DeLong", corrected$analysis, fixed = TRUE)),
    observed = sum(grepl("DeLong", corrected$analysis, fixed = TRUE)), expected = 0L
  ),
  data.table(
    check = "corrected_method_explicit",
    pass = all(grepl("500 nonparametric bootstrap", corrected$analysis[mice_rows], fixed = TRUE)) &&
      all(grepl("Rubin-style variance pooling", corrected$analysis[mice_rows], fixed = TRUE)),
    observed = sum(grepl("500 nonparametric bootstrap", corrected$analysis[mice_rows],
                         fixed = TRUE)),
    expected = 2L
  ),
  data.table(
    check = "no_model_or_estimate_recalculation",
    pass = TRUE, observed = "deterministic text correction only",
    expected = "no refit or numerical modification"
  )
), fill = TRUE)
if (!all(qc$pass)) stop("WP10-D2 method-label correction failed QC.", call. = FALSE)

atomic_fwrite(corrected, OUTPUT_TABLE)
atomic_fwrite(qc, QC_FILE)
atomic_write_lines(c(
  "# WP10-D2 method-label correction",
  "",
  "- The original bridge incorrectly described the MICE-pooled within-imputation",
  "  variance as paired DeLong.",
  "- The authoritative upstream script uses 500 nonparametric bootstrap resamples",
  "  within each imputation, then W + (1 + 1/m)B Rubin-style variance pooling.",
  "- No estimate, confidence interval, cohort, label, clustering result, or model",
  "  was changed or recomputed.",
  "- The original bridge is retained as superseded provenance.",
  "",
  paste0("Original bridge SHA-256: `", EXPECTED_SHA256[["original_bridge"]], "`"),
  paste0("Authoritative upstream R SHA-256: `", EXPECTED_SHA256[["mice_source_script"]], "`")
), LOG_FILE)

output_files <- c(OUTPUT_TABLE, QC_FILE, LOG_FILE, INPUT_MANIFEST)
output_manifest <- data.table(
  path = normalizePath(output_files, winslash = "/", mustWork = TRUE),
  sha256 = vapply(output_files, sha256_file, character(1)),
  bytes = file.info(output_files)$size
)
atomic_fwrite(output_manifest, OUTPUT_MANIFEST)

manifest_check <- fread(OUTPUT_MANIFEST)
manifest_check[, observed_sha256 := vapply(path, sha256_file, character(1))]
if (!all(toupper(manifest_check$sha256) == toupper(manifest_check$observed_sha256))) {
  stop("WP10-D2 corrected output manifest failed recheck.", call. = FALSE)
}

atomic_write_lines(c(
  "RUN COMPLETED",
  "work_package=WP10_Domain2_method_label_correction",
  "statistics_or_models_changed=FALSE",
  "supersedes_table=Table_D1_D2_bridge.csv",
  "authoritative_table=Table_D1_D2_bridge_method_corrected_v2.csv",
  "within_imputation_variance=nonparametric_bootstrap_B500",
  "pooling=W_plus_1_plus_1_over_m_times_B",
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("output_manifest_sha256=", sha256_file(OUTPUT_MANIFEST))
), COMPLETION)

cat("WP10-D2 method label corrected; all numerical fields retained.\n")
print(corrected)
