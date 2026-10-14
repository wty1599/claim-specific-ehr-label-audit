## Domain 1 targeted post-processing: KDIGO-stratified feature profiles.
##
## This script does not refit clustering or modify locked source results.
## It reports median (IQR) values from MICE completed dataset #1, which is
## the dataset used to generate the locked plain-k-means K2 labels.

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
})

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
ANALYSIS_ROOT <- file.path(
  PROJECT_ROOT,
  "domain1_renal_identity_and_method_closure_20260717"
)
TABLE_DIR <- file.path(ANALYSIS_ROOT, "tables")
LOG_DIR <- file.path(ANALYSIS_ROOT, "logs")
PROVENANCE_DIR <- file.path(ANALYSIS_ROOT, "provenance")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROVENANCE_DIR, recursive = TRUE, showWarnings = FALSE)

INPUT <- c(
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  mice = file.path(PROJECT_ROOT, "output", "model", "mice_primary.rds"),
  standardized_matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  baseline = file.path(PROJECT_ROOT, "data", "baseline_covars.csv")
)

if (any(!file.exists(INPUT))) {
  stop("Missing locked input(s): ", paste(names(INPUT)[!file.exists(INPUT)], collapse = ", "))
}

PRIMARY_IMPUTATION <- 1L
EXPECTED_N <- 20049L
EXPECTED_C1 <- 3992L
EXPECTED_COUNTS <- c(
  stage0_C1 = 969L,
  stage0_C2 = 8087L,
  stage23_C1 = 2244L,
  all_C1 = 3992L
)

FEATURES <- c(
  "creatinine_max",
  "bun_max",
  "lactate_max",
  "aniongap_max",
  "bicarbonate_min",
  "ph_min",
  "urine_output_24h_ml"
)

FEATURE_LABEL <- c(
  creatinine_max = "Maximum creatinine",
  bun_max = "Maximum blood urea nitrogen",
  lactate_max = "Maximum lactate",
  aniongap_max = "Maximum anion gap",
  bicarbonate_min = "Minimum bicarbonate",
  ph_min = "Minimum pH",
  urine_output_24h_ml = "Urine output, 0-24 h"
)

FEATURE_UNIT <- c(
  creatinine_max = "mg/dL",
  bun_max = "mg/dL",
  lactate_max = "mmol/L",
  aniongap_max = "mmol/L",
  bicarbonate_min = "mmol/L",
  ph_min = "unitless",
  urine_output_24h_ml = "mL"
)

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

format_number <- function(x, variable) {
  digits <- switch(
    variable,
    creatinine_max = 2L,
    bun_max = 1L,
    lactate_max = 1L,
    aniongap_max = 1L,
    bicarbonate_min = 1L,
    ph_min = 3L,
    urine_output_24h_ml = 1L,
    2L
  )
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

labels <- as.data.table(readRDS(INPUT[["labels"]]))
imp <- readRDS(INPUT[["mice"]])
X_locked <- as.data.table(readRDS(INPUT[["standardized_matrix"]]))
baseline <- fread(INPUT[["baseline"]], showProgress = FALSE)

if (!inherits(imp, "mids")) stop("mice_primary.rds is not a mids object")
if (imp$m < PRIMARY_IMPUTATION) stop("MICE completed dataset #1 is unavailable")
if (!all(c("stay_id", "cluster_k2") %in% names(labels))) stop("Malformed label file")
if (!all(c("stay_id", "aki_stage_0_24h") %in% names(baseline))) stop("Malformed baseline file")
if (anyDuplicated(labels$stay_id) || anyDuplicated(baseline$stay_id)) stop("Duplicate stay_id")
if (nrow(labels) != EXPECTED_N) stop("Locked label denominator changed")
if (sum(labels$cluster_k2 == 1L) != EXPECTED_C1) stop("Locked C1 count changed")

completed <- as.data.table(mice::complete(imp, PRIMARY_IMPUTATION))
if (nrow(completed) != EXPECTED_N) stop("Completed MICE denominator changed")
if (length(setdiff(FEATURES, names(completed)))) {
  stop("MICE object lacks required features: ", paste(setdiff(FEATURES, names(completed)), collapse = ", "))
}

locked_features <- setdiff(names(X_locked), "stay_id")
if (!identical(locked_features, names(imp$data))) stop("Locked feature order differs from MICE data")
if (!identical(as.integer(labels$stay_id), as.integer(X_locked$stay_id))) {
  stop("Locked labels and standardized matrix are not row-aligned")
}
X_rebuilt <- scale(as.matrix(completed[, ..locked_features]))
max_scaled_difference <- max(abs(X_rebuilt - as.matrix(X_locked[, ..locked_features])))
if (!is.finite(max_scaled_difference) || max_scaled_difference > 1e-8) {
  stop("MICE #1 does not reproduce the locked standardized matrix")
}

stage <- baseline$aki_stage_0_24h[match(labels$stay_id, baseline$stay_id)]
if (anyNA(stage) || !all(stage %in% 0:3)) stop("Invalid or missing 0-24 h KDIGO stage")

analysis_data <- cbind(
  labels[, .(stay_id, cluster_k2)],
  data.table(kdigo_stage_0_24h = as.integer(stage)),
  completed[, ..FEATURES]
)

GROUPS <- list(
  stage0_C1 = which(analysis_data$kdigo_stage_0_24h == 0L & analysis_data$cluster_k2 == 1L),
  stage0_C2 = which(analysis_data$kdigo_stage_0_24h == 0L & analysis_data$cluster_k2 == 2L),
  stage23_C1 = which(analysis_data$kdigo_stage_0_24h >= 2L & analysis_data$cluster_k2 == 1L),
  all_C1 = which(analysis_data$cluster_k2 == 1L)
)

GROUP_LABEL <- c(
  stage0_C1 = "KDIGO stage 0 and C1",
  stage0_C2 = "KDIGO stage 0 and C2",
  stage23_C1 = "KDIGO stage 2/3 and C1",
  all_C1 = "All C1"
)

observed_counts <- vapply(GROUPS, length, integer(1))
if (!identical(observed_counts, EXPECTED_COUNTS)) {
  stop(
    "Prespecified group counts changed: ",
    paste(names(observed_counts), observed_counts, sep = "=", collapse = "; ")
  )
}

profile_long <- rbindlist(lapply(names(GROUPS), function(group_id) {
  z <- analysis_data[GROUPS[[group_id]]]
  rbindlist(lapply(FEATURES, function(variable) {
    q <- quantile(z[[variable]], probs = c(0.25, 0.50, 0.75), names = FALSE, type = 7)
    data.table(
      group_id = group_id,
      group = unname(GROUP_LABEL[[group_id]]),
      n = nrow(z),
      variable = variable,
      measure = unname(FEATURE_LABEL[[variable]]),
      unit = unname(FEATURE_UNIT[[variable]]),
      q1 = q[[1L]],
      median = q[[2L]],
      q3 = q[[3L]],
      median_iqr = paste0(
        format_number(q[[2L]], variable), " (",
        format_number(q[[1L]], variable), "-",
        format_number(q[[3L]], variable), ")"
      ),
      source = "MICE completed dataset #1 used for locked plain-k-means labels"
    )
  }))
}))

profile_long[, group_id := factor(group_id, levels = names(GROUPS))]
profile_long[, feature_order := match(variable, FEATURES)]
setorder(profile_long, group_id, feature_order)
profile_long[, group_id := as.character(group_id)]
profile_long[, feature_order := NULL]

profile_wide <- dcast(
  profile_long,
  variable + measure + unit ~ group_id,
  value.var = "median_iqr"
)
setcolorder(profile_wide, c("variable", "measure", "unit", names(GROUPS)))
profile_wide[, feature_order := match(variable, FEATURES)]
setorder(profile_wide, feature_order)
profile_wide[, feature_order := NULL]
setnames(
  profile_wide,
  old = names(GROUPS),
  new = paste0(unname(GROUP_LABEL[names(GROUPS)]), " (n=", observed_counts, ")")
)

qc <- data.table(
  check = c(
    "locked_total_n",
    "locked_C1_n",
    paste0("group_n_", names(GROUPS)),
    "mice_primary_imputation",
    "max_abs_difference_rebuilt_vs_locked_scaled_matrix",
    "missing_profile_values"
  ),
  observed = c(
    nrow(labels),
    sum(labels$cluster_k2 == 1L),
    observed_counts,
    PRIMARY_IMPUTATION,
    max_scaled_difference,
    sum(!is.finite(as.matrix(analysis_data[, ..FEATURES])))
  ),
  expected = c(
    EXPECTED_N,
    EXPECTED_C1,
    EXPECTED_COUNTS,
    PRIMARY_IMPUTATION,
    0,
    0
  ),
  pass = c(
    nrow(labels) == EXPECTED_N,
    sum(labels$cluster_k2 == 1L) == EXPECTED_C1,
    observed_counts == EXPECTED_COUNTS,
    TRUE,
    max_scaled_difference <= 1e-8,
    all(is.finite(as.matrix(analysis_data[, ..FEATURES])))
  )
)

input_info <- file.info(INPUT)
manifest <- data.table(
  input_name = names(INPUT),
  path = unname(INPUT),
  size_bytes = input_info$size,
  modified_time = format(input_info$mtime, "%Y-%m-%d %H:%M:%S"),
  md5 = unname(tools::md5sum(INPUT))
)

out_long <- file.path(TABLE_DIR, "30_Table_D1_KDIGO_stratified_feature_profiles_long.csv")
out_wide <- file.path(TABLE_DIR, "30_Table_D1_KDIGO_stratified_feature_profiles_wide.csv")
out_qc <- file.path(LOG_DIR, "30_KDIGO_stratified_feature_profiles_QC.csv")
out_manifest <- file.path(PROVENANCE_DIR, "30_KDIGO_stratified_feature_profiles_input_manifest.csv")
out_log <- file.path(LOG_DIR, "30_KDIGO_stratified_feature_profiles_log.txt")
out_ok <- file.path(ANALYSIS_ROOT, "30_run_completed.ok")

write_csv_atomic(profile_long, out_long)
write_csv_atomic(profile_wide, out_wide)
write_csv_atomic(qc, out_qc)
write_csv_atomic(manifest, out_manifest)

log_lines <- c(
  "Domain 1 KDIGO-stratified feature profile post-processing",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Values: median (Q1-Q3), quantile type 7",
  "Source: locked MICE completed dataset #1 and locked plain-k-means K2 labels",
  paste0("Locked denominator: ", EXPECTED_N, "; C1: ", EXPECTED_C1),
  paste0("Groups: ", paste(names(observed_counts), observed_counts, sep = "=", collapse = "; ")),
  paste0("Maximum rebuilt-vs-locked standardized matrix difference: ", format(max_scaled_difference, scientific = TRUE)),
  paste0("Long table: ", out_long),
  paste0("Wide table: ", out_wide),
  paste0("All QC passed: ", all(qc$pass))
)
writeLines(log_lines, out_log, useBytes = TRUE)
writeLines(c(log_lines, paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))), out_ok, useBytes = TRUE)

cat(paste(log_lines, collapse = "\n"), "\n")
print(profile_wide)
