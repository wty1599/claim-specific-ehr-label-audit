## =====================================================================
## 34_domain1_union_stratification_and_buncr.R
## Targeted post-lock audit of partially independent renal/acid-base signal.
##
## Stratification uses raw clinical units; OOF models use the frozen
## standardized MICE-imputation-1 matrix and exact saved WP1 outer folds.
## BUN/creatinine is BUN_max / creatinine_max and is explicitly not a
## same-draw laboratory ratio.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "pROC", "digest", "mice")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), winslash = "/", mustWork = TRUE)
CLOSURE_ROOT <- normalizePath(
  file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"),
  winslash = "/", mustWork = TRUE
)
OUT_ROOT <- file.path(CLOSURE_ROOT, "union_stratification_and_buncr_20260718")
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"),
  logs = file.path(OUT_ROOT, "logs"),
  provenance = file.path(OUT_ROOT, "provenance")
)
invisible(lapply(c(list(root = OUT_ROOT), DIRS), dir.create, recursive = TRUE, showWarnings = FALSE))

SCRIPT_PATH <- normalizePath(
  file.path(CLOSURE_ROOT, "scripts", "34_domain1_union_stratification_and_buncr.R"),
  winslash = "/", mustWork = TRUE
)
INPUT <- list(
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  mice = file.path(PROJECT_ROOT, "output", "model", "mice_primary.rds"),
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  baseline = file.path(PROJECT_ROOT, "data", "baseline_covars.csv"),
  wp1_oof = file.path(CLOSURE_ROOT, "tables", "29_Table_D1_identity_oof_predictions.csv"),
  wp1_summary = file.path(CLOSURE_ROOT, "tables", "Table_D1_identity_models_oof.csv")
)
if (any(!file.exists(unlist(INPUT)))) {
  stop("Missing input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))
}

EXPECTED_N <- 20049L
EXPECTED_C1 <- 3992L
EXPECTED_FOLDS <- 10L
BOOT_B <- 2000L
BOOT_SEED <- 2026071834L

OUTPUT <- list(
  definitions = file.path(DIRS$tables, "34_Table_D1_stratification_definitions.csv"),
  auc = file.path(DIRS$tables, "34_Table_D1_union_stratified_oof_auc.csv"),
  predictions = file.path(DIRS$tables, "34_Table_D1_union_stratified_oof_predictions.csv"),
  denominator = file.path(DIRS$tables, "34_Table_D1_union_strata_denominators.csv"),
  buncr = file.path(DIRS$tables, "34_Table_D1_bunmax_creatinine_max_ratio.csv"),
  fold_audit = file.path(DIRS$logs, "34_union_stratified_fold_audit.csv"),
  qc = file.path(DIRS$logs, "34_union_buncr_QC.csv"),
  report = file.path(DIRS$logs, "34_union_buncr_report.md"),
  input_hash = file.path(DIRS$provenance, "34_input_sha256.csv"),
  output_hash = file.path(DIRS$provenance, "34_output_sha256_manifest.csv"),
  session = file.path(DIRS$provenance, "34_sessionInfo.txt"),
  completed = file.path(OUT_ROOT, "34_run_completed.ok")
)
existing <- unlist(OUTPUT)[file.exists(unlist(OUTPUT))]
if (length(existing)) stop("Refusing to overwrite existing output(s):\n", paste(existing, collapse = "\n"))

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

capture_glm <- function(formula, data) {
  warnings <- character()
  fit <- withCallingHandlers(
    glm(formula, data = data, family = binomial()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

run_subset_oof <- function(data, fold_id, subset_mask, analysis_id, subset_label,
                           role, rule, predictors, full_reference_auc) {
  idx <- which(subset_mask)
  d <- data[idx]
  folds <- fold_id[idx]
  if (!nrow(d) || uniqueN(d$C1) != 2L) stop("Invalid subset outcome support: ", analysis_id)
  probability <- rep(NA_real_, nrow(d))
  audits <- vector("list", EXPECTED_FOLDS)
  formula <- reformulate(predictors, response = "C1")

  for (fold in seq_len(EXPECTED_FOLDS)) {
    train_idx <- which(folds != fold)
    test_idx <- which(folds == fold)
    if (!length(test_idx) || uniqueN(d$C1[train_idx]) != 2L || uniqueN(d$C1[test_idx]) != 2L) {
      stop("Fold lacks both classes: ", analysis_id, " fold ", fold)
    }
    captured <- capture_glm(formula, d[train_idx])
    p_test <- as.numeric(predict(captured$fit, newdata = d[test_idx], type = "response"))
    if (any(!is.finite(p_test))) stop("Non-finite prediction: ", analysis_id, " fold ", fold)
    probability[test_idx] <- p_test
    audits[[fold]] <- data.table(
      analysis_id = analysis_id,
      fold = fold,
      n_train = length(train_idx),
      n_test = length(test_idx),
      train_c1_n = sum(d$C1[train_idx]),
      test_c1_n = sum(d$C1[test_idx]),
      converged = isTRUE(captured$fit$converged),
      warning_n = length(captured$warnings),
      warning_text = paste(captured$warnings, collapse = " | ")
    )
  }

  if (anyNA(probability)) stop("Incomplete OOF predictions: ", analysis_id)
  roc_obj <- pROC::roc(d$C1, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  list(
    summary = data.table(
      analysis_id = analysis_id,
      analysis_role = role,
      subset_label = subset_label,
      subset_rule = rule,
      predictors = paste(predictors, collapse = " + "),
      n = nrow(d),
      c1_n = sum(d$C1),
      c2_n = sum(d$C1 == 0L),
      c1_prevalence = mean(d$C1),
      auc = as.numeric(pROC::auc(roc_obj)),
      auc_ci_low = ci[1L],
      auc_ci_high = ci[3L],
      full_cohort_reference_auc = full_reference_auc,
      descriptive_auc_difference_from_full = as.numeric(pROC::auc(roc_obj)) - full_reference_auc,
      inference_note = "Subset and full-cohort AUCs have different denominators; difference is descriptive only"
    ),
    predictions = data.table(
      stay_id = d$stay_id,
      C1 = d$C1,
      fold = folds,
      analysis_id = analysis_id,
      oof_probability = probability
    ),
    audit = rbindlist(audits)
  )
}

bootstrap_median <- function(x, seed, B = BOOT_B) {
  if (!length(x) || any(!is.finite(x))) stop("Invalid ratio vector")
  set.seed(seed)
  n <- length(x)
  estimates <- replicate(B, stats::median(x[sample.int(n, n, replace = TRUE)]))
  q <- as.numeric(stats::quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7))
  list(low = q[1L], high = q[2L], success_n = sum(is.finite(estimates)))
}

cat("Reading frozen labels/matrix/folds and raw clinical-unit data...\n")
labels <- as.data.table(readRDS(INPUT$labels))
X <- as.data.table(readRDS(INPUT$matrix))
imp <- readRDS(INPUT$mice)
completed_imp1 <- as.data.table(mice::complete(imp, 1L))
raw <- fread(
  INPUT$final_full,
  select = c("stay_id", "gender", "creatinine_max", "bun_max", "bicarbonate_min"),
  showProgress = FALSE
)
baseline <- fread(INPUT$baseline, select = c("stay_id", "aki_stage_0_24h"), showProgress = FALSE)
wp1_oof <- fread(INPUT$wp1_oof, select = c("stay_id", "fold", "C1", "model_id"), showProgress = FALSE)
wp1_summary <- fread(INPUT$wp1_summary, showProgress = FALSE)

if (!identical(labels$stay_id, X$stay_id)) stop("Frozen label/matrix row-order mismatch")
if (nrow(labels) != EXPECTED_N || sum(labels$cluster_k2 == 1L) != EXPECTED_C1) stop("Frozen denominator mismatch")
if (nrow(completed_imp1) != EXPECTED_N) stop("MICE imputation-1 denominator mismatch")
if (anyDuplicated(raw$stay_id) || anyDuplicated(baseline$stay_id)) stop("Duplicate raw/baseline stay_id")

fold_check <- wp1_oof[, .(unique_fold_n = uniqueN(fold), unique_c1_n = uniqueN(C1)), by = stay_id]
if (nrow(fold_check) != EXPECTED_N || any(fold_check$unique_fold_n != 1L) || any(fold_check$unique_c1_n != 1L)) {
  stop("WP1 fold map is not invariant across models")
}
fold_map <- unique(wp1_oof[, .(stay_id, fold, C1)])
setkey(fold_map, stay_id)
fold_idx <- match(labels$stay_id, fold_map$stay_id)
raw_idx <- match(labels$stay_id, raw$stay_id)
base_idx <- match(labels$stay_id, baseline$stay_id)
if (anyNA(fold_idx) || anyNA(raw_idx) || anyNA(base_idx)) stop("Incomplete stay_id alignment")

C1 <- as.integer(labels$cluster_k2 == 1L)
fold_id <- as.integer(fold_map$fold[fold_idx])
if (!identical(C1, as.integer(fold_map$C1[fold_idx]))) stop("C1 mismatch against WP1 fold table")
raw <- raw[raw_idx]
baseline <- baseline[base_idx]
if (!identical(raw$stay_id, labels$stay_id) || !identical(baseline$stay_id, labels$stay_id)) stop("Alignment failure")
if (any(!raw$gender %in% c("M", "F"))) stop("Unexpected or missing gender")

requested_features <- c("creatinine_max", "bicarbonate_min", "aniongap_max", "ph_min", "pco2_max")
if (length(setdiff(requested_features, names(X)))) stop("Frozen matrix lacks requested feature")
if (any(!is.finite(as.matrix(X[, ..requested_features])))) stop("Non-finite frozen predictor")
raw_numeric <- c("creatinine_max", "bun_max", "bicarbonate_min")
if (any(raw$creatinine_max <= 0, na.rm = TRUE) || any(raw$bun_max < 0, na.rm = TRUE)) {
  stop("Invalid observed raw clinical-unit value")
}
completed_numeric <- c("creatinine_max", "bun_max", "bicarbonate_min")
if (length(setdiff(completed_numeric, names(completed_imp1))) ||
    any(!is.finite(as.matrix(completed_imp1[, ..completed_numeric]))) ||
    any(completed_imp1$creatinine_max <= 0) || any(completed_imp1$bun_max < 0)) {
  stop("Invalid completed-imputation clinical-unit value")
}

features33 <- setdiff(names(X), "stay_id")
if (!setequal(features33, names(completed_imp1))) stop("MICE/matrix feature mismatch")
restandardized_imp1 <- scale(as.matrix(completed_imp1[, ..features33]))
matrix_reproduction_max_abs_error <- max(abs(restandardized_imp1 - as.matrix(X[, ..features33])))
if (!is.finite(matrix_reproduction_max_abs_error) || matrix_reproduction_max_abs_error > 1e-12) {
  stop("MICE imputation-1 does not reproduce the frozen standardized matrix")
}

analysis_data <- data.table(stay_id = labels$stay_id, C1 = C1)
analysis_data <- cbind(analysis_data, X[, ..requested_features])

## No project-local prespecification was found for these normal-range cutoffs.
## Primary and sensitivity definitions are frozen here before inspecting results.
definitions <- data.table(
  definition_id = c(
    "hco3_no_low_primary", "hco3_no_low_sensitivity",
    "creatinine_no_high_primary", "creatinine_no_high_sensitivity", "creatinine_no_high_unisex"
  ),
  role = c("primary", "sensitivity", "primary", "sensitivity", "sensitivity"),
  variable = c("bicarbonate_min", "bicarbonate_min", "creatinine_max", "creatinine_max", "creatinine_max"),
  rule = c(
    "bicarbonate_min >= 23 mmol/L",
    "bicarbonate_min >= 22 mmol/L",
    "creatinine_max <= 1.35 mg/dL for M and <= 1.04 mg/dL for F",
    "creatinine_max <= 1.30 mg/dL for M and <= 0.95 mg/dL for F",
    "creatinine_max <= 1.20 mg/dL for all sexes"
  ),
  rationale = c(
    "absence of low serum bicarbonate using common chemistry lower reference limit",
    "alternative lower reference limit used in blood-gas references",
    "absence of high creatinine using sex-specific Mayo Clinic upper reference limits",
    "stricter sex-specific MedlinePlus upper reference limits",
    "simple sex-neutral sensitivity threshold"
  ),
  source = c(
    "https://medlineplus.gov/ency/article/003469.htm",
    "https://medlineplus.gov/lab-tests/arterial-blood-gas-abg-test/",
    "https://www.mayoclinic.org/tests-procedures/creatinine-test/about/pac-20384646",
    "https://medlineplus.gov/ency/article/003475.htm",
    "sensitivity definition; not a laboratory-specific reference range"
  ),
  prespecification_status = "fixed in script 34 before result inspection; not prespecified in the earlier project plan"
)

masks <- list(
  hco3_no_low_primary = completed_imp1$bicarbonate_min >= 23,
  hco3_no_low_sensitivity = completed_imp1$bicarbonate_min >= 22,
  creatinine_no_high_primary = (raw$gender == "M" & completed_imp1$creatinine_max <= 1.35) |
    (raw$gender == "F" & completed_imp1$creatinine_max <= 1.04),
  creatinine_no_high_sensitivity = (raw$gender == "M" & completed_imp1$creatinine_max <= 1.30) |
    (raw$gender == "F" & completed_imp1$creatinine_max <= 0.95),
  creatinine_no_high_unisex = completed_imp1$creatinine_max <= 1.20
)

creatinine_reference_auc <- wp1_summary[model_id == "creatinine", auc]
acid_reference_auc <- wp1_summary[model_id == "acid_base", auc]
if (length(creatinine_reference_auc) != 1L || length(acid_reference_auc) != 1L) stop("WP1 reference AUC missing")

specs <- list(
  list(id = "creatinine_in_hco3_ge23", subset = "No low HCO3: minimum >=23 mmol/L", role = "primary",
       rule = definitions[definition_id == "hco3_no_low_primary", rule], mask = masks$hco3_no_low_primary,
       predictors = "creatinine_max", reference = creatinine_reference_auc),
  list(id = "creatinine_in_hco3_ge22", subset = "No low HCO3: minimum >=22 mmol/L", role = "sensitivity",
       rule = definitions[definition_id == "hco3_no_low_sensitivity", rule], mask = masks$hco3_no_low_sensitivity,
       predictors = "creatinine_max", reference = creatinine_reference_auc),
  list(id = "acid_base_in_normal_creatinine_mayo", subset = "No high creatinine: sex-specific Mayo ULN", role = "primary",
       rule = definitions[definition_id == "creatinine_no_high_primary", rule], mask = masks$creatinine_no_high_primary,
       predictors = c("bicarbonate_min", "aniongap_max", "ph_min", "pco2_max"), reference = acid_reference_auc),
  list(id = "acid_base_in_normal_creatinine_medline", subset = "No high creatinine: stricter sex-specific MedlinePlus ULN", role = "sensitivity",
       rule = definitions[definition_id == "creatinine_no_high_sensitivity", rule], mask = masks$creatinine_no_high_sensitivity,
       predictors = c("bicarbonate_min", "aniongap_max", "ph_min", "pco2_max"), reference = acid_reference_auc),
  list(id = "acid_base_in_creatinine_le1p2", subset = "No high creatinine: <=1.20 mg/dL", role = "sensitivity",
       rule = definitions[definition_id == "creatinine_no_high_unisex", rule], mask = masks$creatinine_no_high_unisex,
       predictors = c("bicarbonate_min", "aniongap_max", "ph_min", "pco2_max"), reference = acid_reference_auc)
)

cat("Running same-fold subset OOF models...\n")
results <- lapply(specs, function(s) {
  run_subset_oof(
    analysis_data, fold_id, s$mask, s$id, s$subset, s$role,
    s$rule, s$predictors, s$reference
  )
})
auc_table <- rbindlist(lapply(results, `[[`, "summary"))
predictions <- rbindlist(lapply(results, `[[`, "predictions"))
fold_audit <- rbindlist(lapply(results, `[[`, "audit"))
denominator <- rbindlist(lapply(seq_along(specs), function(i) {
  s <- specs[[i]]
  idx <- which(s$mask)
  data.table(
    analysis_id = s$id,
    analysis_role = s$role,
    subset_rule = s$rule,
    n = length(idx),
    c1_n = sum(C1[idx]),
    c2_n = sum(C1[idx] == 0L),
    c1_prevalence = mean(C1[idx]),
    excluded_n = EXPECTED_N - length(idx)
  )
}))

## BUNmax/creatinine_max is based on separate 24-h extrema, not a same-draw ratio.
ratio_imp1 <- completed_imp1$bun_max / completed_imp1$creatinine_max
ratio_observed <- raw$bun_max / raw$creatinine_max
observed_complete <- is.finite(ratio_observed)
stage <- as.integer(baseline$aki_stage_0_24h)
ratio_groups <- list(
  all_C1 = C1 == 1L,
  all_C2 = C1 == 0L,
  stage0_C1 = C1 == 1L & stage == 0L,
  stage2_3_C1 = C1 == 1L & stage %in% c(2L, 3L)
)
ratio_labels <- c(
  all_C1 = "All C1",
  all_C2 = "All C2",
  stage0_C1 = "KDIGO stage 0 C1",
  stage2_3_C1 = "KDIGO stage 2-3 C1"
)
expected_group_n <- c(all_C1 = 3992L, all_C2 = 16057L, stage0_C1 = 969L, stage2_3_C1 = 2244L)

buncr_imp1 <- rbindlist(lapply(seq_along(ratio_groups), function(i) {
  id <- names(ratio_groups)[i]
  x <- ratio_imp1[ratio_groups[[i]]]
  boot <- bootstrap_median(x, BOOT_SEED + i)
  data.table(
    analysis_version = "MICE_imputation_1_primary",
    group_id = id,
    group_label = unname(ratio_labels[id]),
    n = length(x),
    expected_n = unname(expected_group_n[id]),
    bun_creatinine_ratio_definition = "bun_max_mg_dL / creatinine_max_mg_dL",
    median = median(x),
    q1 = as.numeric(quantile(x, 0.25, names = FALSE)),
    q3 = as.numeric(quantile(x, 0.75, names = FALSE)),
    percentile_bootstrap_median_ci_low = boot$low,
    percentile_bootstrap_median_ci_high = boot$high,
    bootstrap_success_n = boot$success_n,
    bootstrap_requested_n = BOOT_B,
    interpretation_limit = "BUN and creatinine maxima may occur at different times; not a same-draw BUN/Cr ratio"
  )
}))
buncr_observed <- rbindlist(lapply(seq_along(ratio_groups), function(i) {
  id <- names(ratio_groups)[i]
  keep <- ratio_groups[[i]] & observed_complete
  x <- ratio_observed[keep]
  boot <- bootstrap_median(x, BOOT_SEED + 100L + i)
  data.table(
    analysis_version = "observed_complete_case_sensitivity",
    group_id = id,
    group_label = unname(ratio_labels[id]),
    n = length(x),
    expected_n = unname(expected_group_n[id]),
    bun_creatinine_ratio_definition = "observed bun_max_mg_dL / observed creatinine_max_mg_dL",
    median = median(x),
    q1 = as.numeric(quantile(x, 0.25, names = FALSE)),
    q3 = as.numeric(quantile(x, 0.75, names = FALSE)),
    percentile_bootstrap_median_ci_low = boot$low,
    percentile_bootstrap_median_ci_high = boot$high,
    bootstrap_success_n = boot$success_n,
    bootstrap_requested_n = BOOT_B,
    interpretation_limit = "Complete-case sensitivity; BUN and creatinine maxima may occur at different times"
  )
}))
buncr <- rbindlist(list(buncr_imp1, buncr_observed), use.names = TRUE)

qc <- rbindlist(list(
  data.table(
    check = "MICE_imp1_reproduces_frozen_standardized_matrix",
    observed = sprintf("max absolute error=%.3g", matrix_reproduction_max_abs_error),
    required = "<=1e-12",
    pass = matrix_reproduction_max_abs_error <= 1e-12
  ),
  data.table(
    check = "locked_denominator_and_C1",
    observed = paste0(nrow(labels), " rows; C1=", sum(C1)),
    required = paste0(EXPECTED_N, " rows; C1=", EXPECTED_C1),
    pass = nrow(labels) == EXPECTED_N && sum(C1) == EXPECTED_C1
  ),
  data.table(
    check = "same_fold_map_invariant",
    observed = paste0(sum(fold_check$unique_fold_n == 1L), "/", EXPECTED_N),
    required = paste0(EXPECTED_N, "/", EXPECTED_N),
    pass = all(fold_check$unique_fold_n == 1L)
  ),
  data.table(
    check = "subset_models_converged_without_warning",
    observed = paste0(sum(fold_audit$converged & fold_audit$warning_n == 0L), "/", nrow(fold_audit), " folds"),
    required = paste0(nrow(fold_audit), "/", nrow(fold_audit), " folds"),
    pass = all(fold_audit$converged) && sum(fold_audit$warning_n) == 0L
  ),
  data.table(
    check = "all_subset_predictions_complete",
    observed = paste0(nrow(predictions), " predictions for ", nrow(auc_table), " analyses"),
    required = paste0(sum(denominator$n), " predictions"),
    pass = nrow(predictions) == sum(denominator$n) && all(is.finite(predictions$oof_probability))
  ),
  data.table(
    check = "ratio_group_counts_match_locked_KDIGO_cross_tab",
    observed = paste(paste0(buncr_imp1$group_id, "=", buncr_imp1$n), collapse = "; "),
    required = paste(paste0(names(expected_group_n), "=", expected_group_n), collapse = "; "),
    pass = identical(as.integer(buncr_imp1$n), as.integer(expected_group_n[buncr_imp1$group_id]))
  ),
  data.table(
    check = "ratio_bootstrap_complete",
    observed = paste0(sum(buncr$bootstrap_success_n), "/", BOOT_B * nrow(buncr)),
    required = paste0(BOOT_B * nrow(buncr), "/", BOOT_B * nrow(buncr)),
    pass = all(buncr$bootstrap_success_n == BOOT_B)
  )
))
if (!all(qc$pass)) stop("QC failure; no completion marker will be written")

write_csv_atomic(definitions, OUTPUT$definitions)
write_csv_atomic(auc_table, OUTPUT$auc)
write_csv_atomic(predictions, OUTPUT$predictions)
write_csv_atomic(denominator, OUTPUT$denominator)
write_csv_atomic(buncr, OUTPUT$buncr)
write_csv_atomic(fold_audit, OUTPUT$fold_audit)
write_csv_atomic(qc, OUTPUT$qc)

input_hash <- data.table(
  input = names(c(INPUT, script = SCRIPT_PATH)),
  path = unname(c(INPUT, script = SCRIPT_PATH))
)
input_hash[, sha256 := vapply(path, sha256, character(1))]
write_csv_atomic(input_hash, OUTPUT$input_hash)

report <- c(
  "# Domain 1 conditional-axis stratification and BUNmax/creatinine_max audit",
  "",
  "## Threshold governance",
  "",
  "No project-local prespecified normal-range thresholds were found. Primary and sensitivity definitions were frozen in script 34 before inspecting the results.",
  "",
  "## Stratified OOF AUC",
  "",
  paste(capture.output(print(auc_table[, .(analysis_id, analysis_role, n, c1_n, auc, auc_ci_low, auc_ci_high)])), collapse = "\n"),
  "",
  "## BUNmax/creatinine_max summary",
  "",
  paste(capture.output(print(buncr[, .(analysis_version, group_label, n, median, q1, q3, percentile_bootstrap_median_ci_low, percentile_bootstrap_median_ci_high)])), collapse = "\n"),
  "",
  "## Mandatory interpretation limits",
  "",
  "1. Subset AUCs reconstruct the frozen algorithmic C1 label; they do not validate a biological subtype.",
  "2. Persistence of both subset AUCs supports partially independent predictive information, but does not by itself establish disconnected biological subgroups or a logical OR rule.",
  "3. K=2 Euclidean nearest-centroid assignment has a single linear decision boundary in the frozen feature space.",
  "4. BUNmax/creatinine_max uses separate extrema that may occur at different times and must not be interpreted as a same-draw prerenal index."
)
write_lines_atomic(report, OUTPUT$report)
write_lines_atomic(capture.output(sessionInfo()), OUTPUT$session)

manifest_paths <- unlist(OUTPUT[c(
  "definitions", "auc", "predictions", "denominator", "buncr",
  "fold_audit", "qc", "report", "input_hash", "session"
)])
output_hash <- data.table(path = manifest_paths)
output_hash[, `:=`(
  bytes = file.info(path)$size,
  sha256 = vapply(path, sha256, character(1))
)]
write_csv_atomic(output_hash, OUTPUT$output_hash)
write_lines_atomic(c(
  "status=completed",
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", sha256(SCRIPT_PATH)),
  paste0("output_manifest_sha256=", sha256(OUTPUT$output_hash)),
  paste0("qc_pass=", sum(qc$pass), "/", nrow(qc))
), OUTPUT$completed)

cat("Completed. Output directory:\n", OUT_ROOT, "\n", sep = "")
print(auc_table[, .(analysis_id, analysis_role, n, c1_n, auc, auc_ci_low, auc_ci_high)])
print(buncr[, .(analysis_version, group_label, n, median, q1, q3, percentile_bootstrap_median_ci_low, percentile_bootstrap_median_ci_high)])
