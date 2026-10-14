## =====================================================================
## 33_domain1_acid_base_axis_falsification.R
## Targeted falsification audit for the empirical C1 identity statement.
##
## Inputs are frozen WP1 artifacts:
## - MICE-imputation-1 plain-k-means labels
## - locked standardized 33-feature matrix
## - exact outer-fold assignments saved by WP1
##
## New analyses:
## 1. Same-fold OOF AUC for bicarbonate alone and nested acid-base models.
## 2. Pearson correlation matrix for pH, pCO2, bicarbonate, and anion gap.
## 3. Pearson correlations of creatinine/BUN with bicarbonate.
##
## This script writes only to a new isolated output directory.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "pROC", "digest")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), winslash = "/", mustWork = TRUE)
CLOSURE_ROOT <- normalizePath(
  file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"),
  winslash = "/", mustWork = TRUE
)
OUT_ROOT <- file.path(CLOSURE_ROOT, "acid_base_axis_falsification_20260718")
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"),
  logs = file.path(OUT_ROOT, "logs"),
  provenance = file.path(OUT_ROOT, "provenance")
)
invisible(lapply(c(list(root = OUT_ROOT), DIRS), dir.create, recursive = TRUE, showWarnings = FALSE))

SCRIPT_PATH <- normalizePath(
  file.path(CLOSURE_ROOT, "scripts", "33_domain1_acid_base_axis_falsification.R"),
  winslash = "/", mustWork = TRUE
)
INPUT <- list(
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  wp1_oof = file.path(CLOSURE_ROOT, "tables", "29_Table_D1_identity_oof_predictions.csv"),
  wp1_summary = file.path(CLOSURE_ROOT, "tables", "Table_D1_identity_models_oof.csv")
)
if (any(!file.exists(unlist(INPUT)))) {
  stop("Missing input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))
}

EXPECTED_N <- 20049L
EXPECTED_C1 <- 3992L
EXPECTED_FOLDS <- 10L
REPRO_TOLERANCE <- 1e-12

OUTPUT <- list(
  auc = file.path(DIRS$tables, "33_Table_D1_acid_base_nested_oof_auc.csv"),
  predictions = file.path(DIRS$tables, "33_Table_D1_acid_base_nested_oof_predictions.csv"),
  correlation_wide = file.path(DIRS$tables, "33_Table_D1_acid_base_pearson_matrix.csv"),
  correlation_long = file.path(DIRS$tables, "33_Table_D1_acid_base_pearson_long.csv"),
  renal_hco3 = file.path(DIRS$tables, "33_Table_D1_renal_hco3_pearson.csv"),
  fold_audit = file.path(DIRS$logs, "33_acid_base_fold_audit.csv"),
  qc = file.path(DIRS$logs, "33_acid_base_axis_QC.csv"),
  report = file.path(DIRS$logs, "33_acid_base_axis_report.md"),
  input_hash = file.path(DIRS$provenance, "33_input_sha256.csv"),
  output_hash = file.path(DIRS$provenance, "33_output_sha256_manifest.csv"),
  session = file.path(DIRS$provenance, "33_sessionInfo.txt"),
  completed = file.path(OUT_ROOT, "33_run_completed.ok")
)
existing <- unlist(OUTPUT)[file.exists(unlist(OUTPUT))]
if (length(existing)) {
  stop("Refusing to overwrite existing output(s):\n", paste(existing, collapse = "\n"))
}

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

capture_fit <- function(formula, data) {
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

run_oof_model <- function(data, fold_id, model_id, model_label, predictors) {
  probability <- rep(NA_real_, nrow(data))
  audits <- vector("list", EXPECTED_FOLDS)
  formula <- reformulate(predictors, response = "C1")

  for (fold in seq_len(EXPECTED_FOLDS)) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    captured <- capture_fit(formula, data[train_idx])
    p_test <- as.numeric(predict(captured$fit, newdata = data[test_idx], type = "response"))
    if (any(!is.finite(p_test))) stop("Non-finite OOF probability: ", model_id, " fold ", fold)
    probability[test_idx] <- p_test
    audits[[fold]] <- data.table(
      model_id = model_id,
      fold = fold,
      n_train = length(train_idx),
      n_test = length(test_idx),
      train_c1_n = sum(data$C1[train_idx]),
      test_c1_n = sum(data$C1[test_idx]),
      converged = isTRUE(captured$fit$converged),
      warning_n = length(captured$warnings),
      warning_text = paste(captured$warnings, collapse = " | ")
    )
  }

  if (anyNA(probability)) stop("Incomplete OOF predictions: ", model_id)
  roc_obj <- pROC::roc(data$C1, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  list(
    summary = data.table(
      model_id = model_id,
      model_label = model_label,
      predictors = paste(predictors, collapse = " + "),
      predictor_n = length(predictors),
      evaluation = "same locked 10-fold stratified OOF as WP1",
      n = nrow(data),
      c1_n = sum(data$C1),
      auc = as.numeric(pROC::auc(roc_obj)),
      auc_ci_low = auc_ci[1L],
      auc_ci_high = auc_ci[3L],
      auc_ci_method = "DeLong CI on one OOF prediction per patient"
    ),
    predictions = data.table(
      stay_id = data$stay_id,
      C1 = data$C1,
      fold = fold_id,
      model_id = model_id,
      oof_probability = probability
    ),
    audit = rbindlist(audits)
  )
}

cat("Reading frozen labels, matrix, and exact WP1 folds...\n")
labels <- as.data.table(readRDS(INPUT$labels))
X <- as.data.table(readRDS(INPUT$matrix))
wp1_oof <- fread(INPUT$wp1_oof, showProgress = FALSE)
wp1_summary <- fread(INPUT$wp1_summary, showProgress = FALSE)

if (!identical(labels$stay_id, X$stay_id)) stop("Label/matrix row-order mismatch")
if (nrow(labels) != EXPECTED_N || sum(labels$cluster_k2 == 1L) != EXPECTED_C1) {
  stop("Frozen denominator or C1 count mismatch")
}

required_features <- c(
  "bicarbonate_min", "aniongap_max", "ph_min", "pco2_max",
  "creatinine_max", "bun_max"
)
if (length(setdiff(required_features, names(X)))) {
  stop("Missing frozen feature(s): ", paste(setdiff(required_features, names(X)), collapse = ", "))
}

fold_consistency <- wp1_oof[, .(
  unique_fold_n = uniqueN(fold),
  unique_c1_n = uniqueN(C1)
), by = stay_id]
if (nrow(fold_consistency) != EXPECTED_N || any(fold_consistency$unique_fold_n != 1L) ||
    any(fold_consistency$unique_c1_n != 1L)) {
  stop("WP1 fold/C1 assignments are not invariant across saved models")
}
fold_map <- unique(wp1_oof[, .(stay_id, fold, C1)])
setkey(fold_map, stay_id)
fold_idx <- match(labels$stay_id, fold_map$stay_id)
if (anyNA(fold_idx)) stop("Frozen stay absent from WP1 fold map")
fold_id <- as.integer(fold_map$fold[fold_idx])
C1 <- as.integer(labels$cluster_k2 == 1L)
if (!identical(C1, as.integer(fold_map$C1[fold_idx]))) stop("C1 mismatch against WP1 OOF output")
if (!identical(sort(unique(fold_id)), seq_len(EXPECTED_FOLDS))) stop("Unexpected fold IDs")

analysis_data <- data.table(stay_id = labels$stay_id, C1 = C1)
analysis_data <- cbind(analysis_data, X[, ..required_features])
if (any(!is.finite(as.matrix(analysis_data[, ..required_features])))) stop("Non-finite frozen feature value")

models <- list(
  list(
    id = "hco3",
    label = "Bicarbonate minimum",
    predictors = c("bicarbonate_min")
  ),
  list(
    id = "hco3_ag",
    label = "Bicarbonate + anion gap",
    predictors = c("bicarbonate_min", "aniongap_max")
  ),
  list(
    id = "hco3_ag_ph",
    label = "Bicarbonate + anion gap + pH",
    predictors = c("bicarbonate_min", "aniongap_max", "ph_min")
  ),
  list(
    id = "hco3_ag_ph_pco2",
    label = "Bicarbonate + anion gap + pH + pCO2",
    predictors = c("bicarbonate_min", "aniongap_max", "ph_min", "pco2_max")
  )
)

cat("Running four same-fold nested OOF logistic models...\n")
results <- lapply(models, function(m) {
  run_oof_model(analysis_data, fold_id, m$id, m$label, m$predictors)
})
auc_table <- rbindlist(lapply(results, `[[`, "summary"))
auc_table[, `:=`(
  previous_model_id = shift(model_id),
  delta_auc_vs_previous = auc - shift(auc),
  delta_auc_vs_hco3 = auc - first(auc)
)]
predictions <- rbindlist(lapply(results, `[[`, "predictions"))
fold_audit <- rbindlist(lapply(results, `[[`, "audit"))

locked_acid_auc <- wp1_summary[model_id == "acid_base", auc]
if (length(locked_acid_auc) != 1L || !is.finite(locked_acid_auc)) stop("Locked WP1 acid-base AUC unavailable")
reproduced_acid_auc <- auc_table[model_id == "hco3_ag_ph_pco2", auc]
reproduction_abs_error <- abs(reproduced_acid_auc - locked_acid_auc)

acid_vars <- c("ph_min", "pco2_max", "bicarbonate_min", "aniongap_max")
acid_labels <- c(
  ph_min = "pH minimum",
  pco2_max = "pCO2 maximum",
  bicarbonate_min = "Bicarbonate minimum",
  aniongap_max = "Anion gap maximum"
)
acid_matrix <- as.matrix(X[, ..acid_vars])
acid_cor <- stats::cor(acid_matrix, method = "pearson")
correlation_wide <- data.table(variable = unname(acid_labels[rownames(acid_cor)]), as.data.table(acid_cor))
setnames(correlation_wide, acid_vars, unname(acid_labels[acid_vars]))
correlation_long <- as.data.table(as.table(acid_cor))
if (ncol(correlation_long) != 3L) stop("Unexpected long correlation-table shape")
setnames(correlation_long, c("variable_1", "variable_2", "pearson_r"))
correlation_long[, `:=`(
  variable_1_label = unname(acid_labels[as.character(variable_1)]),
  variable_2_label = unname(acid_labels[as.character(variable_2)]),
  n = EXPECTED_N,
  matrix_source = "locked standardized MICE-imputation-1 33-feature matrix"
)]
setcolorder(correlation_long, c(
  "variable_1", "variable_1_label", "variable_2", "variable_2_label",
  "pearson_r", "n", "matrix_source"
))

renal_hco3 <- rbindlist(lapply(c("creatinine_max", "bun_max"), function(v) {
  data.table(
    renal_variable = v,
    acid_base_variable = "bicarbonate_min",
    pearson_r = stats::cor(X[[v]], X[["bicarbonate_min"]], method = "pearson"),
    n = EXPECTED_N,
    matrix_source = "locked standardized MICE-imputation-1 33-feature matrix"
  )
}))

qc <- rbindlist(list(
  data.table(
    check = "locked_denominator_and_C1",
    observed = paste0(nrow(labels), " rows; C1=", sum(C1)),
    required = paste0(EXPECTED_N, " rows; C1=", EXPECTED_C1),
    pass = nrow(labels) == EXPECTED_N && sum(C1) == EXPECTED_C1
  ),
  data.table(
    check = "same_fold_assignments_invariant_across_WP1_models",
    observed = paste0(sum(fold_consistency$unique_fold_n == 1L), "/", EXPECTED_N),
    required = paste0(EXPECTED_N, "/", EXPECTED_N),
    pass = all(fold_consistency$unique_fold_n == 1L)
  ),
  data.table(
    check = "all_new_models_converged_without_warning",
    observed = paste0(sum(fold_audit$converged & fold_audit$warning_n == 0L), "/", nrow(fold_audit), " folds"),
    required = paste0(nrow(fold_audit), "/", nrow(fold_audit), " folds"),
    pass = all(fold_audit$converged) && sum(fold_audit$warning_n) == 0L
  ),
  data.table(
    check = "four_variable_AUC_reproduces_locked_WP1_acid_base_AUC",
    observed = sprintf("new=%.15f; locked=%.15f; abs_error=%.3g", reproduced_acid_auc, locked_acid_auc, reproduction_abs_error),
    required = paste0("absolute error <= ", REPRO_TOLERANCE),
    pass = reproduction_abs_error <= REPRO_TOLERANCE
  ),
  data.table(
    check = "all_requested_correlations_finite",
    observed = paste0(sum(is.finite(c(acid_cor, renal_hco3$pearson_r))), "/", length(c(acid_cor, renal_hco3$pearson_r))),
    required = paste0(length(c(acid_cor, renal_hco3$pearson_r)), "/", length(c(acid_cor, renal_hco3$pearson_r))),
    pass = all(is.finite(c(acid_cor, renal_hco3$pearson_r)))
  )
))
if (!all(qc$pass)) stop("QC failure; no completion marker will be written")

write_csv_atomic(auc_table, OUTPUT$auc)
write_csv_atomic(predictions, OUTPUT$predictions)
write_csv_atomic(correlation_wide, OUTPUT$correlation_wide)
write_csv_atomic(correlation_long, OUTPUT$correlation_long)
write_csv_atomic(renal_hco3, OUTPUT$renal_hco3)
write_csv_atomic(fold_audit, OUTPUT$fold_audit)
write_csv_atomic(qc, OUTPUT$qc)

input_hash <- data.table(
  input = names(c(INPUT, script = SCRIPT_PATH)),
  path = unname(c(INPUT, script = SCRIPT_PATH))
)
input_hash[, sha256 := vapply(path, sha256, character(1))]
write_csv_atomic(input_hash, OUTPUT$input_hash)

ag_hco3 <- acid_cor["aniongap_max", "bicarbonate_min"]
creat_hco3 <- renal_hco3[renal_variable == "creatinine_max", pearson_r]
bun_hco3 <- renal_hco3[renal_variable == "bun_max", pearson_r]
report <- c(
  "# Domain 1 acid-base axis falsification audit",
  "",
  paste0("- Locked n=", EXPECTED_N, "; C1=", EXPECTED_C1, "."),
  "- The locked strict acid-base block already excluded lactate.",
  "- All AUCs use the exact saved WP1 outer folds and logistic OOF procedure.",
  paste0("- Bicarbonate-only OOF AUC: ", sprintf("%.6f", auc_table[model_id == "hco3", auc]), "."),
  paste0("- Full four-variable OOF AUC: ", sprintf("%.6f", reproduced_acid_auc), "."),
  paste0("- Anion gap vs bicarbonate Pearson r: ", sprintf("%.6f", ag_hco3), "."),
  paste0("- Creatinine vs bicarbonate Pearson r: ", sprintf("%.6f", creat_hco3), "."),
  paste0("- BUN vs bicarbonate Pearson r: ", sprintf("%.6f", bun_hco3), "."),
  "",
  "## Nested AUC table",
  "",
  paste(capture.output(print(auc_table[, .(model_id, predictors, auc, auc_ci_low, auc_ci_high, delta_auc_vs_previous)])), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "These results reconstruct an algorithmic label from variables used to create that label. They characterize dimensionality and redundancy; they do not establish an independent biological subtype or causal mechanism."
)
write_lines_atomic(report, OUTPUT$report)
write_lines_atomic(capture.output(sessionInfo()), OUTPUT$session)

manifest_paths <- unlist(OUTPUT[c(
  "auc", "predictions", "correlation_wide", "correlation_long", "renal_hco3",
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
print(auc_table[, .(model_id, predictors, auc, auc_ci_low, auc_ci_high, delta_auc_vs_previous)])
print(correlation_wide)
print(renal_hco3)
