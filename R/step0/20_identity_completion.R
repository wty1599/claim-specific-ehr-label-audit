## =====================================================================
## 29_domain1_renal_identity_completion.R
## Formal completion of the empirical Domain 1 clinical-identity audit.
##
## Outcome: locked MICE-imputation-1 plain-k-means C1 membership.
## Evaluation: identical 10-fold stratified OOF predictions for every model.
## Primary classification threshold: training-fold prevalence matching.
## Sensitivity threshold: training-fold Youden index.
## Pairwise inference: paired patient bootstrap, B = 2,000.
##
## This script reads locked source artifacts and writes only to the isolated
## closure directory. It does not modify source data or earlier results.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "pROC", "mclust", "ggplot2", "glmnet", "digest")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
  library(ggplot2)
})

WORK_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
if (!nzchar(WORK_ROOT)) stop("Set EHR_AUDIT_WORK_ROOT to the authorized workspace.")
PROJECT_ROOT <- normalizePath(WORK_ROOT, winslash = "/", mustWork = TRUE)
OUT_ROOT <- normalizePath(
  file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"),
  winslash = "/", mustWork = TRUE
)
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"),
  figures = file.path(OUT_ROOT, "figures"),
  logs = file.path(OUT_ROOT, "logs"),
  checkpoints = file.path(OUT_ROOT, "checkpoints"),
  provenance = file.path(OUT_ROOT, "provenance")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

SCRIPT_PATH <- normalizePath(
  file.path(OUT_ROOT, "scripts", "29_domain1_renal_identity_completion.R"),
  winslash = "/", mustWork = TRUE
)
RUN_MODE <- "formal"
RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid())
ALLOW_OVERWRITE <- identical(tolower(Sys.getenv("D1_ALLOW_OVERWRITE", "false")), "true")
RUN_IN_PROGRESS <- file.path(OUT_ROOT, "29_run_in_progress.flag")
RUN_COMPLETED <- file.path(OUT_ROOT, "29_run_completed.ok")
OUTPUT_MANIFEST <- file.path(DIRS$provenance, "29_WP1_output_sha256_manifest.csv")

INPUT <- list(
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  baseline = file.path(PROJECT_ROOT, "data", "baseline_covars.csv"),
  registry = file.path(DIRS$provenance, "D1_feature_domain_registry.csv"),
  prior_script_28 = file.path(PROJECT_ROOT, "domain1_next_phase_20260717", "scripts", "28_domain1_renal_identity_audit.R"),
  prior_oof_28 = file.path(PROJECT_ROOT, "domain1_next_phase_20260717", "tables", "28_Table_D1_renal_identity_oof_predictions.csv"),
  prior_auc_28 = file.path(PROJECT_ROOT, "domain1_next_phase_20260717", "tables", "28_Table_D1_renal_identity_auc.csv")
)
if (any(!file.exists(unlist(INPUT)))) {
  stop("Missing required input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))
}

EXPECTED_N <- 20049L
EXPECTED_C1 <- 3992L
CV_FOLDS <- 10L
CV_SEED <- 2026071728L
BOOT_B <- 2000L
BOOT_SEED <- 202607291L
PRIMARY_THRESHOLD <- "closest_achievable_training_prevalence_match_under_ties"
SENSITIVITY_THRESHOLD <- "training_youden"

FINAL_OUTPUTS <- c(
  file.path(DIRS$tables, "Table_D1_identity_models_oof.csv"),
  file.path(DIRS$tables, "Table_D1_identity_pairwise_auc.csv"),
  file.path(DIRS$tables, "Table_D1_identity_threshold_metrics.csv"),
  file.path(DIRS$tables, "Table_D1_identity_discrimination_retention.csv"),
  file.path(DIRS$tables, "29_Table_D1_identity_oof_predictions.csv"),
  file.path(DIRS$logs, "29_cross_validation_fold_audit.csv"),
  file.path(DIRS$logs, "29_same_fold_reproduction_audit.csv"),
  file.path(DIRS$logs, "29_WP1_QC_checks.csv"),
  file.path(DIRS$logs, "29_domain1_renal_identity_completion_log.md"),
  file.path(DIRS$logs, "29_sessionInfo.txt"),
  file.path(DIRS$provenance, "29_WP1_input_checksums.csv"),
  file.path(DIRS$provenance, "29_WP1_run_metadata.csv"),
  OUTPUT_MANIFEST,
  file.path(DIRS$figures, paste0("Figure_D1_identity_model_comparison.", c("pdf", "png", "tiff"))),
  RUN_COMPLETED
)
existing_outputs <- FINAL_OUTPUTS[file.exists(FINAL_OUTPUTS)]
if (length(existing_outputs) && !ALLOW_OVERWRITE) {
  stop(
    "Refusing to overwrite existing WP1 output(s). Archive them or set D1_ALLOW_OVERWRITE=true explicitly:\n",
    paste(existing_outputs, collapse = "\n")
  )
}
if (file.exists(RUN_IN_PROGRESS)) {
  stop("Existing in-progress marker requires manual review: ", RUN_IN_PROGRESS)
}

write_csv_atomic <- function(x, path) {
  if (file.exists(path) && !ALLOW_OVERWRITE) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  if (file.exists(path) && !ALLOW_OVERWRITE) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
  invisible(path)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)
clamp_probability <- function(p) pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)

START_HASH_PATHS <- c(INPUT, current_script_29 = SCRIPT_PATH)
START_HASHES <- vapply(unlist(START_HASH_PATHS, use.names = FALSE), sha256, character(1))
names(START_HASHES) <- names(START_HASH_PATHS)

write_lines_atomic(c(
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script=", SCRIPT_PATH),
  paste0("script_sha256_at_start=", START_HASHES[["current_script_29"]]),
  "internal_provenance_reference_not_distributed=true"
), RUN_IN_PROGRESS)

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (value in sort(unique(y))) {
    idx <- sample(which(y == value), replace = FALSE)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

prevalence_threshold <- function(probability, prevalence) {
  probability <- as.numeric(probability)
  if (!length(probability) || any(!is.finite(probability))) stop("Invalid prevalence-threshold probabilities")
  if (!is.finite(prevalence) || prevalence < 0 || prevalence > 1) stop("Invalid target prevalence")
  n <- length(probability)
  target_n <- round(prevalence * n)
  groups <- data.table(probability = probability)[, .(boundary_tie_n = .N), by = probability]
  setorder(groups, -probability)
  groups[, achieved_positive_n := cumsum(boundary_tie_n)]
  candidates <- rbindlist(list(
    data.table(
      threshold = Inf,
      achieved_positive_n = 0L,
      boundary_tie_n = 0L,
      threshold_location = "above_maximum"
    ),
    groups[, .(
      threshold = probability,
      achieved_positive_n,
      boundary_tie_n,
      threshold_location = "observed_probability"
    )]
  ))
  candidates[, `:=`(
    achieved_prevalence = achieved_positive_n / n,
    absolute_count_deviation = abs(achieved_positive_n - target_n),
    absolute_prevalence_deviation = abs(achieved_positive_n / n - prevalence),
    over_target = achieved_positive_n > target_n
  )]
  # Deterministic tie-break: closest count, then conservative (not above target), then higher threshold.
  chosen <- candidates[order(absolute_count_deviation, over_target, -threshold)][1L]
  list(
    threshold = chosen$threshold,
    target_positive_n = target_n,
    achieved_positive_n = chosen$achieved_positive_n,
    target_prevalence = prevalence,
    achieved_prevalence = chosen$achieved_prevalence,
    absolute_prevalence_deviation = chosen$absolute_prevalence_deviation,
    absolute_count_deviation = chosen$absolute_count_deviation,
    exact_target_count_achievable = chosen$achieved_positive_n == target_n,
    boundary_tie_n = chosen$boundary_tie_n,
    threshold_location = chosen$threshold_location,
    selection_rule = "closest achievable prevalence; ties prefer not exceeding target, then higher threshold"
  )
}

youden_threshold <- function(y, probability) {
  warnings <- character()
  error_message <- ""
  out <- tryCatch(
    withCallingHandlers({
      roc_obj <- pROC::roc(y, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
      pROC::coords(
        roc_obj, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE
      )
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  thresholds <- if (is.null(out)) numeric() else as.numeric(out$threshold)
  finite_thresholds <- thresholds[is.finite(thresholds)]
  status <- if (nzchar(error_message)) {
    "error"
  } else if (!length(finite_thresholds)) {
    "nonfinite"
  } else if (length(warnings)) {
    "warning"
  } else {
    "ok"
  }
  list(
    threshold = if (length(finite_thresholds)) max(finite_thresholds) else NA_real_,
    status = status,
    candidate_n = length(finite_thresholds),
    selection_rule = "highest threshold among tied maximum-Youden candidates (specificity-prioritizing deterministic tie-break)",
    message = paste(c(error_message, unique(warnings)), collapse = " | ")
  )
}

average_precision <- function(y, probability) {
  y <- as.integer(y)
  positives <- sum(y == 1L)
  if (!positives || positives == length(y)) return(NA_real_)
  if (length(probability) != length(y) || any(!is.finite(probability))) return(NA_real_)
  groups <- data.table(y = y, probability = as.numeric(probability))[
    , .(n = .N, positives = sum(y == 1L)), by = probability
  ]
  setorder(groups, -probability)
  groups[, `:=`(
    cumulative_n = cumsum(n),
    cumulative_positives = cumsum(positives)
  )]
  groups[, `:=`(
    precision = cumulative_positives / cumulative_n,
    recall = cumulative_positives / sum(positives)
  )]
  recall_increment <- c(groups$recall[1L], diff(groups$recall))
  sum(recall_increment * groups$precision)
}

binary_metrics <- function(y, predicted_class) {
  y <- as.integer(y)
  predicted_class <- as.integer(predicted_class)
  tp <- sum(y == 1L & predicted_class == 1L)
  tn <- sum(y == 0L & predicted_class == 0L)
  fp <- sum(y == 0L & predicted_class == 1L)
  fn <- sum(y == 1L & predicted_class == 0L)
  n <- length(y)
  observed_agreement <- (tp + tn) / n
  py <- mean(y == 1L)
  pp <- mean(predicted_class == 1L)
  expected_agreement <- py * pp + (1 - py) * (1 - pp)
  kappa <- if (abs(1 - expected_agreement) < 1e-12) NA_real_ else {
    (observed_agreement - expected_agreement) / (1 - expected_agreement)
  }
  sensitivity <- safe_div(tp, tp + fn)
  specificity <- safe_div(tn, tn + fp)
  data.table(
    tp = tp, tn = tn, fp = fp, fn = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = safe_div(tp, tp + fp),
    npv = safe_div(tn, tn + fn),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    accuracy = observed_agreement,
    kappa = kappa,
    ari = mclust::adjustedRandIndex(y, predicted_class),
    predicted_c1_prevalence = pp
  )
}

calibration_metrics <- function(y, probability) {
  lp <- qlogis(clamp_probability(probability))
  audit_glm <- function(formula, offset = NULL) {
    warnings <- character()
    error_message <- ""
    args <- list(formula = formula, family = binomial())
    if (!is.null(offset)) args$offset <- offset
    fit <- tryCatch(
      withCallingHandlers(
        do.call(glm, args),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    coefficient_finite <- !is.null(fit) && all(is.finite(coef(fit)))
    status <- if (nzchar(error_message)) {
      "error"
    } else if (is.null(fit) || !coefficient_finite) {
      "nonfinite"
    } else if (!isTRUE(fit$converged)) {
      "not_converged"
    } else if (length(warnings)) {
      "warning"
    } else {
      "ok"
    }
    list(
      fit = fit,
      status = status,
      message = paste(c(error_message, unique(warnings)), collapse = " | ")
    )
  }
  intercept_audit <- audit_glm(y ~ 1, offset = lp)
  slope_audit <- audit_glm(y ~ lp)
  data.table(
    calibration_intercept = if (is.null(intercept_audit$fit)) NA_real_ else unname(coef(intercept_audit$fit)[1L]),
    calibration_slope = if (is.null(slope_audit$fit) || length(coef(slope_audit$fit)) < 2L) {
      NA_real_
    } else {
      unname(coef(slope_audit$fit)[2L])
    },
    calibration_intercept_status = intercept_audit$status,
    calibration_intercept_message = intercept_audit$message,
    calibration_slope_status = slope_audit$status,
    calibration_slope_message = slope_audit$message,
    calibration_status = if (intercept_audit$status == "ok" && slope_audit$status == "ok") {
      "ok"
    } else if (intercept_audit$status %in% c("ok", "warning") &&
               slope_audit$status %in% c("ok", "warning")) {
      "warning_finite"
    } else {
      "failed_audit"
    }
  )
}

fast_auc <- function(y, probability) {
  y <- as.integer(y)
  n1 <- sum(y == 1L)
  n0 <- length(y) - n1
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

inner_stratified_folds <- function(y, k, seed) make_stratified_folds(y, k, seed)

fit_predict_fold <- function(data, x_all33, model, train_idx, test_idx, fold) {
  warnings <- character()
  if (identical(model$engine, "ridge")) {
    x_train <- x_all33[train_idx, , drop = FALSE]
    x_test <- x_all33[test_idx, , drop = FALSE]
    y_train <- data$C1[train_idx]
    inner_fold <- inner_stratified_folds(y_train, 10L, CV_SEED + 1000L + fold)
    fit <- withCallingHandlers(
      glmnet::cv.glmnet(
        x = x_train, y = y_train, family = "binomial", alpha = 0,
        standardize = FALSE, nfolds = 10L, foldid = inner_fold,
        type.measure = "deviance", keep = FALSE, parallel = FALSE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    p_train <- as.numeric(predict(fit, newx = x_train, s = "lambda.min", type = "response"))
    p_test <- as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
    tuning <- data.table(lambda_min = fit$lambda.min, lambda_1se = fit$lambda.1se)
  } else {
    formula <- reformulate(model$predictors, response = "C1")
    fit <- withCallingHandlers(
      glm(formula, data = data[train_idx], family = binomial()),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    p_train <- as.numeric(predict(fit, newdata = data[train_idx], type = "response"))
    p_test <- as.numeric(predict(fit, newdata = data[test_idx], type = "response"))
    tuning <- data.table(lambda_min = NA_real_, lambda_1se = NA_real_)
  }
  if (any(!is.finite(c(p_train, p_test)))) stop("Non-finite probability: ", model$id, " fold ", fold)
  train_prevalence <- mean(data$C1[train_idx])
  threshold_primary <- prevalence_threshold(p_train, train_prevalence)
  threshold_youden <- youden_threshold(data$C1[train_idx], p_train)
  list(
    probability = p_test,
    class_primary = as.integer(p_test >= threshold_primary$threshold),
    class_youden = if (is.finite(threshold_youden$threshold)) {
      as.integer(p_test >= threshold_youden$threshold)
    } else {
      rep(NA_integer_, length(p_test))
    },
    audit = cbind(data.table(
      model_id = model$id,
      fold = fold,
      n_train = length(train_idx),
      n_test = length(test_idx),
      train_c1_prevalence = train_prevalence,
      test_c1_prevalence = mean(data$C1[test_idx]),
      threshold_primary = threshold_primary$threshold,
      target_train_positive_n_primary = threshold_primary$target_positive_n,
      achieved_train_positive_n_primary = threshold_primary$achieved_positive_n,
      target_train_prevalence_primary = threshold_primary$target_prevalence,
      achieved_train_prevalence_primary = threshold_primary$achieved_prevalence,
      absolute_train_prevalence_deviation_primary = threshold_primary$absolute_prevalence_deviation,
      absolute_train_count_deviation_primary = threshold_primary$absolute_count_deviation,
      exact_target_count_achievable_primary = threshold_primary$exact_target_count_achievable,
      primary_boundary_tie_n = threshold_primary$boundary_tie_n,
      primary_threshold_location = threshold_primary$threshold_location,
      primary_threshold_selection_rule = threshold_primary$selection_rule,
      threshold_youden = threshold_youden$threshold,
      youden_status = threshold_youden$status,
      youden_candidate_n = threshold_youden$candidate_n,
      youden_selection_rule = threshold_youden$selection_rule,
      youden_message = threshold_youden$message,
      warning_count = length(warnings),
      warning_text = paste(unique(warnings), collapse = " | ")
    ), tuning)
  )
}

cross_validated_model <- function(data, x_all33, fold_id, model) {
  n <- nrow(data)
  probability <- rep(NA_real_, n)
  class_primary <- rep(NA_integer_, n)
  class_youden <- rep(NA_integer_, n)
  audits <- vector("list", CV_FOLDS)
  for (fold in seq_len(CV_FOLDS)) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    result <- fit_predict_fold(data, x_all33, model, train_idx, test_idx, fold)
    probability[test_idx] <- result$probability
    class_primary[test_idx] <- result$class_primary
    class_youden[test_idx] <- result$class_youden
    audits[[fold]] <- result$audit
  }
  if (anyNA(probability) || anyNA(class_primary) || anyNA(class_youden)) {
    stop("Incomplete OOF prediction: ", model$id)
  }
  roc_obj <- pROC::roc(data$C1, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  cal <- calibration_metrics(data$C1, probability)
  primary_metrics <- binary_metrics(data$C1, class_primary)
  youden_metrics <- binary_metrics(data$C1, class_youden)
  fold_audit <- rbindlist(audits, fill = TRUE)
  summary <- cbind(data.table(
    model_id = model$id,
    model_label = model$label,
    predictors = paste(model$predictors, collapse = " + "),
    engine = model$engine,
    evaluation = "10-fold stratified out-of-fold",
    n = n,
    c1_n = sum(data$C1),
    c1_prevalence = mean(data$C1),
    auc = as.numeric(pROC::auc(roc_obj)),
    auc_ci_low = auc_ci[1L],
    auc_ci_high = auc_ci[3L],
    auc_ci_method = "DeLong CI on one OOF prediction per patient",
    pr_auc_average_precision = average_precision(data$C1, probability),
    pr_auc_method = "tie-invariant grouped average precision (step-function precision-recall area)",
    brier = mean((probability - data$C1)^2),
    fold_warning_n = sum(fold_audit$warning_count),
    youden_failed_fold_n = sum(fold_audit$youden_status != "ok"),
    primary_exact_match_fold_n = sum(fold_audit$exact_target_count_achievable_primary),
    primary_max_abs_prevalence_deviation = max(fold_audit$absolute_train_prevalence_deviation_primary)
  ), cal, data.table(threshold_rule = PRIMARY_THRESHOLD), primary_metrics)
  threshold_summary <- rbindlist(list(
    cbind(data.table(model_id = model$id, model_label = model$label, threshold_rule = PRIMARY_THRESHOLD), primary_metrics),
    cbind(data.table(model_id = model$id, model_label = model$label, threshold_rule = SENSITIVITY_THRESHOLD), youden_metrics)
  ), fill = TRUE)
  list(
    summary = summary,
    threshold_summary = threshold_summary,
    prediction = data.table(
      subject_id = data$subject_id, stay_id = data$stay_id, fold = fold_id, C1 = data$C1,
      model_id = model$id, oof_probability = probability,
      oof_class_prevalence = class_primary, oof_class_youden = class_youden
    ),
    fold_audit = fold_audit
  )
}

export_plot <- function(plot, stem, width = 8.0, height = 5.6) {
  pdf_file <- file.path(DIRS$figures, paste0(stem, ".pdf"))
  png_file <- file.path(DIRS$figures, paste0(stem, ".png"))
  tiff_file <- file.path(DIRS$figures, paste0(stem, ".tiff"))
  existing <- c(pdf_file, png_file, tiff_file)[file.exists(c(pdf_file, png_file, tiff_file))]
  if (length(existing) && !ALLOW_OVERWRITE) stop("Refusing to overwrite figure(s): ", paste(existing, collapse = ", "))
  if (length(existing)) unlink(existing)
  grDevices::cairo_pdf(pdf_file, width = width, height = height, family = "sans")
  print(plot)
  grDevices::dev.off()
  grDevices::png(png_file, width = width, height = height, units = "in", res = 300, type = "cairo")
  print(plot)
  grDevices::dev.off()
  grDevices::tiff(tiff_file, width = width, height = height, units = "in", res = 600, compression = "lzw")
  print(plot)
  grDevices::dev.off()
  invisible(c(pdf_file, png_file, tiff_file))
}

cat("Reading locked empirical inputs...\n")
labels <- as.data.table(readRDS(INPUT$labels))
X <- as.data.table(readRDS(INPUT$matrix))
final_full <- fread(INPUT$final_full, showProgress = FALSE)
baseline <- fread(INPUT$baseline, showProgress = FALSE)
registry <- fread(INPUT$registry)
prior_oof <- fread(INPUT$prior_oof_28, showProgress = FALSE)
prior_auc <- fread(INPUT$prior_auc_28, showProgress = FALSE)

if (!all(c("stay_id", "cluster_k2") %in% names(labels))) stop("Malformed locked labels")
if (!"stay_id" %in% names(X)) stop("Malformed locked matrix")
if (!identical(labels$stay_id, X$stay_id)) stop("Locked label/matrix row order mismatch")
if (nrow(labels) != EXPECTED_N || sum(labels$cluster_k2 == 1L) != EXPECTED_C1) {
  stop("Locked label denominator/count mismatch")
}
features33 <- setdiff(names(X), "stay_id")
if (length(features33) != 33L || !setequal(features33, registry$feature)) stop("33-feature registry mismatch")

required_final <- c(
  "subject_id", "stay_id", "sofa_score", "respiration", "coagulation", "liver",
  "cardiovascular", "cns", "renal", "mortality_30d", "make30"
)
if (length(setdiff(required_final, names(final_full)))) {
  stop("final_full missing: ", paste(setdiff(required_final, names(final_full)), collapse = ", "))
}
if (!"aki_stage_0_24h" %in% names(baseline)) stop("baseline missing aki_stage_0_24h")
if (anyDuplicated(final_full$stay_id) || anyDuplicated(baseline$stay_id)) stop("Duplicate stay_id in empirical inputs")

idx_final <- match(labels$stay_id, final_full$stay_id)
idx_baseline <- match(labels$stay_id, baseline$stay_id)
if (anyNA(idx_final) || anyNA(idx_baseline)) stop("Locked stay_id absent from empirical input")
clinical <- final_full[idx_final]
kdigo <- baseline$aki_stage_0_24h[idx_baseline]

analysis_data <- data.table(
  subject_id = clinical$subject_id,
  stay_id = labels$stay_id,
  C1 = as.integer(labels$cluster_k2 == 1L),
  kdigo_0_24h = as.numeric(kdigo),
  sofa_renal = as.numeric(clinical$renal),
  sofa_nonrenal = as.numeric(clinical$respiration + clinical$coagulation + clinical$liver + clinical$cardiovascular + clinical$cns),
  sofa_total = as.numeric(clinical$sofa_score),
  mortality_30d = as.numeric(clinical$mortality_30d),
  make30 = as.numeric(clinical$make30)
)
analysis_data <- cbind(analysis_data, X[, ..features33])
if (anyNA(analysis_data$subject_id) || anyDuplicated(analysis_data$subject_id)) {
  stop("Formal patient-level bootstrap requires exactly one analysis row per unique subject_id")
}
if (uniqueN(analysis_data$subject_id) != nrow(analysis_data)) stop("Patient denominator mismatch")
analysis_numeric <- setdiff(names(analysis_data), c("subject_id", "stay_id"))
if (any(!is.finite(as.matrix(analysis_data[, ..analysis_numeric])))) stop("Non-finite analysis value")
if (any(analysis_data$sofa_total != analysis_data$sofa_renal + analysis_data$sofa_nonrenal)) {
  stop("SOFA component sum mismatch")
}

renal_core <- registry[strict_renal_core == 1L, feature]
acid_base <- registry[strict_acid_base == 1L, feature]
renal_acid <- union(renal_core, acid_base)
expanded_renal_acid <- registry[expanded_renal_acid_metabolic == 1L, feature]

if (!identical(sort(renal_core), sort(c("creatinine_max", "bun_max", "urine_output_24h_ml")))) {
  stop("Unexpected strict renal-core registry")
}
if (!identical(sort(acid_base), sort(c("ph_min", "pco2_max", "bicarbonate_min", "aniongap_max")))) {
  stop("Unexpected strict acid-base registry")
}

models <- list(
  list(id = "creatinine", label = "Creatinine maximum", predictors = "creatinine_max", engine = "logistic"),
  list(id = "bun", label = "Blood urea nitrogen maximum", predictors = "bun_max", engine = "logistic"),
  list(id = "creatinine_bun", label = "Creatinine + BUN", predictors = c("creatinine_max", "bun_max"), engine = "logistic"),
  list(id = "kdigo_0_24h", label = "KDIGO stage (0-24 h)", predictors = "kdigo_0_24h", engine = "logistic"),
  list(id = "sofa_nonrenal", label = "Non-renal SOFA", predictors = "sofa_nonrenal", engine = "logistic"),
  list(id = "sofa_renal", label = "SOFA renal component", predictors = "sofa_renal", engine = "logistic"),
  list(id = "sofa_total", label = "Total SOFA", predictors = "sofa_total", engine = "logistic"),
  list(id = "renal_core", label = "Renal core block", predictors = renal_core, engine = "logistic"),
  list(id = "acid_base", label = "Acid-base block", predictors = acid_base, engine = "logistic"),
  list(id = "renal_acid_base", label = "Renal + acid-base block", predictors = renal_acid, engine = "logistic"),
  list(id = "renal_acid_base_expanded", label = "Expanded renal/acid-base/metabolic block", predictors = expanded_renal_acid, engine = "logistic"),
  list(id = "all33_ridge", label = "All 33 features (ridge upper bound)", predictors = features33, engine = "ridge")
)

fold_id <- make_stratified_folds(analysis_data$C1, CV_FOLDS, CV_SEED)
x_all33 <- as.matrix(analysis_data[, ..features33])

cat("Running same-fold OOF models...\n")
model_results <- vector("list", length(models))
names(model_results) <- vapply(models, `[[`, character(1), "id")
checkpoint_path <- file.path(DIRS$checkpoints, paste0("29_wp1_oof_models_checkpoint_", RUN_ID, ".rds"))
for (i in seq_along(models)) {
  model <- models[[i]]
  cat(sprintf("  [%02d/%02d] %s\n", i, length(models), model$id))
  model_results[[model$id]] <- cross_validated_model(analysis_data, x_all33, fold_id, model)
  checkpoint_tmp <- paste0(checkpoint_path, ".tmp")
  saveRDS(model_results, checkpoint_tmp)
  if (file.exists(checkpoint_path)) unlink(checkpoint_path)
  if (!file.rename(checkpoint_tmp, checkpoint_path)) stop("Could not save checkpoint: ", checkpoint_path)
}

model_summary <- rbindlist(lapply(model_results, `[[`, "summary"), fill = TRUE)
threshold_summary <- rbindlist(lapply(model_results, `[[`, "threshold_summary"), fill = TRUE)
predictions <- rbindlist(lapply(model_results, `[[`, "prediction"), fill = TRUE)
fold_audit <- rbindlist(lapply(model_results, `[[`, "fold_audit"), fill = TRUE)

reproduction_models <- c("creatinine", "bun", "creatinine_bun", "kdigo_0_24h", "sofa_nonrenal")
if (anyDuplicated(prior_oof[, .(stay_id, model_id)])) stop("Prior OOF table has duplicate stay/model rows")
current_repro <- predictions[model_id %in% reproduction_models, .(
  stay_id, model_id, current_fold = fold, current_probability = oof_probability
)]
prior_repro <- prior_oof[model_id %in% reproduction_models, .(
  stay_id, model_id, prior_fold = fold, prior_probability = oof_probability
)]
repro_merged <- merge(current_repro, prior_repro, by = c("stay_id", "model_id"), all = TRUE)
if (nrow(repro_merged) != EXPECTED_N * length(reproduction_models) ||
    anyNA(repro_merged$current_fold) || anyNA(repro_merged$prior_fold)) {
  stop("Prior/current OOF reproduction merge is incomplete")
}
same_fold_reproduction <- repro_merged[, .(
  n = .N,
  fold_match_n = sum(current_fold == prior_fold),
  fold_match_rate = mean(current_fold == prior_fold),
  max_abs_probability_difference = max(abs(current_probability - prior_probability)),
  mean_abs_probability_difference = mean(abs(current_probability - prior_probability))
), by = model_id]
same_fold_reproduction <- merge(
  same_fold_reproduction,
  model_summary[model_id %in% reproduction_models, .(model_id, current_auc = auc)],
  by = "model_id", all.x = TRUE
)
same_fold_reproduction <- merge(
  same_fold_reproduction,
  prior_auc[model_id %in% reproduction_models, .(model_id, prior_auc = auc)],
  by = "model_id", all.x = TRUE
)
same_fold_reproduction[, `:=`(
  auc_difference = current_auc - prior_auc,
  reproduction_pass = fold_match_rate == 1 &
    max_abs_probability_difference <= 1e-12 & abs(current_auc - prior_auc) <= 1e-12
)]

all33_auc <- model_summary[model_id == "all33_ridge", auc]
retention_models <- c("creatinine", "bun", "creatinine_bun", "sofa_renal", "sofa_nonrenal", "renal_acid_base")
retention <- model_summary[model_id %in% retention_models, .(
  model_id, model_label, auc, all33_auc = all33_auc,
  discrimination_retention = (auc - 0.5) / (all33_auc - 0.5),
  interpretation = "descriptive only; no pass/fail threshold"
)]

required_pairs <- data.table(
  model_1 = c("creatinine_bun", "creatinine_bun", "creatinine_bun", "renal_acid_base", "sofa_total", "sofa_renal"),
  model_2 = c("sofa_nonrenal", "sofa_total", "sofa_renal", "all33_ridge", "sofa_nonrenal", "sofa_nonrenal")
)
pred_wide <- dcast(predictions, subject_id + stay_id + C1 ~ model_id, value.var = "oof_probability")
pred_wide[, locked_row_order := match(stay_id, analysis_data$stay_id)]
setorder(pred_wide, locked_row_order)
pred_wide[, locked_row_order := NULL]
if (!identical(pred_wide$stay_id, analysis_data$stay_id)) stop("Wide prediction row-order mismatch")
if (anyDuplicated(pred_wide$subject_id) || uniqueN(pred_wide$subject_id) != nrow(pred_wide)) {
  stop("Patient-level prediction table is not one row per subject")
}

pair_labels <- paste(required_pairs$model_1, required_pairs$model_2, sep = "__minus__")
bootstrap_delta <- matrix(NA_real_, nrow = BOOT_B, ncol = nrow(required_pairs), dimnames = list(NULL, pair_labels))
set.seed(BOOT_SEED)
cat("Running paired patient bootstrap (B=", BOOT_B, ")...\n", sep = "")
for (b in seq_len(BOOT_B)) {
  idx <- sample.int(nrow(pred_wide), nrow(pred_wide), replace = TRUE)
  yb <- pred_wide$C1[idx]
  if (length(unique(yb)) == 2L) {
    needed_models <- unique(c(required_pairs$model_1, required_pairs$model_2))
    auc_values <- vapply(needed_models, function(id) fast_auc(yb, pred_wide[[id]][idx]), numeric(1))
    for (j in seq_len(nrow(required_pairs))) {
      bootstrap_delta[b, j] <- auc_values[required_pairs$model_1[j]] - auc_values[required_pairs$model_2[j]]
    }
  }
  if (b %% 100L == 0L) cat(sprintf("  bootstrap %d/%d\n", b, BOOT_B))
}

pairwise <- rbindlist(lapply(seq_len(nrow(required_pairs)), function(j) {
  model_1 <- required_pairs$model_1[j]
  model_2 <- required_pairs$model_2[j]
  values <- bootstrap_delta[, j]
  success <- sum(is.finite(values))
  point <- model_summary[model_id == model_1, auc] - model_summary[model_id == model_2, auc]
  data.table(
    model_1 = model_1,
    model_1_label = model_summary[model_id == model_1, model_label],
    model_2 = model_2,
    model_2_label = model_summary[model_id == model_2, model_label],
    delta_auc_model1_minus_model2 = point,
    delta_auc_ci_low = if (success) unname(quantile(values, 0.025, na.rm = TRUE, type = 6)) else NA_real_,
    delta_auc_ci_high = if (success) unname(quantile(values, 0.975, na.rm = TRUE, type = 6)) else NA_real_,
    bootstrap_B = BOOT_B,
    bootstrap_success_count = success,
    bootstrap_failure_count = BOOT_B - success,
    method = "paired subject-level nonparametric bootstrap (one analysis stay per unique subject); percentile 95% CI"
  )
}))

END_HASHES_PRE_OUTPUT <- vapply(unlist(START_HASH_PATHS, use.names = FALSE), sha256, character(1))
names(END_HASHES_PRE_OUTPUT) <- names(START_HASH_PATHS)
model_required_numeric <- c(
  "auc", "auc_ci_low", "auc_ci_high", "pr_auc_average_precision", "brier",
  "calibration_intercept", "calibration_slope", "sensitivity", "specificity", "ppv", "npv",
  "balanced_accuracy", "accuracy", "kappa", "ari", "predicted_c1_prevalence"
)
threshold_required_numeric <- c(
  "tp", "tn", "fp", "fn", "sensitivity", "specificity", "ppv", "npv",
  "balanced_accuracy", "accuracy", "kappa", "ari",
  "predicted_c1_prevalence"
)
pairwise_required_numeric <- c(
  "delta_auc_model1_minus_model2", "delta_auc_ci_low", "delta_auc_ci_high"
)

qc_checks <- rbindlist(list(
  data.table(
    check = "input_and_script_hashes_unchanged_during_computation",
    observed = paste0(sum(START_HASHES == END_HASHES_PRE_OUTPUT), "/", length(START_HASHES), " unchanged"),
    required = paste0(length(START_HASHES), "/", length(START_HASHES), " unchanged"),
    pass = identical(unname(START_HASHES), unname(END_HASHES_PRE_OUTPUT))
  ),
  data.table(
    check = "locked_denominator_and_unique_patient",
    observed = paste0(nrow(analysis_data), " rows; ", uniqueN(analysis_data$subject_id), " unique subjects"),
    required = paste0(EXPECTED_N, " rows and unique subjects"),
    pass = nrow(analysis_data) == EXPECTED_N && uniqueN(analysis_data$subject_id) == EXPECTED_N
  ),
  data.table(
    check = "same_fold_probability_reproduction_of_script28",
    observed = paste0(sum(same_fold_reproduction$reproduction_pass), "/", nrow(same_fold_reproduction), " models"),
    required = paste0(nrow(same_fold_reproduction), "/", nrow(same_fold_reproduction), " models"),
    pass = all(same_fold_reproduction$reproduction_pass)
  ),
  data.table(
    check = "model_fit_fold_warnings",
    observed = as.character(sum(fold_audit$warning_count)),
    required = "0",
    pass = sum(fold_audit$warning_count) == 0L
  ),
  data.table(
    check = "youden_threshold_audit",
    observed = paste0(sum(fold_audit$youden_status == "ok"), "/", nrow(fold_audit), " folds ok"),
    required = paste0(nrow(fold_audit), "/", nrow(fold_audit), " folds ok"),
    pass = all(fold_audit$youden_status == "ok")
  ),
  data.table(
    check = "calibration_fit_audit",
    observed = paste0(
      sum(model_summary$calibration_status == "ok"), " ok; ",
      sum(model_summary$calibration_status == "warning_finite"), " finite warning; ",
      sum(model_summary$calibration_status == "failed_audit"), " failed"
    ),
    required = "0 failed; finite warnings retained and reported",
    pass = all(model_summary$calibration_status %in% c("ok", "warning_finite"))
  ),
  data.table(
    check = "complete_oof_predictions_and_metrics",
    observed = paste0("prediction rows=", nrow(predictions), "; finite AP=", sum(is.finite(model_summary$pr_auc_average_precision))),
    required = paste0(EXPECTED_N * length(models), " rows; ", length(models), " finite AP"),
    pass = nrow(predictions) == EXPECTED_N * length(models) &&
      all(is.finite(predictions$oof_probability)) && all(is.finite(model_summary$pr_auc_average_precision))
  ),
  data.table(
    check = "all_required_model_threshold_and_contrast_metrics_finite",
    observed = paste0(
      "model=", sum(is.finite(as.matrix(model_summary[, ..model_required_numeric]))), "/", nrow(model_summary) * length(model_required_numeric),
      "; threshold=", sum(is.finite(as.matrix(threshold_summary[, ..threshold_required_numeric]))), "/", nrow(threshold_summary) * length(threshold_required_numeric),
      "; contrast=", sum(is.finite(as.matrix(pairwise[, ..pairwise_required_numeric]))), "/", nrow(pairwise) * length(pairwise_required_numeric)
    ),
    required = "all finite",
    pass = all(is.finite(as.matrix(model_summary[, ..model_required_numeric]))) &&
      all(is.finite(as.matrix(threshold_summary[, ..threshold_required_numeric]))) &&
      all(is.finite(as.matrix(pairwise[, ..pairwise_required_numeric]))) &&
      all(is.finite(retention$discrimination_retention))
  ),
  data.table(
    check = "paired_bootstrap_success",
    observed = paste0(min(pairwise$bootstrap_success_count), "-", max(pairwise$bootstrap_success_count), "/", BOOT_B),
    required = paste0(BOOT_B, "/", BOOT_B, " for every contrast"),
    pass = all(pairwise$bootstrap_success_count == BOOT_B)
  )
), fill = TRUE)
if (any(!qc_checks$pass)) {
  failed <- qc_checks[pass == FALSE, paste0(check, ": observed ", observed, "; required ", required)]
  failure_dir <- file.path(DIRS$logs, "failed_runs")
  dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
  failure_csv <- file.path(failure_dir, paste0("29_WP1_QC_failure_", RUN_ID, ".csv"))
  failure_txt <- file.path(failure_dir, paste0("29_WP1_failure_", RUN_ID, ".txt"))
  fwrite(qc_checks, failure_csv, bom = TRUE)
  writeLines(c(
    paste0("run_id=", RUN_ID),
    paste0("failed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "WP1 QC gate failed before fixed-name output creation.",
    failed
  ), failure_txt, useBytes = TRUE)
  stop("WP1 QC gate failed before final output creation:\n", paste(failed, collapse = "\n"))
}

write_csv_atomic(model_summary, file.path(DIRS$tables, "Table_D1_identity_models_oof.csv"))
write_csv_atomic(pairwise, file.path(DIRS$tables, "Table_D1_identity_pairwise_auc.csv"))
write_csv_atomic(threshold_summary, file.path(DIRS$tables, "Table_D1_identity_threshold_metrics.csv"))
write_csv_atomic(retention, file.path(DIRS$tables, "Table_D1_identity_discrimination_retention.csv"))
write_csv_atomic(predictions, file.path(DIRS$tables, "29_Table_D1_identity_oof_predictions.csv"))
write_csv_atomic(fold_audit, file.path(DIRS$logs, "29_cross_validation_fold_audit.csv"))
write_csv_atomic(same_fold_reproduction, file.path(DIRS$logs, "29_same_fold_reproduction_audit.csv"))
write_csv_atomic(qc_checks, file.path(DIRS$logs, "29_WP1_QC_checks.csv"))

plot_data <- copy(model_summary)
plot_data[, model_label := factor(model_label, levels = rev(model_label[order(auc)]))]
plot <- ggplot(plot_data, aes(x = auc, y = model_label)) +
  geom_errorbar(aes(xmin = auc_ci_low, xmax = auc_ci_high), orientation = "y", width = 0, linewidth = 0.5, colour = "grey55") +
  geom_point(size = 2.4, colour = "#2E2E2E") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey70", linewidth = 0.45) +
  labs(
    x = "Out-of-fold AUC for locked C1 membership",
    y = NULL,
    title = "Reconstruction of the locked empirical C1 label",
    subtitle = "All models used identical folds; all-33 ridge is a reconstruction upper bound"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey35"),
    axis.title.x = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.margin = margin(8, 14, 8, 8)
  )
export_plot(plot, "Figure_D1_identity_model_comparison")

provenance <- data.table(
  source = names(START_HASH_PATHS),
  path = unlist(START_HASH_PATHS, use.names = FALSE),
  sha256 = unname(START_HASHES),
  read_only = TRUE,
  run_id = RUN_ID,
  run_mode = RUN_MODE
)
write_csv_atomic(provenance, file.path(DIRS$provenance, "29_WP1_input_checksums.csv"))

run_metadata <- data.table(
  field = c(
    "run_id", "run_mode", "script_path", "script_sha256", "internal_reference_sha256",
    "n_analysis_rows", "n_unique_subjects", "n_c1", "cv_folds", "cv_seed",
    "primary_threshold", "sensitivity_threshold", "bootstrap_B", "bootstrap_seed",
    "bootstrap_unit", "checkpoint_path", "allow_overwrite"
  ),
  value = c(
    RUN_ID, RUN_MODE, SCRIPT_PATH, START_HASHES[["current_script_29"]], NA_character_,
    nrow(analysis_data), uniqueN(analysis_data$subject_id), sum(analysis_data$C1), CV_FOLDS, CV_SEED,
    PRIMARY_THRESHOLD, SENSITIVITY_THRESHOLD, BOOT_B, BOOT_SEED,
    "unique subject_id (one analysis stay per patient)", checkpoint_path, ALLOW_OVERWRITE
  )
)
write_csv_atomic(run_metadata, file.path(DIRS$provenance, "29_WP1_run_metadata.csv"))

key <- model_summary[model_id %in% c(
  "creatinine_bun", "sofa_renal", "sofa_nonrenal", "sofa_total",
  "renal_core", "acid_base", "renal_acid_base", "renal_acid_base_expanded", "all33_ridge"
), .(model_id, auc, auc_ci_low, auc_ci_high, pr_auc_average_precision, brier,
     calibration_intercept, calibration_slope, balanced_accuracy, kappa, ari,
     predicted_c1_prevalence)]

log_lines <- c(
  "# 29 Domain 1 renal-identity completion log",
  "",
  paste0("- Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- Run ID: ", RUN_ID, "; mode: ", RUN_MODE, "."),
  paste0("- Locked denominator: ", EXPECTED_N, "; C1: ", EXPECTED_C1, "."),
  paste0("- Unique subjects: ", uniqueN(analysis_data$subject_id), "; one analysis stay per patient: TRUE."),
  paste0("- Cross-validation: ", CV_FOLDS, " stratified folds; seed ", CV_SEED, "."),
  paste0("- Primary threshold: closest achievable training-fold prevalence under tied probabilities; exact matching is not claimed when ties prevent it."),
  paste0("- Sensitivity threshold: training-fold Youden."),
  paste0("- Pairwise bootstrap: B=", BOOT_B, "; seed ", BOOT_SEED, "; resampling unit is unique subject_id."),
  "- All-33 model: nested ridge logistic regression; lambda selected inside each outer training fold.",
  "- PR-AUC is reported as tie-invariant grouped average precision (step-function precision-recall area).",
  "- Calibration warnings are retained verbatim; a finite converged fit with a numerical warning is distinguished from error, non-finite output, or non-convergence.",
  "- Strict domain registry: renal core 3 variables; acid-base 4 variables.",
  "- Expanded sensitivity block additionally includes lactate_max and potassium_max; it is not relabeled as strict acid-base.",
  "- Outcome is algorithmic C1 membership, not mortality, MAKE30, or a causal endpoint.",
  "- All formal QC gates passed before fixed-name outputs and completion marker were written.",
  "",
  "## QC gates",
  "",
  paste(capture.output(print(qc_checks)), collapse = "\n"),
  "",
  "## Same-fold reproduction of script 28",
  "",
  paste(capture.output(print(same_fold_reproduction)), collapse = "\n"),
  "",
  "## Key results",
  "",
  paste(capture.output(print(key)), collapse = "\n"),
  "",
  "## Required paired AUC contrasts",
  "",
  paste(capture.output(print(pairwise)), collapse = "\n")
)
write_lines_atomic(log_lines, file.path(DIRS$logs, "29_domain1_renal_identity_completion_log.md"))
write_lines_atomic(capture.output(sessionInfo()), file.path(DIRS$logs, "29_sessionInfo.txt"))

manifest_targets <- setdiff(FINAL_OUTPUTS, c(OUTPUT_MANIFEST, RUN_COMPLETED))
missing_manifest_targets <- manifest_targets[!file.exists(manifest_targets)]
if (length(missing_manifest_targets)) {
  stop("Cannot create output manifest; missing fixed output(s):\n", paste(missing_manifest_targets, collapse = "\n"))
}
output_manifest <- data.table(
  file = basename(manifest_targets),
  path = manifest_targets,
  bytes = file.info(manifest_targets)$size,
  sha256 = vapply(manifest_targets, sha256, character(1)),
  run_id = RUN_ID,
  run_mode = RUN_MODE
)
write_csv_atomic(output_manifest, OUTPUT_MANIFEST)

completion_lines <- c(
  "RUN COMPLETED",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", START_HASHES[["current_script_29"]]),
  "internal_provenance_reference_not_distributed=true",
  paste0("output_manifest_sha256=", sha256(OUTPUT_MANIFEST))
)
completion_tmp <- paste0(RUN_COMPLETED, ".", RUN_ID, ".tmp")
writeLines(completion_lines, completion_tmp, useBytes = TRUE)
if (!file.rename(completion_tmp, RUN_COMPLETED)) stop("Could not atomically create completion marker")
if (!file.exists(RUN_IN_PROGRESS) || !unlink(RUN_IN_PROGRESS)) {
  stop("Completion marker exists but in-progress marker could not be cleared: ", RUN_IN_PROGRESS)
}

cat("WP1 completed.\n")
print(key)
print(pairwise)
