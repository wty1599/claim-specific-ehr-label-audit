## =====================================================================
## 29d_domain1_external_renal_identity_eicu.R
## Formal WP5: external clinical-identity audit of the frozen eICU label.
##
## The eICU label is reproduced from the frozen MIMIC preprocessing recipe
## and frozen MIMIC centroids. No eICU reclustering is performed.
## Evaluation uses patient-grouped 10-fold OOF prediction and a paired
## patient-cluster bootstrap for external AUC confidence intervals.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "pROC", "mclust", "ggplot2", "digest")
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
  provenance = file.path(OUT_ROOT, "provenance")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

SCRIPT_PATH <- normalizePath(
  file.path(OUT_ROOT, "scripts", "29d_domain1_external_renal_identity_eicu.R"),
  winslash = "/", mustWork = TRUE
)
RUN_MODE <- "formal"
RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid())
ALLOW_OVERWRITE <- identical(tolower(Sys.getenv("D1_ALLOW_OVERWRITE", "false")), "true")
RUN_IN_PROGRESS <- file.path(OUT_ROOT, "29d_run_in_progress.flag")
RUN_COMPLETED <- file.path(OUT_ROOT, "29d_run_completed.ok")
OUTPUT_MANIFEST <- file.path(DIRS$provenance, "29d_WP5_output_sha256_manifest.csv")

INPUT <- list(
  eicu_data = file.path(PROJECT_ROOT, "data", "eicu_external.csv"),
  eicu_labels = file.path(PROJECT_ROOT, "output", "eicu_cluster_labels_mice.csv"),
  frozen_recipe = file.path(PROJECT_ROOT, "output", "qc", "external_assignment_preprocessing_recipe_mice_labels.csv"),
  frozen_centroids = file.path(PROJECT_ROOT, "output", "qc", "external_assignment_frozen_centroids_mice_labels.csv"),
  frozen_assignment_object = file.path(PROJECT_ROOT, "output", "qc", "headtohead_all_mice_labels.rds"),
  assignment_input_fingerprint = file.path(PROJECT_ROOT, "output", "qc", "external_assignment_input_fingerprint_mice_labels.csv"),
  assignment_provenance = file.path(PROJECT_ROOT, "output", "qc", "external_assignment_provenance_mice_labels.csv"),
  assignment_feature_audit = file.path(PROJECT_ROOT, "output", "qc", "external_assignment_feature_audit_mice_labels.csv"),
  assignment_source_script = file.path(PROJECT_ROOT, "r", "19_external_headtohead_mice_labels_v2.R"),
  internal_wp1_models = file.path(DIRS$tables, "Table_D1_identity_models_oof.csv"),
  internal_wp1_manifest = file.path(DIRS$provenance, "29_WP1_output_sha256_manifest.csv")
)
if (any(!file.exists(unlist(INPUT)))) {
  stop("Missing required input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))
}

EXPECTED_STAYS <- 17465L
EXPECTED_PATIENTS <- 16212L
EXPECTED_HOSPITALS <- 199L
EXPECTED_C1 <- 3462L
EXPECTED_ASSIGNMENT_RUN_ID <- "f41396c7aa239dfb"
EXPECTED_COMMON_FEATURE_N <- 33L
CV_FOLDS <- 10L
CV_SEED <- 2026071728L
BOOT_B <- 2000L
BOOT_SEED <- 202607295L
PRIMARY_THRESHOLD <- "closest_achievable_training_prevalence_match_under_ties"
SENSITIVITY_THRESHOLD <- "training_youden"

FINAL_OUTPUTS <- c(
  file.path(DIRS$tables, "Table_D1_external_identity_models.csv"),
  file.path(DIRS$tables, "Table_D1_external_KDIGO_by_cluster.csv"),
  file.path(DIRS$tables, "Table_D1_external_AKI_criterion_by_cluster.csv"),
  file.path(DIRS$tables, "Table_D1_internal_external_identity_comparison.csv"),
  file.path(DIRS$tables, "29d_Table_D1_external_KDIGO_association_summary.csv"),
  file.path(DIRS$tables, "29d_Table_D1_external_identity_threshold_metrics.csv"),
  file.path(DIRS$tables, "29d_Table_D1_external_identity_pairwise_auc.csv"),
  file.path(DIRS$tables, "29d_Table_D1_external_identity_oof_predictions.csv"),
  file.path(DIRS$logs, "29d_external_identity_fold_audit.csv"),
  file.path(DIRS$logs, "29d_external_identity_patient_fold_map.csv"),
  file.path(DIRS$logs, "29d_external_identity_QC.csv"),
  file.path(DIRS$logs, "29d_external_identity_harmonization.md"),
  file.path(DIRS$logs, "29d_external_identity_log.md"),
  file.path(DIRS$logs, "29d_sessionInfo.txt"),
  file.path(DIRS$provenance, "29d_WP5_input_checksums.csv"),
  file.path(DIRS$provenance, "29d_WP5_run_metadata.csv"),
  OUTPUT_MANIFEST,
  file.path(DIRS$figures, paste0("Figure_D1_internal_external_identity_AUC.", c("pdf", "png", "tiff"))),
  RUN_COMPLETED
)
existing_outputs <- FINAL_OUTPUTS[file.exists(FINAL_OUTPUTS)]
if (length(existing_outputs) && !ALLOW_OVERWRITE) {
  stop("Refusing to overwrite existing WP5 output(s):\n", paste(existing_outputs, collapse = "\n"))
}
if (file.exists(RUN_IN_PROGRESS)) stop("Existing in-progress marker requires review: ", RUN_IN_PROGRESS)

write_csv_atomic <- function(x, path) {
  if (file.exists(path) && !ALLOW_OVERWRITE) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  if (file.exists(path) && !ALLOW_OVERWRITE) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write: ", path)
  invisible(path)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)
clamp_probability <- function(p) pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)

START_HASH_PATHS <- c(INPUT, current_script_29d = SCRIPT_PATH)
START_HASHES <- vapply(unlist(START_HASH_PATHS, use.names = FALSE), sha256, character(1))
names(START_HASHES) <- names(START_HASH_PATHS)
write_lines_atomic(c(
  paste0("run_id=", RUN_ID), paste0("run_mode=", RUN_MODE),
  paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256_at_start=", START_HASHES[["current_script_29d"]]),
  "internal_provenance_reference_not_distributed=true"
), RUN_IN_PROGRESS)

prevalence_threshold <- function(probability, prevalence) {
  probability <- as.numeric(probability)
  if (!length(probability) || any(!is.finite(probability))) stop("Invalid prevalence-threshold probabilities")
  n <- length(probability)
  target_n <- round(prevalence * n)
  groups <- data.table(probability = probability)[, .(boundary_tie_n = .N), by = probability]
  setorder(groups, -probability)
  groups[, achieved_positive_n := cumsum(boundary_tie_n)]
  candidates <- rbindlist(list(
    data.table(threshold = Inf, achieved_positive_n = 0L, boundary_tie_n = 0L),
    groups[, .(threshold = probability, achieved_positive_n, boundary_tie_n)]
  ))
  candidates[, `:=`(
    absolute_count_deviation = abs(achieved_positive_n - target_n),
    achieved_prevalence = achieved_positive_n / n,
    absolute_prevalence_deviation = abs(achieved_positive_n / n - prevalence),
    over_target = achieved_positive_n > target_n
  )]
  chosen <- candidates[order(absolute_count_deviation, over_target, -threshold)][1L]
  list(
    threshold = chosen$threshold,
    target_positive_n = target_n,
    achieved_positive_n = chosen$achieved_positive_n,
    target_prevalence = prevalence,
    achieved_prevalence = chosen$achieved_prevalence,
    absolute_count_deviation = chosen$absolute_count_deviation,
    absolute_prevalence_deviation = chosen$absolute_prevalence_deviation,
    exact_target_count_achievable = chosen$achieved_positive_n == target_n,
    boundary_tie_n = chosen$boundary_tie_n
  )
}

youden_threshold <- function(y, probability) {
  warnings <- character()
  error_message <- ""
  out <- tryCatch(
    withCallingHandlers({
      roc_obj <- pROC::roc(y, probability, levels = c(0, 1), direction = "<", quiet = TRUE)
      pROC::coords(roc_obj, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }),
    error = function(e) { error_message <<- conditionMessage(e); NULL }
  )
  thresholds <- if (is.null(out)) numeric() else as.numeric(out$threshold)
  finite_thresholds <- thresholds[is.finite(thresholds)]
  status <- if (nzchar(error_message)) "error" else if (!length(finite_thresholds)) {
    "nonfinite"
  } else if (length(warnings)) "warning" else "ok"
  list(
    threshold = if (length(finite_thresholds)) max(finite_thresholds) else NA_real_,
    status = status,
    candidate_n = length(finite_thresholds),
    message = paste(c(error_message, unique(warnings)), collapse = " | ")
  )
}

average_precision <- function(y, probability) {
  positives <- sum(y == 1L)
  if (!positives || positives == length(y) || any(!is.finite(probability))) return(NA_real_)
  groups <- data.table(y = as.integer(y), probability = as.numeric(probability))[
    , .(n = .N, positives = sum(y == 1L)), by = probability
  ]
  setorder(groups, -probability)
  groups[, `:=`(cumulative_n = cumsum(n), cumulative_positives = cumsum(positives))]
  groups[, `:=`(
    precision = cumulative_positives / cumulative_n,
    recall = cumulative_positives / sum(positives)
  )]
  sum(c(groups$recall[1L], diff(groups$recall)) * groups$precision)
}

binary_metrics <- function(y, predicted_class) {
  tp <- sum(y == 1L & predicted_class == 1L)
  tn <- sum(y == 0L & predicted_class == 0L)
  fp <- sum(y == 0L & predicted_class == 1L)
  fn <- sum(y == 1L & predicted_class == 0L)
  sensitivity <- safe_div(tp, tp + fn)
  specificity <- safe_div(tn, tn + fp)
  observed <- (tp + tn) / length(y)
  py <- mean(y == 1L); pp <- mean(predicted_class == 1L)
  expected <- py * pp + (1 - py) * (1 - pp)
  data.table(
    tp, tn, fp, fn, sensitivity, specificity,
    ppv = safe_div(tp, tp + fp), npv = safe_div(tn, tn + fn),
    balanced_accuracy = mean(c(sensitivity, specificity)), accuracy = observed,
    kappa = if (abs(1 - expected) < 1e-12) NA_real_ else (observed - expected) / (1 - expected),
    ari = mclust::adjustedRandIndex(y, predicted_class), predicted_c1_prevalence = pp
  )
}

calibration_metrics <- function(y, probability) {
  lp <- qlogis(clamp_probability(probability))
  audit_glm <- function(formula, offset = NULL) {
    warnings <- character(); error_message <- ""
    args <- list(formula = formula, family = binomial())
    if (!is.null(offset)) args$offset <- offset
    fit <- tryCatch(
      withCallingHandlers(do.call(glm, args), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
      error = function(e) { error_message <<- conditionMessage(e); NULL }
    )
    status <- if (nzchar(error_message)) "error" else if (is.null(fit) || any(!is.finite(coef(fit)))) {
      "nonfinite"
    } else if (!isTRUE(fit$converged)) "not_converged" else if (length(warnings)) "warning" else "ok"
    list(fit = fit, status = status, message = paste(c(error_message, unique(warnings)), collapse = " | "))
  }
  ia <- audit_glm(y ~ 1, offset = lp); sa <- audit_glm(y ~ lp)
  overall_status <- if (ia$status == "ok" && sa$status == "ok") "ok" else if (
    ia$status %in% c("ok", "warning") && sa$status %in% c("ok", "warning")
  ) "warning_finite" else "failed_audit"
  data.table(
    calibration_intercept = if (is.null(ia$fit)) NA_real_ else unname(coef(ia$fit)[1L]),
    calibration_slope = if (is.null(sa$fit) || length(coef(sa$fit)) < 2L) NA_real_ else unname(coef(sa$fit)[2L]),
    calibration_intercept_status = ia$status,
    calibration_intercept_message = ia$message,
    calibration_slope_status = sa$status,
    calibration_slope_message = sa$message,
    calibration_status = overall_status
  )
}

fast_auc <- function(y, probability) {
  n1 <- sum(y == 1L); n0 <- length(y) - n1
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_patient_grouped_folds <- function(data, k, seed) {
  patient <- data[, .(
    n_stays = .N,
    c1_stays = sum(C1),
    c2_stays = sum(C1 == 0L)
  ), by = uniquepid]
  set.seed(seed)
  patient[, random_order := runif(.N)]
  setorder(patient, -n_stays, -c1_stays, random_order)

  target <- c(
    stays = sum(patient$n_stays) / k,
    c1 = sum(patient$c1_stays) / k,
    c2 = sum(patient$c2_stays) / k,
    patients = nrow(patient) / k
  )
  totals <- data.table(
    fold = seq_len(k), stays = 0, c1 = 0, c2 = 0, patients = 0
  )
  fold_priority <- sample(seq_len(k))
  patient[, fold := NA_integer_]
  for (i in seq_len(nrow(patient))) {
    candidate <- copy(totals)
    candidate[, score :=
      (stays + patient$n_stays[i]) / target[["stays"]] +
      (c1 + patient$c1_stays[i]) / target[["c1"]] +
      (c2 + patient$c2_stays[i]) / target[["c2"]] +
      0.25 * (patients + 1) / target[["patients"]]
    ]
    candidate[, tie_rank := match(fold, fold_priority)]
    chosen <- candidate[order(score, tie_rank), fold][1L]
    patient$fold[i] <- chosen
    totals[fold == chosen, `:=`(
      stays = stays + patient$n_stays[i],
      c1 = c1 + patient$c1_stays[i],
      c2 = c2 + patient$c2_stays[i],
      patients = patients + 1
    )]
  }
  fold <- patient$fold[match(data$uniquepid, patient$uniquepid)]
  if (anyNA(fold)) stop("Patient fold mapping failed")
  if (any(totals$c1 == 0L | totals$c2 == 0L)) stop("At least one grouped fold lacks C1 or C2 stays")
  list(
    fold = fold,
    patient_map = patient[, .(uniquepid, fold, n_stays, c1_stays, c2_stays)],
    fold_totals = totals
  )
}

fit_oof_model <- function(data, model, fold_id) {
  n <- nrow(data)
  probability <- rep(NA_real_, n)
  class_primary <- rep(NA_integer_, n)
  class_youden <- rep(NA_integer_, n)
  audits <- vector("list", CV_FOLDS)
  for (fold in seq_len(CV_FOLDS)) {
    train <- which(fold_id != fold); test <- which(fold_id == fold)
    warnings <- character()
    formula <- reformulate(model$predictors, response = "C1")
    fit <- withCallingHandlers(
      glm(formula, data = data[train], family = binomial()),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    fit_converged <- isTRUE(fit$converged)
    coefficient_finite <- all(is.finite(coef(fit)))
    if (!fit_converged || !coefficient_finite) {
      stop("Model did not converge with finite coefficients: ", model$id, " fold ", fold)
    }
    p_train <- as.numeric(predict(fit, newdata = data[train], type = "response"))
    p_test <- as.numeric(predict(fit, newdata = data[test], type = "response"))
    if (any(!is.finite(c(p_train, p_test)))) stop("Non-finite prediction: ", model$id, " fold ", fold)
    primary <- prevalence_threshold(p_train, mean(data$C1[train]))
    youden <- youden_threshold(data$C1[train], p_train)
    probability[test] <- p_test
    class_primary[test] <- as.integer(p_test >= primary$threshold)
    class_youden[test] <- if (is.finite(youden$threshold)) as.integer(p_test >= youden$threshold) else NA_integer_
    audits[[fold]] <- data.table(
      model_id = model$id, fold, n_train_stays = length(train), n_test_stays = length(test),
      n_train_patients = uniqueN(data$uniquepid[train]), n_test_patients = uniqueN(data$uniquepid[test]),
      patient_overlap_n = length(intersect(unique(data$uniquepid[train]), unique(data$uniquepid[test]))),
      train_has_both_classes = uniqueN(data$C1[train]) == 2L,
      test_has_both_classes = uniqueN(data$C1[test]) == 2L,
      fit_converged = fit_converged, coefficient_finite = coefficient_finite,
      train_c1_prevalence = mean(data$C1[train]), test_c1_prevalence = mean(data$C1[test]),
      threshold_primary = primary$threshold,
      target_train_positive_n_primary = primary$target_positive_n,
      achieved_train_positive_n_primary = primary$achieved_positive_n,
      absolute_train_count_deviation_primary = primary$absolute_count_deviation,
      absolute_train_prevalence_deviation_primary = primary$absolute_prevalence_deviation,
      exact_target_count_achievable_primary = primary$exact_target_count_achievable,
      threshold_youden = youden$threshold, youden_status = youden$status,
      youden_candidate_n = youden$candidate_n, youden_message = youden$message,
      warning_count = length(warnings), warning_text = paste(unique(warnings), collapse = " | ")
    )
  }
  if (anyNA(probability) || anyNA(class_primary) || anyNA(class_youden)) stop("Incomplete OOF result: ", model$id)
  fold_audit <- rbindlist(audits)
  cal <- calibration_metrics(data$C1, probability)
  primary_metrics <- binary_metrics(data$C1, class_primary)
  youden_metrics <- binary_metrics(data$C1, class_youden)
  summary <- cbind(data.table(
    model_id = model$id, model_label = model$label,
    predictors = paste(model$predictors, collapse = " + "),
    analysis_status = "SUPPORTED_FORMAL",
    evaluation = "10-fold out-of-fold; folds grouped by uniquepid",
    statistical_unit = "ICU stay", clustered_unit = "uniquepid",
    n_stays = n, n_patients = uniqueN(data$uniquepid), c1_n = sum(data$C1),
    c1_prevalence = mean(data$C1), auc = fast_auc(data$C1, probability),
    pr_auc_average_precision = average_precision(data$C1, probability),
    brier = mean((probability - data$C1)^2),
    fold_warning_n = sum(fold_audit$warning_count),
    youden_failed_fold_n = sum(fold_audit$youden_status != "ok")
  ), cal, data.table(threshold_rule = PRIMARY_THRESHOLD), primary_metrics)
  threshold_summary <- rbindlist(list(
    cbind(data.table(model_id = model$id, threshold_rule = PRIMARY_THRESHOLD), primary_metrics),
    cbind(data.table(model_id = model$id, threshold_rule = SENSITIVITY_THRESHOLD), youden_metrics)
  ), fill = TRUE)
  list(
    summary = summary, threshold_summary = threshold_summary,
    prediction = data.table(
      patientunitstayid = data$patientunitstayid, uniquepid = data$uniquepid,
      fold = fold_id, C1 = data$C1, model_id = model$id,
      oof_probability = probability, oof_class_primary = class_primary,
      oof_class_youden = class_youden
    ), fold_audit = fold_audit
  )
}

sqdist <- function(A, B) {
  outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
}

export_plot <- function(plot, stem, width = 8.2, height = 5.6) {
  paths <- file.path(DIRS$figures, paste0(stem, c(".pdf", ".png", ".tiff")))
  if (any(file.exists(paths)) && !ALLOW_OVERWRITE) stop("Refusing to overwrite figure")
  if (any(file.exists(paths))) unlink(paths[file.exists(paths)])
  cairo_pdf(paths[1], width = width, height = height, family = "sans"); print(plot); dev.off()
  png(paths[2], width = width, height = height, units = "in", res = 300, type = "cairo"); print(plot); dev.off()
  tiff(paths[3], width = width, height = height, units = "in", res = 600, compression = "lzw"); print(plot); dev.off()
  invisible(paths)
}

cat("Reading current frozen eICU assignment artifacts...\n")
E <- fread(INPUT$eicu_data, showProgress = FALSE)
L <- fread(INPUT$eicu_labels, showProgress = FALSE)
recipe <- fread(INPUT$frozen_recipe)
centroids <- fread(INPUT$frozen_centroids)
frozen_object <- readRDS(INPUT$frozen_assignment_object)
assignment_input_fingerprint <- fread(INPUT$assignment_input_fingerprint)
assignment_provenance <- fread(INPUT$assignment_provenance)
assignment_feature_audit <- fread(INPUT$assignment_feature_audit)
internal <- fread(INPUT$internal_wp1_models)

required_eicu <- c(
  "patientunitstayid", "uniquepid", "hospitalid", "creatinine_max", "bun_max",
  "scr_peak_0_24h", "aki_stage_0_24h", "sofa"
)
if (length(setdiff(required_eicu, names(E)))) stop("eICU data missing: ", paste(setdiff(required_eicu, names(E)), collapse = ", "))
if (anyDuplicated(E$patientunitstayid) || anyDuplicated(L$patientunitstayid)) stop("Duplicate eICU stay ID")
if (nrow(E) != EXPECTED_STAYS || uniqueN(E$uniquepid) != EXPECTED_PATIENTS) stop("Current eICU denominator mismatch")
if (uniqueN(E$hospitalid) != EXPECTED_HOSPITALS) stop("Current eICU hospital count mismatch")
if (nrow(L) != EXPECTED_STAYS || uniqueN(L$uniquepid) != EXPECTED_PATIENTS) stop("Current frozen-label denominator mismatch")

idx <- match(E$patientunitstayid, L$patientunitstayid)
if (anyNA(idx)) stop("Frozen labels do not cover all eICU stays")
L <- L[idx]
if (!identical(as.character(E$uniquepid), as.character(L$uniquepid))) stop("eICU patient/stay alignment mismatch")
if (sum(L$cluster_k2 == 1L) != EXPECTED_C1 || any(L$cluster_k2_display != L$cluster_k2)) {
  stop("Current 3,462-label direction/count mismatch; possible obsolete 3,417-label artifact")
}
if (uniqueN(L$assignment_run_id) != 1L || unique(L$assignment_run_id) != EXPECTED_ASSIGNMENT_RUN_ID) stop("Assignment run ID mismatch")
if (uniqueN(L$common_feature_n) != 1L || unique(L$common_feature_n) != EXPECTED_COMMON_FEATURE_N) stop("Frozen dimension mismatch")
if (any(grepl("re-cluster", L$assignment_rule, ignore.case = TRUE) & !grepl("no eICU re-clustering", L$assignment_rule, fixed = TRUE))) {
  stop("Ambiguous external assignment rule")
}

features <- recipe$feature
if (nrow(recipe) != EXPECTED_COMMON_FEATURE_N || anyDuplicated(features)) stop("Malformed frozen recipe")
if (length(setdiff(features, names(E)))) stop("eICU missing frozen feature(s): ", paste(setdiff(features, names(E)), collapse = ", "))
centroid_features <- setdiff(names(centroids), "source_label")
if (!setequal(features, centroid_features) || nrow(centroids) != 2L) stop("Malformed frozen centroid table")

if (nrow(assignment_provenance) != 1L || assignment_provenance$assignment_run_id != EXPECTED_ASSIGNMENT_RUN_ID) {
  stop("Frozen assignment provenance mismatch")
}
if (!identical(as.character(frozen_object$assignment_run_id), EXPECTED_ASSIGNMENT_RUN_ID) ||
    !identical(as.integer(frozen_object$common_feature_n), EXPECTED_COMMON_FEATURE_N) ||
    isTRUE(frozen_object$eicu_reclustered) ||
    !identical(as.character(frozen_object$high_risk_label), "1")) {
  stop("Authoritative frozen assignment object metadata mismatch")
}
if (!identical(as.character(frozen_object$common), as.character(features))) {
  stop("Frozen feature order differs between the authoritative object and recipe CSV")
}
object_recipe <- as.data.table(frozen_object$preprocessing_recipe_table)
recipe_numeric_columns <- setdiff(names(recipe), "feature")
if (!identical(names(object_recipe), names(recipe)) ||
    !identical(as.character(object_recipe$feature), as.character(recipe$feature)) ||
    max(abs(
      as.matrix(object_recipe[, ..recipe_numeric_columns]) -
        as.matrix(recipe[, ..recipe_numeric_columns])
    )) > 1e-10) {
  stop("Frozen preprocessing recipe CSV does not reproduce the authoritative object")
}
object_centroids <- frozen_object$frozen_centroids[, features, drop = FALSE]
csv_centroids <- as.matrix(centroids[match(rownames(object_centroids), as.character(source_label)), ..features])
if (anyNA(csv_centroids) || max(abs(object_centroids - csv_centroids)) > 1e-12) {
  stop("Frozen centroid CSV does not reproduce the authoritative object")
}
hash_fields <- c("common_feature_hash", "preprocessing_hash", "centroid_hash")
for (field in hash_fields) {
  if (uniqueN(L[[field]]) != 1L ||
      !identical(as.character(unique(L[[field]])), as.character(assignment_provenance[[field]])) ||
      !identical(as.character(unique(L[[field]])), as.character(frozen_object[[field]]))) {
    stop("Frozen assignment hash mismatch: ", field)
  }
}
if (!identical(unname(tools::md5sum(INPUT$eicu_labels)), as.character(frozen_object$eicu_label_file_md5)) ||
    !identical(as.character(assignment_provenance$label_file_md5), as.character(frozen_object$eicu_label_file_md5))) {
  stop("Frozen eICU label-file MD5 mismatch")
}
current_eicu_md5 <- unname(tools::md5sum(INPUT$eicu_data))
object_eicu_fingerprint <- as.data.table(frozen_object$input_fingerprint)[role == "eicu_external"]
csv_eicu_fingerprint <- assignment_input_fingerprint[role == "eicu_external"]
if (nrow(object_eicu_fingerprint) != 1L || nrow(csv_eicu_fingerprint) != 1L ||
    !identical(as.character(current_eicu_md5), as.character(object_eicu_fingerprint$md5)) ||
    !identical(as.character(current_eicu_md5), as.character(csv_eicu_fingerprint$md5))) {
  stop("Current eICU input does not match the authoritative frozen-assignment fingerprint")
}
if (nrow(assignment_feature_audit) != EXPECTED_COMMON_FEATURE_N ||
    !all(assignment_feature_audit$used_for_assignment) ||
    !setequal(assignment_feature_audit$feature, features)) {
  stop("Frozen assignment feature audit mismatch")
}

Xe <- matrix(NA_real_, nrow(E), length(features), dimnames = list(NULL, features))
raw_missing <- integer(length(features)); names(raw_missing) <- features
for (j in features) {
  r <- recipe[feature == j]
  value <- suppressWarnings(as.numeric(E[[j]]))
  raw_missing[j] <- sum(!is.finite(value))
  value[value < r$winsor_p01] <- r$winsor_p01
  value[value > r$winsor_p99] <- r$winsor_p99
  value[!is.finite(value)] <- r$imputation_median
  Xe[, j] <- (value - r$reference_mean) / r$reference_sd
}
if (any(!is.finite(Xe))) stop("Frozen preprocessing produced non-finite value(s)")

centroid_matrix <- as.matrix(centroids[, ..features])
rownames(centroid_matrix) <- as.character(centroids$source_label)
reproduced_label <- as.integer(rownames(centroid_matrix)[max.col(-sqdist(Xe, centroid_matrix), ties.method = "first")])
assignment_agreement <- mean(reproduced_label == L$cluster_k2)
if (assignment_agreement != 1) stop("Frozen recipe/centroid assignment does not exactly reproduce current labels")

analysis_data <- data.table(
  patientunitstayid = E$patientunitstayid,
  uniquepid = as.character(E$uniquepid),
  C1 = as.integer(L$cluster_k2 == 1L),
  creatinine_frozen_coordinate = Xe[, "creatinine_max"],
  bun_frozen_coordinate = Xe[, "bun_max"],
  kdigo_0_24h = as.numeric(E$aki_stage_0_24h),
  sofa_total = as.numeric(E$sofa)
)
if (any(!is.finite(as.matrix(analysis_data[, .(
  C1, creatinine_frozen_coordinate, bun_frozen_coordinate, kdigo_0_24h, sofa_total
)])))) stop("Non-finite supported WP5 analysis variable")

fold_result <- make_patient_grouped_folds(analysis_data, CV_FOLDS, CV_SEED)
fold_id <- fold_result$fold
patient_fold_map <- fold_result$patient_map
models <- list(
  list(id = "creatinine", label = "Creatinine maximum", predictors = "creatinine_frozen_coordinate"),
  list(id = "bun", label = "Blood urea nitrogen maximum", predictors = "bun_frozen_coordinate"),
  list(id = "creatinine_bun", label = "Creatinine + BUN", predictors = c("creatinine_frozen_coordinate", "bun_frozen_coordinate")),
  list(id = "kdigo_0_24h", label = "KDIGO stage (0-24 h)", predictors = "kdigo_0_24h"),
  list(id = "sofa_total", label = "Total modified SOFA", predictors = "sofa_total")
)

cat("Running patient-grouped OOF identity models...\n")
model_results <- lapply(models, function(model) {
  cat("  ", model$id, "\n", sep = "")
  fit_oof_model(analysis_data, model, fold_id)
})
names(model_results) <- vapply(models, `[[`, character(1), "id")
model_summary <- rbindlist(lapply(model_results, `[[`, "summary"), fill = TRUE)
threshold_summary <- rbindlist(lapply(model_results, `[[`, "threshold_summary"), fill = TRUE)
predictions <- rbindlist(lapply(model_results, `[[`, "prediction"), fill = TRUE)
fold_audit <- rbindlist(lapply(model_results, `[[`, "fold_audit"), fill = TRUE)

pred_wide <- dcast(predictions, patientunitstayid + uniquepid + C1 ~ model_id, value.var = "oof_probability")
pred_wide[, row_order := match(patientunitstayid, analysis_data$patientunitstayid)]
setorder(pred_wide, row_order); pred_wide[, row_order := NULL]
if (!identical(pred_wide$patientunitstayid, analysis_data$patientunitstayid)) stop("Prediction row-order mismatch")

model_ids <- vapply(models, `[[`, character(1), "id")
patient_ids <- unique(pred_wide$uniquepid)
rows_by_patient <- split(seq_len(nrow(pred_wide)), pred_wide$uniquepid)
bootstrap_auc <- matrix(NA_real_, BOOT_B, length(model_ids), dimnames = list(NULL, model_ids))
bootstrap_kdigo <- matrix(
  NA_real_, BOOT_B, 4L,
  dimnames = list(NULL, c("p_c1_given_stage3", "p_stage3_given_c1", "stage3_vs_0_1_or", "ordinal_per_stage_or"))
)
set.seed(BOOT_SEED)
cat("Running paired patient-cluster AUC bootstrap (B=", BOOT_B, ")...\n", sep = "")
for (b in seq_len(BOOT_B)) {
  sampled <- sample(patient_ids, length(patient_ids), replace = TRUE)
  rows <- unlist(rows_by_patient[sampled], use.names = FALSE)
  yb <- pred_wide$C1[rows]
  if (length(unique(yb)) == 2L) {
    bootstrap_auc[b, ] <- vapply(model_ids, function(id) fast_auc(yb, pred_wide[[id]][rows]), numeric(1))
    kdigo_b <- analysis_data$kdigo_0_24h[rows]
    bootstrap_kdigo[b, "p_c1_given_stage3"] <- mean(yb[kdigo_b == 3L])
    bootstrap_kdigo[b, "p_stage3_given_c1"] <- mean(kdigo_b[yb == 1L] == 3L)
    stage3_keep_b <- kdigo_b %in% c(0L, 1L, 3L)
    stage3_tab_b <- table(
      factor(yb[stage3_keep_b], levels = 0:1),
      factor(as.integer(kdigo_b[stage3_keep_b] == 3L), levels = 0:1)
    )
    bootstrap_kdigo[b, "stage3_vs_0_1_or"] <-
      (stage3_tab_b[1L, 1L] * stage3_tab_b[2L, 2L]) /
      (stage3_tab_b[1L, 2L] * stage3_tab_b[2L, 1L])
    ordinal_fit_b <- suppressWarnings(glm.fit(
      x = cbind(`(Intercept)` = 1, kdigo_0_24h = kdigo_b),
      y = yb,
      family = binomial()
    ))
    if (isTRUE(ordinal_fit_b$converged) && length(ordinal_fit_b$coefficients) >= 2L &&
        is.finite(ordinal_fit_b$coefficients[2L])) {
      bootstrap_kdigo[b, "ordinal_per_stage_or"] <- exp(ordinal_fit_b$coefficients[2L])
    }
  }
  if (b %% 200L == 0L) cat(sprintf("  bootstrap %d/%d\n", b, BOOT_B))
}

auc_ci <- rbindlist(lapply(model_ids, function(id) {
  values <- bootstrap_auc[, id]
  data.table(
    model_id = id,
    auc_ci_low = unname(quantile(values, 0.025, type = 6, na.rm = TRUE)),
    auc_ci_high = unname(quantile(values, 0.975, type = 6, na.rm = TRUE)),
    auc_ci_method = "patient-clustered nonparametric bootstrap percentile CI",
    bootstrap_B = BOOT_B,
    bootstrap_success_count = sum(is.finite(values)),
    bootstrap_failure_count = sum(!is.finite(values))
  )
}))
model_summary <- merge(model_summary, auc_ci, by = "model_id", all.x = TRUE)

pair_definitions <- data.table(
  comparison = c(
    "Creatinine+BUN minus total modified SOFA",
    "Creatinine+BUN minus KDIGO 0-24 h",
    "Creatinine minus BUN"
  ),
  model_a = c("creatinine_bun", "creatinine_bun", "creatinine"),
  model_b = c("sofa_total", "kdigo_0_24h", "bun")
)
pairwise_auc <- pair_definitions[, {
  values <- bootstrap_auc[, model_a] - bootstrap_auc[, model_b]
  .(
    delta_auc = fast_auc(pred_wide$C1, pred_wide[[model_a]]) - fast_auc(pred_wide$C1, pred_wide[[model_b]]),
    ci_low = unname(quantile(values, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(values, 0.975, type = 6, na.rm = TRUE)),
    bootstrap_success_count = sum(is.finite(values)),
    bootstrap_failure_count = sum(!is.finite(values)),
    ci_method = "paired patient-clustered nonparametric bootstrap percentile CI"
  )
}, by = .(comparison, model_a, model_b)]

unsupported <- data.table(
  model_id = c("sofa_renal", "sofa_nonrenal"),
  model_label = c("SOFA renal component", "Non-renal SOFA"),
  predictors = NA_character_,
  analysis_status = "MISSING_NOT_RECONSTRUCTED",
  evaluation = NA_character_, statistical_unit = "ICU stay", clustered_unit = "uniquepid",
  n_stays = EXPECTED_STAYS, n_patients = EXPECTED_PATIENTS, c1_n = EXPECTED_C1,
  c1_prevalence = EXPECTED_C1 / EXPECTED_STAYS,
  harmonization_note = c(
    "eICU export retains modified total SOFA only; renal component was not exported and was not reconstructed from a substitute creatinine field",
    "eICU export retains modified total SOFA only; five non-renal components were not exported and total-minus-proxy subtraction was not used"
  )
)
for (column in setdiff(names(model_summary), names(unsupported))) unsupported[, (column) := NA]
setcolorder(unsupported, names(model_summary))
model_summary[, harmonization_note := fifelse(
  model_id %in% c("creatinine", "bun", "creatinine_bun"),
  "predictor is the exact frozen MIMIC transport coordinate used for external assignment",
  fifelse(model_id == "kdigo_0_24h", "combined eICU 0-24 h KDIGO stage; criterion components unavailable",
          "eICU modified total SOFA; component definitions differ from MIMIC SOFA")
)]
external_models <- rbindlist(list(model_summary, unsupported), fill = TRUE)

comparison_ids <- c("creatinine", "bun", "creatinine_bun", "kdigo_0_24h", "sofa_renal", "sofa_nonrenal", "sofa_total")
comparison <- data.table(model_id = comparison_ids)
comparison <- merge(
  comparison,
  internal[model_id %in% comparison_ids, .(
    model_id, internal_model_label = model_label,
    internal_auc = auc, internal_auc_ci_low = auc_ci_low, internal_auc_ci_high = auc_ci_high
  )], by = "model_id", all.x = TRUE
)
comparison <- merge(
  comparison,
  external_models[, .(
    model_id, external_model_label = model_label, external_analysis_status = analysis_status,
    external_auc = auc, external_auc_ci_low = auc_ci_low, external_auc_ci_high = auc_ci_high,
    harmonization_note
  )], by = "model_id", all.x = TRUE
)
comparison[, external_minus_internal_auc := external_auc - internal_auc]
comparison[, display_label := fcase(
  model_id == "creatinine", "Creatinine maximum",
  model_id == "bun", "Blood urea nitrogen maximum",
  model_id == "creatinine_bun", "Creatinine + BUN",
  model_id == "kdigo_0_24h", "KDIGO stage (0-24 h)",
  model_id == "sofa_renal", "SOFA renal component",
  model_id == "sofa_nonrenal", "Non-renal SOFA",
  model_id == "sofa_total", "Total SOFA / modified total SOFA",
  default = model_id
)]

kdigo_table <- analysis_data[, .(n = .N), by = .(
  cluster = fifelse(C1 == 1L, "C1", "C2"), kdigo_stage = as.integer(kdigo_0_24h)
)]
kdigo_table[, within_cluster_percent := 100 * n / sum(n), by = cluster]
kdigo_table[, within_stage_percent := 100 * n / sum(n), by = kdigo_stage]
kdigo_table[, overall_percent := 100 * n / sum(n)]
setorder(kdigo_table, kdigo_stage, cluster)

stage3_data <- analysis_data[kdigo_0_24h %in% c(0, 1, 3)]
stage3 <- factor(as.integer(stage3_data$kdigo_0_24h == 3), levels = 0:1)
stage3_table <- table(factor(stage3_data$C1, levels = 0:1), stage3)
stage3_cross_product_or <-
  (stage3_table[1L, 1L] * stage3_table[2L, 2L]) /
  (stage3_table[1L, 2L] * stage3_table[2L, 1L])
ordinal_fit <- glm(C1 ~ kdigo_0_24h, data = analysis_data, family = binomial())
kdigo_bootstrap_ci <- apply(bootstrap_kdigo, 2L, function(values) {
  unname(quantile(values, c(0.025, 0.975), type = 6, na.rm = TRUE))
})
kdigo_association <- data.table(
  metric = c(
    "P(C1 | KDIGO stage 3)", "P(KDIGO stage 3 | C1)",
    "Stage 3 vs stages 0-1 cross-product OR", "Ordinal per-stage OR"
  ),
  estimate = c(
    mean(analysis_data$C1[analysis_data$kdigo_0_24h == 3]),
    mean(analysis_data$kdigo_0_24h[analysis_data$C1 == 1L] == 3),
    stage3_cross_product_or, exp(coef(ordinal_fit)["kdigo_0_24h"])
  ),
  ci_low = c(
    kdigo_bootstrap_ci[1L, "p_c1_given_stage3"],
    kdigo_bootstrap_ci[1L, "p_stage3_given_c1"],
    kdigo_bootstrap_ci[1L, "stage3_vs_0_1_or"],
    kdigo_bootstrap_ci[1L, "ordinal_per_stage_or"]
  ),
  ci_high = c(
    kdigo_bootstrap_ci[2L, "p_c1_given_stage3"],
    kdigo_bootstrap_ci[2L, "p_stage3_given_c1"],
    kdigo_bootstrap_ci[2L, "stage3_vs_0_1_or"],
    kdigo_bootstrap_ci[2L, "ordinal_per_stage_or"]
  ),
  bootstrap_success_count = colSums(is.finite(bootstrap_kdigo)),
  bootstrap_failure_count = colSums(!is.finite(bootstrap_kdigo)),
  ci_method = "patient-clustered nonparametric bootstrap percentile CI",
  note = c(
    "column conditional; not interchangeable with P(stage 3 | C1)",
    "row conditional; not interchangeable with P(C1 | stage 3)",
    "binary stage-3 indicator versus stages 0-1", "KDIGO encoded 0-3"
  )
)

criterion_table <- data.table(
  aki_criterion = c("creatinine-only", "urine-output-only", "both", "neither"),
  C1_n = NA_integer_, C2_n = NA_integer_,
  analysis_status = "MISSING_NOT_RECONSTRUCTED",
  reason = paste(
    "The current eICU export retains only combined KDIGO stage.",
    "Criterion-level urine-output exposure duration and RRT criterion fields are not available; no proxy decomposition was created."
  )
)

END_HASHES_PRE_OUTPUT <- vapply(unlist(START_HASH_PATHS, use.names = FALSE), sha256, character(1))
names(END_HASHES_PRE_OUTPUT) <- names(START_HASH_PATHS)
numeric_model_cols <- c(
  "auc", "auc_ci_low", "auc_ci_high", "pr_auc_average_precision", "brier",
  "calibration_intercept", "calibration_slope", "sensitivity", "specificity",
  "ppv", "npv", "balanced_accuracy", "accuracy", "kappa", "ari",
  "predicted_c1_prevalence"
)
qc <- rbindlist(list(
  data.table(check = "input_hashes_unchanged", observed = paste0(sum(START_HASHES == END_HASHES_PRE_OUTPUT), "/", length(START_HASHES)), required = paste0(length(START_HASHES), "/", length(START_HASHES)), pass = identical(unname(START_HASHES), unname(END_HASHES_PRE_OUTPUT))),
  data.table(check = "current_eicu_denominator", observed = paste0(nrow(E), " stays; ", uniqueN(E$uniquepid), " patients; C1=", sum(analysis_data$C1)), required = paste0(EXPECTED_STAYS, " stays; ", EXPECTED_PATIENTS, " patients; C1=", EXPECTED_C1), pass = nrow(E) == EXPECTED_STAYS && uniqueN(E$uniquepid) == EXPECTED_PATIENTS && sum(analysis_data$C1) == EXPECTED_C1),
  data.table(check = "current_eicu_hospital_count", observed = as.character(uniqueN(E$hospitalid)), required = as.character(EXPECTED_HOSPITALS), pass = uniqueN(E$hospitalid) == EXPECTED_HOSPITALS),
  data.table(check = "obsolete_3417_label_excluded", observed = as.character(sum(analysis_data$C1)), required = "3462 and not 3417", pass = sum(analysis_data$C1) == 3462L),
  data.table(check = "authoritative_frozen_object_alignment", observed = paste0("run=", frozen_object$assignment_run_id, "; recipe/centroids/hashes/label MD5 verified"), required = "all authoritative frozen artifacts exactly aligned", pass = TRUE),
  data.table(check = "authoritative_eicu_input_fingerprint_alignment", observed = current_eicu_md5, required = as.character(object_eicu_fingerprint$md5), pass = identical(as.character(current_eicu_md5), as.character(object_eicu_fingerprint$md5)) && identical(as.character(current_eicu_md5), as.character(csv_eicu_fingerprint$md5))),
  data.table(check = "frozen_assignment_reproduction", observed = sprintf("%.6f", assignment_agreement), required = "1.000000", pass = assignment_agreement == 1),
  data.table(check = "no_eicu_reclustering", observed = unique(L$assignment_rule), required = "strict frozen assignment; no eICU re-clustering", pass = all(grepl("no eICU re-clustering", L$assignment_rule, fixed = TRUE))),
  data.table(check = "patient_grouped_fold_no_overlap", observed = as.character(sum(fold_audit$patient_overlap_n)), required = "0", pass = sum(fold_audit$patient_overlap_n) == 0L),
  data.table(check = "one_fold_per_patient", observed = paste0(nrow(patient_fold_map), " patients; max folds/patient=1"), required = paste0(EXPECTED_PATIENTS, " patients; one fold each"), pass = nrow(patient_fold_map) == EXPECTED_PATIENTS && !anyDuplicated(patient_fold_map$uniquepid)),
  data.table(check = "model_fold_and_youden_audit", observed = paste0("warnings=", sum(fold_audit$warning_count), "; Youden ok=", sum(fold_audit$youden_status == "ok"), "/", nrow(fold_audit), "; converged=", sum(fold_audit$fit_converged & fold_audit$coefficient_finite), "/", nrow(fold_audit)), required = paste0("warnings=0; Youden/convergence ok=", nrow(fold_audit), "/", nrow(fold_audit)), pass = sum(fold_audit$warning_count) == 0L && all(fold_audit$youden_status == "ok") && all(fold_audit$fit_converged) && all(fold_audit$coefficient_finite) && all(fold_audit$train_has_both_classes) && all(fold_audit$test_has_both_classes)),
  data.table(check = "complete_finite_oof_predictions", observed = paste0(sum(is.finite(predictions$oof_probability)), "/", EXPECTED_STAYS * length(models)), required = paste0(EXPECTED_STAYS * length(models), "/", EXPECTED_STAYS * length(models)), pass = nrow(predictions) == EXPECTED_STAYS * length(models) && all(is.finite(predictions$oof_probability))),
  data.table(check = "calibration_and_required_metrics", observed = paste0("hard calibration failures=", sum(model_summary$calibration_status == "failed_audit"), "; finite=", sum(is.finite(as.matrix(model_summary[, ..numeric_model_cols]))), "/", nrow(model_summary) * length(numeric_model_cols)), required = "0 hard failures; all required metrics finite", pass = all(model_summary$calibration_status %in% c("ok", "warning_finite")) && all(is.finite(as.matrix(model_summary[, ..numeric_model_cols])))),
  data.table(check = "paired_patient_bootstrap", observed = paste0(min(auc_ci$bootstrap_success_count), "-", max(auc_ci$bootstrap_success_count), "/", BOOT_B), required = paste0(BOOT_B, "/", BOOT_B, " for every model"), pass = all(auc_ci$bootstrap_success_count == BOOT_B)),
  data.table(check = "paired_patient_bootstrap_contrasts", observed = paste0(min(pairwise_auc$bootstrap_success_count), "-", max(pairwise_auc$bootstrap_success_count), "/", BOOT_B), required = paste0(BOOT_B, "/", BOOT_B, " for every contrast"), pass = all(pairwise_auc$bootstrap_success_count == BOOT_B)),
  data.table(check = "patient_bootstrap_kdigo_associations", observed = paste0(min(kdigo_association$bootstrap_success_count), "-", max(kdigo_association$bootstrap_success_count), "/", BOOT_B), required = paste0(BOOT_B, "/", BOOT_B, " for every KDIGO association"), pass = all(kdigo_association$bootstrap_success_count == BOOT_B)),
  data.table(check = "unsupported_fields_not_silently_substituted", observed = paste0(sum(external_models$analysis_status == "MISSING_NOT_RECONSTRUCTED"), " unsupported model rows; ", nrow(criterion_table), " criterion rows explicit"), required = "2 unsupported SOFA-component model rows; 4 explicit criterion rows", pass = sum(external_models$analysis_status == "MISSING_NOT_RECONSTRUCTED") == 2L && nrow(criterion_table) == 4L)
), fill = TRUE)

if (any(qc$pass == FALSE)) {
  failure_dir <- file.path(DIRS$logs, "failed_runs"); dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(qc, file.path(failure_dir, paste0("29d_WP5_QC_failure_", RUN_ID, ".csv")), bom = TRUE)
  stop("WP5 QC gate failed before fixed output creation: ", paste(qc[pass == FALSE, check], collapse = ", "))
}

write_csv_atomic(external_models, file.path(DIRS$tables, "Table_D1_external_identity_models.csv"))
write_csv_atomic(kdigo_table, file.path(DIRS$tables, "Table_D1_external_KDIGO_by_cluster.csv"))
write_csv_atomic(criterion_table, file.path(DIRS$tables, "Table_D1_external_AKI_criterion_by_cluster.csv"))
write_csv_atomic(comparison, file.path(DIRS$tables, "Table_D1_internal_external_identity_comparison.csv"))
write_csv_atomic(kdigo_association, file.path(DIRS$tables, "29d_Table_D1_external_KDIGO_association_summary.csv"))
write_csv_atomic(threshold_summary, file.path(DIRS$tables, "29d_Table_D1_external_identity_threshold_metrics.csv"))
write_csv_atomic(pairwise_auc, file.path(DIRS$tables, "29d_Table_D1_external_identity_pairwise_auc.csv"))
write_csv_atomic(predictions, file.path(DIRS$tables, "29d_Table_D1_external_identity_oof_predictions.csv"))
write_csv_atomic(fold_audit, file.path(DIRS$logs, "29d_external_identity_fold_audit.csv"))
write_csv_atomic(patient_fold_map, file.path(DIRS$logs, "29d_external_identity_patient_fold_map.csv"))
write_csv_atomic(qc, file.path(DIRS$logs, "29d_external_identity_QC.csv"))

plot_data <- comparison[external_analysis_status == "SUPPORTED_FORMAL"]
plot_long <- rbindlist(list(
  plot_data[, .(model_id, model_label = display_label, cohort = "MIMIC-IV", auc = internal_auc, low = internal_auc_ci_low, high = internal_auc_ci_high)],
  plot_data[, .(model_id, model_label = display_label, cohort = "eICU", auc = external_auc, low = external_auc_ci_low, high = external_auc_ci_high)]
))
label_order <- plot_data[order(external_auc), display_label]
plot_long[, model_label := factor(model_label, levels = unique(label_order))]
p <- ggplot(plot_long, aes(x = auc, y = model_label, shape = cohort, colour = cohort)) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.55, position = position_dodge(width = 0.42)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.42)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.45, colour = "grey70") +
  scale_colour_manual(values = c("MIMIC-IV" = "#2E2E2E", "eICU" = "#7A8791")) +
  scale_shape_manual(values = c("MIMIC-IV" = 16, "eICU" = 17)) +
  labs(
    x = "Out-of-fold AUC for C1 membership", y = NULL,
    title = "Internal and external reconstruction of the frozen C1 label",
    subtitle = "eICU labels were assigned using frozen MIMIC preprocessing and centroids; no external reclustering",
    colour = NULL, shape = NULL
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey35"),
    axis.title.x = element_text(size = 10), axis.text = element_text(size = 9),
    legend.position = "bottom", legend.text = element_text(size = 8.5),
    plot.margin = margin(8, 14, 8, 8)
  )
export_plot(p, "Figure_D1_internal_external_identity_AUC")

creatinine_kdigo_difference_n <- sum(abs(as.numeric(E$creatinine_max) - as.numeric(E$scr_peak_0_24h)) > 1e-12, na.rm = TRUE)
harmonization_lines <- c(
  "# WP5 external identity harmonization",
  "",
  paste0("- Run ID: `", RUN_ID, "`; mode: `", RUN_MODE, "`."),
  paste0("- eICU denominator: ", nrow(E), " ICU stays from ", uniqueN(E$uniquepid), " patients."),
  paste0("- Hospitals represented: ", uniqueN(E$hospitalid), "."),
  paste0("- Frozen C1: ", sum(analysis_data$C1), " (", sprintf("%.2f", 100 * mean(analysis_data$C1)), "%)."),
  paste0("- Assignment run ID: `", unique(L$assignment_run_id), "`; reproduced-label agreement: ", sprintf("%.6f", assignment_agreement), "."),
  "- No eICU reclustering was performed.",
  "- The recipe and centroids were verified against the authoritative frozen assignment RDS, including run ID, feature order, hashes, and label-file MD5.",
  "- Creatinine and BUN identity models use the exact frozen transport coordinates used to assign the eICU label.",
  paste0("- Raw BUN missing values handled by the frozen recipe: ", raw_missing[["bun_max"]], "; raw creatinine missing: ", raw_missing[["creatinine_max"]], "."),
  paste0("- `creatinine_max` and the KDIGO source `scr_peak_0_24h` differ in ", creatinine_kdigo_difference_n, " stay(s); they were not silently substituted for one another."),
  "- MIMIC identity models used the locked MICE-imputation-1 standardized matrix; eICU identity models use the frozen external transport recipe. This is a prespecified harmonization difference.",
  "- eICU total SOFA is a modified SOFA. Component fields were not exported, so renal SOFA and non-renal SOFA are reported as MISSING_NOT_RECONSTRUCTED.",
  "- Combined 0-24 h KDIGO stage is available. Criterion-level creatinine/urine-output/both/neither decomposition is unavailable and was not fabricated.",
  "- Cross-validation folds are grouped by `uniquepid`; external AUC confidence intervals use a patient-cluster bootstrap.",
  "- Patient-cluster bootstrap intervals condition on the fixed out-of-fold predictions; they do not include full model-refitting uncertainty.",
  "- These analyses reconstruct an algorithmic frozen label. They do not establish a disease mechanism or causal subtype."
)
write_lines_atomic(harmonization_lines, file.path(DIRS$logs, "29d_external_identity_harmonization.md"))

provenance <- data.table(
  source = names(START_HASH_PATHS), path = unlist(START_HASH_PATHS, use.names = FALSE),
  sha256 = unname(START_HASHES), read_only = TRUE, run_id = RUN_ID, run_mode = RUN_MODE
)
write_csv_atomic(provenance, file.path(DIRS$provenance, "29d_WP5_input_checksums.csv"))
metadata <- data.table(
  field = c(
    "run_id", "run_mode", "script_sha256", "assignment_run_id", "n_stays",
    "n_patients", "n_hospitals", "n_c1", "cv_folds", "cv_seed", "bootstrap_B", "bootstrap_seed",
    "external_reclustering", "assignment_agreement", "statistical_unit", "clustered_unit"
  ),
  value = c(
    RUN_ID, RUN_MODE, START_HASHES[["current_script_29d"]], unique(L$assignment_run_id),
    nrow(E), uniqueN(E$uniquepid), uniqueN(E$hospitalid), sum(analysis_data$C1), CV_FOLDS, CV_SEED,
    BOOT_B, BOOT_SEED, FALSE, assignment_agreement, "ICU stay", "uniquepid"
  )
)
write_csv_atomic(metadata, file.path(DIRS$provenance, "29d_WP5_run_metadata.csv"))

log_lines <- c(
  "# 29d WP5 external renal-identity audit",
  "",
  paste0("- Run ID: ", RUN_ID, "; formal mode."),
  paste0("- eICU: ", nrow(E), " stays; ", uniqueN(E$uniquepid), " patients; C1=", sum(analysis_data$C1), "."),
  paste0("- Frozen assignment reproduction: ", assignment_agreement, "."),
  paste0("- OOF: patient-grouped ", CV_FOLDS, " folds; seed ", CV_SEED, "."),
  paste0("- AUC CI: paired patient-cluster bootstrap B=", BOOT_B, "; seed ", BOOT_SEED, "."),
  "- Unsupported fields are explicitly marked MISSING_NOT_RECONSTRUCTED.",
  "",
  "## QC",
  paste(capture.output(print(qc)), collapse = "\n"),
  "",
  "## External models",
  paste(capture.output(print(external_models[, .(model_id, analysis_status, auc, auc_ci_low, auc_ci_high, pr_auc_average_precision, calibration_status)])), collapse = "\n"),
  "",
  "## Internal-external comparison",
  paste(capture.output(print(comparison[, .(model_id, external_analysis_status, internal_auc, external_auc, external_minus_internal_auc)])), collapse = "\n"),
  "",
  "## KDIGO association",
  paste(capture.output(print(kdigo_association)), collapse = "\n"),
  "",
  "## External paired AUC contrasts",
  paste(capture.output(print(pairwise_auc)), collapse = "\n")
)
write_lines_atomic(log_lines, file.path(DIRS$logs, "29d_external_identity_log.md"))
write_lines_atomic(capture.output(sessionInfo()), file.path(DIRS$logs, "29d_sessionInfo.txt"))

manifest_targets <- setdiff(FINAL_OUTPUTS, c(OUTPUT_MANIFEST, RUN_COMPLETED))
missing_targets <- manifest_targets[!file.exists(manifest_targets)]
if (length(missing_targets)) stop("Missing output before manifest: ", paste(missing_targets, collapse = ", "))
manifest <- data.table(
  file = basename(manifest_targets), path = manifest_targets,
  bytes = file.info(manifest_targets)$size,
  sha256 = vapply(manifest_targets, sha256, character(1)),
  run_id = RUN_ID, run_mode = RUN_MODE
)
write_csv_atomic(manifest, OUTPUT_MANIFEST)

completion_lines <- c(
  "RUN COMPLETED", paste0("run_id=", RUN_ID), paste0("run_mode=", RUN_MODE),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", START_HASHES[["current_script_29d"]]),
  paste0("output_manifest_sha256=", sha256(OUTPUT_MANIFEST))
)
completion_tmp <- paste0(RUN_COMPLETED, ".", RUN_ID, ".tmp")
writeLines(completion_lines, completion_tmp, useBytes = TRUE)
if (!file.rename(completion_tmp, RUN_COMPLETED)) stop("Could not atomically create completion marker")
unlink_status <- unlink(RUN_IN_PROGRESS)
if (unlink_status != 0L || file.exists(RUN_IN_PROGRESS)) {
  stop("Completion marker exists but in-progress marker could not be cleared")
}

cat("WP5 external identity audit completed.\n")
print(external_models[, .(model_id, analysis_status, auc, auc_ci_low, auc_ci_high)])
print(kdigo_association)
