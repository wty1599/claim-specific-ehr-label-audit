#!/usr/bin/env Rscript

## Domain 1 WP1 (first tranche): renal-identity audit
##
## Question: how well can the locked empirical C1 label be reproduced from
## creatinine, BUN, their combination, or 0-24 h KDIGO stage, compared with
## non-renal SOFA? This predicts the algorithmic C1 label, not a clinical
## outcome. It does not alter or replace the locked Domain 1 discreteness rule.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(pROC)
  library(mclust)
  library(ggplot2)
  library(digest)
})

## -------------------------------------------------------------------------
## Paths and immutable analysis settings
## -------------------------------------------------------------------------

project_root <- Sys.getenv(
  "MIMIC_IV_ROOT",
  unset = ""
)
if (!nzchar(project_root)) stop("Set MIMIC_IV_ROOT to the authorized local project root.")
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

analysis_root <- Sys.getenv(
  "D1_NEXT_PHASE_ROOT",
  unset = file.path(project_root, "domain1_next_phase_20260717")
)
analysis_root <- normalizePath(analysis_root, winslash = "/", mustWork = FALSE)

dirs <- list(
  scripts = file.path(analysis_root, "scripts"),
  tables = file.path(analysis_root, "tables"),
  figures = file.path(analysis_root, "figures"),
  logs = file.path(analysis_root, "logs"),
  checkpoints = file.path(analysis_root, "checkpoints"),
  manuscript = file.path(analysis_root, "manuscript_updates"),
  provenance = file.path(analysis_root, "provenance")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

input <- list(
  final_full = file.path(project_root, "data", "final_full.csv"),
  baseline_covars = file.path(project_root, "data", "baseline_covars.csv"),
  mice_primary = file.path(project_root, "output", "model", "mice_primary.rds"),
  labels_primary = file.path(project_root, "output", "model", "labels_primary_mice.rds"),
  standardized_matrix = file.path(project_root, "output", "model", "X_primary33_std_mice.rds"),
  label_builder = file.path(project_root, "r", "23c_save_mice_matrix_labels.R"),
  kdigo_builder = file.path(project_root, "sql", "build_baseline_covars_v2.sql")
)

missing_inputs <- names(input)[!file.exists(unlist(input))]
if (length(missing_inputs)) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

CV_FOLDS <- 10L
CV_SEED <- 2026071728L
PRIMARY_IMPUTATION <- 1L
CI_METHOD <- "DeLong CI on one out-of-fold prediction per patient"
CLASSIFICATION_THRESHOLD <- "fold-specific Youden threshold estimated in training fold"
EXPECTED_N <- 20049L
EXPECTED_C1_N <- 3992L
EXPECTED_C2_N <- 16057L

cat(
  "Domain 1 renal-identity audit\n",
  "  project root: ", project_root, "\n",
  "  output root: ", analysis_root, "\n",
  "  primary imputation: ", PRIMARY_IMPUTATION, "\n",
  "  CV: stratified ", CV_FOLDS, "-fold; seed ", CV_SEED, "\n",
  sep = ""
)

## -------------------------------------------------------------------------
## Utility functions
## -------------------------------------------------------------------------

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically write: ", path)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically write: ", path)
  invisible(path)
}

safe_divide <- function(num, den) ifelse(den > 0, num / den, NA_real_)

winsor <- function(v, p = c(0.01, 0.99)) {
  q <- quantile(v, p, na.rm = TRUE, names = FALSE)
  v[v < q[1L]] <- q[1L]
  v[v > q[2L]] <- q[2L]
  v
}

parse_sex <- function(g) as.integer(toupper(substr(as.character(g), 1L, 1L)) == "M")

script_path <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(hit)) return(NA_character_)
  normalizePath(sub("^--file=", "", hit[1L]), winslash = "/", mustWork = TRUE)
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

make_stratified_folds <- function(y, k, seed) {
  stopifnot(all(y %in% c(0L, 1L)))
  set.seed(seed)
  fold <- integer(length(y))
  for (level in c(0L, 1L)) {
    idx <- which(y == level)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

binary_agreement_metrics <- function(y, pred) {
  y <- as.integer(y)
  pred <- as.integer(pred)
  tp <- sum(y == 1L & pred == 1L)
  tn <- sum(y == 0L & pred == 0L)
  fp <- sum(y == 0L & pred == 1L)
  fn <- sum(y == 1L & pred == 0L)
  n <- length(y)
  accuracy <- (tp + tn) / n
  p_actual_pos <- (tp + fn) / n
  p_pred_pos <- (tp + fp) / n
  p_expected <- p_actual_pos * p_pred_pos + (1 - p_actual_pos) * (1 - p_pred_pos)
  kappa <- if (abs(1 - p_expected) < 1e-12) NA_real_ else (accuracy - p_expected) / (1 - p_expected)

  data.table(
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn,
    sensitivity = safe_divide(tp, tp + fn),
    specificity = safe_divide(tn, tn + fp),
    ppv = safe_divide(tp, tp + fp),
    npv = safe_divide(tn, tn + fn),
    accuracy = accuracy,
    kappa = kappa,
    ari = mclust::adjustedRandIndex(y, pred)
  )
}

cross_validated_logistic <- function(data, formula, model_id, model_label, fold_id) {
  n <- nrow(data)
  oof_probability <- rep(NA_real_, n)
  oof_class <- rep(NA_integer_, n)
  fold_rows <- vector("list", max(fold_id))

  for (fold in sort(unique(fold_id))) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    warning_text <- character()

    fit <- withCallingHandlers(
      glm(formula, data = data[train_idx], family = binomial()),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    p_train <- predict(fit, newdata = data[train_idx], type = "response")
    p_test <- predict(fit, newdata = data[test_idx], type = "response")
    if (any(!is.finite(p_train)) || any(!is.finite(p_test))) {
      stop("Non-finite probability for ", model_id, " fold ", fold)
    }

    roc_train <- pROC::roc(
      response = data$C1[train_idx], predictor = p_train,
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
    threshold_table <- pROC::coords(
      roc_train, x = "best", best.method = "youden",
      ret = "threshold", transpose = FALSE
    )
    threshold <- as.numeric(threshold_table$threshold[1L])
    if (!is.finite(threshold)) threshold <- 0.5

    oof_probability[test_idx] <- p_test
    oof_class[test_idx] <- as.integer(p_test >= threshold)
    fold_rows[[fold]] <- data.table(
      model_id = model_id,
      fold = fold,
      n_train = length(train_idx),
      n_test = length(test_idx),
      train_c1_prevalence = mean(data$C1[train_idx]),
      test_c1_prevalence = mean(data$C1[test_idx]),
      threshold = threshold,
      warning_count = length(warning_text),
      warning_text = paste(unique(warning_text), collapse = " | ")
    )
  }

  if (anyNA(oof_probability) || anyNA(oof_class)) {
    stop("Incomplete out-of-fold predictions for ", model_id)
  }

  roc_object <- pROC::roc(
    response = data$C1, predictor = oof_probability,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  auc_ci <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))

  bounded_probability <- pmin(pmax(oof_probability, 1e-6), 1 - 1e-6)
  logit_probability <- qlogis(bounded_probability)
  calibration_fit <- glm(C1 ~ logit_probability, data = data, family = binomial())
  agreement <- binary_agreement_metrics(data$C1, oof_class)
  fold_audit <- rbindlist(fold_rows)

  summary <- cbind(
    data.table(
      model_id = model_id,
      model_label = model_label,
      predictors = paste(all.vars(formula)[-1L], collapse = " + "),
      evaluation = paste0(CV_FOLDS, "-fold stratified out-of-fold"),
      n = n,
      c1_n = sum(data$C1 == 1L),
      c1_prevalence = mean(data$C1),
      auc = as.numeric(pROC::auc(roc_object)),
      auc_ci_low = auc_ci[1L],
      auc_ci_high = auc_ci[3L],
      auc_ci_method = CI_METHOD,
      brier = mean((oof_probability - data$C1)^2),
      calibration_intercept = unname(coef(calibration_fit)[1L]),
      calibration_slope = unname(coef(calibration_fit)[2L]),
      classification_threshold = CLASSIFICATION_THRESHOLD,
      threshold_median = median(fold_audit$threshold),
      threshold_min = min(fold_audit$threshold),
      threshold_max = max(fold_audit$threshold),
      glm_warning_n = sum(fold_audit$warning_count)
    ),
    agreement
  )

  list(
    summary = summary,
    prediction = data.table(
      stay_id = data$stay_id,
      model_id = model_id,
      fold = fold_id,
      C1 = data$C1,
      oof_probability = oof_probability,
      oof_class = oof_class
    ),
    roc = roc_object,
    fold_audit = fold_audit
  )
}

export_plot <- function(plot, stem, width = 7.2, height = 4.4) {
  pdf_path <- file.path(dirs$figures, paste0(stem, ".pdf"))
  png_path <- file.path(dirs$figures, paste0(stem, ".png"))
  tiff_path <- file.path(dirs$figures, paste0(stem, ".tiff"))

  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "sans")
  print(plot)
  grDevices::dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo")
  print(plot)
  grDevices::dev.off()

  grDevices::tiff(tiff_path, width = width, height = height, units = "in", res = 600, compression = "lzw")
  print(plot)
  grDevices::dev.off()

  c(pdf = pdf_path, png = png_path, tiff = tiff_path)
}

## -------------------------------------------------------------------------
## Recreate the exact row frame used to save MICE imputation #1 labels
## -------------------------------------------------------------------------

final_full <- fread(input$final_full, showProgress = FALSE)
baseline <- fread(input$baseline_covars, showProgress = FALSE)
labels <- as.data.table(readRDS(input$labels_primary))
X_locked <- as.data.table(readRDS(input$standardized_matrix))
imp <- readRDS(input$mice_primary)

if (!inherits(imp, "mids")) stop("mice_primary.rds is not a mids object.")
if (imp$m < PRIMARY_IMPUTATION) stop("Requested imputation is unavailable.")
if (!all(c("stay_id", "cluster_k2") %in% names(labels))) stop("Malformed locked label file.")
if (!"stay_id" %in% names(X_locked)) stop("Malformed locked standardized matrix.")
if (anyDuplicated(labels$stay_id) || anyDuplicated(X_locked$stay_id)) stop("Duplicate locked stay_id.")

features33 <- setdiff(names(X_locked), "stay_id")
if (length(features33) != 33L) stop("Expected 33 locked clustering features.")
if (!identical(features33, names(imp$data))) stop("MICE feature order differs from locked matrix.")

id_col <- "stay_id"
required_final <- c(
  id_col, "age", "gender", "sofa_score", "mortality_30d", "make30",
  features33, "respiration", "coagulation", "liver", "cardiovascular", "cns", "renal",
  "aki_stage_kdigo_7d"
)
missing_final <- setdiff(required_final, names(final_full))
if (length(missing_final)) stop("final_full.csv missing: ", paste(missing_final, collapse = ", "))
if (!"aki_stage_0_24h" %in% names(baseline)) stop("baseline_covars.csv lacks 0-24 h KDIGO stage.")

core <- final_full[, ..required_final]
core[, sex := parse_sex(gender)]
core <- merge(
  core,
  baseline[, .(stay_id, aki_stage_0_24h)],
  by = "stay_id", all.x = TRUE
)
core <- core[complete.cases(core[, .(
  age, sex, sofa_score, aki_stage_0_24h, mortality_30d, make30
)])]

RAW <- copy(core)
for (feature in features33) RAW[[feature]] <- winsor(as.numeric(RAW[[feature]]))
completed <- as.data.table(mice::complete(imp, PRIMARY_IMPUTATION))

if (nrow(RAW) != EXPECTED_N || nrow(completed) != EXPECTED_N) {
  stop("Unexpected empirical denominator: RAW=", nrow(RAW), ", completed=", nrow(completed))
}
if (!identical(as.integer(RAW$stay_id), as.integer(X_locked$stay_id))) {
  stop("Rebuilt MICE row frame does not match locked matrix stay_id order.")
}

X_rebuilt <- scale(as.matrix(completed[, ..features33]))
max_scaled_difference <- max(abs(X_rebuilt - as.matrix(X_locked[, ..features33])))
if (!is.finite(max_scaled_difference) || max_scaled_difference > 1e-8) {
  stop("Completed imputation #1 does not reproduce locked standardized matrix; max difference=", max_scaled_difference)
}

if (!identical(as.integer(labels$stay_id), as.integer(X_locked$stay_id))) {
  stop("Locked labels and matrix are not in identical stay_id order.")
}
label_counts <- table(labels$cluster_k2)
if (!identical(as.integer(label_counts[c("1", "2")]), c(EXPECTED_C1_N, EXPECTED_C2_N))) {
  stop("Locked C1/C2 counts differ from 3,992/16,057.")
}

analysis_data <- data.table(
  stay_id = RAW$stay_id,
  C1 = as.integer(labels$cluster_k2 == 1L),
  creatinine_max = as.numeric(completed$creatinine_max),
  bun_max = as.numeric(completed$bun_max),
  kdigo_stage_0_24h = as.integer(RAW$aki_stage_0_24h),
  kdigo_stage_7d = as.integer(RAW$aki_stage_kdigo_7d),
  sofa_total = as.numeric(RAW$sofa_score),
  sofa_renal = as.numeric(RAW$renal),
  sofa_nonrenal = as.numeric(RAW$respiration + RAW$coagulation + RAW$liver + RAW$cardiovascular + RAW$cns)
)

if (any(!is.finite(as.matrix(analysis_data[, -"stay_id"])))) {
  stop("Non-finite values in analysis dataset.")
}
if (any(analysis_data$sofa_total != analysis_data$sofa_nonrenal + analysis_data$sofa_renal)) {
  stop("SOFA total does not equal renal plus five non-renal components.")
}
if (!all(analysis_data$kdigo_stage_0_24h %in% 0:3)) stop("Invalid 0-24 h KDIGO stage.")

## -------------------------------------------------------------------------
## Same-fold out-of-fold prediction of the locked C1 label
## -------------------------------------------------------------------------

fold_id <- make_stratified_folds(analysis_data$C1, CV_FOLDS, CV_SEED)

model_specs <- list(
  list(id = "creatinine", label = "Creatinine maximum", formula = C1 ~ creatinine_max),
  list(id = "bun", label = "Blood urea nitrogen maximum", formula = C1 ~ bun_max),
  list(id = "creatinine_bun", label = "Creatinine + BUN", formula = C1 ~ creatinine_max + bun_max),
  list(id = "kdigo_0_24h", label = "KDIGO stage (0-24 h)", formula = C1 ~ kdigo_stage_0_24h),
  list(id = "sofa_nonrenal", label = "Non-renal SOFA", formula = C1 ~ sofa_nonrenal)
)

model_results <- lapply(model_specs, function(spec) {
  cross_validated_logistic(
    data = analysis_data,
    formula = spec$formula,
    model_id = spec$id,
    model_label = spec$label,
    fold_id = fold_id
  )
})
names(model_results) <- vapply(model_specs, `[[`, character(1), "id")

auc_table <- rbindlist(lapply(model_results, `[[`, "summary"), fill = TRUE)
prediction_table <- rbindlist(lapply(model_results, `[[`, "prediction"), fill = TRUE)
fold_audit <- rbindlist(lapply(model_results, `[[`, "fold_audit"), fill = TRUE)

pair_ids <- combn(names(model_results), 2L, simplify = FALSE)
contrast_table <- rbindlist(lapply(pair_ids, function(pair) {
  roc_a <- model_results[[pair[1L]]]$roc
  roc_b <- model_results[[pair[2L]]]$roc
  test <- pROC::roc.test(roc_a, roc_b, paired = TRUE, method = "delong", ci = TRUE)
  ci <- as.numeric(test$conf.int)
  data.table(
    model_1 = pair[1L],
    model_1_label = auc_table[model_id == pair[1L], model_label],
    model_2 = pair[2L],
    model_2_label = auc_table[model_id == pair[2L], model_label],
    auc_1 = as.numeric(pROC::auc(roc_a)),
    auc_2 = as.numeric(pROC::auc(roc_b)),
    delta_auc_model1_minus_model2 = as.numeric(pROC::auc(roc_a) - pROC::auc(roc_b)),
    delta_auc_ci_low = ci[1L],
    delta_auc_ci_high = ci[2L],
    p_value = as.numeric(test$p.value),
    method = "paired DeLong test on same out-of-fold predictions"
  )
}))

## -------------------------------------------------------------------------
## 0-24 h KDIGO stage x locked K2 label
## -------------------------------------------------------------------------

kdigo_long <- copy(analysis_data)[, cluster := ifelse(C1 == 1L, "C1", "C2")]
kdigo_counts <- kdigo_long[, .(n = .N), by = .(kdigo_stage_0_24h, cluster)]
kdigo_counts <- merge(
  CJ(kdigo_stage_0_24h = 0:3, cluster = c("C1", "C2"), sorted = TRUE),
  kdigo_counts,
  by = c("kdigo_stage_0_24h", "cluster"), all.x = TRUE
)
kdigo_counts[is.na(n), n := 0L]
kdigo_counts[, row_total := sum(n), by = kdigo_stage_0_24h]
kdigo_counts[, column_total := sum(n), by = cluster]
kdigo_counts[, `:=`(
  row_percent = 100 * n / row_total,
  column_percent = 100 * n / column_total,
  time_window = "ICU admission 0-24 h",
  source_field = "baseline_covars.csv: aki_stage_0_24h"
)]
setcolorder(kdigo_counts, c(
  "time_window", "source_field", "kdigo_stage_0_24h", "cluster",
  "n", "row_total", "row_percent", "column_total", "column_percent"
))

kdigo_matrix <- table(analysis_data$kdigo_stage_0_24h, analysis_data$C1)
chi <- suppressWarnings(chisq.test(kdigo_matrix, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) / (sum(kdigo_matrix) * min(nrow(kdigo_matrix) - 1L, ncol(kdigo_matrix) - 1L)))

kdigo_rule_metrics <- rbindlist(list(
  cbind(
    data.table(rule = "KDIGO stage 2-3", definition = "aki_stage_0_24h >= 2"),
    binary_agreement_metrics(analysis_data$C1, as.integer(analysis_data$kdigo_stage_0_24h >= 2L))
  ),
  cbind(
    data.table(rule = "KDIGO stage 3", definition = "aki_stage_0_24h == 3"),
    binary_agreement_metrics(analysis_data$C1, as.integer(analysis_data$kdigo_stage_0_24h == 3L))
  )
))

or_data <- analysis_data[kdigo_stage_0_24h %in% c(0L, 1L, 3L)]
or_table <- table(
  stage3 = factor(as.integer(or_data$kdigo_stage_0_24h == 3L), levels = 0:1),
  C1 = factor(or_data$C1, levels = 0:1)
)
or_test <- fisher.test(or_table)

kdigo_association <- data.table(
  time_window = "ICU admission 0-24 h",
  n = nrow(analysis_data),
  chi_square = as.numeric(chi$statistic),
  chi_square_df = as.integer(chi$parameter),
  chi_square_p = as.numeric(chi$p.value),
  cramers_v = cramers_v,
  stage3_vs_stage0_1_or_for_C1 = unname(or_test$estimate),
  stage3_vs_stage0_1_or_ci_low = or_test$conf.int[1L],
  stage3_vs_stage0_1_or_ci_high = or_test$conf.int[2L],
  stage3_vs_stage0_1_p = or_test$p.value,
  stage2_excluded_from_or = TRUE
)

## -------------------------------------------------------------------------
## Outputs
## -------------------------------------------------------------------------

write_csv_atomic(auc_table, file.path(dirs$tables, "28_Table_D1_renal_identity_auc.csv"))
write_csv_atomic(kdigo_counts, file.path(dirs$tables, "28_Table_D1_KDIGO_by_cluster.csv"))
write_csv_atomic(kdigo_rule_metrics, file.path(dirs$tables, "28_Table_D1_KDIGO_rule_metrics.csv"))
write_csv_atomic(kdigo_association, file.path(dirs$tables, "28_Table_D1_KDIGO_association_summary.csv"))
write_csv_atomic(contrast_table, file.path(dirs$tables, "28_Table_D1_renal_identity_pairwise_auc_contrasts.csv"))
write_csv_atomic(prediction_table, file.path(dirs$tables, "28_Table_D1_renal_identity_oof_predictions.csv"))
write_csv_atomic(fold_audit, file.path(dirs$logs, "28_cross_validation_fold_audit.csv"))

## Figure contract:
## Core conclusion: compare renal-function summaries with a non-renal severity
## control for reproducing the locked C1 label. This is a quantitative grid.
palette <- c(renal = "#3A5A8C", control = "#555555")
auc_plot_data <- copy(auc_table)
auc_plot_data[, family := ifelse(model_id == "sofa_nonrenal", "Non-renal control", "Renal/KDIGO")]
auc_plot_data[, model_label := factor(model_label, levels = rev(model_label[order(auc)]))]

p_auc <- ggplot(auc_plot_data, aes(x = auc, y = model_label, colour = family)) +
  geom_vline(xintercept = 0.5, colour = "grey70", linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = auc_ci_low, xmax = auc_ci_high),
    orientation = "y", width = 0, linewidth = 0.55
  ) +
  geom_point(size = 2.4) +
  geom_text(
    aes(label = sprintf("%.3f (%.3f-%.3f)", auc, auc_ci_low, auc_ci_high)),
    colour = "#2E2E2E", hjust = -0.05, size = 3.0, family = "sans"
  ) +
  scale_colour_manual(values = c("Renal/KDIGO" = palette[["renal"]], "Non-renal control" = palette[["control"]])) +
  scale_x_continuous(limits = c(0.48, 1.02), breaks = seq(0.5, 1.0, 0.1), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Reproduction of the locked C1 label",
    subtitle = "Ten-fold stratified out-of-fold logistic predictions; bars show 95% DeLong confidence intervals",
    x = "Area under the ROC curve",
    y = NULL,
    colour = NULL
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey35", margin = margin(b = 8)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    plot.margin = margin(8, 80, 8, 8)
  ) +
  coord_cartesian(clip = "off")

heat_data <- copy(kdigo_counts)
heat_data[, cluster := factor(cluster, levels = c("C1", "C2"))]
heat_data[, stage_label := factor(
  paste0("Stage ", kdigo_stage_0_24h),
  levels = rev(paste0("Stage ", 0:3))
)]

p_heat <- ggplot(heat_data, aes(x = cluster, y = stage_label, fill = row_percent)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%s\n%.1f%%", format(n, big.mark = ","), row_percent)), size = 3.2, family = "sans") +
  scale_fill_gradient(low = "#F4F6F8", high = "#3A5A8C", limits = c(0, 100), name = "Row %") +
  labs(
    title = "KDIGO stage and the locked K2 assignment",
    subtitle = "Full KDIGO stage during ICU admission hours 0-24; percentages are within KDIGO stage",
    x = NULL,
    y = "KDIGO stage (0-24 h)"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey35", margin = margin(b = 8)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5),
    plot.margin = margin(8, 8, 8, 8)
  )

auc_paths <- export_plot(p_auc, "28_Figure_D1_renal_identity_auc", width = 7.4, height = 4.5)
heat_paths <- export_plot(p_heat, "28_Figure_D1_KDIGO_cluster_heatmap", width = 6.4, height = 4.6)

script_file <- script_path()
manifest_paths <- c(unlist(input), script = script_file)
manifest_paths <- manifest_paths[!is.na(manifest_paths) & file.exists(manifest_paths)]
input_manifest <- rbindlist(lapply(names(manifest_paths), function(name) {
  path <- manifest_paths[[name]]
  data.table(
    object = name,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = sha256_file(path)
  )
}))
write_csv_atomic(input_manifest, file.path(dirs$provenance, "28_input_and_script_manifest.csv"))

top_auc <- auc_table[which.max(auc)]
nonrenal_auc <- auc_table[model_id == "sofa_nonrenal"]
creat_bun_auc <- auc_table[model_id == "creatinine_bun"]
kdigo_auc <- auc_table[model_id == "kdigo_0_24h"]

log_lines <- c(
  "# 28 Domain 1 renal-identity audit log",
  "",
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope and interpretation boundary",
  "",
  "This analysis estimates how well small renal-function summaries reproduce the locked empirical C1 label. It is not clinical-outcome prediction and does not replace the locked Domain 1 discreteness verdict, which remains inconclusive.",
  "",
  "## Locked inputs and alignment",
  "",
  paste0("- n = ", format(nrow(analysis_data), big.mark = ",")),
  paste0("- C1/C2 = ", sum(analysis_data$C1 == 1L), "/", sum(analysis_data$C1 == 0L)),
  paste0("- Primary MICE completed dataset = #", PRIMARY_IMPUTATION),
  paste0("- Maximum absolute difference when reproducing the locked standardized matrix = ", format(max_scaled_difference, scientific = TRUE)),
  "- Non-renal SOFA was calculated from the five authoritative source components: respiration, coagulation, liver, cardiovascular, and CNS.",
  "- Primary KDIGO predictor was the official full 0-24 h stage (creatinine, urine output, and RRT), aligned with the phenotype window.",
  "- Seven-day KDIGO was available but deliberately excluded from the primary prediction models because it extends beyond the phenotype window.",
  "",
  "## Evaluation",
  "",
  paste0("- Stratified ", CV_FOLDS, "-fold cross-validation; seed ", CV_SEED, "."),
  paste0("- AUC CI: ", CI_METHOD, "."),
  paste0("- Classification threshold: ", CLASSIFICATION_THRESHOLD, "."),
  "- All candidate models used the identical fold assignment and evaluation patients.",
  "",
  "## Key numerical results",
  "",
  paste0("- Highest OOF AUC: ", top_auc$model_label, " = ", sprintf("%.4f (%.4f-%.4f)", top_auc$auc, top_auc$auc_ci_low, top_auc$auc_ci_high), "."),
  paste0("- Creatinine + BUN OOF AUC = ", sprintf("%.4f (%.4f-%.4f)", creat_bun_auc$auc, creat_bun_auc$auc_ci_low, creat_bun_auc$auc_ci_high), "."),
  paste0("- KDIGO 0-24 h OOF AUC = ", sprintf("%.4f (%.4f-%.4f)", kdigo_auc$auc, kdigo_auc$auc_ci_low, kdigo_auc$auc_ci_high), "."),
  paste0("- Non-renal SOFA OOF AUC = ", sprintf("%.4f (%.4f-%.4f)", nonrenal_auc$auc, nonrenal_auc$auc_ci_low, nonrenal_auc$auc_ci_high), "."),
  paste0("- KDIGO x K2 Cramer's V = ", sprintf("%.4f", cramers_v), "."),
  paste0("- Stage 3 versus stages 0-1 OR for C1 = ", sprintf("%.3f (%.3f-%.3f)", unname(or_test$estimate), or_test$conf.int[1L], or_test$conf.int[2L]), "."),
  "",
  "## Missing or deferred WP1 components",
  "",
  "- SOFA renal and total SOFA prediction models: deferred by the user's current tranche, not missing from source data.",
  "- Creatinine-only versus urine-output-only AKI decomposition: component-specific 0-24 h stages were not present in the exported baseline CSV and were not reconstructed in this tranche.",
  "- Manuscript claim upgrade: deferred until the numerical results are reviewed against the prespecified evidence pattern.",
  "",
  "## Generated figures",
  "",
  paste0("- ", paste(auc_paths, collapse = "; ")),
  paste0("- ", paste(heat_paths, collapse = "; "))
)
write_lines_atomic(log_lines, file.path(dirs$logs, "28_domain1_renal_identity_audit_log.md"))

status_lines <- c(
  "# D1 next-phase status",
  "",
  "| Work package | Status | Note |",
  "|---|---|---|",
  "| WP1 renal identity | partially completed | Creatinine, BUN, creatinine+BUN, KDIGO 0-24 h, and non-renal SOFA control completed with same-fold OOF evaluation; remaining WP1 components deferred. |",
  "| WP2 loading contribution | pending | Not run in this tranche. |",
  "| WP3 S1 point-mass fidelity | pending | Not run. |",
  "| WP4 copula controls | blocked | Previous design-v1 Step 0 hard stop; requires a newly locked design before any further run. |",
  "| WP5 formal detection surface | pending | Not run. |",
  "| WP6 evaluation split sensitivity | pending | Not run. |",
  "| WP7 code governance | partially completed | New analysis is isolated, hashed, and uses project-relative configuration. |"
)
write_lines_atomic(status_lines, file.path(analysis_root, "D1_next_phase_status.md"))

discrepancy <- data.table(
  issue_id = c("D1-28-001", "D1-28-002"),
  item = c("AGENTS.md", "KDIGO time window"),
  status = c("missing", "resolved for current tranche"),
  detail = c(
    "No AGENTS.md was found under the project root or supplied context; the handoff and Domain1 context anchor governed this run.",
    "Both 0-24 h and 7-day KDIGO were available. The primary audit used 0-24 h KDIGO to match the phenotype window; 7-day KDIGO was not substituted."
  )
)
write_csv_atomic(discrepancy, file.path(analysis_root, "D1_discrepancy_log.csv"))

claim_map <- data.table(
  claim_id = c("D1-RI-01", "D1-RI-02", "D1-RI-03"),
  claim = c(
    "Creatinine and BUN reproduce the locked C1 label to the observed out-of-fold degree.",
    "The 0-24 h KDIGO stage is associated with locked K2 membership.",
    "Non-renal SOFA provides a same-fold comparator for C1 reproduction."
  ),
  status = "completed_pending_interpretation",
  source_file = c(
    "28_Table_D1_renal_identity_auc.csv",
    "28_Table_D1_KDIGO_by_cluster.csv; 28_Table_D1_KDIGO_association_summary.csv",
    "28_Table_D1_renal_identity_auc.csv"
  ),
  source_object = c("creatinine/bun/creatinine_bun rows", "KDIGO cross-table and association summary", "sofa_nonrenal row"),
  script = basename(script_file),
  run_mode = "formal empirical post-hoc audit",
  allowed_in_main_text = FALSE,
  allowed_in_supplement = TRUE,
  notes = "Do not upgrade the renal-identity hypothesis until numerical review is completed."
)
write_csv_atomic(claim_map, file.path(analysis_root, "D1_claim_evidence_map.csv"))

capture.output(sessionInfo(), file = file.path(dirs$logs, "28_sessionInfo.txt"))
write_lines_atomic("RUN COMPLETED", file.path(analysis_root, "28_run_completed.ok"))

cat("\nCompleted. Main outputs:\n")
cat("  ", file.path(dirs$tables, "28_Table_D1_renal_identity_auc.csv"), "\n", sep = "")
cat("  ", file.path(dirs$tables, "28_Table_D1_KDIGO_by_cluster.csv"), "\n", sep = "")
cat("  ", file.path(dirs$logs, "28_domain1_renal_identity_audit_log.md"), "\n", sep = "")
print(auc_table[, .(model_label, auc, auc_ci_low, auc_ci_high, brier, calibration_intercept, calibration_slope, sensitivity, specificity, ppv, npv, accuracy, kappa, ari)])
print(kdigo_association)
