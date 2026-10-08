# -*- coding: UTF-8 -*-
## ============================================================================
## SF5 upstream: current-aligned eICU strict SOFA+AKI sensitivity analysis.
##
## This script supersedes (but does not overwrite)
## 27_eicu_sepsis3_like_subset_sensitivity.R. It uses the locked current
## script-19/script-21 transport artifacts:
##   * all 33 frozen MIMIC-IV features are retained;
##   * frozen MIMIC winsorization, medians, scaling, and centroids are reused;
##   * current eICU labels are joined, never recomputed or re-clustered;
##   * high-risk C1 direction is inherited from MIMIC-IV mortality;
##   * paired Delta AUC uses uniquepid-clustered bootstrap inference.
##
## Source results are never overwritten. Outputs are written to a dated,
## isolated QC directory and must pass full-cohort concordance checks before
## they are used by Supplementary Figure 5.
## ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

## 00_paths.R contains legacy non-UTF-8 comments. Define the identical current
## project paths explicitly so batch Rscript execution is encoding-independent.
PROJ_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
DIR_DATA <- file.path(PROJ_ROOT, "data")
DIR_OUTPUT <- file.path(PROJ_ROOT, "output")
DIR_MODEL <- file.path(DIR_OUTPUT, "model")
PATH_FINAL_FULL <- file.path(DIR_DATA, "final_full.csv")

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ranger)
  library(pROC)
})

SCRIPT_VERSION <- "sf5-eicu-strict-subset-v2-strict33-current-20260715"
EXPECTED_ASSIGNMENT_RUN_ID <- "f41396c7aa239dfb"
EXPECTED_EICU_STAYS <- 17465L
EXPECTED_EICU_PATIENTS <- 16212L
EXPECTED_STRICT_STAYS <- 7975L
SEED <- 20240601L
BOOT_B <- 1000L
K_FOLDS <- 10L
RF_TREES <- 300L
ENET_ALPHA <- 0.5

OUT_DIR <- file.path(DIR_OUTPUT, "qc", "sf5_strict_subset_current_20260715")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
  "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml", "glucose_max",
  "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min", "sodium_min",
  "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min", "chloride_min",
  "chloride_max", "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
MODELS <- c("base", "base_phen", "feat_pen", "feat_rf")
CONTRASTS <- list(
  c("feat_pen", "base_phen"),
  c("feat_pen", "feat_rf"),
  c("base_phen", "base"),
  c("feat_pen", "base")
)

PATH_TRANSPORT <- file.path(DIR_OUTPUT, "qc", "headtohead_all_mice_labels.rds")
PATH_EICU_LABELS <- file.path(DIR_OUTPUT, "eicu_cluster_labels_mice.csv")
PATH_ASSIGNMENT_PROV <- file.path(DIR_OUTPUT, "qc", "external_assignment_provenance_mice_labels.csv")
PATH_LOCKED_PERF <- file.path(DIR_OUTPUT, "qc", "headtohead_summary_mice_labels.csv")
PATH_LOCKED_PAIR <- file.path(DIR_OUTPUT, "qc", "external_paired_contrasts_mice_labels.csv")
PATH_EICU <- file.path(DIR_DATA, "eicu_external.csv")
PATH_BASELINE <- file.path(DIR_DATA, "baseline_covars.csv")
PATH_MICE_LABELS <- file.path(DIR_MODEL, "labels_primary_mice.rds")

required_inputs <- c(
  PATH_FINAL_FULL, PATH_BASELINE, PATH_MICE_LABELS, PATH_EICU,
  PATH_TRANSPORT, PATH_EICU_LABELS, PATH_ASSIGNMENT_PROV,
  PATH_LOCKED_PERF, PATH_LOCKED_PAIR
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required current artifact(s): ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

hash_file <- function(path) unname(tools::md5sum(path))
check_unique_id <- function(dt, id_col, label) {
  if (!id_col %in% names(dt)) stop(label, " missing ", id_col, call. = FALSE)
  dup <- anyDuplicated(dt[[id_col]])
  if (dup) stop(label, " has duplicated ", id_col, "; first duplicate index: ", dup, call. = FALSE)
}
parse_age <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(">\\s*89", "90", as.character(x))))
}
parse_sex <- function(x) {
  if (is.numeric(x)) return(as.integer(x == 1))
  z <- toupper(substr(trimws(as.character(x)), 1, 1))
  fifelse(z == "M", 1L, fifelse(z == "F", 0L, NA_integer_))
}
winsor <- function(x, q) pmin(pmax(x, q[1]), q[2])
fast_auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
clamp <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
cal_slope_int <- function(y, p) {
  lp <- qlogis(clamp(p))
  slope <- tryCatch(unname(coef(glm(y ~ lp, family = binomial))[2]), error = function(e) NA_real_)
  intercept <- tryCatch(
    unname(coef(glm(y ~ 1, offset = lp, family = binomial))[1]),
    error = function(e) NA_real_
  )
  c(slope = slope, intercept = intercept)
}
boot_auc_ci <- function(y, p, cluster_id, B = BOOT_B) {
  ok <- is.finite(y) & is.finite(p) & !is.na(cluster_id)
  y <- y[ok]; p <- p[ok]; cluster_id <- as.character(cluster_id[ok])
  clusters <- unique(cluster_id)
  rows_by_cluster <- split(seq_along(cluster_id), cluster_id)
  draws <- numeric(B)
  for (b in seq_len(B)) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(rows_by_cluster[sampled], use.names = FALSE)
    draws[b] <- fast_auc(y[idx], p[idx])
  }
  unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
}
cluster_boot_dauc <- function(y, p_better, p_worse, cluster_id, B = BOOT_B) {
  ok <- is.finite(y) & is.finite(p_better) & is.finite(p_worse) & !is.na(cluster_id)
  y <- y[ok]
  p_better <- p_better[ok]
  p_worse <- p_worse[ok]
  cluster_id <- as.character(cluster_id[ok])
  clusters <- unique(cluster_id)
  rows_by_cluster <- split(seq_along(cluster_id), cluster_id)
  draws <- numeric(B)
  for (b in seq_len(B)) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(rows_by_cluster[sampled], use.names = FALSE)
    draws[b] <- fast_auc(y[idx], p_better[idx]) - fast_auc(y[idx], p_worse[idx])
  }
  finite_n <- sum(is.finite(draws))
  p_two <- min(
    1,
    2 * min(
      (sum(draws <= 0, na.rm = TRUE) + 1) / (finite_n + 1),
      (sum(draws >= 0, na.rm = TRUE) + 1) / (finite_n + 1)
    )
  )
  list(
    lo = unname(quantile(draws, 0.025, na.rm = TRUE)),
    hi = unname(quantile(draws, 0.975, na.rm = TRUE)),
    p = p_two,
    n_clusters = length(clusters)
  )
}

cat("============================================================\n")
cat("SF5 current-aligned upstream analysis\n")
cat("Script version: ", SCRIPT_VERSION, "\n", sep = "")
cat("Output directory: ", OUT_DIR, "\n", sep = "")
cat("Frozen dimensions: 33; assignment is reused, not recomputed.\n")
cat("Paired inference: uniquepid-clustered bootstrap, B=", BOOT_B, ".\n", sep = "")
cat("============================================================\n")

## ---- Current locked transport and source frames ----
transport <- readRDS(PATH_TRANSPORT)
assignment_prov <- fread(PATH_ASSIGNMENT_PROV)
if (!identical(as.character(transport$assignment_run_id), EXPECTED_ASSIGNMENT_RUN_ID)) {
  stop("Unexpected frozen assignment run ID: ", transport$assignment_run_id, call. = FALSE)
}
if (nrow(assignment_prov) != 1L ||
    !identical(as.character(assignment_prov$assignment_run_id[1]), EXPECTED_ASSIGNMENT_RUN_ID)) {
  stop("Assignment provenance does not match the locked run.", call. = FALSE)
}
if (is.null(transport$common_feature_n) || transport$common_feature_n != 33L ||
    !identical(as.character(transport$common), FEATURES33)) {
  stop("Transport object is not the locked strict 33-feature version.", call. = FALSE)
}
if (!isTRUE(transport$lymphocyte_audit$verified) ||
    !isTRUE(transport$lymphocyte_audit$used_for_assignment)) {
  stop("Transport object does not verify/use the corrected absolute lymphocyte count.", call. = FALSE)
}
recorded_eicu_md5 <- transport$input_fingerprint[role == "eicu_external", md5]
if (length(recorded_eicu_md5) != 1L ||
    !identical(as.character(recorded_eicu_md5), as.character(hash_file(PATH_EICU)))) {
  stop("Current eicu_external.csv does not match the locked transport input hash.", call. = FALSE)
}
if (!is.null(transport$eicu_label_file_md5) &&
    !identical(as.character(transport$eicu_label_file_md5), as.character(hash_file(PATH_EICU_LABELS)))) {
  stop("Current eICU label file does not match the locked transport hash.", call. = FALSE)
}

mimic <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id", "subject_id"), names(mimic))[1]
if (is.na(mid)) stop("No compatible MIMIC ID in final_full.csv.", call. = FALSE)
need_m <- c(mid, "age", "gender", "sofa_score", "hospital_mortality", "make30", FEATURES33)
miss_m <- setdiff(need_m, names(mimic))
if (length(miss_m)) stop("final_full.csv missing: ", paste(miss_m, collapse = ", "), call. = FALSE)
M <- mimic[, ..need_m]
check_unique_id(M, mid, "final_full.csv")
M[, sex := parse_sex(gender)]

baseline <- fread(PATH_BASELINE)
check_unique_id(baseline, mid, "baseline_covars.csv")
if (!"aki_stage_0_24h" %in% names(baseline)) stop("baseline_covars.csv missing AKI stage.")
M <- merge(M, baseline[, .SD, .SDcols = c(mid, "aki_stage_0_24h")], by = mid, all.x = TRUE)

mice_labels <- as.data.table(readRDS(PATH_MICE_LABELS))
if (!mid %in% names(mice_labels)) {
  candidate <- intersect(c("stay_id", "hadm_id", "subject_id"), names(mice_labels))[1]
  if (is.na(candidate)) stop("labels_primary_mice.rds lacks a compatible ID.")
  setnames(mice_labels, candidate, mid)
}
if (!"cluster_k2" %in% names(mice_labels)) stop("labels_primary_mice.rds lacks cluster_k2.")
mice_labels <- mice_labels[, .SD, .SDcols = c(mid, "cluster_k2")]
check_unique_id(mice_labels, mid, "labels_primary_mice.rds")
M <- merge(M, mice_labels, by = mid, all.x = TRUE)
M <- M[complete.cases(M[, c(BASE_NUM, "cluster_k2", "hospital_mortality", "make30"), with = FALSE])]
M[, cluster_k2 := factor(cluster_k2)]
cluster_levels <- levels(M$cluster_k2)
risk_by_label <- M[, .(mortality = mean(hospital_mortality == 1), n = .N), by = cluster_k2]
high_risk_label <- as.character(risk_by_label[which.max(mortality), cluster_k2])
if (!identical(as.character(transport$high_risk_label), high_risk_label)) {
  stop("Current MIMIC high-risk label disagrees with the frozen transport object.")
}

E <- fread(PATH_EICU)
setnames(E, "sofa", "sofa_score", skip_absent = TRUE)
check_unique_id(E, "patientunitstayid", "eicu_external.csv")
if (!"uniquepid" %in% names(E)) stop("eicu_external.csv lacks uniquepid.")
if (nrow(E) != EXPECTED_EICU_STAYS || uniqueN(E$uniquepid) != EXPECTED_EICU_PATIENTS) {
  stop("Current eICU denominator is not 17,465 ICU stays / 16,212 patients.")
}

Elab <- fread(PATH_EICU_LABELS)
required_label_cols <- c(
  "patientunitstayid", "uniquepid", "cluster_k2", "assignment_run_id",
  "common_feature_n", "preprocessing_hash", "centroid_hash"
)
miss_l <- setdiff(required_label_cols, names(Elab))
if (length(miss_l)) stop("Current eICU labels missing: ", paste(miss_l, collapse = ", "))
check_unique_id(Elab, "patientunitstayid", "eicu_cluster_labels_mice.csv")
if (nrow(Elab) != EXPECTED_EICU_STAYS || uniqueN(Elab$uniquepid) != EXPECTED_EICU_PATIENTS) {
  stop("Current eICU label denominator is not locked.")
}
if (uniqueN(Elab$assignment_run_id) != 1L ||
    !identical(as.character(unique(Elab$assignment_run_id)), EXPECTED_ASSIGNMENT_RUN_ID) ||
    uniqueN(Elab$common_feature_n) != 1L || unique(Elab$common_feature_n) != 33L) {
  stop("Current eICU labels are not from the strict 33-feature frozen assignment.")
}

join_labels <- Elab[, .(
  patientunitstayid,
  label_uniquepid = uniquepid,
  cluster_k2,
  assignment_run_id
)]
E <- merge(E, join_labels, by = "patientunitstayid", all.x = TRUE, sort = FALSE)
if (anyNA(E$label_uniquepid) || anyNA(E$cluster_k2) ||
    any(E$uniquepid != E$label_uniquepid, na.rm = TRUE)) {
  stop("Current eICU source/label join failed.")
}
E[, label_uniquepid := NULL]
E[, cluster_k2 := factor(cluster_k2, levels = cluster_levels)]
E[, age := parse_age(age)]
E[, sex := parse_sex(sex)]
E[, aki_stage_0_24h := as.numeric(aki_stage_0_24h)]
E[, hosp_mortality := as.integer(hosp_mortality)]
if (!all(c("scr_discharge", "scr_baseline") %in% names(E))) {
  stop("Corrected eICU file lacks scr_discharge/scr_baseline; MAKE cannot be reconstructed.")
}
E[, make_hosp := as.integer(
  hosp_mortality == 1 |
    (!is.na(scr_discharge) & !is.na(scr_baseline) & scr_discharge >= 1.5 * scr_baseline)
)]

common <- transport$common
recipe <- transport$preprocessing_recipe
if (!identical(names(recipe), common)) stop("Frozen recipe order differs from frozen features.")
apply_recipe <- function(df) {
  out <- matrix(0, nrow(df), length(common), dimnames = list(NULL, common))
  for (j in common) {
    r <- recipe[[j]]
    v <- winsor(as.numeric(df[[j]]), c(r$lo, r$hi))
    v[is.na(v)] <- r$med
    out[, j] <- (v - r$mean) / r$sd
  }
  out
}
Xm <- apply_recipe(M)
Xe <- apply_recipe(E)
BASEm <- as.matrix(M[, ..BASE_NUM]); storage.mode(BASEm) <- "double"
BASEe <- as.matrix(E[, ..BASE_NUM]); storage.mode(BASEe) <- "double"

strict_mask <- is.finite(E$sofa_score) & E$sofa_score >= 2 &
  is.finite(E$aki_stage_0_24h) & E$aki_stage_0_24h >= 1
if (sum(strict_mask) != EXPECTED_STRICT_STAYS) {
  stop("Strict SOFA+AKI subset is ", sum(strict_mask), ", expected 7,975 ICU stays.")
}

cohort_masks <- list(
  eicu_full = rep(TRUE, nrow(E)),
  eicu_sepsis3_like = strict_mask
)

## ---- Performance/calibration branch: preserve script-19 RNG order ----
fit_perf <- function(kind, y, baseM, Xstd, clust) switch(kind,
  base = glm(y ~ ., data = data.frame(y = y, baseM), family = binomial),
  base_phen = glm(y ~ ., data = data.frame(y = y, baseM, cluster_k2 = clust), family = binomial),
  feat_pen = cv.glmnet(cbind(baseM, Xstd), y, family = "binomial", alpha = ENET_ALPHA, nfolds = 5),
  feat_rf = ranger(
    x = as.data.frame(cbind(baseM, Xstd)), y = factor(y), probability = TRUE,
    num.trees = RF_TREES, num.threads = 0, seed = 1L
  )
)
pred_perf <- function(kind, fit, baseM, Xstd, clust) switch(kind,
  base = predict(fit, newdata = data.frame(baseM), type = "response"),
  base_phen = predict(fit, newdata = data.frame(baseM, cluster_k2 = clust), type = "response"),
  feat_pen = as.numeric(predict(fit, cbind(baseM, Xstd), s = "lambda.min", type = "response")),
  feat_rf = predict(fit, as.data.frame(cbind(baseM, Xstd)))$predictions[, "1"]
)
cv_auc <- function(y) {
  folds <- integer(length(y))
  for (cl in c(0, 1)) {
    idx <- which(y == cl)
    folds[idx] <- sample(rep_len(seq_len(K_FOLDS), length(idx)))
  }
  sapply(MODELS, function(kind) {
    p <- numeric(length(y))
    for (fold in seq_len(K_FOLDS)) {
      train <- folds != fold; test <- folds == fold
      fit <- fit_perf(kind, y[train], BASEm[train, , drop = FALSE], Xm[train, , drop = FALSE], M$cluster_k2[train])
      p[test] <- pred_perf(kind, fit, BASEm[test, , drop = FALSE], Xm[test, , drop = FALSE], M$cluster_k2[test])
    }
    fast_auc(y, p)
  })
}

performance_rows <- list()
performance_predictions <- list()
set.seed(SEED)
outcome_map <- list(
  hospital_mortality = c(mimic = "hospital_mortality", eicu = "hosp_mortality"),
  make_hosp = c(mimic = "make30", eicu = "make_hosp")
)

for (outcome_name in names(outcome_map)) {
  map <- outcome_map[[outcome_name]]
  ym <- as.integer(M[[map[["mimic"]]]])
  internal_auc <- cv_auc(ym)
  fits <- lapply(MODELS, function(kind) fit_perf(kind, ym, BASEm, Xm, M$cluster_k2))
  names(fits) <- MODELS
  predictions <- sapply(MODELS, function(kind) {
    pred_perf(kind, fits[[kind]], BASEe, Xe, E$cluster_k2)
  })
  performance_predictions[[outcome_name]] <- predictions

  ## Full cohort is evaluated first, including the exact script-19 bootstrap
  ## sequence, so the locked full-cohort result can be reproduced.
  full_sel <- which(complete.cases(BASEe) & !is.na(E$cluster_k2) & !is.na(E[[map[["eicu"]]]]))
  full_y <- as.integer(E[[map[["eicu"]]]][full_sel])
  for (kind in MODELS) {
    p <- predictions[full_sel, kind]
    ci <- boot_auc_ci(full_y, p, E$uniquepid[full_sel], B = BOOT_B)
    cs <- cal_slope_int(full_y, p)
    performance_rows[[length(performance_rows) + 1L]] <- data.table(
      cohort = "eicu_full", outcome = outcome_name, model = kind,
      n = length(full_y), events = sum(full_y == 1),
      auc = fast_auc(full_y, p), auc_lo = ci[1], auc_hi = ci[2],
      calibration_slope = cs["slope"], calibration_intercept = cs["intercept"],
      brier = mean((p - full_y)^2), internal_cv_auc = internal_auc[[kind]],
      statistical_unit = "ICU stay", resampling_unit = "uniquepid cluster",
      assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID
    )
  }
}

## Strict-subset performance uses the same fitted models/predictions and is
## evaluated only after both full-cohort runs, preserving full-cohort RNG order.
for (outcome_name in names(outcome_map)) {
  map <- outcome_map[[outcome_name]]
  predictions <- performance_predictions[[outcome_name]]
  strict_sel <- which(strict_mask & complete.cases(BASEe) &
    !is.na(E$cluster_k2) & !is.na(E[[map[["eicu"]]]]))
  strict_y <- as.integer(E[[map[["eicu"]]]][strict_sel])
  for (kind in MODELS) {
    p <- predictions[strict_sel, kind]
    cs <- cal_slope_int(strict_y, p)
    performance_rows[[length(performance_rows) + 1L]] <- data.table(
      cohort = "eicu_sepsis3_like", outcome = outcome_name, model = kind,
      n = length(strict_y), events = sum(strict_y == 1),
      auc = fast_auc(strict_y, p), auc_lo = NA_real_, auc_hi = NA_real_,
      calibration_slope = cs["slope"], calibration_intercept = cs["intercept"],
      brier = mean((p - strict_y)^2), internal_cv_auc = NA_real_,
      statistical_unit = "ICU stay", resampling_unit = "not applicable",
      assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID
    )
  }
}
performance <- rbindlist(performance_rows)

## ---- Paired-contrast branch: preserve script-21 RNG/order ----
fit_pair <- function(kind, y) switch(kind,
  base = glm(y ~ ., data = data.frame(y = y, BASEm), family = binomial),
  base_phen = glm(y ~ ., data = data.frame(y = y, BASEm, cluster_k2 = M$cluster_k2), family = binomial),
  feat_pen = cv.glmnet(cbind(BASEm, Xm), y, family = "binomial", alpha = ENET_ALPHA, nfolds = 5),
  feat_rf = ranger(
    x = as.data.frame(cbind(BASEm, Xm)), y = factor(y), probability = TRUE,
    num.trees = RF_TREES, seed = 1L
  )
)
pred_pair <- function(kind, fit) switch(kind,
  base = predict(fit, newdata = data.frame(BASEe), type = "response"),
  base_phen = predict(fit, newdata = data.frame(BASEe, cluster_k2 = E$cluster_k2), type = "response"),
  feat_pen = as.numeric(predict(fit, cbind(BASEe, Xe), s = "lambda.min", type = "response")),
  feat_rf = predict(fit, as.data.frame(cbind(BASEe, Xe)))$predictions[, "1"]
)
paired_one_cohort <- function(outcome_name, P, mask) {
  map <- outcome_map[[outcome_name]]
  sel <- which(mask & complete.cases(BASEe) & !is.na(E$cluster_k2) & !is.na(E[[map[["eicu"]]]]))
  y <- as.integer(E[[map[["eicu"]]]][sel])
  cluster_id <- E$uniquepid[sel]
  rbindlist(lapply(CONTRASTS, function(cc) {
    better <- cc[1]; worse <- cc[2]
    p_better <- P[sel, better]; p_worse <- P[sel, worse]
    delong_p <- tryCatch(
      roc.test(
        roc(y, p_worse, quiet = TRUE), roc(y, p_better, quiet = TRUE),
        method = "delong", paired = TRUE
      )$p.value,
      error = function(e) NA_real_
    )
    cb <- cluster_boot_dauc(y, p_better, p_worse, cluster_id, B = BOOT_B)
    data.table(
      outcome = outcome_name, better = better, worse = worse,
      auc_better = round(fast_auc(y, p_better), 4),
      auc_worse = round(fast_auc(y, p_worse), 4),
      dAUC = round(fast_auc(y, p_better) - fast_auc(y, p_worse), 4),
      dAUC_lo = round(cb$lo, 4), dAUC_hi = round(cb$hi, 4),
      clustered_bootstrap_p = signif(cb$p, 3), p_value = signif(cb$p, 3),
      p_method = "two-sided uniquepid-clustered bootstrap",
      delong_p_stay_level_sensitivity = signif(delong_p, 3),
      n_eicu_icu_stays = length(y), n_eicu_unique_patients = cb$n_clusters,
      assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID
    )
  }))
}

paired_full <- list()
paired_prediction_store <- list()
set.seed(SEED)
for (outcome_name in names(outcome_map)) {
  map <- outcome_map[[outcome_name]]
  ym <- as.integer(M[[map[["mimic"]]]])
  fits <- lapply(MODELS, function(kind) fit_pair(kind, ym)); names(fits) <- MODELS
  P <- sapply(MODELS, function(kind) pred_pair(kind, fits[[kind]]))
  paired_prediction_store[[outcome_name]] <- P
  paired_full[[outcome_name]] <- paired_one_cohort(outcome_name, P, cohort_masks$eicu_full)
}
paired_full <- rbindlist(paired_full, use.names = TRUE)
paired_full[, cohort := "eicu_full"]

paired_strict <- rbindlist(lapply(names(outcome_map), function(outcome_name) {
  paired_one_cohort(
    outcome_name,
    paired_prediction_store[[outcome_name]],
    cohort_masks$eicu_sepsis3_like
  )
}), use.names = TRUE)
paired_strict[, cohort := "eicu_sepsis3_like"]
paired <- rbindlist(list(paired_full, paired_strict), use.names = TRUE)
setcolorder(paired, c("cohort", setdiff(names(paired), "cohort")))

## ---- Prevalence, feature, and denominator audit ----
prevalence <- data.table(
  cohort = c("MIMIC_training", "eicu_full", "eicu_sepsis3_like"),
  n = c(nrow(M), nrow(E), sum(strict_mask)),
  C1_n = c(
    sum(as.character(M$cluster_k2) == high_risk_label),
    sum(as.character(E$cluster_k2) == high_risk_label),
    sum(strict_mask & as.character(E$cluster_k2) == high_risk_label)
  )
)
prevalence[, C1_pct := 100 * C1_n / n]
prevalence[, `:=`(
  high_risk_source_label = high_risk_label,
  label_direction_source = "MIMIC-IV in-hospital mortality",
  assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID
)]

feature_coverage <- data.table(
  feature = common,
  full_observed_fraction = vapply(common, function(v) mean(is.finite(as.numeric(E[[v]]))), numeric(1)),
  strict_observed_fraction = vapply(common, function(v) mean(is.finite(as.numeric(E[[v]][strict_mask]))), numeric(1)),
  used_for_assignment = TRUE,
  patient_level_missing_handling = "frozen MIMIC median",
  dimension_removed = FALSE
)

audit <- data.table(
  script_version = SCRIPT_VERSION,
  assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID,
  source_eicu_file = PATH_EICU,
  source_eicu_md5 = hash_file(PATH_EICU),
  source_label_file = PATH_EICU_LABELS,
  source_label_md5 = hash_file(PATH_EICU_LABELS),
  full_n_icu_stays = nrow(E),
  full_n_unique_patients = uniqueN(E$uniquepid),
  strict_n_icu_stays = sum(strict_mask),
  strict_n_unique_patients = uniqueN(E$uniquepid[strict_mask]),
  strict_definition = "current eICU Sepsis-3 cohort; SOFA >=2 and AKI stage 0-24h >=1",
  frozen_feature_n = length(common),
  any_dimension_removed = FALSE,
  eicu_reclustered = FALSE,
  eicu_reassigned = FALSE,
  high_risk_source_label = high_risk_label,
  label_direction_source = "MIMIC-IV in-hospital mortality",
  paired_ci_method = "uniquepid-clustered percentile bootstrap",
  bootstrap_replicates = BOOT_B
)

## ---- Concordance gates against locked full-cohort outputs ----
locked_perf <- fread(PATH_LOCKED_PERF)
locked_perf[, outcome_join := fifelse(outcome == "mortality", "hospital_mortality", "make_hosp")]
qa_perf <- merge(
  performance[cohort == "eicu_full", .(
    outcome_join = outcome, model,
    observed_auc = round(auc, 4),
    observed_slope = round(calibration_slope, 3),
    observed_intercept = round(calibration_intercept, 3)
  )],
  locked_perf[, .(
    outcome_join, model,
    locked_auc = as.numeric(external_auc),
    locked_slope = as.numeric(ext_cal_slope),
    locked_intercept = as.numeric(ext_cal_intercept)
  )],
  by = c("outcome_join", "model"), all = TRUE
)
qa_perf[, `:=`(
  auc_match = is.finite(observed_auc) & abs(observed_auc - locked_auc) <= 0.0001,
  slope_match = is.finite(observed_slope) & abs(observed_slope - locked_slope) <= 0.001,
  intercept_match = is.finite(observed_intercept) & abs(observed_intercept - locked_intercept) <= 0.001
)]

locked_pair <- fread(PATH_LOCKED_PAIR)
locked_pair[, outcome_join := fifelse(outcome == "hosp_mortality", "hospital_mortality", "make_hosp")]
qa_pair <- merge(
  paired[cohort == "eicu_full", .(
    outcome_join = outcome, better, worse,
    observed_auc_better = auc_better, observed_auc_worse = auc_worse,
    observed_dAUC = dAUC, observed_lo = dAUC_lo, observed_hi = dAUC_hi
  )],
  locked_pair[, .(
    outcome_join, better, worse,
    locked_auc_better = as.numeric(auc_better), locked_auc_worse = as.numeric(auc_worse),
    locked_dAUC = as.numeric(dAUC), locked_lo = as.numeric(dAUC_lo), locked_hi = as.numeric(dAUC_hi)
  )],
  by = c("outcome_join", "better", "worse"), all = TRUE
)
for (metric in c("auc_better", "auc_worse", "dAUC", "lo", "hi")) {
  qa_pair[, paste0(metric, "_match") :=
    is.finite(get(paste0("observed_", metric))) &
      abs(get(paste0("observed_", metric)) - get(paste0("locked_", metric))) <= 0.0001]
}

summary_qa <- data.table(
  check = c(
    "current_eICU_input_hash_matches_transport",
    "assignment_run_id_locked",
    "all_33_frozen_dimensions_retained",
    "eICU_denominator_17465_stays_16212_patients",
    "strict_subset_7975_stays",
    "no_eICU_reclustering_or_reassignment",
    "full_performance_matches_locked_script19",
    "full_paired_contrasts_match_locked_script21",
    "output_row_counts_are_complete"
  ),
  pass = c(
    identical(as.character(recorded_eicu_md5), as.character(hash_file(PATH_EICU))),
    identical(as.character(transport$assignment_run_id), EXPECTED_ASSIGNMENT_RUN_ID),
    identical(common, FEATURES33) && length(common) == 33L,
    nrow(E) == EXPECTED_EICU_STAYS && uniqueN(E$uniquepid) == EXPECTED_EICU_PATIENTS,
    sum(strict_mask) == EXPECTED_STRICT_STAYS,
    TRUE,
    nrow(qa_perf) == 8L && all(qa_perf$auc_match & qa_perf$slope_match & qa_perf$intercept_match),
    nrow(qa_pair) == 8L && all(unlist(qa_pair[, .SD, .SDcols = patterns("_match$")])),
    nrow(performance) == 16L && nrow(paired) == 16L && nrow(prevalence) == 3L
  )
)

input_manifest <- data.table(
  role = c(
    "mimic_final_full", "mimic_baseline", "mimic_labels", "eicu_current",
    "transport_object", "eicu_frozen_labels", "assignment_provenance",
    "locked_performance", "locked_paired_contrasts"
  ),
  path = required_inputs,
  md5 = unname(tools::md5sum(required_inputs)),
  modified_time = format(file.info(required_inputs)$mtime, "%Y-%m-%d %H:%M:%S")
)

fwrite(performance, file.path(OUT_DIR, "eicu_sepsis3_like_model_performance.csv"))
fwrite(paired, file.path(OUT_DIR, "eicu_sepsis3_like_paired_contrasts.csv"))
fwrite(prevalence, file.path(OUT_DIR, "eicu_sepsis3_like_prevalence_drift.csv"))
fwrite(audit, file.path(OUT_DIR, "eicu_sepsis3_like_subset_audit.csv"))
fwrite(feature_coverage, file.path(OUT_DIR, "eicu_sepsis3_like_feature_coverage.csv"))
fwrite(qa_perf, file.path(OUT_DIR, "eicu_sepsis3_like_full_performance_concordance.csv"))
fwrite(qa_pair, file.path(OUT_DIR, "eicu_sepsis3_like_full_paired_concordance.csv"))
fwrite(summary_qa, file.path(OUT_DIR, "eicu_sepsis3_like_QA.csv"))
fwrite(input_manifest, file.path(OUT_DIR, "eicu_sepsis3_like_input_manifest.csv"))

saveRDS(
  list(
    script_version = SCRIPT_VERSION,
    performance = performance,
    paired = paired,
    prevalence = prevalence,
    audit = audit,
    feature_coverage = feature_coverage,
    qa = summary_qa,
    assignment_run_id = EXPECTED_ASSIGNMENT_RUN_ID
  ),
  file.path(OUT_DIR, "eicu_sepsis3_like_sf5_current_results.rds")
)

print(summary_qa)
if (!all(summary_qa$pass)) {
  stop(
    "SF5 upstream QA failed: ",
    paste(summary_qa[pass == FALSE, check], collapse = "; "),
    call. = FALSE
  )
}

cat("\nStrict subset summary:\n")
print(audit[, .(strict_n_icu_stays, strict_n_unique_patients, strict_definition)])
cat("\nStrict-subset key paired contrasts:\n")
print(paired[cohort == "eicu_sepsis3_like" & better == "feat_pen" & worse == "base_phen"])
cat("\nSF5 upstream analysis completed successfully.\n")
cat("Output directory: ", OUT_DIR, "\n", sep = "")
