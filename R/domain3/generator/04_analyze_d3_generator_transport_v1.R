## =====================================================================
## Outcome-stage analysis for the already frozen regenerated generator.
## This is the first script in the chain permitted to read eICU outcomes.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
for (pkg in c("data.table", "glmnet", "ranger")) require_namespace(pkg)
library(data.table)

if (!file.exists(PATH_GENERATOR) ||
    !file.exists(PATH_EXTERNAL_LABELS) ||
    !file.exists(PATH_EXTERNAL_FREEZE)) {
  stop("Generator or frozen external-label checkpoint is missing.")
}
object <- readRDS(PATH_GENERATOR)
freeze <- fread(PATH_EXTERNAL_FREEZE)
if (hash_file(PATH_GENERATOR) != freeze$generator_file_md5[1] ||
    hash_file(PATH_EXTERNAL_LABELS) != freeze$label_file_md5[1] ||
    hash_file(PATH_EICU) != freeze$eicu_input_md5[1]) {
  stop("A frozen D3 artifact changed before outcome analysis.")
}
if (isTRUE(freeze$external_outcomes_read[1]) ||
    isTRUE(freeze$eicu_reclustered[1])) {
  stop("External label freeze provenance is invalid.")
}

set.seed(20240601L)

## ---------- Helpers ----------
fit_model <- function(kind, y, base, X, cluster) {
  switch(
    kind,
    base = stats::glm(
      y ~ ., data = data.frame(y = y, base), family = stats::binomial()
    ),
    base_gen = stats::glm(
      y ~ .,
      data = data.frame(y = y, base, cluster_k2 = factor(cluster)),
      family = stats::binomial()
    ),
    feat_pen = glmnet::cv.glmnet(
      cbind(base, X), y, family = "binomial", alpha = ENET_ALPHA,
      nfolds = 5
    ),
    feat_rf = ranger::ranger(
      x = as.data.frame(cbind(base, X)),
      y = factor(y),
      probability = TRUE,
      num.trees = RF_TREES,
      num.threads = 0,
      seed = 1L
    ),
    stop("Unknown model kind: ", kind)
  )
}

predict_model <- function(kind, fit, base, X, cluster) {
  switch(
    kind,
    base = as.numeric(stats::predict(
      fit, newdata = data.frame(base), type = "response"
    )),
    base_gen = as.numeric(stats::predict(
      fit,
      newdata = data.frame(base, cluster_k2 = factor(
        cluster, levels = levels(fit$model$cluster_k2)
      )),
      type = "response"
    )),
    feat_pen = as.numeric(stats::predict(
      fit, newx = cbind(base, X), s = "lambda.min", type = "response"
    )),
    feat_rf = as.numeric(stats::predict(
      fit, data = as.data.frame(cbind(base, X))
    )$predictions[, "1"]),
    stop("Unknown model kind: ", kind)
  )
}

make_stratified_folds <- function(y, k) {
  folds <- integer(length(y))
  for (class_value in c(0L, 1L)) {
    ii <- which(y == class_value)
    folds[ii] <- sample(rep_len(seq_len(k), length(ii)))
  }
  folds
}

internal_cv_predictions <- function(y, base, X, cluster, model_kinds) {
  folds <- make_stratified_folds(y, MODEL_CV_FOLDS)
  out <- matrix(NA_real_, nrow = length(y), ncol = length(model_kinds))
  colnames(out) <- model_kinds
  ## Preserve the iteration order of the locked head-to-head script:
  ## model outer loop, fold inner loop.
  for (kind in model_kinds) {
    for (f in seq_len(MODEL_CV_FOLDS)) {
      train <- folds != f
      test <- folds == f
      fit <- fit_model(kind, y[train], base[train, , drop = FALSE],
                       X[train, , drop = FALSE], cluster[train])
      out[test, kind] <- predict_model(
        kind, fit, base[test, , drop = FALSE], X[test, , drop = FALSE],
        cluster[test]
      )
    }
  }
  if (anyNA(out)) stop("Internal cross-validated predictions are incomplete.")
  out
}

bootstrap_auc_ci <- function(y, p, cluster_id, B) {
  cluster_id <- as.character(cluster_id)
  clusters <- unique(cluster_id)
  rows_by_cluster <- split(seq_along(cluster_id), cluster_id)
  values <- vapply(seq_len(B), function(b) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    ii <- unlist(rows_by_cluster[sampled], use.names = FALSE)
    fast_auc(y[ii], p[ii])
  }, numeric(1))
  ci <- percentile_interval(values)
  c(low = ci[1], high = ci[2], mcse = sd(values, na.rm = TRUE) /
      sqrt(sum(is.finite(values))))
}

bootstrap_dauc <- function(y, p_a, p_b, cluster_id, B, seed) {
  rows <- cluster_bootstrap_indices(cluster_id, B, seed)
  values <- vapply(rows, function(ii) {
    fast_auc(y[ii], p_a[ii]) - fast_auc(y[ii], p_b[ii])
  }, numeric(1))
  ci <- percentile_interval(values)
  c(
    delta_auc = fast_auc(y, p_a) - fast_auc(y, p_b),
    low = ci[1], high = ci[2],
    bootstrap_mean = mean(values, na.rm = TRUE),
    mcse = sd(values, na.rm = TRUE) / sqrt(sum(is.finite(values)))
  )
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  if (!n) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}

## ---------- MIMIC development data with the same frozen generator ----------
if (hash_file(PATH_BASELINE) != unname(LOCKED_MD5["baseline_covars"])) {
  stop("Locked baseline-covariate fingerprint mismatch.")
}
M <- fread(
  PATH_FINAL_FULL,
  select = c(
    "stay_id", "subject_id", "age", "gender", "sofa_score",
    "hospital_mortality", "make30", FEATURES33
  ),
  showProgress = FALSE
)
bc <- fread(
  PATH_BASELINE,
  select = c("stay_id", "aki_stage_0_24h"),
  showProgress = FALSE
)
if (anyDuplicated(M$stay_id) || anyDuplicated(bc$stay_id)) {
  stop("Duplicate MIMIC IDs in D3 outcome-stage inputs.")
}
M <- merge(M, bc, by = "stay_id", all.x = TRUE)
M[, sex := parse_sex(gender)]
M <- M[complete.cases(M[, c(
  BASE_NUM, "hospital_mortality", "make30"
), with = FALSE])]
if (nrow(M) != EXPECTED_MIMIC_N) {
  stop("MIMIC D3 model denominator differs from 20,049.")
}
source_assignment <- as.data.table(
  predict_ehr_audit_k2(object, M, id_col = "stay_id")
)
if (!identical(source_assignment$stay_id, M$stay_id)) {
  stop("Source generator assignment is not row-aligned.")
}
M[, `:=`(
  cluster_k2 = factor(
    source_assignment$assigned_cluster, levels = SOURCE_LABEL_LEVELS
  ),
  generator_missing_feature_n = source_assignment$missing_feature_n,
  generator_absolute_margin = source_assignment$absolute_margin
)]
Mprep <- apply_frozen_recipe(M, object$recipe, FEATURES33)
Xm <- Mprep$X
BASEm <- as.matrix(M[, ..BASE_NUM])
storage.mode(BASEm) <- "double"

## ---------- eICU data merged to frozen labels ----------
E <- fread(PATH_EICU, showProgress = FALSE)
labels <- fread(PATH_EXTERNAL_LABELS, showProgress = FALSE)
if (hash_file(PATH_EXTERNAL_LABELS) != freeze$label_file_md5[1]) {
  stop("Frozen external label hash changed during outcome-stage loading.")
}
E <- merge(
  E, labels[, .(
    patientunitstayid, assigned_cluster, distance_c1, distance_c2,
    signed_margin_c1, absolute_margin, missing_feature_n,
    clipped_feature_n, generator_version, generator_content_hash
  )],
  by = "patientunitstayid", all.x = FALSE, all.y = TRUE
)
if (nrow(E) != EXPECTED_EICU_STAYS ||
    uniqueN(E$uniquepid) != EXPECTED_EICU_PATIENTS ||
    uniqueN(E$hospitalid) != EXPECTED_EICU_HOSPITALS) {
  stop("eICU outcome-stage denominator mismatch.")
}
E[, age := parse_age(age)]
E[, sex := parse_sex(sex)]
setnames(E, "sofa", "sofa_score", skip_absent = TRUE)
E[, aki_stage_0_24h := as.numeric(aki_stage_0_24h)]
E[, hosp_mortality := as.integer(hosp_mortality)]
E[, make_hosp := as.integer(
  hosp_mortality == 1L |
    (!is.na(scr_discharge) & !is.na(scr_baseline) &
       scr_discharge >= 1.5 * scr_baseline)
)]
E[, cluster_k2 := factor(
  as.character(assigned_cluster), levels = SOURCE_LABEL_LEVELS
)]
Eprep <- apply_frozen_recipe(E, object$recipe, FEATURES33)
Xe <- Eprep$X
BASEe <- as.matrix(E[, ..BASE_NUM])
storage.mode(BASEe) <- "double"
external_complete_base <- stats::complete.cases(BASEe) & !is.na(E$cluster_k2)

## ---------- Prevalence and descriptive outcomes ----------
mimic_c1_n <- sum(as.character(M$cluster_k2) == SOURCE_HIGHER_RISK_LABEL)
eicu_c1_n <- sum(as.character(E$cluster_k2) == SOURCE_HIGHER_RISK_LABEL)
mimic_ci <- wilson_interval(mimic_c1_n, nrow(M))
eicu_ci <- wilson_interval(eicu_c1_n, nrow(E))
prevalence <- data.table(
  cohort = c("MIMIC-IV", "eICU"),
  label_object = c("regenerated_generator", "regenerated_generator"),
  n = c(nrow(M), nrow(E)),
  c1_n = c(mimic_c1_n, eicu_c1_n),
  c1_prevalence = c(mimic_c1_n / nrow(M), eicu_c1_n / nrow(E)),
  ci_low = c(mimic_ci[1], eicu_ci[1]),
  ci_high = c(mimic_ci[2], eicu_ci[2])
)
prevalence[, prevalence_difference_eicu_minus_mimic :=
  c(NA_real_, c1_prevalence[cohort == "eICU"] -
      c1_prevalence[cohort == "MIMIC-IV"])]
prevalence[, prevalence_ratio_eicu_over_mimic :=
  c(NA_real_, c1_prevalence[cohort == "eICU"] /
      c1_prevalence[cohort == "MIMIC-IV"])]

cluster_outcomes <- rbindlist(list(
  M[, .(
    n = .N,
    hospital_mortality_n = sum(hospital_mortality == 1L),
    hospital_mortality = mean(hospital_mortality),
    kidney_composite_n = sum(make30 == 1L),
    kidney_composite = mean(make30)
  ), by = .(cluster = as.character(cluster_k2))][
    , `:=`(cohort = "MIMIC-IV", kidney_endpoint = "MAKE30")
  ],
  E[, .(
    n = .N,
    hospital_mortality_n = sum(hosp_mortality == 1L, na.rm = TRUE),
    hospital_mortality = mean(hosp_mortality, na.rm = TRUE),
    kidney_composite_n = sum(make_hosp == 1L, na.rm = TRUE),
    kidney_composite = mean(make_hosp, na.rm = TRUE)
  ), by = .(cluster = as.character(cluster_k2))][
    , `:=`(cohort = "eICU", kidney_endpoint = "in-hospital MAKE approximation")
  ]
), use.names = TRUE)
setcolorder(cluster_outcomes, c(
  "cohort", "cluster", "n", "hospital_mortality_n",
  "hospital_mortality", "kidney_composite_n", "kidney_composite",
  "kidney_endpoint"
))

## ---------- Model fitting and external predictions ----------
MODEL_KINDS <- c("base", "base_gen", "feat_pen", "feat_rf")
OUTCOMES <- list(
  mortality = list(mimic = "hospital_mortality", eicu = "hosp_mortality"),
  make = list(mimic = "make30", eicu = "make_hosp")
)
model_summary_rows <- list()
prediction_rows <- list()
fit_store <- list()
row_counter <- 0L

for (outcome_name in names(OUTCOMES)) {
  spec <- OUTCOMES[[outcome_name]]
  y_m <- as.integer(M[[spec$mimic]])
  cv_pred <- internal_cv_predictions(
    y_m, BASEm, Xm, M$cluster_k2, MODEL_KINDS
  )
  fits <- lapply(MODEL_KINDS, function(kind) {
    fit_model(kind, y_m, BASEm, Xm, M$cluster_k2)
  })
  names(fits) <- MODEL_KINDS
  fit_store[[outcome_name]] <- fits

  sel <- which(external_complete_base & !is.na(E[[spec$eicu]]))
  y_e <- as.integer(E[[spec$eicu]][sel])
  for (kind in MODEL_KINDS) {
    p_e <- predict_model(
      kind, fits[[kind]], BASEe[sel, , drop = FALSE],
      Xe[sel, , drop = FALSE], E$cluster_k2[sel]
    )
    auc_ci <- bootstrap_auc_ci(
      y_e, p_e, E$uniquepid[sel], MODEL_BOOT_B
    )
    cal <- calibration_slope_intercept(y_e, p_e)
    internal_auc <- fast_auc(y_m, cv_pred[, kind])
    external_auc <- fast_auc(y_e, p_e)
    row_counter <- row_counter + 1L
    model_summary_rows[[row_counter]] <- data.table(
      outcome = outcome_name,
      model = kind,
      n_external_stays = length(sel),
      n_external_patients = uniqueN(E$uniquepid[sel]),
      n_external_events = sum(y_e == 1L),
      internal_cv_auc = internal_auc,
      external_auc = external_auc,
      external_auc_ci_low = auc_ci["low"],
      external_auc_ci_high = auc_ci["high"],
      external_auc_bootstrap_mcse = auc_ci["mcse"],
      transport_drop_internal_minus_external = internal_auc - external_auc,
      calibration_slope = cal["slope"],
      calibration_intercept = cal["intercept"],
      brier_score = mean((p_e - y_e)^2)
    )
    prediction_rows[[row_counter]] <- data.table(
      patientunitstayid = E$patientunitstayid[sel],
      uniquepid = E$uniquepid[sel],
      hospitalid = E$hospitalid[sel],
      cluster = as.character(E$cluster_k2[sel]),
      absolute_margin = E$absolute_margin[sel],
      outcome = outcome_name,
      outcome_value = y_e,
      model = kind,
      prediction = p_e
    )
  }
}
model_summary <- rbindlist(model_summary_rows)
predictions <- rbindlist(prediction_rows)

## ---------- Paired external AUC contrasts ----------
contrast_specs <- list(
  c("base_gen", "base"),
  c("feat_pen", "base_gen"),
  c("feat_rf", "base_gen"),
  c("feat_pen", "base")
)
contrast_rows <- list()
counter <- 0L
for (outcome_name in names(OUTCOMES)) {
  P <- dcast(
    predictions[outcome == outcome_name],
    patientunitstayid + uniquepid + outcome_value ~ model,
    value.var = "prediction"
  )
  for (pair in contrast_specs) {
    counter <- counter + 1L
    ans <- bootstrap_dauc(
      P$outcome_value, P[[pair[1]]], P[[pair[2]]], P$uniquepid,
      MODEL_BOOT_B, MASTER_SEED + 5000L + counter
    )
    contrast_rows[[counter]] <- data.table(
      outcome = outcome_name,
      model_a = pair[1],
      model_b = pair[2],
      contrast = paste0(pair[1], "_minus_", pair[2]),
      delta_auc = ans["delta_auc"],
      ci_low = ans["low"],
      ci_high = ans["high"],
      bootstrap_mean = ans["bootstrap_mean"],
      bootstrap_mcse = ans["mcse"],
      bootstrap_n = MODEL_BOOT_B
    )
  }
}
paired_contrasts <- rbindlist(contrast_rows)

## ---------- First-recorded-stay sensitivity ----------
setorder(E, uniquepid, patientunitstayid)
first_ids <- E[, .SD[1], by = uniquepid]$patientunitstayid
first_stay_summary <- rbindlist(lapply(
  names(OUTCOMES),
  function(outcome_name) {
    D <- predictions[
      outcome == outcome_name & patientunitstayid %in% first_ids
    ]
    D[, {
      cal <- calibration_slope_intercept(outcome_value, prediction)
      .(
        n_stays = .N,
        n_patients = uniqueN(uniquepid),
        n_events = sum(outcome_value == 1L),
        auc = fast_auc(outcome_value, prediction),
        calibration_slope = cal["slope"],
        calibration_intercept = cal["intercept"],
        brier_score = mean((prediction - outcome_value)^2)
      )
    }, by = .(outcome, model)]
  }
))

first_prevalence <- E[
  patientunitstayid %in% first_ids,
  .(
    n = .N,
    c1_n = sum(as.character(cluster_k2) == SOURCE_HIGHER_RISK_LABEL),
    c1_prevalence = mean(
      as.character(cluster_k2) == SOURCE_HIGHER_RISK_LABEL
    )
  )
]
first_prevalence[, selection_rule :=
  "minimum patientunitstayid per uniquepid; deterministic first-recorded-stay surrogate"]

## ---------- Hospital-aware prevalence and model heterogeneity ----------
hospital_prevalence <- E[, .(
  n_stays = .N,
  n_patients = uniqueN(uniquepid),
  c1_n = sum(as.character(cluster_k2) == SOURCE_HIGHER_RISK_LABEL),
  c1_prevalence = mean(
    as.character(cluster_k2) == SOURCE_HIGHER_RISK_LABEL
  )
), by = hospitalid]
hospital_prevalence[, c("ci_low", "ci_high") := {
  ci <- wilson_interval(c1_n, n_stays)
  list(ci[1], ci[2])
}, by = hospitalid]

hospital_prevalence_summary <- hospital_prevalence[, .(
  hospitals = .N,
  total_stays = sum(n_stays),
  median_prevalence = median(c1_prevalence),
  q25_prevalence = quantile(c1_prevalence, 0.25),
  q75_prevalence = quantile(c1_prevalence, 0.75),
  min_prevalence = min(c1_prevalence),
  max_prevalence = max(c1_prevalence),
  weighted_prevalence = sum(c1_n) / sum(n_stays)
)]

hospital_auc <- predictions[, {
  events <- sum(outcome_value == 1L)
  nonevents <- sum(outcome_value == 0L)
  .(
    n = .N,
    events = events,
    nonevents = nonevents,
    estimable = events > 0L && nonevents > 0L,
    auc = if (events > 0L && nonevents > 0L) {
      fast_auc(outcome_value, prediction)
    } else {
      NA_real_
    }
  )
}, by = .(hospitalid, outcome, model)]

hospital_auc_summary <- hospital_auc[, .(
  hospitals_total = .N,
  hospitals_estimable = sum(estimable),
  median_auc = median(auc, na.rm = TRUE),
  q25_auc = quantile(auc, 0.25, na.rm = TRUE),
  q75_auc = quantile(auc, 0.75, na.rm = TRUE),
  min_auc = min(auc, na.rm = TRUE),
  max_auc = max(auc, na.rm = TRUE)
), by = .(outcome, model)]

## Cluster-hospital bootstrap for pooled prevalence and AUC.
hospital_ids <- unique(E$hospitalid)
hospital_rows_e <- split(seq_len(nrow(E)), E$hospitalid)
set.seed(MASTER_SEED + 7000L)
hospital_boot_prevalence <- numeric(HOSPITAL_BOOT_B)
for (b in seq_len(HOSPITAL_BOOT_B)) {
  sampled <- sample(hospital_ids, length(hospital_ids), replace = TRUE)
  ii <- unlist(hospital_rows_e[as.character(sampled)], use.names = FALSE)
  hospital_boot_prevalence[b] <- mean(
    as.character(E$cluster_k2[ii]) == SOURCE_HIGHER_RISK_LABEL
  )
}
hp_ci <- percentile_interval(hospital_boot_prevalence)
hospital_prevalence_bootstrap <- data.table(
  estimate = mean(as.character(E$cluster_k2) == SOURCE_HIGHER_RISK_LABEL),
  bootstrap_mean = mean(hospital_boot_prevalence),
  bootstrap_sd = sd(hospital_boot_prevalence),
  bootstrap_mcse = sd(hospital_boot_prevalence) /
    sqrt(HOSPITAL_BOOT_B),
  ci_low = hp_ci[1],
  ci_high = hp_ci[2],
  bootstrap_n = HOSPITAL_BOOT_B,
  resampling_unit = "hospital"
)

hospital_boot_auc_rows <- list()
counter <- 0L
for (outcome_name in names(OUTCOMES)) {
  for (kind in MODEL_KINDS) {
    counter <- counter + 1L
    D <- predictions[outcome == outcome_name & model == kind]
    rows_h <- split(seq_len(nrow(D)), D$hospitalid)
    hospitals_d <- names(rows_h)
    set.seed(MASTER_SEED + 8000L + counter)
    vals <- numeric(HOSPITAL_BOOT_B)
    for (b in seq_len(HOSPITAL_BOOT_B)) {
      sampled <- sample(hospitals_d, length(hospitals_d), replace = TRUE)
      ii <- unlist(rows_h[sampled], use.names = FALSE)
      vals[b] <- fast_auc(D$outcome_value[ii], D$prediction[ii])
    }
    ci <- percentile_interval(vals)
    hospital_boot_auc_rows[[counter]] <- data.table(
      outcome = outcome_name,
      model = kind,
      estimate = fast_auc(D$outcome_value, D$prediction),
      bootstrap_mean = mean(vals, na.rm = TRUE),
      bootstrap_sd = sd(vals, na.rm = TRUE),
      bootstrap_mcse = sd(vals, na.rm = TRUE) /
        sqrt(sum(is.finite(vals))),
      ci_low = ci[1],
      ci_high = ci[2],
      bootstrap_n = sum(is.finite(vals)),
      resampling_unit = "hospital"
    )
  }
}
hospital_auc_bootstrap <- rbindlist(hospital_boot_auc_rows)

## ---------- Margin sensitivity ----------
E[, margin_quintile := cut(
  absolute_margin,
  breaks = unique(quantile(
    absolute_margin, probs = seq(0, 1, 0.2), na.rm = TRUE
  )),
  include.lowest = TRUE,
  ordered_result = TRUE
)]
margin_summary <- E[, .(
  n = .N,
  c1_prevalence = mean(
    as.character(cluster_k2) == SOURCE_HIGHER_RISK_LABEL
  ),
  hospital_mortality = mean(hosp_mortality, na.rm = TRUE),
  make_hosp = mean(make_hosp, na.rm = TRUE),
  mean_missing_features = mean(missing_feature_n),
  mean_clipped_features = mean(clipped_feature_n),
  margin_min = min(absolute_margin),
  margin_max = max(absolute_margin)
), by = margin_quintile]

## ---------- Missingness transport qualification ----------
missingness_band <- function(x) {
  factor(
    fifelse(x == 0L, "0",
      fifelse(x <= 2L, "1-2", fifelse(x <= 5L, "3-5", "6+"))
    ),
    levels = c("0", "1-2", "3-5", "6+")
  )
}
missingness_transport <- rbindlist(list(
  M[, .(
    n = .N,
    proportion = .N / nrow(M)
  ), by = .(
    missingness_stratum = missingness_band(generator_missing_feature_n)
  )][, cohort := "MIMIC-IV"],
  E[, .(
    n = .N,
    proportion = .N / nrow(E)
  ), by = .(
    missingness_stratum = missingness_band(missing_feature_n)
  )][, cohort := "eICU"]
))
internal_missingness_fidelity <- fread(file.path(
  OUT_TABLES, "internal_fidelity_by_missingness.csv"
))
missingness_transport <- merge(
  missingness_transport,
  internal_missingness_fidelity[, .(
    missingness_stratum, internal_agreement = agreement,
    internal_ari = ari, internal_c1_recall = recall_c1
  )],
  by = "missingness_stratum", all.x = TRUE
)
missingness_transport[, interpretation := fifelse(
  cohort == "eICU" & missingness_stratum == "6+",
  paste(
    "External patients are concentrated in a source-data stratum with",
    "lower internal surrogate fidelity; external semantic fidelity is unknown."
  ),
  "Descriptive missingness-distribution comparison."
)]

## ---------- D3 component alerts ----------
## Only the baseline + regenerated-label model is a generator-related
## outcome-transport object. Raw-feature models are contextual comparators.
generator_transport <- model_summary[model == "base_gen"]
internal_cv <- fread(file.path(
  OUT_TABLES, "internal_fidelity_crossfit_summary.csv"
))
internal_cv_ari <- internal_cv[metric == "ari", mean]
alert_components <- rbindlist(list(
  data.table(
    component = "prevalence_drift",
    outcome = NA_character_,
    observed = abs(
      prevalence[cohort == "eICU", c1_prevalence] -
        prevalence[cohort == "MIMIC-IV", c1_prevalence]
    ),
    lower_gate = NA_real_,
    upper_gate = D3_PREVALENCE_DRIFT_GATE,
    alert = abs(
      prevalence[cohort == "eICU", c1_prevalence] -
        prevalence[cohort == "MIMIC-IV", c1_prevalence]
    ) > D3_PREVALENCE_DRIFT_GATE,
    evaluability = "evaluated_same_generator_both_cohorts",
    drives_claim_state = TRUE
  ),
  data.table(
    component = "same_distribution_assignment_ari",
    outcome = NA_character_,
    observed = internal_cv_ari,
    lower_gate = D3_ASSIGNMENT_ARI_GATE,
    upper_gate = 1,
    alert = internal_cv_ari < D3_ASSIGNMENT_ARI_GATE,
    evaluability = "internal_reference_qualification_not_external",
    drives_claim_state = FALSE
  ),
  data.table(
    component = "external_semantic_assignment_ari",
    outcome = NA_character_,
    observed = NA_real_,
    lower_gate = D3_ASSIGNMENT_ARI_GATE,
    upper_gate = 1,
    alert = NA,
    evaluability = "REFERENCE_ABSENT",
    drives_claim_state = FALSE
  ),
  generator_transport[, .(
    component = "calibration_slope",
    outcome,
    observed = calibration_slope,
    lower_gate = D3_CALIBRATION_SLOPE_RANGE[1],
    upper_gate = D3_CALIBRATION_SLOPE_RANGE[2],
    alert = calibration_slope < D3_CALIBRATION_SLOPE_RANGE[1] |
      calibration_slope > D3_CALIBRATION_SLOPE_RANGE[2],
    evaluability = fifelse(
      outcome == "mortality",
      "evaluated_same_generator_same_endpoint",
      "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY"
    ),
    drives_claim_state = outcome == "mortality"
  )],
  generator_transport[, .(
    component = "absolute_calibration_intercept",
    outcome,
    observed = abs(calibration_intercept),
    lower_gate = 0,
    upper_gate = D3_CALIBRATION_INTERCEPT_GATE,
    alert = abs(calibration_intercept) >
      D3_CALIBRATION_INTERCEPT_GATE,
    evaluability = fifelse(
      outcome == "mortality",
      "evaluated_same_generator_same_endpoint",
      "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY"
    ),
    drives_claim_state = outcome == "mortality"
  )],
  data.table(
    component = "exact_original_assignment_transport",
    outcome = NA_character_,
    observed = NA_real_,
    lower_gate = NA_real_,
    upper_gate = NA_real_,
    alert = NA,
    evaluability = "GENERATOR_ABSENT",
    drives_claim_state = FALSE
  )
), use.names = TRUE, fill = TRUE)

d3_state <- data.table(
  claim = c(
    "exact original discovery-generator transport",
    "regenerated generator internal fidelity",
    "regenerated generator technical external application",
    "regenerated generator external semantic assignment fidelity",
    "regenerated generator prevalence transport",
    "regenerated generator mortality-model transport",
    "regenerated generator kidney-composite transport",
    "raw-feature outcome-model transport"
  ),
  object = c(
    "original MICE-imputation-1 discovery workflow",
    rep(GENERATOR_NAME, 6),
    "raw-feature prediction models"
  ),
  state = c(
    "GENERATOR_ABSENT",
    "EVALUATED_FIDELITY_QUANTIFIED",
    "EVALUATED_TECHNICAL_APPLICATION",
    "REFERENCE_ABSENT",
    if (alert_components[
      component == "prevalence_drift", alert
    ]) "EVALUATED_WITH_TRANSPORT_ALERT" else "EVALUATED_NO_ALERT",
    if (any(alert_components[
      drives_claim_state & outcome == "mortality" & !is.na(alert), alert
    ])) "EVALUATED_WITH_TRANSPORT_ALERT" else "EVALUATED_NO_ALERT",
    "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY",
    "CONTEXT_ONLY_NOT_GENERATOR_TRANSPORT"
  ),
  interpretation = c(
    "No frozen inductive new-patient mapper was instantiated by the original workflow.",
    "Agreement and ARI are reported without a retrospective scientific pass threshold.",
    "Frozen labels were generated without eICU re-clustering or outcome access.",
    "No independent eICU reference partition exists for external assignment ARI.",
    "Both cohort prevalences were generated by the same frozen object.",
    paste(
      "The baseline-plus-generated-label mortality model used the same",
      "generator in both cohorts; prespecified calibration gates determine the state."
    ),
    paste(
      "The external in-hospital kidney composite is not endpoint-equivalent",
      "to development MAKE30 and is reported descriptively."
    ),
    paste(
      "Raw-EN and Raw-RF evaluate physiology-model transport and do not",
      "determine the regenerated-generator claim state."
    )
  )
)

## ---------- Regression comparison with current authoritative D3 table ----------
reference_path <- file.path(
  PROJECT_ROOT, "output", "qc", "headtohead_summary_mice_labels.csv"
)
regression_comparison <- data.table()
if (file.exists(reference_path)) {
  ref <- fread(reference_path)
  ref_key <- ref[, .(
    outcome, model,
    ref_internal_cv_auc = internal_cv_auc,
    ref_external_auc = external_auc,
    ref_calibration_slope = ext_cal_slope,
    ref_calibration_intercept = ext_cal_intercept
  )]
  model_summary[, reference_model := fifelse(
    model == "base_gen", "base_phen", model
  )]
  regression_comparison <- merge(
    model_summary,
    ref_key,
    by.x = c("outcome", "reference_model"),
    by.y = c("outcome", "model"),
    all = TRUE
  )
  regression_comparison[, `:=`(
    diff_internal_cv_auc = internal_cv_auc - ref_internal_cv_auc,
    diff_external_auc = external_auc - ref_external_auc,
    diff_calibration_slope = calibration_slope - ref_calibration_slope,
    diff_calibration_intercept =
      calibration_intercept - ref_calibration_intercept
  )]
}

## ---------- Write outputs ----------
write_replace_csv(
  prevalence, file.path(OUT_TABLES, "d3_generator_v1_prevalence.csv")
)
write_replace_csv(
  cluster_outcomes,
  file.path(OUT_TABLES, "d3_generator_v1_cluster_outcomes.csv")
)
write_replace_csv(
  model_summary,
  file.path(OUT_TABLES, "d3_generator_v1_model_performance.csv")
)
write_replace_csv(
  model_summary[, .(
    outcome, model, calibration_slope, calibration_intercept, brier_score
  )],
  file.path(OUT_TABLES, "d3_generator_v1_calibration.csv")
)
write_replace_csv(
  paired_contrasts,
  file.path(OUT_TABLES, "d3_generator_v1_paired_auc_contrasts.csv")
)
write_replace_csv(
  first_stay_summary,
  file.path(OUT_TABLES, "d3_generator_v1_first_stay_sensitivity.csv")
)
write_replace_csv(
  first_prevalence,
  file.path(OUT_TABLES, "d3_generator_v1_first_stay_prevalence.csv")
)
write_replace_csv(
  hospital_prevalence,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_prevalence.csv")
)
write_replace_csv(
  hospital_prevalence_summary,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_prevalence_summary.csv")
)
write_replace_csv(
  hospital_prevalence_bootstrap,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_prevalence_bootstrap.csv")
)
write_replace_csv(
  hospital_auc,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_auc.csv")
)
write_replace_csv(
  hospital_auc_summary,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_auc_summary.csv")
)
write_replace_csv(
  hospital_auc_bootstrap,
  file.path(OUT_TABLES, "d3_generator_v1_hospital_auc_bootstrap.csv")
)
write_replace_csv(
  margin_summary,
  file.path(OUT_TABLES, "d3_generator_v1_margin_sensitivity.csv")
)
write_replace_csv(
  missingness_transport,
  file.path(OUT_TABLES, "d3_generator_v1_missingness_transport.csv")
)
write_replace_csv(
  alert_components,
  file.path(OUT_TABLES, "d3_generator_v1_alert_components.csv")
)
write_replace_csv(
  d3_state, file.path(OUT_TABLES, "d3_generator_v1_claim_states.csv")
)
write_replace_csv(
  regression_comparison,
  file.path(OUT_QC, "d3_generator_v1_headtohead_regression_comparison.csv")
)
write_replace_csv(
  predictions,
  file.path(OUT_PRIVATE, "d3_generator_v1_external_predictions_local.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_LOGS, "d3_generator_v1_transport_sessionInfo.txt"),
  useBytes = TRUE
)

analysis_provenance <- data.table(
  item = c(
    "generator_file_md5", "external_label_file_md5", "external_label_frozen_utc",
    "external_outcomes_read", "outcomes_read_after_label_freeze",
    "eicu_reclustered", "mimic_n", "eicu_stays", "eicu_patients",
    "eicu_hospitals", "source_label_object", "external_label_object",
    "external_assignment_reference"
  ),
  value = c(
    hash_file(PATH_GENERATOR), hash_file(PATH_EXTERNAL_LABELS),
    freeze$frozen_utc[1], "true", "true", "false", nrow(M), nrow(E),
    uniqueN(E$uniquepid), uniqueN(E$hospitalid),
    GENERATOR_NAME, GENERATOR_NAME, "absent"
  )
)
write_replace_csv(
  analysis_provenance,
  file.path(OUT_QC, "d3_generator_v1_outcome_analysis_provenance.csv")
)

cat("D3 regenerated-generator transport analysis completed.\n")
cat("\nPrevalence:\n")
print(prevalence)
cat("\nModel summary:\n")
print(model_summary)
cat("\nD3 claim-specific states:\n")
print(d3_state)
