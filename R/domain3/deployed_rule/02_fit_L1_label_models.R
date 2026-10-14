options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
base <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT")
stopifnot(nzchar(output_root), nzchar(base), nzchar(repo))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(output_root, repo)
run_root <- file.path(output_root, "domain3_L1")
aggregate_out <- file.path(run_root, "aggregate_outputs")
private_out <- file.path(run_root, "private_local")
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                               serialize = FALSE))
stopifnot(sha(file.path(repo, "R/domain3/deployed_rule/L1_REANALYSIS_LOCK.md")) ==
  "3BC5D702ED013192FCB3657176B112FFB3A06F88A9509DE987242A6B355545F9")
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))
paths <- c(
  mimic = file.path(base, "final_full.csv"),
  baseline = file.path(base, "baseline_covars.csv"),
  eicu = file.path(base, "eicu_external.csv"),
  mimic_labels = file.path(private_out, "mimic_L0_L1_L2_labels_local.csv"),
  eicu_labels = file.path(private_out, "eicu_L1_L2_labels_local.csv"),
  old_predictions = file.path(output_root,
    "domain3_deployable_generator_v1_20260728/private_local",
    "d3_generator_v1_external_predictions_local.csv"),
  old_performance = file.path(output_root,
    "domain3_deployable_generator_v1_20260728/tables",
    "d3_generator_v1_model_performance.csv")
)
stopifnot(all(file.exists(paths)))
expected_source_md5 <- c(
  mimic = "a3c6582cc68922787b432632f3aff8d7",
  baseline = "972e5b37ab179908d926c2cd34618129",
  eicu = "daa5589d74bdf7c09394863cec456280"
)
observed_md5 <- tolower(unname(tools::md5sum(paths[names(expected_source_md5)])))
stopifnot(identical(observed_md5, unname(expected_source_md5)))
stopifnot(sha(paths[["old_predictions"]]) ==
            "DE04917CF28434AF3A356B464736A7B79E15E8DF0BD0D1872EF74E800CECA019",
          sha(paths[["old_performance"]]) ==
            "002F19E72A5AB7842C718AAC167D8050873BE83BCDAABD020F832D8B840961BC")

base_num <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
outcomes <- list(
  mortality = list(mimic = "hospital_mortality", eicu = "hosp_mortality",
                   expected_n = 17284L, expected_events = 2940L),
  make = list(mimic = "make30", eicu = "make_hosp",
              expected_n = 17333L, expected_events = 5714L)
)

M <- fread(paths[["mimic"]],
           select = c("stay_id", "age", "gender", "sofa_score",
                      "hospital_mortality", "make30"), showProgress = FALSE)
bc <- fread(paths[["baseline"]],
            select = c("stay_id", "aki_stage_0_24h"), showProgress = FALSE)
stopifnot(!anyDuplicated(M$stay_id), !anyDuplicated(bc$stay_id))
M <- merge(M, bc, by = "stay_id", all.x = TRUE)
M[, sex := parse_sex(gender)]
M <- M[complete.cases(M[, c(base_num, "hospital_mortality", "make30"),
                        with = FALSE])]
ML <- fread(paths[["mimic_labels"]], showProgress = FALSE)
idx_m <- match(M$stay_id, ML$stay_id)
stopifnot(nrow(M) == 20049L, !anyNA(idx_m),
          !anyDuplicated(ML$stay_id))
M[, `:=`(l1_k2 = factor(ML$l1_k2[idx_m], levels = 1:2),
         l2_k2 = factor(ML$l2_k2[idx_m], levels = 1:2))]
stopifnot(sum(M$l1_k2 == "1") == 3902L,
          sum(M$l2_k2 == "1") == 3976L)
BASEm <- as.matrix(M[, ..base_num])
storage.mode(BASEm) <- "double"

E <- fread(paths[["eicu"]], showProgress = FALSE)
EL <- fread(paths[["eicu_labels"]], showProgress = FALSE)
idx_e <- match(E$patientunitstayid, EL$patientunitstayid)
stopifnot(nrow(E) == 17465L, !anyNA(idx_e),
          !anyDuplicated(E$patientunitstayid),
          !anyDuplicated(EL$patientunitstayid),
          identical(as.integer(E$hospitalid), as.integer(EL$hospitalid[idx_e])),
          as.character(E$uniquepid) == as.character(EL$uniquepid[idx_e]))
E[, `:=`(l1_k2 = factor(EL$l1_k2[idx_e], levels = 1:2),
         l2_k2 = factor(EL$l2_k2[idx_e], levels = 1:2))]
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
BASEe <- as.matrix(E[, ..base_num])
storage.mode(BASEe) <- "double"
external_complete_base <- stats::complete.cases(BASEe) & !is.na(E$l1_k2)
stopifnot(identical(is.na(E$l1_k2), is.na(E$l2_k2)))

fit_label_glm <- function(y, base_matrix, label) {
  dat <- data.frame(y = y, base_matrix, cluster_k2 = factor(label, levels = c("1", "2")))
  fit <- stats::glm(y ~ ., data = dat, family = stats::binomial())
  if (!isTRUE(fit$converged) || anyNA(stats::coef(fit))) {
    stop("Base+K2 logistic fit did not converge with finite coefficients.")
  }
  fit
}
predict_label_glm <- function(fit, base_matrix, label) {
  newdata <- data.frame(
    base_matrix,
    cluster_k2 = factor(label, levels = levels(fit$model$cluster_k2))
  )
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}
make_stratified_folds <- function(y, k) {
  folds <- integer(length(y))
  for (class_value in c(0L, 1L)) {
    ii <- which(y == class_value)
    folds[ii] <- sample(rep_len(seq_len(k), length(ii)))
  }
  folds
}
internal_cv <- function(y, base_matrix, l1, l2, seed) {
  set.seed(seed)
  folds <- make_stratified_folds(y, 10L)
  predicted <- matrix(NA_real_, nrow = length(y), ncol = 2L,
                      dimnames = list(NULL, c("L1", "L2")))
  for (rule in c("L1", "L2")) {
    label <- if (rule == "L1") l1 else l2
    for (f in 1:10) {
      train <- folds != f
      test <- folds == f
      fit <- fit_label_glm(y[train], base_matrix[train, , drop = FALSE],
                           label[train])
      predicted[test, rule] <- predict_label_glm(
        fit, base_matrix[test, , drop = FALSE], label[test])
    }
  }
  stopifnot(all(is.finite(predicted)))
  c(L1 = fast_auc(y, predicted[, "L1"]),
    L2 = fast_auc(y, predicted[, "L2"]))
}

old_predictions <- fread(paths[["old_predictions"]], showProgress = FALSE)
old_performance <- fread(paths[["old_performance"]], showProgress = FALSE)
stopifnot(all(c("patientunitstayid", "uniquepid", "outcome",
                "outcome_value", "model", "prediction", "cluster") %in%
                names(old_predictions)))
model_rows <- list()
cv_rows <- list()
old_check_rows <- list()
prediction_rows <- list()
coefficient_rows <- list()
fits <- list()

for (outcome_index in seq_along(outcomes)) {
  outcome_name <- names(outcomes)[outcome_index]
  spec <- outcomes[[outcome_index]]
  y_m <- as.integer(M[[spec$mimic]])
  y_e_all <- as.integer(E[[spec$eicu]])
  sel <- which(external_complete_base & !is.na(y_e_all))
  y_e <- y_e_all[sel]
  stopifnot(length(sel) == spec$expected_n,
            sum(y_e == 1L) == spec$expected_events,
            all(y_m %in% 0:1), all(y_e %in% 0:1))

  archived <- old_predictions[outcome == outcome_name &
                                model %in% c("base_gen", "feat_pen")]
  for (kind in c("base_gen", "feat_pen")) {
    D <- archived[model == kind]
    jj <- match(E$patientunitstayid[sel], D$patientunitstayid)
    stopifnot(nrow(D) == length(sel), !anyDuplicated(D$patientunitstayid),
              !anyNA(jj), identical(as.character(E$uniquepid[sel]),
                                     as.character(D$uniquepid[jj])),
              identical(y_e, as.integer(D$outcome_value[jj])))
    if (kind == "base_gen") {
      stopifnot(identical(as.character(E$l2_k2[sel]),
                          as.character(D$cluster[jj])))
    }
    p_archive <- as.numeric(D$prediction[jj])
    ref <- old_performance[outcome == outcome_name & model == kind]
    stopifnot(nrow(ref) == 1L, ref$n_external_stays == length(sel),
              ref$n_external_events == sum(y_e == 1L))
    cal <- calibration_slope_intercept(y_e, p_archive)
    obs_auc <- fast_auc(y_e, p_archive)
    stopifnot(isTRUE(all.equal(obs_auc, ref$external_auc, tolerance = 1e-12)),
              isTRUE(all.equal(unname(cal["slope"]),
                               ref$calibration_slope, tolerance = 1e-12)),
              isTRUE(all.equal(unname(cal["intercept"]),
                               ref$calibration_intercept, tolerance = 1e-12)))
    old_check_rows[[length(old_check_rows) + 1L]] <- data.table(
      outcome = outcome_name, model = kind, n = length(sel),
      auc_reproduced = obs_auc, slope_reproduced = unname(cal["slope"]),
      intercept_reproduced = unname(cal["intercept"]),
      matched_archived_point_estimates = TRUE
    )
  }

  cv <- internal_cv(y_m, BASEm, M$l1_k2, M$l2_k2,
                    seed = 20260728L + 9200L + outcome_index)
  cv_rows[[outcome_index]] <- data.table(
    outcome = outcome_name, rule = names(cv), internal_cv_auc = unname(cv),
    folds = 10L, fold_seed = 20260728L + 9200L + outcome_index,
    caveat = "same new folds for L1 and L2; not byte-identical to historical folds"
  )
  fit <- fit_label_glm(y_m, BASEm, M$l1_k2)
  fits[[outcome_name]] <- fit
  p_e <- predict_label_glm(fit, BASEe[sel, , drop = FALSE], E$l1_k2[sel])
  stopifnot(length(p_e) == length(sel), all(is.finite(p_e)),
            all(p_e >= 0), all(p_e <= 1))
  cal <- calibration_slope_intercept(y_e, p_e)
  model_rows[[outcome_index]] <- data.table(
    outcome = outcome_name, model = "base_l1_k2",
    n_external_stays = length(sel),
    n_external_patients = uniqueN(E$uniquepid[sel]),
    n_external_events = sum(y_e == 1L),
    internal_cv_auc = unname(cv["L1"]),
    external_auc = fast_auc(y_e, p_e),
    calibration_slope = unname(cal["slope"]),
    calibration_intercept = unname(cal["intercept"]),
    brier_score = mean((p_e - y_e)^2),
    glm_converged = isTRUE(fit$converged)
  )
  prediction_rows[[outcome_index]] <- data.table(
    patientunitstayid = E$patientunitstayid[sel],
    uniquepid = E$uniquepid[sel],
    hospitalid = E$hospitalid[sel],
    l1_k2 = as.integer(as.character(E$l1_k2[sel])),
    outcome = outcome_name,
    outcome_value = y_e,
    prediction = p_e
  )
  coefficient_rows[[outcome_index]] <- data.table(
    outcome = outcome_name, term = names(stats::coef(fit)),
    coefficient = as.numeric(stats::coef(fit))
  )
}

model_summary <- rbindlist(model_rows)
cv_summary <- rbindlist(cv_rows)
old_check <- rbindlist(old_check_rows)
new_predictions <- rbindlist(prediction_rows)
coefficients <- rbindlist(coefficient_rows)
fwrite(model_summary,
       file.path(aggregate_out, "L1_base_k2_external_point_estimates.csv"))
fwrite(cv_summary, file.path(aggregate_out, "L1_L2_paired_internal_cv.csv"))
fwrite(old_check, file.path(aggregate_out, "archived_prediction_reproduction.csv"))
fwrite(coefficients, file.path(aggregate_out, "L1_base_k2_coefficients.csv"))
fwrite(new_predictions, file.path(private_out, "L1_external_predictions_local.csv"))
saveRDS(fits, file.path(private_out, "L1_fitted_glm_models_local.rds"), version = 3)
fwrite(data.table(role = names(paths), path = unname(paths),
                  sha256 = vapply(paths, sha, character(1))),
       file.path(aggregate_out, "model_input_manifest.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(aggregate_out, "model_sessionInfo.txt"))
cat("L1_MODEL_PASS: two Base+K2 models fitted; archived Base+K2 and Raw-EN point estimates reproduced.\n")
