options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/source_correction/00_correction_paths.R"))
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                                    serialize = FALSE))
stopifnot(identical(sha(file.path(spec_root, "D3_EXTERNAL_TEMPERATURE_CORRECTION_SPEC.md")),
                    "E6AFB6FD2634ABCA03C77FC06DDCB52B099468335F8D984790AFCA661BAD56B8"))
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))

model_path <- file.path(source_root, "private_local/L1_fitted_glm_models_local.rds")
old_pred_path <- file.path(source_root,
                           "private_local/L1_external_predictions_local.csv")
extract_path <- file.path(base, "eicu_external.csv")
labels_path <- file.path(private, "eicu_corrected_L1_labels_local.csv")
stopifnot(all(file.exists(c(model_path, old_pred_path, extract_path,
                             labels_path))),
          identical(sha(model_path),
                    "1CE13067F4864BAADE880B0F201286D125110F189CB4E01AF6489E8A4BECAEA1"))

fits <- readRDS(model_path)
old <- fread(old_pred_path, showProgress = FALSE)
E <- fread(extract_path, showProgress = FALSE)
L <- fread(labels_path, showProgress = FALSE)
stopifnot(nrow(E) == 17465L, nrow(L) == nrow(E),
          !anyDuplicated(E$patientunitstayid),
          !anyDuplicated(L$patientunitstayid))
li <- match(E$patientunitstayid, L$patientunitstayid)
stopifnot(!anyNA(li))
L <- L[li]
stopifnot(identical(E$patientunitstayid, L$patientunitstayid))

E[, age := parse_age(age)]
E[, sex := parse_sex(sex)]
setnames(E, "sofa", "sofa_score", skip_absent = TRUE)
E[, aki_stage_0_24h := as.numeric(aki_stage_0_24h)]
base_num <- c("age", "sex", "sofa_score", "aki_stage_0_24h")

predict_fit <- function(fit, base_matrix, labels) {
  newdata <- data.frame(base_matrix,
                        cluster_k2 = factor(labels,
                          levels = levels(fit$model$cluster_k2)))
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}

pred_rows <- list()
point_rows <- list()
qc_rows <- list()
for (outcome_name in c("mortality", "make")) {
  oldD <- old[outcome == outcome_name]
  stopifnot(!anyDuplicated(oldD$patientunitstayid),
            nrow(oldD) == if (outcome_name == "mortality") 17284L else 17333L)
  ii <- match(oldD$patientunitstayid, E$patientunitstayid)
  stopifnot(!anyNA(ii),
            identical(as.character(oldD$uniquepid),
                      as.character(E$uniquepid[ii])))
  base_matrix <- as.matrix(E[ii, ..base_num])
  storage.mode(base_matrix) <- "double"
  stopifnot(all(is.finite(base_matrix)))
  fit <- fits[[outcome_name]]
  p_old <- predict_fit(fit, base_matrix, L$old_l1_k2[ii])
  p_new <- predict_fit(fit, base_matrix, L$corrected_l1_k2[ii])
  stopifnot(all(is.finite(p_old)), all(is.finite(p_new)),
            max(abs(p_old - oldD$prediction)) < 1e-12)
  unchanged <- L$old_l1_k2[ii] == L$corrected_l1_k2[ii]
  stopifnot(max(abs(p_old[unchanged] - p_new[unchanged])) < 1e-12)

  y <- as.integer(oldD$outcome_value)
  expected_events <- if (outcome_name == "mortality") 2940L else 5714L
  stopifnot(sum(y == 1L) == expected_events)
  cal <- calibration_slope_intercept(y, p_new)
  point_rows[[outcome_name]] <- data.table(
    outcome = outcome_name, model = "base_l1_k2",
    n_external_stays = length(ii),
    n_external_patients = uniqueN(oldD$uniquepid),
    n_external_events = sum(y == 1L),
    external_auc = fast_auc(y, p_new),
    calibration_slope = unname(cal["slope"]),
    calibration_intercept = unname(cal["intercept"]),
    brier_score = mean((p_new - y)^2))
  qc_rows[[outcome_name]] <- data.table(
    outcome = outcome_name, stays = length(ii),
    changed_label_stays = sum(!unchanged),
    changed_prediction_stays = sum(abs(p_old - p_new) > 1e-12),
    old_prediction_max_abs_reproduction_error =
      max(abs(p_old - oldD$prediction)))
  pred_rows[[outcome_name]] <- data.table(
    patientunitstayid = E$patientunitstayid[ii],
    uniquepid = E$uniquepid[ii],
    hospitalid = E$hospitalid[ii],
    l1_k2 = L$corrected_l1_k2[ii],
    outcome = outcome_name,
    outcome_value = y,
    prediction = p_new)
}

fwrite(rbindlist(pred_rows),
       file.path(private, "corrected_external_predictions_local.csv"))
fwrite(rbindlist(point_rows),
       file.path(aggregate, "L1_base_k2_external_point_estimates.csv"))
fwrite(rbindlist(qc_rows),
       file.path(aggregate, "corrected_prediction_reproduction_qc.csv"))
print(rbindlist(point_rows))
print(rbindlist(qc_rows))
cat("CORRECTED_PREDICTIONS_COMPLETE_WITH_ARCHIVE_REPRODUCTION\n")
