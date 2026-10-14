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
old_root <- file.path(output_root, "domain3_deployable_generator_v1_20260728")
old_pred_path <- file.path(old_root, "private_local",
                           "d3_generator_v1_external_predictions_local.csv")
old_perf_path <- file.path(old_root, "tables",
                           "d3_generator_v1_model_performance.csv")
old_cal_path <- file.path(old_root, "tables",
                          "d3_generator_v1_base_gen_calibration_bootstrap.csv")
old_pair_path <- file.path(old_root, "tables",
                           "d3_generator_v1_paired_auc_contrasts.csv")
old_first_path <- file.path(old_root, "tables",
                            "d3_generator_v1_first_stay_sensitivity.csv")
new_pred_path <- file.path(private_out, "L1_external_predictions_local.csv")
eicu_labels_path <- file.path(private_out, "eicu_L1_L2_labels_local.csv")
new_point_path <- file.path(aggregate_out,
                            "L1_base_k2_external_point_estimates.csv")
stopifnot(all(file.exists(c(old_pred_path, old_perf_path, old_cal_path,
                            old_pair_path, old_first_path, new_pred_path,
                            eicu_labels_path, new_point_path))))
expected_sha <- c(
  old_pred = "DE04917CF28434AF3A356B464736A7B79E15E8DF0BD0D1872EF74E800CECA019",
  old_perf = "002F19E72A5AB7842C718AAC167D8050873BE83BCDAABD020F832D8B840961BC",
  old_cal = "9137AC08F7E6642F9BD939C54243AAF650F8C9448D0E9F713BC32F1E1EAD2254",
  old_pair = "C0042E53B372395DB98ABDBA3EA4B0AFC4D33084E338DEC9EFF69C1F15D4862F"
)
actual_sha <- vapply(c(old_pred = old_pred_path, old_perf = old_perf_path,
                       old_cal = old_cal_path, old_pair = old_pair_path),
                     sha, character(1))
stopifnot(identical(unname(actual_sha), unname(expected_sha)))

old_pred <- fread(old_pred_path, showProgress = FALSE)
old_perf <- fread(old_perf_path, showProgress = FALSE)
old_cal <- fread(old_cal_path, showProgress = FALSE)
old_pair <- fread(old_pair_path, showProgress = FALSE)
old_first <- fread(old_first_path, showProgress = FALSE)
new_pred <- fread(new_pred_path, showProgress = FALSE)
all_labels <- fread(eicu_labels_path, showProgress = FALSE)
stopifnot(nrow(all_labels) == 17465L,
          !anyDuplicated(all_labels$patientunitstayid))

bootstrap_auc <- function(y, p, patient, seed, B = 1000L) {
  ii_list <- cluster_bootstrap_indices(patient, B, seed)
  values <- vapply(ii_list, function(ii) fast_auc(y[ii], p[ii]), numeric(1))
  stopifnot(all(is.finite(values)))
  ci <- percentile_interval(values)
  c(low = ci[1], high = ci[2], bootstrap_mean = mean(values),
    mcse = sd(values) / sqrt(B), valid = B)
}
bootstrap_cal <- function(y, p, ii_list) {
  values <- vapply(ii_list, function(ii) {
    calibration_slope_intercept(y[ii], p[ii])
  }, numeric(2))
  stopifnot(identical(rownames(values), c("slope", "intercept")))
  rbindlist(lapply(c("slope", "intercept"), function(term) {
    v <- as.numeric(values[term, ])
    valid <- is.finite(v)
    stopifnot(sum(valid) > 0L)
    ci <- percentile_interval(v[valid])
    data.table(parameter = term,
               bootstrap_mean = mean(v[valid]),
               bootstrap_sd = sd(v[valid]),
               bootstrap_mcse = sd(v[valid]) / sqrt(sum(valid)),
               ci_low = ci[1], ci_high = ci[2],
               bootstrap_requested = length(v),
               bootstrap_valid = sum(valid),
               bootstrap_failed = sum(!valid))
  }))
}

calibration_rows <- list()
performance_rows <- list()
paired_rows <- list()
first_rows <- list()
archive_qc_rows <- list()
outcomes <- c("mortality", "make")
first_ids <- all_labels[order(uniquepid, patientunitstayid),
                        .SD[1], by = uniquepid]$patientunitstayid
stopifnot(length(first_ids) == uniqueN(all_labels$uniquepid))

for (i in seq_along(outcomes)) {
  outcome_name <- outcomes[i]
  oldD <- old_pred[outcome == outcome_name & model == "base_gen"]
  rawD <- old_pred[outcome == outcome_name & model == "feat_pen"]
  newD <- new_pred[outcome == outcome_name]
  stopifnot(nrow(oldD) == nrow(newD), !anyDuplicated(oldD$patientunitstayid),
            !anyDuplicated(newD$patientunitstayid))
  ii_new <- match(oldD$patientunitstayid, newD$patientunitstayid)
  ii_raw <- match(oldD$patientunitstayid, rawD$patientunitstayid)
  stopifnot(!anyNA(ii_new), !anyNA(ii_raw))
  newD <- newD[ii_new]
  rawD <- rawD[ii_raw]
  stopifnot(identical(oldD$patientunitstayid, newD$patientunitstayid),
            identical(oldD$patientunitstayid, rawD$patientunitstayid),
            identical(as.character(oldD$uniquepid),
                      as.character(newD$uniquepid)),
            identical(as.character(oldD$uniquepid),
                      as.character(rawD$uniquepid)),
            identical(as.integer(oldD$outcome_value),
                      as.integer(newD$outcome_value)),
            identical(as.integer(oldD$outcome_value),
                      as.integer(rawD$outcome_value)))
  y <- as.integer(oldD$outcome_value)
  patient <- oldD$uniquepid
  p_l1 <- as.numeric(newD$prediction)
  p_l2 <- as.numeric(oldD$prediction)
  p_raw <- as.numeric(rawD$prediction)
  point <- fread(new_point_path)[outcome == outcome_name]
  ref <- old_perf[outcome == outcome_name & model == "base_gen"]
  stopifnot(nrow(point) == 1L, nrow(ref) == 1L,
            isTRUE(all.equal(fast_auc(y, p_l1), point$external_auc,
                             tolerance = 1e-12)),
            isTRUE(all.equal(fast_auc(y, p_l2), ref$external_auc,
                             tolerance = 1e-12)))

  auc_seed <- 20260728L + 9000L + i
  auc_ci <- bootstrap_auc(y, p_l1, patient, auc_seed)
  cal_point <- calibration_slope_intercept(y, p_l1)
  cal_index <- match(outcome_name, sort(outcomes))
  cal_seed <- 20260728L + 9100L + cal_index
  ii_cal <- cluster_bootstrap_indices(patient, 1000L, cal_seed)
  old_cal_recomputed <- bootstrap_cal(y, p_l2, ii_cal)
  new_cal_recomputed <- bootstrap_cal(y, p_l1, ii_cal)
  old_cal_archive <- old_cal[outcome == outcome_name & model == "base_gen"]
  for (j in seq_len(nrow(old_cal_recomputed))) {
    z <- old_cal_recomputed[j]
    archived <- old_cal_archive[parameter == z$parameter]
    stopifnot(nrow(archived) == 1L,
              isTRUE(all.equal(z$ci_low, archived$ci_low,
                               tolerance = 1e-10)),
              isTRUE(all.equal(z$ci_high, archived$ci_high,
                               tolerance = 1e-10)),
              identical(z$bootstrap_valid, archived$bootstrap_valid))
  }
  archive_qc_rows[[length(archive_qc_rows) + 1L]] <- data.table(
    outcome = outcome_name,
    item = "old L2 calibration bootstrap intervals",
    exact_numeric_reproduction = TRUE,
    old_n = nrow(oldD)
  )
  new_cal_recomputed[, `:=`(
    outcome = outcome_name, model = "base_l1_k2",
    estimate = c(unname(cal_point["slope"]),
                 unname(cal_point["intercept"])),
    bootstrap_seed = cal_seed,
    resampling_unit = "uniquepid"
  )]
  calibration_rows[[i]] <- new_cal_recomputed
  performance_rows[[i]] <- data.table(
    outcome = outcome_name, model = "base_l1_k2",
    n_external_stays = nrow(newD),
    n_external_patients = uniqueN(patient),
    n_external_events = sum(y == 1L),
    external_auc = fast_auc(y, p_l1),
    external_auc_ci_low = unname(auc_ci["low"]),
    external_auc_ci_high = unname(auc_ci["high"]),
    external_auc_bootstrap_mcse = unname(auc_ci["mcse"]),
    auc_bootstrap_seed = auc_seed,
    auc_bootstrap_n = 1000L,
    calibration_slope = unname(cal_point["slope"]),
    calibration_intercept = unname(cal_point["intercept"]),
    brier_score = mean((p_l1 - y)^2)
  )

  pair_seed <- 20260728L + 5000L +
    if (outcome_name == "mortality") 2L else 6L
  ii_pair <- cluster_bootstrap_indices(patient, 1000L, pair_seed)
  pair_old <- vapply(ii_pair, function(ii) {
    fast_auc(y[ii], p_raw[ii]) - fast_auc(y[ii], p_l2[ii])
  }, numeric(1))
  pair_new <- vapply(ii_pair, function(ii) {
    fast_auc(y[ii], p_raw[ii]) - fast_auc(y[ii], p_l1[ii])
  }, numeric(1))
  stopifnot(all(is.finite(pair_old)), all(is.finite(pair_new)))
  archived_pair <- old_pair[outcome == outcome_name &
                              model_a == "feat_pen" & model_b == "base_gen"]
  old_ci <- percentile_interval(pair_old)
  stopifnot(nrow(archived_pair) == 1L,
            isTRUE(all.equal(old_ci[1], archived_pair$ci_low,
                             tolerance = 1e-12)),
            isTRUE(all.equal(old_ci[2], archived_pair$ci_high,
                             tolerance = 1e-12)),
            isTRUE(all.equal(fast_auc(y, p_raw) - fast_auc(y, p_l2),
                             archived_pair$delta_auc, tolerance = 1e-12)))
  new_ci <- percentile_interval(pair_new)
  paired_rows[[i]] <- data.table(
    outcome = outcome_name,
    contrast = "Raw-EN minus Base+L1 K2",
    delta_auc = fast_auc(y, p_raw) - fast_auc(y, p_l1),
    ci_low = new_ci[1], ci_high = new_ci[2],
    bootstrap_mean = mean(pair_new),
    bootstrap_mcse = sd(pair_new) / sqrt(length(pair_new)),
    bootstrap_n = length(pair_new),
    bootstrap_seed = pair_seed,
    old_L2_contrast_reproduced = TRUE
  )
  archive_qc_rows[[length(archive_qc_rows) + 1L]] <- data.table(
    outcome = outcome_name,
    item = "old L2 Raw-EN paired AUC interval",
    exact_numeric_reproduction = TRUE,
    old_n = nrow(oldD)
  )

  first <- newD[patientunitstayid %in% first_ids]
  first_y <- as.integer(first$outcome_value)
  first_cal <- calibration_slope_intercept(first_y, first$prediction)
  first_rows[[i]] <- data.table(
    outcome = outcome_name, model = "base_l1_k2",
    n_stays = nrow(first), n_patients = uniqueN(first$uniquepid),
    n_events = sum(first_y == 1L),
    auc = fast_auc(first_y, first$prediction),
    calibration_slope = unname(first_cal["slope"]),
    calibration_intercept = unname(first_cal["intercept"]),
    brier_score = mean((first$prediction - first_y)^2)
  )
  old_firstD <- oldD[patientunitstayid %in% first_ids]
  old_first_ref <- old_first[outcome == outcome_name & model == "base_gen"]
  stopifnot(nrow(old_first_ref) == 1L,
            nrow(old_firstD) == old_first_ref$n_stays,
            isTRUE(all.equal(fast_auc(as.integer(old_firstD$outcome_value),
                                      old_firstD$prediction),
                             old_first_ref$auc, tolerance = 1e-12)))
}

calibration <- rbindlist(calibration_rows)
performance <- rbindlist(performance_rows)
paired <- rbindlist(paired_rows)
first_stay <- rbindlist(first_rows)
archive_qc <- rbindlist(archive_qc_rows)
fwrite(calibration,
       file.path(aggregate_out, "L1_base_k2_calibration_bootstrap.csv"))
fwrite(performance,
       file.path(aggregate_out, "L1_base_k2_external_performance.csv"))
fwrite(paired,
       file.path(aggregate_out, "L1_RawEN_paired_auc_contrast.csv"))
fwrite(first_stay,
       file.path(aggregate_out, "L1_first_recorded_stay_sensitivity.csv"))
fwrite(archive_qc,
       file.path(aggregate_out, "archived_interval_reproduction.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(aggregate_out, "bootstrap_sessionInfo.txt"))
cat("L1_BOOTSTRAP_PASS: old calibration and paired contrast intervals reproduced; new L1 intervals complete.\n")
