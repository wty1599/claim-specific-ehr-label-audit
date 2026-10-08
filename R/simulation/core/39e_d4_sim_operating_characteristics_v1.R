if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))

rep_dt <- fread(file.path(SIM_TABLE_DIR, "D4_sim_repeat_level.csv"))
support <- fread(file.path(SIM_TABLE_DIR, "D4_sim_support_metrics_by_K2.csv"))
rep_dt[is.na(expected_state) | !nzchar(expected_state), expected_state := NA_character_]
app <- rep_dt[!is.na(expected_state)]

scenario_oc <- app[, {
  n_requested <- .N
  n_failed <- sum(failure)
  valid <- !failure & !is.na(correct_classification)
  n_valid <- sum(valid)
  correct <- sum(correct_classification[valid])
  ci <- sim_wilson(correct, n_valid)
  correct_failure_wrong <- correct
  ci_fw <- sim_wilson(correct_failure_wrong, n_requested)
  list(
    n_requested = n_requested, n_valid = n_valid, n_failed = n_failed,
    correct = correct, rate = correct / n_valid,
    wilson_low = ci[["low"]], wilson_high = ci[["high"]],
    failure_as_wrong_rate = correct_failure_wrong / n_requested,
    failure_as_wrong_wilson_low = ci_fw[["low"]],
    failure_as_wrong_wilson_high = ci_fw[["high"]]
  )
}, by = .(scenario, expected_state)]
scenario_oc[, operating_characteristic := fifelse(
  expected_state == "SUPPORT_INADEQUATE", "Sensitivity", "Specificity"
)]

overall <- app[, {
  valid <- !failure & !is.na(correct_classification)
  x <- sum(correct_classification[valid])
  n <- sum(valid)
  ci <- sim_wilson(x, n)
  ci_fw <- sim_wilson(x, .N)
  list(
    n_requested = .N, n_valid = n, n_failed = sum(failure), correct = x,
    rate = x / n, wilson_low = ci[["low"]], wilson_high = ci[["high"]],
    failure_as_wrong_rate = x / .N,
    failure_as_wrong_wilson_low = ci_fw[["low"]],
    failure_as_wrong_wilson_high = ci_fw[["high"]]
  )
}]
overall[, `:=`(
  scenario = "All applicable D4 controls",
  expected_state = "scenario-specific",
  operating_characteristic = "Correct classification"
)]
setcolorder(overall, names(scenario_oc))

metric_cols <- c(
  "treatment_prevalence", "outside_overlap_proportion", "ess_ratio",
  "minimum_cell", "extreme_weight_proportion", "ps_auc"
)
support[, (metric_cols) := lapply(.SD, as.numeric), .SDcols = metric_cols]
metric_long <- melt(
  support[phenotype %in% c("C1 higher-risk", "C2 lower-risk")],
  id.vars = c("scenario", "repeat_id", "phenotype"),
  measure.vars = metric_cols, variable.name = "metric", value.name = "value"
)
metric_summary <- metric_long[, {
  x <- value[is.finite(value)]
  list(
    n = length(x), mean = sim_safe_mean(x), sd = sim_safe_sd(x),
    p025 = sim_q(x, 0.025), p975 = sim_q(x, 0.975),
    mcse = sim_safe_sd(x) / sqrt(length(x))
  )
}, by = .(scenario, phenotype, metric)]

failure_audit <- rep_dt[, .(
  n_requested = .N, n_failed = sum(failure),
  failure_rate = mean(failure),
  ledger_qc_failure_n = sum(!ledger_qc_pass, na.rm = TRUE),
  ps_convergence_failure_n = sum(!ps_model_converged, na.rm = TRUE)
), by = scenario]
failure_reasons <- rep_dt[failure == TRUE, .N, by = .(
  scenario, failure_stage, failure_reason
)]

sim_atomic_fwrite(rbindlist(list(scenario_oc, overall), fill = TRUE), file.path(
  SIM_TABLE_DIR, "D4_sim_operating_characteristics.csv"
))
sim_atomic_fwrite(metric_summary, file.path(
  SIM_TABLE_DIR, "D4_sim_continuous_metrics_mcse.csv"
))
sim_atomic_fwrite(failure_audit, file.path(
  SIM_TABLE_DIR, "D4_sim_failure_audit.csv"
))
sim_atomic_fwrite(failure_reasons, file.path(
  SIM_TABLE_DIR, "D4_sim_failure_reasons.csv"
))
