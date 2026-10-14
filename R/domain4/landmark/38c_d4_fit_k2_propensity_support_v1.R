source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

d4_require_stage("38b_riskset")
stage <- "38c_propensity_support"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

risk <- as.data.table(readRDS(file.path(
  D4_PRIVATE_DIR, "D4_landmark_patient_ledger_private.rds"
)))
fit_result <- d4_fit_support(risk, imputation_mode = "median_mode")

if (!fit_result$model_qc$converged ||
    fit_result$model_qc$nonfinite_coefficient_count > 0L ||
    !fit_result$model_qc$cluster_k2_in_formula ||
    fit_result$model_qc$forbidden_term_in_formula ||
    fit_result$model_qc$effect_model_fitted) {
  stop("D4 propensity model failed locked structural QC.", call. = FALSE)
}

d4_atomic_fwrite(
  fit_result$support,
  file.path(D4_TABLE_DIR, "D4_support_metrics_by_K2.csv")
)
d4_atomic_fwrite(
  fit_result$cells,
  file.path(D4_TABLE_DIR, "D4_treatment_cells_and_events_by_K2.csv")
)
d4_atomic_fwrite(
  fit_result$imputation,
  file.path(D4_TABLE_DIR, "D4_propensity_formula_and_imputation.csv")
)
d4_atomic_fwrite(
  fit_result$model_qc,
  file.path(D4_LOG_DIR, "D4_propensity_model_QC.csv")
)
d4_atomic_save_rds(
  fit_result$data[, .(
    stay_id, cluster_k2, first_event_state, event_time_tie,
    rrt_first_24_72, mortality_30d, ps, ipw_ate,
    outside_overlap_within_k2
  )],
  file.path(D4_PRIVATE_DIR, "D4_propensity_patient_audit_private.rds")
)
d4_mark_stage(stage)
d4_msg("Completed ", stage)

