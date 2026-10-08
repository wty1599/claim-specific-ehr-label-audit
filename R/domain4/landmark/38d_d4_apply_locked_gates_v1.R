source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

d4_require_stage("38c_propensity_support")
stage <- "38d_locked_gates"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

support <- fread(file.path(D4_TABLE_DIR, "D4_support_metrics_by_K2.csv"))
cells <- fread(file.path(
  D4_TABLE_DIR, "D4_treatment_cells_and_events_by_K2.csv"
))
model_qc <- fread(file.path(D4_LOG_DIR, "D4_propensity_model_QC.csv"))

state_result <- d4_apply_state(support, cells)
if (!nrow(state_result$gate_table) ||
    !nrow(state_result$support_alerts)) {
  stop("D4 gate application lacks required rows.", call. = FALSE)
}

claim_state <- data.table(
  domain = "D4 Actionability / HTE identifiability",
  audited_claim = "label-specific observed RRT allocation support",
  state = state_result$state,
  primary_gate = D4_PRIMARY_GATE,
  triggers = state_result$triggers,
  interface_complete = TRUE,
  effect_model_fitted = FALSE,
  causal_effect_estimated = FALSE,
  interpretation = if (state_result$state == "SUPPORT_INADEQUATE") {
    paste(
      "The hour-24 landmark analysis was evaluable, but prespecified",
      "support criteria did not permit label-specific effect modeling."
    )
  } else {
    paste(
      "The hour-24 landmark analysis met support criteria; exchangeability",
      "and treatment effects remain untested."
    )
  }
)

d4_atomic_fwrite(
  state_result$gate_table,
  file.path(D4_TABLE_DIR, "D4_gate_A_to_D_by_K2.csv")
)
d4_atomic_fwrite(
  state_result$support_alerts,
  file.path(D4_TABLE_DIR, "D4_structural_support_alerts_by_K2.csv")
)
d4_atomic_fwrite(
  claim_state,
  file.path(D4_TABLE_DIR, "D4_claim_state.csv")
)
d4_atomic_write_lines(
  c(
    paste0("state=", claim_state$state),
    paste0("triggers=", claim_state$triggers),
    paste0("interface_complete=", claim_state$interface_complete),
    "effect_model_fitted=FALSE",
    paste0("model_converged=", model_qc$converged)
  ),
  file.path(D4_LOG_DIR, "D4_claim_state_summary.txt")
)
d4_mark_stage(stage)
d4_msg("Completed ", stage, "; state=", claim_state$state)

