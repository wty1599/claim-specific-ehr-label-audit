source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

d4_require_stage("38f_release")
stage <- "38g_independent_qc"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

flow <- fread(file.path(D4_TABLE_DIR, "D4_24h_landmark_flow.csv"))
events <- fread(file.path(
  D4_TABLE_DIR, "D4_24_72_first_event_by_K2.csv"
))
support <- fread(file.path(D4_TABLE_DIR, "D4_support_metrics_by_K2.csv"))
cells <- fread(file.path(
  D4_TABLE_DIR, "D4_treatment_cells_and_events_by_K2.csv"
))
gates_saved <- fread(file.path(D4_TABLE_DIR, "D4_gate_A_to_D_by_K2.csv"))
alerts_saved <- fread(file.path(
  D4_TABLE_DIR, "D4_structural_support_alerts_by_K2.csv"
))
state_saved <- fread(file.path(D4_TABLE_DIR, "D4_claim_state.csv"))
model_qc <- fread(file.path(D4_LOG_DIR, "D4_propensity_model_QC.csv"))
source_qc <- fread(file.path(D4_TABLE_DIR, "D4_source_and_hash_audit.csv"))
risk_qc <- fread(file.path(D4_LOG_DIR, "D4_24h_landmark_QC.csv"))

recomputed <- d4_apply_state(support, cells)
key_gate <- c(
  "gate", "phenotype", "n_treated", "n_control",
  "events_treated", "events_control", "pass"
)
setorderv(gates_saved, c("gate", "phenotype"))
setorderv(recomputed$gate_table, c("gate", "phenotype"))
gate_equal <- identical(
  gates_saved[, ..key_gate],
  recomputed$gate_table[, ..key_gate]
)
key_alert <- c(
  "phenotype", "outside_overlap_alert", "ess_ratio_alert",
  "minimum_cell_alert", "extreme_weight_alert"
)
setorderv(alerts_saved, "phenotype")
setorderv(recomputed$support_alerts, "phenotype")
alert_equal <- identical(
  alerts_saved[, ..key_alert],
  recomputed$support_alerts[, ..key_alert]
)

script_files <- file.path(
  D4_SCRIPT_DIR,
  c(
    "38_d4_landmark_support_common_v1.R",
    "38a_d4_landmark_source_audit_v1.R",
    "38b_d4_build_24h_opportunity_riskset_v1.R",
    "38c_d4_fit_k2_propensity_support_v1.R",
    "38d_d4_apply_locked_gates_v1.R",
    "38e_d4_landmark_sensitivity_v1.R",
    "38f_d4_release_and_manuscript_tables_v1.R"
  )
)
script_text <- paste(
  unlist(lapply(script_files, readLines, warn = FALSE)),
  collapse = "\n"
)
prohibited_pattern <- paste(
  c(
    "effect_model_fitted\\s*=\\s*TRUE",
    "rrt_first_24_72\\s*\\*\\s*cluster_k2",
    "cluster_k2\\s*\\*\\s*rrt_first_24_72",
    "EValue::",
    "evalues?\\s*\\("
  ),
  collapse = "|"
)
prohibited_code_found <- grepl(
  prohibited_pattern, script_text, ignore.case = FALSE, perl = TRUE
)

final_risk_n <- tail(flow$n_remaining, 1)
qc <- data.table(
  check = c(
    "all_source_QC_pass",
    "all_riskset_QC_pass",
    "event_ledger_reconciles",
    "treatment_cells_reconcile",
    "K2_in_propensity_formula",
    "no_forbidden_propensity_term",
    "propensity_converged",
    "finite_coefficients",
    "effect_model_not_fitted",
    "gate_recomputation_matches",
    "support_alert_recomputation_matches",
    "claim_state_recomputation_matches",
    "no_prohibited_effect_code",
    "formal_state_not_interface_incomplete"
  ),
  pass = c(
    all(source_qc$pass),
    all(risk_qc$pass),
    sum(events$n) == final_risk_n,
    sum(cells$n) == model_qc$n_model,
    model_qc$cluster_k2_in_formula,
    !model_qc$forbidden_term_in_formula,
    model_qc$converged,
    model_qc$nonfinite_coefficient_count == 0L,
    !model_qc$effect_model_fitted && !state_saved$effect_model_fitted,
    gate_equal,
    alert_equal,
    identical(
      c(state_saved$state, state_saved$triggers),
      c(recomputed$state, recomputed$triggers)
    ),
    !prohibited_code_found,
    state_saved$state != "INTERFACE_INCOMPLETE"
  ),
  observed = c(
    paste(sum(source_qc$pass), nrow(source_qc), sep = "/"),
    paste(sum(risk_qc$pass), nrow(risk_qc), sep = "/"),
    paste(sum(events$n), final_risk_n, sep = "/"),
    paste(sum(cells$n), model_qc$n_model, sep = "/"),
    as.character(model_qc$cluster_k2_in_formula),
    as.character(model_qc$forbidden_term_in_formula),
    as.character(model_qc$converged),
    as.character(model_qc$nonfinite_coefficient_count),
    as.character(state_saved$effect_model_fitted),
    as.character(gate_equal),
    as.character(alert_equal),
    paste(
      paste(state_saved$state, state_saved$triggers, sep = ": "),
      paste(recomputed$state, recomputed$triggers, sep = ": "),
      sep = " / "
    ),
    as.character(prohibited_code_found),
    state_saved$state
  ),
  expected = c(
    "all", "all", "equal", "equal", "TRUE", "FALSE", "TRUE", "0",
    "FALSE", "TRUE", "TRUE", "equal", "FALSE",
    "SUPPORT_INADEQUATE or SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT"
  )
)
d4_atomic_fwrite(
  qc,
  file.path(D4_LOG_DIR, "D4_independent_QC.csv")
)
if (any(!qc$pass)) {
  d4_atomic_write_lines(
    c(
      "D4 FORMAL QC FAILURE - DO NOT USE",
      paste0("failed=", paste(qc[pass == FALSE, check], collapse = "; "))
    ),
    file.path(D4_OUTPUT_ROOT, "D4_QC_FAILURE_DO_NOT_USE.txt")
  )
  stop(
    "Independent D4 QC failed: ",
    paste(qc[pass == FALSE, check], collapse = "; "),
    call. = FALSE
  )
}

files <- list.files(D4_OUTPUT_ROOT, recursive = TRUE, full.names = TRUE)
files <- files[
  !file.info(files)$isdir &
    !grepl(
      "D4_output_sha256_manifest.csv$|D4_FORMAL_RUN_COMPLETED.ok$",
      files
    )
]
output_manifest <- rbindlist(lapply(files, function(pth) {
  data.table(
    path = sub(
      paste0("^", gsub("\\\\", "/", normalizePath(D4_OUTPUT_ROOT, winslash = "/")), "/?"),
      "",
      normalizePath(pth, winslash = "/", mustWork = TRUE)
    ),
    bytes = file.info(pth)$size,
    sha256 = d4_sha256(pth),
    public_release = !grepl("private_not_for_release", pth, fixed = TRUE)
  )
}))
d4_atomic_fwrite(
  output_manifest,
  file.path(D4_PROVENANCE_DIR, "D4_output_sha256_manifest.csv")
)
d4_mark_stage(stage)
d4_atomic_write_lines(
  c(
    "D4 formal run completed and passed independent QC.",
    paste0("run_mode=", D4_RUN_MODE),
    paste0("state=", state_saved$state),
    paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  file.path(D4_OUTPUT_ROOT, "D4_FORMAL_RUN_COMPLETED.ok")
)
d4_msg("Completed ", stage, "; release state=", state_saved$state)
