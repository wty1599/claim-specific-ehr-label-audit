if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))
if (!exists("d4_apply_state")) sim_source_empirical_engine()

rep_dt <- fread(file.path(SIM_TABLE_DIR, "D4_sim_repeat_level.csv"))
support <- fread(file.path(SIM_TABLE_DIR, "D4_sim_support_metrics_by_K2.csv"))
gates <- fread(file.path(SIM_TABLE_DIR, "D4_sim_gate_A_to_D_by_K2.csv"))
equiv <- fread(file.path(SIM_LOG_DIR, "D4_empirical_ledger_equivalence.csv"))

make_boundary_support <- function(
    outside = 0.10, ess = 0.25, minimum_cell = 30L, extreme = 0.01) {
  data.table(
    phenotype = c("C1 higher-risk", "C2 lower-risk"), n = 10000L,
    n_treated = 100L, n_control = 9900L, treatment_prevalence = 0.01,
    common_support_low = 0.01, common_support_high = 0.99,
    outside_overlap_proportion = outside, ate_ipw_ess = ess * 10000,
    ess_ratio = ess, minimum_cell = as.integer(minimum_cell),
    extreme_weight_proportion = extreme, ps_auc = 0.50
  )
}
make_boundary_cells <- function(n_treated = 100L, n_control = 9900L,
                                events_treated = 10L, events_control = 10L) {
  rbindlist(lapply(c("1", "2"), function(k) data.table(
    cluster_k2 = k,
    phenotype = if (k == "1") "C1 higher-risk" else "C2 lower-risk",
    rrt_first_24_72 = c(1L, 0L), n = c(n_treated, n_control),
    mortality_30d_events = c(events_treated, events_control),
    mortality_30d_missing = 0L
  )))
}
boundary_cases <- list(
  exact_thresholds = list(
    support = make_boundary_support(), cells = make_boundary_cells(),
    expected = "SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT"
  ),
  outside_above = list(
    support = make_boundary_support(outside = 0.1000001),
    cells = make_boundary_cells(), expected = "SUPPORT_INADEQUATE"
  ),
  ess_below = list(
    support = make_boundary_support(ess = 0.2499999),
    cells = make_boundary_cells(), expected = "SUPPORT_INADEQUATE"
  ),
  minimum_cell_below = list(
    support = make_boundary_support(minimum_cell = 29L),
    cells = make_boundary_cells(), expected = "SUPPORT_INADEQUATE"
  ),
  extreme_weight_above = list(
    support = make_boundary_support(extreme = 0.0100001),
    cells = make_boundary_cells(), expected = "SUPPORT_INADEQUATE"
  ),
  gate_fraction_below = list(
    support = make_boundary_support(),
    cells = make_boundary_cells(n_treated = 99L, n_control = 9901L),
    expected = "SUPPORT_INADEQUATE"
  ),
  gate_events_below = list(
    support = make_boundary_support(),
    cells = make_boundary_cells(events_treated = 9L),
    expected = "SUPPORT_INADEQUATE"
  )
)
boundary_qc <- rbindlist(lapply(names(boundary_cases), function(nm) {
  z <- boundary_cases[[nm]]
  observed <- d4_apply_state(z$support, z$cells)$state
  data.table(
    check = paste0("threshold_boundary_", nm), pass = observed == z$expected,
    detail = paste("observed", observed, "expected", z$expected)
  )
}))

phen <- support[phenotype %in% c("C1 higher-risk", "C2 lower-risk")]
support_check <- phen[, .(
  support_failure = any(
    outside_overlap_proportion > 0.10 |
      ess_ratio < 0.25 |
      minimum_cell < 30 |
      extreme_weight_proportion > 0.01
  )
), by = .(scenario, repeat_id)]
gate_check <- gates[gate == "B", .(
  gate_b_failure = any(!pass)
), by = .(scenario, repeat_id)]
ind <- merge(support_check, gate_check, by = c("scenario", "repeat_id"), all = TRUE)
ind[, independently_reconstructed_state := fifelse(
  support_failure | gate_b_failure,
  "SUPPORT_INADEQUATE", "SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT"
)]
cmp <- merge(
  rep_dt[failure == FALSE, .(scenario, repeat_id, observed_state)],
  ind, by = c("scenario", "repeat_id"), all.x = TRUE
)

expected_rows <- sum(SIM_REPS)
skip_private_qc <- identical(
  Sys.getenv("D4_SKIP_PRIVATE_EQUIVALENCE_QC", unset = "1"), "1"
)
equiv_pass_or_disclosed_skip <-
  (!anyNA(equiv$pass) && all(equiv$pass)) ||
  (skip_private_qc && all(equiv$check == "empirical_authority_files_available"))
qcs <- data.table(
  check = c(
    "repeat_count", "empirical_ledger_equivalence",
    "independent_state_reconstruction", "all_successful_ledgers_pass_qc",
    "all_successful_ps_models_converged", "ledger_stress_created_ties",
    "failure_stage_audit", "g4l_expected_state_defined",
    "g6l_expected_state_defined"
  ),
  pass = c(
    nrow(rep_dt) == expected_rows,
    equiv_pass_or_disclosed_skip,
    all(cmp$observed_state == cmp$independently_reconstructed_state),
    all(rep_dt[failure == FALSE, ledger_qc_pass]),
    all(rep_dt[failure == FALSE, ps_model_converged]),
    all(rep_dt[scenario == SIM_SCENARIOS[["D4LQ"]] & failure == FALSE, n_ties > 0]),
    all(rep_dt[failure == TRUE, failure_stage] %in% c(
      "generation", "landmark_ledger", "propensity_and_imputation",
      "support_metrics_and_state"
    )),
    identical(unname(SIM_EXPECTED_STATE[[SIM_SCENARIOS[["G4L"]]]]),
              "SUPPORT_INADEQUATE"),
    identical(unname(SIM_EXPECTED_STATE[[SIM_SCENARIOS[["G6L"]]]]),
              "SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT")
  ),
  detail = c(
    paste(nrow(rep_dt), "of", expected_rows),
    paste0(
      if (skip_private_qc && !all(equiv$pass)) "EXPLICIT PUBLIC-RELEASE SKIP | " else "",
      paste(equiv$check, equiv$detail, sep = ": ", collapse = " | ")
    ),
    paste(sum(cmp$observed_state == cmp$independently_reconstructed_state),
          "of", nrow(cmp)),
    paste(sum(rep_dt[failure == FALSE, ledger_qc_pass]), "successful repeats"),
    paste(sum(rep_dt[failure == FALSE, ps_model_converged]), "successful repeats"),
    paste("minimum ties per stress repeat:",
          min(rep_dt[scenario == SIM_SCENARIOS[["D4LQ"]] & failure == FALSE, n_ties])),
    paste(sum(rep_dt$failure), "failures with explicit stage attribution"),
    "G4-L is the support-failure positive control",
    "G6-L is the adequate-support negative control"
  )
)
qcs <- rbind(qcs, boundary_qc, fill = TRUE)

forbidden <- list.files(
  SIM_TABLE_DIR, recursive = TRUE, full.names = FALSE,
  pattern = "(hte|interaction|evalue|treatment_effect|risk_difference)",
  ignore.case = TRUE
)
qcs <- rbind(qcs, data.table(
  check = "no_prohibited_effect_output",
  pass = length(forbidden) == 0L,
  detail = if (length(forbidden)) paste(forbidden, collapse = "; ") else "none"
))
sim_atomic_fwrite(qcs, file.path(SIM_LOG_DIR, "D4_sim_independent_QC.csv"))
if (any(!qcs$pass)) {
  stop("Independent D4 simulation QC failed: ",
       paste(qcs[pass == FALSE, check], collapse = ", "))
}
sim_atomic_write_lines(c(
  paste0("run_mode=", SIM_RUN_MODE),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "independent_qc=PASS", "effect_model_fitted=FALSE",
  "causal_effect_estimated=FALSE"
), file.path(SIM_OUTPUT_ROOT, "D4_SIMULATION_COMPLETED.ok"))
