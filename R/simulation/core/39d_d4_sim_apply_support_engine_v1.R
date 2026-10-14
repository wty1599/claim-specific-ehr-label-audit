if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))
source(file.path(SIM_SCRIPT_DIR, "39b_d4_sim_generate_longitudinal_controls_v1.R"))
source(file.path(SIM_SCRIPT_DIR, "39c_d4_sim_build_landmark_ledger_v1.R"))
if (!exists("d4_fit_support")) sim_source_empirical_engine()

equiv <- run_d4_empirical_ledger_equivalence()
sim_atomic_fwrite(equiv, file.path(SIM_LOG_DIR, "D4_empirical_ledger_equivalence.csv"))
if (anyNA(equiv$pass) || any(!equiv$pass)) {
  skip_private_qc <- identical(
    Sys.getenv("D4_SKIP_PRIVATE_EQUIVALENCE_QC", unset = "1"), "1"
  )
  if (!skip_private_qc) {
    stop("Empirical D4 full-chain equivalence failed: ",
         paste(equiv[pass != TRUE | is.na(pass), check], collapse = ", "))
  }
  sim_msg(
    "Restricted empirical ledger files are unavailable; the public run records ",
    "but does not recompute that one equivalence QC. Set ",
    "D4_SKIP_PRIVATE_EQUIVALENCE_QC=0 when restricted inputs are mounted."
  )
}

run_one_d4_sim_repeat <- function(scenario, repeat_id) {
  stage <- "generation"
  tryCatch({
    raw <- simulate_d4_longitudinal_cohort(scenario, repeat_id, SIM_N)
    stage <- "landmark_ledger"
    led <- build_d4_landmark_ledger_sim(raw)
    risk <- led$risk
    if (!nrow(risk)) stop("Empty landmark risk set.")

    stage <- "propensity_and_imputation"
    fit <- d4_fit_support(risk, imputation_mode = "median_mode")
    stage <- "support_metrics_and_state"
    state <- d4_apply_state(fit$support, fit$cells)
    expected <- unname(SIM_EXPECTED_STATE[[scenario]])

    repeat_row <- data.table(
      scenario, repeat_id, seed = sim_seed(scenario, repeat_id),
      n_generated = nrow(raw), n_risk_set = nrow(risk),
      n_c1 = sum(risk$cluster_k2 == 1L),
      n_c2 = sum(risk$cluster_k2 == 2L),
      n_rrt_first = sum(risk$rrt_first_24_72),
      n_ties = sum(risk$event_time_tie),
      ledger_qc_pass = all(led$qc$pass),
      ps_model_converged = isTRUE(fit$model_qc$converged),
      expected_state = expected,
      observed_state = state$state,
      correct_classification = if (is.na(expected)) NA else state$state == expected,
      triggers = state$triggers,
      failure = FALSE,
      failure_stage = NA_character_, failure_reason = NA_character_
    )
    support <- copy(fit$support)[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    gates <- copy(state$gate_table)[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    alerts <- copy(state$support_alerts)[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    flow <- copy(led$flow)[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    events <- risk[, .(n = .N), by = .(cluster_k2, first_event_state)]
    events[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    true_prob <- risk[, .(
      true_rrt_probability_mean = mean(true_rrt_probability),
      true_rrt_probability_p025 = quantile(true_rrt_probability, 0.025),
      true_rrt_probability_p975 = quantile(true_rrt_probability, 0.975)
    ), by = cluster_k2]
    true_prob[, `:=`(scenario = scenario, repeat_id = repeat_id)]
    list(repeat_result = repeat_row, support = support, gates = gates, alerts = alerts,
         flow = flow, events = events, true_prob = true_prob)
  }, error = function(e) {
    list(
      repeat_result = data.table(
        scenario, repeat_id, seed = sim_seed(scenario, repeat_id),
        n_generated = SIM_N, n_risk_set = NA_integer_, n_c1 = NA_integer_,
        n_c2 = NA_integer_, n_rrt_first = NA_integer_, n_ties = NA_integer_,
        ledger_qc_pass = NA, ps_model_converged = NA,
        expected_state = unname(SIM_EXPECTED_STATE[[scenario]]),
        observed_state = NA_character_, correct_classification = NA,
        triggers = NA_character_, failure = TRUE,
        failure_stage = stage,
        failure_reason = conditionMessage(e)
      ),
      support = data.table(), gates = data.table(), alerts = data.table(),
      flow = data.table(), events = data.table(), true_prob = data.table()
    )
  })
}

all_out <- list()
for (key in names(SIM_SCENARIOS)) {
  scenario <- SIM_SCENARIOS[[key]]
  n_rep <- SIM_REPS[[key]]
  sim_msg("Running ", scenario, " (", n_rep, " repetitions)")
  z <- vector("list", n_rep)
  for (i in seq_len(n_rep)) {
    z[[i]] <- run_one_d4_sim_repeat(scenario, i)
    if (i %% max(1L, floor(n_rep / 10L)) == 0L || i == n_rep) {
      sim_msg("  ", scenario, ": ", i, "/", n_rep)
    }
  }
  checkpoint <- file.path(SIM_CHECKPOINT_DIR, paste0(key, "_completed.rds"))
  saveRDS(z, checkpoint, compress = "xz")
  all_out[[key]] <- z
}

bind_component <- function(name) {
  rbindlist(lapply(all_out, function(sc) {
    rbindlist(lapply(sc, `[[`, name), fill = TRUE)
  }), fill = TRUE)
}
sim_atomic_fwrite(bind_component("repeat_result"), file.path(
  SIM_TABLE_DIR, "D4_sim_repeat_level.csv"
))
sim_atomic_fwrite(bind_component("support"), file.path(
  SIM_TABLE_DIR, "D4_sim_support_metrics_by_K2.csv"
))
sim_atomic_fwrite(bind_component("gates"), file.path(
  SIM_TABLE_DIR, "D4_sim_gate_A_to_D_by_K2.csv"
))
sim_atomic_fwrite(bind_component("alerts"), file.path(
  SIM_TABLE_DIR, "D4_sim_structural_support_alerts_by_K2.csv"
))
sim_atomic_fwrite(bind_component("flow"), file.path(
  SIM_TABLE_DIR, "D4_sim_landmark_flow_by_repeat.csv"
))
sim_atomic_fwrite(bind_component("events"), file.path(
  SIM_TABLE_DIR, "D4_sim_first_event_counts_by_repeat.csv"
))
sim_atomic_fwrite(bind_component("true_prob"), file.path(
  SIM_TABLE_DIR, "D4_sim_true_treatment_probability_by_K2.csv"
))
sim_msg("Simulation engine stage completed.")
