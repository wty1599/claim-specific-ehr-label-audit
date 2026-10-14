if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))

sim_candidate_time <- function(x, lower = 24, upper = 72) {
  fifelse(!is.na(x) & x >= lower & x < upper, x, Inf)
}

build_d4_landmark_ledger_sim <- function(dt) {
  d0 <- copy(as.data.table(dt))
  d1 <- d0[is.na(first_rrt_hours) | first_rrt_hours >= 0]
  d2 <- d1[is.na(first_rrt_hours) | first_rrt_hours >= 24]
  d3 <- d2[is.na(death_time_hours) | death_time_hours >= 24]
  risk <- d3[outtime_hours >= 24]

  flow <- data.table(
    step = c(
      "Locked cohort", "Exclude pre-ICU RRT",
      "Exclude RRT from ICU admission to <24h",
      "Exclude death before 24h", "Exclude ICU exit before 24h",
      "Final hour-24 opportunity risk set"
    ),
    n_remaining = c(nrow(d0), nrow(d1), nrow(d2), nrow(d3), nrow(risk), nrow(risk)),
    n_removed_at_step = c(
      0L, nrow(d0) - nrow(d1), nrow(d1) - nrow(d2),
      nrow(d2) - nrow(d3), nrow(d3) - nrow(risk), 0L
    )
  )

  risk[, rrt_candidate := sim_candidate_time(first_rrt_hours)]
  risk[, death_candidate := sim_candidate_time(death_time_hours)]
  risk[, exit_candidate := sim_candidate_time(outtime_hours)]
  risk[, admin_candidate := 72]
  risk[, first_event_time := pmin(
    rrt_candidate, death_candidate, exit_candidate, admin_candidate
  )]
  tol <- 1e-8
  risk[, event_time_tie_n :=
    as.integer(abs(rrt_candidate - first_event_time) <= tol) +
    as.integer(abs(death_candidate - first_event_time) <= tol) +
    as.integer(abs(exit_candidate - first_event_time) <= tol) +
    as.integer(abs(admin_candidate - first_event_time) <= tol)]
  risk[, event_time_tie := event_time_tie_n > 1L]
  risk[, first_event_state := fcase(
    abs(rrt_candidate - first_event_time) <= tol, "RRT_FIRST_24_72",
    abs(death_candidate - first_event_time) <= tol, "DEATH_FIRST_24_72",
    abs(exit_candidate - first_event_time) <= tol, "ICU_EXIT_FIRST_24_72",
    default = "EVENT_FREE_AT_72"
  )]
  risk[, rrt_first_24_72 := as.integer(first_event_state == "RRT_FIRST_24_72")]

  levels_expected <- c(
    "RRT_FIRST_24_72", "DEATH_FIRST_24_72",
    "ICU_EXIT_FIRST_24_72", "EVENT_FREE_AT_72"
  )
  qc <- data.table(
    check = c(
      "all_outtime_at_or_after_24h", "all_alive_at_24h",
      "all_rrt_free_before_24h", "four_states_exhaustive",
      "binary_treatment_matches_state", "locked_labels_retained",
      "tie_priority_rrt"
    ),
    pass = c(
      all(risk$outtime_hours >= 24),
      all(is.na(risk$death_time_hours) | risk$death_time_hours >= 24),
      all(is.na(risk$first_rrt_hours) | risk$first_rrt_hours >= 24),
      !anyNA(risk$first_event_state) &&
        all(risk$first_event_state %in% levels_expected),
      all(risk$rrt_first_24_72 ==
            as.integer(risk$first_event_state == "RRT_FIRST_24_72")),
      all(risk$cluster_k2 %in% c(1L, 2L)),
      all(risk[event_time_tie & is.finite(rrt_candidate) &
                 abs(rrt_candidate - first_event_time) <= tol,
               first_event_state == "RRT_FIRST_24_72"])
    )
  )
  if (any(!qc$pass)) stop("Landmark/ledger QC failed.", call. = FALSE)
  list(risk = risk, flow = flow, qc = qc)
}

run_d4_empirical_ledger_equivalence <- function() {
  formal_root <- file.path(
    SIM_PROJECT_ROOT, "outputs", "domain4_landmark"
  )
  private_source <- file.path(
    formal_root, "private_not_for_release", "D4_locked_source_private.rds"
  )
  private_ledger <- file.path(
    formal_root, "private_not_for_release", "D4_landmark_patient_ledger_private.rds"
  )
  table_path <- function(x) file.path(formal_root, "tables", x)
  required <- c(
    private_source, private_ledger, table_path("D4_24h_landmark_flow.csv"),
    table_path("D4_24_72_first_event_by_K2.csv"),
    table_path("D4_support_metrics_by_K2.csv"),
    table_path("D4_gate_A_to_D_by_K2.csv"),
    table_path("D4_claim_state.csv")
  )
  if (!all(file.exists(required))) {
    return(data.table(
      check = "empirical_authority_files_available", pass = FALSE,
      detail = paste(required[!file.exists(required)], collapse = "; ")
    ))
  }

  compare_tables <- function(a, b, keys, tolerance = 0) {
    a <- copy(as.data.table(a)); b <- copy(as.data.table(b))
    setorderv(a, keys); setorderv(b, keys)
    common <- intersect(names(a), names(b))
    a <- a[, ..common]; b <- b[, ..common]
    isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
  }
  row <- function(check, pass, detail) {
    data.table(check = check, pass = isTRUE(pass), detail = detail)
  }

  x <- build_d4_landmark_ledger_sim(as.data.table(readRDS(private_source)))
  reference_ledger <- as.data.table(readRDS(private_ledger))
  ledger_cols <- c(
    "stay_id", "cluster_k2", "first_event_state", "first_event_time",
    "event_time_tie", "rrt_first_24_72"
  )
  patient_match <- compare_tables(
    x$risk[, ..ledger_cols], reference_ledger[, ..ledger_cols], "stay_id", 1e-12
  )
  flow_match <- compare_tables(
    x$flow[, .(step, n_remaining, n_removed_at_step)],
    fread(table_path("D4_24h_landmark_flow.csv")), "step"
  )

  observed_events <- x$risk[, .(n = .N), by = .(cluster_k2, first_event_state)]
  observed_events[, phenotype := fifelse(
    cluster_k2 == 1L, "C1 higher-risk", "C2 lower-risk"
  )]
  observed_events[, phenotype_n := sum(n), by = cluster_k2]
  observed_events[, proportion := n / phenotype_n]
  setcolorder(observed_events, c(
    "cluster_k2", "phenotype", "first_event_state", "n",
    "phenotype_n", "proportion"
  ))
  event_match <- compare_tables(
    observed_events, fread(table_path("D4_24_72_first_event_by_K2.csv")),
    c("cluster_k2", "first_event_state"), 1e-12
  )

  fit <- d4_fit_support(x$risk, imputation_mode = "median_mode")
  state <- d4_apply_state(fit$support, fit$cells)
  support_match <- compare_tables(
    fit$support, fread(table_path("D4_support_metrics_by_K2.csv")),
    "phenotype", 1e-12
  )
  gate_match <- compare_tables(
    state$gate_table, fread(table_path("D4_gate_A_to_D_by_K2.csv")),
    c("gate", "phenotype"), 1e-12
  )
  claim_ref <- fread(table_path("D4_claim_state.csv"))
  claim_match <- identical(state$state, claim_ref$state[1]) &&
    identical(state$triggers, claim_ref$triggers[1])

  rbindlist(list(
    row("empirical_authority_files_available", TRUE, "all authority files found"),
    row("empirical_flow_exact", flow_match, "locked landmark flow"),
    row("empirical_patient_ledger_exact", patient_match,
        "stay_id, K2, first state/time, tie, treatment"),
    row("empirical_event_counts_exact", event_match, "K2 first-event table"),
    row("empirical_support_metrics_exact", support_match,
        "portable engine vs locked support metrics"),
    row("empirical_gate_table_exact", gate_match,
        "portable engine vs locked Gate A-D table"),
    row("empirical_claim_state_exact", claim_match,
        "portable engine state and trigger string")
  ))
}
