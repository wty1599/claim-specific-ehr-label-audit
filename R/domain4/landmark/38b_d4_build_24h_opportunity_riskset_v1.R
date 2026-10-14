source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

d4_require_stage("38a_source_audit")
stage <- "38b_riskset"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

source_path <- file.path(D4_PRIVATE_DIR, "D4_locked_source_private.rds")
dt <- as.data.table(readRDS(source_path))

d0 <- copy(dt)
d1 <- d0[is.na(first_rrt_hours) | first_rrt_hours >= 0]
d2 <- d1[is.na(first_rrt_hours) | first_rrt_hours >= D4_LANDMARK_H]
d3 <- d2[is.na(death_time_hours) | death_time_hours >= D4_LANDMARK_H]
risk <- d3[outtime_hours >= D4_LANDMARK_H]

flow <- data.table(
  step = c(
    "Locked cohort",
    "Exclude pre-ICU RRT",
    "Exclude RRT from ICU admission to <24h",
    "Exclude death before 24h",
    "Exclude ICU exit before 24h",
    "Final hour-24 opportunity risk set"
  ),
  n_remaining = c(
    nrow(d0), nrow(d1), nrow(d2), nrow(d3), nrow(risk), nrow(risk)
  ),
  n_removed_at_step = c(
    0L, nrow(d0) - nrow(d1), nrow(d1) - nrow(d2),
    nrow(d2) - nrow(d3), nrow(d3) - nrow(risk), 0L
  )
)

candidate_time <- function(x, lower, upper) {
  fifelse(!is.na(x) & x >= lower & x < upper, x, Inf)
}
risk[, rrt_candidate := candidate_time(
  first_rrt_hours, D4_LANDMARK_H, D4_WINDOW_END_H
)]
risk[, death_candidate := candidate_time(
  death_time_hours, D4_LANDMARK_H, D4_WINDOW_END_H
)]
risk[, exit_candidate := candidate_time(
  outtime_hours, D4_LANDMARK_H, D4_WINDOW_END_H
)]
risk[, admin_candidate := D4_WINDOW_END_H]
risk[, first_event_time := pmin(
  rrt_candidate, death_candidate, exit_candidate, admin_candidate
)]
risk[, event_time_tie_n :=
  as.integer(abs(rrt_candidate - first_event_time) <= D4_TIE_TOL) +
  as.integer(abs(death_candidate - first_event_time) <= D4_TIE_TOL) +
  as.integer(abs(exit_candidate - first_event_time) <= D4_TIE_TOL) +
  as.integer(abs(admin_candidate - first_event_time) <= D4_TIE_TOL)]
risk[, event_time_tie := event_time_tie_n > 1L]

# Frozen exact-tie priority: RRT, then death, then ICU exit, then hour 72.
risk[, first_event_state := fcase(
  abs(rrt_candidate - first_event_time) <= D4_TIE_TOL, "RRT_FIRST_24_72",
  abs(death_candidate - first_event_time) <= D4_TIE_TOL, "DEATH_FIRST_24_72",
  abs(exit_candidate - first_event_time) <= D4_TIE_TOL, "ICU_EXIT_FIRST_24_72",
  default = "EVENT_FREE_AT_72"
)]
risk[, rrt_first_24_72 := as.integer(first_event_state == "RRT_FIRST_24_72")]

state_levels <- c(
  "RRT_FIRST_24_72", "DEATH_FIRST_24_72",
  "ICU_EXIT_FIRST_24_72", "EVENT_FREE_AT_72"
)
if (anyNA(risk$first_event_state) ||
    !all(risk$first_event_state %chin% state_levels)) {
  stop("Invalid first-event state.", call. = FALSE)
}
if (sum(table(risk$first_event_state)) != nrow(risk)) {
  stop("First-event states do not reconcile to the risk set.", call. = FALSE)
}
if (any(risk[
  first_event_state == "RRT_FIRST_24_72",
  first_rrt_hours > pmin(death_candidate, exit_candidate, D4_WINDOW_END_H)
])) {
  stop("RRT after a competing/precluding event was coded as treatment.",
       call. = FALSE)
}

event_by_k2 <- risk[, .(n = .N), by = .(cluster_k2, first_event_state)]
event_by_k2 <- merge(
  event_by_k2,
  risk[, .(phenotype_n = .N), by = cluster_k2],
  by = "cluster_k2",
  all.x = TRUE,
  sort = FALSE
)
event_by_k2[, proportion := n / phenotype_n]
event_by_k2[, phenotype := fifelse(
  cluster_k2 == 1L, "C1 higher-risk", "C2 lower-risk"
)]
setcolorder(
  event_by_k2,
  c(
    "cluster_k2", "phenotype", "first_event_state",
    "n", "phenotype_n", "proportion"
  )
)

tie_audit <- risk[event_time_tie == TRUE, .(
  stay_id, cluster_k2, first_event_time, event_time_tie_n,
  rrt_candidate, death_candidate, exit_candidate, admin_candidate,
  assigned_state = first_event_state
)]
tie_summary <- data.table(
  n_risk_set = nrow(risk),
  n_tied = nrow(tie_audit),
  tied_proportion = nrow(tie_audit) / nrow(risk),
  primary_tie_priority = "RRT > death > ICU exit > hour 72"
)

risk_qc <- data.table(
  check = c(
    "risk_set_nonempty",
    "all_outtime_at_or_after_24h",
    "all_alive_at_24h",
    "all_rrt_free_before_24h",
    "four_states_exhaustive",
    "binary_treatment_matches_state",
    "locked_labels_retained"
  ),
  pass = c(
    nrow(risk) > 0L,
    all(risk$outtime_hours >= D4_LANDMARK_H),
    all(is.na(risk$death_time_hours) |
          risk$death_time_hours >= D4_LANDMARK_H),
    all(is.na(risk$first_rrt_hours) |
          risk$first_rrt_hours >= D4_LANDMARK_H),
    sum(event_by_k2$n) == nrow(risk),
    all(risk$rrt_first_24_72 ==
          as.integer(risk$first_event_state == "RRT_FIRST_24_72")),
    all(risk$cluster_k2 %in% c(1L, 2L))
  ),
  observed = c(
    nrow(risk),
    min(risk$outtime_hours),
    sum(!is.na(risk$death_time_hours) &
          risk$death_time_hours < D4_LANDMARK_H),
    sum(!is.na(risk$first_rrt_hours) &
          risk$first_rrt_hours < D4_LANDMARK_H),
    sum(event_by_k2$n),
    sum(risk$rrt_first_24_72),
    paste(risk[, .N, by = cluster_k2][order(cluster_k2), N],
          collapse = "/")
  ),
  expected = c(
    ">0", ">=24", "0", "0", as.character(nrow(risk)),
    as.character(sum(risk$first_event_state == "RRT_FIRST_24_72")),
    "two locked groups"
  )
)
if (any(!risk_qc$pass)) stop("D4 risk-set QC failed.", call. = FALSE)

d4_atomic_fwrite(flow, file.path(D4_TABLE_DIR, "D4_24h_landmark_flow.csv"))
d4_atomic_fwrite(
  event_by_k2,
  file.path(D4_TABLE_DIR, "D4_24_72_first_event_by_K2.csv")
)
d4_atomic_fwrite(
  tie_summary,
  file.path(D4_TABLE_DIR, "D4_event_time_tie_summary.csv")
)
d4_atomic_fwrite(
  risk_qc,
  file.path(D4_LOG_DIR, "D4_24h_landmark_QC.csv")
)
d4_atomic_save_rds(
  risk,
  file.path(D4_PRIVATE_DIR, "D4_landmark_patient_ledger_private.rds")
)
if (nrow(tie_audit)) {
  d4_atomic_fwrite(
    tie_audit,
    file.path(D4_PRIVATE_DIR, "D4_event_time_ties_private.csv")
  )
}
d4_mark_stage(stage)
d4_msg("Completed ", stage, "; n=", nrow(risk))
