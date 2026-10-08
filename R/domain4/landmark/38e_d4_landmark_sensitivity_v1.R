source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

d4_require_stage("38d_locked_gates")
stage <- "38e_sensitivity"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

risk <- as.data.table(readRDS(file.path(
  D4_PRIVATE_DIR, "D4_landmark_patient_ledger_private.rds"
)))

variants <- list(
  primary_recomputed = list(
    data = copy(risk),
    imputation = "median_mode",
    definition = "outtime>=24; frozen tie priority; median/mode"
  ),
  strict_outtime_gt24 = list(
    data = risk[outtime_hours > D4_LANDMARK_H],
    imputation = "median_mode",
    definition = "outtime>24; frozen tie priority; median/mode"
  ),
  exclude_exact_event_ties = list(
    data = risk[event_time_tie == FALSE],
    imputation = "median_mode",
    definition = "outtime>=24; exact event ties excluded; median/mode"
  ),
  complete_case_covariates = list(
    data = copy(risk),
    imputation = "complete_case",
    definition = "outtime>=24; frozen tie priority; complete covariates"
  )
)

variant_rows <- list()
for (nm in names(variants)) {
  v <- variants[[nm]]
  fit <- tryCatch(
    d4_fit_support(v$data, imputation_mode = v$imputation),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    variant_rows[[nm]] <- data.table(
      variant = nm,
      definition = v$definition,
      status = "MODEL_FAILURE",
      error = conditionMessage(fit),
      n = nrow(v$data),
      treated_n = sum(v$data$rrt_first_24_72),
      state = "INTERFACE_INCOMPLETE",
      triggers = conditionMessage(fit)
    )
    next
  }
  state <- d4_apply_state(fit$support, fit$cells)
  variant_rows[[nm]] <- data.table(
    variant = nm,
    definition = v$definition,
    status = "COMPLETED",
    error = "",
    n = nrow(fit$data),
    treated_n = sum(fit$data$rrt_first_24_72),
    c1_n = fit$support[phenotype == "C1 higher-risk", n],
    c1_treated_n = fit$support[phenotype == "C1 higher-risk", n_treated],
    c2_n = fit$support[phenotype == "C2 lower-risk", n],
    c2_treated_n = fit$support[phenotype == "C2 lower-risk", n_treated],
    c1_outside_overlap =
      fit$support[phenotype == "C1 higher-risk", outside_overlap_proportion],
    c2_outside_overlap =
      fit$support[phenotype == "C2 lower-risk", outside_overlap_proportion],
    c1_ess_ratio = fit$support[phenotype == "C1 higher-risk", ess_ratio],
    c2_ess_ratio = fit$support[phenotype == "C2 lower-risk", ess_ratio],
    state = state$state,
    triggers = state$triggers
  )
}
sensitivity <- rbindlist(variant_rows, fill = TRUE, use.names = TRUE)

primary_state <- sensitivity[
  variant == "primary_recomputed" & status == "COMPLETED", state
]
if (length(primary_state) != 1L) {
  stop("Primary sensitivity recomputation failed.", call. = FALSE)
}
sensitivity[, state_matches_primary := state == primary_state]

d4_atomic_fwrite(
  sensitivity,
  file.path(
    D4_TABLE_DIR, "D4_boundary_tie_and_complete_case_sensitivity.csv"
  )
)
d4_mark_stage(stage)
d4_msg("Completed ", stage)

