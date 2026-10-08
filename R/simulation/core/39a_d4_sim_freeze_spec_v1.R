source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))

sim_require_files(c(
  SIM_SPEC, SIM_EMPIRICAL_COMMON, SIM_SUPPORT_ENGINE, SIM_EXPECTED_HASHES
))
sim_init_dirs()
sim_source_empirical_engine()

expected_gates <- data.table(
  gate = c("A", "B", "C", "D"),
  treated_fraction = c(0.005, 0.010, 0.020, 0.010),
  event_threshold = c(5L, 10L, 10L, 20L)
)
observed_gates <- rbindlist(lapply(names(D4_GATES), function(g) {
  data.table(
    gate = g, treated_fraction = D4_GATES[[g]]$frac,
    event_threshold = D4_GATES[[g]]$events
  )
}))
if (!identical(expected_gates, observed_gates)) {
  stop("Empirical D4 gates differ from the frozen simulation specification.")
}
expected_support <- c(
  outside_overlap_max = 0.10, ess_ratio_min = 0.25,
  minimum_cell_min = 30, extreme_weight_cut = 10,
  extreme_weight_prop_max = 0.01
)
observed_support <- unlist(D4_SUPPORT_THRESHOLDS, use.names = TRUE)
if (!isTRUE(all.equal(expected_support, observed_support, tolerance = 0))) {
  stop("Empirical D4 support thresholds differ from the frozen specification.")
}

config <- data.table(
  parameter = c(
    "run_mode", "base_seed", "cohort_n", "g4l_repetitions",
    "g6l_repetitions", "ledger_qc_repetitions", "landmark_hour",
    "window_end_hour", "primary_gate", names(expected_support)
  ),
  value = as.character(c(
    SIM_RUN_MODE, SIM_BASE_SEED, SIM_N, SIM_REPS[["G4L"]],
    SIM_REPS[["G6L"]], SIM_REPS[["D4LQ"]], 24, 72, "B",
    unname(expected_support)
  ))
)
sim_atomic_fwrite(config, file.path(SIM_PROVENANCE_DIR, "D4_sim_frozen_config.csv"))
sim_atomic_fwrite(observed_gates, file.path(
  SIM_PROVENANCE_DIR, "D4_sim_locked_gate_table.csv"
))

input_files <- c(
  frozen_spec = SIM_SPEC,
  empirical_engine_authority = SIM_EMPIRICAL_COMMON,
  portable_support_engine = SIM_SUPPORT_ENGINE,
  simulation_common = file.path(SIM_SCRIPT_DIR, "39_d4_sim_common_v1.R"),
  freeze = file.path(SIM_SCRIPT_DIR, "39a_d4_sim_freeze_spec_v1.R"),
  generator = file.path(SIM_SCRIPT_DIR, "39b_d4_sim_generate_longitudinal_controls_v1.R"),
  ledger = file.path(SIM_SCRIPT_DIR, "39c_d4_sim_build_landmark_ledger_v1.R"),
  engine = file.path(SIM_SCRIPT_DIR, "39d_d4_sim_apply_support_engine_v1.R"),
  operating_characteristics = file.path(SIM_SCRIPT_DIR, "39e_d4_sim_operating_characteristics_v1.R"),
  figures_tables = file.path(SIM_SCRIPT_DIR, "39f_d4_sim_figures_tables_v1.R"),
  independent_qc = file.path(SIM_SCRIPT_DIR, "39g_d4_sim_independent_qc_v1.R"),
  runner = file.path(SIM_SCRIPT_DIR, "39_run_d4_longitudinal_simulation_v1.R")
)
sim_require_files(input_files)
hashes <- data.table(
  artifact = names(input_files), path = unname(input_files),
  sha256 = vapply(input_files, sim_sha256, character(1))
)
expected_hashes <- fread(SIM_EXPECTED_HASHES)
d4_require_columns(expected_hashes, c("artifact", "sha256"),
                   "Frozen expected source hashes")
hash_check <- merge(
  hashes[, .(artifact, observed_sha256 = sha256)],
  expected_hashes[, .(artifact, expected_sha256 = sha256)],
  by = "artifact", all = TRUE
)
hash_check[, pass := !is.na(observed_sha256) & !is.na(expected_sha256) &
             observed_sha256 == expected_sha256]
if (any(!hash_check$pass)) {
  stop("Frozen source-hash qualification failed: ",
       paste(hash_check[pass == FALSE, artifact], collapse = ", "))
}
sim_atomic_fwrite(hashes, file.path(SIM_PROVENANCE_DIR, "D4_sim_input_hashes.csv"))
sim_atomic_fwrite(hash_check, file.path(
  SIM_PROVENANCE_DIR, "D4_sim_source_hash_qualification.csv"
))
sim_msg("Frozen specification and empirical rule equivalence passed.")
