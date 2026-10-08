if (!exists("D3R_ROOT")) source(file.path(D3R_SCRIPT_DIR, "41_d3r_sim_common_v1.R"))
d3r_init_dirs()
if (!file.exists(D3R_SPEC)) stop("Frozen specification is missing: ", D3R_SPEC)

authority_main <- file.path(D3R_SCRIPT_DIR, "00_main_S1_S5_formal.R")

files <- c(
  frozen_spec = D3R_SPEC,
  legacy_simulation_authority = authority_main,
  empirical_generator_config = file.path(D3R_PROJECT_ROOT, "R", "domain3", "generator", "00_d3_generator_config_v1.R"),
  empirical_generator_utils = file.path(D3R_PROJECT_ROOT, "R", "domain3", "generator", "00_d3_generator_utils_v1.R"),
  common = file.path(D3R_SCRIPT_DIR, "41_d3r_sim_common_v1.R"),
  freeze = file.path(D3R_SCRIPT_DIR, "41a_d3r_sim_freeze_spec_v1.R"),
  execute = file.path(D3R_SCRIPT_DIR, "41b_d3r_sim_execute_v1.R"),
  postprocess = file.path(D3R_SCRIPT_DIR, "41c_d3r_sim_postprocess_v1.R"),
  qc = file.path(D3R_SCRIPT_DIR, "41d_d3r_sim_independent_qc_v1.R"),
  runner = file.path(D3R_SCRIPT_DIR, "41_run_d3_regenerated_generator_simulation_v1.R")
)
if (any(!file.exists(files))) stop("Missing source file(s): ", paste(files[!file.exists(files)], collapse = "; "))
observed <- data.table(artifact = names(files), path = unname(files),
                       sha256 = vapply(files, d3r_sha256, character(1)))
d3r_atomic_fwrite(observed, file.path(D3R_PROVENANCE_DIR, "D3_sim_input_hashes.csv"), overwrite = TRUE)

if (D3R_RUN_MODE == "formal") {
  if (!file.exists(D3R_EXPECTED_HASHES)) stop("Formal run requires frozen expected hashes: ", D3R_EXPECTED_HASHES)
  expected <- fread(D3R_EXPECTED_HASHES)
  check <- merge(
    observed[, .(artifact, observed_sha256 = sha256)],
    expected[, .(artifact, expected_sha256 = sha256)], by = "artifact", all = TRUE
  )
  check[, pass := !is.na(observed_sha256) & !is.na(expected_sha256) & observed_sha256 == expected_sha256]
  d3r_atomic_fwrite(check, file.path(D3R_PROVENANCE_DIR, "D3_sim_source_hash_qualification.csv"), overwrite = TRUE)
  if (any(!check$pass)) stop("Formal source-hash qualification failed: ", paste(check[pass == FALSE, artifact], collapse = ", "))
}

config <- data.table(
  parameter = c("run_mode", "base_seed", "n_train", "n_external", "repetitions_per_scenario",
                "mice_m", "mice_maxit", "primary_imputation", "cv_folds", "workers",
                "prevalence_gate", "assignment_ari_gate", "calibration_slope_low",
                "calibration_slope_high", "calibration_intercept_gate"),
  value = as.character(c(D3R_RUN_MODE, D3R_BASE_SEED, D3R_N_TRAIN, D3R_N_EXTERNAL,
                         D3R_N_REP, D3R_MICE_M, D3R_MICE_MAXIT,
                         D3R_PRIMARY_IMPUTATION, D3R_CV_FOLDS, D3R_WORKERS,
                         D3R_PREVALENCE_GATE, D3R_ASSIGNMENT_ARI_GATE,
                         D3R_CALIBRATION_SLOPE_RANGE, D3R_CALIBRATION_INTERCEPT_GATE))
)
d3r_atomic_fwrite(config, file.path(D3R_PROVENANCE_DIR, "D3_sim_frozen_config.csv"), overwrite = TRUE)
d3r_atomic_fwrite(data.table(scenario = names(D3R_EXPECTED_ALERT), expected_alert = unname(D3R_EXPECTED_ALERT)),
                  file.path(D3R_PROVENANCE_DIR, "D3_sim_truth_map.csv"), overwrite = TRUE)
d3r_msg("Specification and provenance lock completed.")
