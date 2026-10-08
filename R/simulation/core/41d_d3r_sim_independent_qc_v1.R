if (!exists("D3R_ROOT")) source(file.path(D3R_SCRIPT_DIR, "41_d3r_sim_common_v1.R"))
rep_dt <- fread(file.path(D3R_TABLE_DIR, "D3_sim_repeat_level.csv"))
truth <- fread(file.path(D3R_PROVENANCE_DIR, "D3_sim_truth_map.csv"))
expected_rows <- length(D3R_SCENARIOS) * D3R_N_REP

recalc <- copy(rep_dt[failure == FALSE])
recalc[, prevalence_alert_recalc := abs(prevalence_drift) > D3R_PREVALENCE_GATE]
recalc[, assignment_alert_recalc := external_reference_ari < D3R_ASSIGNMENT_ARI_GATE]
recalc[, calibration_alert_recalc := abs(calibration_intercept) > D3R_CALIBRATION_INTERCEPT_GATE |
         calibration_slope < D3R_CALIBRATION_SLOPE_RANGE[1] |
         calibration_slope > D3R_CALIBRATION_SLOPE_RANGE[2]]
recalc[, combined_alert_recalc := prevalence_alert_recalc | assignment_alert_recalc | calibration_alert_recalc]

qc <- data.table(
  check = c(
    "planned_repeat_count", "truth_map_exact", "threshold_flags_reproduce",
    "technical_assignment_complete", "batch_invariance", "mice_configuration",
    "formal_source_hash_qualification", "no_manuscript_output"
  ),
  pass = c(
    nrow(rep_dt) == expected_rows && uniqueN(rep_dt[, .(scenario, repeat_id)]) == expected_rows,
    identical(sort(truth$scenario), sort(names(D3R_EXPECTED_ALERT))) &&
      all(truth$expected_alert[match(names(D3R_EXPECTED_ALERT), truth$scenario)] == unname(D3R_EXPECTED_ALERT)),
    all(recalc$prevalence_alert == recalc$prevalence_alert_recalc) &&
      all(recalc$assignment_alert == recalc$assignment_alert_recalc) &&
      all(recalc$calibration_alert == recalc$calibration_alert_recalc) &&
      all(recalc$combined_alert == recalc$combined_alert_recalc),
    all(recalc$external_assignment_complete),
    all(recalc$batch_invariance_pass),
    D3R_MICE_M == 5L && D3R_MICE_MAXIT == 5L && D3R_PRIMARY_IMPUTATION == 1L,
    if (D3R_RUN_MODE == "formal") {
      q <- fread(file.path(D3R_PROVENANCE_DIR, "D3_sim_source_hash_qualification.csv")); all(q$pass)
    } else TRUE,
    !any(grepl("manuscript|supplement", list.files(D3R_OUTPUT_ROOT, recursive = TRUE), ignore.case = TRUE))
  ),
  detail = c(
    paste(nrow(rep_dt), "of", expected_rows),
    paste(truth$scenario, truth$expected_alert, collapse = "; "),
    paste(nrow(recalc), "successful repeats independently reconstructed"),
    paste(sum(recalc$external_assignment_complete), "of", nrow(recalc)),
    paste(sum(recalc$batch_invariance_pass), "of", nrow(recalc)),
    "m=5; maxit=5; primary=1",
    if (D3R_RUN_MODE == "formal") "formal manifest checked" else "not required for smoke",
    "output directory contains simulation artifacts only"
  )
)
d3r_atomic_fwrite(qc, file.path(D3R_LOG_DIR, "D3_sim_independent_QC.csv"), overwrite = TRUE)
if (any(!qc$pass)) stop("Independent D3 simulation QC failed: ", paste(qc[pass == FALSE, check], collapse = ", "))
d3r_atomic_write_lines(c(
  paste0("run_mode=", D3R_RUN_MODE),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "independent_qc=PASS", "manuscript_modified=FALSE"
), file.path(D3R_OUTPUT_ROOT, "D3_SIMULATION_COMPLETED.ok"), overwrite = TRUE)
