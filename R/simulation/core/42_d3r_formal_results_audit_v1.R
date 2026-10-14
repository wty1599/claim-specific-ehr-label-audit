options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
})

root <- normalizePath(
  file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "outputs", "simulation", "D3_G4abc_formal"),
  winslash = "/",
  mustWork = FALSE
)
formal_dir <- root
table_dir <- file.path(formal_dir, "tables")
audit_dir <- file.path(formal_dir, "audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

repeat_file <- file.path(table_dir, "D3_sim_repeat_level.csv")
if (!file.exists(repeat_file)) stop("Formal repeat-level file is missing.")

dt <- fread(repeat_file)

truth <- data.table(
  scenario = c(
    "G3R_same_distribution_no_alert",
    "G3R_covariate_drift_alert",
    "D3CQ_calibration_stress_control"
  ),
  expected_alert_audit = c(FALSE, TRUE, TRUE)
)

prev_gate <- 0.10
ari_gate <- 0.80
slope_low <- 0.80
slope_high <- 1.20
intercept_gate <- 0.20

wilson <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, center - half), high = min(1, center + half))
}

if (nrow(dt) != 300L) stop("Expected 300 formal repeat rows; found ", nrow(dt), ".")
if (uniqueN(dt, by = c("scenario", "repeat_id")) != 300L) {
  stop("Duplicate scenario-repeat keys detected.")
}

dt <- merge(dt, truth, by = "scenario", all.x = TRUE, sort = FALSE)
if (anyNA(dt$expected_alert_audit)) stop("Unknown scenario in repeat-level file.")

dt[, prevalence_alert_audit := abs(prevalence_drift) > prev_gate]
dt[, assignment_alert_audit := external_reference_ari < ari_gate]
dt[, calibration_alert_audit :=
     calibration_slope < slope_low |
     calibration_slope > slope_high |
     abs(calibration_intercept) > intercept_gate]
dt[, combined_alert_audit :=
     prevalence_alert_audit | assignment_alert_audit | calibration_alert_audit]
dt[, correct_audit := combined_alert_audit == expected_alert_audit]

flag_match <-
  all(dt$prevalence_alert == dt$prevalence_alert_audit) &&
  all(dt$assignment_alert == dt$assignment_alert_audit) &&
  all(dt$calibration_alert == dt$calibration_alert_audit) &&
  all(dt$combined_alert == dt$combined_alert_audit) &&
  all(dt$correct_classification == dt$correct_audit)

scenario_summary <- dt[, {
  n <- .N
  correct <- sum(correct_audit)
  ci <- wilson(correct, n)
  .(
    expected_alert = unique(expected_alert_audit),
    n_rep = n,
    failed_n = sum(failure),
    correct_n = correct,
    correct_rate = correct / n,
    wilson_low = ci[["low"]],
    wilson_high = ci[["high"]],
    prevalence_alert_n = sum(prevalence_alert_audit),
    assignment_alert_n = sum(assignment_alert_audit),
    calibration_alert_n = sum(calibration_alert_audit),
    combined_alert_n = sum(combined_alert_audit)
  )
}, by = scenario]

direction_summary <- dt[, {
  n <- .N
  correct <- sum(correct_audit)
  ci <- wilson(correct, n)
  .(
    n_rep = n,
    correct_n = correct,
    rate = correct / n,
    wilson_low = ci[["low"]],
    wilson_high = ci[["high"]]
  )
}, by = .(direction = ifelse(expected_alert_audit, "sensitivity", "specificity"))]

overall_ci <- wilson(sum(dt$correct_audit), nrow(dt))
overall <- data.table(
  direction = "overall_correct_classification",
  n_rep = nrow(dt),
  correct_n = sum(dt$correct_audit),
  rate = mean(dt$correct_audit),
  wilson_low = overall_ci[["low"]],
  wilson_high = overall_ci[["high"]]
)
direction_summary <- rbind(direction_summary, overall)

key_metrics <- c(
  "internal_agreement", "internal_ari", "crossfit_agreement", "crossfit_ari",
  "mice_ari_mean", "mice_ari_min", "external_reference_agreement",
  "external_reference_ari", "derivation_c1_prevalence",
  "external_c1_prevalence", "prevalence_drift", "calibration_slope",
  "calibration_intercept", "external_auc", "external_brier",
  "external_missing_feature_mean", "external_clipped_feature_mean",
  "external_margin_median"
)

metric_summary <- rbindlist(lapply(key_metrics, function(metric_name) {
  dt[, {
    x <- get(metric_name)
    .(
      metric = metric_name,
      n_rep = sum(is.finite(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      mcse = sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))),
      empirical_p025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
      empirical_p975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
      minimum = min(x, na.rm = TRUE),
      maximum = max(x, na.rm = TRUE)
    )
  }, by = scenario]
}))

source_hash_qc <- fread(file.path(
  formal_dir, "provenance", "D3_sim_source_hash_qualification.csv"
))

qc <- data.table(
  check = c(
    "formal_completion_marker",
    "three_scenarios_100_repeats_each",
    "no_failed_repeats",
    "all_assignments_complete",
    "all_batch_invariance_pass",
    "source_hashes_match",
    "independent_flags_match_saved_flags",
    "all_expected_classifications_reproduced"
  ),
  pass = c(
    file.exists(file.path(formal_dir, "D3_SIMULATION_COMPLETED.ok")),
    nrow(scenario_summary) == 3L && all(scenario_summary$n_rep == 100L),
    !any(dt$failure),
    all(dt$external_assignment_complete),
    all(dt$batch_invariance_pass),
    all(source_hash_qc$pass),
    flag_match,
    all(dt$correct_audit)
  )
)

fwrite(scenario_summary, file.path(audit_dir, "D3_formal_scenario_audit.csv"))
fwrite(direction_summary, file.path(audit_dir, "D3_formal_operating_characteristics_audit.csv"))
fwrite(metric_summary, file.path(audit_dir, "D3_formal_metric_summary_audit.csv"))
fwrite(qc, file.path(audit_dir, "D3_formal_second_pass_QC.csv"))

report <- c(
  "# D3 regenerated-generator formal simulation: independent result audit",
  "",
  "This audit reads the locked repeat-level output and independently reconstructs",
  "the pre-specified D3 component alerts. It does not overwrite source results or",
  "modify the manuscript.",
  "",
  "## Scenario-level result",
  "",
  paste(capture.output(print(scenario_summary)), collapse = "\n"),
  "",
  "## Direction-level operating characteristics",
  "",
  paste(capture.output(print(direction_summary)), collapse = "\n"),
  "",
  "## Independent QC",
  "",
  paste(capture.output(print(qc)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "The simulation-only complete-feature external reference quantifies generator",
  "fidelity under a known DGP. It does not create an empirical eICU reference label",
  "and therefore does not remove the empirical REFERENCE_ABSENT boundary."
)
writeLines(report, file.path(audit_dir, "D3_FORMAL_INDEPENDENT_AUDIT.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo_audit.txt"), useBytes = TRUE)

if (any(!qc$pass)) {
  stop("Independent formal audit failed: ", paste(qc[pass == FALSE, check], collapse = ", "))
}

cat("Independent D3 formal audit: PASS\n")
