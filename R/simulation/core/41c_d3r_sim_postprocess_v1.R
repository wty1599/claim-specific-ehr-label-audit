if (!exists("D3R_ROOT")) source(file.path(D3R_SCRIPT_DIR, "41_d3r_sim_common_v1.R"))
rep_dt <- fread(file.path(D3R_TABLE_DIR, "D3_sim_repeat_level.csv"))
ok <- rep_dt[failure == FALSE]

metrics <- c(
  "internal_agreement", "internal_ari", "internal_prevalence_difference",
  "crossfit_agreement", "crossfit_ari", "mice_ari_mean", "mice_ari_min",
  "external_reference_agreement", "external_reference_ari",
  "derivation_c1_prevalence", "external_c1_prevalence", "prevalence_drift",
  "calibration_slope", "calibration_intercept", "external_auc", "external_brier",
  "external_missing_feature_mean", "external_clipped_feature_mean",
  "external_margin_median"
)
continuous <- rbindlist(lapply(metrics, function(m) {
  ok[, {
    x <- get(m); x <- x[is.finite(x)]; n <- length(x)
    .(metric = m, n_rep = n, mean = if (n) mean(x) else NA_real_,
      sd = if (n > 1) sd(x) else NA_real_, mcse = if (n > 1) sd(x) / sqrt(n) else NA_real_,
      empirical_p025 = if (n) quantile(x, 0.025, names = FALSE) else NA_real_,
      empirical_p975 = if (n) quantile(x, 0.975, names = FALSE) else NA_real_)
  }, by = scenario]
}))
setcolorder(continuous, c("scenario", "metric", setdiff(names(continuous), c("scenario", "metric"))))
d3r_atomic_fwrite(continuous, file.path(D3R_TABLE_DIR, "D3_sim_continuous_metrics_mcse.csv"), overwrite = TRUE)

scenario_oc <- rep_dt[, {
  valid <- !failure
  n <- sum(valid); correct <- sum(correct_classification[valid], na.rm = TRUE)
  ci <- d3r_wilson(correct, n)
  n_all <- .N; correct_all <- sum(ifelse(failure, FALSE, correct_classification), na.rm = TRUE)
  ci_all <- d3r_wilson(correct_all, n_all)
  .(planned_n = .N, failed_n = sum(failure), analyzed_n = n, correct_n = correct,
    rate = correct / n, wilson_low = ci[1], wilson_high = ci[2],
    correct_failure_as_error_n = correct_all, failure_as_error_rate = correct_all / n_all,
    failure_as_error_wilson_low = ci_all[1], failure_as_error_wilson_high = ci_all[2])
}, by = .(scenario, expected_alert)]

direction_oc <- rep_dt[, {
  valid <- !failure
  n <- sum(valid); correct <- sum(correct_classification[valid], na.rm = TRUE)
  ci <- d3r_wilson(correct, n)
  .(n_rep = n, correct_n = correct, rate = correct / n,
    wilson_low = ci[1], wilson_high = ci[2], failed_n = sum(failure))
}, by = .(truth_direction = ifelse(expected_alert, "positive_sensitivity", "negative_specificity"))]

overall <- rep_dt[, {
  valid <- !failure; n <- sum(valid); correct <- sum(correct_classification[valid], na.rm = TRUE)
  ci <- d3r_wilson(correct, n)
  .(truth_direction = "overall_correct_classification", n_rep = n,
    correct_n = correct, rate = correct / n, wilson_low = ci[1], wilson_high = ci[2],
    failed_n = sum(failure))
}]
oc <- rbind(direction_oc, overall, fill = TRUE)
d3r_atomic_fwrite(scenario_oc, file.path(D3R_TABLE_DIR, "D3_sim_operating_characteristics_by_scenario.csv"), overwrite = TRUE)
d3r_atomic_fwrite(oc, file.path(D3R_TABLE_DIR, "D3_sim_operating_characteristics.csv"), overwrite = TRUE)

components <- ok[, .(
  n_rep = .N,
  prevalence_alert_rate = mean(prevalence_alert),
  assignment_alert_rate = mean(assignment_alert),
  calibration_alert_rate = mean(calibration_alert),
  combined_alert_rate = mean(combined_alert)
), by = scenario]
d3r_atomic_fwrite(components, file.path(D3R_TABLE_DIR, "D3_sim_component_alert_rates.csv"), overwrite = TRUE)

fail <- rep_dt[, .(planned_n = .N, failed_n = sum(failure), failure_rate = mean(failure)), by = scenario]
fail_reason <- rep_dt[failure == TRUE, .N, by = .(scenario, failure_stage, failure_reason)]
d3r_atomic_fwrite(fail, file.path(D3R_TABLE_DIR, "D3_sim_failure_audit.csv"), overwrite = TRUE)
d3r_atomic_fwrite(fail_reason, file.path(D3R_TABLE_DIR, "D3_sim_failure_reasons.csv"), overwrite = TRUE)

suppressPackageStartupMessages(library(ggplot2))
plot_dt <- melt(ok, id.vars = c("scenario", "repeat_id"),
                measure.vars = c("internal_ari", "external_reference_ari", "calibration_slope"),
                variable.name = "metric", value.name = "value")
plot_dt[, metric := factor(metric, levels = c("internal_ari", "external_reference_ari", "calibration_slope"),
                           labels = c("Internal fidelity ARI", "External reference ARI", "Calibration slope"))]
p <- ggplot(plot_dt, aes(x = scenario, y = value)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35, fill = "grey90") +
  geom_jitter(width = 0.10, alpha = 0.30, size = 0.8, color = "#333333") +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(D3R_FIGURE_DIR, "D3_sim_diagnostic_summary.pdf"), p, width = 7.2, height = 8.2, device = cairo_pdf)
ggsave(file.path(D3R_FIGURE_DIR, "D3_sim_diagnostic_summary.png"), p, width = 7.2, height = 8.2, dpi = 300)

fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "NA")
report <- c(
  "# Regenerated-generator D3 simulation results",
  "",
  paste0("Run mode: ", D3R_RUN_MODE, "; derivation n=", D3R_N_TRAIN,
         "; external n=", D3R_N_EXTERNAL, "; repetitions per scenario=", D3R_N_REP, "."),
  "",
  "## Operating characteristics",
  "",
  paste(capture.output(print(oc)), collapse = "\n"),
  "",
  "## Scenario-specific classification",
  "",
  paste(capture.output(print(scenario_oc)), collapse = "\n"),
  "",
  "## Component alert rates",
  "",
  paste(capture.output(print(components)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "The simulated complete-feature external assignment is a known-DGP reference only. It does not remove the empirical eICU REFERENCE_ABSENT state. The run evaluates the regenerated surrogate interface and the Baseline-plus-generated-K2 calibration engine; it does not replace or overwrite the legacy G3 results."
)
d3r_atomic_write_lines(report, file.path(D3R_OUTPUT_ROOT, "D3_SIMULATION_RESULTS.md"), overwrite = TRUE)
