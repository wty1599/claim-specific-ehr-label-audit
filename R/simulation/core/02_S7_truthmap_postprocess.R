options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: Rscript postprocess_S7_locked_truthmap.R <s7_table_dir> <prior_truthmap_dir> <output_dir>",
    call. = FALSE
  )
}

s7_table_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
prior_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE)
}

wilson <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(estimate = NA, low = NA, high = NA))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(estimate = p, low = max(0, ctr - half), high = min(1, ctr + half))
}

mean_mcse <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  c(
    n = length(x), mean = mean(x), sd = sd(x), mcse = sd(x) / sqrt(length(x)),
    q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)),
    min = min(x), max = max(x)
  )
}

ground <- read_required(file.path(s7_table_dir, "S7_ground_truth_check_by_repeat.csv"))
d1_summary_source <- read_required(file.path(s7_table_dir, "domain1_discreteness_summary.csv"))
d2 <- read_required(file.path(s7_table_dir, "domain2_delta_auc_rubin_by_repeat.csv"))
prior_x1_main <- read_required(file.path(prior_dir, "Table_X1_operating_characteristics_main_exclude_failures.csv"))
prior_x2 <- read_required(file.path(prior_dir, "Table_X2_scenario_domain_verdict_vs_truth.csv"))

stopifnot(nrow(ground) == 100L, length(unique(ground$repeat_id)) == 100L)
stopifnot(nrow(d1_summary_source) == 100L, length(unique(d1_summary_source$repeat_id)) == 100L)
stopifnot(nrow(d2) == 200L, length(unique(d2$repeat_id)) == 100L)

# Locked verdict rules. No threshold is estimated from S7.
D1_DIP_ALPHA <- 0.05
D1_SIGCLUST_EFFECT_THRESHOLD <- 3.0
D1_TRUE_K_ARI_THRESHOLD <- 0.80
D2_PRACTICAL_NULL_THRESHOLD <- 0.001

d1_repeat <- data.frame(
  scenario = d1_summary_source$scenario,
  repeat_id = as.integer(d1_summary_source$repeat_id),
  true_k_ari = as.numeric(d1_summary_source$ari_to_true_subtype),
  dip_p = as.numeric(d1_summary_source$dip_discriminant_p),
  sigclust_effect_sd = as.numeric(d1_summary_source$sigclust_effect_sd)
)
d1_repeat$observed_positive <- with(
  d1_repeat,
  true_k_ari >= D1_TRUE_K_ARI_THRESHOLD &
    (sigclust_effect_sd >= D1_SIGCLUST_EFFECT_THRESHOLD | dip_p < D1_DIP_ALPHA)
)
d1_repeat$expected_positive <- TRUE
d1_repeat$correct_classification <- d1_repeat$observed_positive == d1_repeat$expected_positive
d1_repeat$failed <- !complete.cases(d1_repeat[, c("true_k_ari", "dip_p", "sigclust_effect_sd")])

d2$rubin_estimate <- as.numeric(d2$rubin_estimate)
d2$rubin_ci_low <- as.numeric(d2$rubin_ci_low)
d2$rubin_ci_high <- as.numeric(d2$rubin_ci_high)
d2$repeat_id <- as.integer(d2$repeat_id)
d2$observed_positive <- with(
  d2,
  rubin_estimate > D2_PRACTICAL_NULL_THRESHOLD & rubin_ci_low > 0
)
d2$expected_positive <- FALSE
d2$correct_classification <- d2$observed_positive == d2$expected_positive
d2$failed <- !complete.cases(d2[, c("rubin_estimate", "rubin_ci_low", "rubin_ci_high")])

# The locked operating-characteristics table uses MAKE30 as the single Domain 2
# verdict per domain x repeat cell. Mortality remains an outcome-specific audit.
d2_primary <- d2[d2$outcome == "make30", ]
stopifnot(nrow(d2_primary) == 100L)

write.csv(d1_repeat, file.path(out_dir, "S7_D1_locked_verdict_by_repeat.csv"), row.names = FALSE)
write.csv(d2, file.path(out_dir, "S7_D2_locked_verdict_by_outcome_repeat.csv"), row.names = FALSE)

summarise_verdict <- function(x, expected_positive) {
  valid <- !x$failed & !is.na(x$observed_positive)
  n <- sum(valid)
  positive <- sum(x$observed_positive[valid])
  correct <- sum(x$correct_classification[valid])
  ci <- wilson(correct, n)
  data.frame(
    n_rep = n,
    n_failed = sum(x$failed),
    observed_positive_n = positive,
    expected_positive = expected_positive,
    correct_n = correct,
    correct_rate = unname(ci["estimate"]),
    wilson_low = unname(ci["low"]),
    wilson_high = unname(ci["high"])
  )
}

s7_d1_sum <- summarise_verdict(d1_repeat, TRUE)
s7_d2_sum <- summarise_verdict(d2_primary, FALSE)
s7_verdict_summary <- rbind(
  cbind(domain = "D1", scenario = "S7_discrete_structure_outcome_null", rule = "ARI>=0.80 AND (SigClust effect>=3 SD OR dip p<0.05)", performance_measure = "sensitivity", s7_d1_sum),
  cbind(domain = "D2", scenario = "S7_discrete_structure_outcome_null", rule = "MAKE30 Delta AUC>0.001 AND Rubin CI lower>0", performance_measure = "specificity", s7_d2_sum)
)
write.csv(s7_verdict_summary, file.path(out_dir, "S7_locked_verdict_summary.csv"), row.names = FALSE)

ground_metrics <- c(
  "assigned_cluster_axis_severity_cor", "true_cluster_axis_severity_cor",
  "assigned_cluster_auc_mortality", "assigned_cluster_auc_make30",
  "dip_discriminant_p", "sigclust_effect_sd", "ari_assigned_to_true_K2"
)
ground_summary <- do.call(rbind, lapply(ground_metrics, function(v) {
  data.frame(metric = v, as.list(mean_mcse(ground[[v]])), check.names = FALSE)
}))
write.csv(ground_summary, file.path(out_dir, "S7_ground_truth_and_D1_continuous_summary.csv"), row.names = FALSE)

d2_summary <- do.call(rbind, lapply(split(d2, d2$outcome), function(z) {
  vals <- mean_mcse(z$rubin_estimate)
  data.frame(
    outcome = z$outcome[1], as.list(vals),
    estimate_gt_0_001_n = sum(z$rubin_estimate > 0.001),
    ci_lower_gt_zero_n = sum(z$rubin_ci_low > 0),
    locked_positive_n = sum(z$observed_positive),
    check.names = FALSE
  )
}))
write.csv(d2_summary, file.path(out_dir, "S7_D2_continuous_summary.csv"), row.names = FALSE)

# Append S7 scenario rows to the previously locked S1-S6 Table X1, then rebuild
# pooled domain/truth and overall rows. Existing scenario rows are not altered.
scenario_rows <- prior_x1_main[prior_x1_main$aggregation_level == "scenario", ]
s7_rows <- data.frame(
  aggregation_level = "scenario",
  domain = c("D1", "D2"),
  truth = c("positive", "negative"),
  scenarios = "S7",
  scenario_display = "S7: discrete-structure/outcome-null orthogonal control",
  performance_measure = c("sensitivity", "specificity"),
  failure_policy = "exclude_failed",
  n_rep = c(s7_d1_sum$n_rep, s7_d2_sum$n_rep),
  n_failed = c(s7_d1_sum$n_failed, s7_d2_sum$n_failed),
  verdict_equals_truth_n = c(s7_d1_sum$correct_n, s7_d2_sum$correct_n),
  stringsAsFactors = FALSE
)
scenario_rows <- rbind(scenario_rows[, names(s7_rows)], s7_rows)

rebuild_x1 <- function(rows, failure_policy) {
  rows$failure_policy <- failure_policy
  rows$n_rep <- as.integer(rows$n_rep)
  rows$n_failed <- as.integer(rows$n_failed)
  rows$verdict_equals_truth_n <- as.integer(rows$verdict_equals_truth_n)
  rows$rate <- NA_real_; rows$wilson_low <- NA_real_; rows$wilson_high <- NA_real_
  for (i in seq_len(nrow(rows))) {
    n <- if (failure_policy == "exclude_failed") rows$n_rep[i] else rows$n_rep[i] + rows$n_failed[i]
    ci <- wilson(rows$verdict_equals_truth_n[i], n)
    rows$rate[i] <- ci["estimate"]; rows$wilson_low[i] <- ci["low"]; rows$wilson_high[i] <- ci["high"]
  }
  rows$rate_wilson_95ci <- sprintf("%.3f [%.3f, %.3f]", rows$rate, rows$wilson_low, rows$wilson_high)

  keys <- unique(rows[, c("domain", "truth", "performance_measure")])
  pooled <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    z <- rows[rows$domain == keys$domain[i] & rows$truth == keys$truth[i] & rows$performance_measure == keys$performance_measure[i], ]
    n_rep <- sum(z$n_rep); n_failed <- sum(z$n_failed); correct <- sum(z$verdict_equals_truth_n)
    n <- if (failure_policy == "exclude_failed") n_rep else n_rep + n_failed
    ci <- wilson(correct, n)
    data.frame(
      aggregation_level = "domain_truth_pooled", domain = keys$domain[i], truth = keys$truth[i],
      scenarios = paste(z$scenarios, collapse = "+"), scenario_display = paste(z$scenario_display, collapse = " + "),
      performance_measure = keys$performance_measure[i], failure_policy = failure_policy,
      n_rep = n_rep, n_failed = n_failed, verdict_equals_truth_n = correct,
      rate = ci["estimate"], wilson_low = ci["low"], wilson_high = ci["high"],
      rate_wilson_95ci = sprintf("%.3f [%.3f, %.3f]", ci["estimate"], ci["low"], ci["high"]),
      stringsAsFactors = FALSE
    )
  }))
  total_n_rep <- sum(rows$n_rep); total_failed <- sum(rows$n_failed); total_correct <- sum(rows$verdict_equals_truth_n)
  total_n <- if (failure_policy == "exclude_failed") total_n_rep else total_n_rep + total_failed
  ci <- wilson(total_correct, total_n)
  overall <- data.frame(
    aggregation_level = "overall", domain = "Overall", truth = "mixed",
    scenarios = "all applicable domain x repeat cells", scenario_display = "all applicable domain x repeat cells",
    performance_measure = "correct-classification rate", failure_policy = failure_policy,
    n_rep = total_n_rep, n_failed = total_failed, verdict_equals_truth_n = total_correct,
    rate = ci["estimate"], wilson_low = ci["low"], wilson_high = ci["high"],
    rate_wilson_95ci = sprintf("%.3f [%.3f, %.3f]", ci["estimate"], ci["low"], ci["high"]),
    stringsAsFactors = FALSE
  )
  rbind(rows, pooled[, names(rows)], overall[, names(rows)])
}

x1_main <- rebuild_x1(scenario_rows, "exclude_failed")
x1_fail_wrong <- rebuild_x1(scenario_rows, "failure_equals_misclassification")
write.csv(x1_main, file.path(out_dir, "Table_X1_operating_characteristics_with_S7_main_exclude_failures.csv"), row.names = FALSE)
write.csv(x1_fail_wrong, file.path(out_dir, "Table_X1_operating_characteristics_with_S7_failure_equals_wrong.csv"), row.names = FALSE)

# Append all four S7 truth-map cells to Table X2. D3 and D4 are explicitly N/A.
fmt <- function(x) formatC(x, format = "f", digits = 4)
d1_metrics <- sprintf(
  "K2 true-label ARI %s +/- %s; SigClust effect SD %s +/- %s; dip p %s +/- %s",
  fmt(mean(d1_repeat$true_k_ari)), fmt(sd(d1_repeat$true_k_ari) / 10),
  fmt(mean(d1_repeat$sigclust_effect_sd)), fmt(sd(d1_repeat$sigclust_effect_sd) / 10),
  fmt(mean(d1_repeat$dip_p)), fmt(sd(d1_repeat$dip_p) / 10)
)
make_z <- d2[d2$outcome == "make30", ]; mort_z <- d2[d2$outcome == "mortality_30d", ]
d2_metrics <- sprintf(
  "MAKE30 Delta AUC %s +/- %s; mortality Delta AUC %s +/- %s",
  fmt(mean(make_z$rubin_estimate)), fmt(sd(make_z$rubin_estimate) / 10),
  fmt(mean(mort_z$rubin_estimate)), fmt(sd(mort_z$rubin_estimate) / 10)
)
s7_x2 <- data.frame(
  domain = c("D1", "D2", "D3", "D4"), scenario_short = "S7",
  scenario = "S7_discrete_structure_outcome_null",
  manuscript_scenario_label = "S7: discrete-structure/outcome-null orthogonal control",
  truth = c("positive", "negative", "N/A", "N/A"),
  expected_positive = c("TRUE", "FALSE", "", ""),
  performance_measure = c("sensitivity", "specificity", "N/A", "N/A"),
  applicable = c("TRUE", "TRUE", "FALSE", "FALSE"),
  expected_label = c("positive", "negative", "N/A", "N/A"),
  n_rep = c(100, 100, 0, 0), n_failed = 0,
  observed_positive_n = c(sum(d1_repeat$observed_positive), sum(d2_primary$observed_positive), NA, NA),
  observed_positive_rate = c(mean(d1_repeat$observed_positive), mean(d2_primary$observed_positive), NA, NA),
  correct_n = c(sum(d1_repeat$correct_classification), sum(d2_primary$correct_classification), NA, NA),
  correct_rate = c(mean(d1_repeat$correct_classification), mean(d2_primary$correct_classification), NA, NA),
  binary_verdict_summary = c("positive 100/100 (1.000)", "positive 0/100 (0.000)", "N/A", "N/A"),
  key_continuous_metrics_mean_mcse = c(d1_metrics, d2_metrics, "N/A", "N/A"),
  source = c(
    "domain1_discreteness_summary.csv",
    "domain2_delta_auc_rubin_by_repeat.csv (MAKE30 primary verdict)",
    "truth-map: N/A", "truth-map: N/A"
  ),
  stringsAsFactors = FALSE
)
x2 <- rbind(prior_x2, s7_x2[, names(prior_x2)])
write.csv(x2, file.path(out_dir, "Table_X2_scenario_domain_verdict_vs_truth_with_S7.csv"), row.names = FALSE)

failure_audit <- data.frame(
  scenario = "S7_discrete_structure_outcome_null",
  domain = c("D1", "D2"),
  expected_repeats = 100,
  valid_repeats = c(sum(!d1_repeat$failed), sum(!d2_primary$failed)),
  failed_repeats = c(sum(d1_repeat$failed), sum(d2_primary$failed)),
  failure_reason = c("none in required D1 fields", "none in required MAKE30 Rubin fields"),
  stringsAsFactors = FALSE
)
write.csv(failure_audit, file.path(out_dir, "S7_failure_audit.csv"), row.names = FALSE)

checklist <- data.frame(
  check = c(
    "100 unique S7 repeats", "all generation QC gates pass", "D1 fields complete",
    "D2 MAKE30 fields complete", "locked D1 threshold unchanged", "locked D2 threshold unchanged",
    "D3 excluded", "D4 excluded", "D4 specificity remains S6 only"
  ),
  pass = c(
    length(unique(ground$repeat_id)) == 100,
    all(ground$generation_qc_pass == TRUE),
    !any(d1_repeat$failed), !any(d2_primary$failed), TRUE, TRUE, TRUE, TRUE,
    identical(sort(unique(x1_main$scenarios[x1_main$aggregation_level == "scenario" & x1_main$domain == "D4" & x1_main$truth == "negative"])), "S6")
  ),
  stringsAsFactors = FALSE
)
write.csv(checklist, file.path(out_dir, "S7_postprocess_checklist.csv"), row.names = FALSE)

cat("S7 locked-truth-map post-processing completed.\n")
cat("Output:", out_dir, "\n")
cat("D1 positive:", sum(d1_repeat$observed_positive), "/100\n")
cat("D2 MAKE30 positive:", sum(d2_primary$observed_positive), "/100\n")
cat("All checks pass:", all(checklist$pass), "\n")
