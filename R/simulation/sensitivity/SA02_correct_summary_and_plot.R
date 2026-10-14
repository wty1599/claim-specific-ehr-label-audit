#!/usr/bin/env Rscript

# Sensitivity 02 post-processing only (base R; no extra packages required).
# It reads completed repeat-level Rubin results, joins oracle accuracy from
# the task manifest, and writes corrected outputs without changing sources.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) if (!is.null(x$ofile)) as.character(x$ofile) else NA_character_,
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[length(frame_files)], winslash = "/")))
  }

  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
analysis_override <- Sys.getenv("SENSITIVITY02_DIR", unset = "")
analysis_dir <- if (nzchar(analysis_override)) {
  normalizePath(analysis_override, winslash = "/", mustWork = TRUE)
} else {
  script_dir
}

source_table_dir <- file.path(analysis_dir, "output", "tables")
rubin_file <- file.path(source_table_dir, "domain2_delta_auc_rubin_by_repeat.csv")
manifest_file <- file.path(source_table_dir, "sensitivity_task_manifest.csv")

if (!file.exists(rubin_file)) {
  stop("Repeat-level Rubin file not found: ", rubin_file, call. = FALSE)
}
if (!file.exists(manifest_file)) {
  stop(
    "Task manifest not found: ", manifest_file,
    ". Oracle accuracy is not stored in the Rubin file.",
    call. = FALSE
  )
}

output_root <- file.path(analysis_dir, "output", "corrected_domain2_v2")
output_table_dir <- file.path(output_root, "tables")
output_figure_dir <- file.path(output_root, "figures")
output_log_dir <- file.path(output_root, "logs")
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_log_dir, recursive = TRUE, showWarnings = FALSE)

practical_null_threshold <- 0.001
practical_null_thresholds <- c(0.001, 0.002, 0.005)

wilson_interval <- function(successes, n, conf_level = 0.95) {
  if (n <= 0L) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf_level) / 2)
  p <- successes / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  pmax(0, pmin(1, c(center - half, center + half)))
}

rubin <- read.csv(
  rubin_file,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN"),
  check.names = FALSE
)
manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN"),
  check.names = FALSE
)

required_rubin_columns <- c(
  "repeat_id", "outcome", "rubin_estimate", "rubin_ci_low", "rubin_ci_high",
  "m_used", "variance_method", "domain2_label_source"
)
required_manifest_columns <- c(
  "rep_id", "local_rep_id", "sensitivity_id",
  "S5_ORACLE_LABEL_ACCURACY_value", "scenario"
)

missing_rubin_columns <- setdiff(required_rubin_columns, names(rubin))
missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_rubin_columns)) {
  stop("Rubin file is missing: ", paste(missing_rubin_columns, collapse = ", "), call. = FALSE)
}
if (length(missing_manifest_columns)) {
  stop("Manifest is missing: ", paste(missing_manifest_columns, collapse = ", "), call. = FALSE)
}

rubin$repeat_id <- as.character(rubin$repeat_id)
manifest$rep_id <- as.character(manifest$rep_id)

manifest_map <- unique(data.frame(
  repeat_id = manifest$rep_id,
  local_repeat_id = as.integer(manifest$local_rep_id),
  sensitivity_id = manifest$sensitivity_id,
  oracle_accuracy = as.numeric(manifest$S5_ORACLE_LABEL_ACCURACY_value),
  manifest_scenario = manifest$scenario,
  stringsAsFactors = FALSE
))

if (anyDuplicated(manifest_map$repeat_id)) {
  stop("Manifest rep_id is not unique; the accuracy mapping is ambiguous.", call. = FALSE)
}

repeat_level <- merge(
  rubin,
  manifest_map,
  by = "repeat_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(repeat_level$oracle_accuracy)) {
  missing_ids <- unique(repeat_level$repeat_id[is.na(repeat_level$oracle_accuracy)])
  stop(
    "Oracle accuracy could not be mapped for repeat_id(s): ",
    paste(head(missing_ids, 10), collapse = ", "),
    call. = FALSE
  )
}

key <- paste(
  repeat_level$sensitivity_id,
  repeat_level$oracle_accuracy,
  repeat_level$outcome,
  repeat_level$local_repeat_id,
  sep = "|"
)
if (anyDuplicated(key)) {
  stop("Expected one Rubin result per setting, outcome, and repeat.", call. = FALSE)
}

repeat_level$outcome_label <- ifelse(
  repeat_level$outcome == "make30",
  "MAKE30",
  ifelse(
    repeat_level$outcome == "mortality_30d",
    "30-day mortality",
    repeat_level$outcome
  )
)

group_key <- interaction(
  repeat_level$sensitivity_id,
  repeat_level$oracle_accuracy,
  repeat_level$outcome,
  drop = TRUE,
  lex.order = TRUE
)
groups <- split(seq_len(nrow(repeat_level)), group_key)

summarize_group <- function(idx) {
  x <- repeat_level[idx, , drop = FALSE]
  n <- nrow(x)
  mean_est <- mean(x$rubin_estimate, na.rm = TRUE)
  sd_est <- sd(x$rubin_estimate, na.rm = TRUE)
  mcse <- sd_est / sqrt(n)
  crit <- qt(0.975, df = n - 1L)
  ci_low <- mean_est - crit * mcse
  ci_high <- mean_est + crit * mcse
  practical_null <- abs(mean_est) < practical_null_threshold

  classification <- if (practical_null) {
    "practically_null"
  } else if (ci_low > 0) {
    "positive_incremental_value"
  } else if (ci_high < 0) {
    "negative_incremental_value"
  } else {
    "mean_ci_includes_zero"
  }

  expected_role <- if (isTRUE(all.equal(x$oracle_accuracy[1], 0.50))) {
    "chance_level_null_control"
  } else {
    "dose_response_positive_control"
  }
  expected_classification <- if (expected_role == "chance_level_null_control") {
    "practically_null"
  } else {
    "positive_incremental_value"
  }

  positive_n <- sum(x$rubin_ci_low > 0, na.rm = TRUE)
  negative_n <- sum(x$rubin_ci_high < 0, na.rm = TRUE)
  positive_ci <- wilson_interval(positive_n, n)
  negative_ci <- wilson_interval(negative_n, n)
  positive_rate <- positive_n / n

  study_level_interpretation <- if (expected_role == "chance_level_null_control") {
    "chance_level_practically_null"
  } else if (x$oracle_accuracy[1] <= 0.60) {
    "weak_non_null_signal_low_study_level_detectability"
  } else if (positive_rate < 0.80) {
    "non_null_signal_incomplete_study_level_detectability"
  } else if (positive_rate < 0.95) {
    "non_null_signal_reliably_detectable"
  } else {
    "non_null_signal_consistently_detectable"
  }

  data.frame(
    sensitivity_id = x$sensitivity_id[1],
    oracle_accuracy = x$oracle_accuracy[1],
    outcome = x$outcome[1],
    outcome_label = x$outcome_label[1],
    n_repeats = n,
    m_used_min = min(x$m_used, na.rm = TRUE),
    m_used_max = max(x$m_used, na.rm = TRUE),
    mean_delta_auc = mean_est,
    median_delta_auc = median(x$rubin_estimate, na.rm = TRUE),
    sd_across_repeats = sd_est,
    mcse_mean = mcse,
    mean_ci_low = ci_low,
    mean_ci_high = ci_high,
    empirical_p025 = unname(quantile(x$rubin_estimate, 0.025, na.rm = TRUE)),
    empirical_p975 = unname(quantile(x$rubin_estimate, 0.975, na.rm = TRUE)),
    repeat_delta_auc_positive_proportion = mean(x$rubin_estimate > 0, na.rm = TRUE),
    repeat_rubin_ci_positive_n = positive_n,
    repeat_rubin_ci_positive_proportion = positive_rate,
    repeat_rubin_ci_positive_ci_low = positive_ci[1],
    repeat_rubin_ci_positive_ci_high = positive_ci[2],
    repeat_rubin_ci_negative_n = negative_n,
    repeat_rubin_ci_negative_proportion = negative_n / n,
    repeat_rubin_ci_negative_ci_low = negative_ci[1],
    repeat_rubin_ci_negative_ci_high = negative_ci[2],
    repeat_rubin_ci_includes_zero_proportion = mean(
      x$rubin_ci_low <= 0 & x$rubin_ci_high >= 0,
      na.rm = TRUE
    ),
    practical_null_threshold = practical_null_threshold,
    practical_null = practical_null,
    classification = classification,
    expected_role = expected_role,
    expected_classification = expected_classification,
    correct_classification = identical(classification, expected_classification),
    study_level_interpretation = study_level_interpretation,
    variance_method = paste(sort(unique(x$variance_method)), collapse = ";"),
    domain2_label_source = paste(sort(unique(x$domain2_label_source)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

setting_summary <- do.call(rbind, lapply(groups, summarize_group))
setting_summary <- setting_summary[
  order(setting_summary$oracle_accuracy, setting_summary$outcome),
  ,
  drop = FALSE
]
rownames(setting_summary) <- NULL

if (nrow(setting_summary) != 12L) {
  stop("Corrected summary must contain exactly 12 rows; found ", nrow(setting_summary), call. = FALSE)
}
if (any(setting_summary$n_repeats != 100L)) {
  stop("Each accuracy-outcome setting must contain 100 repeats.", call. = FALSE)
}
if (any(setting_summary$m_used_min != 5L | setting_summary$m_used_max != 5L)) {
  stop("Not all repeat-level estimates used five imputations.", call. = FALSE)
}

summary_file <- file.path(
  output_table_dir,
  "domain2_delta_auc_setting_summary_corrected_12rows.csv"
)
source_data_file <- file.path(
  output_table_dir,
  "domain2_delta_auc_repeat_level_with_oracle_accuracy.csv"
)
write.csv(setting_summary, summary_file, row.names = FALSE, na = "")
write.csv(repeat_level, source_data_file, row.names = FALSE, na = "")

# Practical-null threshold sensitivity. This changes interpretation only and
# does not alter any fitted model or repeat-level estimate.
threshold_sensitivity <- do.call(rbind, lapply(practical_null_thresholds, function(threshold) {
  out <- setting_summary[, c(
    "sensitivity_id", "oracle_accuracy", "outcome", "outcome_label",
    "n_repeats", "mean_delta_auc", "mean_ci_low", "mean_ci_high",
    "repeat_rubin_ci_positive_proportion"
  )]
  out$practical_null_threshold <- threshold
  out$practical_null <- abs(out$mean_delta_auc) < threshold
  out$threshold_classification <- ifelse(
    out$practical_null,
    "practically_null",
    ifelse(
      out$mean_ci_low > 0,
      "positive_incremental_value",
      ifelse(out$mean_ci_high < 0, "negative_incremental_value", "mean_ci_includes_zero")
    )
  )
  out
}))
threshold_sensitivity <- threshold_sensitivity[
  order(
    threshold_sensitivity$practical_null_threshold,
    threshold_sensitivity$oracle_accuracy,
    threshold_sensitivity$outcome
  ),
  , drop = FALSE
]
write.csv(
  threshold_sensitivity,
  file.path(output_table_dir, "domain2_practical_null_threshold_sensitivity.csv"),
  row.names = FALSE, na = ""
)

# Seed-reconstructed oracle-label fidelity QC. Individual latent/oracle labels
# were not retained in the original result object. The original deterministic
# seed contract is therefore replayed without regenerating observed features,
# missingness, MICE datasets, outcomes, or Domain 1-4 models.
global_seed <- 20250101L
scenario_index <- 5L
latent_prevalence_target <- 0.35

reconstruct_oracle_qc <- function(manifest_row, cohort, n, cohort_offset) {
  global_repeat_id <- as.integer(manifest_row$rep_id)
  target_accuracy <- as.numeric(manifest_row$S5_ORACLE_LABEL_ACCURACY_value)
  seed_base <- global_seed + global_repeat_id * 1000L + scenario_index * 100000L
  dataset_seed <- seed_base + cohort_offset

  set.seed(dataset_seed)
  latent_z <- rnorm(n)
  latent_state <- as.integer(latent_z > qnorm(1 - latent_prevalence_target))

  set.seed(dataset_seed + 909L)
  flipped <- rbinom(n, size = 1L, prob = 1 - target_accuracy) == 1L
  oracle_state <- ifelse(flipped, 1L - latent_state, latent_state)

  positive_n <- sum(latent_state == 1L)
  negative_n <- sum(latent_state == 0L)
  data.frame(
    sensitivity_id = manifest_row$sensitivity_id,
    global_repeat_id = global_repeat_id,
    local_repeat_id = as.integer(manifest_row$local_rep_id),
    cohort = cohort,
    n = n,
    target_accuracy = target_accuracy,
    realized_accuracy = mean(oracle_state == latent_state),
    target_minus_realized_accuracy = target_accuracy - mean(oracle_state == latent_state),
    oracle_label_positive_prevalence = mean(oracle_state == 1L),
    latent_state_positive_prevalence = mean(latent_state == 1L),
    sensitivity_vs_latent_state = if (positive_n) {
      mean(oracle_state[latent_state == 1L] == 1L)
    } else NA_real_,
    specificity_vs_latent_state = if (negative_n) {
      mean(oracle_state[latent_state == 0L] == 0L)
    } else NA_real_,
    reconstruction_source = "deterministic_seed_contract_label_only_no_model_rerun",
    stringsAsFactors = FALSE
  )
}

manifest_unique <- manifest[!duplicated(manifest$rep_id), , drop = FALSE]
oracle_qc_repeat <- do.call(rbind, lapply(seq_len(nrow(manifest_unique)), function(i) {
  row <- manifest_unique[i, , drop = FALSE]
  rbind(
    reconstruct_oracle_qc(row, "derivation", as.integer(row$N_TRAIN_value), 1L),
    reconstruct_oracle_qc(row, "external", as.integer(row$N_EXTERNAL_value), 2L)
  )
}))
rownames(oracle_qc_repeat) <- NULL

qc_group <- interaction(
  oracle_qc_repeat$sensitivity_id,
  oracle_qc_repeat$target_accuracy,
  oracle_qc_repeat$cohort,
  drop = TRUE,
  lex.order = TRUE
)
oracle_qc_summary <- do.call(rbind, lapply(split(seq_len(nrow(oracle_qc_repeat)), qc_group), function(idx) {
  d <- oracle_qc_repeat[idx, , drop = FALSE]
  summarize_metric <- function(x, prefix) {
    n <- sum(!is.na(x))
    s <- sd(x, na.rm = TRUE)
    setNames(
      c(
        mean(x, na.rm = TRUE), s,
        unname(quantile(x, 0.025, na.rm = TRUE, type = 8)),
        unname(quantile(x, 0.975, na.rm = TRUE, type = 8)),
        s / sqrt(n)
      ),
      paste0(prefix, c("_mean", "_sd", "_empirical_p025", "_empirical_p975", "_mcse"))
    )
  }
  metrics <- c(
    summarize_metric(d$realized_accuracy, "realized_accuracy"),
    summarize_metric(d$target_minus_realized_accuracy, "target_minus_realized_accuracy"),
    summarize_metric(d$oracle_label_positive_prevalence, "oracle_label_positive_prevalence"),
    summarize_metric(d$latent_state_positive_prevalence, "latent_state_positive_prevalence"),
    summarize_metric(d$sensitivity_vs_latent_state, "sensitivity_vs_latent_state"),
    summarize_metric(d$specificity_vs_latent_state, "specificity_vs_latent_state")
  )
  data.frame(
    sensitivity_id = d$sensitivity_id[1],
    target_accuracy = d$target_accuracy[1],
    cohort = d$cohort[1],
    n_repeats = nrow(d),
    as.list(metrics),
    reconstruction_source = d$reconstruction_source[1],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(oracle_qc_summary) <- NULL
oracle_qc_summary <- oracle_qc_summary[
  order(oracle_qc_summary$target_accuracy, oracle_qc_summary$cohort),
  , drop = FALSE
]

write.csv(
  oracle_qc_repeat,
  file.path(output_table_dir, "oracle_label_fidelity_qc_by_repeat_seed_reconstructed.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  oracle_qc_summary,
  file.path(output_table_dir, "oracle_label_fidelity_qc_setting_summary_12rows.csv"),
  row.names = FALSE, na = ""
)

deprecated_outputs <- data.frame(
  file = c(
    file.path(source_table_dir, "simulation_four_domain_aggregate_summary.csv"),
    file.path(source_table_dir, "domain2_delta_auc_mean_ci_flags.csv")
  ),
  status = "deprecated_grouping_error",
  reason = c(
    "sens_local_rep_id was retained in grouping; nonmissing_n=1 and MC-SE is unavailable",
    "setting-level grouping retained repeat identifiers and does not produce the required 12-row summary"
  ),
  allowed_for_final_analysis = FALSE,
  replacement = c(summary_file, summary_file),
  stringsAsFactors = FALSE
)
write.csv(
  deprecated_outputs,
  file.path(output_table_dir, "deprecated_source_outputs_do_not_use.csv"),
  row.names = FALSE, na = ""
)

canonical_outputs <- data.frame(
  role = c(
    "repeat_level_delta_auc", "setting_level_delta_auc", "operating_characteristics",
    "threshold_sensitivity", "oracle_fidelity_repeat", "oracle_fidelity_summary"
  ),
  file = c(
    source_data_file,
    summary_file,
    summary_file,
    file.path(output_table_dir, "domain2_practical_null_threshold_sensitivity.csv"),
    file.path(output_table_dir, "oracle_label_fidelity_qc_by_repeat_seed_reconstructed.csv"),
    file.path(output_table_dir, "oracle_label_fidelity_qc_setting_summary_12rows.csv")
  ),
  stringsAsFactors = FALSE
)
write.csv(
  canonical_outputs,
  file.path(output_table_dir, "sensitivity02_canonical_output_manifest.csv"),
  row.names = FALSE, na = ""
)

scope_note <- c(
  "Sensitivity 02 scope and interpretation",
  "",
  "This analysis is an internal Domain 2 metric-level positive control.",
  "The noisy oracle label is not an ordinary k-means K=2 phenotype label.",
  "It is generated by symmetrically corrupting an unobserved binary latent state",
  "at the prespecified accuracy and is fixed across all five imputed datasets.",
  "",
  "The external-cohort rows in the oracle fidelity QC verify only that the",
  "deterministic label-generation mechanism achieved its target properties.",
  "They do not constitute external predictive validation or an external Delta-AUC",
  "dose-response analysis. Existing Domain 3 outputs must not be interpreted as",
  "an external oracle-label positive control.",
  "",
  "A separate Sensitivity 02b would be required for external oracle-label",
  "incremental discrimination, calibration, Brier score, or IDI analyses."
)
writeLines(
  scope_note,
  file.path(output_log_dir, "sensitivity02_scope_and_interpretation.txt")
)

# Figure contract:
# Claim: Domain 2 responds monotonically to genuine label-specific prognostic
# information, but reliable study-level detection emerges only at stronger
# oracle-label accuracy (approximately 0.80-0.90).
# Evidence: mean Delta AUC with Monte Carlo mean CIs and repeat-level Rubin-CI
# positive rates with binomial Monte Carlo CIs.

palette <- c(
  "MAKE30" = "#2F6B9A",
  "30-day mortality" = "#C65A46"
)
panel_order <- c("MAKE30", "30-day mortality")

draw_delta_panel <- function(label, panel_tag) {
  d <- setting_summary[setting_summary$outcome_label == label, , drop = FALSE]
  d <- d[order(d$oracle_accuracy), , drop = FALSE]
  col <- unname(palette[label])
  y_limits <- range(
    c(setting_summary$mean_ci_low, setting_summary$mean_ci_high,
      -practical_null_threshold, practical_null_threshold),
    finite = TRUE
  )
  y_padding <- diff(y_limits) * 0.08
  y_limits <- y_limits + c(-y_padding, y_padding)

  plot(
    d$oracle_accuracy,
    d$mean_delta_auc,
    type = "n",
    xlim = c(0.48, 1.02),
    ylim = y_limits,
    xaxt = "n",
    xlab = "",
    ylab = if (label == panel_order[1]) {
      "Mean Delta AUC\n(Raw-EN + oracle label minus Raw-EN)"
    } else "",
    main = label,
    cex.main = 0.96,
    cex.lab = 0.82,
    cex.axis = 0.78,
    bty = "l"
  )
  rect(
    par("usr")[1], -practical_null_threshold,
    par("usr")[2], practical_null_threshold,
    col = adjustcolor("#BDBDBD", alpha.f = 0.34), border = NA
  )
  abline(h = 0, lty = 2, lwd = 0.7, col = "#4D4D4D")
  axis(
    1, at = sort(unique(setting_summary$oracle_accuracy)),
    labels = sprintf("%.2f", sort(unique(setting_summary$oracle_accuracy))),
    cex.axis = 0.78
  )
  segments(
    d$oracle_accuracy, d$mean_ci_low,
    d$oracle_accuracy, d$mean_ci_high,
    col = col, lwd = 1.05
  )
  segments(
    d$oracle_accuracy - 0.007, d$mean_ci_low,
    d$oracle_accuracy + 0.007, d$mean_ci_low,
    col = col, lwd = 1.05
  )
  segments(
    d$oracle_accuracy - 0.007, d$mean_ci_high,
    d$oracle_accuracy + 0.007, d$mean_ci_high,
    col = col, lwd = 1.05
  )
  points(
    d$oracle_accuracy, d$mean_delta_auc,
    pch = 21, bg = "white", col = col, lwd = 1.15, cex = 1.05
  )
  mtext(panel_tag, side = 3, adj = 0, line = 0.35, font = 2, cex = 0.95)
}

draw_detection_panel <- function(label, panel_tag) {
  d <- setting_summary[setting_summary$outcome_label == label, , drop = FALSE]
  d <- d[order(d$oracle_accuracy), , drop = FALSE]
  col <- unname(palette[label])

  plot(
    d$oracle_accuracy,
    d$repeat_rubin_ci_positive_proportion,
    type = "n",
    xlim = c(0.48, 1.02),
    ylim = c(0, 1.05),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = if (label == panel_order[1]) {
      "Repeat-level positive detection rate"
    } else "",
    main = label,
    cex.main = 0.96,
    cex.lab = 0.82,
    cex.axis = 0.78,
    bty = "l"
  )
  abline(h = c(0.80, 0.90), lty = c(3, 2), lwd = 0.6, col = c("#A0A0A0", "#707070"))
  axis(
    1, at = sort(unique(setting_summary$oracle_accuracy)),
    labels = sprintf("%.2f", sort(unique(setting_summary$oracle_accuracy))),
    cex.axis = 0.78
  )
  axis(2, at = seq(0, 1, 0.2), labels = sprintf("%d%%", seq(0, 100, 20)), cex.axis = 0.78)
  segments(
    d$oracle_accuracy, d$repeat_rubin_ci_positive_ci_low,
    d$oracle_accuracy, d$repeat_rubin_ci_positive_ci_high,
    col = col, lwd = 1.05
  )
  segments(
    d$oracle_accuracy - 0.007, d$repeat_rubin_ci_positive_ci_low,
    d$oracle_accuracy + 0.007, d$repeat_rubin_ci_positive_ci_low,
    col = col, lwd = 1.05
  )
  segments(
    d$oracle_accuracy - 0.007, d$repeat_rubin_ci_positive_ci_high,
    d$oracle_accuracy + 0.007, d$repeat_rubin_ci_positive_ci_high,
    col = col, lwd = 1.05
  )
  points(
    d$oracle_accuracy, d$repeat_rubin_ci_positive_proportion,
    pch = 21, bg = "white", col = col, lwd = 1.15, cex = 1.05
  )
  mtext(panel_tag, side = 3, adj = 0, line = 0.35, font = 2, cex = 0.95)
}

draw_figure <- function(mode = c("combined", "delta", "detection")) {
  mode <- match.arg(mode)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if (mode == "combined") {
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
  } else {
    layout(matrix(1:2, nrow = 1, byrow = TRUE))
  }
  par(
    mar = c(3.5, 4.2, 2.2, 0.9),
    oma = c(2.2, 0.5, 0.2, 0.2),
    mgp = c(2.25, 0.60, 0),
    tcl = -0.25,
    family = "sans",
    las = 1,
    xaxs = "r",
    yaxs = "r"
  )

  if (mode %in% c("combined", "delta")) {
    draw_delta_panel(panel_order[1], "a")
    draw_delta_panel(panel_order[2], "b")
  }
  if (mode == "combined") {
    draw_detection_panel(panel_order[1], "c")
    draw_detection_panel(panel_order[2], "d")
  } else if (mode == "detection") {
    draw_detection_panel(panel_order[1], "a")
    draw_detection_panel(panel_order[2], "b")
  }

  mtext("Oracle-label accuracy", side = 1, outer = TRUE, line = 0.55, cex = 0.82)
}

width_mm <- 183
width_in <- width_mm / 25.4

export_figure <- function(stem, mode, height_mm) {
  height_in <- height_mm / 25.4
  pdf_file <- file.path(output_figure_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_figure_dir, paste0(stem, ".png"))
  cairo_pdf(pdf_file, width = width_in, height = height_in, family = "sans", onefile = TRUE)
  draw_figure(mode)
  dev.off()
  png(
    png_file, width = width_in, height = height_in, units = "in",
    res = 600, type = "cairo", bg = "white"
  )
  draw_figure(mode)
  dev.off()
  c(pdf = pdf_file, png = png_file)
}

delta_files <- export_figure(
  "Figure_Sensitivity02_oracle_label_delta_auc_point_range_v2",
  "delta", 94
)
detection_files <- export_figure(
  "Figure_Sensitivity02_oracle_label_detection_rate_point_range_v2",
  "detection", 94
)
combined_files <- export_figure(
  "Figure_Sensitivity02_oracle_label_delta_auc_and_detection_v2",
  "combined", 154
)

figure_legend <- paste(
  "Sensitivity analysis 02: metric-level positive control for Domain 2.",
  "Panels a-b show the across-repeat mean change in AUC after adding the noisy oracle label",
  "to the raw-feature elastic-net model; error bars are 95% Monte Carlo confidence intervals",
  "for the across-repeat mean, not empirical repeat distributions or within-study Rubin intervals.",
  "The gray band denotes the prespecified practical-null region (absolute Delta AUC < 0.001).",
  "Panels c-d show the proportion of simulated studies whose repeat-level Rubin 95% confidence",
  "interval was entirely above zero; error bars are Wilson 95% binomial intervals.",
  "Horizontal reference lines in panels c-d mark detection rates of 80% and 90%.",
  "The oracle label was generated by symmetrically corrupting the latent binary state at the",
  "prespecified accuracy and was held fixed across the five imputed datasets.",
  "Each point summarizes 100 Monte Carlo repeats with five imputations per repeat."
)
writeLines(
  figure_legend,
  file.path(output_figure_dir, "Figure_Sensitivity02_oracle_label_legend_v2.txt")
)

monotonicity <- do.call(rbind, lapply(panel_order, function(label) {
  d <- setting_summary[setting_summary$outcome_label == label, , drop = FALSE]
  d <- d[order(d$oracle_accuracy), , drop = FALSE]
  data.frame(
    outcome_label = label,
    spearman_rho = cor(d$oracle_accuracy, d$mean_delta_auc, method = "spearman"),
    strictly_increasing = all(diff(d$mean_delta_auc) > 0),
    stringsAsFactors = FALSE
  )
}))

qc_lines <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Source Rubin rows: ", nrow(rubin)),
  paste0("Manifest rows: ", nrow(manifest_map)),
  paste0("Joined repeat-level rows: ", nrow(repeat_level)),
  paste0("Corrected setting-level rows: ", nrow(setting_summary)),
  paste0("Practical-null threshold: ", practical_null_threshold),
  paste0("All settings have 100 repeats: ", all(setting_summary$n_repeats == 100L)),
  paste0(
    "All repeat estimates used m=5: ",
    all(setting_summary$m_used_min == 5L & setting_summary$m_used_max == 5L)
  ),
  paste0("Oracle fidelity repeat-level rows: ", nrow(oracle_qc_repeat)),
  paste0("Oracle fidelity setting-level rows: ", nrow(oracle_qc_summary)),
  paste0(
    "Maximum absolute mean target-realized difference: ",
    sprintf("%.6f", max(abs(oracle_qc_summary$target_minus_realized_accuracy_mean)))
  ),
  "Original aggregate and mean-CI flag files are deprecated because of grouping errors.",
  "Monotonicity by outcome:",
  paste0(
    "  ", monotonicity$outcome_label,
    ": Spearman rho=", sprintf("%.3f", monotonicity$spearman_rho),
    "; strictly increasing=", monotonicity$strictly_increasing
  ),
  paste0("Summary: ", summary_file),
  paste0("Delta-AUC figure PDF: ", delta_files["pdf"]),
  paste0("Detection-rate figure PDF: ", detection_files["pdf"]),
  paste0("Combined figure PDF: ", combined_files["pdf"]),
  paste0("Combined figure PNG: ", combined_files["png"])
)
writeLines(qc_lines, file.path(output_log_dir, "sensitivity02_correction_qc.txt"), useBytes = TRUE)

cat(paste(qc_lines, collapse = "\n"), "\n")
