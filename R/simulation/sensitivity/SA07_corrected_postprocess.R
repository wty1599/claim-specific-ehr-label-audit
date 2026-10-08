#!/usr/bin/env Rscript

# Corrected post-processing for Sensitivity 07 and 07b.
# Base R only. Source results are read-only; all outputs go to new folders.

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

base_dir <- Sys.getenv("SIMULATION_RESULTS_DIR", unset = get_script_dir())
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)

dir07 <- file.path(base_dir, "SA07_clustering_gmm_ward_scan")
dir07b <- file.path(base_dir, "SA07b_gmm_ward_transport")

if (!dir.exists(dir07)) stop("Sensitivity 07 directory not found: ", dir07, call. = FALSE)
if (!dir.exists(dir07b)) stop("Sensitivity 07b directory not found: ", dir07b, call. = FALSE)

out07 <- file.path(dir07, "output", "postprocessing_corrected")
out07b <- file.path(dir07b, "output", "postprocessing_corrected")
for (p in c(
  file.path(out07, "tables"), file.path(out07, "logs"),
  file.path(out07b, "tables"), file.path(out07b, "logs")
)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"), check.names = FALSE)
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_sd <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE)
safe_quantile <- function(x, p) {
  if (!any(!is.na(x))) return(NA_real_)
  unname(quantile(x, probs = p, na.rm = TRUE, type = 8))
}

metric_bounds <- function(metric) {
  m <- tolower(metric)
  if (grepl("ari", m, fixed = TRUE)) return(c(-1, 1))
  if (grepl("abs.*prevalence|prevalence.*abs", m)) return(c(0, 1))
  if (grepl("prevalence.*drift|signed.*prevalence|prevalence.*difference", m)) {
    return(c(-1, 1))
  }
  if (grepl("agreement|proportion|auc$|dip.*p|identifiable|rate", m)) {
    return(c(0, 1))
  }
  c(-Inf, Inf)
}

clip_metric <- function(x, metric) {
  b <- metric_bounds(metric)
  pmin(b[2], pmax(b[1], x))
}

summarize_values <- function(x, metric) {
  n <- sum(!is.na(x))
  s <- safe_sd(x)
  data.frame(
    nonmissing_n = n,
    mean_estimate = safe_mean(x),
    sd_across_repeats = s,
    empirical_p025 = clip_metric(safe_quantile(x, 0.025), metric),
    empirical_p975 = clip_metric(safe_quantile(x, 0.975), metric),
    mcse_mean = if (n >= 2L) s / sqrt(n) else NA_real_,
    stringsAsFactors = FALSE
  )
}

summarize_long <- function(data, group_columns, value_column = "estimate") {
  group_key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), group_key)
  rows <- lapply(groups, function(idx) {
    d <- data[idx, , drop = FALSE]
    metric <- if ("metric" %in% names(d)) d$metric[1] else "estimate"
    cbind(
      d[1, group_columns, drop = FALSE],
      n_repeats = length(unique(d$repeat_id)),
      summarize_values(as.numeric(d[[value_column]]), metric)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

wilson_interval <- function(successes, n, conf = 0.95) {
  if (is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- successes / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  pmax(0, pmin(1, c(center - half, center + half)))
}

# -------------------------------------------------------------------------
# Sensitivity 07: corrected aggregate and MC-SE
# -------------------------------------------------------------------------

tab07 <- file.path(dir07, "output", "tables")
master07 <- read_csv(file.path(tab07, "simulation_four_domain_master_summary.csv"))
grid07 <- read_csv(file.path(tab07, "sensitivity_grid.csv"))
centroids07 <- read_csv(file.path(tab07, "primary_imputation1_centroids.csv"))

algorithm_map07 <- setNames(
  tolower(grid07$CLUSTER_ALGORITHM_value),
  as.character(grid07$sens_setting_index)
)
id_map07 <- setNames(grid07$sensitivity_id, as.character(grid07$sens_setting_index))
master07$sensitivity_id <- unname(id_map07[as.character(master07$sens_sens_setting_index)])
master07$requested_algorithm <- unname(
  algorithm_map07[as.character(master07$sens_sens_setting_index)]
)
status_by_id07 <- tapply(
  centroids07$cluster_algorithm_status,
  centroids07$sensitivity_id,
  function(x) all(x == "ok")
)
master07$algorithm_valid <- unname(status_by_id07[master07$sensitivity_id])
master07$algorithm <- ifelse(
  master07$algorithm_valid,
  master07$requested_algorithm,
  "kmeans_fallback"
)
if (anyNA(master07$requested_algorithm) || anyNA(master07$algorithm_valid)) {
  stop("Could not map requested/effective algorithm in Sensitivity 07.", call. = FALSE)
}
master07$eligible_for_formal_inference <- master07$algorithm_valid

# All original 07 labels are K2 except the explicit true-K K3 ARI diagnostic.
master07$cluster_k <- ifelse(
  grepl("ari_k3|k3_", master07$metric, ignore.case = TRUE),
  3L,
  2L
)

aggregate07 <- summarize_long(
  master07,
  c(
    "scenario", "requested_algorithm", "algorithm", "algorithm_valid",
    "eligible_for_formal_inference", "cluster_k", "domain", "metric"
  )
)
aggregate07 <- aggregate07[
  order(aggregate07$scenario, aggregate07$algorithm, aggregate07$cluster_k,
        aggregate07$domain, aggregate07$metric),
  , drop = FALSE
]

formal_metric <- with(master07,
  (domain == "Domain1_Discreteness" & metric %in% c(
    "dip_discriminant_p", "sigclust_effect_sd",
    "ari_to_true_subtype", "ari_k3_to_true_subtype"
  )) |
  (domain == "Domain2_IncrementalValue" & grepl("decisive_delta_auc", metric)) |
  (domain == "Domain3_Transportability" & metric ==
    "C1_prevalence_drift_validation_minus_derivation")
)
formal07 <- summarize_long(
  master07[formal_metric, , drop = FALSE],
  c(
    "scenario", "requested_algorithm", "algorithm", "algorithm_valid",
    "eligible_for_formal_inference", "cluster_k", "domain", "metric"
  )
)
formal07 <- formal07[
  order(formal07$scenario, formal07$algorithm, formal07$cluster_k,
        formal07$domain, formal07$metric),
  , drop = FALSE
]

# MICE ARI is retained as a separate formal metric family.
mice07 <- read_csv(file.path(tab07, "mice_label_stability_ari_summary.csv"))
mice07$sensitivity_id <- unname(id_map07[as.character(mice07$sens_sens_setting_index)])
mice07$requested_algorithm <- unname(
  algorithm_map07[as.character(mice07$sens_sens_setting_index)]
)
mice07$algorithm_valid <- unname(status_by_id07[mice07$sensitivity_id])
mice07$algorithm <- ifelse(
  mice07$algorithm_valid,
  mice07$requested_algorithm,
  "kmeans_fallback"
)
mice_metric_columns <- c(
  "ari_mean_vs_primary", "ari_min_vs_primary", "ari_prop_ge_0_90"
)
mice_long <- do.call(rbind, lapply(mice_metric_columns, function(metric) {
  data.frame(
    scenario = mice07$scenario,
    repeat_id = mice07$repeat_id,
    requested_algorithm = mice07$requested_algorithm,
    algorithm = mice07$algorithm,
    algorithm_valid = mice07$algorithm_valid,
    eligible_for_formal_inference = mice07$algorithm_valid,
    cluster_k = 2L,
    domain = "MICE_LabelStability",
    metric = metric,
    estimate = as.numeric(mice07[[metric]]),
    stringsAsFactors = FALSE
  )
}))
mice_summary07 <- summarize_long(
  mice_long,
  c(
    "scenario", "requested_algorithm", "algorithm", "algorithm_valid",
    "eligible_for_formal_inference", "cluster_k", "domain", "metric"
  )
)

# Audit whether requested algorithms actually ran.
audit_key07 <- interaction(
  centroids07$sensitivity_id,
  centroids07$sens_CLUSTER_ALGORITHM_value,
  centroids07$cluster_algorithm_status,
  drop = TRUE,
  lex.order = TRUE
)
algorithm_audit07 <- do.call(rbind, lapply(split(seq_len(nrow(centroids07)), audit_key07), function(idx) {
  d <- centroids07[idx, , drop = FALSE]
  data.frame(
    sensitivity_id = d$sensitivity_id[1],
    requested_algorithm = d$sens_CLUSTER_ALGORITHM_value[1],
    recorded_algorithm = d$cluster_algorithm[1],
    execution_status = d$cluster_algorithm_status[1],
    n_rows = nrow(d),
    n_repeats = length(unique(d$repeat_id)),
    n_scenarios = length(unique(d$scenario)),
    valid_requested_algorithm = identical(d$cluster_algorithm_status[1], "ok"),
    stringsAsFactors = FALSE
  )
}))
rownames(algorithm_audit07) <- NULL

write.csv(
  aggregate07,
  file.path(out07, "tables", "sensitivity07_aggregate_summary_corrected.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  formal07,
  file.path(out07, "tables", "sensitivity07_formal_metrics_summary_corrected.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  mice_summary07,
  file.path(out07, "tables", "sensitivity07_mice_ari_summary_corrected.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  algorithm_audit07,
  file.path(out07, "tables", "sensitivity07_algorithm_execution_audit.csv"),
  row.names = FALSE, na = ""
)

# -------------------------------------------------------------------------
# Sensitivity 07b: empirical intervals and corrected operating definitions
# -------------------------------------------------------------------------

tab07b <- file.path(dir07b, "output", "tables")
assign07b <- read_csv(file.path(tab07b, "domain3_assignment_transportability_summary.csv"))
grid07b <- read_csv(file.path(tab07b, "sensitivity_grid.csv"))
centroids07b <- read_csv(file.path(tab07b, "primary_imputation1_centroids.csv"))
algorithm_map07b <- setNames(
  tolower(grid07b$CLUSTER_ALGORITHM_value),
  grid07b$sensitivity_id
)
assign07b$requested_algorithm <- unname(algorithm_map07b[assign07b$sensitivity_id])
status_by_id07b <- tapply(
  centroids07b$cluster_algorithm_status,
  centroids07b$sensitivity_id,
  function(x) all(x == "ok")
)
assign07b$algorithm_valid <- unname(status_by_id07b[assign07b$sensitivity_id])
assign07b$algorithm <- ifelse(
  assign07b$algorithm_valid,
  assign07b$requested_algorithm,
  "kmeans_fallback"
)
assign07b$eligible_for_formal_inference <- assign07b$algorithm_valid
if (anyNA(assign07b$requested_algorithm) || anyNA(assign07b$algorithm_valid)) {
  stop("Could not map requested/effective algorithm in Sensitivity 07b.", call. = FALSE)
}

transport_metrics <- c(
  frozen_vs_reclustered_ari = "assignment_ari",
  frozen_vs_reclustered_matched_agreement = "matched_agreement",
  mean_abs_cluster_prevalence_difference = "mean_abs_prevalence_difference",
  max_abs_cluster_prevalence_difference = "max_abs_prevalence_difference"
)

transport_long <- do.call(rbind, lapply(names(transport_metrics), function(column) {
  data.frame(
    sensitivity_id = assign07b$sensitivity_id,
    requested_algorithm = assign07b$requested_algorithm,
    algorithm = assign07b$algorithm,
    algorithm_valid = assign07b$algorithm_valid,
    eligible_for_formal_inference = assign07b$eligible_for_formal_inference,
    scenario = assign07b$scenario,
    repeat_id = assign07b$repeat_id,
    cluster_k = as.integer(assign07b$cluster_k),
    domain = "Domain3_AssignmentTransportability",
    metric = unname(transport_metrics[column]),
    estimate = as.numeric(assign07b[[column]]),
    stringsAsFactors = FALSE
  )
}))

transport_summary07b <- summarize_long(
  transport_long,
  c(
    "sensitivity_id", "requested_algorithm", "algorithm", "algorithm_valid",
    "eligible_for_formal_inference", "scenario", "cluster_k", "domain", "metric"
  )
)
transport_summary07b <- transport_summary07b[
  order(transport_summary07b$algorithm, transport_summary07b$scenario,
        transport_summary07b$cluster_k, transport_summary07b$metric),
  , drop = FALSE
]

# Prespecified repeat-level assignment failure rule.
assign07b$assignment_failure <- with(assign07b,
  frozen_vs_reclustered_ari < 0.80 |
  frozen_vs_reclustered_matched_agreement < 0.90 |
  max_abs_cluster_prevalence_difference > 0.10
)
assign07b$evaluation_role <- ifelse(
  assign07b$scenario == "S1_pure_severity_continuum",
  "partition_reproducibility_stress_test",
  ifelse(
    assign07b$scenario == "S2_true_discrete_subtypes",
    "expected_assignment_success",
    "expected_assignment_failure"
  )
)
assign07b$included_in_sensitivity_specificity <-
  assign07b$scenario %in% c("S2_true_discrete_subtypes", "S3_transport_drift")
assign07b$correct_classification <- ifelse(
  assign07b$scenario == "S2_true_discrete_subtypes",
  !assign07b$assignment_failure,
  ifelse(
    assign07b$scenario == "S3_transport_drift",
    assign07b$assignment_failure,
    NA
  )
)

op_key <- interaction(
  assign07b$sensitivity_id, assign07b$requested_algorithm, assign07b$algorithm,
  assign07b$algorithm_valid, assign07b$scenario,
  assign07b$cluster_k, drop = TRUE, lex.order = TRUE
)
operating07b <- do.call(rbind, lapply(split(seq_len(nrow(assign07b)), op_key), function(idx) {
  d <- assign07b[idx, , drop = FALSE]
  n <- nrow(d)
  failures <- sum(d$assignment_failure)
  failure_ci <- wilson_interval(failures, n)
  correct <- d$correct_classification[!is.na(d$correct_classification)]
  correct_ci <- if (length(correct)) wilson_interval(sum(correct), length(correct)) else c(NA, NA)
  data.frame(
    sensitivity_id = d$sensitivity_id[1],
    requested_algorithm = d$requested_algorithm[1],
    algorithm = d$algorithm[1],
    algorithm_valid = d$algorithm_valid[1],
    eligible_for_formal_inference = d$eligible_for_formal_inference[1],
    scenario = d$scenario[1],
    cluster_k = d$cluster_k[1],
    evaluation_role = d$evaluation_role[1],
    included_in_sensitivity_specificity = d$included_in_sensitivity_specificity[1],
    n_repeats = n,
    assignment_failure_rate = failures / n,
    assignment_failure_ci_low = failure_ci[1],
    assignment_failure_ci_high = failure_ci[2],
    assignment_success_rate = 1 - failures / n,
    correct_classification_rate = if (length(correct)) mean(correct) else NA_real_,
    correct_classification_ci_low = correct_ci[1],
    correct_classification_ci_high = correct_ci[2],
    stringsAsFactors = FALSE
  )
}))
rownames(operating07b) <- NULL
operating07b <- operating07b[
  order(operating07b$algorithm, operating07b$scenario, operating07b$cluster_k),
  , drop = FALSE
]

ss_key <- interaction(
  assign07b$requested_algorithm, assign07b$algorithm,
  assign07b$algorithm_valid, assign07b$cluster_k,
  drop = TRUE, lex.order = TRUE
)
sensitivity_specificity07b <- do.call(rbind, lapply(split(seq_len(nrow(assign07b)), ss_key), function(idx) {
  d <- assign07b[idx, , drop = FALSE]
  s2 <- d[d$scenario == "S2_true_discrete_subtypes", , drop = FALSE]
  s3 <- d[d$scenario == "S3_transport_drift", , drop = FALSE]
  sensitivity <- mean(s3$assignment_failure)
  specificity <- mean(!s2$assignment_failure)
  sens_ci <- wilson_interval(sum(s3$assignment_failure), nrow(s3))
  spec_ci <- wilson_interval(sum(!s2$assignment_failure), nrow(s2))
  data.frame(
    requested_algorithm = d$requested_algorithm[1],
    algorithm = d$algorithm[1],
    algorithm_valid = d$algorithm_valid[1],
    eligible_for_formal_inference = d$algorithm_valid[1],
    cluster_k = d$cluster_k[1],
    sensitivity_for_assignment_failure = sensitivity,
    sensitivity_ci_low = sens_ci[1],
    sensitivity_ci_high = sens_ci[2],
    specificity_for_assignment_success = specificity,
    specificity_ci_low = spec_ci[1],
    specificity_ci_high = spec_ci[2],
    balanced_accuracy = mean(c(sensitivity, specificity)),
    n_s2 = nrow(s2),
    n_s3 = nrow(s3),
    s1_excluded = TRUE,
    stringsAsFactors = FALSE
  )
}))
rownames(sensitivity_specificity07b) <- NULL

audit_key07b <- interaction(
  centroids07b$sensitivity_id,
  centroids07b$sens_CLUSTER_ALGORITHM_value,
  centroids07b$cluster_algorithm_status,
  drop = TRUE,
  lex.order = TRUE
)
algorithm_audit07b <- do.call(rbind, lapply(split(seq_len(nrow(centroids07b)), audit_key07b), function(idx) {
  d <- centroids07b[idx, , drop = FALSE]
  data.frame(
    sensitivity_id = d$sensitivity_id[1],
    requested_algorithm = d$sens_CLUSTER_ALGORITHM_value[1],
    recorded_algorithm = d$cluster_algorithm[1],
    execution_status = d$cluster_algorithm_status[1],
    n_rows = nrow(d),
    n_repeats = length(unique(d$repeat_id)),
    n_scenarios = length(unique(d$scenario)),
    valid_requested_algorithm = identical(d$cluster_algorithm_status[1], "ok"),
    stringsAsFactors = FALSE
  )
}))
rownames(algorithm_audit07b) <- NULL

write.csv(
  transport_summary07b,
  file.path(out07b, "tables", "sensitivity07b_transport_metrics_empirical_intervals.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  operating07b,
  file.path(out07b, "tables", "sensitivity07b_operating_characteristics_corrected.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  sensitivity_specificity07b,
  file.path(out07b, "tables", "sensitivity07b_sensitivity_specificity_corrected.csv"),
  row.names = FALSE, na = ""
)
write.csv(
  algorithm_audit07b,
  file.path(out07b, "tables", "sensitivity07b_algorithm_execution_audit.csv"),
  row.names = FALSE, na = ""
)

gmm_valid07 <- any(
  algorithm_audit07$requested_algorithm == "gmm" &
  algorithm_audit07$valid_requested_algorithm
)
gmm_valid07b <- any(
  algorithm_audit07b$requested_algorithm == "gmm" &
  algorithm_audit07b$valid_requested_algorithm
)

direct_status <- data.frame(
  requested_analysis = "GMM_vs_Ward_direct_label_agreement",
  individual_labels_available = FALSE,
  gmm_arm_valid_in_07 = gmm_valid07,
  gmm_arm_valid_in_07b = gmm_valid07b,
  postprocessing_possible = FALSE,
  required_action = paste(
    "Lightweight targeted rerun of GMM and Ward labels on paired datasets;",
    "use mclust::Mclust inside each multisession worker and disable k-means fallback."
  ),
  reason = paste(
    "Individual labels were not saved, and the recorded GMM arm fell back to k-means",
    "after mclustBIC was unavailable in workers."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  direct_status,
  file.path(out07b, "tables", "sensitivity07_direct_algorithm_agreement_status.csv"),
  row.names = FALSE, na = ""
)

targeted_rerun_spec <- c(
  "Sensitivity 07 targeted paired-label rerun specification",
  "",
  "Scope: S1, S2, and S3; K=2 and K=3; GMM and Ward on the same generated datasets.",
  "Do not rerun Domain 2, Domain 3 outcome models, or Domain 4.",
  "Use the same primary completed dataset and frozen preprocessing within each repeat.",
  "Call mclust::Mclust explicitly inside each multisession worker.",
  "Declare future.packages including mclust when using future_lapply.",
  "Treat any GMM error as a failed repeat; do not fall back to k-means.",
  "Save one long CSV with columns:",
  "scenario, repeat_id, cluster_k, patient_id, algorithm, cluster, algorithm_status",
  "Optionally append standardized features prefixed feature__ for cluster-profile correlations.",
  "Expected filename: sensitivity07_targeted_pair_labels.csv",
  "Then run postprocess_sensitivity07_direct_algorithm_agreement.R."
)
writeLines(
  targeted_rerun_spec,
  file.path(out07b, "logs", "sensitivity07_targeted_paired_label_rerun_spec.txt")
)

qc07 <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Original 07 master rows: ", nrow(master07)),
  paste0("Corrected 07 aggregate rows: ", nrow(aggregate07)),
  paste0("Formal 07 summary rows: ", nrow(formal07)),
  paste0("All nonmissing formal groups have n <= 50: ", all(formal07$nonmissing_n <= 50)),
  paste0("GMM requested algorithm valid in 07: ", gmm_valid07),
  "SigClust-like p-values are excluded; only effect size is retained."
)
writeLines(qc07, file.path(out07, "logs", "sensitivity07_postprocessing_qc.txt"))

qc07b <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("07b repeat-level assignment rows: ", nrow(assign07b)),
  paste0("07b empirical interval rows: ", nrow(transport_summary07b)),
  paste0("07b corrected operating rows: ", nrow(operating07b)),
  "S1 is a partition reproducibility stress test and is excluded from sensitivity/specificity.",
  "S2 defines expected assignment success; S3 defines expected assignment failure.",
  paste0("GMM requested algorithm valid in 07b: ", gmm_valid07b),
  "ARI intervals are bounded to [-1, 1]; agreement/proportion intervals are bounded to [0, 1]."
)
writeLines(qc07b, file.path(out07b, "logs", "sensitivity07b_postprocessing_qc.txt"))

cat(paste(c(qc07, "", qc07b), collapse = "\n"), "\n")
