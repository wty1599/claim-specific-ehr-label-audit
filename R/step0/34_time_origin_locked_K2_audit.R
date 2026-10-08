# Time-origin audit for the locked MIMIC-IV K=2 partition.
# This is a post hoc semantic-composition analysis and does not regenerate labels.

options(stringsAsFactors = FALSE)

work_root <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT", unset = "")
if (!nzchar(work_root) || !nzchar(output_root)) {
  stop("Set EHR_AUDIT_WORK_ROOT and EHR_AUDIT_OUTPUT_ROOT")
}
input_data <- file.path(work_root, "data", "final_full.csv")
input_labels <- file.path(work_root, "output", "model", "labels_primary_mice.rds")
output_dir <- file.path(output_root, "time_origin_audit")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_data), file.exists(input_labels))

dat <- read.csv(
  input_data,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(stay_id = "numeric")
)
labels <- as.data.frame(readRDS(input_labels))

required_data <- c("stay_id", "icu_intime", "sepsis_onset_time")
required_labels <- c("stay_id", "cluster_k2")
stopifnot(
  all(required_data %in% names(dat)),
  all(required_labels %in% names(labels)),
  nrow(dat) == 20049L,
  nrow(labels) == 20049L,
  !anyDuplicated(dat$stay_id),
  !anyDuplicated(labels$stay_id)
)

audit <- merge(
  dat[required_data],
  labels[required_labels],
  by = "stay_id",
  all = FALSE,
  sort = FALSE
)
stopifnot(nrow(audit) == 20049L, !anyNA(audit$cluster_k2))

parse_datetime <- function(x) {
  as.POSIXct(
    x,
    tz = "UTC",
    tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%OS")
  )
}

audit$icu_time <- parse_datetime(audit$icu_intime)
audit$onset_time <- parse_datetime(audit$sepsis_onset_time)
stopifnot(!anyNA(audit$icu_time), !anyNA(audit$onset_time))

audit$onset_minus_icu_hours <- as.numeric(
  difftime(audit$onset_time, audit$icu_time, units = "hours")
)
audit$time_origin <- ifelse(
  audit$onset_minus_icu_hours < 0,
  "Before ICU admission",
  ifelse(
    audit$onset_minus_icu_hours > 0,
    "After ICU admission",
    "Same timestamp"
  )
)
audit$phenotype <- ifelse(
  audit$cluster_k2 == 1,
  "C1 higher-risk",
  "C2 lower-risk"
)

stopifnot(
  sum(audit$time_origin == "Before ICU admission") == 11025L,
  sum(audit$time_origin == "After ICU admission") == 9020L,
  sum(audit$time_origin == "Same timestamp") == 4L,
  sum(audit$phenotype == "C1 higher-risk") == 3992L,
  sum(audit$phenotype == "C2 lower-risk") == 16057L
)

wilson_ci <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(low = center - half, high = center + half)
}

time_levels <- c(
  "Before ICU admission",
  "After ICU admission",
  "Same timestamp"
)

crosstab_rows <- lapply(time_levels, function(level) {
  idx <- audit$time_origin == level
  n_total <- sum(idx)
  n_c1 <- sum(idx & audit$cluster_k2 == 1)
  n_c2 <- sum(idx & audit$cluster_k2 == 2)
  ci <- wilson_ci(n_c1, n_total)
  data.frame(
    time_origin = level,
    C1_n = n_c1,
    C2_n = n_c2,
    total_n = n_total,
    C1_prevalence = n_c1 / n_total,
    C1_prevalence_wilson_low = ci["low"],
    C1_prevalence_wilson_high = ci["high"],
    stringsAsFactors = FALSE
  )
})
crosstab <- do.call(rbind, crosstab_rows)

total_ci <- wilson_ci(sum(audit$cluster_k2 == 1), nrow(audit))
crosstab <- rbind(
  crosstab,
  data.frame(
    time_origin = "Overall",
    C1_n = sum(audit$cluster_k2 == 1),
    C2_n = sum(audit$cluster_k2 == 2),
    total_n = nrow(audit),
    C1_prevalence = mean(audit$cluster_k2 == 1),
    C1_prevalence_wilson_low = total_ci["low"],
    C1_prevalence_wilson_high = total_ci["high"],
    stringsAsFactors = FALSE
  )
)
rownames(crosstab) <- NULL

two_by_two <- subset(audit, time_origin != "Same timestamp")
tab <- table(
  before_icu = two_by_two$time_origin == "Before ICU admission",
  C1 = two_by_two$cluster_k2 == 1
)

a <- unname(tab["TRUE", "TRUE"])
b <- unname(tab["TRUE", "FALSE"])
c_count <- unname(tab["FALSE", "TRUE"])
d <- unname(tab["FALSE", "FALSE"])

p_before <- a / (a + b)
p_after <- c_count / (c_count + d)
risk_ratio <- p_before / p_after
se_log_rr <- sqrt(
  1 / a - 1 / (a + b) +
    1 / c_count - 1 / (c_count + d)
)
rr_ci <- exp(log(risk_ratio) + c(-1, 1) * qnorm(0.975) * se_log_rr)

odds_ratio <- (a * d) / (b * c_count)
se_log_or <- sqrt(1 / a + 1 / b + 1 / c_count + 1 / d)
or_ci <- exp(log(odds_ratio) + c(-1, 1) * qnorm(0.975) * se_log_or)

risk_difference <- p_before - p_after
se_rd <- sqrt(
  p_before * (1 - p_before) / (a + b) +
    p_after * (1 - p_after) / (c_count + d)
)
rd_ci <- risk_difference + c(-1, 1) * qnorm(0.975) * se_rd

chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
cramers_v <- sqrt(unname(chi$statistic) / sum(tab))

rank_auc <- function(outcome, score) {
  outcome <- as.logical(outcome)
  ranks <- rank(score, ties.method = "average")
  n_positive <- sum(outcome)
  n_negative <- sum(!outcome)
  (
    sum(ranks[outcome]) - n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

binary_time_origin_auc <- rank_auc(
  audit$cluster_k2 == 1,
  as.numeric(audit$time_origin == "Before ICU admission")
)
continuous_onset_offset_auc <- rank_auc(
  audit$cluster_k2 == 1,
  -audit$onset_minus_icu_hours
)

effect_estimates <- data.frame(
  estimand = c(
    "C1 prevalence before ICU onset",
    "C1 prevalence after ICU onset",
    "C1 prevalence risk difference: before minus after",
    "C1 prevalence risk ratio: before versus after",
    "C1 prevalence odds ratio: before versus after",
    "Pearson chi-square",
    "Cramer's V",
    "AUC for C1 using before-versus-after onset only",
    "AUC for C1 using earlier continuous onset offset"
  ),
  estimate = c(
    p_before,
    p_after,
    risk_difference,
    risk_ratio,
    odds_ratio,
    unname(chi$statistic),
    cramers_v,
    binary_time_origin_auc,
    continuous_onset_offset_auc
  ),
  ci_low = c(NA, NA, rd_ci[1], rr_ci[1], or_ci[1], NA, NA, NA, NA),
  ci_high = c(NA, NA, rd_ci[2], rr_ci[2], or_ci[2], NA, NA, NA, NA),
  p_value = c(NA, NA, NA, NA, NA, chi$p.value, NA, NA, NA),
  stringsAsFactors = FALSE
)

offset_summary <- do.call(
  rbind,
  lapply(c("C1 higher-risk", "C2 lower-risk"), function(group) {
    x <- audit$onset_minus_icu_hours[audit$phenotype == group]
    data.frame(
      phenotype = group,
      n = length(x),
      median_onset_minus_icu_hours = median(x),
      q1_onset_minus_icu_hours = unname(quantile(x, 0.25)),
      q3_onset_minus_icu_hours = unname(quantile(x, 0.75)),
      before_icu_n = sum(x < 0),
      before_icu_proportion = mean(x < 0),
      stringsAsFactors = FALSE
    )
  })
)

audit$onset_offset_bin <- cut(
  audit$onset_minus_icu_hours,
  breaks = c(-Inf, -12, -6, 0, 6, 12, Inf),
  right = FALSE,
  labels = c(
    "<-12 h",
    "-12 to <-6 h",
    "-6 to <0 h",
    "0 to <6 h",
    "6 to <12 h",
    ">=12 h"
  )
)

offset_bin_summary <- do.call(
  rbind,
  lapply(levels(audit$onset_offset_bin), function(level) {
    idx <- audit$onset_offset_bin == level
    n_total <- sum(idx)
    n_c1 <- sum(idx & audit$cluster_k2 == 1)
    ci <- wilson_ci(n_c1, n_total)
    data.frame(
      onset_minus_icu_bin = level,
      C1_n = n_c1,
      C2_n = n_total - n_c1,
      total_n = n_total,
      C1_prevalence = n_c1 / n_total,
      C1_prevalence_wilson_low = ci["low"],
      C1_prevalence_wilson_high = ci["high"],
      stringsAsFactors = FALSE
    )
  })
)
rownames(offset_bin_summary) <- NULL

hash_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  unname(tools::md5sum(path))
}

provenance <- data.frame(
  item = c(
    "analysis_status",
    "analysis_scope",
    "input_data",
    "input_data_hash",
    "input_labels",
    "input_labels_hash",
    "label_definition",
    "same_timestamp_handling",
    "R_version"
  ),
  value = c(
    "post hoc completed; not yet incorporated into the manuscript",
    "full locked MIMIC-IV cohort; no risk-set restriction",
    input_data,
    hash_file(input_data),
    input_labels,
    hash_file(input_labels),
    "cluster_k2=1 is C1 higher-risk; cluster_k2=2 is C2 lower-risk",
    "reported in the descriptive table; excluded from 2x2 effect estimates",
    R.version.string
  ),
  stringsAsFactors = FALSE
)

write.csv(
  crosstab,
  file.path(output_dir, "Table_time_origin_by_locked_K2.csv"),
  row.names = FALSE
)
write.csv(
  effect_estimates,
  file.path(output_dir, "Table_time_origin_K2_effect_estimates.csv"),
  row.names = FALSE
)
write.csv(
  offset_summary,
  file.path(output_dir, "Table_onset_offset_hours_by_locked_K2.csv"),
  row.names = FALSE
)
write.csv(
  offset_bin_summary,
  file.path(output_dir, "Table_onset_offset_bins_by_locked_K2.csv"),
  row.names = FALSE
)
write.csv(
  provenance,
  file.path(output_dir, "time_origin_audit_provenance.csv"),
  row.names = FALSE
)

cat("Time-origin audit completed.\n")
print(crosstab)
print(effect_estimates)
print(offset_summary)
print(offset_bin_summary)
cat("Output:", output_dir, "\n")
