## =====================================================================
## 29a_domain1_aki_definition_alignment.R
## KDIGO alignment for the locked empirical K2 labels.
##
## The project-authoritative 0-24 h export contains the combined KDIGO
## stage but not separate creatinine and urine-output criterion fields.
## This script therefore does not reconstruct creatinine-only, urine-only,
## both, or neither groups. Their absence is reported explicitly.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "ggplot2", "digest")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
  library(ggplot2)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), winslash = "/", mustWork = TRUE)
OUT_ROOT <- normalizePath(file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"), winslash = "/", mustWork = TRUE)
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"), figures = file.path(OUT_ROOT, "figures"),
  logs = file.path(OUT_ROOT, "logs"), provenance = file.path(OUT_ROOT, "provenance")
)
INPUT <- list(
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  baseline = file.path(PROJECT_ROOT, "data", "baseline_covars.csv"),
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  kdigo_sql = file.path(PROJECT_ROOT, "sql", "build_baseline_covars_v2.sql")
)
if (any(!file.exists(unlist(INPUT)))) stop("Missing input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}
write_lines_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)

binary_rule_metrics <- function(stage, c1, rule, label, window) {
  pred <- as.integer(rule(stage))
  tp <- sum(pred == 1L & c1 == 1L)
  tn <- sum(pred == 0L & c1 == 0L)
  fp <- sum(pred == 1L & c1 == 0L)
  fn <- sum(pred == 0L & c1 == 1L)
  data.table(
    time_window = window, rule = label, n = length(c1),
    tp = tp, tn = tn, fp = fp, fn = fn,
    sensitivity_P_rule_given_C1 = safe_div(tp, tp + fn),
    specificity_P_not_rule_given_C2 = safe_div(tn, tn + fp),
    ppv_P_C1_given_rule = safe_div(tp, tp + fp),
    npv_P_C2_given_not_rule = safe_div(tn, tn + fn)
  )
}

or_from_glm <- function(fit, term) {
  beta <- unname(coef(fit)[term])
  se <- sqrt(diag(vcov(fit)))[term]
  data.table(or = exp(beta), or_ci_low = exp(beta - 1.96 * se), or_ci_high = exp(beta + 1.96 * se))
}

summarize_window <- function(stage, c1, window, source_field) {
  if (anyNA(stage) || !all(stage %in% 0:3)) stop("Invalid stage values for ", window)
  dt <- data.table(stage = as.integer(stage), C1 = as.integer(c1), cluster = ifelse(c1 == 1L, "C1", "C2"))
  counts <- dt[, .(n = .N), by = .(stage, cluster)]
  counts <- merge(CJ(stage = 0:3, cluster = c("C1", "C2"), sorted = TRUE), counts, by = c("stage", "cluster"), all.x = TRUE)
  counts[is.na(n), n := 0L]
  counts[, row_total := sum(n), by = stage]
  counts[, column_total := sum(n), by = cluster]
  counts[, `:=`(
    row_percent_P_cluster_given_stage = 100 * n / row_total,
    column_percent_P_stage_given_cluster = 100 * n / column_total,
    time_window = window,
    source_field = source_field
  )]
  setcolorder(counts, c(
    "time_window", "source_field", "stage", "cluster", "n", "row_total",
    "row_percent_P_cluster_given_stage", "column_total", "column_percent_P_stage_given_cluster"
  ))

  matrix <- table(dt$stage, dt$C1)
  chi <- suppressWarnings(chisq.test(matrix, correct = FALSE))
  cramer_v <- sqrt(as.numeric(chi$statistic) / (nrow(dt) * min(nrow(matrix) - 1L, ncol(matrix) - 1L)))
  trend <- prop.trend.test(
    x = vapply(0:3, function(s) sum(dt$C1 == 1L & dt$stage == s), numeric(1)),
    n = vapply(0:3, function(s) sum(dt$stage == s), numeric(1)),
    score = 0:3
  )
  ordinal_fit <- glm(C1 ~ stage, data = dt, family = binomial())
  ordinal_or <- or_from_glm(ordinal_fit, "stage")

  stage3_subset <- dt[stage %in% c(0L, 1L, 3L)]
  stage3_table <- table(
    stage3 = factor(as.integer(stage3_subset$stage == 3L), levels = 0:1),
    C1 = factor(stage3_subset$C1, levels = 0:1)
  )
  stage3_test <- fisher.test(stage3_table)

  association <- cbind(data.table(
    time_window = window,
    source_field = source_field,
    n = nrow(dt),
    chi_square = as.numeric(chi$statistic),
    chi_square_df = unname(chi$parameter),
    chi_square_p = chi$p.value,
    cramer_v = cramer_v,
    trend_chi_square = as.numeric(trend$statistic),
    trend_df = unname(trend$parameter),
    trend_p = trend$p.value,
    ordinal_stage_or_interpretation = "odds ratio for C1 membership per one-stage increase"
  ), ordinal_or[, .(
    ordinal_stage_or = or,
    ordinal_stage_or_ci_low = or_ci_low,
    ordinal_stage_or_ci_high = or_ci_high
  )], data.table(
    stage3_vs_stage01_or = unname(stage3_test$estimate),
    stage3_vs_stage01_or_ci_low = stage3_test$conf.int[1L],
    stage3_vs_stage01_or_ci_high = stage3_test$conf.int[2L],
    stage3_vs_stage01_p = stage3_test$p.value,
    stage2_excluded_from_stage3_or = TRUE,
    stage3_or_method = "Fisher exact conditional maximum-likelihood estimate"
  ))

  metrics <- rbindlist(list(
    binary_rule_metrics(dt$stage, dt$C1, function(x) x >= 2L, "KDIGO stage 2-3", window),
    binary_rule_metrics(dt$stage, dt$C1, function(x) x == 3L, "KDIGO stage 3", window)
  ))
  list(counts = counts, association = association, metrics = metrics)
}

export_plot <- function(plot, stem, width = 7.2, height = 4.8) {
  pdf_file <- file.path(DIRS$figures, paste0(stem, ".pdf"))
  png_file <- file.path(DIRS$figures, paste0(stem, ".png"))
  tiff_file <- file.path(DIRS$figures, paste0(stem, ".tiff"))
  cairo_pdf(pdf_file, width = width, height = height, family = "sans")
  print(plot)
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 300, type = "cairo")
  print(plot)
  dev.off()
  tiff(tiff_file, width = width, height = height, units = "in", res = 600, compression = "lzw")
  print(plot)
  dev.off()
}

labels <- as.data.table(readRDS(INPUT$labels))
baseline <- fread(INPUT$baseline, showProgress = FALSE)
final_full <- fread(INPUT$final_full, showProgress = FALSE)
if (nrow(labels) != 20049L || sum(labels$cluster_k2 == 1L) != 3992L) stop("Locked labels changed")
if (anyDuplicated(labels$stay_id) || anyDuplicated(baseline$stay_id) || anyDuplicated(final_full$stay_id)) stop("Duplicate stay_id")
idx24 <- match(labels$stay_id, baseline$stay_id)
idx7 <- match(labels$stay_id, final_full$stay_id)
if (anyNA(idx24) || anyNA(idx7)) stop("Unmatched stay_id")
if (!"aki_stage_0_24h" %in% names(baseline)) stop("Missing authoritative 0-24 h KDIGO stage")
if (!"aki_stage_kdigo_7d" %in% names(final_full)) stop("Missing 7-day combined KDIGO stage")

c1 <- as.integer(labels$cluster_k2 == 1L)
result24 <- summarize_window(
  baseline$aki_stage_0_24h[idx24], c1,
  "ICU admission 0-24 h", "baseline_covars.csv:aki_stage_0_24h"
)
result7 <- summarize_window(
  final_full$aki_stage_kdigo_7d[idx7], c1,
  "ICU admission through day 7", "final_full.csv:aki_stage_kdigo_7d"
)

criteria_status <- data.table(
  requested_analysis = c("creatinine-only AKI", "urine-output-only AKI", "both criteria", "neither criterion"),
  status = "MISSING_NOT_RECONSTRUCTED",
  reason = paste0(
    "The authoritative 0-24 h export retains only combined aki_stage_0_24h; ",
    "component criterion stages are not present in the locked CSV/RDS inputs."
  ),
  governing_rule = "Do not reconstruct criteria without an existing QC-approved field.",
  denominator_used = 0L
)

write_csv_atomic(criteria_status, file.path(DIRS$tables, "Table_D1_AKI_criterion_by_cluster.csv"))
write_csv_atomic(result24$counts, file.path(DIRS$tables, "Table_D1_KDIGO_24h_by_cluster.csv"))
write_csv_atomic(result7$counts, file.path(DIRS$tables, "Table_D1_KDIGO_7d_by_cluster.csv"))
write_csv_atomic(rbindlist(list(result24$metrics, result7$metrics), fill = TRUE), file.path(DIRS$tables, "Table_D1_KDIGO_classification_metrics.csv"))
write_csv_atomic(rbindlist(list(result24$association, result7$association), fill = TRUE), file.path(DIRS$tables, "29a_Table_D1_KDIGO_association_summary.csv"))

plot24 <- copy(result24$counts)
plot24[, row_fraction := row_percent_P_cluster_given_stage / 100]
plot24[, label := sprintf("%s\n%d (%.1f%%)", cluster, n, row_percent_P_cluster_given_stage)]
p <- ggplot(plot24, aes(x = cluster, y = factor(stage, levels = 3:0), fill = row_fraction)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.2, family = "sans") +
  scale_fill_gradient(
    low = "#F2F2F2", high = "#3A5A8C", limits = c(0, 1),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = NULL, y = "KDIGO stage (0-24 h)", fill = "Within-stage\nproportion",
    title = "KDIGO stage and locked K2 membership",
    subtitle = "Percentages are P(cluster | KDIGO stage); time window is ICU admission 0-24 h"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey35"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5)
  )
export_plot(p, "Figure_D1_KDIGO_cluster_heatmap")

write_lines_atomic(c(
  "NOT GENERATED",
  "The AKI criterion heatmap was not generated because authoritative exported 0-24 h creatinine-only and urine-output-only fields are missing.",
  "No criteria were reconstructed from raw data."
), file.path(DIRS$figures, "Figure_D1_AKI_criterion_cluster_heatmap_NOT_GENERATED.txt"))

log_lines <- c(
  "# 29a AKI definition alignment log",
  "",
  paste0("- Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Locked labels: 20,049; C1=3,992; C2=16,057.",
  "- Main KDIGO window: ICU admission 0-24 h, combined creatinine/urine-output/RRT stage from baseline_covars.csv.",
  "- Sensitivity KDIGO window: ICU admission through day 7 from final_full.csv.",
  "- The two time windows are reported separately.",
  "- AKI criterion decomposition is missing and was not reconstructed.",
  "",
  "## 0-24 h association",
  paste(capture.output(print(result24$association)), collapse = "\n"),
  "",
  "## 7-day association",
  paste(capture.output(print(result7$association)), collapse = "\n"),
  "",
  "## Classification metrics",
  paste(capture.output(print(rbindlist(list(result24$metrics, result7$metrics)))), collapse = "\n")
)
write_lines_atomic(log_lines, file.path(DIRS$logs, "29a_domain1_aki_alignment_log.md"))
capture.output(sessionInfo(), file = file.path(DIRS$logs, "29a_sessionInfo.txt"))
write_lines_atomic("RUN COMPLETED WITH PRESPECIFIED MISSING AKI-CRITERION DECOMPOSITION", file.path(OUT_ROOT, "29a_run_completed.ok"))

cat("WP2 supported analyses completed; criterion decomposition remains missing.\n")
print(result24$counts)
print(result24$metrics)
print(result24$association)
