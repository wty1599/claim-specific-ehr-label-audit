# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 8 - SA-02 oracle-label accuracy dose response.
##
## Presentation-only script. It reads the corrected canonical 12-row summary.
## The label is a fixed noisy oracle hidden-state label, not the empirical K=2
## phenotype. No source result or post-processing classification is changed.
## ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE),
                        error = function(e) NA_character_)
if (is.na(script_file)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) script_file <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
}
if (is.na(script_file)) stop("Run this script with Rscript or source() so its package-relative path can be resolved.")
WORK_ROOT <- normalizePath(file.path(dirname(script_file), "..", "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(WORK_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")
OUTPUT_DIR <- file.path(WORK_ROOT, "figures", "supplement")
SOURCE_DIR <- file.path(WORK_ROOT, "data", "figure_source", "supplement")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

SUMMARY_FILE <- file.path(SOURCE_DIR, "SF8_SA02_corrected_12row_summary.csv")
FIDELITY_FILE <- file.path(SOURCE_DIR, "SF8_SA02_oracle_fidelity_summary.csv")

required <- c(SUMMARY_FILE, FIDELITY_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF8 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 150
DPI <- 600
BASE_FAMILY <- "sans"

COLORS <- c(
  structure = "#2E4B6E",
  alert = "#C77A30",
  neutral_focus = "#333333",
  neutral_main = "#4A4A4A",
  neutral_line = "grey55",
  reference = "grey70",
  C1 = "#B4553F",
  C2 = "#3E6E93",
  pale_grey = "#F2F2F2"
)

FONT_PT <- c(
  panel_title = 11,
  panel_tag = 11,
  subtitle = 8.5,
  axis_title = 9,
  axis_text = 8,
  annotation = 8,
  legend = 8
)
activate_supplement_style()

cat("============================================================\n")
cat("Supplementary Figure 8 - presentation-only\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Label source: fixed noisy oracle hidden state (not empirical K=2).\n")
cat("============================================================\n")

d <- fread(SUMMARY_FILE, encoding = "UTF-8")
fidelity <- fread(FIDELITY_FILE, encoding = "UTF-8")

stopifnot(
  all(c("oracle_accuracy", "outcome", "outcome_label", "n_repeats", "m_used_min",
        "m_used_max", "mean_delta_auc", "mean_ci_low", "mean_ci_high",
        "repeat_rubin_ci_positive_proportion", "repeat_rubin_ci_positive_ci_low",
        "repeat_rubin_ci_positive_ci_high", "practical_null_threshold",
        "classification", "domain2_label_source") %in% names(d)),
  nrow(d) == 12L,
  nrow(fidelity) == 12L,
  all(d$n_repeats == 100L),
  all(d$m_used_min == 5L & d$m_used_max == 5L)
)

accuracy_levels <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
stopifnot(identical(sort(unique(d$oracle_accuracy)), accuracy_levels))
d[, signal_status := factor(
  fifelse(practical_null, "Chance-level practical null", "Incremental signal detected"),
  levels = c("Chance-level practical null", "Incremental signal detected")
)]
d[, detection_status := factor(
  fifelse(repeat_rubin_ci_positive_proportion == 0, "No positive detections", "Positive detections observed"),
  levels = c("No positive detections", "Positive detections observed")
)]

theme_pub <- theme_supplement() +
  theme(
    plot.title = element_text(size = FONT_PT["panel_title"], face = "bold", hjust = 0),
    plot.subtitle = element_text(size = FONT_PT["subtitle"], colour = COLORS["neutral_main"], hjust = 0),
    axis.title = element_text(size = FONT_PT["axis_title"]),
    axis.text = element_text(size = FONT_PT["axis_text"], colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = FONT_PT["legend"]),
    legend.position = "bottom",
    plot.margin = margin(5, 7, 5, 7)
  )

signal_cols <- c(
  "Chance-level practical null" = unname(COLORS["neutral_focus"]),
  "Incremental signal detected" = unname(COLORS["alert"])
)
detection_cols <- c(
  "No positive detections" = unname(COLORS["neutral_focus"]),
  "Positive detections observed" = unname(COLORS["alert"])
)

make_delta_panel <- function(outcome_code, title, show_y = TRUE) {
  z <- d[outcome == outcome_code][order(oracle_accuracy)]
  null_gate <- unique(z$practical_null_threshold)
  stopifnot(length(null_gate) == 1L)
  ggplot(z, aes(oracle_accuracy, mean_delta_auc, group = 1)) +
    annotate(
      "rect", xmin = -Inf, xmax = Inf, ymin = -null_gate, ymax = null_gate,
      fill = "grey85", alpha = 0.55
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
    geom_line(colour = COLORS["neutral_line"], linewidth = 0.55) +
    geom_errorbar(
      aes(ymin = mean_ci_low, ymax = mean_ci_high, colour = signal_status),
      width = 0.012, linewidth = 0.5
    ) +
    geom_point(aes(colour = signal_status), size = 2.3) +
    scale_colour_manual(values = signal_cols, guide = "none") +
    scale_x_continuous(breaks = accuracy_levels, limits = c(0.48, 1.02), labels = number_format(accuracy = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    labs(
      title = title,
      x = "Fixed noisy-oracle accuracy",
      y = if (show_y) "Delta AUC" else NULL
    ) +
    theme_pub + theme(legend.position = "none")
}

make_detection_panel <- function(outcome_code, title, show_y = TRUE) {
  z <- d[outcome == outcome_code][order(oracle_accuracy)]
  ggplot(z, aes(oracle_accuracy, repeat_rubin_ci_positive_proportion, group = 1)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
    geom_line(colour = COLORS["neutral_line"], linewidth = 0.55) +
    geom_errorbar(
      aes(
        ymin = repeat_rubin_ci_positive_ci_low,
        ymax = repeat_rubin_ci_positive_ci_high,
        colour = detection_status
      ),
      width = 0.012, linewidth = 0.5
    ) +
    geom_point(aes(colour = detection_status), size = 2.3) +
    scale_colour_manual(values = detection_cols, guide = "none") +
    scale_x_continuous(breaks = accuracy_levels, limits = c(0.48, 1.02), labels = number_format(accuracy = 0.1)) +
    scale_y_continuous(
      limits = c(0, 1.04), breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title = title,
      x = "Fixed noisy-oracle accuracy",
      y = if (show_y) "Positive detection rate" else NULL
    ) +
    theme_pub + theme(legend.position = "none")
}

pA <- make_delta_panel("make30", "Delta AUC: renal composite", TRUE)
pB <- make_delta_panel("mortality_30d", "Delta AUC: 30-day mortality", FALSE)
pC <- make_detection_panel("make30", "Detection rate: renal composite", TRUE)
pD <- make_detection_panel("mortality_30d", "Detection rate: 30-day mortality", FALSE)
pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C")
pD <- tag_panel(pD, "D")

fig <- (pA | pB) / (pC | pD) +
  plot_layout() +
  plot_annotation() &
  theme(
    plot.caption = element_text(size = FONT_PT["annotation"], colour = COLORS["neutral_main"], hjust = 0),
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "D2  Noisy-oracle label-accuracy sensitivity")

RENDER_AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(RENDER_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(d, file.path(RENDER_AUDIT_DIR, "SF8_label_accuracy_12row_summary.csv"))
fwrite(fidelity, file.path(RENDER_AUDIT_DIR, "SF8_SA02_oracle_fidelity_summary.csv"))
fwrite(data.table(
  input = file.path("figures", "supplement_source_data", basename(required)),
  md5 = unname(tools::md5sum(required))
), file.path(RENDER_AUDIT_DIR, "SF8_input_manifest.csv"))

stem_name <- "Supplementary_Figure_S8_oracle_accuracy_dose_response_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

monotonic <- d[order(oracle_accuracy), .(
  strictly_increasing_delta_auc = all(diff(mean_delta_auc) > 0),
  spearman_rho = cor(oracle_accuracy, mean_delta_auc, method = "spearman")
), by = outcome]
outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
legend_text <- paste(
  "Supplementary Figure S8. Noisy-oracle label-accuracy sensitivity.",
  "Panels A and B show mean delta AUC and 95% intervals across 100 repeats for the renal composite and 30-day mortality.",
  "Panels C and D show the corresponding positive-detection rates with Wilson intervals.",
  "The gray band is the prespecified practical-null range for delta AUC; the dashed 80% line is a detection-performance reference, not a single-repeat decision threshold.",
  "Dark points at accuracy 0.5 denote chance-level practical-null results in A-B and no positive detections in C-D; orange denotes a non-null or detected signal."
)
writeLines(legend_text, paste0(stem, "_legend.txt"), useBytes = TRUE)
qa <- data.table(
  check = c(
    "all_required_inputs_exist", "canonical_rows_12", "six_accuracy_levels",
    "one_hundred_repeats_each", "five_imputations_each",
    "chance_level_practical_null_both_outcomes", "delta_auc_strictly_increasing_both_outcomes",
    "spearman_rho_one_both_outcomes", "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)), nrow(d) == 12L,
    identical(sort(unique(d$oracle_accuracy)), accuracy_levels), all(d$n_repeats == 100L),
    all(d$m_used_min == 5L & d$m_used_max == 5L),
    all(d[oracle_accuracy == 0.5, practical_null]),
    all(monotonic$strictly_increasing_delta_auc), all(monotonic$spearman_rho == 1),
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(OUTPUT_DIR, "Supplementary_Figure_S8_render_QA.csv"))
if (!all(qa$pass)) stop("SF8 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("SF8 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")
