# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 10 - SA-05 missingness mechanism x intensity.
##
## Presentation-only script. Domain 1 and Domain 2 verdicts come from the
## locked framework-aligned postprocessing tables. No threshold is recomputed.
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

REALIZED_FILE <- file.path(SOURCE_DIR, "SF10_SA05_realized_missingness_summary.csv")
D1_FILE <- file.path(SOURCE_DIR, "SF10_SA05_D1_framework_rates.csv")
D2_NULL_FILE <- file.path(SOURCE_DIR, "SF10_SA05_D2_null_false_positive_rates.csv")
D2_POS_FILE <- file.path(SOURCE_DIR, "SF10_SA05_D2_positive_control_detection.csv")
QC_FILE <- file.path(SOURCE_DIR, "SF10_SA05_critical_QC.csv")
FALLBACK_FILE <- file.path(SOURCE_DIR, "SF10_SA05_kmeans_fallback_audit.csv")

required <- c(
  REALIZED_FILE, D1_FILE, D2_NULL_FILE, D2_POS_FILE, QC_FILE, FALLBACK_FILE
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF10 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 210
DPI <- 600
BASE_FAMILY <- "sans"

COLORS <- c(
  structure = "#2E4B6E",
  alert = "#C77A30",
  neutral_focus = "#333333",
  neutral_main = "#4A4A4A",
  neutral_mid = "#77838E",
  neutral_light = "#A7A7A7",
  reference = "grey70",
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
cat("Supplementary Figure 10 - presentation-only\n")
cat("Core conclusion: the four-domain verdicts remain stable across\n")
cat("MCAR/MAR/MNAR missingness at 0%, 15%, and 30%.\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Amber marks a setting with at least one classification error.\n")
cat("============================================================\n")

realized <- fread(REALIZED_FILE, encoding = "UTF-8")
d1 <- fread(D1_FILE, encoding = "UTF-8")
d2 <- rbindlist(list(
  fread(D2_NULL_FILE, encoding = "UTF-8"),
  fread(D2_POS_FILE, encoding = "UTF-8")
), use.names = TRUE, fill = TRUE)
qc <- fread(QC_FILE, encoding = "UTF-8")
fallback <- fread(FALLBACK_FILE, encoding = "UTF-8")

stopifnot(
  all(c("scenario", "missingness_mechanism", "target_missingness", "repeats",
        "realized_mean", "realized_q025", "realized_q975") %in% names(realized)),
  all(c("scenario", "missingness_mechanism", "target_missingness", "metric",
        "numerator", "denominator", "proportion", "wilson_low", "wilson_high") %in% names(d1)),
  all(c("scenario", "missingness_mechanism", "target_missingness", "outcome", "metric",
        "numerator", "denominator", "proportion", "wilson_low", "wilson_high") %in% names(d2)),
  nrow(realized) == 27L,
  nrow(d1) == 18L,
  nrow(d2) == 54L,
  all(qc[critical == TRUE, pass] == TRUE),
  nrow(fallback) == 54L
)

fallback[, fallback_fits := mean_imputation_fallback_rate * repeats * 5]
fallback_fits_n <- sum(fallback$fallback_fits)
fallback_fits_total <- sum(fallback$repeats) * 5

scenario_labels <- c(
  "S1_pure_severity_continuum" = "Continuous\ngenerator",
  "S2_true_discrete_subtypes" = "Three-component\ndiscrete control",
  "S5_domain2_positive_control" = "Noisy-oracle\ngenerator"
)
mechanism_levels <- c("MCAR", "MAR", "MNAR")
mechanism_cols <- c(
  "MCAR" = unname(COLORS["structure"]),
  "MAR" = unname(COLORS["neutral_mid"]),
  "MNAR" = unname(COLORS["neutral_main"])
)
mechanism_shapes <- c("MCAR" = 16, "MAR" = 17, "MNAR" = 15)
mechanism_lines <- c("MCAR" = "solid", "MAR" = "dashed", "MNAR" = "dotted")

for (z in list(realized, d1, d2)) {
  stopifnot(all(z$missingness_mechanism %in% mechanism_levels))
}
realized[, mechanism := factor(missingness_mechanism, levels = mechanism_levels)]
d1[, mechanism := factor(missingness_mechanism, levels = mechanism_levels)]
d2[, mechanism := factor(missingness_mechanism, levels = mechanism_levels)]
realized[, scenario_display := factor(scenario_labels[scenario], levels = unname(scenario_labels))]
d1[, scenario_display := factor(scenario_labels[scenario], levels = unname(scenario_labels))]
d2[, scenario_display := factor(scenario_labels[scenario], levels = unname(scenario_labels))]
realized[, scenario_short := factor(
  scenario,
  levels = c("S1_pure_severity_continuum", "S2_true_discrete_subtypes", "S5_domain2_positive_control"),
  labels = c("Continuous\ncontrol", "Three-component\ndiscrete", "Noisy-oracle\ncontrol")
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
    strip.text = element_text(size = FONT_PT["subtitle"], face = "bold"),
    strip.background = element_rect(fill = COLORS["pale_grey"], colour = NA),
    panel.spacing = grid::unit(4, "mm"),
    plot.margin = margin(5, 7, 5, 7)
  )

common_mechanism_scales <- list(
  scale_colour_manual(values = mechanism_cols, limits = mechanism_levels, drop = FALSE),
  scale_shape_manual(values = mechanism_shapes, limits = mechanism_levels, drop = FALSE),
  scale_linetype_manual(values = mechanism_lines, limits = mechanism_levels, drop = FALSE)
)

pA <- ggplot(
  realized,
  aes(target_missingness, realized_mean, colour = mechanism, shape = mechanism, linetype = mechanism)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(aes(ymin = realized_q025, ymax = realized_q975), width = 0.008, linewidth = 0.42) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 2.0) +
  facet_wrap(~ scenario_short, nrow = 1, labeller = labeller(
    scenario_short = function(x) ifelse(x == "Three-component\ndiscrete", "Discrete\ncontrol", x)
  )) +
  common_mechanism_scales +
  scale_x_continuous(breaks = c(0, 0.15, 0.30), labels = percent_format(accuracy = 1), limits = c(-0.015, 0.315)) +
  scale_y_continuous(breaks = c(0, 0.15, 0.30), labels = percent_format(accuracy = 1), limits = c(-0.015, 0.315)) +
  labs(
    title = "Missingness fidelity",
    x = "Target missingness", y = "Realized missingness"
  ) +
  theme_pub +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

d1_plot <- d1[metric %in% c("D1_joint_false_positive_rate", "D1_joint_detection_rate")]
d1_plot[, endpoint := factor(
  metric,
  levels = c("D1_joint_false_positive_rate", "D1_joint_detection_rate"),
  labels = c("False-positive rate\n(continuous control)", "Detection rate\n(discrete control)")
)]
d1_plot[, any_error := fifelse(metric == "D1_joint_false_positive_rate", numerator > 0, numerator < denominator)]
ring_label <- "Any unexpected decision"
zero_note <- data.table(
  endpoint = factor("False-positive rate\n(continuous control)", levels = levels(d1_plot$endpoint)),
  x = 0.15, y = 0.12, label = "0/50 in every setting"
)

pB <- ggplot(
  d1_plot,
  aes(target_missingness, proportion, colour = mechanism, shape = mechanism, linetype = mechanism)
) +
  geom_errorbar(aes(ymin = wilson_low, ymax = wilson_high), width = 0.008, linewidth = 0.42) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 2.0) +
  geom_point(
    data = d1_plot[any_error == TRUE], aes(fill = ring_label),
    colour = COLORS["alert"], shape = 21, stroke = 0.8, size = 2.6
  ) +
  geom_text(
    data = zero_note, aes(x = x, y = y, label = label), inherit.aes = FALSE,
    size = FONT_PT["annotation"] / ggplot2::.pt, colour = COLORS["neutral_main"], family = BASE_FAMILY
  ) +
  facet_wrap(~ endpoint, nrow = 1) +
  common_mechanism_scales +
  scale_x_continuous(breaks = c(0, 0.15, 0.30), labels = percent_format(accuracy = 1), limits = c(-0.015, 0.315)) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title = "D1 control decisions",
    x = "Target missingness", y = "Decision rate"
  ) +
  scale_fill_manual(values = setNames("white", ring_label), limits = ring_label, name = NULL) +
  theme_pub + theme(legend.position = "none")

d2_null <- d2[
  metric == "D2_joint_false_positive_rate" &
    scenario %in% c("S1_pure_severity_continuum", "S2_true_discrete_subtypes")
]
d2_null[, outcome_display := factor(outcome, levels = c("mortality_30d", "make30"), labels = c("Mortality", "Renal composite"))]
d2_null[, any_error := numerator > 0]

pC <- ggplot(
  d2_null,
  aes(target_missingness, proportion, colour = mechanism, shape = mechanism, linetype = mechanism)
) +
  geom_errorbar(aes(ymin = wilson_low, ymax = wilson_high), width = 0.008, linewidth = 0.42) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 1.9) +
  geom_point(
    data = d2_null[any_error == TRUE], aes(fill = ring_label),
    colour = COLORS["alert"], shape = 21, stroke = 0.8, size = 2.5
  ) +
  facet_grid(outcome_display ~ scenario_display) +
  common_mechanism_scales +
  scale_x_continuous(breaks = c(0, 0.15, 0.30), labels = percent_format(accuracy = 1), limits = c(-0.015, 0.315)) +
  scale_y_continuous(limits = c(-0.005, 0.145), breaks = c(0, 0.05, 0.10), labels = percent_format(accuracy = 1)) +
  labs(
    title = "D2 null false positives",
    x = "Target missingness", y = "False-positive rate"
  ) +
  scale_fill_manual(values = setNames("white", ring_label), name = NULL) +
  theme_pub + theme(legend.position = "none")

d2_pos <- d2[
  metric == "D2_joint_detection_rate" & scenario == "S5_domain2_positive_control"
]
d2_pos[, outcome_display := factor(outcome, levels = c("mortality_30d", "make30"), labels = c("Mortality", "Renal composite"))]
d2_pos[, any_error := numerator < denominator]

pD <- ggplot(
  d2_pos,
  aes(target_missingness, proportion, colour = mechanism, shape = mechanism, linetype = mechanism)
) +
  geom_errorbar(aes(ymin = wilson_low, ymax = wilson_high), width = 0.008, linewidth = 0.42) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 2.0) +
  geom_point(
    data = d2_pos[any_error == TRUE], aes(fill = ring_label),
    colour = COLORS["alert"], shape = 21, stroke = 0.8, size = 2.7
  ) +
  facet_wrap(~ outcome_display, nrow = 1) +
  common_mechanism_scales +
  scale_x_continuous(breaks = c(0, 0.15, 0.30), labels = percent_format(accuracy = 1), limits = c(-0.015, 0.315)) +
  scale_y_continuous(limits = c(0.86, 1.01), breaks = c(0.90, 0.95, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title = "D2 positive-control detection",
    x = "Target missingness", y = "Detection rate"
  ) +
  scale_fill_manual(values = setNames("white", ring_label), name = NULL) +
  theme_pub +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = COLORS["alert"], size = 2.7)))

layout_design <- "
AB
CD
EE
"

pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C")
pD <- tag_panel(pD, "D")

fig <- pA + pB + pC + pD + guide_area() +
  plot_layout(
    design = layout_design,
    guides = "collect",
    heights = c(0.70, 1.36, 0.13)
  ) +
  plot_annotation() &
  theme(
    legend.box = "horizontal",
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "Sensitivity  Missingness mechanism and intensity")

RENDER_AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(RENDER_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(realized, file.path(RENDER_AUDIT_DIR, "SF10_SA05_realized_missingness_summary.csv"))
fwrite(d1_plot, file.path(RENDER_AUDIT_DIR, "SF10_SA05_D1_framework_rates.csv"))
fwrite(d2_null, file.path(RENDER_AUDIT_DIR, "SF10_SA05_D2_null_false_positive_rates.csv"))
fwrite(d2_pos, file.path(RENDER_AUDIT_DIR, "SF10_SA05_D2_positive_control_detection.csv"))
fwrite(qc, file.path(RENDER_AUDIT_DIR, "SF10_SA05_critical_QC.csv"))
fwrite(fallback, file.path(RENDER_AUDIT_DIR, "SF10_SA05_kmeans_fallback_audit.csv"))
fwrite(data.table(
  input = file.path("figures", "supplement_source_data", basename(required)),
  md5 = unname(tools::md5sum(required))
),
       file.path(RENDER_AUDIT_DIR, "SF10_input_manifest.csv"))

stem_name <- "Supplementary_Figure_S10_missingness_robustness_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
legend_text <- paste(
  "Supplementary Figure S10. Missingness-mechanism sensitivity.",
  "Panel A compares realized and target missingness.",
  "Panels B-D show D1 and D2 control behavior across MCAR, MAR, and MNAR settings.",
  "Orange rings identify settings in which at least one repeat differed from expected control behavior; they do not denote crossing a 5% threshold.",
  "The continuous-control panel in B represents 0/50 false-positive decisions in every setting."
)
writeLines(legend_text, paste0(stem, "_legend.txt"), useBytes = TRUE)
qa_out <- data.table(
  check = c(
    "all_required_inputs_exist", "fifty_repeats_per_setting",
    "all_critical_QC_passed",
    "realized_missingness_fidelity", "D1_false_positive_zero_all_settings",
    "D1_detection_one_all_settings", "D2_null_false_positive_at_most_four_percent",
    "S5_detection_at_least_98_percent",
    "fallback_250_of_13500", "fallback_excludes_K2_and_S2_K3",
    "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)),
    all(d1$denominator == 50L) && all(d2$denominator == 50L),
    all(qc[critical == TRUE, pass] == TRUE),
    max(abs(realized$realized_mean - realized$target_missingness)) < 0.015,
    all(d1[metric == "D1_joint_false_positive_rate", proportion] == 0),
    all(d1[metric == "D1_joint_detection_rate", proportion] == 1),
    max(d2_null$proportion) <= 0.04,
    min(d2_pos$proportion) >= 0.98,
    fallback_fits_n == 250 && fallback_fits_total == 13500,
    all(fallback[k == 2, mean_imputation_fallback_rate] == 0) &&
      all(fallback[scenario == "S2_true_discrete_subtypes", mean_imputation_fallback_rate] == 0),
    all(file.exists(outputs))
  )
)
fwrite(qa_out, file.path(OUTPUT_DIR, "Supplementary_Figure_S10_render_QA.csv"))
if (!all(qa_out$pass)) stop("SF10 QA failed: ", paste(qa_out[pass == FALSE, check], collapse = "; "))

cat("SF10 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")
