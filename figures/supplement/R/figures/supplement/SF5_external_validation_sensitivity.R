# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 5 - external-validation sensitivity and calibration.
##
## Presentation-only script. It compares the current simulated Sepsis-3 eICU cohort
## (n=17,465) with the prespecified stricter SOFA+AKI subset (n=7,975). Model-
## specific complete-case counts remain those in the locked input files.
## The pending one-stay analysis is not included.
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
PACKAGE_ROOT <- normalizePath(file.path(dirname(script_file), "..", "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(PACKAGE_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")
OUTPUT_DIR <- file.path(PACKAGE_ROOT, "figures", "supplement")
SOURCE_DIR <- file.path(PACKAGE_ROOT, "data", "figure_source", "supplement")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

PERF_FILE <- file.path(SOURCE_DIR, "SF5_external_model_performance.csv")
PAIR_FILE <- file.path(SOURCE_DIR, "SF5_key_paired_dAUC.csv")
PREV_FILE <- file.path(SOURCE_DIR, "SF5_prevalence_sensitivity.csv")
AUDIT_FILE <- file.path(SOURCE_DIR, "SF5_subset_audit.csv")
UPSTREAM_QA_FILE <- file.path(SOURCE_DIR, "SF5_upstream_QA.csv")

required <- c(PERF_FILE, PAIR_FILE, PREV_FILE, AUDIT_FILE, UPSTREAM_QA_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF5 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 180
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
cat("Supplementary Figure 5 - presentation-only\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Cohorts: eICU n=17,465; stricter SOFA+AKI subset n=7,975.\n")
cat("Pending one-stay sensitivity is excluded.\n")
cat("============================================================\n")

perf <- fread(PERF_FILE, encoding = "UTF-8")
pair <- fread(PAIR_FILE, encoding = "UTF-8")
prev <- fread(PREV_FILE, encoding = "UTF-8")
audit <- fread(AUDIT_FILE, encoding = "UTF-8")
upstream_qa <- fread(UPSTREAM_QA_FILE, encoding = "UTF-8")

stopifnot(
  all(c("cohort", "outcome", "model", "n", "events", "auc",
        "calibration_slope", "calibration_intercept") %in% names(perf)),
  all(c("cohort", "outcome", "better", "worse", "dAUC", "dAUC_lo", "dAUC_hi") %in% names(pair)),
  all(c("cohort", "n", "C1_pct") %in% names(prev)),
  all(c("check", "pass") %in% names(upstream_qa)),
  all(upstream_qa$pass),
  prev[cohort == "eicu_full", n] == 17465L,
  prev[cohort == "eicu_sepsis3_like", n] == 7975L
)

cohort_levels <- c("eicu_full", "eicu_sepsis3_like")
cohort_labels <- c(
  "eicu_full" = "Full eICU",
  "eicu_sepsis3_like" = "Strict SOFA+AKI"
)
cohort_shapes <- c(
  "Full eICU" = 16,
  "Strict SOFA+AKI" = 1
)
cohort_lines <- c(
  "Full eICU" = "solid",
  "Strict SOFA+AKI" = "22"
)

model_labels <- c(
  "base" = "Base",
  "base_phen" = "Base + K2",
  "feat_pen" = "Raw-EN",
  "feat_rf" = "Raw-RF"
)
model_levels <- c("base", "base_phen", "feat_pen", "feat_rf")
outcome_labels <- c(
  "hospital_mortality" = "In-hospital mortality",
  "make_hosp" = "eICU composite"
)

perf[, cohort_display := factor(cohort, levels = cohort_levels, labels = cohort_labels[cohort_levels])]
perf[, model_display := factor(model, levels = model_levels, labels = model_labels[model_levels])]
perf[, outcome_display := factor(outcome, levels = names(outcome_labels), labels = outcome_labels)]

theme_pub <- theme_supplement() +
  theme(
    plot.title = element_text(size = FONT_PT["panel_title"], face = "bold", hjust = 0),
    plot.subtitle = element_text(size = FONT_PT["subtitle"], colour = COLORS["neutral_main"], hjust = 0),
    axis.title = element_text(size = FONT_PT["axis_title"]),
    axis.text = element_text(size = FONT_PT["axis_text"], colour = "black"),
    strip.text = element_text(size = FONT_PT["annotation"], face = "bold"),
    strip.background = element_rect(fill = COLORS["pale_grey"], colour = NA),
    legend.title = element_text(size = FONT_PT["legend"], face = "bold"),
    legend.text = element_text(size = FONT_PT["legend"]),
    legend.position = "bottom",
    panel.spacing = unit(3, "mm"),
    plot.margin = margin(5, 7, 5, 7)
  )

prev_plot <- prev[cohort %chin% c("MIMIC_training", cohort_levels)]
prev_plot[, cohort_display := factor(
  cohort,
  levels = c("MIMIC_training", "eicu_full", "eicu_sepsis3_like"),
  labels = c("MIMIC-IV", "Full eICU", "Strict SOFA+AKI")
)]
prev_plot[, count_label := sprintf("%.1f%%", C1_pct)]

pA <- ggplot(prev_plot, aes(cohort_display, C1_pct)) +
  geom_col(width = 0.58, fill = COLORS["C1"]) +
  geom_text(aes(label = count_label), vjust = -0.25, size = 3.1,
            family = BASE_FAMILY, colour = COLORS["neutral_focus"], lineheight = 0.95) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 40, by = 10),
                     labels = label_percent(scale = 1)) +
  scale_x_discrete(labels = c("MIMIC-IV"="MIMIC-IV", "Full eICU"="Full eICU",
                              "Strict SOFA+AKI"="Strict\nSOFA+AKI")) +
  labs(
    title = "C1 prevalence",
    x = NULL, y = "C1 prevalence"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(lineheight = 0.95), legend.position = "none")

pB <- ggplot(perf, aes(auc, model_display, shape = cohort_display)) +
  geom_point(colour = COLORS["neutral_focus"], size = 2.0,
             position = position_dodge(width = 0.48)) +
  facet_wrap(~outcome_display, nrow = 1) +
  scale_shape_manual(values = cohort_shapes) +
  scale_x_continuous(
    limits = c(0.69, 0.82), breaks = c(0.70, 0.75, 0.80),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "External AUC",
    x = "External AUC", y = NULL, shape = "External cohort definition"
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

pair_key <- pair[better == "feat_pen" & worse == "base_phen"]
stopifnot(nrow(pair_key) == 4L)
pair_key[, cohort_display := factor(cohort, levels = cohort_levels, labels = cohort_labels[cohort_levels])]
pair_key[, outcome_display := factor(outcome, levels = names(outcome_labels), labels = outcome_labels)]

pC <- ggplot(pair_key, aes(dAUC, outcome_display, shape = cohort_display)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.5) +
  geom_errorbar(aes(xmin = dAUC_lo, xmax = dAUC_hi), orientation = "y",
                width = 0.16, linewidth = 0.25, colour = COLORS["neutral_main"],
                position = position_dodge(width = 0.45)) +
  geom_point(size = 2.1, colour = COLORS["neutral_focus"],
             position = position_dodge(width = 0.45)) +
  scale_shape_manual(values = cohort_shapes) +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title = "Raw-EN versus Base + K2",
    x = "Delta AUC", y = NULL, shape = "External cohort definition"
  ) +
  theme_pub +
  theme(legend.position = "none")

cal <- melt(
  perf,
  id.vars = c("cohort_display", "outcome_display", "model_display"),
  measure.vars = c("calibration_slope", "calibration_intercept"),
  variable.name = "metric", value.name = "value"
)
cal[, metric_display := factor(
  metric,
  levels = c("calibration_slope", "calibration_intercept"),
  labels = c("Calibration slope", "Calibration intercept")
)]
cal[, alert_flag := fifelse(
  metric == "calibration_slope",
  value < 0.80 | value > 1.20,
  abs(value) > 0.20
)]
cal[, status := factor(
  fifelse(alert_flag, "Alert", "Within range"),
  levels = c("Within range", "Alert")
)]
cal[, facet_label := factor(
  paste(
    fifelse(metric == "calibration_slope", "Slope", "Intercept"),
    fifelse(as.character(outcome_display) == "In-hospital mortality", "mortality", "eICU composite"),
    sep = ":\n"
  ),
  levels = c(
    "Slope:\nmortality", "Slope:\neICU composite",
    "Intercept:\nmortality", "Intercept:\neICU composite"
  )
)]

ref_lines <- unique(cal[, .(
  facet_label,
  reference = fifelse(metric == "calibration_slope", 1, 0)
)])

pD <- ggplot(cal, aes(value, model_display, colour = status, shape = cohort_display)) +
  geom_vline(data = ref_lines, aes(xintercept = reference),
             linetype = "dashed", colour = COLORS["reference"], linewidth = 0.5,
             inherit.aes = FALSE) +
  # Preserve the original axis ranges without drawing the inapplicable line.
  geom_blank(data = data.table(value = c(0, 1)), aes(x = value), inherit.aes = FALSE) +
  geom_point(size = 1.9, position = position_dodge(width = 0.48)) +
  facet_wrap(~facet_label, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c(
    "Within range" = unname(COLORS["neutral_focus"]),
    "Alert" = unname(COLORS["alert"])
  )) +
  scale_shape_manual(values = cohort_shapes) +
  labs(
    title = "Calibration sensitivity",
    x = "Calibration estimate", y = NULL,
    colour = "Calibration", shape = "External cohort definition"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    legend.box = "vertical",
    legend.box.just = "left"
  ) +
  guides(
    shape = guide_legend(order = 1, nrow = 1, override.aes = list(colour = COLORS["neutral_focus"])),
    colour = guide_legend(order = 2, nrow = 1, override.aes = list(shape = 16))
  )

pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C")
pD <- tag_panel(pD, "D")

legend_data <- data.table(x = c(0.02,0.31,0.62,0.81), y = .5,
  label = c("Full eICU", "Strict SOFA+AKI", "Within range", "Calibration alert"),
  shape = c(16,1,16,16), colour = c(rep(COLORS[["neutral_focus"]],3),COLORS[["alert"]]))
p_legend <- ggplot(legend_data, aes(x,y)) +
  geom_point(aes(shape=I(shape),colour=I(colour)),size=2) +
  geom_text(aes(label=label),nudge_x=.02,hjust=0,
            size=8/ggplot2::.pt,family=FONT_FAMILY,colour=PAL[["ink"]]) +
  coord_cartesian(xlim=c(0,1.02),ylim=c(0,1),clip="off") +
  theme_void() + theme(plot.margin=margin(0,8,0,8))
top <- (pA | pB) + plot_layout(widths=c(.85,1.5))
fig <- free(top) / pC / pD / p_legend +
  plot_layout(heights = c(1.05, .67, 1.68, .15)) +
  plot_annotation(
    caption = "B-D: point shape denotes the external cohort definition; orange denotes a calibration alert."
  ) &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "D3  Reconstructed-surrogate application sensitivity")

RENDER_AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(RENDER_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(perf, file.path(RENDER_AUDIT_DIR, "SF5_external_model_performance.csv"))
fwrite(pair_key, file.path(RENDER_AUDIT_DIR, "SF5_key_paired_dAUC.csv"))
fwrite(prev_plot, file.path(RENDER_AUDIT_DIR, "SF5_prevalence_sensitivity.csv"))
fwrite(audit, file.path(RENDER_AUDIT_DIR, "SF5_subset_audit.csv"))
fwrite(upstream_qa, file.path(RENDER_AUDIT_DIR, "SF5_upstream_QA.csv"))
fwrite(cal, file.path(RENDER_AUDIT_DIR, "SF5_calibration_long.csv"))
fwrite(data.table(
  input = file.path("figures", "supplement_source_data", basename(required)),
  md5 = unname(tools::md5sum(required))
), file.path(RENDER_AUDIT_DIR, "SF5_input_manifest.csv"))

stem_name <- "Supplementary_Figure_S5_surrogate_application_sensitivity_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
qa <- data.table(
  check = c(
    "all_required_inputs_exist", "full_eICU_n_17465", "strict_subset_n_7975",
    "upstream_current_version_QA_passed",
    "four_models_two_outcomes_two_external_definitions", "four_key_paired_contrasts",
    "one_stay_not_included", "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)), prev[cohort == "eicu_full", n] == 17465L,
    prev[cohort == "eicu_sepsis3_like", n] == 7975L,
    all(upstream_qa$pass),
    nrow(perf) == 16L, nrow(pair_key) == 4L,
    !any(grepl("one.stay", c(PERF_FILE, PAIR_FILE, PREV_FILE, AUDIT_FILE), ignore.case = TRUE)),
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(OUTPUT_DIR, "Supplementary_Figure_S5_render_QA.csv"))
if (!all(qa$pass)) stop("SF5 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("SF5 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")

legend_text <- paste0(
  "Supplementary Figure S5. Reconstructed-surrogate application sensitivity. ",
  "Analyses compare the full eICU cohort with the stricter SOFA-plus-AKI subset. ",
  "Cohort definition is encoded by point shape. Panels B and C use neutral marks because they compare model performance; ",
  "in panel D, orange identifies estimates outside the prespecified calibration ranges. ",
  "The four visible combinations arise from two independent encodings rather than a four-level categorical state."
)
writeLines(legend_text, file.path(OUTPUT_DIR, paste0(stem_name, "_legend.txt")), useBytes = TRUE)
