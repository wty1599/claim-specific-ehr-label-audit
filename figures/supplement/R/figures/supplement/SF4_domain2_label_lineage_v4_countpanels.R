# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 4 - Domain 2 label-lineage analyses.
##
## Presentation-only script. It reads locked formal outputs and does not rerun
## imputation, clustering, prediction, bootstrap, or Rubin pooling.
## ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", ".."), winslash = "/")
source(file.path(PROJECT_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")
source(file.path(SCRIPT_DIR, "00_count_panel_revision_helpers_v4.R"), encoding = "UTF-8")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "figures", "supplement")
SOURCE_DIR <- file.path(PROJECT_ROOT, "data", "figure_source", "supplement")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
QA_DIR <- file.path(PROJECT_ROOT, "qa")
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

AUC_FILE <- file.path(SOURCE_DIR, "SF4_crossfitted_auc_ladder.csv")
SUMMARY_FILE <- file.path(SOURCE_DIR, "SF4_crossfit_formal_summary.csv")
DAUC_FILE <- file.path(SOURCE_DIR, "SF4_locked_imputation_foldwise_delta_auc.csv")
FOLD_LOG_FILE <- file.path(SOURCE_DIR, "SF4_crossfit_fold_log.csv")
LOCKED_FILE <- file.path(SOURCE_DIR, "SF4_locked_variants.csv")

required <- c(AUC_FILE, SUMMARY_FILE, DAUC_FILE, FOLD_LOG_FILE, LOCKED_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF4 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 158
DPI <- 600
BASE_FAMILY <- "sans"

COLORS <- c(
  structure = "#2E4B6E",
  indigo = "#3A5A8C",
  neutral_focus = "#333333",
  neutral_main = "#4A4A4A",
  neutral_mid = "#7A7A7A",
  neutral_light = "#C5C5C5",
  neutral_line = "grey55",
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
cat("Supplementary Figure 4 - presentation-only label-lineage audit\n")
cat("Core conclusion: fixed-label reuse, imputation-specific, and strictly fold-wise\n")
cat("analyses are qualitatively concordant near zero beyond Raw-EN.\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("============================================================\n")

auc <- fread(AUC_FILE, encoding = "UTF-8")
summary <- fread(SUMMARY_FILE, encoding = "UTF-8")
dauc <- fread(DAUC_FILE, encoding = "UTF-8")
fold_log <- fread(FOLD_LOG_FILE, encoding = "UTF-8")
locked <- fread(LOCKED_FILE, encoding = "UTF-8", check.names = FALSE)
locked <- locked[clustering == "kmeans_k2_mice"]

# Public aggregate input may already contain the two locked-label display rows.
# Rebuild those rows from the locked authority table below so they are included
# exactly once and retain the same provenance as the original script.
if ("analysis" %in% names(dauc) && any(dauc$analysis == "Locked label")) {
  dauc <- dauc[analysis != "Locked label"]
}

stopifnot(
  all(c("outcome", "outcome_label", "model", "AUC") %in% names(auc)),
  all(c("outcome", "K_folds", "auc_base", "auc_base_phen", "auc_feat",
        "auc_feat_lab", "decisive_dAUC", "dAUC_lo", "dAUC_hi") %in% names(summary)),
  all(c("outcome", "analysis", "d", "lo", "hi", "row_lab", "txt") %in% names(dauc)),
  nrow(auc) == 8L,
  nrow(summary) == 2L,
  nrow(dauc) == 4L,
  nrow(fold_log) == 20L,
  nrow(locked) == 1L,
  all(summary$K_folds == 10L),
  all(fold_log$K_folds == 10L),
  all(fold_log$MICE_M == 5L),
  all(fold_log$MICE_MAXIT == 10L),
  all(as.character(fold_log$QUICK_RUN) == "FALSE")
)

model_levels <- c("Baseline", "Baseline + K2", "Raw33 elastic-net", "Raw33 elastic-net + K2")
model_labels <- c("Base", "Base + K2", "Raw-EN", "Raw-EN + K2")
auc[, model_display := factor(model, levels = model_levels, labels = model_labels)]
auc[, outcome_display := factor(
  outcome_label,
  levels = c("30-day mortality", "MAKE-30"),
  labels = c("30-day mortality", "Renal composite")
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
    plot.margin = margin(5, 8, 5, 7)
  )

pA <- ggplot(auc, aes(AUC, model_display)) +
  geom_point(
    colour = COLORS["neutral_focus"], shape = 16, size = 2.4
  ) +
  geom_text(
    aes(
      x = AUC + 0.005,
      label = sprintf("%.3f", AUC)
    ),
    hjust = 0, size = 2.82, colour = COLORS["neutral_focus"],
    family = BASE_FAMILY
  ) +
  facet_wrap(~ outcome_display, nrow = 1) +
  scale_x_continuous(
    limits = c(0.62, 0.84), breaks = seq(0.65, 0.80, 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Strict fold-wise model discrimination",
    x = "Cross-fitted AUC", y = NULL
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank())

dauc[, analysis_display := fifelse(
  analysis == "MICE Rubin-pooled", "Imputation-specific", "Strict fold-wise"
)]
dauc[, outcome_display := fifelse(
  outcome == "MAKE-30", "Renal composite", "30-day mortality"
)]

locked_dauc <- data.table(
  outcome = c("30-day mortality", "MAKE-30"),
  analysis = "Locked label",
  d = c(locked$mort_dAUC, locked$make_dAUC),
  lo = NA_real_, hi = NA_real_,
  row_lab = NA_character_, txt = NA_character_,
  analysis_display = "Fixed-label reuse",
  outcome_display = c("30-day mortality", "Renal composite")
)
dauc <- rbindlist(list(locked_dauc, dauc), fill = TRUE)
dauc[, row_display := paste(outcome_display, analysis_display, sep = " | ")]
row_levels <- c(
  "Renal composite | Strict fold-wise",
  "Renal composite | Imputation-specific",
  "Renal composite | Fixed-label reuse",
  "30-day mortality | Strict fold-wise",
  "30-day mortality | Imputation-specific",
  "30-day mortality | Fixed-label reuse"
)
dauc[, row_id := match(row_display, row_levels)]
stopifnot(!anyNA(dauc$row_id))
dauc[, `:=`(d_scaled = d * 1e4, lo_scaled = lo * 1e4, hi_scaled = hi * 1e4)]
# One decimal in the display unit retains all five source decimal places.
for (v in c("d", "lo", "hi")) {
  scaled <- dauc[[paste0(v, "_scaled")]]
  stopifnot(isTRUE(all.equal(scaled / 1e4, dauc[[v]], tolerance = 1e-15)),
            isTRUE(all.equal(round(scaled, 1) / 1e4, dauc[[v]], tolerance = 1e-15)))
}
dauc[, method_row := factor(analysis_display, levels = c(
  "Strict fold-wise", "Imputation-specific", "Fixed-label reuse"
))]
dauc[, outcome_display := factor(outcome_display,
  levels = c("30-day mortality", "Renal composite"))]
dauc[, value_label := fifelse(
  analysis_display == "Fixed-label reuse",
  sprintf("%+.1f  [-]", d_scaled),
  sprintf("%+.1f  [%+.1f to %+.1f]", d_scaled, lo_scaled, hi_scaled)
)]
dauc[, uncertainty_status := fifelse(
  is.na(lo) | is.na(hi), "Point estimate only", "Interval available"
)]
stopifnot(
  sum(dauc$uncertainty_status == "Point estimate only") == 2L,
  all(dauc[uncertainty_status == "Point estimate only", analysis_display] == "Fixed-label reuse")
)

analysis_cols <- c(
  "Fixed-label reuse" = unname(COLORS["indigo"]),
  "Imputation-specific" = unname(COLORS["neutral_mid"]),
  "Strict fold-wise" = unname(COLORS["neutral_focus"])
)

pB_forest <- ggplot(dauc, aes(d_scaled, method_row, colour = analysis_display)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(
    data = dauc[!is.na(lo_scaled)],
    aes(xmin = lo_scaled, xmax = hi_scaled),
    orientation = "y", width = 0.16, linewidth = 0.25
  ) +
  geom_point(aes(shape = uncertainty_status), size = 2.3, stroke = .65) +
  scale_colour_manual(values = analysis_cols) +
  scale_shape_manual(values = c("Interval available" = 16, "Point estimate only" = 1)) +
  scale_x_continuous(limits = c(-3.0, 3.6), breaks = seq(-3, 3, 1)) +
  scale_y_discrete(expand = expansion(add = .5)) +
  facet_wrap(~ outcome_display, ncol = 1) +
  labs(
    title = "Incremental AUC: Raw-EN + K2 vs Raw-EN",
    x = "\u0394AUC (10^-4 units)", y = NULL
  ) +
  theme_pub +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
  guides(
    colour = "none",
    shape = guide_legend(title = NULL, nrow = 1, override.aes = list(colour = COLORS["neutral_focus"]))
  )

pB_values <- ggplot(dauc, aes(x = 0, y = method_row, label = value_label)) +
  geom_text(
    hjust = 0, vjust = 0.5, colour = COLORS["neutral_focus"],
    size = FONT_PT["annotation"] / ggplot2::.pt, family = BASE_FAMILY
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = .5)) +
  facet_wrap(~ outcome_display, ncol = 1, labeller = function(x) {
    x[] <- lapply(x, function(v) rep(" ", length(v)))
    x
  }) +
  labs(title = "Estimate [95% CI] (10^-4 units)") +
  theme_void(base_family = BASE_FAMILY) +
  theme(
    plot.title = element_text(
      size = FONT_PT["annotation"], face = "bold", hjust = 0,
      colour = COLORS["neutral_focus"], margin = margin(b = 5)
    ),
    strip.text = element_text(size = FONT_PT["subtitle"], face = "bold",
                              margin = margin(t = 2.2, b = 2.2)),
    strip.background = element_blank(),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(5, 4, 5, 2)
  )

pB_forest <- tag_panel(pB_forest, "B")
pB <- pB_forest + pB_values + plot_layout(widths = c(1.65, 1.05))

pA <- tag_panel(pA, "A")

fig <- pA / pB +
  plot_layout(heights = c(1.04, 1), guides = "keep") +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "D2  Comparator ladder and label-lineage sensitivity")

pA_single <- tag_panel(
  pA +
    facet_wrap(~ outcome_display, ncol = 1) +
    scale_x_continuous(
      limits = c(0.62, 0.86), breaks = seq(0.65, 0.85, 0.05),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = "Fold-wise discrimination") +
    theme(axis.text.x = element_text(size = FONT_PT["axis_text"])),
  "A"
)
pB_single <- ggplot(dauc, aes(d_scaled, method_row, colour = analysis_display)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(
    data = dauc[!is.na(lo_scaled)],
    aes(xmin = lo_scaled, xmax = hi_scaled),
    orientation = "y", width = 0.16, linewidth = 0.25
  ) +
  geom_point(aes(shape = uncertainty_status), size = 2.3, stroke = .65) +
  scale_colour_manual(values = analysis_cols, guide = "none") +
  scale_shape_manual(values = c("Interval available" = 16, "Point estimate only" = 1)) +
  scale_x_continuous(limits = c(-3.0, 3.6), breaks = seq(-3, 3, 1)) +
  scale_y_discrete(expand = expansion(add = .5)) +
  facet_wrap(~ outcome_display, ncol = 1) +
  labs(title = "Raw-EN + K2 vs Raw-EN", x = "\u0394AUC (10^-4 units)", y = NULL) +
  theme_pub +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(title = NULL, nrow = 2, override.aes = list(colour = COLORS["neutral_focus"])))
pB_single <- tag_panel(pB_single, "B")

fig_single <- pA_single / pB_single +
  plot_layout(heights = c(1.08, .92))
fig_single <- add_supplement_header(fig_single, "D2  Label lineage")

fwrite(data.table(input = required, md5 = unname(tools::md5sum(required))),
       file.path(QA_DIR, "layout_1_2_4_6_S4_inputs.csv"))

stem_name <- "Supplementary_Figure_S4_domain2_label_lineage_v5_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
outputs <- save_count_panel_versions(
  fig, stem_name, OUTPUT_DIR, single_plot = fig_single,
  main_width_mm = WIDTH_MM, main_height_mm = HEIGHT_MM,
  single_width_mm = 89, single_height_mm = 205, dpi = DPI
)

qa <- data.table(
  check = c(
    "all_required_inputs_exist", "formal_ten_fold_run", "five_imputations",
    "twenty_completed_fold_logs", "two_outcomes", "four_models_per_outcome",
    "locked_imputation_foldwise_dauc_rows", "all_decisive_estimates_practically_null",
    "two_point_only_rows", "pdf_svg_png_main_and_single_written"
  ),
  pass = c(
    all(file.exists(required)), all(summary$K_folds == 10L), all(fold_log$MICE_M == 5L),
    nrow(fold_log) == 20L, uniqueN(auc$outcome) == 2L,
    all(auc[, .N, by = outcome]$N == 4L), nrow(dauc) == 6L,
    all(abs(dauc$d) < 0.001),
    sum(dauc$uncertainty_status == "Point estimate only") == 2L,
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(QA_DIR, "layout_1_2_4_6_S4_render_QA.csv"))
if (!all(qa$pass)) stop("SF4 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

fwrite(
  dauc[, .(outcome_display, analysis_display, d, lo, hi, d_scaled, lo_scaled,
           hi_scaled, value_label, uncertainty_status)],
  file.path(QA_DIR, "layout_1_2_4_6_S4_display_inverse_check.csv")
)
legend_text <- paste0(
  "Supplementary Figure S4. Comparator ladder and label-lineage sensitivity. ",
  "(A) Strict fold-wise AUC estimates for Base, Base + K2, Raw-EN, and Raw-EN + K2. ",
  "The Raw-EN and Raw-EN + K2 estimates overlap exactly at the displayed precision in the strict fold-wise analysis. ",
  "(B) Delta AUC (Raw-EN + K2 minus Raw-EN) under fixed-label reuse, imputation-specific regeneration, ",
  "and strict fold-wise regeneration, grouped by outcome. The forest axis and numeric column both use units of 10^-4 AUC. ",
  "Filled points have available analysis-specific 95% intervals. Imputation-specific intervals are conditional on computed predictions; ",
  "fold-wise regeneration provides leakage control and is not equivalent to propagating all pipeline uncertainty. ",
  "Open points and a dash in the interval position denote the two fixed-label full-cohort estimates, ",
  "for which no full-pipeline interval was available. The primary comparison is conditional on the full-cohort label vector. ",
  "Raw-EN, elastic-net model containing the source features."
)
writeLines(legend_text, file.path(QA_DIR, "layout_1_2_4_6_S4_caption.txt"), useBytes = TRUE)

cat("SF4 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")
