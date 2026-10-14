# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 1 - MICE convergence, distributional diagnostics,
## and across-imputation label stability.
##
## Presentation-only script. It reads locked output CSV files, selects the six
## variables with the highest observed missingness for compact display, and
## does not refit MICE or clustering models. Historical outputs are untouched.
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

CONVERGENCE_FILE <- file.path(SOURCE_DIR, "SF1_selected_convergence_data.csv")
DENSITY_FILE <- file.path(SOURCE_DIR, "SF1_density_curves_public.csv")
MISSINGNESS_FILE <- file.path(SOURCE_DIR, "SF1_missingness_summary_public.csv")
ARI_FILE <- file.path(SOURCE_DIR, "SF1_pairwise_ARI_data.csv")

required <- c(CONVERGENCE_FILE, DENSITY_FILE, MISSINGNESS_FILE, ARI_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF1 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 190
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
cat("Supplementary Figure 1 - presentation-only\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Selection rule: six variables with the greatest observed missingness.\n")
cat("============================================================\n")

conv <- fread(CONVERGENCE_FILE, encoding = "UTF-8")
dens <- fread(DENSITY_FILE, encoding = "UTF-8")
missing_audit <- fread(MISSINGNESS_FILE, encoding = "UTF-8")
ari <- fread(ARI_FILE, encoding = "UTF-8")

stopifnot(
  all(c("variable", "iteration", "imputation", "value", "metric", "variable_label") %in% names(conv)),
  all(c("variable_label", "source_display", "value", "density") %in% names(dens)),
  all(c("variable", "variable_label", "observed_n", "missing_n", "missing_pct") %in% names(missing_audit)),
  "pairwise_ARI" %in% names(ari),
  nrow(ari) == choose(5, 2)
)

N_ANALYSIS <- 20049L
stopifnot(all(missing_audit$observed_n + missing_audit$missing_n == N_ANALYSIS),
          all(abs(missing_audit$missing_pct -
                    100 * missing_audit$missing_n / N_ANALYSIS) < 1e-8))
setorder(missing_audit, -missing_pct, variable)
selected <- missing_audit[1:6, variable]

conv_plot <- conv[variable %chin% selected & metric == "Mean"]
dens_plot <- dens[variable_label %chin% missing_audit[variable %chin% selected,
                                                     variable_label]]

conv_plot[, display_label := factor(
  variable_label,
  levels = missing_audit[variable %chin% selected, variable_label]
)]
dens_plot[, display_label := factor(
  variable_label,
  levels = missing_audit[variable %chin% selected, variable_label]
)]

chain_types <- c(
  "Chain 1" = "solid", "Chain 2" = "22", "Chain 3" = "42",
  "Chain 4" = "F1", "Chain 5" = "12"
)
chain_shapes <- c(
  "Chain 1" = 16, "Chain 2" = 17, "Chain 3" = 15,
  "Chain 4" = 18, "Chain 5" = 3
)

theme_pub <- theme_supplement() +
  theme(
    plot.title = element_text(size = FONT_PT["panel_title"], face = "bold", hjust = 0),
    plot.subtitle = element_text(size = FONT_PT["subtitle"], colour = COLORS["neutral_main"], hjust = 0),
    axis.title = element_text(size = FONT_PT["axis_title"]),
    axis.text = element_text(size = FONT_PT["axis_text"], colour = "black"),
    strip.text = element_text(size = FONT_PT["annotation"], face = "bold", colour = "black"),
    strip.background = element_rect(fill = COLORS["pale_grey"], colour = NA),
    legend.title = element_text(size = FONT_PT["legend"], face = "bold"),
    legend.text = element_text(size = FONT_PT["legend"]),
    legend.position = "bottom",
    legend.key.width = unit(8, "mm"),
    panel.spacing = unit(3.0, "mm"),
    plot.margin = margin(5, 7, 5, 7)
  )

pA <- ggplot(
  conv_plot,
  aes(iteration, value, group = imputation, linetype = imputation, shape = imputation)
) +
  geom_line(colour = COLORS["neutral_main"], linewidth = 0.35) +
  geom_point(colour = COLORS["neutral_focus"], size = 1.0) +
  facet_wrap(~display_label, scales = "free_y", ncol = 3) +
  scale_linetype_manual(values = chain_types) +
  scale_shape_manual(values = chain_shapes) +
  scale_x_continuous(breaks = seq(1, 10, by = 3)) +
  labs(
    title = "MICE chain convergence for the six most incomplete features",
    x = "Iteration", y = "Chain mean (native units)",
    linetype = "Imputation chain", shape = "Imputation chain"
  ) +
  theme_pub

pB <- ggplot(dens_plot, aes(value, density, colour = source_display,
                            linewidth = source_display,
                            linetype = source_display)) +
  geom_line(na.rm = TRUE, key_glyph = draw_key_path) +
  facet_wrap(~display_label, scales = "free", ncol = 3) +
  scale_colour_manual(values = c(
    "Observed" = unname(COLORS["neutral_focus"]),
    "Pooled imputed" = unname(COLORS["structure"])
  )) +
  scale_linewidth_manual(values = c("Observed" = 0.65, "Pooled imputed" = 0.55)) +
  scale_linetype_manual(values = c("Observed" = "solid", "Pooled imputed" = "22")) +
  labs(
    title = "Observed versus pooled-imputed distributions",
    x = "Observed or imputed value (native units)", y = "Density",
    colour = NULL, linewidth = NULL, linetype = NULL
  ) +
  theme_pub

ari_mean <- mean(ari$pairwise_ARI)
ari_range <- range(ari$pairwise_ARI)
ari_summary <- data.table(
  statistic = "Pairwise ARI across five completed imputations",
  n_pairs = nrow(ari), mean = ari_mean,
  minimum = ari_range[1], maximum = ari_range[2], reference = 0.80
)

pA_main <- tag_panel(pA, "A")
pB_main <- tag_panel(pB, "B")
pA_single <- tag_panel(
  pA +
    facet_wrap(~display_label, scales = "free_y", ncol = 2) +
    labs(title = "MICE chain convergence") +
    guides(linetype = guide_legend(nrow = 3, byrow = TRUE)) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.key.width = unit(5, "mm")
    ),
  "A"
)
pB_single <- tag_panel(
  pB +
    facet_wrap(~display_label, scales = "free", ncol = 2) +
    labs(title = "Observed versus imputed distributions") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)),
  "B"
)

fig <- pA_main / pB_main +
  plot_layout(heights = c(1, 1)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(
  fig,
  "Multiple-imputation diagnostics"
)

fig_single <- pA_single / pB_single +
  plot_layout(heights = c(1, 1)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig_single <- add_supplement_header(
  fig_single,
  "Multiple-imputation diagnostics"
)

fwrite(ari_summary, file.path(QA_DIR, "layout_1_2_4_6_S1_ARI_summary.csv"))
fwrite(data.table(
  input = c(CONVERGENCE_FILE, DENSITY_FILE, MISSINGNESS_FILE, ARI_FILE),
  md5 = unname(tools::md5sum(c(CONVERGENCE_FILE, DENSITY_FILE, MISSINGNESS_FILE, ARI_FILE)))
), file.path(QA_DIR, "layout_1_2_4_6_S1_inputs.csv"))

stem_name <- "Supplementary_Figure_S1_MICE_diagnostics_v5_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
outputs <- save_count_panel_versions(
  fig, stem_name, OUTPUT_DIR, single_plot = fig_single,
  main_width_mm = WIDTH_MM, main_height_mm = HEIGHT_MM,
  single_width_mm = 89, single_height_mm = 220, dpi = DPI
)

qa <- data.table(
  check = c(
    "all_required_inputs_exist", "five_chains_present", "ten_iterations_present",
    "six_variables_selected", "ten_pairwise_ARI_values", "all_ARI_at_least_0.80",
    "pairwise_ARI_summary_matches_source", "pdf_svg_png_main_and_single_written"
  ),
  pass = c(
    all(file.exists(required)), uniqueN(conv_plot$imputation) == 5L,
    max(conv_plot$iteration) == 10L, length(selected) == 6L,
    nrow(ari) == 10L, all(ari$pairwise_ARI >= 0.80),
    abs(ari_mean - mean(ari$pairwise_ARI)) < 1e-12 &&
      all(abs(ari_range - range(ari$pairwise_ARI)) < 1e-12),
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(QA_DIR, "layout_1_2_4_6_S1_render_QA.csv"))
if (!all(qa$pass)) stop("SF1 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

legend_text <- paste0(
  "Supplementary Figure S1. Multiple-imputation diagnostics and label stability. ",
  "(A) Chain means across 10 iterations for the six variables with the greatest observed missingness. ",
  "(B) Observed and pooled-imputed distributions for the same variables. ",
  "Solid lines show observed values; dashed lines show pooled-imputed values. ",
  "Lymphocytes denotes the absolute count, not the percentage. ",
  "Across the five completed imputations, the 10 pairwise K=2 label comparisons had mean ARI ",
  sprintf("%.3f (range %.3f-%.3f); reference 0.80. ", ari_mean, ari_range[1], ari_range[2]),
  "The complete value ranges and original density estimation are retained. ",
  "MICE, multiple imputation by chained equations."
)
writeLines(legend_text, file.path(QA_DIR, "layout_1_2_4_6_S1_caption.txt"), useBytes = TRUE)

cat("SF1 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")
