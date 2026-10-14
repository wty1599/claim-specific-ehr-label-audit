# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 11 - SA-06/SA-07 algorithm robustness.
##
## Presentation-only script. SA-06 contrasts plain and consensus k-means.
## SA-07 uses the audited formal GMM/Ward agreement summaries and the optimal-
## rematching postprocessing. The six affected convenience-label rows are not
## read. No clustering or simulation is rerun.
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
OUTPUT_DIR <- file.path(PROJECT_ROOT, "figures", "supplement")
SOURCE_DIR <- file.path(PROJECT_ROOT, "data", "figure_source", "supplement")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

PAC_FILE <- file.path(SOURCE_DIR, "SF11_SA06_PAC_summary.csv")
SELECTED_FILE <- file.path(SOURCE_DIR, "SF11_SA06_selected_K_distribution.csv")
AGREEMENT_FILE <- file.path(SOURCE_DIR, "SF11_SA06_SA07_direct_algorithm_agreement.csv")
RECOVERY_FILE <- file.path(SOURCE_DIR, "SF11_SA06_SA07_latent_subtype_recovery.csv")
INTEGRITY_FILE <- file.path(SOURCE_DIR, "SF11_SA06_SA07_integrity_checks.csv")

required <- c(PAC_FILE, SELECTED_FILE, AGREEMENT_FILE, RECOVERY_FILE, INTEGRITY_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF11 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 194
DPI <- 600
BASE_FAMILY <- "sans"

COLORS <- c(
  structure = "#2E4B6E",
  indigo = "#3A5A8C",
  neutral_focus = "#333333",
  neutral_main = "#4A4A4A",
  neutral_mid = "#737373",
  neutral_light = "#A0A0A0",
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
cat("Supplementary Figure 11 - presentation-only\n")
cat("Core conclusion: strong discrete structure is algorithm independent,\n")
cat("whereas partitions of a severity continuum are algorithm dependent.\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("SA-07 uses optimal rematching; convenience-label rows are excluded.\n")
cat("============================================================\n")

a_data <- fread(PAC_FILE, encoding = "UTF-8")
b_data <- fread(SELECTED_FILE, encoding = "UTF-8")
c_data <- fread(AGREEMENT_FILE, encoding = "UTF-8")
d_data <- fread(RECOVERY_FILE, encoding = "UTF-8")
integrity <- fread(INTEGRITY_FILE, encoding = "UTF-8")

scenario_publication_labels <- c(
  "G1 continuum" = "Continuous generator",
  "L1 legacy three-component" = "Three-component discrete control",
  "L2 legacy drift" = "Transport-drift control"
)
for (dt in list(a_data, b_data, c_data, d_data)) {
  stopifnot(all(as.character(dt$scenario_display) %in% names(scenario_publication_labels)))
  dt[, scenario_display := unname(scenario_publication_labels[as.character(scenario_display)])]
}

a_data[, scenario_display := factor(scenario_display, levels = unique(scenario_display))]
a_data[, K_display := factor(K_display, levels = c("K=2", "K=3", "K=4"))]
b_data[, scenario_display := factor(scenario_display, levels = unique(scenario_display))]
b_data[, method_display := factor(method_display, levels = c("Minimum PAC", "Silhouette", "PAC delta area"))]
b_data[, selected_K_display := factor(selected_K_display, levels = c("K=2", "K=3", "K=4"))]
b_small <- b_data[proportion > 0 & proportion < 0.05]
b_small[, small_label := scales::percent(proportion, accuracy = 1)]
c_data[, scenario_display := factor(scenario_display, levels = unique(scenario_display))]
c_data[, K_display := factor(K_display, levels = c("K=2", "K=3"))]
c_data[, comparison := factor(comparison, levels = c("Plain vs consensus", "GMM vs Ward"))]
d_data[, scenario_display := factor(scenario_display, levels = unique(scenario_display))]
d_data[, K_display := factor(K_display, levels = c("K=2", "K=3"))]
d_data[, algorithm_display := fifelse(
  algorithm_display == "Ward.D2", "Ward linkage", algorithm_display
)]
d_data[, algorithm_display := factor(
  algorithm_display,
  levels = c("Plain k-means", "Consensus k-means", "GMM", "Ward linkage")
)]

stopifnot(
  nrow(a_data) == 9L,
  nrow(b_data) == 12L,
  nrow(c_data) == 12L,
  nrow(d_data) == 24L,
  all(a_data$n == 100L),
  all(b_data[, sum(n), by = .(scenario, method)]$V1 == 100L),
  all(c_data$n == 100L), all(d_data$n == 100L)
)

k_cols <- c("K=2" = unname(COLORS["neutral_focus"]), "K=3" = unname(COLORS["neutral_mid"]), "K=4" = unname(COLORS["neutral_light"]))
k_shapes <- c("K=2" = 16, "K=3" = 17, "K=4" = 15)
algorithm_cols <- c(
  "Plain k-means" = unname(COLORS["neutral_focus"]),
  "Consensus k-means" = unname(COLORS["indigo"]),
  "GMM" = unname(COLORS["neutral_mid"]),
  "Ward linkage" = unname(COLORS["neutral_light"])
)
algorithm_shapes <- c("Plain k-means" = 16, "Consensus k-means" = 17, "GMM" = 15, "Ward linkage" = 18)
scenario_x_labels <- c(
  "Continuous generator" = "Continuous\ngenerator",
  "Three-component discrete control" = "Discrete\ncontrol",
  "Transport-drift control" = "Transport-drift\ncontrol"
)
scenario_y_labels <- c(
  "Continuous generator" = "Continuous\ngenerator",
  "Three-component discrete control" = "Discrete\ncontrol",
  "Transport-drift control" = "Transport-drift\ncontrol"
)

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

pd_k <- position_dodge(width = 0.5)

pA <- ggplot(a_data, aes(scenario_display, mean, colour = K_display, shape = K_display)) +
  geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.14, linewidth = 0.48, position = pd_k) +
  geom_point(size = 2.25, position = pd_k) +
  scale_colour_manual(values = k_cols) +
  scale_shape_manual(values = k_shapes) +
  scale_x_discrete(labels = scenario_x_labels) +
  scale_y_continuous(limits = c(-0.02, 0.57), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = "PAC ambiguity",
    x = NULL, y = "PAC"
  ) +
  theme_pub +
  labs(colour = "K (A)", shape = "K (A)") +
  theme(legend.position = "bottom", legend.title = element_text(size = FONT_PT["legend"], face = "bold"))

pB <- ggplot(b_data, aes(proportion, method_display, fill = selected_K_display)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.22) +
  geom_text(
    aes(label = fifelse(proportion >= 0.05, proportion_label, "")),
    position = position_stack(vjust = 0.5),
    colour = "white", size = FONT_PT["annotation"] / ggplot2::.pt, family = BASE_FAMILY
  ) +
  geom_segment(
    data = b_small,
    aes(x = proportion / 2, xend = 0.12,
        y = as.numeric(method_display) - 0.35, yend = as.numeric(method_display) - 0.57),
    inherit.aes = FALSE, linewidth = 0.3, colour = COLORS["neutral_main"]
  ) +
  geom_text(
    data = b_small,
    aes(x = 0.14, y = as.numeric(method_display) - 0.57, label = small_label),
    inherit.aes = FALSE, hjust = 0, colour = COLORS["neutral_focus"],
    size = FONT_PT["annotation"] / ggplot2::.pt, family = BASE_FAMILY
  ) +
  facet_wrap(
    ~ scenario_display, ncol = 1, scales = "free_y",
    labeller = labeller(scenario_display = label_wrap_gen(width = 25))
  ) +
  scale_fill_manual(values = k_cols, drop = FALSE) +
  scale_y_discrete(expand = expansion(add = c(0.85, 0.55))) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Selected K across criteria",
    x = "Proportion of repeats", y = NULL, fill = "Selected K (B)"
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = FONT_PT["legend"], face = "bold"),
    panel.spacing.y = grid::unit(1.8, "mm"),
    strip.text.x = element_text(size = FONT_PT["annotation"], margin = margin(2, 2, 2, 2))
  )

pC <- ggplot(c_data, aes(mean, scenario_display, colour = K_display, shape = K_display)) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(aes(xmin = p025, xmax = p975), orientation = "y", width = 0.14, linewidth = 0.48, position = pd_k) +
  geom_point(size = 2.25, position = pd_k) +
  facet_wrap(~ comparison, nrow = 1, labeller = as_labeller(c(
    "Plain vs consensus" = "Plain vs\nconsensus",
    "GMM vs Ward" = "GMM vs\nWard"
  ))) +
  scale_colour_manual(values = k_cols[c("K=2", "K=3")]) +
  scale_shape_manual(values = k_shapes[c("K=2", "K=3")]) +
  scale_y_discrete(labels = scenario_y_labels) +
  scale_x_continuous(limits = c(-0.03, 1.03), breaks = c(0, 0.5, 1.0)) +
  labs(
    title = "Algorithm agreement",
    x = "Agreement (ARI)", y = NULL
  ) +
  theme_pub +
  labs(colour = "K (C)", shape = "K (C)") +
  theme(legend.position = "bottom", legend.title = element_text(size = FONT_PT["legend"], face = "bold"))

d_plot <- d_data[as.character(scenario_display) != "Continuous generator"]
no_truth <- data.table(
  scenario_display = factor("Continuous generator", levels = levels(d_data$scenario_display)),
  K_display = factor(c("K=2", "K=3"), levels = levels(d_data$K_display)),
  x = 0.50,
  label = "No latent\nsubtype truth"
)

pD <- ggplot(d_plot, aes(mean, scenario_display, colour = algorithm_display, shape = algorithm_display)) +
  annotate("segment", x = 0.8, xend = 0.8, y = 1.5, yend = 3.5,
           linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = p025, xmax = p975), orientation = "y", width = 0.12, linewidth = 0.45,
    position = position_dodge(width = 0.78)
  ) +
  geom_point(size = 2.05, position = position_dodge(width = 0.78)) +
  geom_text(
    data = no_truth,
    aes(x = x, y = scenario_display, label = label),
    inherit.aes = FALSE, lineheight = 0.92,
    size = FONT_PT["annotation"] / ggplot2::.pt,
    colour = COLORS["neutral_main"], family = BASE_FAMILY
  ) +
  facet_wrap(~ K_display, nrow = 1) +
  scale_colour_manual(values = algorithm_cols) +
  scale_shape_manual(values = algorithm_shapes) +
  scale_y_discrete(labels = scenario_y_labels, drop = FALSE) +
  scale_x_continuous(limits = c(-0.03, 1.03), breaks = c(0, 0.5, 1.0)) +
  labs(
    title = "Latent-label recovery",
    x = "ARI versus latent subtype", y = NULL,
    colour = "Algorithm (D)", shape = "Algorithm (D)"
  ) +
  theme_pub +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.title.position = "top",
        legend.title = element_text(size = FONT_PT["legend"], face = "bold"),
        legend.key.width = grid::unit(3, "mm"), legend.spacing.x = grid::unit(1, "mm"))

pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C")
pD <- tag_panel(pD, "D")

fig <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "keep", heights = c(1.12, 1)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "Sensitivity  Cluster-number and algorithm robustness")

RENDER_AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(RENDER_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(a_data, file.path(RENDER_AUDIT_DIR, "SF11_SA06_PAC_summary.csv"))
fwrite(b_data, file.path(RENDER_AUDIT_DIR, "SF11_SA06_selected_K_distribution.csv"))
fwrite(c_data, file.path(RENDER_AUDIT_DIR, "SF11_SA06_SA07_direct_algorithm_agreement.csv"))
fwrite(d_data, file.path(RENDER_AUDIT_DIR, "SF11_SA06_SA07_latent_subtype_recovery.csv"))
fwrite(integrity, file.path(RENDER_AUDIT_DIR, "SF11_SA06_SA07_integrity_checks.csv"))
fwrite(data.table(input = required, md5 = unname(tools::md5sum(required))),
       file.path(RENDER_AUDIT_DIR, "SF11_input_manifest.csv"))

stem_name <- "Supplementary_Figure_S11_algorithm_robustness_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
qa <- data.table(
  check = c(
    "all_required_inputs_exist", "all_integrity_checks_passed",
    "six_convenience_label_rows_fixed_and_excluded",
    "formal_300_tasks_SA06", "formal_300_tasks_SA07", "no_task_failures_SA06",
    "no_repeat_failures_SA07", "SA07_rematching_validation_complete",
    "one_hundred_repeats_each", "S2_S3_K3_PAC_zero",
    "S2_S3_minimum_PAC_K3_rates_22_24_percent",
    "S2_S3_silhouette_and_delta_area_select_K3",
    "S2_S3_K3_direct_agreement_one", "S2_S3_K3_latent_recovery_one",
    "continuous_generator_has_no_latent_truth", "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)), all(integrity$resolved_pass),
    integrity[
      analysis == 7L & check == "affected_presaved_matched_label_rows",
      observed == 6 & passed == "FIXED_BY_POSTPROCESSING"
    ],
    integrity[analysis == 6L & check == "formal_unique_tasks", observed] == 300,
    integrity[analysis == 7L & check == "formal_unique_tasks", observed] == 300,
    integrity[analysis == 6L & check == "task_failure_rows", observed] == 0,
    integrity[analysis == 7L & check == "repeat_failure_rows", observed] == 0,
    integrity[analysis == 7L & check == "direct_agreement_validation_rows", observed] == 1800,
    all(a_data$n == 100L) && all(c_data$n == 100L) && all(d_data$n == 100L) &&
      all(b_data[, sum(n), by = .(scenario, method)]$V1 == 100L),
    all(a_data[scenario %in% c("S2_true_discrete_subtypes", "S3_transport_drift") & k == 3, mean] == 0),
    b_data[scenario == "S2_true_discrete_subtypes" & method == "consensus_PAC_minimum" & selected_k == 3, proportion] == 0.22 &&
      b_data[scenario == "S3_transport_drift" & method == "consensus_PAC_minimum" & selected_k == 3, proportion] == 0.24,
    all(b_data[
      scenario %in% c("S2_true_discrete_subtypes", "S3_transport_drift") &
        method %in% c("plain_kmeans_silhouette", "PAC_delta_area") & selected_k == 3,
      proportion
    ] == 1),
    all(c_data[scenario %in% c("S2_true_discrete_subtypes", "S3_transport_drift") & k == 3, mean] >= 0.9999),
    all(d_data[scenario %in% c("S2_true_discrete_subtypes", "S3_transport_drift") & k == 3, mean] >= 0.9999),
    all(abs(d_data[scenario == "S1_pure_severity_continuum", mean]) < 0.001),
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(OUTPUT_DIR, "Supplementary_Figure_S11_render_QA.csv"))
if (!all(qa$pass)) stop("SF11 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("SF11 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")

legend_text <- paste0(
  "Supplementary Figure S11. Cluster-number and algorithm robustness. ",
  "Panels A-C summarize consensus ambiguity, selected K, and direct agreement among clustering procedures. ",
  "Panel D compares recovered labels with planted latent labels only for the two generators that contain latent classes. ",
  "The continuous generator has no latent subtype truth, so no recovery ARI is defined or plotted for that condition. ",
  "Ward linkage denotes Ward's minimum-variance hierarchical method."
)
writeLines(legend_text, file.path(OUTPUT_DIR, paste0(stem_name, "_legend.txt")), useBytes = TRUE)
