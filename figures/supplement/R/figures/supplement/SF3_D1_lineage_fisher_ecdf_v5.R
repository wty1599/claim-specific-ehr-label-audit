#!/usr/bin/env Rscript

# Supplementary Figure S3. Panel A uses a separate numeric column; Panel B displays the exact
# empirical cumulative distribution of the decisive repeat-level Fisher-axis
# BH-adjusted dip-test P value. No analysis is rerun and no value is transformed.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})

SCRIPT_ARG <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(SCRIPT_ARG)) stop("Run this script with Rscript.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", SCRIPT_ARG[[1]]), winslash = "/", mustWork = TRUE)
REPO_ROOT <- normalizePath(file.path(dirname(SCRIPT_PATH), "..", "..", ".."), winslash = "/")
source(file.path(REPO_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")

DIR_FIG <- file.path(REPO_ROOT, "figures", "supplement")
DIR_SOURCE <- file.path(REPO_ROOT, "data", "figure_source", "supplement")
DIR_PROV <- file.path(DIR_FIG, "render_audit")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV, recursive = TRUE, showWarnings = FALSE)

PALETTE <- c(
  neutral_dark = "#2E2E2E",
  neutral_mid = "#8A8A8A",
  neutral_light = "#D8D8D8",
  indigo = "#3A5A8C",
  deep_indigo = "#1F3A5F",
  alert_amber = "#C77A30",
  grid = "grey85",
  reference = "grey70"
)
activate_supplement_palette()
FONT <- c(panel_tag = 11, panel_title = 9.8, axis_title = 9, axis_text = 8,
          annotation = 8, legend = 8, footnote = 8)
BASE_FAMILY <- FONT_FAMILY
DPI <- 600L
DIP_ALPHA <- 0.05
EXPECTED_PANEL_A_SHA256 <- "45f0d7b43a833175f870884af5f8c14f6f6703d38c36f63ffe19f2735aec2c99"

cat("Figure contract: all source values are unchanged; only layout is revised.\n")
cat("Distribution design: exact ECDF steps on the complete 0-1 P-value scale; no smoothing, jitter, or transformation.\n")
cat("Palette:\n"); print(PALETTE)
cat("Font sizes (pt):\n"); print(FONT)
cat("Export: full-width 180 x 145 mm; vector PDF/SVG and 600-dpi PNG.\n")

INPUT <- c(
  panel_a = file.path(DIR_SOURCE, "S3_panelA_source_data.csv"),
  panel_b = file.path(REPO_ROOT, "..", "main", "source_data", "fig23",
                      "Figure3_panelB_fisher_axis_dip_values_v5.csv"),
  count_qa = file.path(DIR_SOURCE, "S3_panelB_count_consistency_v4.csv")
)
if (!all(file.exists(INPUT))) {
  stop("Missing input(s): ", paste(names(INPUT)[!file.exists(INPUT)], collapse = "; "))
}
if (digest(INPUT[["panel_a"]], algo = "sha256", file = TRUE, serialize = FALSE) !=
    EXPECTED_PANEL_A_SHA256) {
  stop("Panel-A source SHA-256 changed; refusing to render.")
}

panel_a <- fread(INPUT[["panel_a"]])
panel_b <- fread(INPUT[["panel_b"]])
count_qa <- fread(INPUT[["count_qa"]])
if (nrow(panel_b) != 900L || anyNA(panel_b$fisher_axis_dip_p_bh6) ||
    !all(count_qa$count_match & count_qa$denominator_match & count_qa$label_match)) {
  stop("Panel-B repeat-level or count-consistency audit failed.")
}
if (!all(panel_b$shape_alert == (panel_b$fisher_axis_dip_p_bh6 < DIP_ALPHA))) {
  stop("Panel-B threshold decisions do not match the plotted continuous values.")
}
if (any(panel_b$fisher_axis_dip_p_bh6 < 0 | panel_b$fisher_axis_dip_p_bh6 > 1)) {
  stop("Panel-B P values fall outside [0, 1].")
}
if (any(panel_b$fisher_axis_dip_p_bh6 == DIP_ALPHA)) {
  stop("A Panel-B P value lies exactly on the strict decision boundary.")
}

# Panel A setup is unchanged from SF3_D1_lineage_scenario_ids.R.
gate <- panel_a[evidence_type == "Reference"]
if (nrow(gate) != 1L || is.na(gate$estimate)) {
  stop("Panel-A must contain exactly one non-missing N1 reference row.")
}
panel_a[, evidence := factor(evidence, levels = c("Locked full-cohort K2", "Lloyd-aligned N1 q95"))]

condition_levels <- c(
  "N1 continuous control", "Empirical Lloyd refits", "G2 discrete control"
)
panel_b[, condition := factor(condition, levels = condition_levels)]
panel_b[, split := factor(split, levels = c("1k", "5k", "10k"))]

condition_colours <- setNames(c(PAL[["reference"]], PAL[["ink"]],
                                SUPP_COLORS[["alert"]]), condition_levels)
condition_linetypes <- c(
  "N1 continuous control" = "dashed",
  "Empirical Lloyd refits" = "solid",
  "G2 discrete control" = "dotdash"
)
condition_display <- c(
  "N1 continuous control" = "Continuous copula\nreference",
  "Empirical Lloyd refits" = "Empirical Lloyd\nrefits",
  "G2 discrete control" = "Outcome-independent\ndiscrete control"
)
legend_data <- data.table(
  split = factor("1k", levels = c("1k", "5k", "10k")),
  condition = factor(condition_levels, levels = condition_levels),
  y = c(0.73, 0.56, 0.39),
  x0 = 0.16,
  x1 = 0.27,
  x_text = 0.30,
  label = unname(condition_display[condition_levels])
)

# Build exact step coordinates from the locked values. Duplicated x=0 or x=1
# coordinates are intentional and preserve boundary point masses.
ecdf_data <- panel_b[, {
  tab <- data.table(p_value = sort(unique(fisher_axis_dip_p_bh6)))
  tab[, cumulative_proportion := vapply(
    p_value,
    function(z) mean(fisher_axis_dip_p_bh6 <= z),
    numeric(1)
  )]
  rbind(
    data.table(p_value = 0, cumulative_proportion = 0),
    tab,
    data.table(p_value = 1, cumulative_proportion = 1)
  )
}, by = .(condition, split)]

gate_qa <- panel_b[, .(
  source_n = .N,
  source_alert_n = sum(fisher_axis_dip_p_bh6 < DIP_ALPHA),
  ecdf_height_at_gate = mean(fisher_axis_dip_p_bh6 < DIP_ALPHA)
), by = .(condition, split)]
gate_qa <- merge(
  gate_qa,
  count_qa[, .(
    condition,
    split,
    archived_n = observed_n,
    archived_alert_n = observed_shape_alert_count,
    archived_count_label = observed_count_label
  )],
  by = c("condition", "split"),
  all.x = TRUE,
  sort = FALSE
)
gate_qa[, `:=`(
  denominator_match = source_n == archived_n,
  count_match = source_alert_n == archived_alert_n,
  height_match = abs(ecdf_height_at_gate - archived_alert_n / archived_n) < 1e-12
)]
if (nrow(gate_qa) != 9L || !all(gate_qa$denominator_match & gate_qa$count_match & gate_qa$height_match)) {
  stop("Panel-B ECDF height at P=0.05 does not reproduce the archived 9-cell count audit.")
}

theme_pub <- function() {
  theme_supplement() +
    theme(
      plot.title = element_text(size = FONT[["panel_title"]], face = "bold", hjust = 0, margin = margin(b = 7)),
      plot.subtitle = element_text(size = FONT[["annotation"]], colour = PALETTE[["neutral_mid"]], margin = margin(b = 7)),
      axis.title = element_text(size = FONT[["axis_title"]]),
      axis.text = element_text(size = FONT[["axis_text"]], colour = PALETTE[["neutral_dark"]]),
      axis.line = element_line(linewidth = 0.16, colour = PAL[["muted"]]),
      axis.ticks = element_line(linewidth = 0.16, colour = PAL[["muted"]]),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 12, 8, 10)
    )
}

# The numeric column is independent of the scientific x scale.
panel_a[, value_label := ifelse(is.na(ci_low), sprintf("%.3f", estimate),
                                sprintf("%.3f [%.3f-%.3f]", estimate, ci_low, ci_high))]
p_a <- ggplot(panel_a, aes(x = estimate, y = evidence, colour = evidence_type)) +
  geom_vline(xintercept = gate$estimate, linetype = "dashed", linewidth = 0.45, colour = PALETTE[["reference"]]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.12, linewidth = 0.25, na.rm = TRUE) +
  geom_point(size = 2.3, stroke = 0.5) +
  scale_colour_manual(values = c("Reference" = PAL[["ink"]], "Separation alert" = PALETTE[["alert_amber"]])) +
  scale_x_continuous(limits = c(2.67, 3.63), breaks = c(2.7, 3.0, 3.3, 3.6), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Pooled-within Mahalanobis\nseparation",
    y = NULL
  ) +
  scale_y_discrete(labels = c(
    "Locked full-cohort K2" = "Fixed full-cohort K2",
    "Lloyd-aligned N1 q95" = "Copula-reference q95"
  )) +
  theme_pub() + theme(plot.margin = margin(6, 8, 5, 8), panel.grid.major.y = element_blank())

p_a_values <- ggplot(panel_a, aes(x = 0, y = evidence, label = value_label)) +
  geom_text(hjust = 0, size = 8 / ggplot2::.pt, family = FONT_FAMILY,
            colour = PAL[["ink"]]) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  labs(title = "Separation [95% CI]", x = " \n ", y = NULL) +
  theme_pub() +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.grid = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        plot.margin = margin(6, 10, 5, 0))

p_b <- ggplot(
  ecdf_data,
  aes(
    x = p_value,
    y = cumulative_proportion,
    colour = condition,
    linetype = condition,
    group = condition
  )
) +
  geom_vline(
    xintercept = DIP_ALPHA, linetype = "dashed", linewidth = 0.45,
    colour = PALETTE[["reference"]]
  ) +
  geom_step(direction = "hv", linewidth = 0.72, lineend = "butt") +
  facet_wrap(vars(split), nrow = 1, scales = "fixed",
             labeller = as_labeller(c(`1k`="1,000",`5k`="5,000",`10k`="10,025"))) +
  scale_colour_manual(values = condition_colours, labels = condition_display, name = NULL) +
  scale_linetype_manual(values = condition_linetypes, labels = condition_display, name = NULL) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(add = c(0.025, 0.025))
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(add = c(0.025, 0.025))
  ) +
  labs(
    title = "Held-out Fisher-axis shape evidence",
    subtitle = "Dashed vertical lines: adjusted P = 0.05",
    x = "BH-adjusted dip-test P value",
    y = "Cumulative proportion of repeats"
  ) +
  theme_pub() +
  theme(
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = FONT[["axis_text"]], face = "bold", hjust = 0,
      colour = PALETTE[["neutral_dark"]], margin = margin(t = 2, b = 2)
    ),
    panel.spacing.x = unit(4, "mm"),
    legend.position = "bottom", legend.key.width = unit(10, "mm"),
    legend.text = element_text(size = 8), legend.margin = margin(t = 4)
  ) + guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

p_a <- tag_panel(p_a, "A", "Separation relative to the copula reference")
p_b <- tag_panel(p_b, "B")

assemble <- function(width_mm, height_mm, single_column = FALSE) {
  header_mm <- 5.7
  panel_a_mm <- 55.0
  body_mm <- height_mm - header_mm
  if (body_mm <= panel_a_mm) stop("Canvas is too short for Panel A.")
  panel_a_plot <- p_a
  panel_b_plot <- p_b
  header_label <- "D1  Separation and distributional shape"
  if (isTRUE(single_column)) {
    # The scientific text is unchanged; line wrapping and the shorter running
    # header prevent clipping at the 86-mm submission width.
    panel_a_plot <- p_a +
      labs(title = "Separation relative to the\ncopula reference") +
      theme(
        plot.title = element_text(
          size = 8, face = "bold", hjust = 0,
          margin = margin(b = 5, l = 12)
        )
      )
    panel_b_plot <- p_b +
      labs(title = "Held-out Fisher-axis\nshape evidence") +
      theme(
        plot.title = element_text(
          size = 8, face = "bold", hjust = 0,
          margin = margin(b = 5, l = 12)
        )
      )
    header_label <- "D1  Diagnostic lineage"
  }
  top <- (panel_a_plot | p_a_values) + plot_layout(widths = c(2.4, 1))
  body <- (free(top) / panel_b_plot) +
    plot_layout(heights = c(panel_a_mm, body_mm - panel_a_mm)) +
    plot_annotation(
      theme = theme(
        text = element_text(family = "sans"),
        plot.tag = element_text(size = FONT[["panel_tag"]], face = "bold", colour = PALETTE[["neutral_dark"]]),
        plot.tag.position = c(0.002, 0.995)
      )
    )
  add_supplement_header(
    body, header_label,
    fraction = header_mm / height_mm
  )
}

single_stem <- "Supplementary_Figure_S3_D1_lineage_ecdf_v5"
fig <- assemble(180, 145)
save_supplement_figure(fig, single_stem, DIR_FIG, 180, 145, DPI)
single_paths <- file.path(DIR_FIG, paste0(single_stem, c(".pdf", ".svg", ".png")))

ecdf_source_path <- file.path(DIR_PROV, "S3_panelB_ecdf_curve_v5.csv")
gate_qa_path <- file.path(DIR_PROV, "S3_panelB_ecdf_gate_consistency_v5.csv")
legend_path <- file.path(DIR_FIG, "Supplementary_Figure_S3_legend_v5.txt")
design_path <- file.path(DIR_FIG, "Supplementary_Figure_S3_design_decision_v5.md")
session_path <- file.path(DIR_PROV, "S3_v5_sessionInfo.txt")
manifest_path <- file.path(DIR_PROV, "S3_v5_output_sha256_manifest.csv")
completion_path <- file.path(DIR_PROV, "S3_v5_completed.ok")

fwrite(ecdf_data, ecdf_source_path)
fwrite(gate_qa, gate_qa_path)

legend <- paste(
  "Supplementary Figure S3. D1 diagnostic lineage and held-out shape evidence.",
  "(A) Pooled-within Mahalanobis separation for the fixed full-cohort K=2 partition",
  "relative to the sample-size-aligned 95th percentile of the empirical-margin",
  "Gaussian-copula reference.",
  "(B) Empirical cumulative distributions across 100 confirmatory repeats of the BH-adjusted",
  "Hartigan dip-test P value computed on the held-out Fisher discriminant score. The Fisher direction",
  "was estimated in training data as the ridge-stabilized pooled-within covariance",
  "inverse multiplied by the difference between the two Lloyd centroids, then fixed",
  "and applied to held-out observations. Multiplicity adjustment was across the Fisher",
  "axis and PC1-PC5; the displayed and decision-relevant value is the adjusted Fisher-axis",
  "P value, rather than the minimum across projections. The vertical dashed line is the",
  "prespecified shape-alert threshold (adjusted P < 0.05).",
  "The continuous control uses an empirical-margin Gaussian-copula reference; empirical Lloyd refits",
  "repeat Lloyd K=2 training and apply the fitted discriminant direction to held-out observations;",
  "the outcome-independent discrete control contains injected",
  "two-component structure. The cumulative proportion at the 0.05 threshold is the alert rate.",
  "Shape alerts occurred in 0/100 continuous-control repeats, 0/100 empirical",
  "Lloyd refits, and 100/100 discrete-control repeats at each evaluation sample size (1,000, 5,000,",
  "and 10,025). All plotted values come from the archived",
  "repeat-level analysis; no statistic was recomputed for presentation.",
  sep = " "
)
writeLines(legend, legend_path, useBytes = TRUE)

writeLines(c(
  "# Supplementary Figure S3 Panel B design decision",
  "",
  "Panel B uses exact empirical cumulative distributions, one for each of the three",
  "conditions at each evaluation sample size. ECDF steps retain every formal repeat",
  "without kernel smoothing, jitter, coordinate truncation, or transformation. This is",
  "preferable to boxplots because the locked values contain point masses at the 0 and 1",
  "boundaries. Colour and line type redundantly distinguish the three conditions in",
  "grayscale. The vertical gate maps the continuous statistic directly to the archived",
  "shape-alert counts.",
  "",
  "The plotted quantity is named precisely: it is the BH-adjusted Hartigan dip-test",
  "P value on the held-out Fisher discriminant score, not an unrecorded Fisher score",
  "magnitude and not the minimum adjusted P value across six projections. The D1-R copula reference does not",
  "follow a near-uniform null distribution in this locked run: its values concentrate",
  "near 1. This is retained as observed. Panel A values and scientific coordinates are",
  "unchanged; labels occupy a separate column on a 180-mm-wide canvas."
), design_path, useBytes = TRUE)

writeLines(capture.output(sessionInfo()), session_path, useBytes = TRUE)
all_paths <- c(single_paths, ecdf_source_path, gate_qa_path, legend_path, design_path,
               INPUT[["panel_b"]], INPUT[["count_qa"]], session_path)
manifest <- data.table(
  path = normalizePath(all_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(all_paths)$size,
  sha256 = vapply(all_paths, digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE)
)
fwrite(manifest, manifest_path)

writeLines(c(
  "status=completed",
  "analysis_rerun=FALSE",
  "panel_a_layout_changed=TRUE",
  "source_rows=900",
  "count_cells_matched=9/9",
  "panel_b_geometry=ecdf_step",
  "p_value_transform=none",
  "canvas_mm=180x145",
  "png_dpi=600",
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("plot_runtime=", R.version.string),
  paste0("script_sha256=", digest(SCRIPT_PATH, algo = "sha256", file = TRUE, serialize = FALSE)),
  paste0("manifest_sha256=", digest(manifest_path, algo = "sha256", file = TRUE, serialize = FALSE))
), completion_path, useBytes = TRUE)

cat("S3 v5 exported without rerunning analysis:\n")
cat(paste0("- ", single_paths), sep = "\n")
cat("\nLegend: ", legend_path, "\n", sep = "")
cat("Design record: ", design_path, "\n", sep = "")
cat("Gate consistency: ", gate_qa_path, "\n", sep = "")
