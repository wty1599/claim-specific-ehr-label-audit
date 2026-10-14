# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 9 - SA-03b subtype separation, PAC, and assignment
## reproducibility.
##
## Presentation-only script. It summarizes formal repeat-level outputs using
## means and empirical 2.5th-97.5th percentiles. Gap statistics are excluded.
## The assignment panel concerns independent same-distribution replicates and
## must not be described as transportability.
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

RECOVERY_FILE <- file.path(SOURCE_DIR, "SF9_SA03b_true_label_recovery_summary.csv")
PAC_FILE <- file.path(SOURCE_DIR, "SF9_SA03b_PAC_summary.csv")
ASSIGN_FILE <- file.path(SOURCE_DIR, "SF9_SA03b_same_distribution_assignment_summary.csv")

required <- c(RECOVERY_FILE, PAC_FILE, ASSIGN_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF9 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 158
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
  K2 = "#2E2E2E",
  K3 = "#666666",
  K4 = "#9A9A9A",
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
cat("Supplementary Figure 9 - presentation-only\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Terminology: same-distribution assignment reproducibility (not transportability).\n")
cat("Gap statistics are intentionally excluded.\n")
cat("============================================================\n")

recovery_summary <- fread(RECOVERY_FILE, encoding = "UTF-8")
pac_summary <- fread(PAC_FILE, encoding = "UTF-8")
assign_summary <- fread(ASSIGN_FILE, encoding = "UTF-8")

delta_levels <- c(0, 0.75, 1.5, 2.25, 3.2, 4)
for (x in list(recovery_summary, pac_summary, assign_summary)) {
  stopifnot(all(c("delta", "K", "mean", "p025", "p975", "n") %in% names(x)))
}

stopifnot(
  all(recovery_summary$n == 100L),
  all(pac_summary$n == 100L),
  all(assign_summary$n == 100L)
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
    plot.margin = margin(5, 7, 5, 7)
  )

k_cols <- c("K=2" = unname(COLORS["K2"]), "K=3" = unname(COLORS["K3"]), "K=4" = unname(COLORS["K4"]))
k_shapes <- c("K=2" = 16, "K=3" = 17, "K=4" = 15)
k_lines <- c("K=2" = "solid", "K=3" = "dashed", "K=4" = "dotted")
legend_labels <- c("K=2" = "K=2 (A-C)", "K=3" = "K=3 (A-C)", "K=4" = "K=4 (B only)")

base_layers <- list(
  scale_x_continuous(breaks = delta_levels, labels = number_format(accuracy = 0.01)),
  scale_colour_manual(
    values = k_cols,
    name = "K",
    labels = legend_labels,
    limits = names(k_cols),
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(
        shape = unname(k_shapes),
        linetype = unname(k_lines)
      ),
      nrow = 1
    )
  ),
  scale_shape_manual(values = k_shapes, limits = names(k_shapes), drop = FALSE, guide = "none"),
  scale_linetype_manual(values = k_lines, limits = names(k_lines), drop = FALSE, guide = "none")
)

pA <- ggplot(recovery_summary, aes(delta, mean, group = K, colour = K, shape = K, linetype = K)) +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.05, linewidth = 0.45) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2.1) +
  base_layers +
  scale_y_continuous(limits = c(-0.03, 1.03), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Latent-structure recovery",
    x = "Injected separation (delta)", y = "ARI versus latent truth"
  ) +
  theme_pub +
  theme(legend.position = "none") + guides(colour = "none")

pB <- ggplot(pac_summary, aes(delta, mean, group = K, colour = K, shape = K, linetype = K)) +
  geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.05, linewidth = 0.45) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2.1) +
  base_layers +
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.7, 0.1)) +
  labs(
    title = "Consensus ambiguity across K",
    x = "Injected separation (delta)", y = "PAC"
  ) +
  theme_pub +
  theme(legend.position = "bottom")

pC <- ggplot(assign_summary, aes(delta, mean, group = K, colour = K, shape = K, linetype = K)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = COLORS["reference"], linewidth = 0.45) +
  geom_errorbar(aes(ymin = p025, ymax = p975), width = 0.05, linewidth = 0.45) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2.1) +
  base_layers +
  scale_y_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Same-distribution assignment reproducibility",
    x = "Injected separation (delta)", y = "Fixed assignment vs reclustering (ARI)"
  ) +
  theme_pub +
  theme(legend.position = "none") + guides(colour = "none")

pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C")

fig <- pA + pB + pC + guide_area() +
  plot_layout(design = "AB\nCC\nDD", heights = c(1, 1.05, 0.12), guides = "collect") +
  plot_annotation() &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(fig, "D1  Separation and same-distribution reproducibility")

RENDER_AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(RENDER_AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(recovery_summary, file.path(RENDER_AUDIT_DIR, "SF9_SA03b_true_label_recovery_summary.csv"))
fwrite(pac_summary, file.path(RENDER_AUDIT_DIR, "SF9_SA03b_PAC_summary.csv"))
fwrite(assign_summary, file.path(RENDER_AUDIT_DIR, "SF9_SA03b_same_distribution_assignment_summary.csv"))
fwrite(data.table(
  input = required,
  md5 = unname(tools::md5sum(required))
), file.path(RENDER_AUDIT_DIR, "SF9_input_manifest.csv"))

stem_name <- "Supplementary_Figure_S9_separation_reproducibility_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
legend_text <- paste(
  "Supplementary Figure S9. Separation and same-distribution reproducibility.",
  "Panels A and B show latent-label recovery and consensus ambiguity across the prespecified injected-separation grid; delta = 3.2 is an intentional design point.",
  "Panel C shows assignment reproducibility between fixed assignment and same-distribution reclustering.",
  "Dashed lines mark the prespecified 0.90 latent-recovery and 0.80 assignment-reproducibility references.",
  "All K=2 series are displayed on the same untruncated scale."
)
writeLines(legend_text, paste0(stem, "_legend.txt"), useBytes = TRUE)
qa <- data.table(
  check = c(
    "all_required_inputs_exist", "one_hundred_repeats_per_cell",
    "six_separation_levels", "recovery_summary_complete", "PAC_summary_complete",
    "assignment_summary_complete", "Gap_not_used",
    "assignment_terminology_same_distribution", "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)),
    all(recovery_summary$n == 100L) && all(pac_summary$n == 100L) && all(assign_summary$n == 100L),
    identical(sort(unique(recovery_summary$delta)), delta_levels), nrow(recovery_summary) == 12L,
    nrow(pac_summary) == 18L, nrow(assign_summary) == 12L,
    !any(grepl("gap", c(names(recovery_summary), names(pac_summary), names(assign_summary)), ignore.case = TRUE)),
    TRUE, all(file.exists(outputs))
  )
)
fwrite(qa, file.path(OUTPUT_DIR, "Supplementary_Figure_S9_render_QA.csv"))
if (!all(qa$pass)) stop("SF9 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("SF9 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")
