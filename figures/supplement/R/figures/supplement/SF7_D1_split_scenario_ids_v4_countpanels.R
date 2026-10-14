#!/usr/bin/env Rscript

# Publication figure for the split-allocation sensitivity results.
# Presentation only: this script reads the source tables and does
# not recompute any statistic or modify source results.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(digest)
  library(svglite)
})

SCRIPT_ARG <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(SCRIPT_ARG)) stop("Run this script with Rscript.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", SCRIPT_ARG[[1]]), winslash = "/", mustWork = TRUE)
REPO_ROOT <- normalizePath(file.path(dirname(SCRIPT_PATH), "..", "..", ".."), winslash = "/")
source(file.path(REPO_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")
source(file.path(dirname(SCRIPT_PATH), "00_count_panel_revision_helpers_v4.R"), encoding = "UTF-8")
DIR_TABLE <- file.path(REPO_ROOT, "data", "figure_source", "supplement")
DIR_FIG <- file.path(REPO_ROOT, "figures", "supplement")
DIR_PROV <- file.path(DIR_FIG, "render_audit")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PROV, recursive = TRUE, showWarnings = FALSE)

COL <- c(
  discrete = "#3A5A8C",
  inconclusive = "#C77A30",
  no_discrete = "#8A8A8A",
  ink = "#2E2E2E",
  line = "grey55",
  reference = "grey70",
  paper = "white"
)
FS <- c(panel = 11, axis_title = 10, axis_text = 9, annotation = 8.5, legend = 8.5, tag = 11)
COL <- c(
  discrete = PAL[["raw_en"]],
  inconclusive = FIGURE_SPEC$states$S1$stroke,
  no_discrete = PAL[["muted"]],
  ink = PAL[["ink"]],
  line = PAL[["muted"]],
  reference = PAL[["reference"]],
  paper = PAL[["paper"]]
)
FS <- c(panel = 9.8, axis_title = 9, axis_text = 8, annotation = 8, legend = 8, tag = 11)
FIG_WIDTH_IN <- 180 / 25.4
FIG_HEIGHT_IN <- 165 / 25.4

cat("Palette:\n")
print(COL)
cat("Typography (pt):\n")
print(FS)
cat(sprintf("Figure size: %.2f x %.2f in; TIFF/PNG 600 dpi\n", FIG_WIDTH_IN, FIG_HEIGHT_IN))

IN <- c(
  gate = file.path(DIR_TABLE, "SF7_gate_split_source.csv"),
  states = file.path(DIR_TABLE, "SF7_state_split_source.csv"),
  empirical = file.path(DIR_TABLE, "SF7_empirical_split_source.csv"),
  transitions = file.path(DIR_TABLE, "SF7_transition_stability_source.csv")
)
if (!all(file.exists(IN))) stop("Missing source files: ", paste(IN[!file.exists(IN)], collapse = ", "))

OUT <- c(
  manifest = file.path(DIR_PROV, "S7_v4_output_sha256_manifest.csv"),
  completion = file.path(DIR_PROV, "S7_v4_completed.ok")
)

gate <- fread(IN[["gate"]])
states <- fread(IN[["states"]])
emp <- fread(IN[["empirical"]])
trans <- fread(IN[["transitions"]])

stopifnot(nrow(gate) == 3L, uniqueN(gate$split_id) == 3L)
stopifnot(nrow(states) == 27L, all(states$n == 100L))
stopifnot(nrow(emp) == 3L, all(emp$final_state == "Inconclusive"))

split_levels <- gate[order(split_order)]$split_id
split_labels <- setNames(
  sprintf("%s /\n%s", comma(gate[order(split_order)]$n_train), comma(gate[order(split_order)]$n_eval)),
  split_levels
)
eval_labels <- setNames(comma(gate[order(split_order)]$n_eval), split_levels)
gate[, split_id := factor(split_id, levels = split_levels)]
states[, split_id := factor(split_id, levels = split_levels)]
emp[, split_id := factor(split_id, levels = split_levels)]

state_levels <- c("No discrete evidence", "Inconclusive", "Discrete evidence")
states[, state_label := factor(state_label, levels = state_levels)]
stage_levels <- c("Continuous copula\nreference", "Outcome-independent\ndiscrete control", "Empirical K2")
states[, stage_label := fifelse(
  stage == "N1_independent_control", "Continuous copula\nreference",
  fifelse(stage == "S7_positive_control", "Outcome-independent\ndiscrete control", "Empirical K2")
)]
states[, stage_label := factor(stage_label, levels = stage_levels)]
states_plot <- states[as.character(stage_label) != "Empirical K2"]
empirical_state_summary <- states[as.character(stage_label) == "Empirical K2"]
states_plot[, stage_label := factor(
  as.character(stage_label), levels = c("Continuous copula\nreference", "Outcome-independent\ndiscrete control")
)]
stopifnot(
  nrow(states_plot) == 18L,
  nrow(empirical_state_summary) == 9L,
  all(empirical_state_summary[state_label == "Inconclusive", count] == 100L),
  all(empirical_state_summary[state_label != "Inconclusive", count] == 0L)
)

theme_pub <- function() {
  theme_supplement() +
    theme(
      text = element_text(family = FONT_FAMILY, colour = COL[["ink"]]),
      plot.title = element_text(size = FS[["panel"]], face = "bold", hjust = 0, margin = margin(b = 5)),
      axis.title = element_text(size = FS[["axis_title"]]),
      axis.text = element_text(size = FS[["axis_text"]], colour = COL[["ink"]]),
      plot.tag = element_text(size = FS[["tag"]], face = "bold"),
      plot.tag.position = "topleft",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = FS[["legend"]]),
      plot.margin = margin(7, 9, 7, 8)
    )
}

p_a <- ggplot(gate, aes(x = split_id, y = q95, group = 1)) +
  geom_hline(
    yintercept = gate[!is.na(cross_stage_n19049_q95)]$cross_stage_n19049_q95,
    colour = COL[["reference"]], linewidth = 0.55, linetype = "dashed"
  ) +
  geom_line(colour = COL[["line"]], linewidth = 0.65) +
  geom_errorbar(
    aes(ymin = q95_bootstrap_ci_low, ymax = q95_bootstrap_ci_high),
    width = 0.10, linewidth = 0.55, colour = COL[["discrete"]]
  ) +
  geom_point(shape = 21, size = 2.5, stroke = 0.65, fill = COL[["paper"]], colour = COL[["discrete"]]) +
  annotate(
    "text", x = 0.65, y = 2.95,
    label = "Reference q95 at n = 19,049", hjust = 0, vjust = 1,
    size = FS[["annotation"]] / ggplot2::.pt, family = FONT_FAMILY, colour = COL[["line"]]
  ) +
  scale_x_discrete(labels = split_labels) +
  scale_y_continuous(limits = c(0, 3.05), breaks = 0:3, expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Copula-reference threshold",
    x = "Training / evaluation sample size",
    y = "Copula-reference q95 threshold"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8.3), plot.margin = margin(8, 8, 7, 10))

states_small <- states_plot[count > 0 & rate < 0.08]
stopifnot(all(states_small$state_label == "Inconclusive"))
p_b <- ggplot(states_plot, aes(x = split_id, y = rate, fill = state_label)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.25) +
  geom_text(
    data = states_plot[count > 0],
    aes(label = fifelse(rate >= 0.08, sprintf("%d", count), "")),
    position = position_stack(vjust = 0.5),
    size = FS[["annotation"]] / ggplot2::.pt, family = FONT_FAMILY, colour = "white", fontface = "bold"
  ) +
  geom_segment(
    data = states_small,
    aes(x = as.numeric(split_id) + 0.34, xend = as.numeric(split_id) + 0.48,
        y = rate / 2, yend = 0.10), inherit.aes = FALSE,
    linewidth = 0.3, colour = COL[["inconclusive"]]
  ) +
  geom_text(
    data = states_small,
    aes(x = as.numeric(split_id) + 0.48, y = 0.12, label = count),
    inherit.aes = FALSE, size = FS[["annotation"]] / ggplot2::.pt,
    family = FONT_FAMILY, colour = COL[["ink"]]
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  facet_wrap(~stage_label, ncol = 2, drop = TRUE) +
  scale_fill_manual(
    values = c(
      "No discrete evidence" = COL[["no_discrete"]],
      "Inconclusive" = COL[["inconclusive"]],
      "Discrete evidence" = COL[["discrete"]]
    ),
    drop = FALSE
  ) +
  scale_x_discrete(labels = eval_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Three-state decisions",
    x = "Evaluation sample size",
    y = "Repeat-level decision rate"
  ) +
  theme_pub() +
  theme(
    strip.text = element_text(size = 8.2, face = "bold", margin = margin(b = 3)),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 8, lineheight = 0.90),
    legend.position = "bottom",
    legend.key.width = grid::unit(10, "pt"),
    panel.spacing.x = grid::unit(8, "pt")
  )

emp_long <- rbindlist(list(
  emp[, .(
    split_id, series = "Empirical Lloyd refits",
    estimate = separation_mean, low = separation_p025, high = separation_p975
  )],
  emp[, .(
    split_id, series = "Continuous-reference q95",
    estimate = q95, low = q95_bootstrap_ci_low, high = q95_bootstrap_ci_high
  )]
))
emp_long[, series := factor(series, levels = c("Continuous-reference q95", "Empirical Lloyd refits"))]

p_c <- ggplot(emp_long, aes(x = split_id, y = estimate, colour = series, shape = series, group = series)) +
  geom_line(position = position_dodge(width = 0.15), linewidth = 0.55) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    width = 0.08, linewidth = 0.5, position = position_dodge(width = 0.15)
  ) +
  geom_point(size = 2.25, stroke = 0.6, position = position_dodge(width = 0.15)) +
  scale_colour_manual(values = c("Continuous-reference q95" = COL[["discrete"]], "Empirical Lloyd refits" = COL[["inconclusive"]])) +
  scale_shape_manual(values = c("Continuous-reference q95" = 21, "Empirical Lloyd refits" = 19)) +
  scale_x_discrete(labels = split_labels) +
  scale_y_continuous(limits = c(2.60, 3.70), breaks = seq(2.6, 3.6, 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Empirical separation and\nreference threshold",
    x = "Training / evaluation sample size",
    y = "Mahalanobis separation"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8.3), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

p_d <- ggplot(emp, aes(x = split_id, y = ari_mean, group = 1)) +
  geom_hline(yintercept = 1, colour = COL[["reference"]], linetype = "dashed", linewidth = 0.5) +
  geom_line(colour = COL[["line"]], linewidth = 0.65) +
  geom_errorbar(aes(ymin = ari_p025, ymax = ari_p975), width = 0.10, linewidth = 0.55, colour = COL[["ink"]]) +
  geom_point(shape = 21, fill = COL[["paper"]], colour = COL[["ink"]], size = 2.5, stroke = 0.65) +
  scale_x_discrete(labels = split_labels) +
  scale_y_continuous(limits = c(0.90, 1.01), breaks = c(0.90, 0.95, 1.00), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Empirical label reproducibility",
    x = "Training / evaluation sample size",
    y = "ARI versus fixed K2 labels"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 8.3))

p_a <- tag_panel(p_a, "A")
p_c <- tag_panel(p_c, "B")
p_b <- tag_panel(p_b, "C")
p_d <- tag_panel(p_d, "D")

fig <- (p_a | p_c) / (p_b | p_d) +
  plot_layout(widths = c(1, 1), heights = c(0.92, 1.08), guides = "keep") +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = FS[["tag"]], colour = COL[["ink"]]),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) &
  theme(legend.position = "bottom")
fig <- add_supplement_header(fig, "D1  Split-allocation sensitivity")

p_a_single <- p_a +
  labs(title = "Copula-reference threshold") +
  theme(plot.tag.position = "topleft", plot.margin = margin(7, 9, 7, 10))
p_c_single <- p_c +
  labs(title = "Empirical separation and\nreference threshold") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(
    legend.position = "bottom", plot.tag.position = "topleft",
    plot.margin = margin(7, 9, 7, 10)
  )
p_b_single <- p_b +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")
p_d_single <- p_d +
  theme(plot.tag.position = "topleft", plot.margin = margin(7, 9, 7, 10))
fig_single <- p_a_single / p_c_single / p_b_single / p_d_single +
  plot_layout(heights = c(.78, .86, 1.08, .86), guides = "keep")
fig_single <- add_supplement_header(fig_single, "D1  Split-allocation sensitivity")

stem <- "Supplementary_Figure_S7_D1_split_sensitivity_v5_submission"
outputs <- save_count_panel_versions(
  fig, stem, DIR_FIG, single_plot = fig_single,
  main_width_mm = 180, main_height_mm = 178,
  single_width_mm = 89, single_height_mm = 228, dpi = 600
)

fwrite(states_plot, file.path(DIR_PROV, "SF7_panelC_three_state_display_v4.csv"))
fwrite(empirical_state_summary, file.path(DIR_PROV, "SF7_panelC_empirical_legend_v4.csv"))

empirical_counts <- empirical_state_summary[state_label == "Inconclusive", count]
legend_text <- paste0(
  "Supplementary Figure S7. Split-allocation sensitivity. ",
  "(A) Bootstrap-calibrated q95 separation thresholds; the horizontal line is the reference calibrated at n = 19,049. ",
  "(B) Empirical Lloyd-refit Mahalanobis separation relative to the corresponding continuous-control threshold. ",
  "(C) Repeat-level three-state decisions for the continuous copula reference and outcome-independent discrete control. ",
  "At evaluation samples of 1k, 5k, and 10k, the continuous control yielded 95/100, 98/100, and 95/100 no-discrete-evidence decisions, ",
  "with the remaining 5/100, 2/100, and 5/100 inconclusive; the discrete control yielded 100/100 discrete-evidence decisions at each size. ",
  "The empirical K2 result is reported here rather than repeated as three constant bars: ",
  paste(empirical_counts, collapse = "/100, "), "/100 inconclusive across the same three evaluation sizes. ",
  "(D) ARI between empirical refits and fixed K2 labels."
)
writeLines(legend_text, file.path(DIR_FIG, paste0(stem, "_legend.txt")), useBytes = TRUE)

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
manifest_targets <- outputs
manifest <- data.table(
  artifact = sub(paste0("^", stem, "_"), "", tools::file_path_sans_ext(basename(manifest_targets))),
  format = tools::file_ext(manifest_targets),
  path = basename(manifest_targets),
  bytes = as.numeric(file.info(manifest_targets)$size),
  sha256 = vapply(manifest_targets, sha256_file, character(1))
)
fwrite(manifest, OUT[["manifest"]])

completion <- c(
  "FIGURE GENERATION COMPLETED",
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("r_version=", R.version.string),
  paste0("script_sha256=", sha256_file(SCRIPT_PATH)),
  paste0("input_bundle_sha256=", digest(vapply(IN, sha256_file, character(1)), algo = "sha256")),
  paste0("output_manifest_sha256=", sha256_file(OUT[["manifest"]])),
  "source_statistics_modified=NO",
  "figure_role=split_allocation_sensitivity",
  "constant_empirical_state_bars_moved_to_legend=YES"
)
writeLines(completion, OUT[["completion"]], useBytes = TRUE)

qa <- data.table(
  check = c(
    "source_tables_present", "three_split_allocations", "one_hundred_repeats_per_control_cell",
    "empirical_refits_inconclusive_at_all_sizes", "all_figure_outputs_written",
    "legend_written", "manifest_written", "completion_record_written"
  ),
  pass = c(
    all(file.exists(IN)), nrow(gate) == 3L && uniqueN(gate$split_id) == 3L,
    nrow(states) == 27L && all(states$n == 100L),
    nrow(emp) == 3L && all(emp$final_state == "Inconclusive"),
    all(file.exists(outputs)),
    file.exists(file.path(DIR_FIG, paste0(stem, "_legend.txt"))),
    file.exists(OUT[["manifest"]]), file.exists(OUT[["completion"]])
  )
)
fwrite(qa, file.path(DIR_FIG, "Supplementary_Figure_S7_v5_render_QA.csv"))
if (!all(qa$pass)) stop("S7 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("37 figure generation completed.\n")
cat("Outputs:\n", paste(manifest_targets, collapse = "\n"), "\n")
