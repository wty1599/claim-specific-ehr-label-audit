suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/"))
REPO_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", ".."), winslash = "/")
source(file.path(REPO_ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")
source(file.path(SCRIPT_DIR, "00_count_panel_revision_helpers_v4.R"), encoding = "UTF-8")

SOURCE_DIR <- file.path(REPO_ROOT, "data", "figure_source", "supplement")
OUTPUT_DIR <- file.path(REPO_ROOT, "figures", "supplement")
QA_DIR <- file.path(REPO_ROOT, "qa")
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

COL <- c(
  dark = PAL_NEUTRAL[["ink"]], mid = PAL_NEUTRAL[["outline"]],
  light = PAL_NEUTRAL[["reference"]], alert_pale = SUPP_COLORS[["alert_pale"]],
  grid = PAL_NEUTRAL[["grid"]], c1 = PAL_ENTITY[["C1"]], c2 = PAL_ENTITY[["C2"]]
)
event_cols <- c(
  "RRT first" = PAL[["ink"]], "Death first" = "#C594B4",
  "ICU exit first" = PAL[["acidbase"]], "Event-free at 72 h" = PAL[["neutral_fill"]]
)

paths <- file.path(SOURCE_DIR, c(
  "SF6_panelA_first_events_v2.csv",
  "SF6_panelB_treatment_v2.csv",
  "SF6_claim_state_v2.csv"
))
assert_files(paths)

events <- fread(paths[1])
support <- fread(paths[2])[phenotype != "Overall"]
state <- fread(paths[3])
stopifnot(
  nrow(events) == 8L,
  nrow(support) == 2L,
  nrow(state) == 1L,
  setequal(events$phenotype, c("C1 higher-risk", "C2 lower-risk")),
  setequal(events$first_event_state, c(
    "RRT_FIRST_24_72", "DEATH_FIRST_24_72",
    "ICU_EXIT_FIRST_24_72", "EVENT_FREE_AT_72"
  )),
  identical(support$n, c(3151L, 14374L)),
  identical(support$n_treated, c(265L, 49L)),
  all(abs(support$treatment_prevalence - support$n_treated / support$n) < 1e-12),
  identical(state$state, "SUPPORT_INADEQUATE"),
  identical(state$primary_gate, "B")
)
event_totals <- events[, .(event_n = sum(n), denominator = unique(phenotype_n)), by = phenotype]
stopifnot(
  all(event_totals$event_n == event_totals$denominator),
  all(abs(events[, sum(proportion), by = phenotype]$V1 - 1) < 1e-12)
)

events[, event_label := factor(
  fcase(
    first_event_state == "RRT_FIRST_24_72", "RRT first",
    first_event_state == "DEATH_FIRST_24_72", "Death first",
    first_event_state == "ICU_EXIT_FIRST_24_72", "ICU exit first",
    default = "Event-free at 72 h"
  ),
  levels = c("RRT first", "Death first", "ICU exit first", "Event-free at 72 h")
)]
events[, phenotype := factor(
  phenotype, levels = c("C1 higher-risk", "C2 lower-risk")
)]
support[, phenotype := factor(
  phenotype, levels = c("C1 higher-risk", "C2 lower-risk")
)]
support[, rate_label := sprintf(
  "RRT first: %.2f%% (%s/%s)", 100 * treatment_prevalence,
  format(n_treated, big.mark = ",", trim = TRUE),
  format(n, big.mark = ",", trim = TRUE)
)]

rrt_event <- events[first_event_state == "RRT_FIRST_24_72", .(
  phenotype, event_n = n, event_denominator = phenotype_n,
  event_rate = n / phenotype_n
)]
rrt_check <- merge(
  rrt_event,
  support[, .(phenotype, n_treated, n, treatment_prevalence)],
  by = "phenotype", sort = FALSE
)
stopifnot(
  all(rrt_check$event_n == rrt_check$n_treated),
  all(rrt_check$event_denominator == rrt_check$n),
  all(abs(rrt_check$event_rate - rrt_check$treatment_prevalence) < 1e-12)
)

pA <- ggplot(events, aes(phenotype, proportion, fill = event_label)) +
  geom_col(width = .62, colour = "white", linewidth = .25) +
  scale_fill_manual(values = event_cols) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1.10), breaks = seq(0, 1, .25),
    expand = expansion(mult = c(0, 0))
  ) +
  geom_text(
    data = support,
    aes(x = phenotype, y = 1.055, label = rate_label),
    inherit.aes = FALSE, family = FONT_FAMILY,
    size = SUPP_FS[["annotation"]] / ggplot2::.pt,
    colour = COL[["dark"]]
  ) +
  labs(x = NULL, y = "First-event composition", fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_supplement() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(5, "mm")
  )
pA <- tag_panel(pA, "A", "First events from ICU hour 24 to before hour 72")

C <- rbindlist(list(
  support[, .(
    phenotype, metric = "Outside overlap",
    ratio = outside_overlap_proportion / .10
  )],
  support[, .(
    phenotype, metric = "ESS/N",
    ratio = .25 / ess_ratio
  )],
  support[, .(
    phenotype, metric = "Minimum treatment cell",
    ratio = 30 / minimum_cell
  )],
  support[, .(
    phenotype, metric = "Extreme weights",
    ratio = pmax(extreme_weight_proportion / .01, 1e-3)
  )]
))
C[, metric := factor(
  metric,
  levels = rev(c(
    "Outside overlap", "ESS/N",
    "Minimum treatment cell", "Extreme weights"
  ))
)]
C[, ratio_label := sprintf("%.2f", ratio)]
C[, `:=`(
  label_ratio = fifelse(ratio >= 1, ratio * 1.12, ratio / 1.12),
  label_hjust = fifelse(ratio >= 1, 0, 1)
)]
pC <- ggplot(C, aes(ratio, metric, colour = phenotype)) +
  annotate(
    "rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
    fill = COL[["alert_pale"]], alpha = .45, colour = NA
  ) +
  geom_vline(
    xintercept = 1, linetype = "dashed",
    colour = COL[["mid"]], linewidth = .45
  ) +
  geom_point(aes(shape = phenotype), position = position_dodge(width = .48), size = 2.35) +
  geom_text(
    aes(x = label_ratio, label = ratio_label, hjust = label_hjust),
    position = position_dodge(width = .48),
    vjust = 0.5, size = SUPP_FS[["annotation"]] / ggplot2::.pt,
    family = FONT_FAMILY, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "C1 higher-risk" = COL[["c1"]],
    "C2 lower-risk" = COL[["c2"]]
  ), breaks = c("C1 higher-risk", "C2 lower-risk")) +
  scale_shape_manual(values = c("C1 higher-risk" = 16, "C2 lower-risk" = 1),
                     breaks = c("C1 higher-risk", "C2 lower-risk")) +
  scale_x_log10(
    limits = c(.08, 25),
    breaks = c(.1, .25, .5, 1, 2, 5, 10, 25),
    labels = label_number(accuracy = .01)
  ) +
  scale_y_discrete(labels = c(
    "Outside overlap" = "Outside overlap", "ESS/N" = "ESS/N",
    "Minimum treatment cell" = "Minimum treatment-group\npatient count",
    "Extreme weights" = "Extreme weights"
  )) +
  labs(x = "Direction-aligned threshold ratio (log scale; >1 = worse)", y = NULL, colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_supplement() +
  theme(legend.position = "bottom")
pC <- tag_panel(pC, "B", "Treatment-support diagnostics")
pB <- pC

fig <- pA / pC +
  plot_layout(heights = c(.92, 1.08), guides = "keep")
fig <- add_supplement_header(fig, "D4  Landmark first events and treatment-support diagnostics")

support[, rate_label_single := sprintf(
  "RRT first: %.2f%%\n%s/%s", 100 * treatment_prevalence,
  format(n_treated, big.mark = ",", trim = TRUE),
  format(n, big.mark = ",", trim = TRUE)
)]
pA_single <- ggplot(events, aes(phenotype, proportion, fill = event_label)) +
  geom_col(width = .62, colour = "white", linewidth = .25) +
  geom_text(
    data = support,
    aes(x = phenotype, y = 1.075, label = rate_label_single),
    inherit.aes = FALSE, family = FONT_FAMILY,
    size = SUPP_FS[["annotation"]] / ggplot2::.pt,
    colour = COL[["dark"]], lineheight = .95
  ) +
  scale_fill_manual(values = event_cols) +
  scale_x_discrete(labels = c("C1 higher-risk" = "C1", "C2 lower-risk" = "C2")) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1), limits = c(0, 1.14),
    breaks = seq(0, 1, .25), expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "First-event composition", fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_supplement() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom", legend.key.width = unit(4, "mm")
  )
pA_single <- tag_panel(pA_single, "A", "First events: ICU 24 to <72 h")

C_single <- copy(C)
C_single[, `:=`(
  label_ratio = fifelse(ratio > 10, ratio / 1.4, ratio * 1.3),
  label_hjust = fifelse(ratio > 10, 1, 0)
)]
pC_single <- pC + C_single +
  scale_x_log10(
    limits = c(.08, 25), breaks = c(.1, 1, 5, 25),
    labels = label_number(accuracy = .01)
  ) +
  scale_y_discrete(labels = c(
    "Outside overlap" = "Outside overlap", "ESS/N" = "ESS/N",
    "Minimum treatment cell" = "Minimum group\nsize (patients)",
    "Extreme weights" = "Extreme weights"
  )) +
  labs(x = "Threshold ratio\n(log scale; >1 = worse)") +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
pC_single <- tag_panel(pC_single, "B", "Support diagnostics")
fig_single <- pA_single / pC_single +
  plot_layout(heights = c(.94, 1.06), guides = "keep")
fig_single <- add_supplement_header(fig_single, "D4  Events and treatment support")

fwrite(C, file.path(QA_DIR, "layout_1_2_4_6_S6_ratio_check.csv"))
fwrite(rrt_check, file.path(QA_DIR, "layout_1_2_4_6_S6_RRT_check.csv"))
fwrite(data.table(input = paths, md5 = unname(tools::md5sum(paths))),
       file.path(QA_DIR, "layout_1_2_4_6_S6_inputs.csv"))

stem <- "Supplementary_Figure_S6_D4_support_v5_submission"
outputs <- save_count_panel_versions(
  fig, stem, OUTPUT_DIR, single_plot = fig_single,
  main_width_mm = 180, main_height_mm = 150,
  single_width_mm = 89, single_height_mm = 185, dpi = 600
)

legend_text <- paste0(
  "Supplementary Figure S6. Landmark first events and treatment-support diagnostics. ",
  "(A) First events from ICU hour 24 to before hour 72 among patients alive, still in the ICU, ",
  "and RRT-free at hour 24. The RRT-first annotations show the numerators and denominators: ",
  paste(support$rate_label, collapse = "; "), ". ",
  "Death-first, ICU-exit-first, and event-free categories are mutually exclusive with RRT-first. ",
  "(B) Direction-aligned ratios to prespecified thresholds: threshold divided by observed value for ESS/N ",
  "and minimum treatment-group patient count, and observed value divided by threshold for outside overlap ",
  "and extreme weights. The minimum treatment-group patient count is the smaller treated or untreated ",
  "group within each phenotype, not a biological cell count. Values greater than 1 indicate an alert. ",
  "The event-count criterion uses 30-day mortality, separately from death-first events in panel A. ",
  "Treatment support was inadequate, and no treatment-effect or interaction model was fitted. ",
  "RRT, renal replacement therapy; ESS, effective sample size."
)
writeLines(legend_text, file.path(QA_DIR, "layout_1_2_4_6_S6_caption.txt"), useBytes = TRUE)

qa <- data.table(
  check = c(
    "RRT_counts_match_event_ledger", "RRT_denominators_match_event_ledger",
    "RRT_rates_match_event_ledger", "all_main_single_outputs_written"
  ),
  pass = c(
    all(rrt_check$event_n == rrt_check$n_treated),
    all(rrt_check$event_denominator == rrt_check$n),
    all(abs(rrt_check$event_rate - rrt_check$treatment_prevalence) < 1e-12),
    all(file.exists(outputs))
  )
)
fwrite(qa, file.path(QA_DIR, "layout_1_2_4_6_S6_render_QA.csv"))
if (!all(qa$pass)) stop("SF6 v4 QA failed.")
