if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))
suppressPackageStartupMessages(library(ggplot2))

support <- fread(file.path(SIM_TABLE_DIR, "D4_sim_support_metrics_by_K2.csv"))
oc <- fread(file.path(SIM_TABLE_DIR, "D4_sim_operating_characteristics.csv"))

scenario_labels <- c(
  G4L_longitudinal_support_failure = "G4-L: support failure",
  G6L_longitudinal_adequate_support = "G6-L: adequate support",
  D4LQ_ledger_stress_control = "D4-LQ: ledger stress"
)
phenotype_labels <- c(
  "C1 higher-risk" = "C1 higher-risk",
  "C2 lower-risk" = "C2 lower-risk"
)
metric_labels <- c(
  treatment_prevalence = "RRT-first proportion",
  outside_overlap_proportion = "Outside common support",
  ess_ratio = "ATE-IPW ESS / N",
  minimum_cell = "Minimum treatment cell",
  extreme_weight_proportion = "Weights >10"
)

long <- melt(
  support[phenotype %in% names(phenotype_labels)][,
    (names(metric_labels)) := lapply(.SD, as.numeric),
    .SDcols = names(metric_labels)
  ],
  id.vars = c("scenario", "repeat_id", "phenotype"),
  measure.vars = names(metric_labels),
  variable.name = "metric", value.name = "value"
)
long[, scenario_label := factor(
  scenario_labels[scenario], levels = unname(scenario_labels)
)]
long[, phenotype := factor(phenotype, levels = names(phenotype_labels))]
long[, metric_label := factor(metric_labels[metric], levels = unname(metric_labels))]

thresholds <- data.table(
  metric_label = factor(
    unname(metric_labels[c(
      "treatment_prevalence", "outside_overlap_proportion", "ess_ratio",
      "minimum_cell", "extreme_weight_proportion"
    )]),
    levels = unname(metric_labels)
  ),
  threshold = c(0.01, 0.10, 0.25, 30, 0.01)
)

base_theme <- theme_minimal(base_family = "sans", base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0),
    plot.margin = margin(7, 9, 7, 7)
  )

p_support <- ggplot(long, aes(x = scenario_label, y = value, colour = phenotype)) +
  geom_hline(
    data = thresholds, aes(yintercept = threshold),
    inherit.aes = FALSE, linetype = "dashed", colour = "grey55", linewidth = 0.4
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.68), width = 0.54,
    outlier.shape = NA, linewidth = 0.45
  ) +
  stat_summary(
    fun = mean, geom = "point", position = position_dodge(width = 0.68),
    shape = 21, fill = "white", size = 1.7, stroke = 0.45
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(
    "C1 higher-risk" = "#B4553F", "C2 lower-risk" = "#3E6E93"
  )) +
  labs(
    title = "Longitudinal D4 support controls",
    subtitle = "Dashed lines show the locked alert or Gate B thresholds",
    x = NULL, y = NULL, colour = NULL
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

oc_plot <- oc[scenario != "All applicable D4 controls"]
oc_plot[, label := factor(
  scenario_labels[scenario],
  levels = rev(unname(scenario_labels[1:2]))
)]
p_oc <- ggplot(oc_plot, aes(x = rate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70") +
  geom_errorbar(
    aes(xmin = wilson_low, xmax = wilson_high),
    width = 0.12, linewidth = 0.5, colour = "#3A5A8C", orientation = "y"
  ) +
  geom_point(size = 2.1, colour = "#1F3A5F") +
  geom_text(
    aes(label = sprintf("%d/%d", correct, n_valid)),
    hjust = -0.15, size = 3, colour = "#2E2E2E"
  ) +
  scale_x_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Truth-aligned D4 operating characteristics",
    subtitle = "Point estimates and Wilson 95% confidence intervals",
    x = "Correct classification rate", y = NULL
  ) + base_theme

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(SIM_FIGURE_DIR, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(SIM_FIGURE_DIR, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(SIM_FIGURE_DIR, paste0(stem, ".tiff")), plot,
         width = width, height = height, units = "in", dpi = 600,
         compression = "lzw", bg = "white")
}
save_plot(p_support, "Figure_D4L_support_diagnostics", 8.3, 7.2)
save_plot(p_oc, "Figure_D4L_operating_characteristics", 7.2, 3.5)

ademp <- data.table(
  component = c("Aim", "Data-generating mechanisms", "Estimands",
                "Methods", "Performance measures"),
  specification = c(
    "Validate the revised hour-24 D4 support interface under known truth.",
    paste(
      "G4-L phenotype-specific RRT scarcity; G6-L bounded adequate support;",
      "D4-LQ pre-landmark exclusion and exact-tie ledger stress control."
    ),
    paste(
      "Support state, Gate B, phenotype-specific overlap, ESS/N, minimum cell,",
      "and extreme-weight proportion; no treatment-effect estimand."
    ),
    paste(
      "Hour-24 risk set; 24-72-hour first-event ledger; fixed K2;",
      "K2-inclusive propensity model; median/mode covariate imputation."
    ),
    paste(
      "Sensitivity, specificity, Wilson intervals, continuous-metric MCSE,",
      "failure audit, and failure-as-incorrect sensitivity analysis."
    )
  )
)
sim_atomic_fwrite(ademp, file.path(SIM_TABLE_DIR, "D4_sim_ADEMP_table.csv"))
