#!/usr/bin/env Rscript

# Supplementary Figure S12: D3 reconstructed-surrogate controls.
# Presentation only. The display tables come from the confirmatory analysis.

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(arg)) stop("Run with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/"))
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", ".."), winslash = "/")
source(file.path(ROOT, "R", "figures", "00_supplement_theme.R"), encoding = "UTF-8")

SOURCE_DIR <- file.path(ROOT, "data", "figure_source", "supplement")
OUTPUT_DIR <- file.path(ROOT, "figures", "supplement")
AUDIT_DIR <- file.path(OUTPUT_DIR, "render_audit")
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  alerts = file.path(SOURCE_DIR, "SF12_component_alert_rates.csv"),
  metrics = file.path(SOURCE_DIR, "SF12_continuous_metrics.csv"),
  operating = file.path(SOURCE_DIR, "SF12_operating_characteristics.csv")
)
if (!all(file.exists(paths))) stop("Missing S12 source file(s).")

alerts <- fread(paths[["alerts"]])
metrics <- fread(paths[["metrics"]])
operating <- fread(paths[["operating"]])
stopifnot(
  nrow(alerts) == 9L,
  nrow(metrics) == 12L,
  nrow(operating) == 3L,
  all(operating$failed_n == 0L),
  all(operating$planned_n == 100L),
  all(operating$analyzed_n == 100L),
  all(operating$correct_n == 100L),
  setequal(operating$scenario, c(
    "G3R_same_distribution_no_alert",
    "G3R_covariate_drift_alert",
    "D3CQ_calibration_stress_control"
  )),
  setequal(alerts$scenario_label, c(
    "G4a same distribution", "G4b covariate drift", "G4c calibration stress"
  )),
  setequal(alerts$component, c("Prevalence", "Assignment", "Calibration")),
  all(c("scenario", "component", "alert_rate", "scenario_label") %in% names(alerts)),
  all(c("scenario", "metric", "mean", "empirical_p025", "empirical_p975", "scenario_label") %in% names(metrics))
)

scenario_order <- c(
  "G4a same distribution", "G4b covariate drift", "G4c calibration stress"
)
scenario_publication_labels <- c(
  "G4a same distribution" = "Same-distribution baseline",
  "G4b covariate drift" = "Covariate-drift control",
  "G4c calibration stress" = "Calibration-stress control"
)
alerts[, scenario_label := factor(
  unname(scenario_publication_labels[scenario_label]),
  levels = unname(scenario_publication_labels[scenario_order])
)]
metrics[, scenario_label := factor(
  unname(scenario_publication_labels[scenario_label]),
  levels = unname(scenario_publication_labels[scenario_order])
)]
alerts[, component := factor(component, levels = c("Prevalence", "Assignment", "Calibration"))]
baseline_note <- data.table(
  scenario_label = factor(
    "Same-distribution baseline",
    levels = unname(scenario_publication_labels[scenario_order])
  ),
  x = 0.07,
  label = "All components: 0%"
)

pA <- ggplot(alerts, aes(alert_rate, scenario_label, shape = component)) +
  geom_point(size = 2.4, stroke = .6, colour = PAL[["ink"]], fill = "white",
             position = position_dodge(width = 0.60, orientation = "y")) +
  geom_text(
    data = baseline_note,
    aes(x = x, y = scenario_label, label = label),
    inherit.aes = FALSE, hjust = 0,
    size = SUPP_FS[["annotation"]] / ggplot2::.pt,
    family = FONT_FAMILY, colour = PAL[["muted"]]
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_y_discrete(labels = label_wrap_gen(width = 19)) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, .5, 1),
    labels = label_percent(accuracy = 1)
  ) +
  labs(x = "Component alert rate", y = NULL, shape = "Component (A)") +
  theme_supplement() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        panel.grid.major.y = element_blank()) +
  guides(shape = guide_legend(nrow = 1))
pA <- tag_panel(pA, "A", "Component alerts")

prev <- metrics[metric == "prevalence_drift"]
pB <- ggplot(prev, aes(mean, scenario_label)) +
  annotate("rect", xmin = -.10, xmax = .10, ymin = -Inf, ymax = Inf,
           fill = PAL[["neutral_fill"]], colour = NA) +
  geom_vline(
    xintercept = c(-.10, .10), colour = PAL[["reference"]],
    linetype = "dashed", linewidth = .35
  ) +
  geom_errorbar(
    aes(xmin = empirical_p025, xmax = empirical_p975), orientation = "y",
    width = .12, linewidth = .48, colour = PAL[["muted"]]
  ) +
  geom_point(size = 2.35, colour = PAL[["ink"]]) +
  scale_x_continuous(limits = c(-.12, .26), breaks = c(-.10, 0, .10, .20)) +
  labs(x = "C1 prevalence difference", y = NULL) +
  theme_supplement() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
pB <- tag_panel(pB, "B", "Prevalence drift")

cal <- metrics[metric %in% c("calibration_slope", "calibration_intercept")]
make_calibration_panel <- function(dt, limits, breaks, gate, title, show_y = TRUE) {
  ggplot(dt, aes(mean, scenario_label)) +
    annotate(
      "rect", xmin = gate[[1]], xmax = gate[[2]], ymin = -Inf, ymax = Inf,
      fill = PAL[["neutral_fill"]], colour = NA
    ) +
    geom_vline(
      xintercept = gate, colour = PAL[["reference"]],
      linetype = "dashed", linewidth = .35
    ) +
    geom_errorbar(
      aes(xmin = empirical_p025, xmax = empirical_p975), orientation = "y",
      width = .12, linewidth = .48, colour = PAL[["muted"]]
    ) +
    geom_point(size = 2.35, colour = PAL[["ink"]]) +
    facet_wrap(~ metric, labeller = as_labeller(c(
      calibration_slope = "Slope", calibration_intercept = "Intercept"
    ))) +
    scale_y_discrete(labels = label_wrap_gen(width = 19)) +
    scale_x_continuous(limits = limits, breaks = breaks) +
    labs(x = title, y = NULL) +
    theme_supplement() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line(colour = PAL[["muted"]]) else element_blank(),
      axis.line.y = if (show_y) element_line(colour = PAL[["muted"]]) else element_blank()
    )
}

p_slope <- make_calibration_panel(
  cal[metric == "calibration_slope"], c(.70, 1.50), c(.8, 1.0, 1.2, 1.4),
  c(.8, 1.2), "Calibration slope", TRUE
)
p_slope <- tag_panel(p_slope, "C", "Calibration")
p_intercept <- make_calibration_panel(
  cal[metric == "calibration_intercept"], c(-.30, .90), seq(-.2, .8, .2),
  c(-.2, .2), "Calibration intercept", FALSE
) +
  labs(title = " ") +
  theme(
    plot.title = element_text(
      size = SUPP_FS[["panel_title"]], face = "bold", hjust = 0,
      margin = margin(b = 2.5)
    )
  )
pC <- p_slope + p_intercept +
  plot_layout(widths = c(1, 1))

panels <- list(A = pA, B = pB, C_slope = p_slope, C_intercept = p_intercept)
fig <- pA + pB + p_slope + p_intercept + guide_area() +
  plot_layout(design = "AB\nCD\nEE", widths = c(1, 1),
              heights = c(1, 0.92, 0.12), guides = "collect") &
  theme(legend.position = "bottom")
fig <- add_supplement_header(fig, "D3  Reconstructed-surrogate operating characteristics")

stem <- "Supplementary_Figure_S12_D3_reconstructed_generator_validation_v4_submission"
save_supplement_figure(fig, stem, OUTPUT_DIR, width_mm = 180, height_mm = 172)
output_paths <- file.path(OUTPUT_DIR, paste0(stem, c(".pdf", ".png", ".tiff")))

fwrite(alerts, file.path(AUDIT_DIR, "SF12_component_alert_rates.csv"))
fwrite(metrics, file.path(AUDIT_DIR, "SF12_continuous_metrics.csv"))
fwrite(operating, file.path(AUDIT_DIR, "SF12_operating_characteristics.csv"))
legend_text <- paste0(
  "Supplementary Figure S12. Reconstructed-surrogate operating behavior. ",
  "The same-distribution baseline, covariate-drift control, and calibration-stress control evaluate distinct transport alerts. ",
  "Shaded bands are the prespecified non-alert ranges and dashed lines are their boundaries. ",
  "Panel A reports component alert rates; the baseline produced 0% alerts for all three components. ",
  "Panels B and C show means and empirical 2.5th-97.5th percentiles across 100 repeats."
)
writeLines(legend_text, file.path(OUTPUT_DIR, paste0(stem, "_legend.txt")), useBytes = TRUE)

qa <- data.table(
  check = c(
    "source_tables_present", "component_alert_grid_complete", "continuous_metrics_complete",
    "three_control_conditions_complete", "no_failed_repeats", "all_decisions_correct",
    "all_figure_outputs_written", "legend_written"
  ),
  pass = c(
    all(file.exists(paths)), nrow(alerts) == 9L, nrow(metrics) == 12L,
    nrow(operating) == 3L && all(operating$planned_n == 100L) && all(operating$analyzed_n == 100L),
    all(operating$failed_n == 0L), all(operating$correct_n == 100L),
    all(file.exists(output_paths)), file.exists(file.path(OUTPUT_DIR, paste0(stem, "_legend.txt")))
  )
)
fwrite(qa, file.path(OUTPUT_DIR, "Supplementary_Figure_S12_render_QA.csv"))
if (!all(qa$pass)) stop("S12 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))
cat("S12 completed: ", file.path(OUTPUT_DIR, stem), "\n", sep = "")
