## =====================================================================
## Publication outputs for D3 regenerated deployable generator v1.
## Presentation only: reads frozen result tables and changes no estimates.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
for (pkg in c(
  "data.table", "ggplot2", "patchwork", "svglite", "ragg",
  "flextable", "officer", "scales"
)) require_namespace(pkg)
library(data.table)
library(ggplot2)
library(patchwork)

## Figure contract
CORE_CONCLUSION <- paste(
  "The regenerated frozen mapper closely reproduces the locked partition;",
  "same-object mortality transport is evaluable and triggers a calibration",
  "alert, while external semantic assignment fidelity lacks a reference."
)
FIGURE_ARCHETYPE <- "quantitative grid"
WIDTH_MM <- 183
HEIGHT_MM <- 150
DPI <- 600

COLORS <- c(
  c1 = "#B4553F",
  c2 = "#3E6E93",
  neutral_dark = "#2E2E2E",
  neutral_mid = "#777777",
  neutral_light = "#D9D9D9",
  indigo = "#3A5A8C",
  alert = "#C77A30",
  reference = "#B3B3B3"
)

cat("D3 generator figure contract\n")
cat("Conclusion:", CORE_CONCLUSION, "\n")
cat("Archetype:", FIGURE_ARCHETYPE, "\n")
cat("Palette:\n")
print(COLORS)
cat("Panel size:", WIDTH_MM, "x", HEIGHT_MM, "mm\n")

theme_pub <- function() {
  theme_classic(base_size = 9.5, base_family = "sans") +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      strip.text = element_text(size = 9.5, face = "bold"),
      plot.title = element_text(size = 9.8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.8, colour = COLORS["neutral_mid"]),
      plot.tag = element_text(size = 11, face = "bold"),
      plot.margin = margin(6, 8, 6, 6),
      panel.grid.major.x = element_line(
        colour = "#ECECEC", linewidth = 0.3
      ),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35)
    )
}

## ---------- Read locked result tables ----------
boot <- fread(file.path(
  OUT_TABLES, "internal_fidelity_bootstrap_summary.csv"
))
cv <- fread(file.path(
  OUT_TABLES, "internal_fidelity_crossfit_summary.csv"
))
hosp <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_hospital_prevalence.csv"
))
prev <- fread(file.path(OUT_TABLES, "d3_generator_v1_prevalence.csv"))
perf <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_model_performance.csv"
))
alerts <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_alert_components.csv"
))
states <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_claim_states.csv"
))
apparent <- fread(file.path(
  OUT_TABLES, "internal_fidelity_summary.csv"
))
crossfit <- fread(file.path(
  OUT_TABLES, "internal_fidelity_crossfit_summary.csv"
))
cluster_outcomes <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_cluster_outcomes.csv"
))
calibration_bootstrap <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_base_gen_calibration_bootstrap.csv"
))

required_nonempty <- list(
  boot = boot, cv = cv, hosp = hosp, prev = prev, perf = perf,
  alerts = alerts, states = states
)
if (any(vapply(required_nonempty, nrow, integer(1)) == 0L)) {
  stop("At least one required D3 result table is empty.")
}

## ---------- Panel A: internal fidelity ----------
A_boot <- boot[metric %in% c("agreement", "ari"), .(
  metric,
  analysis = "Full-cohort apparent",
  estimate,
  low = ci_low,
  high = ci_high
)]
A_cv <- cv[metric %in% c("agreement", "ari"), .(
  metric,
  analysis = "Repeated held-out",
  estimate = mean,
  low = empirical_p025,
  high = empirical_p975
)]
A <- rbind(A_boot, A_cv)
A[, metric_label := fifelse(
  metric == "agreement", "Label agreement", "Adjusted Rand index"
)]
A[, analysis := factor(
  analysis, levels = c("Repeated held-out", "Full-cohort apparent")
)]

pA <- ggplot(A, aes(
  x = estimate, y = analysis, colour = metric_label, shape = metric_label
)) +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    width = 0.12, linewidth = 0.55, orientation = "y"
  ) +
  geom_point(size = 2.2, stroke = 0.7) +
  facet_wrap(~metric_label, ncol = 1) +
  scale_colour_manual(values = c(
    "Adjusted Rand index" = unname(COLORS["indigo"]),
    "Label agreement" = unname(COLORS["neutral_dark"])
  )) +
  scale_shape_manual(values = c(
    "Adjusted Rand index" = 16,
    "Label agreement" = 17
  )) +
  scale_x_continuous(
    limits = c(0.90, 1.00),
    breaks = c(0.90, 0.95, 1.00),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Internal fidelity",
    subtitle = "Locked-label agreement and ARI",
    x = "Fidelity metric",
    y = NULL
  ) +
  guides(colour = "none", shape = "none") +
  theme_pub()

## ---------- Panel B: hospital prevalence ----------
setorder(hosp, c1_prevalence, hospitalid)
hosp[, hospital_rank := seq_len(.N)]
mimic_prev <- prev[cohort == "MIMIC-IV", c1_prevalence]
eicu_prev <- prev[cohort == "eICU", c1_prevalence]

pB <- ggplot(hosp, aes(hospital_rank, 100 * c1_prevalence)) +
  geom_point(
    aes(size = n_stays), colour = COLORS["neutral_mid"], alpha = 0.70
  ) +
  geom_hline(
    yintercept = 100 * mimic_prev, linetype = "dashed",
    colour = COLORS["c1"], linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 100 * eicu_prev, colour = COLORS["c1"], linewidth = 0.8
  ) +
  annotate(
    "text", x = max(hosp$hospital_rank), y = 100 * eicu_prev,
    label = sprintf("eICU pooled %.1f%%", 100 * eicu_prev),
    hjust = 1, vjust = -0.6, size = 2.8, colour = COLORS["c1"]
  ) +
  annotate(
    "text", x = max(hosp$hospital_rank), y = 100 * mimic_prev,
    label = sprintf("MIMIC-IV %.1f%%", 100 * mimic_prev),
    hjust = 1, vjust = 1.4, size = 2.8, colour = COLORS["c1"]
  ) +
  scale_size_continuous(range = c(0.7, 2.4), guide = "none") +
  labs(
    title = "C1 prevalence across eICU hospitals",
    subtitle = "199 hospitals ordered by generated prevalence",
    x = "Hospital rank",
    y = "C1 prevalence (%)"
  ) +
  theme_pub()

## ---------- Panel C: external discrimination ----------
model_labels <- c(
  base = "Clinical baseline",
  base_gen = "Baseline + generated K2",
  feat_pen = "Raw-penalized",
  feat_rf = "Raw-RF"
)
outcome_labels <- c(
  mortality = "Mortality",
  make = "Kidney composite"
)
outcome_facet_labels <- c(
  mortality = "Mortality",
  make = "Kidney\ncomposite"
)
C <- copy(perf)
C[, model_label := factor(
  model_labels[model],
  levels = rev(unname(model_labels))
)]
C[, outcome_label := outcome_facet_labels[outcome]]
C[, model_group := fifelse(
  model == "feat_pen", "Primary raw-feature comparator", "Other comparator"
)]

pC <- ggplot(C, aes(external_auc, model_label)) +
  geom_errorbar(
    aes(xmin = external_auc_ci_low, xmax = external_auc_ci_high),
    width = 0.14, linewidth = 0.55,
    colour = COLORS["neutral_mid"], orientation = "y"
  ) +
  geom_point(
    aes(fill = model_group), shape = 21, size = 2.2,
    colour = COLORS["neutral_dark"], stroke = 0.5
  ) +
  facet_wrap(~outcome_label, nrow = 1) +
  scale_fill_manual(values = c(
    "Primary raw-feature comparator" = unname(COLORS["indigo"]),
    "Other comparator" = "white"
  )) +
  scale_x_continuous(
    limits = c(0.70, 0.83),
    breaks = c(0.70, 0.76, 0.82),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    title = "External AUC",
    subtitle = "Blue = primary Raw-penalized comparator",
    x = "External AUC",
    y = NULL,
    fill = NULL
  ) +
  theme_pub() +
  guides(fill = "none")

## ---------- Panel D: calibration-based transport alerts ----------
D <- alerts[
  component %in% c(
    "calibration_slope", "absolute_calibration_intercept"
  )
]
D[, label := paste(
  fifelse(outcome == "mortality", "Mortality", "Kidney composite"),
  fifelse(
    component == "calibration_slope", "slope", "|intercept|"
  )
)]
D[, status := fifelse(
  outcome == "make",
  "Descriptive only",
  fifelse(alert, "Transport alert", "Within gate")
)]
D[, label := factor(label, levels = rev(c(
  "Mortality slope", "Kidney composite slope",
  "Mortality |intercept|", "Kidney composite |intercept|"
)))]
calibration_plot <- copy(calibration_bootstrap)
calibration_plot[, component := fifelse(
  parameter == "slope",
  "calibration_slope",
  "absolute_calibration_intercept"
)]
calibration_plot[, `:=`(
  plot_ci_low = ci_low,
  plot_ci_high = ci_high
)]
calibration_plot[parameter == "intercept", `:=`(
  plot_ci_low = fifelse(
    ci_low <= 0 & ci_high >= 0,
    0,
    pmin(abs(ci_low), abs(ci_high))
  ),
  plot_ci_high = pmax(abs(ci_low), abs(ci_high))
)]
D <- merge(
  D,
  calibration_plot[, .(outcome, component, plot_ci_low, plot_ci_high)],
  by = c("outcome", "component"),
  all.x = TRUE,
  sort = FALSE
)

pD <- ggplot(D, aes(observed, label)) +
  geom_segment(
    aes(x = lower_gate, xend = upper_gate, yend = label),
    linewidth = 3.8, colour = "#E5E5E5", lineend = "round"
  ) +
  geom_errorbar(
    aes(xmin = plot_ci_low, xmax = plot_ci_high),
    width = 0.12, linewidth = 0.5,
    colour = COLORS["neutral_mid"], orientation = "y"
  ) +
  geom_point(
    aes(colour = status), size = 2.4
  ) +
  geom_text(
    aes(label = sprintf("%.3f", observed)),
    nudge_y = 0.24, size = 2.8, colour = COLORS["neutral_dark"]
  ) +
  scale_colour_manual(values = c(
    "Transport alert" = unname(COLORS["alert"]),
    "Within gate" = unname(COLORS["indigo"]),
    "Descriptive only" = unname(COLORS["neutral_dark"])
  )) +
  scale_x_continuous(
    limits = c(0, 1.30),
    breaks = c(0, 0.2, 0.8, 1.0, 1.2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title = "Calibration transport gates",
    subtitle = "Amber = claim-driving point alert",
    x = "Observed value",
    y = NULL,
    colour = NULL
  ) +
  theme_pub() +
  guides(colour = "none")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "D3 audit of a regenerated deployable surrogate generator",
    subtitle = "Original generator absent | Regenerated surrogate evaluated",
    caption = paste(
      "Intervals use patient-level bootstrap or repeated held-out splits.",
      "Hospital points show all 199 eICU hospitals.\n",
      paste(
        "The kidney composite is descriptive and not endpoint-equivalent",
        "to MAKE30; the regenerated generator did not replace the locked",
        "partition."
      )
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        family = "sans", face = "bold", size = 11, hjust = 0
      ),
      plot.subtitle = element_text(
        family = "sans", size = 8.2, colour = COLORS["neutral_mid"]
      ),
      plot.caption = element_text(
        family = "sans", size = 7.3, colour = COLORS["neutral_mid"],
        hjust = 0
      ),
      plot.tag = element_text(family = "sans", face = "bold", size = 10.5)
    )
  )

base_file <- file.path(OUT_FIGURES, "Figure_D3_generator_v1")
ggsave(
  paste0(base_file, ".pdf"), fig,
  width = WIDTH_MM / 25.4, height = HEIGHT_MM / 25.4,
  device = grDevices::cairo_pdf, family = "sans"
)
ggsave(
  paste0(base_file, ".png"), fig,
  width = WIDTH_MM / 25.4, height = HEIGHT_MM / 25.4,
  dpi = 300, bg = "white"
)
ragg::agg_tiff(
  paste0(base_file, ".tiff"),
  width = WIDTH_MM / 25.4, height = HEIGHT_MM / 25.4,
  units = "in", res = DPI, compression = "lzw", background = "white"
)
print(fig)
grDevices::dev.off()
svglite::svglite(
  paste0(base_file, ".svg"),
  width = WIDTH_MM / 25.4, height = HEIGHT_MM / 25.4,
  bg = "white"
)
print(fig)
grDevices::dev.off()

## ---------- Word tables ----------
make_ft <- function(x, caption) {
  ft <- flextable::flextable(as.data.frame(x))
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::autofit(ft)
  ft <- flextable::set_caption(ft, caption = caption)
  ft
}

write_docx_table <- function(x, caption, path) {
  doc <- officer::read_docx()
  doc <- flextable::body_add_flextable(doc, make_ft(x, caption))
  print(doc, target = path)
}

lineage_table <- states[, .(
  Claim = claim,
  Object = object,
  State = state,
  Interpretation = interpretation
)]
fidelity_table <- rbindlist(list(
  data.table(
    Analysis = "Full-cohort apparent",
    Agreement = apparent$agreement,
    ARI = apparent$ari,
    Agreement_interval = sprintf(
      "%.4f to %.4f",
      boot[metric == "agreement", ci_low],
      boot[metric == "agreement", ci_high]
    ),
    ARI_interval = sprintf(
      "%.4f to %.4f",
      boot[metric == "ari", ci_low],
      boot[metric == "ari", ci_high]
    )
  ),
  data.table(
    Analysis = "Repeated held-out",
    Agreement = crossfit[metric == "agreement", mean],
    ARI = crossfit[metric == "ari", mean],
    Agreement_interval = sprintf(
      "%.4f to %.4f",
      crossfit[metric == "agreement", empirical_p025],
      crossfit[metric == "agreement", empirical_p975]
    ),
    ARI_interval = sprintf(
      "%.4f to %.4f",
      crossfit[metric == "ari", empirical_p025],
      crossfit[metric == "ari", empirical_p975]
    )
  )
))
external_table <- perf[, .(
  Outcome = outcome_labels[outcome],
  Model = model_labels[model],
  External_AUC = round(external_auc, 3),
  AUC_95CI = sprintf(
    "%.3f to %.3f", external_auc_ci_low, external_auc_ci_high
  ),
  Calibration_slope = round(calibration_slope, 3),
  Calibration_intercept = round(calibration_intercept, 3),
  Brier = round(brier_score, 3)
)]
external_table[, `:=`(
  Calibration_slope_95CI = "",
  Calibration_intercept_95CI = ""
)]
for (outcome_name in unique(perf$outcome)) {
  outcome_label <- unname(outcome_labels[outcome_name])
  slope_interval <- calibration_bootstrap[
    outcome == outcome_name & parameter == "slope"
  ]
  intercept_interval <- calibration_bootstrap[
    outcome == outcome_name & parameter == "intercept"
  ]
  external_table[
    Outcome == outcome_label & Model == model_labels["base_gen"],
    `:=`(
      Calibration_slope_95CI = sprintf(
        "%.3f to %.3f", slope_interval$ci_low, slope_interval$ci_high
      ),
      Calibration_intercept_95CI = sprintf(
        "%.3f to %.3f",
        intercept_interval$ci_low,
        intercept_interval$ci_high
      )
    )
  ]
}

write_docx_table(
  lineage_table,
  "D3 claim-specific object lineage and audit states",
  file.path(OUT_TABLES, "Table_D3_generator_lineage.docx")
)
write_docx_table(
  fidelity_table,
  "Internal fidelity of the regenerated deployable surrogate",
  file.path(OUT_TABLES, "Table_D3_generator_fidelity.docx")
)
write_docx_table(
  external_table,
  "External performance of models using the regenerated generator interface",
  file.path(OUT_TABLES, "Table_D3_generator_external.docx")
)

figure_manifest <- data.table(
  panel = c("A", "B", "C", "D"),
  conclusion = c(
    "Regenerated labels closely reproduce the locked partition.",
            "Observed C1 prevalence varies across hospitals; no shrinkage model is implied.",
    "External discrimination is preserved across model families.",
    paste(
      "The same-object mortality model triggers the prespecified",
      "point-estimate calibration alert; its bootstrap interval is displayed."
    )
  ),
  source_file = c(
    "internal_fidelity_bootstrap_summary.csv and internal_fidelity_crossfit_summary.csv",
    "d3_generator_v1_hospital_prevalence.csv and d3_generator_v1_prevalence.csv",
    "d3_generator_v1_model_performance.csv",
    paste(
      "d3_generator_v1_alert_components.csv and",
      "d3_generator_v1_base_gen_calibration_bootstrap.csv"
    )
  ),
  interval = c(
    "patient-level percentile bootstrap or repeated-split empirical interval",
    "hospital-level distribution with pooled reference lines",
    "patient-cluster percentile bootstrap",
    "patient-cluster percentile bootstrap over prespecified point-estimate gates"
  )
)
write_replace_csv(
  figure_manifest,
  file.path(OUT_QC, "Figure_D3_generator_v1_manifest.csv")
)
writeLines(
  c(
    paste0("Core conclusion: ", CORE_CONCLUSION),
    paste0("Archetype: ", FIGURE_ARCHETYPE),
    paste0("Backend: R only"),
    paste0("Final size: ", WIDTH_MM, " x ", HEIGHT_MM, " mm"),
    "Exports: SVG, PDF, 600-dpi LZW TIFF, 300-dpi PNG preview",
    "Panel A: patient bootstrap and repeated held-out empirical intervals",
    "Panel B: all 199 eICU hospitals",
    "Panel C: patient-cluster bootstrap intervals",
    paste(
      "Panel D: patient-cluster bootstrap intervals displayed over",
      "prespecified point-estimate gates"
    ),
    "Visual QA: pending manual preview inspection"
  ),
  file.path(OUT_QC, "Figure_D3_generator_v1_QA.txt"),
  useBytes = TRUE
)

cat("D3 publication outputs rendered.\n")
cat("Figure base:", base_file, "\n")
