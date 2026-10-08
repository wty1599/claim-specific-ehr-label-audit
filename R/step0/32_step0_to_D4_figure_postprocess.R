#!/usr/bin/env Rscript

## WP10-D4 publication-figure postprocess.
## Presentation only: reads locked WP10-D4 outputs and does not refit any model.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
  library(patchwork)
})

AUDIT_ROOT <- file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "domain1_renal_identity_and_method_closure_20260717")
SCRIPT_PATH <- file.path(AUDIT_ROOT, "scripts", "37b_domain1_to_domain4_bridge_figure_postprocess.R")
RUN_ROOT <- file.path(AUDIT_ROOT, "wp10_domain2_domain4_bridge_20260718", "domain4")
TABLE_DIR <- file.path(RUN_ROOT, "tables")
FIGURE_DIR <- file.path(RUN_ROOT, "figures")
PROVENANCE_DIR <- file.path(RUN_ROOT, "provenance")

SUMMARY_FILE <- file.path(TABLE_DIR, "Table_D1_D4_restricted_overlap.csv")
PATIENT_FILE <- file.path(TABLE_DIR, "37_WP10_D4_propensity_patient_level.csv")
CORE_COMPLETION <- file.path(RUN_ROOT, "37_WP10_D4_run_completed_with_declared_missing.ok")

EXPECTED_SHA256 <- c(
  summary = "368B53336B6BB7BF2A9DD2ED62C0D18EBBA257F6317A8EF0ED5297DB3A1937C4",
  patient = "8B4FA860EA9643543F569DE6DE92201094D1658E5C77757E58277BF90119F9B7"
)

FIGURE_WIDTH_MM <- 183
FIGURE_HEIGHT_MM <- 165
PHENOTYPE_COLOURS <- c(
  "C1 higher-risk" = "#B4553F",
  "C2 lower-risk" = "#3E6E93"
)
ALERT_COLOURS <- c(
  "Within gate" = "#3A5A8C",
  "Alert" = "#C77A30"
)
NEUTRAL_POINT <- "#2E2E2E"
REFERENCE_GREY <- "grey70"

cat("WP10-D4 publication figure postprocess\n")
cat("Phenotype colours:", paste(names(PHENOTYPE_COLOURS), PHENOTYPE_COLOURS, collapse = "; "), "\n")
cat("Gate colours:", paste(names(ALERT_COLOURS), ALERT_COLOURS, collapse = "; "), "\n")
cat("Figure size:", FIGURE_WIDTH_MM, "x", FIGURE_HEIGHT_MM, "mm\n")

for (p in c(SCRIPT_PATH, SUMMARY_FILE, PATIENT_FILE, CORE_COMPLETION)) {
  if (!file.exists(p)) stop("Missing required locked input: ", p, call. = FALSE)
}

sha256_file <- function(path) {
  toupper(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

observed_hash <- c(
  summary = sha256_file(SUMMARY_FILE),
  patient = sha256_file(PATIENT_FILE)
)
if (!identical(observed_hash, EXPECTED_SHA256)) {
  stop(
    "Locked WP10-D4 source hashes do not match the audited inputs.\n",
    paste(names(observed_hash), observed_hash, sep = "=", collapse = "\n"),
    call. = FALSE
  )
}

fig_base <- file.path(FIGURE_DIR, "Figure_D1_D4_restricted_overlap_submission_v2")
figure_paths <- paste0(fig_base, c(".pdf", ".png", ".tiff"))
source_path <- file.path(TABLE_DIR, "37b_Figure_D1_D4_restricted_overlap_submission_v2_source_data.csv")
manifest_path <- file.path(PROVENANCE_DIR, "37b_WP10_D4_figure_output_sha256_manifest.csv")
completion_path <- file.path(RUN_ROOT, "37b_WP10_D4_figure_postprocess_completed.ok")
if (any(file.exists(c(figure_paths, source_path, manifest_path, completion_path)))) {
  stop("Versioned WP10-D4 publication outputs already exist; refusing to overwrite.", call. = FALSE)
}

summary_dt <- fread(SUMMARY_FILE)
patient_dt <- fread(PATIENT_FILE)

required_summary <- c(
  "risk_set_mode", "eligibility", "c1_treatment_prevalence",
  "c2_treatment_prevalence", "ps_c_statistic",
  "outside_overlap_proportion", "ess_ratio",
  "minimum_phenotype_treatment_cell", "extreme_weight_proportion"
)
required_patient <- c(
  "stay_id", "risk_set_mode", "eligibility", "phenotype", "treat_24_72", "ps"
)
if (length(setdiff(required_summary, names(summary_dt)))) {
  stop("Summary input is missing required plotting columns.", call. = FALSE)
}
if (length(setdiff(required_patient, names(patient_dt)))) {
  stop("Patient input is missing required plotting columns.", call. = FALSE)
}
if (nrow(summary_dt) != 4L || any(!is.finite(summary_dt$ps_c_statistic))) {
  stop("Locked four-row summary failed structural validation.", call. = FALSE)
}

population_levels <- c(
  "window_start__full", "window_start__KDIGO_2_3",
  "window_end__full", "window_end__KDIGO_2_3"
)
population_labels <- c(
  "24 h\nFull", "24 h\nKDIGO 2-3", "72 h\nFull", "72 h\nKDIGO 2-3"
)

add_population <- function(dt) {
  dt[, population := factor(
    paste(risk_set_mode, eligibility, sep = "__"),
    levels = population_levels,
    labels = population_labels
  )]
  dt
}

phen_long <- melt(
  copy(summary_dt),
  id.vars = c("risk_set_mode", "eligibility"),
  measure.vars = c("c1_treatment_prevalence", "c2_treatment_prevalence"),
  variable.name = "phenotype", value.name = "estimate"
)
phen_long[, phenotype := factor(
  phenotype,
  levels = c("c1_treatment_prevalence", "c2_treatment_prevalence"),
  labels = c("C1 higher-risk", "C2 lower-risk")
)]
phen_long <- add_population(phen_long)

ps_plot <- patient_dt[risk_set_mode == "window_start"]
ps_plot[, phenotype := factor(
  phenotype, levels = c("C1 higher-risk", "C2 lower-risk")
)]
ps_plot[, treatment := factor(
  treat_24_72, levels = c(0, 1), labels = c("No RRT", "RRT 24-72 h")
)]
ps_plot[, eligibility_label := factor(
  eligibility,
  levels = c("full", "KDIGO_2_3"),
  labels = c("Full cohort", "KDIGO stage 2-3")
)]

gate_long <- rbindlist(list(
  summary_dt[, .(
    risk_set_mode, eligibility, metric = "Outside overlap",
    estimate = outside_overlap_proportion, threshold = 0.10,
    alert_ratio = outside_overlap_proportion / 0.10
  )],
  summary_dt[, .(
    risk_set_mode, eligibility, metric = "ESS/N",
    estimate = ess_ratio, threshold = 0.25,
    alert_ratio = 0.25 / ess_ratio
  )],
  summary_dt[, .(
    risk_set_mode, eligibility, metric = "Minimum cell",
    estimate = minimum_phenotype_treatment_cell, threshold = 30,
    alert_ratio = 30 / minimum_phenotype_treatment_cell
  )],
  summary_dt[, .(
    risk_set_mode, eligibility, metric = "Extreme weights",
    estimate = extreme_weight_proportion, threshold = 0.01,
    alert_ratio = extreme_weight_proportion / 0.01
  )]
))
gate_long[alert_ratio == 0, alert_ratio := 1e-3]
gate_long[, gate_status := factor(
  ifelse(alert_ratio > 1, "Alert", "Within gate"),
  levels = c("Within gate", "Alert")
)]
gate_long[, metric := factor(
  metric,
  levels = c("Outside overlap", "ESS/N", "Minimum cell", "Extreme weights")
)]
gate_long <- add_population(gate_long)

cstat_plot <- add_population(summary_dt[, .(
  risk_set_mode, eligibility, ps_c_statistic
)])

figure_source <- rbindlist(list(
  phen_long[, .(
    panel = "A", risk_set_mode, eligibility,
    group = as.character(phenotype), x = as.character(population),
    estimate, threshold = NA_real_, alert_ratio = NA_real_
  )],
  ps_plot[, .(
    panel = "B", risk_set_mode, eligibility,
    group = paste(phenotype, treatment, sep = ";"), x = as.character(stay_id),
    estimate = ps, threshold = NA_real_, alert_ratio = NA_real_
  )],
  gate_long[, .(
    panel = "C", risk_set_mode, eligibility,
    group = paste(metric, gate_status, sep = ";"), x = as.character(population),
    estimate, threshold, alert_ratio
  )],
  cstat_plot[, .(
    panel = "D", risk_set_mode, eligibility,
    group = "PS C-statistic", x = as.character(population),
    estimate = ps_c_statistic, threshold = NA_real_, alert_ratio = NA_real_
  )]
), fill = TRUE)
fwrite(figure_source, source_path, na = "")

theme_pub <- function() {
  theme_classic(base_family = "sans", base_size = 9.5) +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, colour = "grey35", margin = margin(b = 4)),
      plot.tag = element_text(size = 11, face = "bold"),
      strip.text = element_text(size = 8.5, face = "bold"),
      legend.title = element_text(size = 8.5, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(6, "mm"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin = margin(6, 8, 6, 6)
    )
}

p_a <- ggplot(phen_long, aes(population, estimate, fill = phenotype)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = PHENOTYPE_COLOURS, drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Subsequent RRT prevalence",
    subtitle = "First RRT from 24 h to <72 h",
    x = NULL, y = "Treatment prevalence", fill = "Phenotype"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p_b <- ggplot(ps_plot, aes(ps, colour = phenotype, linetype = treatment)) +
  geom_density(linewidth = 0.65, adjust = 1.1) +
  facet_wrap(~eligibility_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = PHENOTYPE_COLOURS, drop = FALSE) +
  scale_linetype_manual(values = c("No RRT" = "solid", "RRT 24-72 h" = "22")) +
  labs(
    title = "Propensity-score distributions",
    subtitle = "Primary 24-h risk set; descriptive curves",
    x = "Estimated treatment propensity", y = "Density",
    colour = "Phenotype", linetype = "Treatment"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    linetype = guide_legend(nrow = 1, byrow = TRUE, order = 2)
  )

p_c <- ggplot(
  gate_long,
  aes(population, alert_ratio, shape = metric, colour = gate_status, group = metric)
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  geom_point(size = 2.1, stroke = 0.65, position = position_dodge(width = 0.62)) +
  scale_shape_manual(values = c(
    "Outside overlap" = 16, "ESS/N" = 17,
    "Minimum cell" = 15, "Extreme weights" = 3
  )) +
  scale_colour_manual(values = ALERT_COLOURS, drop = FALSE) +
  scale_y_log10() +
  labs(
    title = "Prespecified overlap gates",
    subtitle = "Values above 1 trigger the corresponding alert",
    x = NULL, y = "Alert ratio (log scale)",
    shape = "Metric", colour = "Gate status"
  ) +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    shape = guide_legend(nrow = 2, byrow = TRUE, order = 2)
  )

p_d <- ggplot(cstat_plot, aes(population, ps_c_statistic)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = REFERENCE_GREY, linewidth = 0.4) +
  geom_point(size = 2.1, colour = NEUTRAL_POINT) +
  scale_y_continuous(limits = c(0.5, 1), breaks = c(0.5, 0.75, 1)) +
  labs(
    title = "Treatment-model discrimination",
    subtitle = "Higher C-statistic = stronger treatment selection",
    x = NULL, y = "PS C-statistic"
  ) +
  theme_pub() +
  theme(legend.position = "none")

fig <- ((p_a | p_b) / (p_c | p_d)) +
  plot_layout(widths = c(0.96, 1.04), heights = c(1, 1.08), guides = "keep") +
  plot_annotation(tag_levels = "A")

ggsave(
  figure_paths[1], fig,
  width = FIGURE_WIDTH_MM / 25.4, height = FIGURE_HEIGHT_MM / 25.4,
  device = cairo_pdf, family = "sans", bg = "white"
)
ggsave(
  figure_paths[2], fig,
  width = FIGURE_WIDTH_MM / 25.4, height = FIGURE_HEIGHT_MM / 25.4,
  units = "in", dpi = 600, bg = "white"
)
grDevices::tiff(
  figure_paths[3],
  width = FIGURE_WIDTH_MM / 25.4, height = FIGURE_HEIGHT_MM / 25.4,
  units = "in", res = 600, compression = "lzw", type = "cairo", bg = "white"
)
print(fig)
grDevices::dev.off()

if (!all(file.exists(figure_paths)) || any(file.info(figure_paths)$size <= 0)) {
  stop("At least one publication figure export is missing or empty.", call. = FALSE)
}

output_files <- c(figure_paths, source_path)
manifest <- data.table(
  path = normalizePath(output_files, winslash = "/", mustWork = TRUE),
  sha256 = vapply(output_files, sha256_file, character(1)),
  bytes = file.info(output_files)$size
)
fwrite(manifest, manifest_path)
writeLines(c(
  "RUN COMPLETED",
  "work_package=WP10_Domain4_figure_postprocess",
  "statistics_or_models_changed=FALSE",
  paste0("source_summary_sha256=", observed_hash[["summary"]]),
  paste0("source_patient_sha256=", observed_hash[["patient"]]),
  paste0("figure_width_mm=", FIGURE_WIDTH_MM),
  paste0("figure_height_mm=", FIGURE_HEIGHT_MM),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
), completion_path, useBytes = TRUE)

cat("Publication figure postprocess completed without changing statistics.\n")
print(manifest)
