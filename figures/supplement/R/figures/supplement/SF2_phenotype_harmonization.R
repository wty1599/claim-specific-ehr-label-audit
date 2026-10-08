# -*- coding: UTF-8 -*-
## ============================================================================
## Supplementary Figure 2 - phenotype-profile contrasts and cross-database
## feature harmonization.
##
## Presentation-only script. It reads the locked v7 bootstrap profile table and
## the current strict 33-feature external-assignment audit. It does not change labels,
## preprocessing, assignment, bootstrap estimates, or confidence intervals.
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
QA_DIR <- file.path(PROJECT_ROOT, "qa")
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

PROFILE_FILE <- file.path(SOURCE_DIR, "SF2_profile_shift_bootstrap_v7.csv")
HARMONIZATION_FILE <- file.path(SOURCE_DIR, "SF2_feature_harmonization_audit.csv")

required <- c(PROFILE_FILE, HARMONIZATION_FILE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing SF2 input(s): ", paste(missing, collapse = "; "))

WIDTH_MM <- 180
HEIGHT_MM <- 228
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
cat("Supplementary Figure 2 - presentation-only\n")
cat("Palette:\n"); print(COLORS)
cat("Typography (pt):\n"); print(FONT_PT)
cat(sprintf("Canvas: %d x %d mm; raster export: %d dpi\n", WIDTH_MM, HEIGHT_MM, DPI))
cat("Profile source: v7 bootstrap table. Harmonization source: strict 33-feature assignment audit.\n")
cat("============================================================\n")

profile <- fread(PROFILE_FILE, encoding = "UTF-8")
harm <- fread(HARMONIZATION_FILE, encoding = "UTF-8")

stopifnot(
  all(c("cohort", "feature", "feature_label", "estimate", "lower", "upper",
        "n_clusters_total", "bootstrap_B") %in% names(profile)),
  all(c("feature", "eicu_coverage", "used_for_assignment") %in% names(harm)),
  uniqueN(profile$cohort) == 2L,
  unique(profile$bootstrap_B) == 300L,
  nrow(harm) == 33L,
  sum(harm$used_for_assignment, na.rm = TRUE) == 33L
)

feature_display <- c(
  bilirubin_total_max = "Bilirubin (max)", alt_max = "ALT (max)",
  pao2fio2ratio_min = "PaO2/FiO2 (min)", abs_lymphocytes_min = "Absolute lymphocytes (min)",
  lactate_max = "Lactate (max)", ph_min = "pH (min)", pco2_max = "PaCO2 (max)",
  calcium_min = "Calcium (min)", calcium_max = "Calcium (max)", ptt_max = "PTT (max)",
  inr_max = "INR (max)", temperature_min = "Temperature (min)",
  temperature_max = "Temperature (max)", urine_output_24h_ml = "Urine output (24 h)",
  glucose_max = "Glucose (max)", aniongap_max = "Anion gap (max)",
  potassium_min = "Potassium (min)", potassium_max = "Potassium (max)",
  hemoglobin_min = "Hemoglobin (min)", sodium_min = "Sodium (min)",
  sodium_max = "Sodium (max)", wbc_max = "WBC (max)", platelets_min = "Platelets (min)",
  bicarbonate_min = "HCO3- (min)", chloride_min = "Chloride (min)",
  chloride_max = "Chloride (max)", bun_max = "BUN (max)", creatinine_max = "Creatinine (max)",
  resp_rate_max = "Respiratory rate (max)", gcs_min = "GCS (min)",
  spo2_min = "SpO2 (min)", heart_rate_max = "Heart rate (max)", mbp_min = "MAP (min)"
)
stopifnot(
  all(profile$feature %in% names(feature_display)),
  all(harm$feature %in% names(feature_display)),
  !"ptt_max" %in% profile[cohort == "eICU", feature],
  "ptt_max" %in% harm$feature
)
profile[, feature_label_clean := unname(feature_display[feature])]

mimic_order <- profile[cohort == "MIMIC-IV"][order(estimate), feature]
label_levels <- profile[cohort == "MIMIC-IV"][match(mimic_order, feature), feature_label_clean]
profile[, feature_label_clean := factor(feature_label_clean, levels = label_levels)]
profile[, cohort_display := factor(
  cohort,
  levels = c("MIMIC-IV", "eICU"),
  labels = c("MIMIC-IV derivation", "eICU external")
)]

theme_pub <- theme_supplement() +
  theme(
    plot.title = element_text(size = FONT_PT["panel_title"], face = "bold", hjust = 0),
    plot.subtitle = element_text(size = FONT_PT["subtitle"], colour = COLORS["neutral_main"], hjust = 0),
    axis.title = element_text(size = FONT_PT["axis_title"]),
    axis.text = element_text(size = FONT_PT["axis_text"], colour = "black"),
    legend.title = element_text(size = FONT_PT["legend"], face = "bold"),
    legend.text = element_text(size = FONT_PT["legend"]),
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "pt"),
    strip.text = element_text(size = FONT_PT["annotation"], face = "bold"),
    strip.background = element_rect(fill = COLORS["pale_grey"], colour = NA),
    panel.spacing = unit(3.5, "mm"),
    plot.margin = margin(3, 7, 3, 7)
  )

cohort_cols <- c(
  "MIMIC-IV derivation" = unname(COLORS["neutral_focus"]),
  "eICU external" = unname(COLORS["structure"])
)

pA <- ggplot(profile, aes(estimate, feature_label_clean, colour = cohort_display,
                          shape = cohort_display)) +
  geom_vline(xintercept = 0, colour = COLORS["reference"], linewidth = 0.5, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                linewidth = 0.48, position = position_dodge(width = 0.56)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.56)) +
  scale_colour_manual(values = cohort_cols) +
  scale_shape_manual(values = c("MIMIC-IV derivation" = 16, "eICU external" = 1)) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  labs(
    title = "C1-C2 partition-profile contrasts",
    x = "Standardized mean difference (C1 - C2)", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_pub +
  theme(axis.text.y = element_text(size = 8))

harm[, feature_label_clean := unname(feature_display[feature])]
harm[, coverage_plot := eicu_coverage]
harm[, label_factor := factor(feature_label_clean,
                              levels = harm[order(coverage_plot, feature_label_clean), feature_label_clean])]

pB <- ggplot(harm, aes(coverage_plot, label_factor)) +
  geom_segment(
    aes(x = 0, xend = coverage_plot, y = label_factor, yend = label_factor),
    colour = COLORS["grid"], linewidth = 0.25, inherit.aes = FALSE
  ) +
  geom_point(size = 2.15, colour = COLORS["structure"]) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = percent_format()) +
  labs(
    title = "Observed coverage of the 33 assignment features in eICU",
    x = "Observed eICU coverage", y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

pA <- tag_panel(pA, "A")
pB <- tag_panel(pB, "B")

fig <- pA / pB +
  plot_layout(heights = c(24.2, 33.2)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = FONT_PT["panel_tag"], face = "bold", family = BASE_FAMILY),
    plot.tag.position = c(0, 1)
  )
fig <- add_supplement_header(
  fig,
  "Feature contrasts and external data coverage",
  fraction = 0.04
)

fwrite(data.table(
  input = c(PROFILE_FILE, HARMONIZATION_FILE),
  md5 = unname(tools::md5sum(c(PROFILE_FILE, HARMONIZATION_FILE)))
), file.path(QA_DIR, "layout_1_2_4_6_S2_inputs.csv"))

stem_name <- "Supplementary_Figure_S2_profile_harmonization_v4_submission"
stem <- file.path(OUTPUT_DIR, stem_name)
save_supplement_figure(fig, stem_name, OUTPUT_DIR, WIDTH_MM, HEIGHT_MM, DPI)

outputs <- c(paste0(stem, ".pdf"), paste0(stem, ".png"), paste0(stem, ".tiff"))
qa <- data.table(
  check = c(
    "all_required_inputs_exist", "bootstrap_B_is_300", "MIMIC_profile_features_24",
    "eICU_profile_features_23", "assignment_features_33",
    "all_33_used_for_assignment", "all_33_present_in_eICU", "pdf_png_tiff_written"
  ),
  pass = c(
    all(file.exists(required)), unique(profile$bootstrap_B) == 300L,
    uniqueN(profile[cohort == "MIMIC-IV", feature]) == 24L,
    uniqueN(profile[cohort == "eICU", feature]) == 23L,
    nrow(harm) == 33L, sum(harm$used_for_assignment, na.rm = TRUE) == 33L,
    sum(harm$present_in_eicu, na.rm = TRUE) == 33L, all(file.exists(outputs))
  )
)
fwrite(qa, file.path(QA_DIR, "layout_1_2_4_6_S2_render_QA.csv"))
if (!all(qa$pass)) stop("SF2 QA failed: ", paste(qa[pass == FALSE, check], collapse = "; "))

cat("SF2 completed successfully.\n")
cat("Output stem: ", stem, "\n", sep = "")

legend_text <- paste0(
  "Supplementary Figure S2. Partition-profile contrasts and external feature coverage. ",
  "(A) C1-C2 standardized mean differences for the 24 profile variables in MIMIC-IV and eICU, ",
  "with 300-bootstrap 95% intervals; eICU resampling was clustered by patient. ",
  "Filled circles denote MIMIC-IV; open circles denote eICU. ",
  "The eICU PTT profile contrast was unavailable and is therefore not plotted. ",
  "(B) Observed eICU coverage for all 33 surrogate inputs, including PTT, before MIMIC-median completion. ",
  "The eICU labels were assigned by the reconstructed surrogate; these contrasts do not provide ",
  "an independent semantic reference. PTT, partial thromboplastin time."
)
writeLines(legend_text, file.path(QA_DIR, "layout_1_2_4_6_S2_caption.txt"), useBytes = TRUE)
