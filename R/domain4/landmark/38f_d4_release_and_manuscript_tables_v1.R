source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

suppressPackageStartupMessages({
  library(ggplot2)
})

d4_require_stage("38e_sensitivity")
stage <- "38f_release"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}
d4_msg("Starting ", stage)

support <- fread(file.path(D4_TABLE_DIR, "D4_support_metrics_by_K2.csv"))
alerts <- fread(file.path(
  D4_TABLE_DIR, "D4_structural_support_alerts_by_K2.csv"
))
gates <- fread(file.path(D4_TABLE_DIR, "D4_gate_A_to_D_by_K2.csv"))
state <- fread(file.path(D4_TABLE_DIR, "D4_claim_state.csv"))
flow <- fread(file.path(D4_TABLE_DIR, "D4_24h_landmark_flow.csv"))
events <- fread(file.path(
  D4_TABLE_DIR, "D4_24_72_first_event_by_K2.csv"
))
sensitivity <- fread(file.path(
  D4_TABLE_DIR, "D4_boundary_tie_and_complete_case_sensitivity.csv"
))

primary_groups <- support[phenotype != "Overall"]
figure_source <- rbindlist(list(
  primary_groups[, .(
    phenotype,
    metric = "Outside overlap",
    estimate = outside_overlap_proportion,
    threshold = D4_SUPPORT_THRESHOLDS$outside_overlap_max,
    alert_ratio = outside_overlap_proportion /
      D4_SUPPORT_THRESHOLDS$outside_overlap_max,
    direction = "higher_is_worse"
  )],
  primary_groups[, .(
    phenotype,
    metric = "ESS/N",
    estimate = ess_ratio,
    threshold = D4_SUPPORT_THRESHOLDS$ess_ratio_min,
    alert_ratio = D4_SUPPORT_THRESHOLDS$ess_ratio_min / ess_ratio,
    direction = "lower_is_worse"
  )],
  primary_groups[, .(
    phenotype,
    metric = "Minimum treatment cell",
    estimate = minimum_cell,
    threshold = D4_SUPPORT_THRESHOLDS$minimum_cell_min,
    alert_ratio = D4_SUPPORT_THRESHOLDS$minimum_cell_min / minimum_cell,
    direction = "lower_is_worse"
  )],
  primary_groups[, .(
    phenotype,
    metric = "Extreme weights",
    estimate = extreme_weight_proportion,
    threshold = D4_SUPPORT_THRESHOLDS$extreme_weight_prop_max,
    alert_ratio = fifelse(
      extreme_weight_proportion == 0,
      1e-3,
      extreme_weight_proportion /
        D4_SUPPORT_THRESHOLDS$extreme_weight_prop_max
    ),
    direction = "higher_is_worse"
  )]
))
figure_source[, alert := alert_ratio > 1]
figure_source[, metric := factor(
  metric,
  levels = c(
    "Outside overlap", "ESS/N", "Minimum treatment cell", "Extreme weights"
  )
)]

d4_atomic_fwrite(
  figure_source,
  file.path(D4_TABLE_DIR, "Figure5D_D4_landmark_support_source.csv")
)

palette <- c(
  "C1 higher-risk" = "#B4553F",
  "C2 lower-risk" = "#3E6E93"
)
p <- ggplot(
  figure_source,
  aes(x = alert_ratio, y = metric, color = phenotype, shape = alert)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey70", linewidth = 0.5) +
  geom_segment(
    aes(x = 1, xend = alert_ratio, yend = metric),
    position = position_dodge(width = 0.45),
    linewidth = 0.6,
    alpha = 0.65
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.4,
    stroke = 0.8
  ) +
  scale_x_log10(
    name = "Ratio to prespecified alert threshold (log scale)",
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_y_discrete(name = NULL) +
  scale_color_manual(values = palette, name = NULL) +
  scale_shape_manual(
    values = c(`FALSE` = 1, `TRUE` = 16),
    labels = c(`FALSE` = "Gate passed", `TRUE` = "Alert triggered"),
    name = NULL
  ) +
  labs(
    title = "D4 landmark treatment-support audit",
    subtitle = paste0(
      "Hour-24 risk set; state: ",
      gsub("_", " ", state$state)
    ),
    caption = paste(
      "Values >1 cross a prespecified support alert threshold.",
      "No treatment-effect model was fitted."
    )
  ) +
  theme_classic(base_family = "sans", base_size = 9.5) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    plot.caption = element_text(size = 8, hjust = 0),
    plot.margin = margin(8, 14, 8, 8)
  )

pdf_path <- file.path(
  D4_FIGURE_DIR, "Figure5D_D4_landmark_support_candidate.pdf"
)
png_path <- file.path(
  D4_FIGURE_DIR, "Figure5D_D4_landmark_support_candidate.png"
)
tiff_path <- file.path(
  D4_FIGURE_DIR, "Figure5D_D4_landmark_support_candidate.tiff"
)
if (any(file.exists(c(pdf_path, png_path, tiff_path)))) {
  stop("Refusing to overwrite an existing D4 figure.", call. = FALSE)
}
ggsave(pdf_path, p, width = 6.6, height = 4.3, units = "in",
       device = cairo_pdf)
ggsave(png_path, p, width = 6.6, height = 4.3, units = "in",
       dpi = 300, bg = "white")
ggsave(tiff_path, p, width = 6.6, height = 4.3, units = "in",
       dpi = 600, compression = "lzw", bg = "white")

result_lines <- c(
  "# D4 corrected landmark support result",
  "",
  paste0("- Run mode: `", D4_RUN_MODE, "`"),
  paste0("- State: `", state$state, "`"),
  paste0("- Trigger summary: ", state$triggers),
  paste0("- Corrected hour-24 risk set: ", tail(flow$n_remaining, 1)),
  paste0(
    "- First RRT before death/ICU exit and before hour 72: ",
    sum(events[first_event_state == "RRT_FIRST_24_72", n])
  ),
  "",
  "This result concerns observed treatment-allocation support only.",
  "No treatment effect, HTE, interaction, RD, RR, OR, or E-value was estimated.",
  "",
  "## Phenotype-specific support",
  paste(capture.output(print(primary_groups)), collapse = "\n"),
  "",
  "## Sensitivity states",
  paste(capture.output(print(
    sensitivity[, .(variant, n, treated_n, state, state_matches_primary)]
  )), collapse = "\n")
)
d4_atomic_write_lines(
  result_lines,
  file.path(D4_OUTPUT_ROOT, "D4_FORMAL_RESULT_README.md")
)

public_files <- list.files(
  D4_OUTPUT_ROOT, recursive = TRUE, full.names = TRUE
)
public_files <- public_files[
  !grepl("private_not_for_release", public_files, fixed = TRUE) &
    !file.info(public_files)$isdir
]
release_manifest <- rbindlist(lapply(public_files, function(pth) {
  data.table(
    path = sub(
      paste0("^", gsub("\\\\", "/", normalizePath(D4_OUTPUT_ROOT, winslash = "/")), "/?"),
      "",
      normalizePath(pth, winslash = "/", mustWork = TRUE)
    ),
    bytes = file.info(pth)$size,
    sha256 = d4_sha256(pth)
  )
}))
d4_atomic_fwrite(
  release_manifest,
  file.path(D4_PROVENANCE_DIR, "D4_release_manifest_preQC.csv")
)
capture.output(
  sessionInfo(),
  file = file.path(D4_PROVENANCE_DIR, "D4_sessionInfo.txt")
)
d4_mark_stage(stage)
d4_msg("Completed ", stage)

