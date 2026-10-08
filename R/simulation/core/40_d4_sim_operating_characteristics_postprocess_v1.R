#!/usr/bin/env Rscript

# Presentation-only post-processing for the frozen formal D4 simulation.
# Statistical values are read verbatim from the formal operating-characteristics CSV.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("40_d4_sim_operating_characteristics_postprocess_v1.R", winslash = "/", mustWork = TRUE)
}

root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "output_formal_v2", "tables", "D4_sim_operating_characteristics.csv")
output_dir <- file.path(root, "output_formal_v2", "figures_postprocessed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file))
oc <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
oc <- subset(oc, scenario %in% c(
  "G4L_longitudinal_support_failure",
  "G6L_longitudinal_adequate_support"
))

label_map <- c(
  G4L_longitudinal_support_failure = "G4-L support-failure control",
  G6L_longitudinal_adequate_support = "G6-L adequate-support control"
)
oc$control <- unname(label_map[oc$scenario])
oc$control <- factor(
  oc$control,
  levels = rev(unname(label_map))
)
oc$display <- sprintf(
  "%d/%d  (%.1f%%)",
  oc$correct,
  oc$n_valid,
  100 * oc$rate
)

COL_POINT <- "#3A5A8C"
COL_TEXT <- "#2E2E2E"
COL_GRID <- "grey70"

cat("D4 longitudinal OC figure palette:\n")
print(c(point = COL_POINT, text = COL_TEXT, reference = COL_GRID))
cat("Canvas: 160 x 76 mm; export: PDF/SVG + 600 dpi TIFF/PNG\n")

p <- ggplot(oc, aes(x = rate, y = control)) +
  geom_vline(xintercept = 1, linewidth = 0.45, linetype = "dashed", colour = COL_GRID) +
  geom_errorbar(
    aes(xmin = wilson_low, xmax = wilson_high),
    orientation = "y",
    width = 0,
    linewidth = 0.65,
    colour = COL_POINT
  ) +
  geom_point(size = 2.6, colour = COL_POINT) +
  geom_text(
    aes(x = 1.002, label = display),
    hjust = 0,
    size = 3.2,
    colour = COL_TEXT
  ) +
  annotate(
    "text",
    x = 0.965,
    y = 0.55,
    label = "D4-LQ: 100/100 ledger stress repeats completed\nImplementation check; excluded from performance denominators",
    hjust = 0,
    size = 2.75,
    lineheight = 1.05,
    colour = "grey35"
  ) +
  scale_x_continuous(
    limits = c(0.965, 1.022),
    breaks = c(0.97, 0.98, 0.99, 1.00),
    labels = function(x) sprintf("%.2f", x),
    expand = c(0, 0)
  ) +
  labs(
    title = "Truth-aligned Domain 4 support classification",
    subtitle = "Formal longitudinal controls; Wilson 95% confidence intervals",
    x = "Correct classification rate",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "sans", base_size = 9) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9.5, colour = COL_TEXT, margin = margin(r = 7)),
    axis.text.x = element_text(size = 8.5, colour = COL_TEXT),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 9, colour = "grey35", hjust = 0, margin = margin(b = 9)),
    plot.margin = margin(7, 44, 15, 8),
    panel.grid = element_blank()
  )

base <- file.path(output_dir, "Figure_D4L_operating_characteristics_postprocessed_v1")
ggsave(paste0(base, ".pdf"), p, width = 160, height = 76, units = "mm", device = cairo_pdf)
ggsave(paste0(base, ".svg"), p, width = 160, height = 76, units = "mm", device = svglite::svglite)
ggsave(
  paste0(base, ".tiff"), p,
  width = 160, height = 76, units = "mm", dpi = 600,
  device = ragg::agg_tiff, compression = "lzw"
)
ggsave(
  paste0(base, ".png"), p,
  width = 160, height = 76, units = "mm", dpi = 600,
  device = ragg::agg_png
)

cat("Saved publication exports to: ", output_dir, "\n", sep = "")
