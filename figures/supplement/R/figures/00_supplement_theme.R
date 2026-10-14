# Shared presentation specification for Supplementary Figures S1-S12.
#
# This file reads the same figure_spec.json used by the current main figures.
# Supplementary panels apply an 8-point floor to dense text, as required by the
# submission specification. No study data or statistical decision logic lives
# in this file.

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(grid)
  library(patchwork)
})

ROOT_FROM_THEME <- normalizePath(
  file.path(dirname(sys.frame(1)$ofile), "..", ".."),
  winslash = "/", mustWork = TRUE
)

source(file.path(ROOT_FROM_THEME, "R", "figures", "main_figure_style_spec.R"))
source(file.path(ROOT_FROM_THEME, "R", "figures", "supplementary_style_adapter.R"))
FIGURE_SPEC <- read_main_figure_style_spec(ROOT_FROM_THEME)

PAL <- unlist(FIGURE_SPEC$palette, use.names = TRUE)
MAIN_FS <- unlist(FIGURE_SPEC$sizes_pt, use.names = TRUE)
MAIN_LW <- unlist(FIGURE_SPEC$line_width_pt, use.names = TRUE)
FONT_FAMILY <- FIGURE_SPEC$font$family
options(device = function(...) grDevices::cairo_pdf(file = tempfile(fileext = ".pdf"),
                                                   family = FONT_FAMILY, ...))
PT_TO_MM <- 25.4 / 72

# The main figures and supplementary figures share one palette and type family.
# The supplementary floor is raised to 8 pt without changing relative roles.
SUPP_FS <- supplement_type_scale(MAIN_FS)

SUPP_COLORS <- c(
  structure = PAL[["raw_en"]],
  alert = FIGURE_SPEC$states$S4$stroke,
  alert_pale = FIGURE_SPEC$states$S4$fill,
  neutral_focus = PAL[["ink"]],
  neutral_main = PAL[["muted"]],
  neutral_dark = PAL[["neutral_dark"]],
  neutral_mid = PAL[["muted"]],
  neutral_light = PAL[["reference"]],
  neutral_line = PAL[["muted"]],
  reference = PAL[["reference"]],
  grid = PAL[["grid"]],
  border = PAL[["border"]],
  C1 = PAL[["c1"]],
  C2 = PAL[["c2"]],
  renal = PAL[["renal"]],
  acidbase = PAL[["acidbase"]],
  other = PAL[["other"]],
  baseline = PAL[["base"]],
  base_k2 = PAL[["base_k2"]],
  raw_en = PAL[["raw_en"]],
  raw_rf = PAL[["raw_rf"]],
  indigo = PAL[["raw_en"]],
  K2 = PAL[["ink"]],
  K3 = PAL[["muted"]],
  K4 = PAL[["reference"]],
  pale_grey = PAL[["neutral_fill"]],
  paper = PAL[["paper"]]
)

# Backward-compatible semantic aliases used by the locked presentation scripts.
# They point to the current main-figure specification rather than creating a
# second palette.
PAL_NEUTRAL <- c(
  ink = PAL[["ink"]], text_muted = PAL[["muted"]],
  outline = PAL[["muted"]], reference = PAL[["reference"]],
  grid = PAL[["grid"]], band = PAL[["neutral_fill"]], paper = PAL[["paper"]]
)
PAL_ENTITY <- c(C1 = PAL[["c1"]], C2 = PAL[["c2"]])
PAL_BLOCK <- c(renal = PAL[["renal"]], acidbase = PAL[["acidbase"]], other = PAL[["other"]])
PAL_MODEL <- c(
  baseline = PAL[["base"]], base_k2 = PAL[["base_k2"]],
  raw_en = PAL[["raw_en"]], raw_rf = PAL[["raw_rf"]]
)

activate_supplement_style <- function(
    palette_object = "COLORS",
    font_object = "FONT_PT",
    family_object = "BASE_FAMILY") {
  env <- parent.frame()
  assign(palette_object, SUPP_COLORS, envir = env)
  assign(font_object, c(
    panel_title = SUPP_FS[["panel_title"]],
    panel_tag = SUPP_FS[["panel_letter"]],
    subtitle = SUPP_FS[["subtitle"]],
    axis_title = SUPP_FS[["axis_title"]],
    axis_text = SUPP_FS[["axis_text"]],
    annotation = SUPP_FS[["annotation"]],
    legend = SUPP_FS[["legend"]],
    tag = SUPP_FS[["panel_letter"]],
    panel = SUPP_FS[["panel_title"]]
  ), envir = env)
  assign(family_object, FONT_FAMILY, envir = env)
  invisible(TRUE)
}

activate_supplement_palette <- function(palette_object = "PALETTE") {
  env <- parent.frame()
  assign(palette_object, c(
    neutral_dark = SUPP_COLORS[["neutral_focus"]],
    neutral_mid = SUPP_COLORS[["neutral_main"]],
    neutral_light = SUPP_COLORS[["neutral_light"]],
    grid = SUPP_COLORS[["grid"]],
    reference = SUPP_COLORS[["reference"]],
    indigo = SUPP_COLORS[["indigo"]],
    alert_amber = SUPP_COLORS[["alert"]],
    shade = SUPP_COLORS[["pale_grey"]],
    paper = SUPP_COLORS[["paper"]]
  ), envir = env)
  invisible(TRUE)
}

theme_supplement <- function() {
  theme_classic(base_family = FONT_FAMILY, base_size = SUPP_FS[["axis_text"]]) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = PAL[["ink"]]),
      plot.title = element_text(
        size = SUPP_FS[["panel_title"]], face = "bold", hjust = 0,
        margin = margin(b = 2.5)
      ),
      plot.subtitle = element_text(
        size = SUPP_FS[["subtitle"]], colour = PAL[["muted"]], hjust = 0,
        margin = margin(b = 2.5)
      ),
      plot.caption = element_text(
        size = SUPP_FS[["annotation"]], colour = PAL[["muted"]], hjust = 0
      ),
      plot.tag = element_text(
        size = SUPP_FS[["panel_letter"]], face = "bold", hjust = 0
      ),
      plot.tag.position = c(0, 1),
      axis.title = element_text(size = SUPP_FS[["axis_title"]]),
      axis.text = element_text(size = SUPP_FS[["axis_text"]], colour = PAL[["ink"]]),
      axis.line = element_line(
        colour = PAL[["muted"]], linewidth = MAIN_LW[["axis"]] * PT_TO_MM
      ),
      axis.ticks = element_line(
        colour = PAL[["muted"]], linewidth = MAIN_LW[["axis"]] * PT_TO_MM
      ),
      axis.ticks.length = unit(1.25, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        size = SUPP_FS[["axis_text"]], face = "bold", hjust = 0,
        margin = margin(t = 2.2, b = 2.2)
      ),
      legend.position = "bottom",
      legend.title = element_text(size = SUPP_FS[["legend"]], face = "bold"),
      legend.text = element_text(size = SUPP_FS[["legend"]]),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(4.5, "mm"),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(4, 4.5, 3.5, 4.5, unit = "mm")
    )
}

theme_paper <- theme_supplement

tag_panel <- function(plot, letter, title = NULL) {
  if (!is.null(title)) plot <- plot + labs(title = title)
  plot + labs(tag = letter) +
    theme(
      plot.tag = element_text(
        size = SUPP_FS[["panel_letter"]], face = "bold", hjust = 0
      ),
      plot.tag.position = c(0, 1),
      plot.title = element_text(
        size = SUPP_FS[["panel_title"]], face = "bold", hjust = 0,
        margin = margin(b = 2.5, l = 17)
      )
    )
}

assert_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing figure input(s): ", paste(missing, collapse = "; "))
  invisible(TRUE)
}

supplement_header <- function(label) {
  ggplot() +
    annotate(
      "segment", x = 0, xend = 1, y = .12, yend = .12,
      colour = PAL[["border"]], linewidth = .42
    ) +
    annotate(
      "text", x = .01, y = .64, label = label, hjust = 0, vjust = .5,
      family = FONT_FAMILY, fontface = "bold",
      size = SUPP_FS[["axis_text"]] / ggplot2::.pt,
      colour = PAL[["ink"]]
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FONT_FAMILY)
}

add_supplement_header <- function(plot, label, fraction = 0.052) {
  plot & theme(panel.background = element_rect(fill = PAL[["paper"]], colour = NA),
               plot.background = element_rect(fill = PAL[["paper"]], colour = NA),
               panel.grid = element_blank(), strip.background = element_blank())
}

save_supplement_figure <- function(
    plot, stem, output_dir, width_mm = 180, height_mm, dpi = 600) {
  if (height_mm > 230) stop("Supplementary figure height exceeds 230 mm: ", height_mm)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  message("Figure specification: ", FIGURE_SPEC$version)
  message("Canvas: ", width_mm, " x ", height_mm, " mm; Arial; 600 dpi")
  message("Palette: ", paste(names(PAL), PAL, collapse = "; "))
  message("Supplementary font floor: 8 pt")
  has_ragg <- requireNamespace("ragg", quietly = TRUE)
  png_device <- if (has_ragg) ragg::agg_png else "png"
  tiff_device <- if (has_ragg) ragg::agg_tiff else "tiff"
  message("Raster device: ", if (has_ragg) "ragg" else "grDevices fallback (ragg unavailable)")
  ggsave(
    file.path(output_dir, paste0(stem, ".svg")), plot,
    width = width_in, height = height_in, units = "in",
    device = svglite::svglite, bg = PAL[["paper"]]
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width_in, height = height_in, units = "in",
    device = grDevices::cairo_pdf, family = FONT_FAMILY,
    bg = PAL[["paper"]]
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width_in, height = height_in, units = "in",
    device = png_device, dpi = dpi, bg = PAL[["paper"]]
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".tiff")), plot,
    width = width_in, height = height_in, units = "in",
    device = tiff_device, dpi = dpi, compression = "lzw",
    bg = PAL[["paper"]]
  )
  capture_supplement_render(plot, stem, width_mm, height_mm, parent.frame())
  invisible(file.path(output_dir, stem))
}

reader_scenario_labels <- c(
  S1_pure_severity_continuum = "Continuous generator",
  S7_discrete_structure_outcome_null = "Outcome-null discrete generator",
  S5_domain2_positive_control = "Noisy-oracle generator",
  S5_domain2_incremental_value_positive_control = "Noisy-oracle generator",
  S2_true_discrete_subtypes = "Three-component discrete control",
  S3_external_shift = "Transport-drift control",
  S3_external_prevalence_shift = "Transport-drift control",
  S4_positivity_failure = "Support-failure control",
  S6_domain4_actionability_positive_control = "Adequate-support control"
)

print_supplement_spec <- function(stem, height_mm) {
  message("Figure: ", stem)
  message("Shared main-figure spec: ", FIGURE_SPEC$version)
  message("Canvas: 180 x ", height_mm, " mm")
  message("Font: ", FONT_FAMILY, "; minimum displayed text: 8 pt")
}
