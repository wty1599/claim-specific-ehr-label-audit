# Supplement-specific physical text floor; scientific scales are not altered.
supplement_type_scale <- function(main_sizes) {
  c(panel_letter = main_sizes[["panel_letter"]],
    panel_title = main_sizes[["panel_title"]], subtitle = 8,
    axis_title = 9, axis_text = 8, annotation = 8, legend = 8, small = 8)
}

capture_supplement_render <- function(plot, stem, width_mm, height_mm, env) {
  dest <- file.path(ROOT_FROM_THEME, "qa", "snapshots")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  names <- ls(env, all.names = TRUE)
  objects <- mget(names, envir = env, inherits = FALSE)
  tables <- Filter(is.data.frame, objects)
  panels <- Filter(function(x) inherits(x, "ggplot") && !inherits(x, "patchwork"), objects)
  built <- lapply(panels, function(p) {
    b <- ggplot2::ggplot_build(p)
    list(data = b$data, layout = b$layout$layout,
         geoms = vapply(p$layers, function(l) class(l$geom)[1], character(1)),
         ranges = lapply(b$layout$panel_params, function(x) list(x=x$x.range,y=x$y.range)))
  })
  saveRDS(list(stem=stem, width_mm=width_mm, height_mm=height_mm,
               tables=tables, panels=built), file.path(dest,paste0(stem,".rds")))
  saveRDS(plot, file.path(dest,paste0(stem,"_plot.rds")))
}
