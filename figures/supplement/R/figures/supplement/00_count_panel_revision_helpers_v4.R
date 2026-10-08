# Shared export helper for the count-panel revision package.
# This file changes presentation only. It does not read or modify study data.

save_count_panel_versions <- function(
    plot, stem, output_dir, single_plot = plot,
    main_width_mm = 180, main_height_mm,
    single_width_mm = 89, single_height_mm,
    dpi = 600) {
  if (main_height_mm > 230 || single_height_mm > 230) {
    stop("Requested figure height exceeds 230 mm.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  save_one <- function(plot_object, suffix, width_mm, height_mm) {
    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4
    out_stem <- file.path(output_dir, paste0(stem, "_", suffix))

    svglite::svglite(
      paste0(out_stem, ".svg"), width = width_in, height = height_in,
      bg = PAL[["paper"]]
    )
    print(plot_object)
    grDevices::dev.off()

    grDevices::cairo_pdf(
      paste0(out_stem, ".pdf"), width = width_in, height = height_in,
      family = FONT_FAMILY, bg = PAL[["paper"]]
    )
    print(plot_object)
    grDevices::dev.off()

    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(
        paste0(out_stem, ".png"), width = width_in, height = height_in,
        units = "in", res = dpi, background = PAL[["paper"]]
      )
    } else {
      grDevices::png(
        paste0(out_stem, ".png"), width = width_in, height = height_in,
        units = "in", res = dpi, bg = PAL[["paper"]], type = "cairo"
      )
    }
    print(plot_object)
    grDevices::dev.off()

    paste0(out_stem, c(".pdf", ".svg", ".png"))
  }

  message("Count-panel revision: ", stem)
  message("Main canvas: ", main_width_mm, " x ", main_height_mm, " mm")
  message("Single-column canvas: ", single_width_mm, " x ", single_height_mm, " mm")
  message("Vector PDF/SVG; PNG ", dpi, " dpi; source values unchanged")
  message("Raster device: ", if (requireNamespace("ragg", quietly = TRUE)) "ragg" else "grDevices fallback (ragg unavailable)")

  capture_supplement_render(plot, paste0(stem, "_main"), main_width_mm, main_height_mm, parent.frame())
  c(
    save_one(plot, "main", main_width_mm, main_height_mm),
    save_one(single_plot, "single", single_width_mm, single_height_mm)
  )
}
