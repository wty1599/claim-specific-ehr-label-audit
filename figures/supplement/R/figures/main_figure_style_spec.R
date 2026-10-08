# Read-only transcription of the shared main-figure visual specification.
read_main_figure_style_spec <- function(root) {
  spec <- jsonlite::fromJSON(file.path(root, "figure_spec.json"), simplifyVector = TRUE)
  stopifnot(spec$font$family == "Arial", spec$export$dpi == 600)
  spec
}
