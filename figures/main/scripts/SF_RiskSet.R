script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
source(file.path(script_dir, "common.R"))

filename <- "Figure6_panelA_landmark_flow_v2.csv"
flow_path <- file.path(DATA, "fig5", filename)
manifest <- read.csv(file.path(DATA, "fig5", "source_manifest.csv"),
                     stringsAsFactors = FALSE)
stopifnot(file.exists(flow_path), !anyDuplicated(manifest$source),
          filename %in% manifest$source,
          unname(tools::md5sum(flow_path)) == manifest$md5[match(filename, manifest$source)])
d <- read.csv(flow_path, check.names = FALSE)
d <- d[order(d$order), ]
stopifnot(identical(as.integer(d$n_remaining),
                    c(20049L, 19658L, 19262L, 17525L, 17525L)),
          identical(as.integer(d$n_removed_at_step),
                    c(0L, 391L, 396L, 1737L, 0L)),
          all(d$n_removed_at_step[d$n_removed_at_step > 0L] >= 10L),
          all(head(d$n_remaining, -1L) - tail(d$n_removed_at_step, -1L) ==
                tail(d$n_remaining, -1L)))

main_labels <- c("MIMIC-IV cohort", "No RRT-positive record\nbefore ICU hour 24",
                 "Alive at ICU hour 24", "Hour-24 risk set")
side_labels <- c("RRT-positive record\nbefore ICU hour 24",
                 "Death before hour 24", "ICU exit before hour 24")
main_text <- paste0(main_labels, "\nn = ", fmt_n(d$n_remaining[1:4]))
side_text <- paste0(side_labels, "\nExcluded: ", fmt_n(d$n_removed_at_step[2:4]))
height_mm <- 150
draw <- function() {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = FONT, fontsize = FS)))
  x_main <- 0.29
  x_side <- 0.76
  main_y <- c(0.87, 0.65, 0.43, 0.21)
  side_y <- (main_y[-length(main_y)] + main_y[-1L]) / 2
  main_w <- 0.32
  side_w <- 0.31
  box_h <- 0.125
  line_gp <- gpar(col = INK, fill = INK, lwd = max(0.8, LW * 96 / 25.4))
  arrow_tip <- arrow(type = "closed", length = unit(1.7, "mm"))
  for (i in seq_along(side_y)) {
    grid.segments(x_main, main_y[i] - box_h / 2,
                  x_main, main_y[i + 1L] + box_h / 2,
                  default.units = "npc", arrow = arrow_tip, gp = line_gp)
    grid.segments(x_main, side_y[i], x_side - side_w / 2, side_y[i],
                  default.units = "npc", arrow = arrow_tip, gp = line_gp)
  }
  for (i in seq_along(main_text)) {
    grid.rect(x_main, main_y[i], main_w, box_h, default.units = "npc",
              gp = gpar(col = INK, fill = "white", lwd = 0.8))
    grid.text(main_text[i], x_main, main_y[i], default.units = "npc",
              gp = gpar(fontfamily = FONT, fontsize = FS, col = INK))
  }
  for (i in seq_along(side_text)) {
    grid.rect(x_side, side_y[i], side_w, box_h, default.units = "npc",
              gp = gpar(col = INK, fill = "white", lwd = 0.8))
    grid.text(side_text[i], x_side, side_y[i], default.units = "npc",
              gp = gpar(fontfamily = FONT, fontsize = FS, col = INK))
  }
  popViewport()
}

export_figure(draw, "SF_RiskSet", height_mm)
stopifnot(file.copy(file.path(QA, "SF_RiskSet.png"),
                    file.path(OUT, "SF_RiskSet.png"), overwrite = TRUE))
stopifnot(unname(tools::md5sum(flow_path)) ==
            manifest$md5[match(filename, manifest$source)])
pass("figure5_riskset", c(
  "Sequential exclusions reconcile to the hour-24 risk set.",
  "Early RRT exclusions are combined in the public figure source.",
  "The restricted original flow and all locked patient-level outputs remain unchanged."
))
