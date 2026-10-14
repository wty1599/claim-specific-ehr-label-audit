args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript prepare_Figure3_shape_summary.R repeat_csv output_csv")
}
repeat_rows <- read.csv(args[1L], check.names = FALSE)
required <- c("condition", "n_eval", "repeat_id", "fisher_axis_dip_p_bh6",
              "shape_alert_threshold", "shape_alert")
stopifnot(all(required %in% names(repeat_rows)), nrow(repeat_rows) == 900L,
          !anyDuplicated(repeat_rows[c("condition", "n_eval", "repeat_id")]),
          all(is.finite(repeat_rows$fisher_axis_dip_p_bh6)),
          all(repeat_rows$shape_alert_threshold == 0.05),
          all(repeat_rows$shape_alert ==
                as.integer(repeat_rows$fisher_axis_dip_p_bh6 < 0.05)))
groups <- split(repeat_rows,
                interaction(repeat_rows$condition, repeat_rows$n_eval, drop = TRUE))
summary_rows <- do.call(rbind, lapply(groups, function(z) {
  data.frame(condition = z$condition[1L], n_eval = z$n_eval[1L],
             count = sum(z$shape_alert), n = nrow(z),
             median_p = median(z$fisher_axis_dip_p_bh6),
             shape_alert_threshold = 0.05)
}))
stopifnot(nrow(summary_rows) == 9L, all(summary_rows$n == 100L),
          identical(sort(unique(summary_rows$n_eval)), c(1000L, 5000L, 10025L)))
summary_rows <- summary_rows[order(summary_rows$condition,
                                   summary_rows$n_eval), ]
write.csv(summary_rows, args[2L], row.names = FALSE)
