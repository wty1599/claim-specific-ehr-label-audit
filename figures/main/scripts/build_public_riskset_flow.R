args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript build_public_riskset_flow.R private_flow.csv public_flow.csv")
source <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
target <- args[2]
if (normalizePath(dirname(target), winslash = "/", mustWork = TRUE) == dirname(source)) {
  stop("The publication copy must not overwrite the restricted source.")
}

d <- read.csv(source, check.names = FALSE, stringsAsFactors = FALSE)
expected <- c("Locked cohort", "Exclude pre-ICU RRT",
              "Exclude RRT from ICU admission to <24h",
              "Exclude death before 24h", "Exclude ICU exit before 24h",
              "Final hour-24 opportunity risk set")
stopifnot(identical(d$step, expected),
          all(d$n_remaining[-1L] == d$n_remaining[-nrow(d)] -
                d$n_removed_at_step[-1L]),
          d$n_remaining[1L] == 20049L,
          tail(d$n_remaining, 1L) == 17525L)

public <- data.frame(
  step = c("Locked cohort", "Exclude RRT before ICU hour 24",
           "Exclude death before 24h", "Exclude ICU exit before 24h",
           "Final hour-24 opportunity risk set"),
  n_remaining = d$n_remaining[c(1L, 3L, 4L, 5L, 6L)],
  n_removed_at_step = c(0L, sum(d$n_removed_at_step[2:3]),
                        d$n_removed_at_step[4:5], 0L),
  step_short = c("Locked\ncohort", "No RRT\nbefore 24 h", "Alive at\n24 h",
                 "In ICU at\n24 h", "Hour-24\nrisk set"),
  order = seq_len(5L),
  stringsAsFactors = FALSE
)
public$removed_label <- ifelse(public$n_removed_at_step == 0L, "",
                               paste0("-", format(public$n_removed_at_step,
                                                   big.mark = ",", trim = TRUE)))
stopifnot(all(public$n_removed_at_step[public$n_removed_at_step > 0L] >= 10L),
          all(public$n_remaining[-1L] == public$n_remaining[-nrow(public)] -
                public$n_removed_at_step[-1L]))
write.csv(public, target, row.names = FALSE, na = "")
