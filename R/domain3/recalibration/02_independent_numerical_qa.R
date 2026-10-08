options(stringsAsFactors = FALSE)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
stopifnot(nzchar(output_root))
root <- file.path(output_root, "domain3_recalibration")
private <- file.path(root, "private_local", "joined_holdout_records.csv")
point_path <- file.path(root, "aggregate_outputs",
                        "heldout_performance_point.csv")
split_path <- file.path(root, "aggregate_outputs",
                        "hospital_holdout_split_summary.csv")
interval_path <- file.path(root, "aggregate_outputs",
                           "heldout_performance_intervals.csv")
stopifnot(all(file.exists(c(private, point_path, split_path, interval_path))))
d <- read.csv(private, colClasses = c(patient = "character",
                                     hospital = "character"))
point <- read.csv(point_path)
split <- read.csv(split_path)
interval <- read.csv(interval_path)
u <- d[d$side == "updating", ]
v <- d[d$side == "evaluation", ]
stopifnot(nrow(u) + nrow(v) == 17284L,
          length(intersect(unique(u$hospital), unique(v$hospital))) == 0L,
          length(intersect(unique(u$patient), unique(v$patient))) == 0L,
          nrow(u) == split$stay_n[split$side == "updating"],
          nrow(v) == split$stay_n[split$side == "evaluation"],
          sum(u$y) == split$event_n[split$side == "updating"],
          sum(v$y) == split$event_n[split$side == "evaluation"])
for (model in c("base", "base_l1_k2")) {
  pcol <- if (model == "base") "p_base" else "p_k2"
  lp_up <- qlogis(u[[pcol]])
  lp_ev <- qlogis(v[[pcol]])
  intercept <- coef(glm(u$y ~ 1, offset = lp_up, family = binomial()))[1]
  joint <- coef(glm(u$y ~ lp_up, family = binomial()))
  predictions <- list(
    no_update = v[[pcol]],
    intercept_only = plogis(intercept + lp_ev),
    intercept_slope = plogis(joint[1] + joint[2] * lp_ev))
  for (method in names(predictions)) {
    p <- predictions[[method]]
    z <- point[point$model == model & point$method == method, ]
    stopifnot(nrow(z) == 1L,
              abs(mean((v$y - p)^2) - z$brier) < 1e-10,
              abs(mean(p) - z$mean_prediction) < 1e-10,
              abs(coef(glm(v$y ~ 1, offset = qlogis(p),
                           family = binomial()))[1] - z$citl) < 1e-9,
              abs(coef(glm(v$y ~ qlogis(p),
                           family = binomial()))[2] - z$slope) < 1e-9)
  }
  base_auc <- point$auc[point$model == model &
                          point$method == "no_update"]
  stopifnot(all(abs(point$auc[point$model == model] - base_auc) < 1e-12))
}
primary <- interval[interval$metric == "delta_brier_vs_no_update" &
                      interval$method == "intercept_slope", ]
stopifnot(nrow(primary) == 2L,
          all(primary$valid_replicates == 1000L),
          all(primary$ci_low < 0 & primary$ci_high < 0))
cat("INDEPENDENT_QA_PASS: matched patients/hospitals, split counts,",
    "both-model refits, Brier, CITL, slopes, invariant AUC,",
    "and 1000 paired primary intervals\n")
