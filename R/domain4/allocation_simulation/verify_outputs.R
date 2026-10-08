args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
out_dir <- file.path(dirname(normalizePath(file_arg)), "output")
read_output <- function(name) {
  utils::read.csv(file.path(out_dir, name), stringsAsFactors = FALSE)
}

settings <- read_output("settings.csv")
ledger <- read_output("replicate_ledger.csv")
failures <- read_output("failure_ledger.csv")
counts <- read_output("setting_counts.csv")
events <- read_output("event_summary.csv")
metrics <- read_output("metric_summary.csv")
checks <- character()
check <- function(name, condition) {
  if (!isTRUE(condition)) stop("Output verification failed: ", name)
  checks <<- c(checks, name)
}
near <- function(x, y, tolerance = 1e-10) {
  all(is.finite(x) & is.finite(y) & abs(x - y) <= tolerance)
}

check("31 settings and 15500 generated repeats",
      nrow(settings) == 31L && nrow(ledger) == 15500L &&
        all(table(ledger$setting_id) == 500L))
check("unique and consecutive replicate identifiers",
      !anyDuplicated(paste(ledger$setting_id, ledger$replicate)) &&
        all(vapply(split(ledger$replicate, ledger$setting_id),
                   function(x) identical(sort(x), seq_len(500L)), logical(1))))
check("allocated counts and observed prevalence agree",
      all(ledger$n_treated + ledger$n_control == ledger$N) &&
        near(ledger$observed_prevalence, ledger$n_treated / ledger$N))
check("replicate-level marginal true probabilities calibrated",
      near(ledger$true_mean, ledger$target_p, 1e-9))
structural <- ledger[ledger$scenario == "structural_zero", ]
check("structural-zero region has no treated draws",
      nrow(structural) == 500L &&
        sum(structural$n_treated_zero_region) == 0L &&
        all(structural$n_zero_region > 0L & structural$n_zero_region < structural$N) &&
        counts$zero_region_people[counts$setting_id == "S31"] ==
          sum(structural$n_zero_region))
check("working-model scores remain positive in zero region",
      all(structural$mean_fitted_e_zero_region > 0))
check("absence of generated outcomes and full-D4 verdict",
      !any(grepl("death|mortality|verdict|icu_exit|rrt_time",
                 names(ledger), ignore.case = TRUE)))

for (target in c("ATE", "ATT", "ATO")) {
  ess1 <- ledger[[paste0("ess_treated_", target)]]
  ess0 <- ledger[[paste0("ess_control_", target)]]
  retention1 <- ledger[[paste0("retention_treated_", target)]]
  retention0 <- ledger[[paste0("retention_control_", target)]]
  share1 <- ledger[[paste0("top1_share_treated_", target)]]
  share0 <- ledger[[paste0("top1_share_control_", target)]]
  pooled <- ledger[[paste0("pooled_ess_ratio_", target)]]
  valid <- is.finite(ess1) & is.finite(ess0)
  check(paste0(target, " arm ESS, retention, and top-share bounds"),
        all(ess1[valid] > 0 & ess1[valid] <= ledger$n_treated[valid] + 1e-8) &&
          all(ess0[valid] > 0 & ess0[valid] <= ledger$n_control[valid] + 1e-8) &&
          all(retention1[valid] > 0 & retention1[valid] <= 1 + 1e-10) &&
          all(retention0[valid] > 0 & retention0[valid] <= 1 + 1e-10) &&
          all(share1[valid] > 0 & share1[valid] <= 1) &&
          all(share0[valid] > 0 & share0[valid] <= 1) &&
          all(pooled[valid] > 0 & pooled[valid] <= 1 + 1e-10))
}
check("zero-treatment replicate has no two-arm support metric",
      all(is.na(ledger$pooled_ess_ratio_ATE[ledger$n_treated == 0L])) &&
        all(is.na(ledger$ess_treated_ATE[ledger$n_treated == 0L])))

issues <- ledger$fit_status != "ok" | ledger$support_status != "two_arm" |
  ledger$metric_status != "ok" | ledger$warning_count > 0L
key <- function(d) paste(d$setting_id, d$replicate)
check("failure ledger contains exactly flagged repeats",
      nrow(failures) == sum(issues) &&
        identical(sort(key(failures)), sort(key(ledger[issues, ]))))
for (i in seq_len(nrow(counts))) {
  d <- ledger[ledger$setting_id == counts$setting_id[i], ]
  check(paste0(counts$setting_id[i], " count summary agrees"),
        counts$generated_reps[i] == nrow(d) &&
          counts$two_arm_reps[i] == sum(d$support_status == "two_arm") &&
          counts$fit_ok_reps[i] == sum(d$fit_status == "ok") &&
          counts$complete_metric_reps[i] == sum(d$metric_status == "ok") &&
          counts$issue_ledger_rows[i] == sum(issues[ledger$setting_id == counts$setting_id[i]]))
}

check("event MCSEs use the generated-replicate denominator",
      all(events$generated_reps == 500L) &&
        near(events$probability, events$count / 500) &&
        near(events$mcse,
             sqrt(events$probability * (1 - events$probability) / 500)))
check("metric availability and means match replicate ledger",
      all(vapply(seq_len(nrow(metrics)), function(i) {
        d <- ledger[ledger$setting_id == metrics$setting_id[i], ]
        x <- d[[metrics$metric[i]]]
        x <- x[is.finite(x)]
        length(x) == metrics$available_reps[i] &&
          (length(x) == 0L || abs(mean(x) - metrics$mean[i]) < 1e-9) &&
          (length(x) <= 1L ||
             abs(stats::sd(x) / sqrt(length(x)) - metrics$mcse[i]) < 1e-9)
      }, logical(1))))

utils::write.csv(data.frame(check = checks, passed = TRUE),
                 file.path(out_dir, "verification_results.csv"), row.names = FALSE)
cat(length(checks), "output verification checks passed\n")
