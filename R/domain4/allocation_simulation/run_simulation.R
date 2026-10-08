args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(file_arg))
source(file.path(script_dir, "diagnostics.R"))
mode <- commandArgs(trailingOnly = TRUE)
if (length(mode) != 1L || !mode %in% c("--smoke", "--full")) {
  stop("Use --smoke, inspect its checks, then use --full")
}

out_dir <- file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE)
seed <- 20260924L
reps <- 500L
targets <- c("ATE", "ATT", "ATO")
target_metric_stems <- c("pooled_ess_ratio", "ess_treated", "ess_control",
                         "retention_treated", "retention_control",
                         "top1_share_treated", "top1_share_control", "abs_smd")
target_metrics <- as.vector(outer(target_metric_stems, targets, paste, sep = "_"))
metric_names <- c("observed_prevalence", "unweighted_abs_smd",
                  "overlap_below_treated", "overlap_above_treated",
                  "overlap_below_control", "overlap_above_control",
                  target_metrics, "extreme_ate_weight_fraction", "max_ate_weight",
                  "clipped_fraction", "zero_region_fraction",
                  "mean_fitted_e_zero_region")

settings <- expand.grid(N = c(3151L, 14374L),
                        target_p = c(0.003, 0.01, 0.05, 0.10, 0.50),
                        beta = c(0, 0.6, 1.5),
                        KEEP.OUT.ATTRS = FALSE)
settings$setting_id <- sprintf("S%02d", seq_len(nrow(settings)))
settings$scenario <- "logistic"
settings <- settings[, c("setting_id", "scenario", "N", "target_p", "beta")]
settings <- rbind(settings,
                  data.frame(setting_id = "S31", scenario = "structural_zero",
                             N = 14374L, target_p = 0.05, beta = NA_real_))
stopifnot(nrow(settings) == 31L, !anyDuplicated(settings$setting_id))
utils::write.csv(settings, file.path(out_dir, "settings.csv"), row.names = FALSE,
                 na = "")

new_row <- function(setting, replicate) {
  c(list(setting_id = setting$setting_id, scenario = setting$scenario,
         N = setting$N, target_p = setting$target_p, beta = setting$beta,
         replicate = replicate, n_treated = NA_integer_, n_control = NA_integer_,
         n_zero_region = NA_integer_, n_treated_zero_region = NA_integer_,
         true_mean = NA_real_, n_clipped_low = NA_integer_,
         n_clipped_high = NA_integer_, empty_overlap = NA,
         fit_status = "not_run", support_status = "not_drawn",
         metric_status = "not_run", warning_count = 0L, warning_text = ""),
    as.list(stats::setNames(rep(NA_real_, length(metric_names)), metric_names)))
}

simulate_one <- function(setting, replicate) {
  row <- new_row(setting, replicate)
  n <- setting$N
  x <- stats::rnorm(n)
  is_structural <- setting$scenario == "structural_zero"
  zero_region <- if (is_structural) x > 0 else rep(FALSE, n)
  if (is_structural) {
    row$n_zero_region <- sum(zero_region)
    row$zero_region_fraction <- mean(zero_region)
  }

  calibration <- tryCatch({
    if (is_structural) {
      eligible <- !zero_region
      conditional_p <- setting$target_p / mean(eligible)
      if (conditional_p <= 0 || conditional_p >= 1) {
        stop("Structural-zero target cannot be calibrated")
      }
      e_true <- ifelse(eligible, conditional_p, 0)
    } else {
      alpha <- if (setting$beta == 0) {
        stats::qlogis(setting$target_p)
      } else {
        stats::uniroot(
          function(a) mean(stats::plogis(a + setting$beta * x)) - setting$target_p,
          interval = c(-40, 40), tol = 1e-12)$root
      }
      e_true <- stats::plogis(alpha + setting$beta * x)
    }
    e_true
  }, error = function(err) err)
  if (inherits(calibration, "error")) {
    row$metric_status <- "calibration_error"
    row$warning_text <- conditionMessage(calibration)
    return(row)
  }

  e_true <- calibration
  row$true_mean <- mean(e_true)
  A <- stats::rbinom(n, size = 1L, prob = e_true)
  row$n_treated <- sum(A)
  row$n_control <- n - row$n_treated
  if (is_structural) row$n_treated_zero_region <- sum(A[zero_region])
  row$observed_prevalence <- row$n_treated / n
  row$support_status <- if (row$n_treated == 0L) "no_treated" else
    if (row$n_control == 0L) "no_control" else "two_arm"

  fit_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(stats::glm(A ~ x, family = stats::binomial()),
                        warning = function(w) {
                          fit_warnings <<- c(fit_warnings, conditionMessage(w))
                          invokeRestart("muffleWarning")
                        }),
    error = function(err) err)
  row$warning_count <- length(fit_warnings)
  row$warning_text <- paste(unique(fit_warnings), collapse = " | ")
  if (inherits(fit, "error")) {
    row$fit_status <- "error"
    row$metric_status <- "unavailable_fit"
    row$warning_text <- paste(c(row$warning_text, conditionMessage(fit)),
                              collapse = " | ")
    return(row)
  }
  if (!isTRUE(fit$converged)) {
    row$fit_status <- "nonconverged"
    row$metric_status <- "unavailable_fit"
    return(row)
  }
  raw_e <- fit$fitted.values
  if (length(raw_e) != n || any(!is.finite(raw_e)) ||
      any(raw_e < 0 | raw_e > 1) || any(!is.finite(stats::coef(fit)))) {
    row$fit_status <- "invalid_scores"
    row$metric_status <- "unavailable_fit"
    return(row)
  }
  row$fit_status <- "ok"
  row$n_clipped_low <- sum(raw_e < 1e-6)
  row$n_clipped_high <- sum(raw_e > 1 - 1e-6)
  row$clipped_fraction <- (row$n_clipped_low + row$n_clipped_high) / n
  e <- pmin(pmax(raw_e, 1e-6), 1 - 1e-6)
  if (is_structural) row$mean_fitted_e_zero_region <- mean(e[zero_region])
  if (row$support_status != "two_arm") {
    row$metric_status <- paste0("unavailable_", row$support_status)
    return(row)
  }

  overlap <- empirical_overlap(A, e)
  row$empty_overlap <- as.logical(overlap["empty"])
  row$overlap_below_treated <- overlap["below_treated"]
  row$overlap_above_treated <- overlap["above_treated"]
  row$overlap_below_control <- overlap["below_control"]
  row$overlap_above_control <- overlap["above_control"]
  row$unweighted_abs_smd <- absolute_smd(x, A)
  for (target in targets) {
    w <- target_weights(A, e, target)
    arm1 <- arm_diagnostics(w[A == 1L])
    arm0 <- arm_diagnostics(w[A == 0L])
    row[[paste0("pooled_ess_ratio_", target)]] <- pooled_ess(w) / n
    row[[paste0("ess_treated_", target)]] <- arm1["ess"]
    row[[paste0("ess_control_", target)]] <- arm0["ess"]
    row[[paste0("retention_treated_", target)]] <- arm1["retention"]
    row[[paste0("retention_control_", target)]] <- arm0["retention"]
    row[[paste0("top1_share_treated_", target)]] <- arm1["top1_share"]
    row[[paste0("top1_share_control_", target)]] <- arm0["top1_share"]
    row[[paste0("abs_smd_", target)]] <- absolute_smd(x, A, w)
    if (target == "ATE") {
      row$extreme_ate_weight_fraction <- mean(w > 10)
      row$max_ate_weight <- max(w)
    }
  }
  required_metrics <- c("unweighted_abs_smd", "overlap_below_treated",
                        "overlap_above_treated", "overlap_below_control",
                        "overlap_above_control", target_metrics,
                        "extreme_ate_weight_fraction", "max_ate_weight")
  row$metric_status <- if (all(is.finite(unlist(row[required_metrics])))) {
    "ok"
  } else {
    "partial_metric"
  }
  row
}

run_records <- function(selected_settings, n_reps) {
  chunks <- vector("list", nrow(selected_settings))
  for (i in seq_len(nrow(selected_settings))) {
    setting <- selected_settings[i, ]
    rows <- lapply(seq_len(n_reps), function(r) {
      as.data.frame(simulate_one(setting, r), stringsAsFactors = FALSE)
    })
    chunks[[i]] <- do.call(rbind, rows)
    message(setting$setting_id, ": ", n_reps, " replicates complete")
  }
  rownames_result <- do.call(rbind, chunks)
  rownames(rownames_result) <- NULL
  rownames_result
}

issue_rows <- function(d) {
  d$fit_status != "ok" | d$support_status != "two_arm" |
    d$metric_status != "ok" | d$warning_count > 0L
}

set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
if (mode == "--smoke") {
  smoke_index <- c(which(settings$N == 3151L & settings$target_p == 0.003 &
                           settings$beta == 1.5),
                   which(settings$N == 3151L & settings$target_p == 0.50 &
                           settings$beta == 0), 31L)
  stopifnot(length(smoke_index) == 3L)
  smoke <- run_records(settings[smoke_index, ], 2L)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  replay <- simulate_one(settings[smoke_index[1L], ], 1L)
  replay_frame <- as.data.frame(replay, stringsAsFactors = FALSE)
  rownames(replay_frame) <- NULL
  checks <- data.frame(
    check = c("six_generated", "counts_conserved", "calibrated_true_mean",
              "structural_zero_respected", "seed_replay_identical",
              "failure_ledger_subsets_smoke"),
    passed = c(nrow(smoke) == 6L,
               all(smoke$n_treated + smoke$n_control == smoke$N),
               all(abs(smoke$true_mean - smoke$target_p) < 1e-9),
               all(smoke$n_treated_zero_region[smoke$scenario == "structural_zero"] == 0L),
               isTRUE(all.equal(smoke[1L, , drop = FALSE], replay_frame,
                                tolerance = 0, check.attributes = FALSE)),
               sum(issue_rows(smoke)) <= nrow(smoke)))
  utils::write.csv(smoke, file.path(out_dir, "smoke_ledger.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(smoke[issue_rows(smoke), , drop = FALSE],
                   file.path(out_dir, "smoke_failure_ledger.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(checks, file.path(out_dir, "smoke_checks.csv"),
                   row.names = FALSE)
  if (!all(checks$passed)) stop("Smoke checks failed; inspect smoke_checks.csv")
  cat("Smoke checks passed; smoke draws are excluded from the formal run.\n")
  quit(save = "no")
}

smoke_checks_path <- file.path(out_dir, "smoke_checks.csv")
if (!file.exists(smoke_checks_path) ||
    !all(utils::read.csv(smoke_checks_path)$passed)) {
  stop("The smoke run must pass before the formal run")
}
formal <- run_records(settings, reps)
stopifnot(nrow(formal) == 31L * reps,
          all(table(formal$setting_id) == reps),
          all(!is.na(formal$n_treated) | formal$metric_status == "calibration_error"),
          all(is.na(formal$n_treated_zero_region[formal$scenario != "structural_zero"])),
          all(formal$n_treated_zero_region[formal$scenario == "structural_zero"] == 0L))
utils::write.csv(formal, file.path(out_dir, "replicate_ledger.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(formal[issue_rows(formal), , drop = FALSE],
                 file.path(out_dir, "failure_ledger.csv"),
                 row.names = FALSE, na = "")

summarize_metric <- function(d, name) {
  x <- d[[name]]
  x <- x[is.finite(x)]
  available <- length(x)
  if (available == 0L) {
    return(data.frame(metric = name, generated_reps = nrow(d),
                      available_reps = 0L, mean = NA_real_, mcse = NA_real_,
                      sd = NA_real_, q05 = NA_real_, median = NA_real_,
                      q95 = NA_real_))
  }
  spread <- if (available > 1L) stats::sd(x) else NA_real_
  quantiles <- stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  data.frame(metric = name, generated_reps = nrow(d),
             available_reps = available, mean = mean(x),
             mcse = spread / sqrt(available), sd = spread,
             q05 = quantiles[1L], median = quantiles[2L], q95 = quantiles[3L])
}

event_names <- c("zero_treated", "zero_control", "calibration_failure",
                 "fit_error", "fit_nonconverged", "fit_invalid_scores",
                 "fit_warning", "incomplete_support_metrics")
event_values <- function(d) {
  list(zero_treated = !is.na(d$n_treated) & d$n_treated == 0L,
       zero_control = !is.na(d$n_control) & d$n_control == 0L,
       calibration_failure = d$metric_status == "calibration_error",
       fit_error = d$fit_status == "error",
       fit_nonconverged = d$fit_status == "nonconverged",
       fit_invalid_scores = d$fit_status == "invalid_scores",
       fit_warning = d$warning_count > 0L,
       incomplete_support_metrics = d$metric_status != "ok")
}

metric_chunks <- vector("list", nrow(settings))
event_chunks <- vector("list", nrow(settings))
count_chunks <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  d <- formal[formal$setting_id == settings$setting_id[i], , drop = FALSE]
  metric_chunks[[i]] <- do.call(rbind, lapply(metric_names, function(name) {
    cbind(setting_id = settings$setting_id[i], summarize_metric(d, name))
  }))
  events <- event_values(d)
  event_chunks[[i]] <- do.call(rbind, lapply(event_names, function(name) {
    count <- sum(events[[name]])
    p_hat <- count / reps
    data.frame(setting_id = settings$setting_id[i], event = name,
               count = count, generated_reps = reps, probability = p_hat,
               mcse = sqrt(p_hat * (1 - p_hat) / reps))
  }))
  count_chunks[[i]] <- data.frame(
    setting_id = settings$setting_id[i], generated_reps = nrow(d),
    two_arm_reps = sum(d$support_status == "two_arm"),
    fit_ok_reps = sum(d$fit_status == "ok"),
    complete_metric_reps = sum(d$metric_status == "ok"),
    zero_treated_reps = sum(events$zero_treated),
    zero_control_reps = sum(events$zero_control),
    fit_warning_reps = sum(events$fit_warning),
    issue_ledger_rows = sum(issue_rows(d)),
    zero_region_people = if (settings$scenario[i] == "structural_zero")
      sum(d$n_zero_region) else NA_integer_,
    treated_in_zero_region = if (settings$scenario[i] == "structural_zero")
      sum(d$n_treated_zero_region) else NA_integer_)
}
utils::write.csv(do.call(rbind, metric_chunks),
                 file.path(out_dir, "metric_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(do.call(rbind, event_chunks),
                 file.path(out_dir, "event_summary.csv"), row.names = FALSE, na = "")
utils::write.csv(do.call(rbind, count_chunks),
                 file.path(out_dir, "setting_counts.csv"), row.names = FALSE, na = "")
writeLines(c(paste0("Seed: ", seed), paste0("R: ", R.version.string),
             paste0("Settings: ", nrow(settings)),
             paste0("Generated replicates: ", nrow(formal)),
             "Smoke results excluded from formal summaries.",
             "Structural-zero scenario: X>0 has true probability zero; X<=0 has constant calibrated probability.",
             "No deaths, ICU exits, RRT times, treatment effects, or full-D4 verdict were generated.",
             capture.output(utils::sessionInfo())),
           file.path(out_dir, "run_metadata.txt"))
cat("Formal run complete: ", nrow(formal), " generated replicates; ",
    sum(formal$metric_status == "ok"), " complete metric records.\n", sep = "")
