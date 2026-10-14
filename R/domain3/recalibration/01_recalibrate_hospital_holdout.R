options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(splines)
})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
script <- normalizePath(sub("^--file=", "", arg), winslash = "/")
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
old_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT")
plan <- Sys.getenv("D3_RECAL_ORIGINAL_PLAN")
stopifnot(nzchar(output_root), nzchar(old_root), nzchar(repo), nzchar(plan))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(output_root, repo)
root <- file.path(output_root, "domain3_recalibration")
private <- file.path(root, "private_local")
aggregate <- file.path(root, "aggregate_outputs")
qa <- file.path(root, "qa")
dir.create(private, recursive = TRUE, showWarnings = FALSE)
dir.create(aggregate, recursive = TRUE, showWarnings = FALSE)
dir.create(qa, recursive = TRUE, showWarnings = FALSE)

sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                                    serialize = FALSE))
stopifnot(sha(plan) ==
  "EC95EC5CBAF1D91F4242D79C443FAC1BDCEDA9F5C0FF27F1702D5DEC3D503AB0")
source_dir <- file.path(output_root, "domain3_L1")
new_path <- file.path(source_dir, "private_local",
                      "L1_external_predictions_local.csv")
model_path <- file.path(source_dir, "private_local",
                        "L1_fitted_glm_models_local.rds")
old_path <- file.path(old_root, "domain3_deployable_generator_v1_20260728",
                      "private_local/d3_generator_v1_external_predictions_local.csv")
old_perf_path <- file.path(old_root,
  "domain3_deployable_generator_v1_20260728/tables/d3_generator_v1_model_performance.csv")
new_perf_path <- file.path(source_dir, "aggregate_outputs",
                           "L1_base_k2_external_performance.csv")
expected <- c(
  new = "A6BC6D188A619677A143493AD2DB557DD00B0984E0BCFF79B7C7AAA06CF4D364",
  old = "DE04917CF28434AF3A356B464736A7B79E15E8DF0BD0D1872EF74E800CECA019"
)
stopifnot(file.exists(model_path))
actual <- c(new = sha(new_path), old = sha(old_path))
stopifnot(identical(actual, expected))
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))

old <- fread(old_path, showProgress = FALSE)
new <- fread(new_path, showProgress = FALSE)
old <- old[outcome == "mortality" & model == "base",
           .(patientunitstayid, patient_old = as.character(uniquepid),
             hospital_old = as.character(hospitalid),
             outcome_old = as.integer(outcome_value),
             p_base = as.numeric(prediction))]
new <- new[outcome == "mortality",
           .(patientunitstayid, patient_new = as.character(uniquepid),
             hospital_new = as.character(hospitalid),
             outcome_new = as.integer(outcome_value),
             l1_k2 = as.integer(l1_k2), p_k2 = as.numeric(prediction))]
stopifnot(nrow(old) == 17284L, nrow(new) == 17284L,
          !anyDuplicated(old$patientunitstayid),
          !anyDuplicated(new$patientunitstayid))
d <- merge(old, new, by = "patientunitstayid", all = TRUE, sort = TRUE)
stopifnot(nrow(d) == 17284L, !anyNA(d),
          identical(d$patient_old, d$patient_new),
          identical(d$hospital_old, d$hospital_new),
          identical(d$outcome_old, d$outcome_new),
          all(d$outcome_old %in% 0:1),
          all(d$l1_k2 %in% 1:2),
          all(d$p_base >= 0 & d$p_base <= 1),
          all(d$p_k2 >= 0 & d$p_k2 <= 1),
          uniqueN(d$patient_old) == 16058L,
          uniqueN(d$hospital_old) == 199L,
          sum(d$outcome_old) == 2940L)
setnames(d, c("patient_old", "hospital_old", "outcome_old"),
         c("patient", "hospital", "y"))
d[, c("patient_new", "hospital_new", "outcome_new") := NULL]

old_ref <- fread(old_perf_path)[outcome == "mortality" & model == "base"]
new_ref <- fread(new_perf_path)[outcome == "mortality" & model == "base_l1_k2"]
stopifnot(nrow(old_ref) == 1L, nrow(new_ref) == 1L)
baseline <- rbindlist(lapply(c("base", "base_l1_k2"), function(model) {
  p <- d[[if (model == "base") "p_base" else "p_k2"]]
  ref <- if (model == "base") old_ref else new_ref
  cal <- calibration_slope_intercept(d$y, p)
  calculated <- c(auc = fast_auc(d$y, p),
                  brier = mean((p - d$y)^2),
                  citl = unname(cal["intercept"]),
                  slope = unname(cal["slope"]))
  archived <- c(auc = ref$external_auc, brier = ref$brier_score,
                citl = ref$calibration_intercept,
                slope = ref$calibration_slope)
  stopifnot(isTRUE(all.equal(calculated, archived, tolerance = 1e-9)))
  data.table(model = model, full_eicu_stays = nrow(d),
             full_eicu_patients = uniqueN(d$patient),
             full_eicu_hospitals = uniqueN(d$hospital),
             full_eicu_events = sum(d$y),
             auc = calculated["auc"], brier = calculated["brier"],
             citl = calculated["citl"], slope = calculated["slope"],
             archived_values_reproduced = TRUE)
}))

# Connected hospitals prevent a person with records at two hospitals from
# appearing on both sides of the held-out comparison.
hospitals <- sort(unique(d$hospital))
parent <- seq_along(hospitals)
find_root <- function(i) {
  while (parent[i] != i) {
    parent[i] <<- parent[parent[i]]
    i <- parent[i]
  }
  i
}
join_roots <- function(i, j) {
  ri <- find_root(i)
  rj <- find_root(j)
  if (ri != rj) parent[rj] <<- ri
}
patient_hospital <- unique(d[, .(patient, hospital)])
for (h in split(patient_hospital$hospital, patient_hospital$patient)) {
  if (length(h) > 1L) {
    ids <- match(h, hospitals)
    for (j in ids[-1L]) join_roots(ids[1L], j)
  }
}
components <- vapply(seq_along(hospitals), find_root, integer(1))
component_ids <- sort(unique(components))
component_hospital_n <- tabulate(match(components, component_ids),
                                 nbins = length(component_ids))
stopifnot(length(component_ids) >= 2L)
set.seed(20260926L)
permutation <- sample(seq_along(component_ids))
hospital_cumsum <- cumsum(component_hospital_n[permutation])
candidate <- seq_len(length(permutation) - 1L)
selected_prefix <- candidate[which.min(abs(hospital_cumsum[candidate] -
                                         0.7 * length(hospitals)))]
update_components <- component_ids[permutation[seq_len(selected_prefix)]]
hospital_component <- data.table(hospital = hospitals,
                                 component = component_ids[
                                   match(components, component_ids)])
hospital_component[, side := ifelse(component %in% update_components,
                                     "updating", "evaluation")]
d <- merge(d, hospital_component, by = "hospital", sort = FALSE)
stopifnot(!anyNA(d$side), all(c("updating", "evaluation") %in% d$side))
up <- d[side == "updating"]
ev <- d[side == "evaluation"]
stopifnot(length(intersect(unique(up$hospital), unique(ev$hospital))) == 0L,
          length(intersect(unique(up$patient), unique(ev$patient))) == 0L,
          all(c(0L, 1L) %in% up$y), all(c(0L, 1L) %in% ev$y))
split_summary <- d[, .(
  hospital_n = uniqueN(hospital), component_n = uniqueN(component),
  patient_n = uniqueN(patient), stay_n = .N, event_n = sum(y),
  c1_stay_n = sum(l1_k2 == 1L), c2_stay_n = sum(l1_k2 == 2L),
  repeat_stays = .N - uniqueN(patient),
  base_p_q25 = unname(quantile(p_base, 0.25)),
  base_p_median = median(p_base), base_p_q75 = unname(quantile(p_base, 0.75)),
  k2_p_q25 = unname(quantile(p_k2, 0.25)),
  k2_p_median = median(p_k2), k2_p_q75 = unname(quantile(p_k2, 0.75))
), by = side]
split_summary[, connected_components_all := length(component_ids)]
split_summary[, cross_hospital_patient_n := uniqueN(
  patient_hospital[, .N, by = patient][N > 1L, patient])]

eps <- 1e-6
clamp <- function(p) pmin(pmax(p, eps), 1 - eps)
boundary_counts <- rbindlist(lapply(c("base", "base_l1_k2"), function(model) {
  p <- d[[if (model == "base") "p_base" else "p_k2"]]
  data.table(model = model, p_zero_n = sum(p == 0), p_one_n = sum(p == 1),
             p_clamped_n = sum(p <= eps | p >= 1 - eps))
}))
d[, lp_base := qlogis(clamp(p_base))]
d[, lp_k2 := qlogis(clamp(p_k2))]
up <- d[side == "updating"]
ev <- d[side == "evaluation"]

fit_glm <- function(y, lp, w, offset_only = FALSE) {
  active <- w > 0
  if (!all(c(0L, 1L) %in% y[active])) stop("one outcome class in sample")
  x <- if (offset_only) matrix(1, nrow = length(y), ncol = 1L) else
    cbind(1, lp)
  warning_text <- character()
  z <- withCallingHandlers(
    stats::glm.fit(x = x, y = y, weights = w,
                   offset = if (offset_only) lp else rep(0, length(lp)),
                   family = stats::binomial(),
                   control = stats::glm.control(maxit = 50)),
    warning = function(e) {
      warning_text <<- c(warning_text, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(z$converged) || any(!is.finite(z$coefficients)))
    stop("logistic fit did not converge with finite coefficients")
  list(coef = unname(z$coefficients), warning = paste(unique(warning_text),
                                                     collapse = " | "))
}
weighted_auc <- function(y, p, w) {
  keep <- w > 0
  y <- y[keep]; p <- p[keep]; w <- w[keep]
  positive <- sum(w[y == 1L]); negative <- sum(w[y == 0L])
  if (positive == 0 || negative == 0) return(NA_real_)
  ord <- order(p)
  y <- y[ord]; p <- p[ord]; w <- w[ord]
  tie <- cumsum(c(TRUE, diff(p) != 0))
  yes <- as.vector(rowsum(w * y, tie, reorder = FALSE))
  no <- as.vector(rowsum(w * (1L - y), tie, reorder = FALSE))
  sum(yes * (cumsum(no) - no + 0.5 * no)) / (positive * negative)
}
stopifnot(abs(weighted_auc(d$y, d$p_base, rep(1, nrow(d))) -
              fast_auc(d$y, d$p_base)) < 1e-12)
evaluate <- function(y, original_lp, p, w) {
  predicted_lp <- qlogis(clamp(p))
  citl <- tryCatch(fit_glm(y, predicted_lp, w, TRUE)$coef[1],
                   error = function(e) NA_real_)
  slope <- tryCatch(fit_glm(y, predicted_lp, w, FALSE)$coef[2],
                    error = function(e) NA_real_)
  c(mean_prediction = weighted.mean(p, w),
    observed_rate = weighted.mean(y, w),
    auc = weighted_auc(y, p, w),
    brier = weighted.mean((p - y)^2, w),
    citl = citl, slope = slope)
}
predict_modes <- function(lp, intercept_fit, logistic_fit) {
  list(no_update = plogis(lp),
       intercept_only = plogis(lp + intercept_fit$coef[1]),
       intercept_slope = plogis(logistic_fit$coef[1] +
                                 logistic_fit$coef[2] * lp))
}
stopifnot(max(abs(plogis(0 + 1 * d$lp_base) - plogis(d$lp_base))) == 0,
          all(plogis(d$lp_base) >= 0 & plogis(d$lp_base) <= 1))

point_rows <- list()
subgroup_rows <- list()
coefficient_rows <- list()
curve_rows <- list()
distribution_rows <- list()
model_names <- c("base", "base_l1_k2")
for (model in model_names) {
  lp_col <- if (model == "base") "lp_base" else "lp_k2"
  up_lp <- up[[lp_col]]; ev_lp <- ev[[lp_col]]
  intercept_fit <- fit_glm(up$y, up_lp, rep(1, nrow(up)), TRUE)
  logistic_fit <- fit_glm(up$y, up_lp, rep(1, nrow(up)), FALSE)
  coefficients <- data.table(model = model,
    update_hospitals = uniqueN(up$hospital), update_stays = nrow(up),
    update_events = sum(up$y),
    intercept_only_a0 = intercept_fit$coef[1],
    logistic_a = logistic_fit$coef[1], logistic_b = logistic_fit$coef[2],
    intercept_only_warning = intercept_fit$warning,
    logistic_warning = logistic_fit$warning)
  coefficient_rows[[model]] <- coefficients
  predictions <- predict_modes(ev_lp, intercept_fit, logistic_fit)
  for (method in names(predictions)) {
    p <- predictions[[method]]
    if (method != "no_update" &&
        (method == "intercept_only" || logistic_fit$coef[2] > 0)) {
      stopifnot(abs(fast_auc(ev$y, p) -
                    fast_auc(ev$y, predictions$no_update)) < 1e-12)
    }
    v <- evaluate(ev$y, ev_lp, p, rep(1, nrow(ev)))
    point_rows[[paste(model, method)]] <- data.table(
      model = model, method = method, evaluation_stays = nrow(ev),
      evaluation_patients = uniqueN(ev$patient),
      evaluation_hospitals = uniqueN(ev$hospital),
      evaluation_events = sum(ev$y),
      mean_prediction = v["mean_prediction"],
      observed_rate = v["observed_rate"], auc = v["auc"],
      brier = v["brier"], citl = v["citl"], slope = v["slope"])
    for (group in 1:2) {
      idx <- which(ev$l1_k2 == group)
      sg <- evaluate(ev$y[idx], ev_lp[idx], p[idx], rep(1, length(idx)))
      subgroup_rows[[paste(model, method, group)]] <- data.table(
        model = model, method = method, l1_group = paste0("C", group),
        stays = length(idx), patients = uniqueN(ev$patient[idx]),
        events = sum(ev$y[idx]), brier = sg["brier"],
        citl = sg["citl"], slope = sg["slope"],
        auc = sg["auc"])
    }
    br <- seq(0, 1, by = 0.05)
    hist <- tabulate(findInterval(p, br, rightmost.closed = TRUE),
                     nbins = length(br) - 1L)
    distribution_rows[[paste(model, method)]] <- data.table(
      model = model, method = method, bin_low = head(br, -1L),
      bin_high = tail(br, -1L), stays = hist)
  }
  basis <- splines::ns(ev_lp, df = 3)
  curve_fit <- stats::glm.fit(x = cbind(1, basis), y = ev$y,
                              family = stats::binomial())
  stopifnot(isTRUE(curve_fit$converged),
            all(is.finite(curve_fit$coefficients)))
  grid_lp <- seq(min(ev_lp), max(ev_lp), length.out = 201L)
  smooth_rate <- plogis(drop(cbind(1, predict(basis, grid_lp)) %*%
                               curve_fit$coefficients))
  for (method in names(predictions)) {
    grid_risk <- switch(method,
      no_update = plogis(grid_lp),
      intercept_only = plogis(grid_lp + intercept_fit$coef[1]),
      intercept_slope = plogis(logistic_fit$coef[1] +
                                 logistic_fit$coef[2] * grid_lp))
    curve_rows[[paste(model, method)]] <- data.table(
      model = model, method = method, grid_lp = grid_lp,
      predicted_risk = grid_risk, observed_smooth = smooth_rate)
  }
}
point <- rbindlist(point_rows)
subgroups <- rbindlist(subgroup_rows)
coefficients <- rbindlist(coefficient_rows)
curves <- rbindlist(curve_rows)
distributions <- rbindlist(distribution_rows)
point[, delta_brier_vs_no_update := brier - brier[method == "no_update"],
      by = model]

# The same hospital-component draws are used for both source models and all
# methods. Frequency weights retain a repeatedly drawn component intact.
up_components <- sort(unique(up$component))
ev_components <- sort(unique(ev$component))
up_comp_idx <- match(up$component, up_components)
ev_comp_idx <- match(ev$component, ev_components)
set.seed(20261026L)
B <- 1000L
boot_rows <- vector("list", B * length(model_names) * 3L)
boot_coeff <- vector("list", B * length(model_names))
failure_rows <- list()
row_id <- 0L
coef_id <- 0L
for (b in seq_len(B)) {
  up_draw <- sample.int(length(up_components), length(up_components),
                        replace = TRUE)
  ev_draw <- sample.int(length(ev_components), length(ev_components),
                        replace = TRUE)
  up_w <- tabulate(up_draw, nbins = length(up_components))[up_comp_idx]
  ev_w <- tabulate(ev_draw, nbins = length(ev_components))[ev_comp_idx]
  for (model in model_names) {
    lp_col <- if (model == "base") "lp_base" else "lp_k2"
    result <- tryCatch({
      f0 <- fit_glm(up$y, up[[lp_col]], up_w, TRUE)
      f1 <- fit_glm(up$y, up[[lp_col]], up_w, FALSE)
      predictions <- predict_modes(ev[[lp_col]], f0, f1)
      metrics <- lapply(predictions, function(p) {
        evaluate(ev$y, ev[[lp_col]], p, ev_w)
      })
      list(f0 = f0, f1 = f1, metrics = metrics)
    }, error = function(e) e)
    if (inherits(result, "error")) {
      failure_rows[[length(failure_rows) + 1L]] <- data.table(
        replicate = b, model = model, reason = conditionMessage(result))
      next
    }
    coef_id <- coef_id + 1L
    boot_coeff[[coef_id]] <- data.table(replicate = b, model = model,
      intercept_only_a0 = result$f0$coef[1],
      logistic_a = result$f1$coef[1], logistic_b = result$f1$coef[2],
      intercept_only_warning = result$f0$warning,
      logistic_warning = result$f1$warning)
    ref_brier <- result$metrics$no_update["brier"]
    for (method in names(result$metrics)) {
      row_id <- row_id + 1L
      v <- result$metrics[[method]]
      boot_rows[[row_id]] <- data.table(
        replicate = b, model = model, method = method,
        mean_prediction = v["mean_prediction"],
        observed_rate = v["observed_rate"], auc = v["auc"],
        brier = v["brier"], citl = v["citl"], slope = v["slope"],
        delta_brier_vs_no_update = v["brier"] - ref_brier)
    }
  }
  if (b %% 100L == 0L) message("Bootstrap attempt ", b, "/", B)
}
bootstrap <- rbindlist(boot_rows[seq_len(row_id)])
boot_coefficient <- rbindlist(boot_coeff[seq_len(coef_id)])
failures <- if (length(failure_rows)) rbindlist(failure_rows) else
  data.table(replicate = integer(), model = character(), reason = character())
ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10L) return(c(low = NA_real_, high = NA_real_, n = length(x)))
  c(low = unname(quantile(x, 0.025)),
    high = unname(quantile(x, 0.975)), n = length(x))
}
metric_names <- c("mean_prediction", "observed_rate", "auc", "brier",
                  "citl", "slope", "delta_brier_vs_no_update")
intervals <- rbindlist(lapply(seq_len(nrow(point)), function(i) {
  z <- point[i]
  q <- bootstrap[model == z$model & method == z$method]
  rbindlist(lapply(metric_names, function(metric) {
    bounds <- ci(q[[metric]])
    data.table(model = z$model, method = z$method, metric = metric,
      estimate = z[[metric]], ci_low = bounds["low"],
      ci_high = bounds["high"], valid_replicates = bounds["n"],
      attempted_replicates = B, failed_replicates = B - bounds["n"])
  }))
}))
coefficient_intervals <- rbindlist(lapply(seq_len(nrow(coefficients)),
  function(i) {
    z <- coefficients[i]
    q <- boot_coefficient[model == z$model]
    rbindlist(lapply(c("intercept_only_a0", "logistic_a", "logistic_b"),
      function(term) {
        bounds <- ci(q[[term]])
        data.table(model = z$model, term = term,
          estimate = z[[term]], ci_low = bounds["low"],
          ci_high = bounds["high"], valid_replicates = bounds["n"],
          attempted_replicates = B)
      }))
  }))
stopifnot(nrow(bootstrap) + 3L * nrow(failures) == B * 2L * 3L)

fwrite(baseline, file.path(aggregate, "full_eicu_original_reproduction.csv"))
fwrite(split_summary, file.path(aggregate, "hospital_holdout_split_summary.csv"))
fwrite(boundary_counts, file.path(aggregate, "prediction_boundary_counts.csv"))
fwrite(point, file.path(aggregate, "heldout_performance_point.csv"))
fwrite(intervals, file.path(aggregate, "heldout_performance_intervals.csv"))
fwrite(subgroups, file.path(aggregate, "heldout_L1_subgroups.csv"))
fwrite(coefficients, file.path(aggregate, "updating_parameters_point.csv"))
fwrite(coefficient_intervals,
       file.path(aggregate, "updating_parameters_intervals.csv"))
fwrite(curves, file.path(aggregate, "heldout_calibration_curves.csv"))
fwrite(distributions, file.path(aggregate, "heldout_risk_distributions.csv"))
fwrite(bootstrap, file.path(aggregate, "hospital_bootstrap_aggregate.csv"))
fwrite(failures, file.path(aggregate, "hospital_bootstrap_failures.csv"))
fwrite(d[, .(patientunitstayid, patient, hospital, component, side, l1_k2,
             y, p_base, p_k2)], file.path(private, "joined_holdout_records.csv"))
fwrite(hospital_component, file.path(private, "hospital_split_mapping.csv"))
writeLines(capture.output(sessionInfo()), file.path(qa, "sessionInfo.txt"))
message("D3 recalibration completed: ", nrow(ev),
        " held-out stays; ", nrow(failures), " bootstrap model failures")
