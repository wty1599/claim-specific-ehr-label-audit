options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

formal <- Sys.getenv("D4_FORMAL_OUTPUT_ROOT")
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT")
stopifnot(nzchar(formal), nzchar(repo))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(output_root, repo)
audit_path <- file.path(formal, "private_not_for_release", "D4_propensity_patient_audit_private.rds")
ledger_path <- file.path(formal, "private_not_for_release", "D4_landmark_patient_ledger_private.rds")
impute_path <- file.path(formal, "tables", "D4_propensity_formula_and_imputation.csv")
formal_support_path <- file.path(formal, "tables", "D4_support_metrics_by_K2.csv")
out_dir <- file.path(output_root, "domain4_arm_support", "empirical_v2")
stopifnot(all(file.exists(c(audit_path, ledger_path, impute_path, formal_support_path))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

audit <- as.data.table(readRDS(audit_path))
ledger <- as.data.table(readRDS(ledger_path))
audit[, cluster_k2 := as.integer(as.character(cluster_k2))]
ledger[, cluster_k2 := as.integer(as.character(cluster_k2))]
stopifnot(nrow(audit) == 17525L, nrow(ledger) == nrow(audit),
          uniqueN(audit$stay_id) == nrow(audit),
          uniqueN(ledger$stay_id) == nrow(ledger))

covariates <- c("age", "sofa_score", "aki_stage_0_24h", "sapsii",
                "creatinine_max", "bun_max", "potassium_max",
                "bicarbonate_min", "ph_min", "lactate_max",
                "urine_output_24h_ml", "mbp_min", "heart_rate_max",
                "resp_rate_max", "gcs_min")
needed <- c("stay_id", "cluster_k2", "rrt_first_24_72", "mortality_30d",
            covariates, "gender")
stopifnot(all(needed %in% names(ledger)))
z <- merge(audit, ledger[, ..needed], by = "stay_id", suffixes = c("_audit", "_ledger"), sort = FALSE)
stopifnot(nrow(z) == nrow(audit),
          all(z$cluster_k2_audit == z$cluster_k2_ledger),
          all(z$rrt_first_24_72_audit == z$rrt_first_24_72_ledger),
          !anyNA(z$mortality_30d_audit),
          !anyNA(z$mortality_30d_ledger),
          all(z$mortality_30d_audit == z$mortality_30d_ledger))
setnames(z, c("cluster_k2_audit", "rrt_first_24_72_audit", "mortality_30d_audit"),
         c("cluster_k2", "A", "death30"))
z[, c("cluster_k2_ledger", "rrt_first_24_72_ledger", "mortality_30d_ledger") := NULL]
stopifnot(all(z$A %in% 0:1), all(z$cluster_k2 %in% 1:2),
          all(is.finite(z$ps) & z$ps > 0 & z$ps < 1),
          all(is.finite(z$ipw_ate) & z$ipw_ate > 0))

fill <- fread(impute_path)
missing_rows <- rbindlist(lapply(c(covariates, "gender"), function(v) {
  x <- z[[v]]
  data.table(variable = v, n_missing_before = sum(is.na(x) |
    (is.character(x) & !nzchar(x))))
}))
missing_arm_rows <- rbindlist(lapply(1:2, function(k) {
  rbindlist(lapply(0:1, function(a) {
    za <- z[cluster_k2 == k & A == a]
    rbindlist(lapply(c(covariates, "gender"), function(v) {
      x <- za[[v]]
      data.table(phenotype = if (k == 1L) "C1 higher-risk" else "C2 lower-risk",
                 treatment = a, variable = v, n = nrow(za),
                 n_missing_before = sum(is.na(x) |
                   (is.character(x) & !nzchar(x))))
    }))
  }))
}))
for (v in covariates) {
  fv <- fill[variable == v, fill_value]
  stopifnot(length(fv) == 1L, is.finite(as.numeric(fv)))
  set(z, i = which(is.na(z[[v]])), j = v, value = as.numeric(fv))
}
z[, gender_male := as.numeric(gender == "M")]
stopifnot(!anyNA(z[, c(covariates, "gender_male"), with = FALSE]))
fwrite(missing_rows, file.path(out_dir, "covariate_missing_before_imputation.csv"))
fwrite(missing_arm_rows,
       file.path(out_dir, "covariate_missing_by_arm_before_imputation.csv"))

expected_ate <- fifelse(z$A == 1L, 1 / z$ps, 1 / (1 - z$ps))
stopifnot(max(abs(expected_ate - z$ipw_ate)) < 1e-8)
z[, `:=`(w_ATE = ipw_ate,
         w_ATT = fifelse(A == 1L, 1, ps / (1 - ps)),
         w_ATO = fifelse(A == 1L, 1 - ps, ps))]
weight_names <- c(ATE = "w_ATE", ATT = "w_ATT", ATO = "w_ATO")

ess <- function(w) {
  if (!length(w) || any(!is.finite(w) | w <= 0) || sum(w) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}
weighted_avg <- function(x, w) {
  if (!length(x) || anyNA(x) || any(!is.finite(x)) || any(!is.finite(w) | w < 0) ||
      sum(w) <= 0) return(NA_real_)
  sum(x * w) / sum(w)
}

group_rows <- list()
pooled_rows <- list()
score_rows <- list()
balance_rows <- list()
target_rows <- list()
formal_support <- fread(formal_support_path)
balance_vars <- c(covariates, "gender_male")
target_vars <- c("age", "sofa_score", "creatinine_max", "aki_stage_0_24h")

for (k in 1:2) {
  d <- z[cluster_k2 == k]
  name <- if (k == 1L) "C1 higher-risk" else "C2 lower-risk"
  n1 <- sum(d$A == 1L)
  n0 <- sum(d$A == 0L)
  n <- nrow(d)
  stopifnot(n1 > 0, n0 > 0)
  p <- n1 / n
  original <- formal_support[phenotype == name]
  overlap_low <- max(min(d[A == 1L, ps]), min(d[A == 0L, ps]))
  overlap_high <- min(max(d[A == 1L, ps]), max(d[A == 0L, ps]))
  stopifnot(nrow(original) == 1L,
            abs(overlap_low - original$common_support_low) < 1e-8,
            abs(overlap_high - original$common_support_high) < 1e-8,
            abs(mean(d$ps < overlap_low | d$ps > overlap_high) -
                original$outside_overlap_proportion) < 1e-9)

  w_constant <- fifelse(d$A == 1L, 1 / p, 1 / (1 - p))
  constant_ratio_numeric <- ess(w_constant) / n
  constant_ratio_formula <- 4 * p * (1 - p)
  stopifnot(abs(constant_ratio_numeric - constant_ratio_formula) < 1e-12)
  pooled_rows[[length(pooled_rows) + 1L]] <- data.table(
    phenotype = name, n = n, n_rrt = n1, n_no_rrt = n0, observed_p = p,
    observed_unstabilized_ate_pooled_ess = ess(d$w_ATE),
    observed_unstabilized_ate_pooled_ess_per_n = ess(d$w_ATE) / n,
    constant_score_pooled_ess_per_n_numeric = constant_ratio_numeric,
    constant_score_reference_4p1minus_p = constant_ratio_formula,
    constant_score_reference_is_patient_estimate = FALSE
  )

  for (a in 0:1) {
    da <- d[A == a]
    qa <- as.numeric(quantile(da$ps, probs = c(0, .01, .05, .25, .5,
                                           .75, .95, .99, 1), names = FALSE))
    score_rows[[length(score_rows) + 1L]] <- data.table(
      phenotype = name, treatment = a, n = nrow(da),
      ps_min = qa[1], ps_p01 = qa[2], ps_p05 = qa[3], ps_p25 = qa[4],
      ps_median = qa[5], ps_p75 = qa[6], ps_p95 = qa[7],
      ps_p99 = qa[8], ps_max = qa[9],
      empirical_overlap_low = overlap_low,
      empirical_overlap_high = overlap_high,
      n_below = sum(da$ps < overlap_low),
      n_above = sum(da$ps > overlap_high),
      fraction_outside = mean(da$ps < overlap_low | da$ps > overlap_high)
    )
  }

  for (target in names(weight_names)) {
    wname <- weight_names[[target]]
    w <- d[[wname]]
    arm_ess <- numeric(2)
    for (a in 0:1) {
      da <- d[A == a]
      wa <- da[[wname]]
      ea <- ess(wa)
      arm_ess[a + 1L] <- ea
      shares <- sort(wa / sum(wa), decreasing = TRUE)
      top_n <- ceiling(0.01 * nrow(da))
      stopifnot(is.finite(ea), ea <= nrow(da) + 1e-8,
                abs(sum(shares) - 1) < 1e-12)
      group_rows[[length(group_rows) + 1L]] <- data.table(
        phenotype = name, target = target, treatment = a,
        n = nrow(da), deaths_30d = sum(da$death30 == 1L, na.rm = TRUE),
        survived_30d = sum(da$death30 == 0L, na.rm = TRUE),
        outcome_unknown = sum(is.na(da$death30)),
        weight_sum = sum(wa), weight_max = max(wa),
        arm_ess = ea, arm_ess_per_n = ea / nrow(da),
        maximum_normalized_weight_share = shares[1L],
        top_one_percent_n = top_n,
        top_one_percent_normalized_weight_share = sum(shares[seq_len(top_n)]),
        negative_urine_n = sum(da$urine_output_24h_ml < 0),
        negative_urine_min_ml = if (any(da$urine_output_24h_ml < 0))
          min(da$urine_output_24h_ml[da$urine_output_24h_ml < 0]) else NA_real_,
        negative_urine_max_ml = if (any(da$urine_output_24h_ml < 0))
          max(da$urine_output_24h_ml[da$urine_output_24h_ml < 0]) else NA_real_,
        negative_urine_weight_share =
          sum(wa[da$urine_output_24h_ml < 0]) / sum(wa)
      )
    }
    deff <- (1 / arm_ess[2L] + 1 / arm_ess[1L]) / (1 / n1 + 1 / n0)
    pooled_rows[[length(pooled_rows)]][, paste0("deff_", target) := deff]

    for (v in balance_vars) {
      x1 <- d[A == 1L][[v]]
      x0 <- d[A == 0L][[v]]
      scale_ref <- sqrt((var(x1) + var(x0)) / 2)
      treated_mean <- weighted_avg(x1, d[A == 1L][[wname]])
      control_mean <- weighted_avg(x0, d[A == 0L][[wname]])
      smd <- if (is.finite(scale_ref) && scale_ref > 0) {
        (treated_mean - control_mean) / scale_ref
      } else NA_real_
      balance_rows[[length(balance_rows) + 1L]] <- data.table(
        phenotype = name, target = target, variable = v,
        unweighted_reference_sd = scale_ref,
        treated_mean = treated_mean, untreated_mean = control_mean,
        smd_fixed_unweighted_scale = smd,
        unweighted_treated_mean = mean(x1),
        unweighted_untreated_mean = mean(x0),
        unweighted_smd_fixed_scale = if (is.finite(scale_ref) && scale_ref > 0)
          (mean(x1) - mean(x0)) / scale_ref else NA_real_,
        data_quality_flag = if (v == "urine_output_24h_ml")
          "Negative source values; do not interpret urine balance clinically"
          else ""
      )
    }
    for (v in target_vars) {
      target_rows[[length(target_rows) + 1L]] <- data.table(
        phenotype = name, target = target, variable = v,
        observed_unweighted_mean = mean(d[[v]]),
        observed_weighted_mixture_mean = weighted_avg(d[[v]], w),
        treated_weighted_mean = weighted_avg(d[A == 1L][[v]], d[A == 1L][[wname]]),
        untreated_weighted_mean = weighted_avg(d[A == 0L][[v]], d[A == 0L][[wname]])
      )
    }
  }
}

arm_summary <- rbindlist(group_rows)
pooled_summary <- rbindlist(pooled_rows, fill = TRUE)
scores <- rbindlist(score_rows)
balance <- rbindlist(balance_rows)
target_population <- rbindlist(target_rows)
stopifnot(nrow(arm_summary) == 12L, nrow(pooled_summary) == 2L,
          nrow(scores) == 4L, nrow(balance) == 2L * 3L * length(balance_vars))

# Mathematical identity and invariance tests use synthetic weights, not patient data.
toy_A <- c(1L, 1L, 0L, 0L, 0L)
toy_w <- c(2, 5, 1, 3, 4)
toy_scaled <- toy_w * fifelse(toy_A == 1L, 7, 2)
same_arm_ess <- all(vapply(0:1, function(a)
  abs(ess(toy_w[toy_A == a]) - ess(toy_scaled[toy_A == a])) < 1e-12,
  logical(1)))
same_arm_shares <- all(vapply(0:1, function(a) {
  old <- toy_w[toy_A == a] / sum(toy_w[toy_A == a])
  new <- toy_scaled[toy_A == a] / sum(toy_scaled[toy_A == a])
  max(abs(old - new)) < 1e-12
}, logical(1)))
pooled_changed <- abs(ess(toy_w) - ess(toy_scaled)) > 1e-6
p_toy <- mean(toy_A)
constant_unstabilized <- fifelse(toy_A == 1L, 1 / p_toy, 1 / (1 - p_toy))
constant_stabilized <- constant_unstabilized * fifelse(toy_A == 1L, p_toy, 1 - p_toy)
toy_ps <- c(.2, .8, .3, .6, .4)
toy_att <- fifelse(toy_A == 1L, 1, toy_ps / (1 - toy_ps))
toy_ato <- fifelse(toy_A == 1L, 1 - toy_ps, toy_ps)
toy_smd <- (weighted_avg(c(2, 6), c(1, 3)) -
              weighted_avg(c(1, 3, 5), c(1, 1, 1))) / sqrt(6)
toy_deff <- (1 / ess(c(1, 3)) + 1 / ess(c(1, 1, 1))) / (1 / 2 + 1 / 3)
tests <- data.table(
  check = c("constant_score_identity", "arm_ess_rescaling_invariant",
            "arm_normalized_shares_rescaling_invariant", "pooled_ess_can_change",
            "stabilized_constant_weights_all_one", "stabilized_does_not_add_treated",
            "equal_arm_weights_ess_equals_n", "empty_arm_ess_uncomputable",
            "zero_weight_ess_uncomputable", "negative_weight_ess_uncomputable",
            "nonfinite_weight_ess_uncomputable", "att_formula",
            "ato_formula", "fixed_scale_smd_formula", "deff_formula"),
  pass = c(all(abs(pooled_summary$constant_score_pooled_ess_per_n_numeric -
                     pooled_summary$constant_score_reference_4p1minus_p) < 1e-12),
           same_arm_ess, same_arm_shares, pooled_changed,
           all(abs(constant_stabilized - 1) < 1e-12),
           sum(toy_A) == 2L,
           abs(ess(rep(3, 10)) - 10) < 1e-12,
           is.na(ess(numeric(0))), is.na(ess(c(0, 1))),
           is.na(ess(c(-1, 2))), is.na(ess(c(1, Inf))),
           max(abs(toy_att - c(1, 1, 3/7, 1.5, 2/3))) < 1e-12,
           max(abs(toy_ato - c(.8, .2, .3, .6, .4))) < 1e-12,
           abs(toy_smd - 2 / sqrt(6)) < 1e-12,
           abs(toy_deff - 1.15) < 1e-12)
)
stopifnot(all(tests$pass))

private_dir <- file.path(out_dir, "private_not_for_release")
dir.create(private_dir, showWarnings = FALSE)
fwrite(arm_summary, file.path(private_dir, "arm_ess_concentration_by_target_full.csv"))
public_arm_cols <- c("phenotype", "target", "treatment", "n", "deaths_30d",
                     "survived_30d", "outcome_unknown", "weight_sum",
                     "arm_ess", "arm_ess_per_n")
fwrite(arm_summary[, ..public_arm_cols],
       file.path(out_dir, "arm_ess_concentration_by_target.csv"))
fwrite(pooled_summary, file.path(out_dir, "pooled_ess_constant_reference.csv"))
fwrite(scores, file.path(out_dir, "propensity_overlap_by_arm.csv"))
fwrite(balance, file.path(out_dir, "covariate_balance_fixed_scale.csv"))
fwrite(target_population, file.path(out_dir, "weighted_population_descriptives.csv"))
fwrite(tests, file.path(out_dir, "formula_unit_tests.csv"))
fwrite(data.table(
  source = c(audit_path, ledger_path, impute_path, formal_support_path),
  sha256 = vapply(c(audit_path, ledger_path, impute_path, formal_support_path),
                  function(p) toupper(digest(p, algo = "sha256", file = TRUE,
                                           serialize = FALSE)), character(1)),
  output_contains_individuals = FALSE
), file.path(out_dir, "input_provenance.csv"))

cat("Empirical support diagnostics written.\n")
print(arm_summary[, .(phenotype, target, treatment, n, arm_ess, arm_ess_per_n,
                       maximum_normalized_weight_share,
                       top_one_percent_normalized_weight_share)])
print(pooled_summary)
