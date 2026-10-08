options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(file_arg) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
out <- file.path(root, "outputs")
private <- file.path(root, "private_not_for_release")
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
stopifnot(sha(file.path(root, "01_LOCKED_POSTHOC_SPEC.md")) ==
            "95A3A1FF8ACA21B7992AB7438E90E65C5ABFCDB36F4AD4393A9B429953C7115C")

project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
d4 <- file.path(project, "d4_landmark_support_interface_20260729")
paths <- c(
  risk_ledger = file.path(d4, "output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds"),
  original_ps = file.path(d4, "output_formal_v1/private_not_for_release/D4_propensity_patient_audit_private.rds"),
  original_imputation = file.path(d4, "output_formal_v1/tables/D4_propensity_formula_and_imputation.csv"),
  original_common_script = file.path(d4, "scripts/38_d4_landmark_support_common_v1.R"),
  ablation_labels = file.path(private, "A4_first_seed_labels_private.rds")
)
stopifnot(all(file.exists(paths)))
before_sha <- vapply(paths, sha, character(1))
stopifnot(before_sha[["risk_ledger"]] ==
            "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
          before_sha[["original_ps"]] ==
            "BC6106F70DD96457D4D3BA11915A89714EC2C8C581ABDC5F7DAF0C95FE72DC44",
          before_sha[["original_common_script"]] ==
            "2B6ECC405E15924CB69A77C7ACB3A2E88DE79D2E8DE99D397B0C9759DF1F78DD")

source(paths[["original_common_script"]], local = TRUE)
risk <- as.data.table(readRDS(paths[["risk_ledger"]]))
old <- as.data.table(readRDS(paths[["original_ps"]]))
labels <- as.data.table(readRDS(paths[["ablation_labels"]]))
impute_archive <- fread(paths[["original_imputation"]])
risk[, cluster_k2 := as.integer(as.character(cluster_k2))]
old[, cluster_k2 := as.integer(as.character(cluster_k2))]
stopifnot(nrow(risk) == 17525L, nrow(old) == nrow(risk),
          uniqueN(risk$stay_id) == nrow(risk),
          uniqueN(old$stay_id) == nrow(old),
          uniqueN(labels$stay_id) == nrow(labels),
          all(risk$rrt_first_24_72 %in% 0:1),
          all(risk$mortality_30d %in% 0:1),
          !anyNA(risk$mortality_30d))
idx <- match(risk$stay_id, labels$stay_id)
old_idx <- match(risk$stay_id, old$stay_id)
stopifnot(!anyNA(idx), !anyNA(old_idx),
          all(risk$cluster_k2 == labels$cluster_k2[idx]),
          all(risk$cluster_k2 == old$cluster_k2[old_idx]),
          all(risk$rrt_first_24_72 == old$rrt_first_24_72[old_idx]))
risk[, ablation_group := labels$ablation_group[idx]]
stopifnot(identical(as.integer(risk[cluster_k2 == 1L, .N]), 3151L),
          identical(as.integer(risk[cluster_k2 == 2L, .N]), 14374L),
          identical(as.integer(risk[cluster_k2 == 1L, sum(rrt_first_24_72)]), 265L),
          identical(as.integer(risk[cluster_k2 == 2L, sum(rrt_first_24_72)]), 49L),
          identical(as.integer(risk[cluster_k2 == 1L & rrt_first_24_72 == 1L,
                                   sum(mortality_30d)]), 149L),
          identical(as.integer(risk[cluster_k2 == 1L & rrt_first_24_72 == 0L,
                                   sum(mortality_30d)]), 992L),
          identical(as.integer(risk[cluster_k2 == 2L & rrt_first_24_72 == 1L,
                                   sum(mortality_30d)]), 18L),
          identical(as.integer(risk[cluster_k2 == 2L & rrt_first_24_72 == 0L,
                                   sum(mortality_30d)]), 1762L),
          all(tabulate(risk$ablation_group, 2L) >= 500L))

covars <- setdiff(D4_PS_COVARIATES, "cluster_k2")
stopifnot(length(D4_PS_COVARIATES) == 17L, length(covars) == 16L,
          !any(c("cluster_k2", "ablation_group") %in% covars),
          setequal(covars, c("age", "gender", "sofa_score", "aki_stage_0_24h",
                            "sapsii", "creatinine_max", "bun_max", "potassium_max",
                            "bicarbonate_min", "ph_min", "lactate_max",
                            "urine_output_24h_ml", "mbp_min", "heart_rate_max",
                            "resp_rate_max", "gcs_min")))
imp <- d4_impute_data(risk, covars, "median_mode")
dm <- as.data.table(imp$data)
set(dm, j = "gender", value = factor(dm$gender))
for (v in covars) {
  source_row <- impute_archive[variable == v]
  new_row <- imp$audit[variable == v]
  stopifnot(nrow(source_row) == 1L, nrow(new_row) == 1L,
            source_row$missing_n == new_row$missing_n,
            as.character(source_row$fill_value) == as.character(new_row$fill_value))
}

form <- reformulate(covars, response = "rrt_first_24_72")
stopifnot(!any(c("cluster_k2", "ablation_group", "mortality_30d") %in% all.vars(form)))
fit_warnings <- character()
fit <- withCallingHandlers(
  glm(form, data = dm, family = binomial()),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(isTRUE(fit$converged), all(is.finite(coef(fit))),
          nobs(fit) == nrow(dm))
ps_raw <- as.numeric(predict(fit, type = "response"))
eps <- 1e-6
ps <- pmin(pmax(ps_raw, eps), 1 - eps)
A <- as.integer(dm$rrt_first_24_72)
dm[, `:=`(ps_common = ps,
         w_ATE_common = fifelse(A == 1L, 1 / ps, 1 / (1 - ps)),
         w_ATO_common = fifelse(A == 1L, 1 - ps, ps))]
stopifnot(all(is.finite(dm$ps_common)), all(dm$ps_common > 0 & dm$ps_common < 1),
          all(is.finite(dm$w_ATE_common) & dm$w_ATE_common > 0),
          all(is.finite(dm$w_ATO_common) & dm$w_ATO_common > 0))

model_qc <- data.table(
  formula = paste(deparse(form), collapse = ""),
  n_input = nrow(risk), n_model = nobs(fit),
  n_record_positive = sum(A), converged = fit$converged,
  coefficient_count = length(coef(fit)),
  nonfinite_coefficient_count = sum(!is.finite(coef(fit))),
  warning_count = length(fit_warnings),
  warning_text = paste(unique(fit_warnings), collapse = " | "),
  probability_min_raw = min(ps_raw), probability_max_raw = max(ps_raw),
  clipped_probability_n = sum(ps_raw < eps | ps_raw > 1 - eps),
  pooled_allocation_auc = d4_binary_auc(A, ps),
  any_partition_in_formula = any(c("cluster_k2", "ablation_group") %in% all.vars(form)),
  effect_model_fitted = FALSE
)

ess <- function(w) {
  if (!length(w) || any(!is.finite(w) | w <= 0) || sum(w) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}
wmean <- function(x, w) sum(x * w) / sum(w)
balance_vars <- c(setdiff(covars, "gender"), "gender_male")
dm[, gender_male := as.numeric(gender == "M")]
stopifnot(length(balance_vars) == 16L, all(vapply(dm[, ..balance_vars],
                                                  function(x) all(is.finite(x)), logical(1))))
strict_block <- c("urine_output_24h_ml", "bun_max", "creatinine_max",
                  "ph_min", "pco2_max", "aniongap_max", "bicarbonate_min")
intergroup_vars <- unique(c(strict_block, "lactate_max", "potassium_min",
                            "potassium_max", "chloride_min", "chloride_max"))
full_matrix <- as.data.table(readRDS(file.path(project, "output/model/X_primary33_std_mice.rds")))
matrix_idx <- match(dm$stay_id, full_matrix$stay_id)
stopifnot(!anyNA(matrix_idx), all(intergroup_vars %in% names(full_matrix)))
intergroup_data <- data.table(cluster_k2 = dm$cluster_k2,
                              ablation_group = dm$ablation_group)
for (v in intergroup_vars) {
  set(intergroup_data, j = v, value = full_matrix[[v]][matrix_idx])
}
stopifnot(all(vapply(intergroup_data[, ..intergroup_vars],
                     function(x) all(is.finite(x)), logical(1))))

group_rows <- list()
balance_rows <- list()
premise_rows <- list()
decision_rows <- list()
for (partition in c("K2", "A4")) {
  group_col <- if (partition == "K2") "cluster_k2" else "ablation_group"
  for (g in 1:2) {
    d <- dm[get(group_col) == g]
    positive <- d[rrt_first_24_72 == 1L]
    negative <- d[rrt_first_24_72 == 0L]
    n1 <- nrow(positive); n0 <- nrow(negative)
    if (n1 == 0L || n0 == 0L) {
      stop("Cannot evaluate empty record-status cell in ", partition, " group ", g)
    }
    overlap_low <- max(min(positive$ps_common), min(negative$ps_common))
    overlap_high <- min(max(positive$ps_common), max(negative$ps_common))
    stopifnot(overlap_low <= overlap_high)

    row <- data.table(
      partition = partition, group = paste0(if (partition == "K2") "C" else "A", g),
      n = nrow(d), n_record_positive = n1, n_record_negative = n0,
      record_positive_rate = n1 / nrow(d),
      original_one_percent_rate_criterion = n1 / nrow(d) >= 0.01,
      positive_deaths_30d = sum(positive$mortality_30d == 1L),
      negative_deaths_30d = sum(negative$mortality_30d == 1L),
      positive_outcome_unknown = sum(is.na(positive$mortality_30d)),
      negative_outcome_unknown = sum(is.na(negative$mortality_30d)),
      overlap_low = overlap_low, overlap_high = overlap_high,
      positive_outside_n = sum(positive$ps_common < overlap_low |
                                 positive$ps_common > overlap_high),
      negative_outside_n = sum(negative$ps_common < overlap_low |
                                 negative$ps_common > overlap_high)
    )
    for (target in c("ATE", "ATO")) {
      wcol <- paste0("w_", target, "_common")
      w1 <- positive[[wcol]]; w0 <- negative[[wcol]]
      shares <- sort(w1 / sum(w1), decreasing = TRUE)
      n_top <- ceiling(0.01 * n1)
      set(row, j = paste0(target, "_positive_ESS"), value = ess(w1))
      set(row, j = paste0(target, "_negative_ESS"), value = ess(w0))
      set(row, j = paste0(target, "_positive_ESS_per_n"), value = ess(w1) / n1)
      set(row, j = paste0(target, "_negative_ESS_per_n"), value = ess(w0) / n0)
      set(row, j = paste0(target, "_top_one_percent_positive_n"), value = n_top)
      set(row, j = paste0(target, "_top_one_percent_positive_weight_share"),
          value = sum(shares[seq_len(n_top)]))
      stopifnot(ess(w1) <= n1 + 1e-8, ess(w0) <= n0 + 1e-8)
    }
    group_rows[[length(group_rows) + 1L]] <- row

    for (v in balance_vars) {
      x1 <- positive[[v]]; x0 <- negative[[v]]
      scale <- sqrt((var(x1) + var(x0)) / 2)
      smd <- function(w1, w0) {
        if (!is.finite(scale) || scale <= 0) return(NA_real_)
        (wmean(x1, w1) - wmean(x0, w0)) / scale
      }
      balance_rows[[length(balance_rows) + 1L]] <- data.table(
        partition = partition, group = row$group, variable = v,
        n_record_positive = n1, n_record_negative = n0,
        unweighted_reference_sd = scale,
        SMD_unweighted = smd(rep(1, n1), rep(1, n0)),
        SMD_ATE = smd(positive$w_ATE_common, negative$w_ATE_common),
        SMD_ATO = smd(positive$w_ATO_common, negative$w_ATO_common),
        interpretable = v != "urine_output_24h_ml"
      )
    }
  }

  g1 <- intergroup_data[get(group_col) == 1L]
  g2 <- intergroup_data[get(group_col) == 2L]
  for (v in intergroup_vars) {
    x1 <- g1[[v]]; x2 <- g2[[v]]
    scale <- sqrt((var(x1) + var(x2)) / 2)
    premise_rows[[length(premise_rows) + 1L]] <- data.table(
      partition = partition, variable = v,
      in_strict_removed_block = v %in% strict_block,
      group1_mean = mean(x1), group2_mean = mean(x2),
      unweighted_pooled_sd = scale,
      SMD_group1_minus_group2 = if (is.finite(scale) && scale > 0)
        (mean(x1) - mean(x2)) / scale else NA_real_,
      variable_scale = "standardized completed matrix"
    )
  }
}

groups <- rbindlist(group_rows, fill = TRUE)
balance <- rbindlist(balance_rows)
premise <- rbindlist(premise_rows)
stopifnot(nrow(groups) == 4L, nrow(balance) == 4L * 16L,
          nrow(premise) == 2L * length(intergroup_vars),
          all(is.finite(balance$SMD_unweighted)),
          all(is.finite(balance$SMD_ATE)), all(is.finite(balance$SMD_ATO)))

for (i in seq_len(nrow(groups))) {
  r <- groups[i]
  b <- balance[partition == r$partition & group == r$group & interpretable]
  stopifnot(nrow(b) == 15L)
  violations <- b[abs(SMD_ATO) >= 0.1, .(variable, SMD_ATO)]
  group_pass <- nrow(violations) == 0L &&
    r$positive_deaths_30d >= 10L && r$negative_deaths_30d >= 10L
  decision_rows[[i]] <- data.table(
    partition = r$partition, group = r$group,
    count_SMD_ATE_above_point_one = b[abs(SMD_ATE) > 0.1, .N],
    count_SMD_ATO_above_point_one = b[abs(SMD_ATO) > 0.1, .N],
    count_SMD_ATO_at_or_above_point_one = nrow(violations),
    variables_failing_ATO = paste(sprintf("%s:%.6f", violations$variable,
                                           violations$SMD_ATO), collapse = "; "),
    death_cell_minimum_met = r$positive_deaths_30d >= 10L &&
      r$negative_deaths_30d >= 10L,
    measured_balance_and_event_screen_met = group_pass
  )
}
decisions <- rbindlist(decision_rows)
overall <- decisions[, .(both_groups_meet_screen = all(measured_balance_and_event_screen_met)),
                     by = partition]
overall[, predicted_A4_both_groups_meet := ifelse(partition == "A4",
                                                  TRUE, NA)]
stopifnot(nrow(overall) == 2L, setequal(overall$partition, c("K2", "A4")),
          overall[partition == "A4", predicted_A4_both_groups_meet])

after_sha <- vapply(paths, sha, character(1))
stopifnot(identical(before_sha, after_sha))
manifest <- data.table(role = names(paths), path = unname(paths),
                       sha256_before = unname(before_sha),
                       sha256_after = unname(after_sha), unchanged = TRUE)

saveRDS(dm[, .(stay_id, cluster_k2, ablation_group, rrt_first_24_72,
               mortality_30d, ps_common, w_ATE_common, w_ATO_common)],
        file.path(private, "common_propensity_patient_audit_private.rds"))
fwrite(model_qc, file.path(out, "common_propensity_model_QC.csv"))
fwrite(imp$audit, file.path(out, "common_propensity_imputation.csv"))
fwrite(groups, file.path(out, "parallel_group_support.csv"))
fwrite(balance, file.path(out, "parallel_covariate_balance.csv"))
fwrite(premise, file.path(out, "between_group_physiology_SMD.csv"))
fwrite(decisions, file.path(out, "posthoc_group_screen.csv"))
fwrite(overall, file.path(out, "posthoc_between_group_screen.csv"))
fwrite(manifest, file.path(out, "source_manifest_comparison.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_session_common_score.txt"))
cat("Common score fit once; partition-neutral formula verified; A4 locked prediction=",
    overall[partition == "A4", predicted_A4_both_groups_meet],
    "; observed screen=", overall[partition == "A4", both_groups_meet_screen],
    "\n", sep = "")
