options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
parent <- dirname(root)
out <- file.path(root, "outputs")
private <- file.path(root, "private_not_for_release")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(private, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
stopifnot(sha(file.path(root, "01_REVISION_SPEC.md")) ==
            "5B51EC3B29B0D9E503099145E88D376C933283BF56AE75AC6FDE4B0757FF7984")

project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
d4 <- file.path(project, "d4_landmark_support_interface_20260729")
paths <- c(
  risk_ledger = file.path(d4, "output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds"),
  ablation_labels = file.path(parent, "private_not_for_release/A4_first_seed_labels_private.rds"),
  v1_patient_scores = file.path(parent, "private_not_for_release/common_propensity_patient_audit_private.rds"),
  original_common_script = file.path(d4, "scripts/38_d4_landmark_support_common_v1.R"),
  original_imputation = file.path(d4, "output_formal_v1/tables/D4_propensity_formula_and_imputation.csv"),
  v1_balance = file.path(parent, "outputs/parallel_covariate_balance.csv"),
  v1_groups = file.path(parent, "outputs/parallel_group_support.csv"),
  v1_screen = file.path(parent, "outputs/posthoc_group_screen.csv"),
  v1_between_screen = file.path(parent, "outputs/posthoc_between_group_screen.csv"),
  v1_sofa_audit = file.path(parent, "outputs/sofa_time_boundary_aggregate.csv")
)
stopifnot(all(file.exists(paths)))
before_sha <- vapply(paths, sha, character(1))
expected <- c(
  risk_ledger = "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
  ablation_labels = "0A0E9AFAB6D320B48CBD5C7A13B814A46C5C3F8A84DF47DC3A69763A13795987",
  v1_patient_scores = "CA2B096AE452776BFD243F2116A2036A64F92ACB9BC5BEDCA72BA7B70E459887",
  original_common_script = "2B6ECC405E15924CB69A77C7ACB3A2E88DE79D2E8DE99D397B0C9759DF1F78DD",
  original_imputation = "D85600EC7BFB2E4D6C949BD1E2FF432C0B47FDF76091885349A5FC89B0A3A019",
  v1_balance = "923C6D2511F92750DC69CF3B62E2C72337AB40F13AC37C7B0FB4D2BC10893B59",
  v1_groups = "164113E5F7AE7EA617E6B64D04399A445A4F2FF03CA417EFDFE5049C4E809954",
  v1_screen = "2A770AAE55F0A7363BC75F1FC50184267FD0188FEB4C27AB639D3D84C94D1602",
  v1_between_screen = "1D8B55AFBB810523A66C9E33522C8D4556B1790DFF2EABFE31560AC763498343",
  v1_sofa_audit = "84F59CD655FB47A26BF005C05E2C179D992AA6040D03394BB5AB57769BBDA14A"
)
stopifnot(all(before_sha[names(expected)] == expected))

source(paths[["original_common_script"]], local = TRUE)
risk <- as.data.table(readRDS(paths[["risk_ledger"]]))
labels <- as.data.table(readRDS(paths[["ablation_labels"]]))
risk[, cluster_k2 := as.integer(as.character(cluster_k2))]
stopifnot(nrow(risk) == 17525L, uniqueN(risk$stay_id) == nrow(risk),
          uniqueN(labels$stay_id) == nrow(labels),
          !anyNA(risk$rrt_first_24_72),
          all(risk$rrt_first_24_72 %in% 0:1),
          !anyNA(risk$mortality_30d), all(risk$mortality_30d %in% 0:1))
idx <- match(risk$stay_id, labels$stay_id)
stopifnot(!anyNA(idx), all(risk$cluster_k2 == labels$cluster_k2[idx]),
          identical(as.integer(risk[cluster_k2 == 1L, .N]), 3151L),
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
                                   sum(mortality_30d)]), 1762L))
risk[, ablation_group := as.integer(labels$ablation_group[idx])]
stopifnot(identical(as.integer(risk[ablation_group == 1L, .N]), 5168L),
          identical(as.integer(risk[ablation_group == 2L, .N]), 12357L))

covars <- setdiff(D4_PS_COVARIATES, c("cluster_k2", "sofa_score"))
stopifnot(length(D4_PS_COVARIATES) == 17L, length(covars) == 15L,
          setequal(covars, c("age", "gender", "aki_stage_0_24h", "sapsii",
                            "creatinine_max", "bun_max", "potassium_max",
                            "bicarbonate_min", "ph_min", "lactate_max",
                            "urine_output_24h_ml", "mbp_min", "heart_rate_max",
                            "resp_rate_max", "gcs_min")))
imp <- d4_impute_data(risk, covars, "median_mode")
dm <- as.data.table(imp$data)
set(dm, j = "gender", value = factor(dm$gender))
impute_archive <- fread(paths[["original_imputation"]])
for (v in covars) {
  source_row <- impute_archive[variable == v]
  new_row <- imp$audit[variable == v]
  stopifnot(nrow(source_row) == 1L, nrow(new_row) == 1L,
            source_row$missing_n == new_row$missing_n,
            as.character(source_row$fill_value) == as.character(new_row$fill_value))
}

form <- reformulate(covars, response = "rrt_first_24_72")
stopifnot(!any(c("cluster_k2", "ablation_group", "sofa_score", "mortality_30d") %in%
                 all.vars(form)))
fit_warnings <- character()
fit <- withCallingHandlers(glm(form, data = dm, family = binomial()),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(isTRUE(fit$converged), all(is.finite(coef(fit))), nobs(fit) == nrow(dm))
raw_e <- as.numeric(predict(fit, type = "response"))
eps <- 1e-6
e <- pmin(pmax(raw_e, eps), 1 - eps)
A <- as.integer(dm$rrt_first_24_72)
dm[, `:=`(common_e = e,
         w_ATE = fifelse(A == 1L, 1 / e, 1 / (1 - e)),
         w_ATO = fifelse(A == 1L, 1 - e, e),
         gender_male = as.numeric(gender == "M"))]
stopifnot(all(is.finite(dm$w_ATE) & dm$w_ATE > 0),
          all(is.finite(dm$w_ATO) & dm$w_ATO > 0),
          all(is.finite(dm$common_e) & dm$common_e > 0 & dm$common_e < 1))

qc <- data.table(
  formula = paste(deparse(form), collapse = ""),
  n_input = nrow(risk), n_model = nobs(fit), record_positive_n = sum(A),
  converged = fit$converged,
  coefficient_count = length(coef(fit)),
  nonfinite_coefficient_count = sum(!is.finite(coef(fit))),
  warning_count = length(fit_warnings),
  warning_text = paste(unique(fit_warnings), collapse = " | "),
  probability_min_raw = min(raw_e), probability_max_raw = max(raw_e),
  clipped_probability_n = sum(raw_e < eps | raw_e > 1 - eps),
  record_event_discrimination_auc = d4_binary_auc(A, e),
  partition_or_sofa_in_formula = any(c("cluster_k2", "ablation_group", "sofa_score") %in%
                                       all.vars(form)),
  effect_model_fitted = FALSE
)
coefficients <- data.table(term = names(coef(fit)), coefficient = unname(coef(fit)))

ess <- function(w) {
  if (!length(w) || any(!is.finite(w) | w <= 0) || sum(w) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}
wmean <- function(x, w) sum(x * w) / sum(w)
balance_vars <- c(setdiff(covars, "gender"), "gender_male")
stopifnot(length(balance_vars) == 15L,
          all(vapply(dm[, ..balance_vars], function(x) all(is.finite(x)), logical(1))))

group_rows <- list()
balance_rows <- list()
for (partition_name in c("K2", "A4")) {
  group_col <- if (partition_name == "K2") "cluster_k2" else "ablation_group"
  for (g in 1:2) {
    z <- dm[get(group_col) == g]
    pos <- z[rrt_first_24_72 == 1L]
    neg <- z[rrt_first_24_72 == 0L]
    n1 <- nrow(pos); n0 <- nrow(neg)
    if (n1 == 0L || n0 == 0L) stop("Empty record-status cell")
    low <- max(min(pos$common_e), min(neg$common_e))
    high <- min(max(pos$common_e), max(neg$common_e))
    stopifnot(low <= high)
    group_name <- paste0(if (partition_name == "K2") "C" else "A", g)
    row <- data.table(
      partition = partition_name, group = group_name,
      n = nrow(z), n_record_positive = n1, n_record_negative = n0,
      record_positive_rate = n1 / nrow(z),
      original_one_percent_rate_criterion = n1 / nrow(z) >= 0.01,
      positive_deaths_30d = sum(pos$mortality_30d == 1L),
      negative_deaths_30d = sum(neg$mortality_30d == 1L),
      positive_outcome_unknown = sum(is.na(pos$mortality_30d)),
      negative_outcome_unknown = sum(is.na(neg$mortality_30d)),
      overlap_low = low, overlap_high = high,
      positive_outside_n = sum(pos$common_e < low | pos$common_e > high),
      negative_outside_n = sum(neg$common_e < low | neg$common_e > high)
    )
    for (target in c("ATE", "ATO")) {
      weight <- paste0("w_", target)
      w1 <- pos[[weight]]; w0 <- neg[[weight]]
      top_n <- ceiling(n1 / 100)
      set(row, j = paste0(target, "_positive_ESS"), value = ess(w1))
      set(row, j = paste0(target, "_negative_ESS"), value = ess(w0))
      set(row, j = paste0(target, "_positive_ESS_per_n"), value = ess(w1) / n1)
      set(row, j = paste0(target, "_negative_ESS_per_n"), value = ess(w0) / n0)
      set(row, j = paste0(target, "_positive_max_normalized_weight"),
          value = max(w1 / sum(w1)))
      set(row, j = paste0(target, "_negative_max_normalized_weight"),
          value = max(w0 / sum(w0)))
      set(row, j = paste0(target, "_positive_top_one_percent_n"), value = top_n)
      set(row, j = paste0(target, "_positive_top_one_percent_weight_share"),
          value = sum(sort(w1 / sum(w1), decreasing = TRUE)[seq_len(top_n)]))
      stopifnot(ess(w1) <= n1 + 1e-8, ess(w0) <= n0 + 1e-8)
    }
    group_rows[[length(group_rows) + 1L]] <- row

    for (v in balance_vars) {
      x1 <- pos[[v]]; x0 <- neg[[v]]
      reference_sd <- sqrt((var(x1) + var(x0)) / 2)
      if (!is.finite(reference_sd) || reference_sd <= 0) {
        stop("Undefined balance reference SD: ", partition_name, "/", group_name, "/", v)
      }
      smd <- function(w1, w0) (wmean(x1, w1) - wmean(x0, w0)) / reference_sd
      balance_rows[[length(balance_rows) + 1L]] <- data.table(
        partition = partition_name, group = group_name, variable = v,
        unweighted_reference_sd = reference_sd,
        SMD_unweighted = smd(rep(1, n1), rep(1, n0)),
        SMD_ATE = smd(pos$w_ATE, neg$w_ATE),
        SMD_ATO = smd(pos$w_ATO, neg$w_ATO),
        interpretable = v != "urine_output_24h_ml"
      )
    }
  }
}
groups <- rbindlist(group_rows)
balance <- rbindlist(balance_rows)
stopifnot(nrow(groups) == 4L, nrow(balance) == 60L,
          all(is.finite(balance$SMD_unweighted)),
          all(is.finite(balance$SMD_ATE)), all(is.finite(balance$SMD_ATO)))

screen <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
  r <- groups[i]
  b <- balance[partition == r$partition & group == r$group & interpretable]
  stopifnot(nrow(b) == 14L)
  failure <- b[abs(SMD_ATO) >= 0.1]
  data.table(
    partition = r$partition, group = r$group,
    count_SMD_ATE_above_point_one = b[abs(SMD_ATE) > 0.1, .N],
    count_SMD_ATO_at_or_above_point_one = nrow(failure),
    variables_failing_ATO = paste(sprintf("%s:%.6f", failure$variable,
                                           failure$SMD_ATO), collapse = "; "),
    death_cell_minimum_met = r$positive_deaths_30d >= 10L &&
      r$negative_deaths_30d >= 10L,
    measured_balance_and_event_screen_met = nrow(failure) == 0L &&
      r$positive_deaths_30d >= 10L && r$negative_deaths_30d >= 10L
  )
}))
between <- screen[, .(both_groups_meet_screen = all(measured_balance_and_event_screen_met)),
                  by = partition]
between[, historical_A4_prediction_both_groups_meet :=
          ifelse(partition == "A4", TRUE, NA)]
stopifnot(nrow(between) == 2L, setequal(between$partition, c("K2", "A4")))

v1_balance <- fread(paths[["v1_balance"]])
v1_groups <- fread(paths[["v1_groups"]])
v1_screen <- fread(paths[["v1_screen"]])
compare <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
  r <- groups[i]
  old_b <- v1_balance[partition == r$partition & group == r$group &
                        interpretable == TRUE & variable != "sofa_score"]
  new_b <- balance[partition == r$partition & group == r$group & interpretable]
  old_r <- v1_groups[partition == r$partition & group == r$group]
  old_screen <- v1_screen[partition == r$partition & group == r$group]
  stopifnot(nrow(old_b) == 14L, nrow(new_b) == 14L,
            setequal(old_b$variable, new_b$variable),
            nrow(old_r) == 1L, nrow(old_screen) == 1L,
            old_r$n == r$n, old_r$n_record_positive == r$n_record_positive)
  data.table(
    partition = r$partition, group = r$group, n = r$n,
    record_positive_n = r$n_record_positive,
    v1_ATE_positive_ESS = old_r$ATE_positive_ESS,
    v2_ATE_positive_ESS = r$ATE_positive_ESS,
    v1_ATO_positive_ESS = old_r$ATO_positive_ESS,
    v2_ATO_positive_ESS = r$ATO_positive_ESS,
    v1_common14_ATO_SMD_fail_n = sum(abs(old_b$SMD_ATO) >= 0.1),
    v2_common14_ATO_SMD_fail_n = sum(abs(new_b$SMD_ATO) >= 0.1),
    v1_original15_ATO_SMD_fail_n = old_screen$count_SMD_ATO_at_or_above_point_one
  )
}))

after_sha <- vapply(paths, sha, character(1))
stopifnot(identical(before_sha, after_sha))
manifest <- data.table(role = names(paths), path = unname(paths),
                       sha256_before = unname(before_sha),
                       sha256_after = unname(after_sha), unchanged = TRUE)

saveRDS(dm[, .(stay_id, cluster_k2, ablation_group, rrt_first_24_72,
               mortality_30d, common_e, w_ATE, w_ATO)],
        file.path(private, "no_sofa_common_score_private.rds"))
fwrite(qc, file.path(out, "no_sofa_model_QC.csv"))
fwrite(coefficients, file.path(out, "no_sofa_model_coefficients.csv"))
fwrite(imp$audit, file.path(out, "no_sofa_imputation.csv"))
fwrite(groups, file.path(out, "no_sofa_group_support.csv"))
fwrite(balance, file.path(out, "no_sofa_covariate_balance.csv"))
fwrite(screen, file.path(out, "no_sofa_group_screen.csv"))
fwrite(between, file.path(out, "no_sofa_between_group_screen.csv"))
fwrite(compare, file.path(out, "v1_v2_common14_comparison.csv"))
fwrite(manifest, file.path(out, "source_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "R_session.txt"))
cat("No-SOFA common record-event model fitted once; A4 screen=",
    between[partition == "A4", both_groups_meet_screen], "\n", sep = "")
