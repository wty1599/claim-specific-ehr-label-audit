options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
parent <- dirname(root)
project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
src <- file.path(project, "d4_landmark_support_interface_20260729/output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds")
score_path <- file.path(root, "private_not_for_release/no_sofa_common_score_private.rds")
out <- file.path(root, "outputs")
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
stopifnot(sha(file.path(root, "01_REVISION_SPEC.md")) ==
            "5B51EC3B29B0D9E503099145E88D376C933283BF56AE75AC6FDE4B0757FF7984",
          sha(src) == "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
          sha(file.path(parent, "private_not_for_release/A4_first_seed_labels_private.rds")) ==
            "0A0E9AFAB6D320B48CBD5C7A13B814A46C5C3F8A84DF47DC3A69763A13795987",
          sha(file.path(parent, "private_not_for_release/common_propensity_patient_audit_private.rds")) ==
            "CA2B096AE452776BFD243F2116A2036A64F92ACB9BC5BEDCA72BA7B70E459887")

risk <- as.data.table(readRDS(src))
patient <- as.data.table(readRDS(score_path))
groups <- fread(file.path(out, "no_sofa_group_support.csv"))
balance <- fread(file.path(out, "no_sofa_covariate_balance.csv"))
screen <- fread(file.path(out, "no_sofa_group_screen.csv"))
comparison <- fread(file.path(out, "v1_v2_common14_comparison.csv"))
qc <- fread(file.path(out, "no_sofa_model_QC.csv"))
coefficients <- fread(file.path(out, "no_sofa_model_coefficients.csv"))
imputation <- fread(file.path(out, "no_sofa_imputation.csv"))

stopifnot(nrow(risk) == 17525L, nrow(patient) == nrow(risk),
          uniqueN(risk$stay_id) == nrow(risk), uniqueN(patient$stay_id) == nrow(patient),
          !anyNA(match(risk$stay_id, patient$stay_id)),
          nrow(groups) == 4L, nrow(balance) == 60L, nrow(screen) == 4L,
          nrow(comparison) == 4L, nrow(qc) == 1L,
          !any(grepl("sofa|cluster_k2|ablation_group", coefficients$term, ignore.case = TRUE)),
          !any(balance$variable == "sofa_score"),
          isTRUE(qc$converged), qc$warning_count == 0L,
          qc$partition_or_sofa_in_formula == FALSE,
          qc$effect_model_fitted == FALSE)
idx <- match(risk$stay_id, patient$stay_id)
patient <- patient[idx]
stopifnot(identical(as.numeric(risk$stay_id), as.numeric(patient$stay_id)),
          identical(as.integer(risk$rrt_first_24_72), as.integer(patient$rrt_first_24_72)),
          identical(as.integer(risk$mortality_30d), as.integer(patient$mortality_30d)),
          all(as.integer(as.character(risk$cluster_k2)) == patient$cluster_k2))

e <- patient$common_e
a <- patient$rrt_first_24_72
ate_check <- ifelse(a == 1L, 1 / e, 1 / (1 - e))
ato_check <- ifelse(a == 1L, 1 - e, e)
stopifnot(max(abs(ate_check - patient$w_ATE)) < 1e-9,
          max(abs(ato_check - patient$w_ATO)) < 1e-12,
          all(e > 0 & e < 1), sum(a) == 314L)

model_covars <- c("age", "gender", "aki_stage_0_24h", "sapsii",
                  "creatinine_max", "bun_max", "potassium_max",
                  "bicarbonate_min", "ph_min", "lactate_max",
                  "urine_output_24h_ml", "mbp_min", "heart_rate_max",
                  "resp_rate_max", "gcs_min")
for (v in model_covars) {
  fill_row <- imputation[variable == v]
  stopifnot(nrow(fill_row) == 1L)
  x <- risk[[v]]
  missing <- is.na(x) | (is.character(x) & !nzchar(x))
  stopifnot(sum(missing) == fill_row$missing_n)
  if (v == "gender") {
    x <- as.character(x)
    x[missing] <- as.character(fill_row$fill_value)
  } else {
    x[missing] <- as.numeric(fill_row$fill_value)
  }
  risk[[v]] <- x
}
risk[, gender_male := as.numeric(gender == "M")]
balance_vars <- c(setdiff(model_covars, "gender"), "gender_male")
stopifnot(length(balance_vars) == 15L)

max_group_error <- 0
max_balance_error <- 0
verification <- list()
for (i in seq_len(nrow(groups))) {
  r <- groups[i]
  s <- screen[partition == r$partition & group == r$group]
  col <- if (r$partition == "K2") "cluster_k2" else "ablation_group"
  keep <- patient[[col]] == as.integer(substring(r$group, 2L))
  z <- patient[keep]
  x <- risk[keep]
  positive <- z$rrt_first_24_72 == 1L
  stopifnot(nrow(z) == r$n, sum(positive) == r$n_record_positive,
            sum(!positive) == r$n_record_negative,
            sum(z$mortality_30d[positive] == 1L) == r$positive_deaths_30d,
            sum(z$mortality_30d[!positive] == 1L) == r$negative_deaths_30d)
  low <- max(min(z$common_e[positive]), min(z$common_e[!positive]))
  high <- min(max(z$common_e[positive]), max(z$common_e[!positive]))
  outside_pos <- sum(z$common_e[positive] < low | z$common_e[positive] > high)
  outside_neg <- sum(z$common_e[!positive] < low | z$common_e[!positive] > high)
  stopifnot(abs(low - r$overlap_low) < 1e-12,
            abs(high - r$overlap_high) < 1e-12,
            outside_pos == r$positive_outside_n,
            outside_neg == r$negative_outside_n)
  for (target in c("ATE", "ATO")) {
    w <- z[[paste0("w_", target)]]
    for (arm in c("positive", "negative")) {
      ww <- w[if (arm == "positive") positive else !positive]
      n <- length(ww)
      ess <- sum(ww)^2 / sum(ww^2)
      ratio <- ess / n
      max_share <- max(ww / sum(ww))
      err <- max(abs(c(ess - r[[paste0(target, "_", arm, "_ESS")]],
                       ratio - r[[paste0(target, "_", arm, "_ESS_per_n")]],
                       max_share - r[[paste0(target, "_", arm, "_max_normalized_weight")]])))
      max_group_error <- max(max_group_error, err)
      if (arm == "positive") {
        top_n <- ceiling(n / 100)
        share <- sum(sort(ww / sum(ww), decreasing = TRUE)[seq_len(top_n)])
        max_group_error <- max(max_group_error,
                               abs(share - r[[paste0(target, "_positive_top_one_percent_weight_share")]]))
      }
    }
  }
  failures <- character()
  for (v in balance_vars) {
    b <- balance[partition == r$partition & group == r$group & variable == v]
    stopifnot(nrow(b) == 1L)
    xx <- x[[v]]
    xx1 <- xx[positive]; xx0 <- xx[!positive]
    sd_fixed <- sqrt((stats::var(xx1) + stats::var(xx0)) / 2)
    smd <- function(w) {
      weighted.mean(xx1, w[positive]) / sd_fixed -
        weighted.mean(xx0, w[!positive]) / sd_fixed
    }
    calc <- c(reference = sd_fixed,
              unweighted = smd(rep(1, nrow(z))),
              ATE = smd(z$w_ATE), ATO = smd(z$w_ATO))
    reported <- c(reference = b$unweighted_reference_sd,
                  unweighted = b$SMD_unweighted,
                  ATE = b$SMD_ATE, ATO = b$SMD_ATO)
    max_balance_error <- max(max_balance_error, max(abs(calc - reported)))
    if (v != "urine_output_24h_ml" && abs(calc[["ATO"]]) >= 0.1) {
      failures <- c(failures, v)
    }
  }
  stopifnot(length(failures) == s$count_SMD_ATO_at_or_above_point_one,
            (length(failures) == 0L && min(r$positive_deaths_30d,
                                           r$negative_deaths_30d) >= 10L) ==
              s$measured_balance_and_event_screen_met)
  verification[[i]] <- data.table(partition = r$partition, group = r$group,
                                  n = nrow(z), record_positive_n = sum(positive),
                                  recalculated_ATO_fail_n = length(failures),
                                  result_matched = TRUE)
}
stopifnot(max_group_error < 1e-7, max_balance_error < 1e-9)

v1_hashes <- c(
  balance = "923C6D2511F92750DC69CF3B62E2C72337AB40F13AC37C7B0FB4D2BC10893B59",
  groups = "164113E5F7AE7EA617E6B64D04399A445A4F2FF03CA417EFDFE5049C4E809954",
  screen = "2A770AAE55F0A7363BC75F1FC50184267FD0188FEB4C27AB639D3D84C94D1602",
  between = "1D8B55AFBB810523A66C9E33522C8D4556B1790DFF2EABFE31560AC763498343"
)
v1_paths <- file.path(parent, "outputs",
                      c("parallel_covariate_balance.csv", "parallel_group_support.csv",
                        "posthoc_group_screen.csv", "posthoc_between_group_screen.csv"))
stopifnot(all(vapply(v1_paths, sha, character(1)) == unname(v1_hashes)))

fwrite(rbindlist(verification), file.path(out, "independent_numeric_verification.csv"))
writeLines(c(
  paste0("Verification passed for all 4 partition groups and 60 balance rows."),
  paste0("Maximum absolute aggregate-group error: ", format(max_group_error, digits = 16)),
  paste0("Maximum absolute balance error: ", format(max_balance_error, digits = 16)),
  "Original risk ledger and v1 patient-level and aggregate objects retain the locked hashes.",
  "No effect model was fitted in verification."
), file.path(out, "independent_verification_summary.txt"))
cat("Independent no-SOFA numeric verification passed\n")
