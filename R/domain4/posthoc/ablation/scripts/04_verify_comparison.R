options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
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
risk_path <- file.path(d4, "output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds")
original_path <- file.path(d4, "output_formal_v1/private_not_for_release/D4_propensity_patient_audit_private.rds")
stopifnot(sha(risk_path) == "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
          sha(original_path) == "BC6106F70DD96457D4D3BA11915A89714EC2C8C581ABDC5F7DAF0C95FE72DC44")

risk <- as.data.table(readRDS(risk_path))
score <- as.data.table(readRDS(file.path(private, "common_propensity_patient_audit_private.rds")))
aggregate <- fread(file.path(out, "parallel_group_support.csv"))
balance <- fread(file.path(out, "parallel_covariate_balance.csv"))
screen <- fread(file.path(out, "posthoc_group_screen.csv"))
overall <- fread(file.path(out, "posthoc_between_group_screen.csv"))
qc <- fread(file.path(out, "common_propensity_model_QC.csv"))

stopifnot(nrow(risk) == 17525L, nrow(score) == 17525L,
          uniqueN(score$stay_id) == nrow(score), uniqueN(risk$stay_id) == nrow(risk),
          nrow(aggregate) == 4L, nrow(balance) == 64L, nrow(screen) == 4L,
          nrow(overall) == 2L, nrow(qc) == 1L,
          identical(as.logical(overall[partition == "A4", predicted_A4_both_groups_meet]), TRUE),
          identical(as.logical(overall[partition == "A4", both_groups_meet_screen]), FALSE),
          !qc$any_partition_in_formula, !qc$effect_model_fitted,
          qc$converged, qc$nonfinite_coefficient_count == 0L,
          qc$n_model == nrow(score), qc$n_record_positive == 314L)
setkey(risk, stay_id)
setkey(score, stay_id)
stopifnot(identical(risk$stay_id, score$stay_id),
          identical(as.integer(as.character(risk$cluster_k2)), score$cluster_k2),
          identical(as.integer(risk$rrt_first_24_72), score$rrt_first_24_72),
          identical(as.integer(risk$mortality_30d), score$mortality_30d))

close_enough <- function(actual, expected, label, tol = 1e-9) {
  if (length(actual) != 1L || length(expected) != 1L || !is.finite(actual) ||
      !is.finite(expected) || abs(actual - expected) > tol) {
    stop(sprintf("Verification failed for %s: actual=%s expected=%s",
                 label, actual, expected))
  }
  TRUE
}

e <- score$ps_common
A <- score$rrt_first_24_72
stopifnot(all(e > 0 & e < 1), all(is.finite(e)),
          max(abs(score$w_ATE_common - ifelse(A == 1L, 1 / e, 1 / (1 - e)))) < 1e-9,
          max(abs(score$w_ATO_common - ifelse(A == 1L, 1 - e, e))) < 1e-12)

calc_ess <- function(w) sum(w)^2 / sum(w^2)
checks <- list()
for (partition in c("K2", "A4")) {
  p <- if (partition == "K2") score$cluster_k2 else score$ablation_group
  for (level in 1:2) {
    group <- paste0(if (partition == "K2") "C" else "A", level)
    partition_name <- partition
    group_name <- group
    i <- which(p == level)
    z <- score[i]
    ref <- aggregate[partition == partition_name & group == group_name]
    verdict <- screen[partition == partition_name & group == group_name]
    stopifnot(nrow(ref) == 1L, nrow(verdict) == 1L)
    n1 <- sum(z$rrt_first_24_72 == 1L)
    n0 <- sum(z$rrt_first_24_72 == 0L)
    stopifnot(nrow(z) == ref$n, n1 == ref$n_record_positive,
              n0 == ref$n_record_negative,
              sum(z$mortality_30d[z$rrt_first_24_72 == 1L]) == ref$positive_deaths_30d,
              sum(z$mortality_30d[z$rrt_first_24_72 == 0L]) == ref$negative_deaths_30d)
    close_enough(n1 / nrow(z), ref$record_positive_rate, paste(group, "rate"))
    pos <- z[rrt_first_24_72 == 1L]
    neg <- z[rrt_first_24_72 == 0L]
    low <- max(min(pos$ps_common), min(neg$ps_common))
    high <- min(max(pos$ps_common), max(neg$ps_common))
    stopifnot(sum(pos$ps_common < low | pos$ps_common > high) == ref$positive_outside_n,
              sum(neg$ps_common < low | neg$ps_common > high) == ref$negative_outside_n)
    for (target in c("ATE", "ATO")) {
      weight <- paste0("w_", target, "_common")
      expected_ess1 <- calc_ess(pos[[weight]])
      expected_ess0 <- calc_ess(neg[[weight]])
      close_enough(expected_ess1, ref[[paste0(target, "_positive_ESS")]],
                   paste(group, target, "positive ESS"))
      close_enough(expected_ess0, ref[[paste0(target, "_negative_ESS")]],
                   paste(group, target, "negative ESS"))
      close_enough(expected_ess1 / n1, ref[[paste0(target, "_positive_ESS_per_n")]],
                   paste(group, target, "positive ESS fraction"))
      share <- sort(pos[[weight]] / sum(pos[[weight]]), decreasing = TRUE)
      top_n <- ceiling(n1 / 100)
      close_enough(sum(share[seq_len(top_n)]),
                   ref[[paste0(target, "_top_one_percent_positive_weight_share")]],
                   paste(group, target, "top share"))
    }
    b <- balance[partition == partition_name & group == group_name & interpretable == TRUE]
    stopifnot(nrow(b) == 15L,
              sum(abs(b$SMD_ATO) >= 0.1) == verdict$count_SMD_ATO_at_or_above_point_one,
              identical(as.logical(verdict$measured_balance_and_event_screen_met),
                        all(abs(b$SMD_ATO) < 0.1) &&
                          ref$positive_deaths_30d >= 10L && ref$negative_deaths_30d >= 10L))
    checks[[length(checks) + 1L]] <- data.table(
      partition = partition, group = group, n = nrow(z), positive_n = n1,
      positive_deaths = ref$positive_deaths_30d,
      ATE_positive_ESS = calc_ess(pos$w_ATE_common),
      ATO_positive_ESS = calc_ess(pos$w_ATO_common),
      ATO_failed_covariates = sum(abs(b$SMD_ATO) >= 0.1),
      all_checks_passed = TRUE
    )
  }
}

source(file.path(d4, "scripts/38_d4_landmark_support_common_v1.R"), local = TRUE)
covars <- setdiff(D4_PS_COVARIATES, "cluster_k2")
imp <- d4_impute_data(risk, covars, "median_mode")
check_data <- as.data.table(imp$data)
check_data[, gender_male := as.numeric(as.character(gender) == "M")]
for (partition in c("K2", "A4")) {
  p <- if (partition == "K2") score$cluster_k2 else score$ablation_group
  for (level in 1:2) {
    group <- paste0(if (partition == "K2") "C" else "A", level)
    partition_name <- partition
    group_name <- group
    rows <- which(p == level)
    treated <- rows[A[rows] == 1L]
    untreated <- rows[A[rows] == 0L]
    for (v in c(setdiff(covars, c("gender", "urine_output_24h_ml")), "gender_male")) {
      variable_name <- v
      x1 <- check_data[[v]][treated]
      x0 <- check_data[[v]][untreated]
      sd_ref <- sqrt((var(x1) + var(x0)) / 2)
      expected <- ((sum(x1 * score$w_ATO_common[treated]) /
                      sum(score$w_ATO_common[treated])) -
                     (sum(x0 * score$w_ATO_common[untreated]) /
                      sum(score$w_ATO_common[untreated]))) / sd_ref
      observed <- balance[partition == partition_name & group == group_name & variable == variable_name,
                          SMD_ATO]
      close_enough(expected, observed, paste(partition, group, v, "ATO SMD"))
    }
  }
}

report <- rbindlist(checks)
fwrite(report, file.path(out, "independent_numeric_verification.csv"))
cat("PASS: 4 group rows, both weight families, 60 interpretable SMDs, and source joins.\n")
