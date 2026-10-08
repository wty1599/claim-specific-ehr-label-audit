options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/source_correction/00_correction_paths.R"))
out <- file.path(aggregate, "missing_pattern_transfer")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                                    serialize = FALSE))

spec <- file.path(spec_root, "D3_MISSING_PATTERN_TRANSFER_SPEC_20261008.md")
stopifnot(identical(sha(spec),
  "784593E379B1AD81BB9013D4A19834BF2F6C21852DFBBB69BE78D78EE4E7F1FA"))
param_path <- file.path(recovery_root,
                        "outputs/original_parameters_aggregate_v1.rds")
stopifnot(identical(sha(param_path),
  "70CEAD3058658339D3AE00753F7B21045F0BBCCB74B8E9267FC2F5B906931C99"))
e_path <- file.path(base, "eicu_external.csv")
t_path <- file.path(run_root, "private_local/eicu_temperature_corrected_v1.csv")
m_path <- file.path(base, "final_full.csv")
l_path <- file.path(source_root, "private_local/mimic_L0_L1_L2_labels_local.csv")
stopifnot(all(file.exists(c(e_path, t_path, m_path, l_path))))
stopifnot(identical(sha(e_path),
  "249F28CAC3C900B2F8B187D8CC3F2887230DB7B6BFE54EAB19717377F223B5A5"))
stopifnot(identical(sha(t_path),
  "F45422CC054DAD2493309AE8E225362A713E80D33EC37DE90773B4BA9A3FF827"))
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))

params <- readRDS(param_path)
features <- params$feature_names
stopifnot(length(features) == 33L, !anyDuplicated(features))
recipe <- data.table(
  feature = features,
  winsor_p01 = params$winsor_limits[, 1L],
  winsor_p99 = params$winsor_limits[, 2L],
  imputation_median =
    params$inductive_median_from_pre_imputation_winsorized_observed,
  reference_mean = unname(params$original_center),
  reference_sd = unname(params$original_scale)
)
stopifnot(identical(as.character(recipe$feature), features))
stopifnot(all(recipe$imputation_median >= recipe$winsor_p01 &
              recipe$imputation_median <= recipe$winsor_p99))

M <- fread(m_path, select = c("stay_id", features), showProgress = FALSE)
L <- fread(l_path, select = c("stay_id", "locked_k2", "missing_feature_n"),
           showProgress = FALSE)
stopifnot(nrow(M) == 20049L, nrow(L) == nrow(M),
          !anyDuplicated(M$stay_id), !anyDuplicated(L$stay_id))
L <- L[match(M$stay_id, L$stay_id)]
stopifnot(identical(M$stay_id, L$stay_id))
Mprep <- apply_frozen_recipe(M, recipe, features)
stopifnot(identical(as.integer(Mprep$missing_n),
                    as.integer(L$missing_feature_n)))
complete_i <- which(Mprep$missing_n == 0L)
stopifnot(length(complete_i) == 4689L)
X0 <- Mprep$X[complete_i, , drop = FALSE]
locked <- as.integer(L$locked_k2[complete_i])
baseline <- as.integer(assign_nearest_centroid(
  X0, params$original_centroids)$assigned_cluster)
stopifnot(identical(baseline, locked))
baseline_c1 <- mean(locked == 1L)

E <- fread(e_path, select = c("patientunitstayid", features),
           showProgress = FALSE)
T <- fread(t_path, showProgress = FALSE)
stopifnot(nrow(E) == 17465L, nrow(T) == nrow(E),
          !anyDuplicated(E$patientunitstayid),
          !anyDuplicated(T$patientunitstayid))
T <- T[match(E$patientunitstayid, T$patientunitstayid)]
stopifnot(identical(E$patientunitstayid, T$patientunitstayid))
E[, `:=`(temperature_min = T$temperature_min_corrected,
         temperature_max = T$temperature_max_corrected)]
external_mask <- !is.finite(as.matrix(E[, ..features]))
stopifnot(nrow(external_mask) == 17465L,
          ncol(external_mask) == 33L)
external_missing <- rowSums(external_mask)
audit_labels <- fread(file.path(run_root, "private_local",
                                "eicu_corrected_L1_labels_local.csv"),
                      select = c("patientunitstayid", "corrected_missing_n"))
audit_labels <- audit_labels[match(E$patientunitstayid,
                                   audit_labels$patientunitstayid)]
stopifnot(identical(E$patientunitstayid, audit_labels$patientunitstayid),
          identical(as.integer(external_missing),
                    as.integer(audit_labels$corrected_missing_n)))
observed_external <- fread(file.path(run_root, "aggregate_outputs",
                                     "corrected_assignment_by_missingness.csv"))
stopifnot(sum(observed_external$stays) == 17465L)

band <- function(n) {
  fifelse(n == 0L, "0", fifelse(n <= 2L, "1-2",
          fifelse(n <= 5L, "3-5", "6+")))
}
stand_median <- (recipe$imputation_median - recipe$reference_mean) /
  recipe$reference_sd
n <- length(complete_i)
reps <- 500L
set.seed(20261008L)
results <- vector("list", reps * 5L)
row <- 0L
for (rep in seq_len(reps)) {
  idx <- sample.int(nrow(external_mask), n, replace = TRUE)
  mask <- external_mask[idx, , drop = FALSE]
  X <- X0
  for (j in seq_along(features)) {
    if (any(mask[, j])) X[mask[, j], j] <- stand_median[j]
  }
  assigned <- as.integer(assign_nearest_centroid(
    X, params$original_centroids)$assigned_cluster)
  stopifnot(length(assigned) == n, all(assigned %in% 1:2))
  missing_n <- external_missing[idx]
  groups <- band(missing_n)
  for (stratum in c("all", "0", "1-2", "3-5", "6+")) {
    subset <- if (stratum == "all") rep(TRUE, n) else groups == stratum
    row <- row + 1L
    results[[row]] <- data.table(
      repetition = rep, stratum = stratum, n = sum(subset),
      baseline_c1 = mean(locked[subset] == 1L),
      induced_c1 = mean(assigned[subset] == 1L),
      c1_to_c2 = sum(subset & locked == 1L & assigned == 2L),
      c2_to_c1 = sum(subset & locked == 2L & assigned == 1L),
      flips = sum(subset & locked != assigned)
    )
  }
}
R <- rbindlist(results)
stopifnot(nrow(R) == reps * 5L, all(R$n > 0L),
          all(R[ stratum == "all", n] == n),
          all(R$flips == R$c1_to_c2 + R$c2_to_c1))
R[, `:=`(flip_rate = flips / n,
         c1_to_c2_rate = c1_to_c2 / n,
         c2_to_c1_rate = c2_to_c1 / n,
         c1_change_pp = 100 * (induced_c1 - baseline_c1))]

metric_cols <- c("n", "baseline_c1", "induced_c1", "flip_rate",
                 "c1_to_c2_rate", "c2_to_c1_rate", "c1_change_pp")
summary <- rbindlist(lapply(metric_cols, function(metric) {
  R[, .(metric = metric, mean = mean(get(metric)),
        mcse = sd(get(metric)) / sqrt(.N),
        p025 = as.numeric(quantile(get(metric), .025)),
        p975 = as.numeric(quantile(get(metric), .975))),
    by = stratum]
}))
setorder(summary, stratum, metric)
fwrite(summary, file.path(out, "transfer_metric_summary.csv"))
fwrite(R, file.path(out, "transfer_repetition_aggregates.csv"))
fwrite(data.table(
  complete_mimic_n = n,
  complete_mimic_c1_n = sum(locked == 1L),
  complete_mimic_c1_prevalence = baseline_c1,
  eicu_mask_pool_n = nrow(external_mask),
  repetitions = reps,
  baseline_labels_reproduced = sum(baseline == locked)),
  file.path(out, "transfer_input_and_baseline_checks.csv"))
fwrite(observed_external,
       file.path(out, "observed_corrected_eicu_by_missingness.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(summary[stratum == "all"])
cat("MISSING_PATTERN_TRANSFER_COMPLETE\n")
