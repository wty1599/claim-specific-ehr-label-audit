options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
base <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT")
stopifnot(nzchar(output_root), nzchar(base), nzchar(repo))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(output_root, repo)
run_root <- file.path(output_root, "domain3_L1")
recovery_root <- file.path(output_root, "domain3_parameter_recovery")
aggregate_out <- file.path(run_root, "aggregate_outputs")
private_out <- file.path(run_root, "private_local")
dir.create(aggregate_out, recursive = TRUE, showWarnings = FALSE)
dir.create(private_out, recursive = TRUE, showWarnings = FALSE)

sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                               serialize = FALSE))
lock_path <- file.path(repo, "R/domain3/deployed_rule/L1_REANALYSIS_LOCK.md")
stopifnot(sha(lock_path) ==
  "3BC5D702ED013192FCB3657176B112FFB3A06F88A9509DE987242A6B355545F9")
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))
parameter_path <- file.path(recovery_root, "outputs",
                            "original_parameters_aggregate_v1.rds")
stopifnot(file.exists(parameter_path))
params <- readRDS(parameter_path)
stopifnot(identical(params$object_type,
                    "recovered_original_centroids_with_inductive_median_adapter"),
          identical(params$verification$fixed_label_matches, 20049L),
          identical(unname(params$source_sha256[c("matrix", "labels", "feature_sets", "mice")]),
                    c("9896AD4906207CA910E54EB187507037106BFF212299586141908F387CFE6279",
                      "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
                      "8BC1AABF3ACB88BD6C37F8A857A308ABE92955C65264249D7F31F4827F73C8DF",
                      "B9EF638A9BBD3AA6672D15F788187A70765197B11279A2C1A4A721B0BE87EA80")))
features <- params$feature_names
stopifnot(length(features) == 33L, !anyDuplicated(features))

input_paths <- c(
  mimic = file.path(base, "final_full.csv"),
  eicu = file.path(base, "eicu_external.csv"),
  locked_labels = file.path(output_root, "model/labels_primary_mice.rds"),
  l2_object = file.path(output_root, "domain3_deployable_generator_v1_20260728",
                        "model/ehr_audit_k2_deployable_surrogate_generator_v1.rds"),
  l2_eicu = file.path(output_root, "domain3_deployable_generator_v1_20260728",
                      "labels/eicu_d3_generator_v1_labels.csv")
)
expected_md5 <- c(
  mimic = "a3c6582cc68922787b432632f3aff8d7",
  eicu = "daa5589d74bdf7c09394863cec456280",
  locked_labels = "e741e9cf6bc39bf8075127a76dab5b2a",
  l2_object = "4722d24103151ccc5bcbcc737772f249",
  l2_eicu = "5958bcba237a61ed2689e590a572e3ca"
)
stopifnot(all(file.exists(input_paths)))
actual_md5 <- tolower(unname(tools::md5sum(input_paths)))
names(actual_md5) <- names(input_paths)
stopifnot(identical(actual_md5, expected_md5))
old <- readRDS(input_paths[["l2_object"]])
stopifnot(identical(as.character(old$feature_names), features))

recipe <- data.table(
  feature = features,
  winsor_p01 = params$winsor_limits[, 1L],
  winsor_p99 = params$winsor_limits[, 2L],
  imputation_median = params$inductive_median_from_pre_imputation_winsorized_observed,
  reference_mean = unname(params$original_center),
  reference_sd = unname(params$original_scale)
)
stopifnot(identical(as.character(recipe$feature), features))
missing_band <- function(x) {
  factor(fifelse(x == 0L, "0",
          fifelse(x <= 2L, "1-2", fifelse(x <= 5L, "3-5", "6+"))),
         levels = c("0", "1-2", "3-5", "6+"))
}

# MIMIC: no outcome column is loaded during assignment.
M <- fread(input_paths[["mimic"]], select = c("stay_id", features),
           showProgress = FALSE)
L <- as.data.table(readRDS(input_paths[["locked_labels"]]))
stopifnot(nrow(M) == 20049L, nrow(L) == 20049L,
          !anyDuplicated(M$stay_id), !anyDuplicated(L$stay_id))
ii <- match(L$stay_id, M$stay_id)
stopifnot(!anyNA(ii))
M <- M[ii]
stopifnot(identical(M$stay_id, L$stay_id))
Mprep <- apply_frozen_recipe(M, recipe, features)
L1_M <- as.integer(assign_nearest_centroid(
  Mprep$X, params$original_centroids)$assigned_cluster)
L2_M <- as.integer(predict_ehr_audit_k2(old, M, return_diagnostics = FALSE)$assigned_cluster)
L0_M <- as.integer(L$cluster_k2)
stopifnot(sum(L1_M != L0_M) == 224L,
          sum(L2_M != L0_M) == 282L,
          sum(L1_M == 1L) == 3902L,
          sum(L2_M == 1L) == 3976L)
Mlabels <- data.table(
  stay_id = L$stay_id,
  locked_k2 = L0_M,
  l1_k2 = L1_M,
  l2_k2 = L2_M,
  missing_feature_n = Mprep$missing_n
)
fwrite(Mlabels, file.path(private_out, "mimic_L0_L1_L2_labels_local.csv"))
mimic_missingness <- Mlabels[, .(
  n = .N,
  locked_c1_n = sum(locked_k2 == 1L),
  l1_c1_n = sum(l1_k2 == 1L),
  flip_n = sum(locked_k2 != l1_k2),
  c1_to_c2_n = sum(locked_k2 == 1L & l1_k2 == 2L),
  c2_to_c1_n = sum(locked_k2 == 2L & l1_k2 == 1L)
), by = .(missingness_stratum = missing_band(missing_feature_n))]
mimic_missingness[, `:=`(
  flip_rate_total = flip_n / n,
  c1_to_c2_rate_among_locked_c1 = fifelse(locked_c1_n > 0L,
                                          c1_to_c2_n / locked_c1_n, NA_real_),
  c2_to_c1_rate_among_locked_c2 = fifelse(n - locked_c1_n > 0L,
                                          c2_to_c1_n / (n - locked_c1_n),
                                          NA_real_)
)]
setorder(mimic_missingness, missingness_stratum)
stopifnot(sum(mimic_missingness$flip_n) == 224L)
fwrite(mimic_missingness,
       file.path(aggregate_out, "mimic_L0_to_L1_flip_by_missingness.csv"))

# eICU: only IDs, feature inputs, and feature-derivation provenance are read.
provenance_cols <- c("alc_exact_match_count", "alc_time_matching_verified",
                     "alc_derivation")
E <- fread(input_paths[["eicu"]],
           select = unique(c("patientunitstayid", "uniquepid", "hospitalid",
                             features, provenance_cols)), showProgress = FALSE)
stopifnot(nrow(E) == 17465L, uniqueN(E$uniquepid) == 16212L,
          uniqueN(E$hospitalid) == 199L,
          !anyDuplicated(E$patientunitstayid))
verified_flag <- E$alc_time_matching_verified
if (is.numeric(verified_flag)) verified_flag <- verified_flag == 1
if (!is.logical(verified_flag)) {
  verified_flag <- tolower(trimws(as.character(verified_flag))) %in%
    c("true", "t", "1", "yes", "y")
}
stopifnot(!anyNA(verified_flag), all(verified_flag))
derivation <- paste(unique(na.omit(E$alc_derivation)), collapse = " | ")
stopifnot(grepl("WBC", derivation, ignore.case = TRUE),
          grepl("lymph", derivation, ignore.case = TRUE),
          grepl("labresultoffset", derivation, ignore.case = TRUE))
alc <- suppressWarnings(as.numeric(E$abs_lymphocytes_min))
alc_observed <- is.finite(alc)
match_count <- suppressWarnings(as.numeric(E$alc_exact_match_count))
stopifnot(any(alc_observed), all(alc[alc_observed] >= 0),
          all(is.finite(match_count[alc_observed])),
          all(match_count[alc_observed] >= 1))
coverage <- data.table(
  feature_order = seq_along(features), feature = features,
  observed_n = vapply(features, function(j) sum(is.finite(
    suppressWarnings(as.numeric(E[[j]])))), integer(1))
)
coverage[, observed_proportion := observed_n / nrow(E)]
stopifnot(all(coverage$observed_n > 0L))
fwrite(coverage, file.path(aggregate_out, "eicu_feature_coverage.csv"))

Eprep <- apply_frozen_recipe(E, recipe, features)
L1_E <- as.integer(assign_nearest_centroid(
  Eprep$X, params$original_centroids)$assigned_cluster)
L2_archive <- fread(input_paths[["l2_eicu"]],
                    select = c("patientunitstayid", "assigned_cluster",
                               "missing_feature_n"), showProgress = FALSE)
stopifnot(nrow(L2_archive) == 17465L,
          !anyDuplicated(L2_archive$patientunitstayid))
jj <- match(E$patientunitstayid, L2_archive$patientunitstayid)
stopifnot(!anyNA(jj))
L2_archive <- L2_archive[jj]
stopifnot(identical(E$patientunitstayid, L2_archive$patientunitstayid),
          identical(as.integer(Eprep$missing_n),
                    as.integer(L2_archive$missing_feature_n)),
          sum(L2_archive$assigned_cluster == 1L) == 3462L)
Elabels <- data.table(
  patientunitstayid = E$patientunitstayid,
  uniquepid = E$uniquepid,
  hospitalid = E$hospitalid,
  l1_k2 = L1_E,
  l2_k2 = as.integer(L2_archive$assigned_cluster),
  missing_feature_n = Eprep$missing_n,
  clipped_feature_n = Eprep$clipped_n
)
fwrite(Elabels, file.path(private_out, "eicu_L1_L2_labels_local.csv"))

pooled <- rbindlist(list(
  data.table(cohort = "MIMIC-IV", n = nrow(Mlabels),
             locked_c1_n = sum(Mlabels$locked_k2 == 1L),
             l1_c1_n = sum(Mlabels$l1_k2 == 1L),
             l2_c1_n = sum(Mlabels$l2_k2 == 1L),
             l1_l2_different_n = sum(Mlabels$l1_k2 != Mlabels$l2_k2)),
  data.table(cohort = "eICU", n = nrow(Elabels), locked_c1_n = NA_integer_,
             l1_c1_n = sum(Elabels$l1_k2 == 1L),
             l2_c1_n = sum(Elabels$l2_k2 == 1L),
             l1_l2_different_n = sum(Elabels$l1_k2 != Elabels$l2_k2))
), use.names = TRUE)
pooled[, `:=`(l1_c1_prevalence = l1_c1_n / n,
             l2_c1_prevalence = l2_c1_n / n,
             locked_c1_prevalence = locked_c1_n / n)]
fwrite(pooled, file.path(aggregate_out, "L1_L2_pooled_prevalence.csv"))

eicu_missingness <- Elabels[, .(
  n = .N,
  l1_c1_n = sum(l1_k2 == 1L),
  l2_c1_n = sum(l2_k2 == 1L),
  l1_l2_different_n = sum(l1_k2 != l2_k2),
  l2_c1_to_l1_c2_n = sum(l2_k2 == 1L & l1_k2 == 2L),
  l2_c2_to_l1_c1_n = sum(l2_k2 == 2L & l1_k2 == 1L)
), by = .(missingness_stratum = missing_band(missing_feature_n))]
eicu_missingness[, `:=`(l1_c1_prevalence = l1_c1_n / n,
                        l2_c1_prevalence = l2_c1_n / n,
                        l1_l2_difference_rate = l1_l2_different_n / n)]
setorder(eicu_missingness, missingness_stratum)
fwrite(eicu_missingness,
       file.path(aggregate_out, "eicu_L1_L2_prevalence_by_missingness.csv"))

hospital <- Elabels[, .(
  n_stays = .N,
  l1_c1_n = sum(l1_k2 == 1L),
  l2_c1_n = sum(l2_k2 == 1L)
), by = hospitalid]
hospital[, `:=`(l1_c1_prevalence = l1_c1_n / n_stays,
               l2_c1_prevalence = l2_c1_n / n_stays)]
mimic_l1_prevalence <- pooled[cohort == "MIMIC-IV", l1_c1_prevalence]
mimic_l2_prevalence <- pooled[cohort == "MIMIC-IV", l2_c1_prevalence]
hospital[, `:=`(
  l1_drift_gt10pp = abs(l1_c1_prevalence - mimic_l1_prevalence) > 0.10,
  l2_drift_gt10pp = abs(l2_c1_prevalence - mimic_l2_prevalence) > 0.10
)]
stopifnot(nrow(hospital) == 199L, sum(hospital$n_stays) == 17465L)
fwrite(hospital, file.path(private_out, "hospital_L1_L2_prevalence_local.csv"))
hospital_summary <- data.table(
  rule = c("L1", "L2"),
  hospitals = 199L,
  median_prevalence = c(median(hospital$l1_c1_prevalence),
                        median(hospital$l2_c1_prevalence)),
  q25_prevalence = c(as.numeric(quantile(hospital$l1_c1_prevalence, 0.25)),
                     as.numeric(quantile(hospital$l2_c1_prevalence, 0.25))),
  q75_prevalence = c(as.numeric(quantile(hospital$l1_c1_prevalence, 0.75)),
                     as.numeric(quantile(hospital$l2_c1_prevalence, 0.75))),
  min_prevalence = c(min(hospital$l1_c1_prevalence),
                     min(hospital$l2_c1_prevalence)),
  max_prevalence = c(max(hospital$l1_c1_prevalence),
                     max(hospital$l2_c1_prevalence)),
  beyond_10pp_n = c(sum(hospital$l1_drift_gt10pp),
                    sum(hospital$l2_drift_gt10pp))
)
fwrite(hospital_summary,
       file.path(aggregate_out, "hospital_L1_L2_prevalence_summary.csv"))

input_manifest <- data.table(
  role = names(input_paths), path = unname(input_paths),
  md5 = unname(actual_md5), sha256 = vapply(input_paths, sha, character(1))
)
fwrite(input_manifest, file.path(aggregate_out, "assignment_input_manifest.csv"))
cat("L1_ASSIGNMENT_PASS: MIMIC", nrow(Mlabels), "eICU", nrow(Elabels),
    "eICU_changed", sum(Elabels$l1_k2 != Elabels$l2_k2), "\n")
