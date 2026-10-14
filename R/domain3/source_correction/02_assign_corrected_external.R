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

sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE,
                                    serialize = FALSE))
spec <- file.path(spec_root, "D3_EXTERNAL_TEMPERATURE_CORRECTION_SPEC.md")
expected_spec <- "E6AFB6FD2634ABCA03C77FC06DDCB52B099468335F8D984790AFCA661BAD56B8"
stopifnot(identical(sha(spec), expected_spec))

old_extract_path <- file.path(base, "eicu_external.csv")
temp_path <- file.path(private_dir, "eicu_temperature_corrected_v1.csv")
param_path <- file.path(recovery_root,
                        "outputs/original_parameters_aggregate_v1.rds")
old_labels_path <- file.path(source_root, "private_local/eicu_L1_L2_labels_local.csv")
stopifnot(all(file.exists(c(old_extract_path, temp_path, param_path,
                             old_labels_path))))
stopifnot(identical(sha(old_extract_path),
                    "249F28CAC3C900B2F8B187D8CC3F2887230DB7B6BFE54EAB19717377F223B5A5"))
stopifnot(identical(sha(param_path),
                    "70CEAD3058658339D3AE00753F7B21045F0BBCCB74B8E9267FC2F5B906931C99"))
source(file.path(repo, "R/domain3/generator",
                 "00_d3_generator_utils_v1.R"))

E <- fread(old_extract_path, showProgress = FALSE)
T <- fread(temp_path, showProgress = FALSE)
L <- fread(old_labels_path, showProgress = FALSE)
stopifnot(nrow(E) == 17465L, nrow(T) == nrow(E), nrow(L) == nrow(E),
          !anyDuplicated(E$patientunitstayid),
          !anyDuplicated(T$patientunitstayid),
          !anyDuplicated(L$patientunitstayid))
ti <- match(E$patientunitstayid, T$patientunitstayid)
li <- match(E$patientunitstayid, L$patientunitstayid)
stopifnot(!anyNA(ti), !anyNA(li))
T <- T[ti]
L <- L[li]
stopifnot(identical(E$patientunitstayid, T$patientunitstayid),
          identical(E$patientunitstayid, L$patientunitstayid))
stopifnot(all(T$temperature_source %in%
                c("nurse_c", "nurse_f_converted", "periodic_fallback", "missing")))
stopifnot(all(T$temperature_min_corrected <= T$temperature_max_corrected,
              na.rm = TRUE),
          all(T$temperature_min_corrected >= 25 &
                T$temperature_max_corrected <= 45, na.rm = TRUE),
          all((is.na(T$temperature_min_corrected) &
                 is.na(T$temperature_max_corrected)) ==
                (T$temperature_source == "missing")))

params <- readRDS(param_path)
features <- params$feature_names
stopifnot(length(features) == 33L, !anyDuplicated(features),
          all(features %in% names(E)),
          all(c("temperature_min", "temperature_max") %in% features))
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

Enew <- copy(E)
Enew[, `:=`(temperature_min = T$temperature_min_corrected,
            temperature_max = T$temperature_max_corrected)]
unchanged <- setdiff(names(E), c("temperature_min", "temperature_max"))
stopifnot(identical(E[, ..unchanged], Enew[, ..unchanged]))
prep <- apply_frozen_recipe(Enew, recipe, features)
new_labels <- as.integer(assign_nearest_centroid(
  prep$X, params$original_centroids)$assigned_cluster)
old_labels <- as.integer(L$l1_k2)
stopifnot(all(new_labels %in% 1:2), all(old_labels %in% 1:2))

P <- data.table(patientunitstayid = E$patientunitstayid,
                hospitalid = E$hospitalid,
                old_l1_k2 = old_labels,
                corrected_l1_k2 = new_labels,
                old_missing_n = L$missing_feature_n,
                corrected_missing_n = prep$missing_n,
                temperature_source = T$temperature_source)
fwrite(P, file.path(private_dir, "eicu_corrected_L1_labels_local.csv"))

source_summary <- P[, .(stays = .N,
                        old_c1 = sum(old_l1_k2 == 1L),
                        corrected_c1 = sum(corrected_l1_k2 == 1L),
                        changed = sum(old_l1_k2 != corrected_l1_k2)),
                    by = temperature_source][order(-stays)]
fwrite(source_summary, file.path(aggregate_dir, "temperature_source_summary.csv"))

summary <- data.table(
  stays = nrow(P),
  old_temp_observed = sum(!is.na(E$temperature_min)),
  corrected_temp_observed = sum(!is.na(Enew$temperature_min)),
  old_c1 = sum(old_labels == 1L),
  corrected_c1 = sum(new_labels == 1L),
  c1_to_c2 = sum(old_labels == 1L & new_labels == 2L),
  c2_to_c1 = sum(old_labels == 2L & new_labels == 1L),
  changed = sum(old_labels != new_labels)
)
summary[, `:=`(old_c1_rate = old_c1 / stays,
               corrected_c1_rate = corrected_c1 / stays)]
fwrite(summary, file.path(aggregate_dir, "corrected_assignment_summary.csv"))

P[, missingness_stratum := fifelse(corrected_missing_n == 0L, "0",
  fifelse(corrected_missing_n <= 2L, "1-2",
    fifelse(corrected_missing_n <= 5L, "3-5", "6+")))]
missing_summary <- P[, .(stays = .N,
                         old_c1 = sum(old_l1_k2 == 1L),
                         corrected_c1 = sum(corrected_l1_k2 == 1L),
                         c1_to_c2 = sum(old_l1_k2 == 1L & corrected_l1_k2 == 2L),
                         c2_to_c1 = sum(old_l1_k2 == 2L & corrected_l1_k2 == 1L)),
                     by = missingness_stratum]
fwrite(missing_summary,
       file.path(aggregate_dir, "corrected_assignment_by_missingness.csv"))

hospital <- P[, .(stays = .N,
                  old_c1 = sum(old_l1_k2 == 1L),
                  corrected_c1 = sum(corrected_l1_k2 == 1L)),
              by = hospitalid]
stopifnot(nrow(hospital) == 199L, sum(hospital$stays) == nrow(P))
fwrite(hospital, file.path(private_dir, "corrected_hospital_counts_local.csv"))
hospital[, corrected_prevalence := corrected_c1 / stays]
hospital_summary <- data.table(
  hospitals = nrow(hospital),
  stays = sum(hospital$stays),
  median_prevalence = median(hospital$corrected_prevalence),
  q25_prevalence = as.numeric(quantile(hospital$corrected_prevalence, .25)),
  q75_prevalence = as.numeric(quantile(hospital$corrected_prevalence, .75)))
fwrite(hospital_summary,
       file.path(aggregate_dir, "corrected_hospital_summary.csv"))

input_manifest <- data.table(
  role = c("spec", "old_external_extract", "corrected_temperature",
           "recovered_parameters", "old_l1_labels"),
  sha256 = vapply(c(spec, old_extract_path, temp_path, param_path,
                    old_labels_path), sha, character(1)))
fwrite(input_manifest, file.path(aggregate_dir, "assignment_input_hashes.csv"))

print(summary)
print(source_summary)
print(hospital_summary)
cat("ASSIGNMENT_COMPLETE_WITHOUT_OUTCOMES\n")
