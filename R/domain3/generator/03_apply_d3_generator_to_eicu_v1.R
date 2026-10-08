## =====================================================================
## Apply the frozen regenerated generator to eICU without reading outcomes.
## The external label artifact is frozen and hashed before outcome analysis.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
require_namespace("data.table")
library(data.table)

if (!file.exists(PATH_GENERATOR)) {
  stop("Frozen generator is missing. Run scripts 01 and 02 first.")
}
if (!file.exists(file.path(OUT_QC, "generator_v1_invariance_tests.csv"))) {
  stop("Internal invariance audit is missing. Run script 02 first.")
}
invariance <- fread(file.path(OUT_QC, "generator_v1_invariance_tests.csv"))
if (!all(invariance$pass)) {
  stop("Internal invariance audit did not pass.")
}

object <- readRDS(PATH_GENERATOR)
generator_file_md5 <- hash_file(PATH_GENERATOR)
build_manifest <- fread(file.path(OUT_QC, "generator_v1_build_manifest.csv"))
if (!identical(generator_file_md5, build_manifest$generator_file_md5[1])) {
  stop("Generator file hash changed after construction.")
}
if (!identical(hash_file(PATH_EICU), unname(LOCKED_MD5["eicu_external"]))) {
  stop("eICU input fingerprint differs from the locked configuration.")
}

provenance_cols <- c(
  "alc_exact_match_count", "alc_time_matching_verified", "alc_derivation"
)
assignment_cols <- unique(c(
  "patientunitstayid", "uniquepid", "hospitalid", FEATURES33,
  provenance_cols
))

## Deliberately exclude all outcomes and treatment variables here.
forbidden_cols <- c(
  "hosp_mortality", "icu_mortality", "scr_discharge", "make_hosp",
  "hospitaldischargeoffset", "unitdischargeoffset", "rrt", "make30"
)
if (any(assignment_cols %in% forbidden_cols)) {
  stop("Assignment-stage schema includes a forbidden outcome column.")
}

E <- fread(PATH_EICU, select = assignment_cols, showProgress = FALSE)
if (nrow(E) != EXPECTED_EICU_STAYS) {
  stop("eICU stay denominator mismatch: ", nrow(E))
}
if (uniqueN(E$uniquepid) != EXPECTED_EICU_PATIENTS) {
  stop("eICU patient denominator mismatch: ", uniqueN(E$uniquepid))
}
if (uniqueN(E$hospitalid) != EXPECTED_EICU_HOSPITALS) {
  stop("eICU hospital denominator mismatch: ", uniqueN(E$hospitalid))
}
if (anyDuplicated(E$patientunitstayid)) {
  stop("eICU patientunitstayid is not unique.")
}
validate_feature_schema(E, FEATURES33, "eICU assignment data")

verified_flag <- E$alc_time_matching_verified
if (is.logical(verified_flag)) {
  verified_flag <- verified_flag
} else if (is.numeric(verified_flag)) {
  verified_flag <- verified_flag == 1
} else {
  verified_flag <- tolower(trimws(as.character(verified_flag))) %in%
    c("true", "t", "1", "yes", "y")
}
if (anyNA(verified_flag) || !all(verified_flag)) {
  stop("Absolute lymphocyte time-matching provenance is not verified.")
}
derivation_text <- paste(unique(na.omit(E$alc_derivation)), collapse = " | ")
if (!grepl("WBC", derivation_text, ignore.case = TRUE) ||
    !grepl("lymph", derivation_text, ignore.case = TRUE) ||
    !grepl("labresultoffset", derivation_text, ignore.case = TRUE)) {
  stop("Absolute lymphocyte derivation text is incomplete.")
}
alc <- suppressWarnings(as.numeric(E$abs_lymphocytes_min))
alc_observed <- is.finite(alc)
if (!any(alc_observed) || any(alc[alc_observed] < 0)) {
  stop("Absolute lymphocyte values are absent or implausible.")
}
match_count <- suppressWarnings(as.numeric(E$alc_exact_match_count))
if (any(alc_observed & (!is.finite(match_count) | match_count < 1))) {
  stop("Observed absolute lymphocyte values lack exact-match provenance.")
}

coverage <- data.table(
  feature_order = seq_along(FEATURES33),
  feature = FEATURES33,
  observed_n = vapply(
    FEATURES33,
    function(j) sum(is.finite(suppressWarnings(as.numeric(E[[j]])))),
    integer(1)
  )
)
coverage[, observed_proportion := observed_n / nrow(E)]
if (any(coverage$observed_n == 0L)) {
  stop(
    "At least one frozen feature has zero external coverage: ",
    paste(coverage[observed_n == 0L, feature], collapse = ", ")
  )
}

pred <- predict_ehr_audit_k2(object, E, id_col = "patientunitstayid")
labels <- data.table(
  patientunitstayid = as.integer(E$patientunitstayid),
  uniquepid = as.character(E$uniquepid),
  hospitalid = as.integer(E$hospitalid),
  assigned_cluster = as.integer(as.character(pred$assigned_cluster)),
  distance_c1 = pred$distance_c1,
  distance_c2 = pred$distance_c2,
  signed_margin_c1 = pred$signed_margin_c1,
  absolute_margin = pred$absolute_margin,
  missing_feature_n = pred$missing_feature_n,
  clipped_feature_n = pred$clipped_feature_n,
  generator_version = object$generator_version,
  generator_content_hash = object$generator_content_hash
)
if (anyNA(labels$assigned_cluster) ||
    !all(labels$assigned_cluster %in% c(1L, 2L))) {
  stop("External generator produced invalid labels.")
}
if (labels[assigned_cluster == 1L, .N] != EXPECTED_EICU_C1_N) {
  discrepancy <- labels[, .N, by = assigned_cluster]
  write_replace_csv(
    discrepancy, file.path(OUT_QC, "eicu_generator_v1_count_discrepancy.csv")
  )
  stop("External C1 count differs from the locked regression reference.")
}

legacy <- fread(
  PATH_LEGACY_EICU_LABELS,
  select = c("patientunitstayid", "cluster_k2"),
  showProgress = FALSE
)
if (anyDuplicated(legacy$patientunitstayid)) {
  stop("Legacy eICU label file contains duplicate IDs.")
}
cmp <- merge(
  labels[, .(patientunitstayid, assigned_cluster)],
  legacy,
  by = "patientunitstayid", all = TRUE
)
comparison <- data.table(
  n_new = nrow(labels),
  n_legacy = nrow(legacy),
  n_matched_id = sum(!is.na(cmp$assigned_cluster) & !is.na(cmp$cluster_k2)),
  n_equal_label = sum(
    cmp$assigned_cluster == as.integer(cmp$cluster_k2), na.rm = TRUE
  ),
  n_unequal_label = sum(
    cmp$assigned_cluster != as.integer(cmp$cluster_k2), na.rm = TRUE
  ),
  agreement = mean(
    cmp$assigned_cluster == as.integer(cmp$cluster_k2), na.rm = TRUE
  ),
  legacy_file_md5 = hash_file(PATH_LEGACY_EICU_LABELS),
  expected_legacy_file_md5 = unname(LOCKED_MD5["legacy_eicu_labels"])
)
if (comparison$n_matched_id != EXPECTED_EICU_STAYS ||
    comparison$n_unequal_label != 0L ||
    comparison$legacy_file_md5 != comparison$expected_legacy_file_md5) {
  write_replace_csv(
    comparison,
    file.path(OUT_QC, "eicu_d3_generator_v1_legacy_assignment_comparison.csv")
  )
  stop("New external labels do not exactly reproduce the legacy surrogate.")
}

if (file.exists(PATH_EXTERNAL_LABELS) || file.exists(PATH_EXTERNAL_FREEZE)) {
  stop(
    "Frozen external label output already exists. Refusing to overwrite. ",
    "Delete only after a documented release reset."
  )
}
write_atomic_csv(labels, PATH_EXTERNAL_LABELS)
label_md5 <- hash_file(PATH_EXTERNAL_LABELS)
freeze_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)

fingerprint <- data.table(
  label_path = PATH_EXTERNAL_LABELS,
  label_file_md5 = label_md5,
  generator_path = PATH_GENERATOR,
  generator_file_md5 = generator_file_md5,
  generator_content_hash = object$generator_content_hash,
  eicu_input_path = PATH_EICU,
  eicu_input_md5 = hash_file(PATH_EICU),
  n_icu_stays = nrow(labels),
  n_unique_patients = uniqueN(labels$uniquepid),
  n_hospitals = uniqueN(labels$hospitalid),
  n_c1 = labels[assigned_cluster == 1L, .N],
  n_c2 = labels[assigned_cluster == 2L, .N],
  external_outcomes_read = FALSE,
  eicu_reclustered = FALSE,
  frozen_utc = freeze_time
)
write_atomic_csv(fingerprint, PATH_EXTERNAL_FREEZE)

provenance <- data.table(
  item = c(
    "assignment_rule", "missingness_policy", "distance_metric",
    "feature_n", "alc_derivation", "alc_observed_n",
    "alc_observed_proportion", "outcome_columns_loaded", "eicu_reclustered"
  ),
  value = c(
    object$assignment_rule, object$missingness_policy,
    object$distance_metric, length(FEATURES33), derivation_text,
    sum(alc_observed), mean(alc_observed), "none", "false"
  )
)

write_replace_csv(
  coverage, file.path(OUT_QC, "eicu_d3_generator_v1_feature_coverage.csv")
)
write_replace_csv(
  comparison,
  file.path(OUT_QC, "eicu_d3_generator_v1_legacy_assignment_comparison.csv")
)
write_replace_csv(
  provenance,
  file.path(OUT_QC, "eicu_d3_generator_v1_assignment_provenance.csv")
)
write_replace_csv(
  data.table(
    role = c("generator", "eicu_features", "frozen_external_labels"),
    path = c(PATH_GENERATOR, PATH_EICU, PATH_EXTERNAL_LABELS),
    md5 = c(generator_file_md5, hash_file(PATH_EICU), label_md5)
  ),
  file.path(OUT_QC, "eicu_d3_generator_v1_input_fingerprint.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_LOGS, "eicu_d3_generator_v1_assignment_sessionInfo.txt"),
  useBytes = TRUE
)

cat("External labels frozen without outcome access.\n")
print(fingerprint)
cat("Legacy assignment agreement:", comparison$agreement, "\n")

