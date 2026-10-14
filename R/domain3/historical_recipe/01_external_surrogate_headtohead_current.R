## =====================================================================
## 19_external_headtohead_mice_labels_v2.R
## MIMIC-trained models frozen and transported to eICU, using the MICE primary
## K=2 label as the subphenotype label source.
##
## Changes from older 19_external_headtohead.R:
##   - Reads output/model/labels_primary_mice.rds, not labels_primary.rds.
##   - Stops if MICE labels are missing, lack IDs, or lack cluster_k2.
##   - Writes label provenance and high-risk-label prevalence drift.
##   - Uses high-risk label defined by MIMIC in-hospital mortality rather than
##     assuming the first sorted label is C1.
##   - Requires the upstream SQL-verified 0-24h absolute lymphocyte count,
##     derived from WBC x lymphocyte percentage / 100 at identical offsets.
##   - Requires all 33 frozen MIMIC features for external assignment; patient-
##     level missing values are imputed by the frozen MIMIC medians.
##   - Retains the frozen MIMIC recipe and centroids in the transport object.
##   - Writes feature-level provenance, input hashes, ICU-stay/patient counts,
##     and patient-clustered external AUC confidence intervals.
## =====================================================================
## Resolve 00_paths.R relative to this script, not to getwd(). Calling
## source() with an absolute script path does not change the working directory.
.script_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
.script_dir <- if (!is.null(.script_ofile) && length(.script_ofile) && !is.na(.script_ofile)) {
  dirname(normalizePath(.script_ofile, winslash = "/", mustWork = FALSE))
} else {
  NA_character_
}
.paths_candidates <- unique(c(
  if (!is.na(.script_dir)) file.path(.script_dir, "00_paths.R") else NA_character_,
  file.path(getwd(), "00_paths.R"),
  file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "r", "00_paths.R")
))
.paths_candidates <- .paths_candidates[!is.na(.paths_candidates)]
.paths_file <- .paths_candidates[file.exists(.paths_candidates)][1]
if (is.na(.paths_file)) {
  stop(
    "Cannot locate 00_paths.R. Current working directory: ", getwd(),
    "\nChecked: ", paste(.paths_candidates, collapse = "; "),
    call. = FALSE
  )
}
source(.paths_file, chdir = FALSE)
message("Loaded project paths from: ", normalizePath(.paths_file, winslash = "/"))
rm(.script_ofile, .script_dir, .paths_candidates, .paths_file)
suppressPackageStartupMessages({
  need <- c("data.table", "glmnet", "ranger", "pROC")
  for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(glmnet)
  library(ranger)
  library(pROC)
})
set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SCRIPT_VERSION <- "external-headtohead-mice-v2-strict33-alc-20260713"
EXPECTED_EICU_STAYS <- 17465L
EXPECTED_EICU_PATIENTS <- 16212L
STOP_ON_EICU_COUNT_MISMATCH <- TRUE
LYMPH_FEATURE <- "abs_lymphocytes_min"

K_FOLDS <- 10L
RF_TREES <- 300L
ENET_ALPHA <- 0.5
BOOT_B <- 1000L
DO_DCA <- TRUE
DCA_THRESH <- seq(0.02, 0.60, by = 0.01)

FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
  "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml", "glucose_max",
  "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min", "sodium_min",
  "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min", "chloride_min",
  "chloride_max", "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")

## ---------- helpers ----------
fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

boot_auc_ci <- function(y, p, cluster_id = NULL, B = BOOT_B) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  if (is.null(cluster_id)) cluster_id <- seq_along(y) else cluster_id <- as.character(cluster_id[ok])
  if (anyNA(cluster_id)) stop("Missing cluster ID in external AUC bootstrap.")
  clusters <- unique(cluster_id)
  rows_by_cluster <- split(seq_along(cluster_id), cluster_id)
  a <- numeric(B)
  for (b in seq_len(B)) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    i <- unlist(rows_by_cluster[sampled_clusters], use.names = FALSE)
    a[b] <- fast_auc(y[i], p[i])
  }
  quantile(a, c(.025, .975), na.rm = TRUE)
}

clamp <- function(p, e = 1e-6) pmin(pmax(p, e), 1 - e)

cal_slope_int <- function(y, p) {
  lp <- qlogis(clamp(p))
  s <- tryCatch(coef(glm(y ~ lp, family = binomial))[2], error = function(e) NA_real_)
  i <- tryCatch(coef(glm(y ~ 1, offset = lp, family = binomial))[1], error = function(e) NA_real_)
  c(slope = unname(s), intercept = unname(i))
}

brier <- function(y, p) mean((p - y)^2)

net_benefit <- function(y, p, th) {
  n <- length(y)
  sapply(th, function(t) {
    pos <- p >= t
    TP <- sum(pos & y == 1)
    FP <- sum(pos & y == 0)
    TP / n - FP / n * (t / (1 - t))
  })
}

parse_age <- function(a) {
  if (is.numeric(a)) return(a)
  a <- as.character(a)
  a <- gsub(">\\s*89", "90", a)
  suppressWarnings(as.numeric(a))
}

parse_sex <- function(g) {
  if (is.numeric(g)) return(as.integer(g == 1))
  x <- toupper(substr(trimws(as.character(g)), 1, 1))
  fifelse(x == "M", 1L, fifelse(x == "F", 0L, NA_integer_))
}

winsor <- function(x, q) pmin(pmax(x, q[1]), q[2])

hash_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

hash_object <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf, version = 3)
  unname(tools::md5sum(tf))
}

derive_verified_alc <- function(E) {
  ## The upstream eICU SQL must derive abs_lymphocytes_min from WBC and
  ## lymphocyte percentage measured at the identical labresultoffset. The old
  ## eICU column was only the '-lymphs' percentage and must never be accepted.
  required <- c(
    LYMPH_FEATURE,
    "alc_exact_match_count",
    "alc_time_matching_verified",
    "alc_derivation"
  )
  missing_required <- setdiff(required, names(E))
  if (length(missing_required)) {
    stop(
      "eICU absolute lymphocyte provenance is incomplete. Missing: ",
      paste(missing_required, collapse = ", "),
      "\nRe-run the corrected eicu_06_features.sql, then eicu_07_sofa.sql and eicu_09_dedup.sql before script 19.",
      call. = FALSE
    )
  }

  verified_flag <- E[["alc_time_matching_verified"]]
  if (is.logical(verified_flag)) {
    verified_flag <- verified_flag
  } else if (is.numeric(verified_flag)) {
    verified_flag <- verified_flag == 1
  } else {
    verified_flag <- tolower(trimws(as.character(verified_flag))) %in% c("true", "t", "1", "yes", "y")
  }
  if (anyNA(verified_flag) || !all(verified_flag)) {
    stop("alc_time_matching_verified is not TRUE for every eICU row.", call. = FALSE)
  }

  derivation_values <- unique(na.omit(as.character(E[["alc_derivation"]])))
  derivation_text <- paste(derivation_values, collapse = " | ")
  if (!length(derivation_values) ||
      !grepl("WBC", derivation_text, ignore.case = TRUE) ||
      !grepl("lymph", derivation_text, ignore.case = TRUE) ||
      !grepl("labresultoffset", derivation_text, ignore.case = TRUE)) {
    stop("alc_derivation does not document exact-offset WBC/lymphocyte matching.", call. = FALSE)
  }

  alc <- suppressWarnings(as.numeric(E[[LYMPH_FEATURE]]))
  observed <- is.finite(alc)
  if (!any(observed)) stop("Verified eICU absolute lymphocyte column has no observed values.", call. = FALSE)
  if (any(alc[observed] < 0)) stop("Verified eICU absolute lymphocyte column contains negative values.", call. = FALSE)
  alc_p99 <- unname(stats::quantile(alc[observed], 0.99, na.rm = TRUE))
  alc_median <- stats::median(alc[observed], na.rm = TRUE)
  if (!is.finite(alc_median) || alc_median > 10 || !is.finite(alc_p99) || alc_p99 > 50) {
    stop(
      sprintf("Verified ALC distribution is implausible (median %.3f, p99 %.3f); percentage values may still be present.",
              alc_median, alc_p99),
      call. = FALSE
    )
  }

  match_count <- suppressWarnings(as.numeric(E[["alc_exact_match_count"]]))
  if (any(observed & (!is.finite(match_count) | match_count < 1))) {
    stop("Observed verified ALC values lack a positive exact-match count.", call. = FALSE)
  }

  list(
    data = E,
    verified = TRUE,
    derivation = derivation_text,
    source_columns = paste(required, collapse = ";"),
    coverage = mean(observed),
    observed_n = sum(observed),
    imputed_n = sum(!observed),
    median = alc_median,
    p99 = alc_p99
  )
}

check_unique_id <- function(dt, id_col, name) {
  if (!(id_col %in% names(dt))) stop(name, " missing ID column: ", id_col)
  nd <- anyDuplicated(dt[[id_col]])
  if (nd) stop(name, " contains duplicated ", id_col, " values; first duplicated index = ", nd)
}

load_mice_k2_labels <- function(id_col) {
  label_path <- PATH_MICE_LABELS
  if (!file.exists(label_path)) {
    stop("Missing MICE label file: ", label_path,
         "\nRun 23_mice_primary_pooling_v2_fixed.R or 23c_save_mice_matrix_labels.R first.")
  }
  lab <- as.data.table(readRDS(label_path))
  id_candidates <- intersect(c(id_col, "stay_id", "hadm_id", "subject_id"), names(lab))
  if (!(id_col %in% names(lab))) {
    if (length(id_candidates) == 0) stop("labels_primary_mice.rds lacks an ID column compatible with ", id_col)
    setnames(lab, id_candidates[1], id_col)
  }
  if (!"cluster_k2" %in% names(lab)) stop("labels_primary_mice.rds lacks cluster_k2.")
  lab <- lab[, .SD, .SDcols = c(id_col, "cluster_k2")]
  check_unique_id(lab, id_col, "labels_primary_mice.rds")
  lab[, cluster_k2 := factor(cluster_k2)]
  fwrite(data.table(
    artifact = "labels_primary_mice.rds",
    path = label_path,
    n_labels = nrow(lab),
    columns = paste(names(lab), collapse = ";")
  ), file.path(OUT_DIR, "external_mice_label_provenance.csv"))
  lab
}

sqdist <- function(A, B) {
  aa <- rowSums(A^2)
  bb <- rowSums(B^2)
  outer(aa, bb, "+") - 2 * A %*% t(B)
}

PATH_MICE_LABELS <- file.path(DIR_MODEL, "labels_primary_mice.rds")
PATH_BASELINE_COVARS <- file.path(DIR_DATA, "baseline_covars.csv")
PATH_EICU_EXTERNAL <- file.path(DIR_DATA, "eicu_external.csv")

input_fingerprint <- data.table(
  role = c("mimic_final_full", "mimic_labels", "baseline_covariates", "eicu_external"),
  path = c(PATH_FINAL_FULL, PATH_MICE_LABELS, PATH_BASELINE_COVARS, PATH_EICU_EXTERNAL)
)
input_fingerprint[, exists := file.exists(path)]
input_fingerprint[, size_bytes := ifelse(exists, file.info(path)$size, NA_real_)]
input_fingerprint[, modified_utc := ifelse(exists, format(file.info(path)$mtime, tz = "UTC", usetz = TRUE), NA_character_)]
input_fingerprint[, md5 := vapply(path, hash_file, character(1))]
if (any(!input_fingerprint$exists)) {
  stop("Missing required input(s): ", paste(input_fingerprint[!exists, path], collapse = ", "))
}
fwrite(input_fingerprint, file.path(OUT_DIR, "external_assignment_input_fingerprint_mice_labels.csv"))

## ---------- 1. MIMIC frame ----------
mimic <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id", "subject_id"), names(mimic))[1]
if (is.na(mid)) stop("No stay_id/hadm_id/subject_id found in PATH_FINAL_FULL.")

need_ff <- c(mid, "age", "gender", "sofa_score", "hospital_mortality", "make30", FEATURES33)
miss <- setdiff(need_ff, names(mimic))
if (length(miss)) stop("final_full missing: ", paste(miss, collapse = ", "))
M <- mimic[, ..need_ff]
check_unique_id(M, mid, "PATH_FINAL_FULL")
M[, sex := parse_sex(gender)]

bc <- fread(PATH_BASELINE_COVARS)
check_unique_id(bc, mid, "baseline_covars.csv")
if (!"aki_stage_0_24h" %in% names(bc)) stop("baseline_covars.csv missing aki_stage_0_24h")
M <- merge(M, bc[, .SD, .SDcols = c(mid, "aki_stage_0_24h")], by = mid, all.x = TRUE)

labs <- load_mice_k2_labels(mid)
M <- merge(M, labs, by = mid, all.x = TRUE)
M <- M[complete.cases(M[, c(BASE_NUM, "cluster_k2", "hospital_mortality", "make30"), with = FALSE])]
M[, cluster_k2 := factor(cluster_k2)]
CLUST_LEVELS <- levels(M$cluster_k2)

risk_by_label <- M[, .(mortality = mean(hospital_mortality == 1), n = .N), by = cluster_k2]
high_risk_label <- as.character(risk_by_label[which.max(mortality), cluster_k2])
fwrite(risk_by_label, file.path(OUT_DIR, "external_mice_label_mimic_risk_by_label.csv"))

cat(sprintf("MIMIC n=%d | in-hosp death=%d (%.1f%%) | make30=%d (%.1f%%) | high-risk label=%s\n",
            nrow(M), sum(M$hospital_mortality), 100 * mean(M$hospital_mortality),
            sum(M$make30), 100 * mean(M$make30), high_risk_label))

## ---------- 2. eICU frame ----------
eicu <- fread(PATH_EICU_EXTERNAL)
setnames(eicu, "sofa", "sofa_score", skip_absent = TRUE)
E <- copy(eicu)
if (!"patientunitstayid" %in% names(E)) stop("eicu_external.csv missing patientunitstayid")
if (!"uniquepid" %in% names(E)) stop("eicu_external.csv missing uniquepid; clustered external inference is not possible")
if (anyDuplicated(E$patientunitstayid)) stop("eicu_external.csv contains duplicated patientunitstayid values")

n_eicu_stays <- nrow(E)
n_eicu_patients <- uniqueN(E$uniquepid)
eicu_denominator_audit <- data.table(
  n_icu_stays = n_eicu_stays,
  n_unique_patients = n_eicu_patients,
  repeated_stays = n_eicu_stays - n_eicu_patients,
  statistical_unit = "ICU stay",
  clustered_inference_unit = "uniquepid"
)
fwrite(eicu_denominator_audit, file.path(OUT_DIR, "external_eicu_denominator_audit_mice_labels.csv"))
count_mismatch <- n_eicu_stays != EXPECTED_EICU_STAYS || n_eicu_patients != EXPECTED_EICU_PATIENTS
if (count_mismatch) {
  txt <- sprintf(
    "eICU denominator mismatch: observed %d ICU stays/%d patients; expected %d/%d.",
    n_eicu_stays, n_eicu_patients, EXPECTED_EICU_STAYS, EXPECTED_EICU_PATIENTS
  )
  if (isTRUE(STOP_ON_EICU_COUNT_MISMATCH)) stop(txt) else warning(txt)
}

alc_result <- derive_verified_alc(E)
E <- alc_result$data
E[, age := parse_age(age)]
E[, sex := parse_sex(sex)]
E[, aki_stage_0_24h := as.numeric(aki_stage_0_24h)]
E[, hosp_mortality := as.integer(hosp_mortality)]
if (all(c("scr_discharge", "scr_baseline") %in% names(E))) {
  E[, make_hosp := as.integer(hosp_mortality == 1 |
    (!is.na(scr_discharge) & !is.na(scr_baseline) & scr_discharge >= 1.5 * scr_baseline))]
} else {
  warning("scr_discharge/scr_baseline missing -> make_hosp = in-hosp death only")
  E[, make_hosp := hosp_mortality]
}

## ---------- 3. strict frozen 33-feature assignment + MIMIC-coordinate recipe ----------
missing_mimic_features <- setdiff(FEATURES33, names(M))
missing_eicu_features <- setdiff(FEATURES33, names(E))
if (length(missing_mimic_features)) {
  stop("MIMIC data are missing frozen feature(s): ", paste(missing_mimic_features, collapse = ", "))
}
if (length(missing_eicu_features)) {
  stop("eICU data are missing frozen feature(s): ", paste(missing_eicu_features, collapse = ", "))
}
if (!isTRUE(alc_result$verified)) stop("Verified absolute lymphocyte derivation is required for strict 33-feature assignment.")

## Strict transport preserves all 33 frozen dimensions. Missing patient-level
## values are handled below by the frozen MIMIC medians; low external coverage
## is audited but does not delete a dimension.
common <- FEATURES33
cov_e <- sapply(common, function(v) mean(is.finite(suppressWarnings(as.numeric(E[[v]])))))
zero_coverage <- names(cov_e)[!is.finite(cov_e) | cov_e <= 0]
if (length(zero_coverage)) {
  stop(
    "Strict 33-feature assignment cannot proceed because eICU has zero observed coverage for: ",
    paste(zero_coverage, collapse = ", "),
    "\nRebuild eicu_external.csv with the corrected feature SQL (including PaCO2 mapping).",
    call. = FALSE
  )
}
cat("\n---- strict frozen assignment features (", length(common), ") | eICU coverage ----\n", sep = "")
print(round(sort(cov_e), 3))
cat("Absolute lymphocyte handling:", alc_result$derivation, "\n")
cat("Absolute lymphocyte used for assignment:", LYMPH_FEATURE %in% common, "\n")
cat("Assignment dimensionality: strict ", length(common), " of 33 frozen MIMIC features\n", sep = "")

feature_audit <- data.table(
  feature = FEATURES33,
  present_in_mimic = FEATURES33 %in% names(M),
  present_in_eicu = FEATURES33 %in% names(E),
  eicu_coverage = vapply(FEATURES33, function(v) {
    if (!v %in% names(E)) return(NA_real_)
    mean(is.finite(suppressWarnings(as.numeric(E[[v]]))))
  }, numeric(1))
)
feature_audit[, used_for_assignment := feature %in% common]
feature_audit[, exclusion_reason := fifelse(used_for_assignment, "included in strict frozen 33-feature assignment", "not included")]
feature_audit[, `:=`(derivation = NA_character_, source_columns = NA_character_, time_matching_verified = NA)]
feature_audit[feature == LYMPH_FEATURE, `:=`(
  derivation = alc_result$derivation,
  source_columns = alc_result$source_columns,
  time_matching_verified = isTRUE(alc_result$verified),
  exclusion_reason = "included: verified exact-offset derived absolute count"
)]
fwrite(feature_audit, file.path(OUT_DIR, "external_assignment_feature_audit_mice_labels.csv"))

recipe <- lapply(common, function(j) {
  v <- as.numeric(M[[j]])
  q <- quantile(v, c(.01, .99), na.rm = TRUE)
  vw <- winsor(v, q)
  sdv <- sd(vw, na.rm = TRUE)
  if (!is.finite(sdv) || sdv < 1e-8) sdv <- 1
  list(lo = q[1], hi = q[2], mean = mean(vw, na.rm = TRUE), sd = sdv, med = median(vw, na.rm = TRUE))
})
names(recipe) <- common

recipe_table <- rbindlist(lapply(names(recipe), function(j) {
  r <- recipe[[j]]
  data.table(feature = j, winsor_p01 = r$lo, winsor_p99 = r$hi,
             imputation_median = r$med, reference_mean = r$mean, reference_sd = r$sd)
}))
fwrite(recipe_table, file.path(OUT_DIR, "external_assignment_preprocessing_recipe_mice_labels.csv"))
recipe_hash <- hash_object(recipe_table)
common_feature_hash <- hash_object(common)

apply_recipe <- function(df) {
  out <- matrix(0, nrow(df), length(common), dimnames = list(NULL, common))
  for (j in common) {
    r <- recipe[[j]]
    v <- as.numeric(df[[j]])
    v <- winsor(v, c(r$lo, r$hi))
    v[is.na(v)] <- r$med
    out[, j] <- (v - r$mean) / r$sd
  }
  out
}

Xm <- apply_recipe(M)
Xe <- apply_recipe(E)
BASEm <- as.matrix(M[, ..BASE_NUM]); storage.mode(BASEm) <- "double"
BASEe <- as.matrix(E[, ..BASE_NUM]); storage.mode(BASEe) <- "double"

## ---------- 4. eICU phenotype: nearest MIMIC centroid on common features ----------
make_centroids <- function(Xtr, lab) {
  ks <- sort(unique(lab))
  cent <- t(sapply(ks, function(k) colMeans(Xtr[lab == k, , drop = FALSE])))
  rownames(cent) <- ks
  cent
}

assign_centroid <- function(cent, Xnew) {
  ks <- rownames(cent)
  idx <- max.col(-sqdist(Xnew, cent), ties.method = "first")
  ks[idx]
}
frozen_centroids <- make_centroids(Xm, as.character(M$cluster_k2))
centroid_hash <- hash_object(frozen_centroids)
centroid_table <- as.data.table(frozen_centroids, keep.rownames = "source_label")
fwrite(centroid_table, file.path(OUT_DIR, "external_assignment_frozen_centroids_mice_labels.csv"))
E_k2 <- assign_centroid(frozen_centroids, Xe)
E[, cluster_k2 := factor(E_k2, levels = CLUST_LEVELS)]
cat("\neICU K2 assignment from MICE-label centroids:\n")
print(table(E$cluster_k2))
## ---------- Save eICU MICE-label assignment for Table 1 ----------
assignment_run_id <- substr(hash_object(list(
  script_version = SCRIPT_VERSION,
  input_md5 = input_fingerprint[, .(role, md5)],
  common = common,
  recipe_hash = recipe_hash,
  centroid_hash = centroid_hash,
  high_risk_label = high_risk_label
)), 1, 16)

eicu_label_mice <- data.table(
  patientunitstayid = E$patientunitstayid,
  uniquepid = E$uniquepid,
  cluster_k2 = as.integer(as.character(E$cluster_k2)),
  cluster_k2_display = fifelse(as.character(E$cluster_k2) == high_risk_label, 1L, 2L),
  hosp_mortality = E$hosp_mortality,
  make_hosp = E$make_hosp,
  label_source = "labels_primary_mice.rds",
  label_direction_source = "MIMIC-IV in-hospital mortality only",
  high_risk_source_label = high_risk_label,
  assignment_rule = "strict frozen MIMIC 33-feature preprocessing plus nearest frozen MIMIC centroid; no eICU re-clustering",
  common_feature_n = length(common),
  assignment_dimension_rule = "all 33 frozen MIMIC features retained; patient-level missing values imputed with frozen MIMIC medians",
  common_feature_hash = common_feature_hash,
  preprocessing_hash = recipe_hash,
  centroid_hash = centroid_hash,
  assignment_run_id = assignment_run_id
)

fwrite(
  eicu_label_mice,
  file.path(DIR_OUTPUT, "eicu_cluster_labels_mice.csv")
)

fwrite(
  eicu_label_mice,
  file.path(OUT_DIR, "eicu_cluster_labels_mice.csv")
)

eicu_label_hash <- hash_file(file.path(DIR_OUTPUT, "eicu_cluster_labels_mice.csv"))
assignment_provenance <- data.table(
  script_version = SCRIPT_VERSION,
  assignment_run_id = assignment_run_id,
  label_file = file.path(DIR_OUTPUT, "eicu_cluster_labels_mice.csv"),
  label_file_md5 = eicu_label_hash,
  n_icu_stays = nrow(eicu_label_mice),
  n_unique_patients = uniqueN(eicu_label_mice$uniquepid),
  common_feature_n = length(common),
  assignment_dimension_rule = "all 33 frozen MIMIC features retained; patient-level missing values imputed with frozen MIMIC medians",
  common_feature_hash = common_feature_hash,
  preprocessing_hash = recipe_hash,
  centroid_hash = centroid_hash,
  high_risk_source_label = high_risk_label,
  label_direction_source = "MIMIC-IV in-hospital mortality only",
  eicu_reclustered = FALSE,
  lymphocyte_verified = isTRUE(alc_result$verified),
  lymphocyte_used_for_assignment = LYMPH_FEATURE %in% common,
  lymphocyte_observed_n = alc_result$observed_n,
  lymphocyte_imputed_n = alc_result$imputed_n,
  lymphocyte_observed_coverage = alc_result$coverage,
  lymphocyte_observed_median = alc_result$median,
  lymphocyte_observed_p99 = alc_result$p99
)
fwrite(assignment_provenance, file.path(OUT_DIR, "external_assignment_provenance_mice_labels.csv"))

cat("\nSaved eICU MICE-label assignment for Table 1:\n")
print(eicu_label_mice[, .N, by = cluster_k2])

prev <- data.table(
  cohort = c("MIMIC", "eICU"),
  high_risk_label = high_risk_label,
  high_risk_pct = c(
    round(100 * mean(as.character(M$cluster_k2) == high_risk_label), 1),
    round(100 * mean(as.character(E$cluster_k2) == high_risk_label), 1)
  ),
  n = c(nrow(M), nrow(E))
)
fwrite(prev, file.path(OUT_DIR, "headtohead_mice_prevalence_drift.csv"))

## ---------- 5. model fit / predict on common ----------
fit_model <- function(kind, y, baseM, Xstd, clust) switch(kind,
  base      = glm(y ~ ., data = data.frame(y = y, baseM), family = binomial),
  base_phen = glm(y ~ ., data = data.frame(y = y, baseM, cluster_k2 = clust), family = binomial),
  feat_pen  = cv.glmnet(cbind(baseM, Xstd), y, family = "binomial", alpha = ENET_ALPHA, nfolds = 5),
  feat_rf   = ranger(x = as.data.frame(cbind(baseM, Xstd)), y = factor(y), probability = TRUE,
                     num.trees = RF_TREES, num.threads = 0, seed = 1L)
)

pred_model <- function(kind, mdl, baseM, Xstd, clust) switch(kind,
  base      = predict(mdl, newdata = data.frame(baseM), type = "response"),
  base_phen = predict(mdl, newdata = data.frame(baseM, cluster_k2 = clust), type = "response"),
  feat_pen  = as.numeric(predict(mdl, cbind(baseM, Xstd), s = "lambda.min", type = "response")),
  feat_rf   = predict(mdl, as.data.frame(cbind(baseM, Xstd)))$predictions[, "1"]
)
MODELS <- c("base", "base_phen", "feat_pen", "feat_rf")

cv_auc <- function(y) {
  if (anyNA(y) || !all(y %in% c(0, 1))) stop("Internal outcome must be binary without NA.")
  folds <- integer(length(y))
  for (cl in c(0, 1)) {
    ix <- which(y == cl)
    folds[ix] <- sample(rep_len(seq_len(K_FOLDS), length(ix)))
  }
  sapply(MODELS, function(k) {
    p <- numeric(length(y))
    for (f in seq_len(K_FOLDS)) {
      tr <- folds != f
      te <- folds == f
      m <- fit_model(k, y[tr], BASEm[tr, , drop = FALSE], Xm[tr, , drop = FALSE], M$cluster_k2[tr])
      p[te] <- pred_model(k, m, BASEm[te, , drop = FALSE], Xm[te, , drop = FALSE], M$cluster_k2[te])
    }
    fast_auc(y, p)
  })
}

## ---------- 6. per-outcome transport ----------
OUTCOMES <- list(
  mortality = list(mimic = "hospital_mortality", eicu = "hosp_mortality", clean = TRUE),
  make      = list(mimic = "make30", eicu = "make_hosp", clean = FALSE)
)
okE <- stats::complete.cases(BASEe) & !is.na(E$cluster_k2)
cat(sprintf("eICU scored rows: %d of %d\n", sum(okE), nrow(E)))

run <- function(onm) {
  oc <- OUTCOMES[[onm]]
  ym <- as.integer(M[[oc$mimic]])
  ecol <- oc$eicu
  cat(sprintf("\n=== OUTCOME %s | internal=%s external=%s | same-endpoint=%s ===\n",
              onm, oc$mimic, ecol, oc$clean))
  internal <- cv_auc(ym)
  fits <- lapply(MODELS, function(k) fit_model(k, ym, BASEm, Xm, M$cluster_k2))
  names(fits) <- MODELS
  sel <- which(okE & !is.na(E[[ecol]]))
  ye <- as.integer(E[[ecol]][sel])
  cat(sprintf("  external evaluable rows for %s: %d\n", onm, length(sel)))
  rows <- rbindlist(lapply(MODELS, function(k) {
    pe <- pred_model(k, fits[[k]], BASEe[sel, , drop = FALSE], Xe[sel, , drop = FALSE], E$cluster_k2[sel])
    aci <- boot_auc_ci(ye, pe, cluster_id = E$uniquepid[sel])
    cs <- cal_slope_int(ye, pe)
    data.table(
      outcome = onm,
      model = k,
      label_source = "labels_primary_mice.rds",
      internal_cv_auc = round(internal[[k]], 4),
      external_auc = round(fast_auc(ye, pe), 4),
      ext_auc_lo = round(aci[1], 4),
      ext_auc_hi = round(aci[2], 4),
      transport_drop = round(internal[[k]] - fast_auc(ye, pe), 4),
      ext_cal_slope = round(cs["slope"], 3),
      ext_cal_intercept = round(cs["intercept"], 3),
      ext_brier = round(brier(ye, pe), 4),
      external_ci_resampling_unit = "uniquepid cluster",
      n_external_icu_stays = length(sel),
      n_external_unique_patients = uniqueN(E$uniquepid[sel])
    )
  }))
  dca <- NULL
  if (DO_DCA) {
    nb <- rbindlist(lapply(MODELS, function(k) {
      pe <- pred_model(k, fits[[k]], BASEe[sel, , drop = FALSE], Xe[sel, , drop = FALSE], E$cluster_k2[sel])
      data.table(outcome = onm, model = k, threshold = DCA_THRESH,
                 net_benefit = net_benefit(ye, pe, DCA_THRESH))
    }))
    prev_y <- mean(ye)
    dca <- rbind(
      nb,
      data.table(outcome = onm, model = "treat_all", threshold = DCA_THRESH,
                 net_benefit = prev_y - (1 - prev_y) * (DCA_THRESH / (1 - DCA_THRESH))),
      data.table(outcome = onm, model = "treat_none", threshold = DCA_THRESH, net_benefit = 0)
    )
  }
  list(summary = rows, dca = dca)
}

r_mort <- run("mortality")
r_make <- run("make")
summary_tab <- rbind(r_mort$summary, r_make$summary)
dca_tab <- rbind(r_mort$dca, r_make$dca)

## ---------- 7. readout ----------
cat("\n================ TRANSPORTABILITY HEAD-TO-HEAD (MICE labels) ================\n")
print(summary_tab)
cat("\n===== MICE-label prevalence drift =====\n")
print(prev)

m <- summary_tab[outcome == "mortality"]
bp <- m[model == "base_phen"]
rf <- m[model == "feat_rf"]
pn <- m[model == "feat_pen"]
cat("\n----- Route-1 readout (in-hospital mortality, primary) -----\n")
cat(sprintf("external AUC  : base_phen %.3f | feat_pen %.3f | feat_rf %.3f\n", bp$external_auc, pn$external_auc, rf$external_auc))
cat(sprintf("transport drop: base_phen %.3f | feat_pen %.3f | feat_rf %.3f  (smaller=more robust)\n", bp$transport_drop, pn$transport_drop, rf$transport_drop))
cat(sprintf("cal slope     : base_phen %.2f | feat_pen %.2f | feat_rf %.2f  (closer to 1=better)\n", bp$ext_cal_slope, pn$ext_cal_slope, rf$ext_cal_slope))
route1 <- (bp$transport_drop < rf$transport_drop) && (abs(bp$ext_cal_slope - 1) < abs(rf$ext_cal_slope - 1))
cat(sprintf("\n>>> Route 1 (parsimonious stratifier transports more robustly) SUPPORTED: %s\n", route1))
cat("    Claim is 'degrades less / calibrates better / fewer inputs', NOT necessarily higher external AUC.\n")

fwrite(summary_tab, file.path(OUT_DIR, "headtohead_summary.csv"))
fwrite(summary_tab, file.path(OUT_DIR, "headtohead_summary_mice_labels.csv"))
fwrite(data.table(feature_order = seq_along(common), feature = common),
       file.path(OUT_DIR, "external_assignment_features_mice_labels.csv"))
if (DO_DCA) {
  fwrite(dca_tab, file.path(OUT_DIR, "headtohead_dca.csv"))
  fwrite(dca_tab, file.path(OUT_DIR, "headtohead_dca_mice_labels.csv"))
}
saveRDS(list(
             script_version = SCRIPT_VERSION,
             assignment_run_id = assignment_run_id,
             summary = summary_tab,
             dca = dca_tab,
             common = common,
             common_feature_n = length(common),
             assignment_dimension_rule = "all 33 frozen MIMIC features retained; patient-level missing values imputed with frozen MIMIC medians",
             common_feature_hash = common_feature_hash,
             coverage = cov_e[common],
             feature_audit = feature_audit,
             preprocessing_recipe = recipe,
             preprocessing_recipe_table = recipe_table,
             preprocessing_hash = recipe_hash,
             frozen_centroids = frozen_centroids,
             centroid_hash = centroid_hash,
             input_fingerprint = input_fingerprint,
             eicu_label_file_md5 = eicu_label_hash,
             eicu_denominator = eicu_denominator_audit,
             lymphocyte_audit = list(
               verified = isTRUE(alc_result$verified),
               used_for_assignment = LYMPH_FEATURE %in% common,
               derivation = alc_result$derivation,
               source_columns = alc_result$source_columns,
               observed_n = alc_result$observed_n,
               imputed_n = alc_result$imputed_n,
               coverage = alc_result$coverage,
               median = alc_result$median,
               p99 = alc_result$p99
             ),
             label_source = "labels_primary_mice.rds",
             label_derivation = "plain K-means K=2 on MICE imputation #1",
             prevalence_drift = prev,
             high_risk_label = high_risk_label,
             label_direction_source = "MIMIC-IV in-hospital mortality only",
             eicu_reclustered = FALSE
           ),
        file.path(OUT_DIR, "headtohead_all_mice_labels.rds"))
message("Done. Wrote headtohead_summary_mice_labels.csv (+ dca) to ", OUT_DIR)
