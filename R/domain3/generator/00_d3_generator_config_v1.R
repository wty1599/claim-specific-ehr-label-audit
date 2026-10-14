## =====================================================================
## D3 regenerated deployable surrogate generator v1: locked configuration
## This file contains configuration only. It must not read external outcomes.
## =====================================================================

PROJECT_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
RESTRICTED_DATA_ROOT <- Sys.getenv(
  "EHR_AUDIT_RESTRICTED_DATA_ROOT",
  unset = file.path(PROJECT_ROOT, "data_restricted")
)
RUNTIME_EXPECTED <- "4.6.0"
if (as.character(getRversion()) != RUNTIME_EXPECTED) {
  stop(
    "D3 generator v1 requires R ", RUNTIME_EXPECTED,
    "; current runtime is ", as.character(getRversion()), "."
  )
}
MASTER_SEED <- 20260728L
BOOT_B <- 2000L
CV_REPEATS <- 10L
CV_FOLDS <- 5L
MODEL_BOOT_B <- 1000L
HOSPITAL_BOOT_B <- 1000L
MODEL_CV_FOLDS <- 10L
RF_TREES <- 300L
ENET_ALPHA <- 0.5

D3_PREVALENCE_DRIFT_GATE <- 0.10
D3_ASSIGNMENT_ARI_GATE <- 0.80
D3_CALIBRATION_SLOPE_RANGE <- c(0.80, 1.20)
D3_CALIBRATION_INTERCEPT_GATE <- 0.20

EXPECTED_MIMIC_N <- 20049L
EXPECTED_EICU_STAYS <- 17465L
EXPECTED_EICU_PATIENTS <- 16212L
EXPECTED_EICU_HOSPITALS <- 199L

EXPECTED_APPARENT_AGREEMENT_N <- 19767L
EXPECTED_APPARENT_DISCORDANT_N <- 282L
EXPECTED_APPARENT_ARI_ROUNDED <- 0.936
EXPECTED_EICU_C1_N <- 3462L
EXPECTED_GENERATOR_FILE_MD5 <- "4722d24103151ccc5bcbcc737772f249"
EXPECTED_GENERATOR_CONTENT_HASH <- "aef8e864a31e8e725a99e9e0b6b7c009"
EXPECTED_EXTERNAL_LABEL_MD5 <- "5958bcba237a61ed2689e590a572e3ca"

FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min",
  "abs_lymphocytes_min", "lactate_max", "ph_min", "pco2_max",
  "calcium_min", "calcium_max", "ptt_max", "inr_max",
  "temperature_min", "temperature_max", "urine_output_24h_ml",
  "glucose_max", "aniongap_max", "potassium_min", "potassium_max",
  "hemoglobin_min", "sodium_min", "sodium_max", "wbc_max",
  "platelets_min", "bicarbonate_min", "chloride_min", "chloride_max",
  "bun_max", "creatinine_max", "resp_rate_max", "gcs_min", "spo2_min",
  "heart_rate_max", "mbp_min"
)
stopifnot(length(FEATURES33) == 33L, !anyDuplicated(FEATURES33))

BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
SOURCE_LABEL_LEVELS <- c("1", "2")
SOURCE_HIGHER_RISK_LABEL <- "1"

PATH_FINAL_FULL <- file.path(RESTRICTED_DATA_ROOT, "final_full.csv")
PATH_EICU <- file.path(RESTRICTED_DATA_ROOT, "eicu_external.csv")
PATH_BASELINE <- file.path(RESTRICTED_DATA_ROOT, "baseline_covars.csv")
PATH_MICE <- file.path(RESTRICTED_DATA_ROOT, "model", "mice_primary.rds")
PATH_MICE_X <- file.path(
  RESTRICTED_DATA_ROOT, "model", "X_primary33_std_mice.rds"
)
PATH_LOCKED_LABELS <- file.path(
  RESTRICTED_DATA_ROOT, "model", "labels_primary_mice.rds"
)
PATH_LEGACY_EICU_LABELS <- file.path(
  RESTRICTED_DATA_ROOT, "eicu_cluster_labels_mice.csv"
)
PATH_REFERENCE_RECIPE <- file.path(
  RESTRICTED_DATA_ROOT, "qc",
  "external_assignment_preprocessing_recipe_mice_labels.csv"
)
PATH_REFERENCE_CENTROIDS <- file.path(
  RESTRICTED_DATA_ROOT, "qc",
  "external_assignment_frozen_centroids_mice_labels.csv"
)
PATH_REFERENCE_PROVENANCE <- file.path(
  RESTRICTED_DATA_ROOT, "qc",
  "external_assignment_provenance_mice_labels.csv"
)

LOCKED_MD5 <- c(
  final_full = "a3c6582cc68922787b432632f3aff8d7",
  eicu_external = "daa5589d74bdf7c09394863cec456280",
  locked_labels = "e741e9cf6bc39bf8075127a76dab5b2a",
  mice_primary = "9da5ae19f60efaa10d4b54845d897f4e",
  mice_matrix = "15a6f073e7e6d6c2d611b0b2d04b8f99",
  legacy_eicu_labels = "f208f66249d5219a1a87c5a203089c2e",
  baseline_covars = "972e5b37ab179908d926c2cd34618129",
  reference_recipe = "9943c37b255fc4f8dcc93d9e577e70c6",
  reference_centroids = "9fa3baf829958c9a0ab4ad8a0c6cc3fe"
)

SCRIPT_DIR <- file.path(PROJECT_ROOT, "R", "domain3", "generator")
OUT_ROOT <- file.path(
  PROJECT_ROOT, "outputs", "domain3_deployable_generator_v1"
)
OUT_MODEL <- file.path(OUT_ROOT, "model")
OUT_LABELS <- file.path(OUT_ROOT, "labels")
OUT_TABLES <- file.path(OUT_ROOT, "tables")
OUT_QC <- file.path(OUT_ROOT, "qc")
OUT_FIGURES <- file.path(OUT_ROOT, "figures")
OUT_LOGS <- file.path(OUT_ROOT, "logs")
OUT_RELEASE <- file.path(OUT_ROOT, "release")
OUT_PRIVATE <- file.path(OUT_ROOT, "private_local")

for (d in c(
  OUT_ROOT, OUT_MODEL, OUT_LABELS, OUT_TABLES, OUT_QC, OUT_FIGURES,
  OUT_LOGS, OUT_RELEASE, OUT_PRIVATE
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PATH_GENERATOR <- file.path(
  OUT_MODEL, "ehr_audit_k2_deployable_surrogate_generator_v1.rds"
)
PATH_EXTERNAL_LABELS <- file.path(
  OUT_LABELS, "eicu_d3_generator_v1_labels.csv"
)
PATH_EXTERNAL_FREEZE <- file.path(
  OUT_QC, "eicu_d3_generator_v1_label_fingerprint.csv"
)

GENERATOR_NAME <- "ehr_audit_k2_deployable_surrogate_generator_v1"
GENERATOR_VERSION <- "1.0.0-20260728"
GENERATOR_STATUS <- "post_hoc_regenerated_deployable_surrogate"
DISTANCE_METRIC <- "squared_euclidean"
MISSINGNESS_POLICY <- "frozen_MIMIC_post_winsor_median"

REQUIRED_INPUTS <- c(
  PATH_FINAL_FULL, PATH_EICU, PATH_BASELINE, PATH_MICE, PATH_MICE_X,
  PATH_LOCKED_LABELS, PATH_LEGACY_EICU_LABELS, PATH_REFERENCE_RECIPE,
  PATH_REFERENCE_CENTROIDS, PATH_REFERENCE_PROVENANCE
)

message("Loaded D3 generator configuration: ", GENERATOR_VERSION)
