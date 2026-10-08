## Domain 1 observable split-sample revision configuration (v3)
## Created after the 2026-07-16/17 reviewer stress test.
## This file defines a revised, observable-only rule. It does not overwrite or
## retroactively relabel the original preregistered analysis.

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- Sys.getenv(
  "MIMIC_IV_ROOT",
  unset = Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
)

D1_REVISION_ROOT <- Sys.getenv(
  "D1_REVISION_ROOT",
  unset = file.path(
    PROJECT_ROOT, "analysis_archive", "simulations",
    "Domain1_revision_v3_observable_heldout_20260717"
  )
)

D1_RUN_MODE <- tolower(Sys.getenv("D1_RUN_MODE", unset = "test"))
if (!D1_RUN_MODE %in% c("test", "smoke", "formal")) {
  stop("D1_RUN_MODE must be 'test', 'smoke', or 'formal'.", call. = FALSE)
}
D1_IS_FORMAL <- identical(D1_RUN_MODE, "formal")
D1_PARALLEL_ENABLED <- D1_RUN_MODE %in% c("smoke", "formal")

D1_OUTPUT_ROOT <- file.path(D1_REVISION_ROOT, "output", D1_RUN_MODE)
D1_LOG_DIR <- file.path(D1_REVISION_ROOT, "logs")
D1_PROVENANCE_DIR <- file.path(D1_REVISION_ROOT, "provenance")
dir.create(D1_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(D1_LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(D1_PROVENANCE_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Revised rule: observable evidence only
## ---------------------------------------------------------------------------
## The original 3-SD gate is retained only as a historical comparator. The
## revised verdict requires both:
##   1. shape evidence: BH-adjusted dip-test P < 0.05 on a discriminant axis
##      learned in training data and evaluated in independent/held-out data;
##   2. practical separation: a pooled-within Mahalanobis centroid distance
##      above an independently calibrated S1-continuum threshold.
##
## The rejected post-hoc 5% relative-effect margin remains only in the archived
## v1 prototype under history/ and is not exposed by the active v3 scripts.
## True-label/true-K ARI is recovery-only and must never enter a D1 verdict.

D1_RULE_VERSION <- "heldout_shape_plus_pipeline_matched_mahalanobis_v3_20260717"
D1_CHECKPOINT_SCHEMA_VERSION <- "d1_v3_schema_20260717_2"
D1_DIP_FDR_ALPHA <- 0.05
D1_EVALUATION_N <- 1000L
D1_TRAIN_FRACTION <- 0.50
D1_KMEANS_NSTART <- if (D1_IS_FORMAL) 50L else 10L
D1_MAHALANOBIS_RIDGE_MULTIPLIER <- 1e-4
D1_CALIBRATION_QUANTILE <- 0.95
D1_MICE_CALIBRATION_REP <- if (D1_IS_FORMAL) 200L else 2L
D1_COMPLETE_CALIBRATION_REP <- if (D1_IS_FORMAL) 200L else 5L
D1_MICE_M <- 5L
D1_MICE_MAXIT <- 5L
D1_PRIMARY_IMPUTATION_ID <- 1L
D1_KNOWN_TRUTH_TRAIN_N <- 5000L
D1_KNOWN_TRUTH_REP <- if (D1_IS_FORMAL) 100L else 2L

## Detection-limit scan. The evaluation subset remains fixed at 1,000, so
## increasing the generated cohort size does not mechanically inflate the
## Gaussian-null standardized deviation. Only D1 is run; outcomes, Domain 2-4,
## PAC, Gap, and external validation are deliberately excluded.
D1_N_GRID <- if (D1_IS_FORMAL) {
  c(2000L, 5000L, 20000L, 50000L)
} else {
  c(2000L, 5000L)
}

D1_DELTA_GRID <- if (D1_IS_FORMAL) {
  c(0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.20)
} else {
  c(0, 1.5, 3.0)
}

D1_N_REP <- if (D1_IS_FORMAL) 100L else 2L
D1_EMPIRICAL_SUBSAMPLES <- if (D1_IS_FORMAL) 100L else {
  if (identical(D1_RUN_MODE, "smoke")) 4L else 5L
}
D1_EMPIRICAL_PRIMARY_SPLIT_ID <- 1L
D1_BASE_SEED <- 2026071701L

detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores)) detected_cores <- 2L
D1_N_WORKERS <- as.integer(Sys.getenv(
  "D1_N_WORKERS",
  unset = if (D1_IS_FORMAL) {
    as.character(max(1L, min(4L, detected_cores - 2L)))
  } else if (identical(D1_RUN_MODE, "smoke")) {
    "2"
  } else {
    "1"
  }
))

## ---------------------------------------------------------------------------
## Locked inputs. Historical directories are never searched recursively.
## ---------------------------------------------------------------------------
LOCKED_ROOT <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations", "locked_results_20260712"
)

D1_INPUTS <- list(
  main_d1 = file.path(
    LOCKED_ROOT, "01_extended_main_locked", "tables",
    "domain1_discreteness_summary.csv"
  ),
  s7_d1 = file.path(
    PROJECT_ROOT, "analysis_archive", "simulations", "S7_audit_20260712", "extracted",
    "simulation_v9_scenario_S7_discrete_structure_outcome_null", "output",
    "tables", "domain1_discreteness_summary.csv"
  ),
  sa01_d1 = file.path(
    LOCKED_ROOT, "10_SA01_formal_locked_source_confirmed", "source_inputs",
    "domain1_discreteness_summary.csv"
  ),
  sa03b_diagnostics = file.path(
    LOCKED_ROOT, "04_SA03b_locked", "tables",
    "sensitivity03b_diagnostic_summary.csv"
  ),
  sa03b_recovery = file.path(
    LOCKED_ROOT, "04_SA03b_locked", "tables",
    "sensitivity03b_domain1_by_repeat.csv"
  ),
  sa03b_engine = file.path(
    LOCKED_ROOT, "04_SA03b_locked", "script",
    "03b_subtype_separation_pac_assignment_reproducibility_ARCHIVE.R"
  ),
  main_s1_s2_s7_engine = file.path(
    LOCKED_ROOT, "09_S7_discriminant_validity_locked",
    "08_s7_discrete_structure_outcome_null_v9_FORMAL.R"
  ),
  empirical_matrix = file.path(
    PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"
  ),
  empirical_labels = file.path(
    PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"
  ),
  empirical_old_dip = file.path(
    PROJECT_ROOT, "output", "qc_domain1_mice", "clusterability_dip_v2.csv"
  ),
  empirical_old_sigclust = file.path(
    PROJECT_ROOT, "output", "qc_domain1_mice",
    "clusterability_sigclust_sensitivity.csv"
  )
)

required_packages <- c(
  "data.table", "mclust", "diptest", "ggplot2", "mice", "pROC"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this package.",
    call. = FALSE
  )
}

cat(
  "Domain 1 revision v3 configuration\n",
  "  mode: ", D1_RUN_MODE, "\n",
  "  rule: ", D1_RULE_VERSION, "\n",
  "  checkpoint schema: ", D1_CHECKPOINT_SCHEMA_VERSION, "\n",
  "  dip FDR alpha: ", D1_DIP_FDR_ALPHA, "\n",
  "  primary effect: pipeline-matched independently calibrated Mahalanobis separation\n",
  "  calibration quantile: ", D1_CALIBRATION_QUANTILE, "\n",
  "  fixed evaluation n: ", D1_EVALUATION_N, "\n",
  "  MICE-gate calibration repeats: ", D1_MICE_CALIBRATION_REP, "\n",
  "  complete-data gate calibration repeats: ", D1_COMPLETE_CALIBRATION_REP, "\n",
  "  split-sample train fraction: ", D1_TRAIN_FRACTION, "\n",
  "  parallel code path: ", D1_PARALLEL_ENABLED, "\n",
  "  workers: ", D1_N_WORKERS, "\n",
  sep = ""
)
