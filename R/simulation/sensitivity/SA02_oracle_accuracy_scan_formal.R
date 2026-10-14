
################################################################################
# simulation_four_domain_framework_mice_parallel_v6_manuscript_aligned.R
#
# Semi-synthetic simulation for an EHR-derived EHR audit subphenotype validation
# framework with MICE-based imputation, formal winsorization, frozen centroids,
# and four-domain validation:
#
#   Domain 1: Discreteness / clusterability
#   Domain 2: Incremental information value
#   Domain 3: Transportability
#   Domain 4: Actionability / positivity
#
# This script is independent of real MIMIC-IV or eICU-CRD data.
#
# Major changes versus the first draft:
#   1. Primary imputation uses MICE, m = 5, maxit = 5.
#   2. Each imputed dataset is analyzed separately.
#   3. Main results are summarized as across-imputation means.
#   4. ARI, AUC, Brier, calibration slope, and positivity metrics are summarized
#      across the five imputed datasets. Rubin pooling is not forced for
#      non-regression metrics.
#   5. Formal winsorization is performed after MICE and before standardization.
#   6. Derivation-set winsor bounds, means, SDs, and centroids are exported.
#   7. External validation uses MICE for missing data but retains frozen
#      derivation preprocessing: derivation winsor bounds, derivation mean/SD,
#      and derivation centroids. External cohorts are not re-clustered.
#   8. Scenario 2 has been strengthened as a positive-control setting:
#      stronger module shifts, attenuated continuous severity dominance, and
#      K=3 true-K ARI reporting in addition to the K=2 compressed label.
#   9. Scenario 5 is added as a Domain 2 metric-level positive-control setting.
#      A hidden latent state is weakly expressed by the 33 raw features, directly
#      affects outcomes, and is represented in Domain 2 by a noisy oracle K2
#      label that is deliberately unavailable to the raw-feature model.
#  10. MICE handling is aligned with the empirical manuscript workflow:
#      m=5; completed dataset 1 defines the primary K2 label; all five completed
#      datasets are re-analysed; cross-imputation label stability is quantified
#      by ARI; Domain 2 decisive Delta AUC is Rubin-style pooled.
#  11. Scenario 5's noisy oracle label is generated once in the raw dataset and
#      remains fixed across all five completed datasets.
#  12. Domain 4 explicitly reports empty cells, sparse cells, fitting status,
#      convergence, and quasi-separation; non-identifiable interactions are not
#      silently reduced to NA.
#  13. Across-repeat flags use a confidence interval for the across-repeat mean,
#      while empirical 2.5th/97.5th percentiles are retained only as dispersion.
#  14. Scenario 4 treatment assignment is strengthened to induce near-structural
#      non-positivity in the low-risk phenotype, including empty/sparse cells
#      and truly non-identifiable interaction terms in some repeats.
#  15. safe_median() now always returns double, avoiding data.table grouped
#      aggregation failures when integer columns are all-missing in some groups
#      and non-missing in others.
################################################################################

rm(list = ls())

## =============================================================================
## 0. User settings
## =============================================================================

PROJECT_DIR <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
SIM_ROOT    <- file.path(PROJECT_DIR, "simulation_four_domain_framework_mice_sensitivity_safeplot_02_s5_oracle_accuracy_scan")

OUT_DIR <- file.path(SIM_ROOT, "output")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
LOG_DIR <- file.path(OUT_DIR, "logs")
RAW_DIR <- file.path(OUT_DIR, "raw_results")

dir.create(SIM_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RAW_DIR,  recursive = TRUE, showWarnings = FALSE)

CREATE_STANDARD_FIGURES <- TRUE

GLOBAL_SEED <- 20250101
set.seed(GLOBAL_SEED)

N_TRAIN    <- 5000
N_EXTERNAL <- 3000
P_FEATURES <- 33

N_REP      <- 100
N_REP_TEST <- 2
USE_TEST_MODE <- FALSE
N_ACTIVE_REP <- if (USE_TEST_MODE) N_REP_TEST else N_REP

## MICE is the primary imputation method. Median imputation sensitivity is
## intentionally not run in this version.
IMPUTATION_METHOD <- "mice"
MICE_M     <- 5L
MICE_MAXIT <- 5L
PRIMARY_IMPUTATION_ID <- 1L

if (MICE_M != 5L || PRIMARY_IMPUTATION_ID != 1L) {
  stop("Manuscript-aligned settings require MICE_M=5 and PRIMARY_IMPUTATION_ID=1.",
       call. = FALSE)
}

## Domain 2: within-imputation variance for paired Delta AUC is estimated by
## paired DeLong, then pooled across the five completed datasets using Rubin's
## variance decomposition. This avoids a computationally prohibitive nested
## bootstrap inside the 100-repetition simulation.
D2_RUBIN_CONF_LEVEL <- 0.95

K_MAIN <- 2
K_SENS <- 3

## Runtime controls. For quick tests, reduce GAP_BOOT and SIGCLUST_NULL_B.
K_GAP_MAX <- 10
GAP_BOOT <- 5
GAP_SUBSAMPLE_N <- 1000
SIGCLUST_NULL_B <- 10
SIGCLUST_SUBSAMPLE_N <- 1000
GLMNET_NFOLD <- 5

## Scenario 2 is a positive-control data-generating mechanism with true
## discrete latent subtypes. The following parameters intentionally strengthen
## subtype separation and reduce the dominance of the continuous severity axis.
S2_SUBTYPE_DELTA <- 3.2
S2_SEVERITY_SIGNAL_MULTIPLIER <- 0.35
S2_ORGAN_SIGNAL_MULTIPLIER <- 0.75
S2_NOISE_SD <- 0.45

## Scenario 5 is a Domain 2 metric-level positive-control mechanism. The hidden
## latent state is only weakly expressed by the 33 raw features, directly affects
## outcomes, and is represented in Domain 2 by a noisy oracle K2 label. The raw
## feature model cannot access this label or the hidden latent variables.
S5_LATENT_PREVALENCE <- 0.35
S5_FEATURE_DELTA <- 0.35
S5_LATENT_OUTCOME_LOGOR_MORT <- 1.70
S5_LATENT_OUTCOME_LOGOR_MAKE <- 1.45
S5_SEVERITY_SIGNAL_MULTIPLIER <- 0.65
S5_ORGAN_SIGNAL_MULTIPLIER <- 0.85
S5_NOISE_SD <- 0.70
S5_ORACLE_LABEL_ACCURACY <- 0.90

## Scenario 4 treatment-assignment parameters. These are intentionally extreme
## so that the low-risk phenotype has approximately zero treated patients in
## some repeats, creating empty/sparse treatment-by-phenotype cells and making
## the phenotype-by-treatment interaction genuinely non-identifiable.
S4_STRONG_GAMMA0 <- -8.50
S4_STRONG_GAMMA_SEVERITY <- 3.20
S4_STRONG_GAMMA_RENAL <- 1.40
S4_STRONG_GAMMA_SHOCK <- 1.40

WINSOR_PROBS <- c(0.01, 0.99)

SCENARIOS <- c(
  "S5_domain2_incremental_value_positive_control"
)

## =============================================================================
## 1. Packages and logging
## =============================================================================

required_pkgs <- c(
  "data.table", "ggplot2", "patchwork", "MASS", "cluster",
  "diptest", "mclust", "glmnet", "pROC", "mice", "future.apply" 
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Required packages are not installed:\n  ",
      paste(missing_pkgs, collapse = ", "),
      "\n\nInstall them with:\n  install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))\n"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(MASS)
  library(cluster)
  library(diptest)
  library(mclust)
  library(glmnet)
  library(pROC)
  library(mice)
})


options(future.globals.maxSize = 8 * 1024^3)


future::plan(future::multisession, workers = 14)

setDTthreads(percent = 50)

message("Parallel workers: ", future::nbrOfWorkers(), " (multisession)")

log_file <- file.path(LOG_DIR, paste0("simulation_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Simulation started.")


## =============================================================================
## Sensitivity-analysis grid: sensitivity_02_s5_oracle_accuracy_scan
## =============================================================================
SENSITIVITY_CATEGORY <- "sensitivity_02_s5_oracle_accuracy_scan"
SENSITIVITY_GRID <- data.table(
  sensitivity_id = paste0("oracle_acc_", gsub("\\.", "p", sprintf("%.2f", c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00)))),
  S5_ORACLE_LABEL_ACCURACY_value = c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00),
  N_TRAIN_value = 5000L,
  N_EXTERNAL_value = 3000L
)
if (!"sensitivity_id" %in% names(SENSITIVITY_GRID)) {
  stop("SENSITIVITY_GRID must contain sensitivity_id.", call. = FALSE)
}
SENSITIVITY_GRID <- data.table::as.data.table(SENSITIVITY_GRID)
SENSITIVITY_GRID[, sens_setting_index := .I]
data.table::fwrite(SENSITIVITY_GRID, file.path(TAB_DIR, "sensitivity_grid.csv"))
log_msg("Sensitivity category: ", SENSITIVITY_CATEGORY)
log_msg("Sensitivity settings: ", nrow(SENSITIVITY_GRID))

apply_sensitivity_params <- function(params) {
  ## This function is called inside each parallel worker before data generation.
  ## Assignments are worker-local under multisession and preserve the user's
  ## parallel harness.
  if ("N_TRAIN_value" %in% names(params) && is.finite(as.numeric(params$N_TRAIN_value))) {
    N_TRAIN <<- as.integer(params$N_TRAIN_value)
  }
  if ("N_EXTERNAL_value" %in% names(params) && is.finite(as.numeric(params$N_EXTERNAL_value))) {
    N_EXTERNAL <<- as.integer(params$N_EXTERNAL_value)
  }
  if ("S2_SUBTYPE_DELTA_value" %in% names(params) && is.finite(as.numeric(params$S2_SUBTYPE_DELTA_value))) {
    S2_SUBTYPE_DELTA <<- as.numeric(params$S2_SUBTYPE_DELTA_value)
  }
  if ("S5_ORACLE_LABEL_ACCURACY_value" %in% names(params) && is.finite(as.numeric(params$S5_ORACLE_LABEL_ACCURACY_value))) {
    S5_ORACLE_LABEL_ACCURACY <<- as.numeric(params$S5_ORACLE_LABEL_ACCURACY_value)
  }
  if ("SENS_MISSING_MECHANISM_value" %in% names(params) && !is.na(params$SENS_MISSING_MECHANISM_value)) {
    SENS_MISSING_MECHANISM <<- as.character(params$SENS_MISSING_MECHANISM_value)
  }
  if ("SENS_MISSING_INTENSITY_value" %in% names(params) && is.finite(as.numeric(params$SENS_MISSING_INTENSITY_value))) {
    SENS_MISSING_INTENSITY <<- as.numeric(params$SENS_MISSING_INTENSITY_value)
  }
  if ("CLUSTER_ALGORITHM_value" %in% names(params) && !is.na(params$CLUSTER_ALGORITHM_value)) {
    CLUSTER_ALGORITHM <<- as.character(params$CLUSTER_ALGORITHM_value)
  }
  if ("S6_TRUE_INTERACTION_LOGOR_value" %in% names(params) && is.finite(as.numeric(params$S6_TRUE_INTERACTION_LOGOR_value))) {
    S6_TRUE_INTERACTION_LOGOR <<- as.numeric(params$S6_TRUE_INTERACTION_LOGOR_value)
  }
  invisible(TRUE)
}

attach_sensitivity_metadata <- function(res, params) {
  meta <- as.list(params)
  drop_names <- c("scenario", "rep_id")
  meta <- meta[setdiff(names(meta), drop_names)]
  meta_names <- names(meta)
  for (nm in names(res)) {
    if (data.table::is.data.table(res[[nm]]) && nrow(res[[nm]]) > 0) {
      res[[nm]][, sensitivity_category := SENSITIVITY_CATEGORY]
      for (mn in meta_names) {
        val <- meta[[mn]]
        if (length(val) != 1L) next
        col <- if (mn == "sensitivity_id") "sensitivity_id" else paste0("sens_", mn)
        res[[nm]][, (col) := val]
      }
    }
  }
  res
}

## ---- Robust table binding -----------------------------------------------------
## Prevent final-stage failures from duplicated column names, especially generic
## columns such as k/K produced by gap statistic or K-sensitivity objects.
fix_duplicate_colnames <- function(x) {
  x <- data.table::as.data.table(x)
  if (anyDuplicated(names(x))) {
    data.table::setnames(x, make.unique(names(x), sep = "_dup"))
  }
  x
}

safe_rbindlist <- function(x, fill = TRUE, use.names = TRUE) {
  x <- Filter(function(z) {
    !is.null(z) && nrow(data.table::as.data.table(z)) > 0
  }, x)
  if (length(x) == 0) return(data.table::data.table())
  x <- lapply(x, fix_duplicate_colnames)
  data.table::rbindlist(x, fill = fill, use.names = use.names)
}

log_msg("SIM_ROOT: ", SIM_ROOT)
log_msg("N_ACTIVE_REP: ", N_ACTIVE_REP)
log_msg("IMPUTATION_METHOD: ", IMPUTATION_METHOD, "; MICE_M: ", MICE_M, "; MICE_MAXIT: ", MICE_MAXIT)

## =============================================================================
## 2. Feature metadata and general utilities
## =============================================================================

feature_names <- c(
  "aniongap_max", "creatinine_max", "bun_max", "lactate_max",
  "bicarbonate_min", "ph_min", "urine_output_24h", "mbp_min",
  "pao2fio2_min", "spo2_min", "heart_rate_max", "resp_rate_max",
  "wbc_max", "platelets_min", "bilirubin_max", "inr_max", "ptt_max",
  "potassium_max", "sodium_max", "temperature_max", "gcs_min",
  "hemoglobin_min", "glucose_max", "chloride_max", "alt_max", "ast_max",
  "albumin_min", "calcium_min", "magnesium_max", "vasopressor_dose_max",
  "shock_index_max", "baseexcess_min", "fio2_max"
)

stopifnot(length(feature_names) == P_FEATURES)

baseline_vars <- c("age", "sex", "sofa", "aki_stage")

severity_loading <- c(
  aniongap_max = 0.80, creatinine_max = 0.95, bun_max = 0.90,
  lactate_max = 0.85, bicarbonate_min = -0.85, ph_min = -0.75,
  urine_output_24h = -0.95, mbp_min = -0.70, pao2fio2_min = -0.65,
  spo2_min = -0.35, heart_rate_max = 0.55, resp_rate_max = 0.50,
  wbc_max = 0.35, platelets_min = -0.40, bilirubin_max = 0.45,
  inr_max = 0.45, ptt_max = 0.40, potassium_max = 0.35,
  sodium_max = 0.10, temperature_max = 0.25, gcs_min = -0.60,
  hemoglobin_min = -0.25, glucose_max = 0.25, chloride_max = 0.15,
  alt_max = 0.30, ast_max = 0.35, albumin_min = -0.45,
  calcium_min = -0.20, magnesium_max = 0.20, vasopressor_dose_max = 0.75,
  shock_index_max = 0.70, baseexcess_min = -0.80, fio2_max = 0.55
)

feature_scale <- data.table(
  feature = feature_names,
  mean = c(
    16, 2.0, 38, 2.8, 20, 7.31, 1200, 65, 230, 95, 115, 26, 16, 155,
    2.0, 1.4, 42, 4.6, 140, 37.8, 12, 10, 165, 104, 90, 110,
    2.8, 8.0, 2.2, 0.08, 1.0, -4.0, 0.45
  ),
  sd = c(
    5, 1.2, 25, 2.2, 5, 0.10, 850, 14, 90, 4, 22, 7, 9, 80,
    2.5, 0.5, 20, 0.8, 5, 1.0, 3, 2.0, 80, 6, 300, 500,
    0.7, 0.7, 0.5, 0.15, 0.5, 5.0, 0.25
  )
)

positive_features <- setdiff(feature_names, c("ph_min", "temperature_max", "baseexcess_min"))
bounded_low  <- c(ph_min = 6.70, spo2_min = 50, gcs_min = 3, calcium_min = 4)
bounded_high <- c(ph_min = 7.60, spo2_min = 100, gcs_min = 15, temperature_max = 43, fio2_max = 1)

safe_plogis <- function(x) plogis(pmin(pmax(x, -30), 30))

clamp_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

logit <- function(p) {
  p <- clamp_prob(p)
  log(p / (1 - p))
}

## Formal winsorization. Bounds are estimated in the derivation/training set and
## then applied to validation/test data.
winsor_fit <- function(dt, vars, probs = WINSOR_PROBS) {
  out <- safe_rbindlist(lapply(vars, function(v) {
    x <- dt[[v]]
    qs <- suppressWarnings(stats::quantile(x, probs = probs, na.rm = TRUE, type = 7))
    if (!all(is.finite(qs))) qs <- c(NA_real_, NA_real_)
    data.table(feature = v, winsor_low = as.numeric(qs[1]), winsor_high = as.numeric(qs[2]))
  }))
  out[]
}

winsor_apply <- function(dt, bounds) {
  dt <- copy(as.data.table(dt))
  for (i in seq_len(nrow(bounds))) {
    v <- bounds$feature[i]
    lo <- bounds$winsor_low[i]
    hi <- bounds$winsor_high[i]
    if (v %in% names(dt) && is.finite(lo) && is.finite(hi)) {
      dt[[v]] <- pmin(pmax(dt[[v]], lo), hi)
    }
  }
  dt[]
}

standardize_train <- function(dt, vars) {
  X <- as.matrix(dt[, ..vars])
  mu <- colMeans(X, na.rm = TRUE)
  sig <- apply(X, 2, stats::sd, na.rm = TRUE)
  sig[!is.finite(sig) | sig == 0] <- 1
  Xz <- sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
  list(Xz = Xz, center = mu, scale = sig)
}

standardize_apply <- function(dt, vars, center, scale) {
  X <- as.matrix(dt[, ..vars])
  sweep(sweep(X, 2, center, "-"), 2, scale, "/")
}

assert_no_missing_features <- function(dt, vars, where = "") {
  miss <- colSums(is.na(dt[, ..vars]))
  bad <- names(miss)[miss > 0]
  if (length(bad) > 0) {
    stop(
      paste0("Remaining missing values after MICE at ", where, ": ",
             paste(bad, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## =============================================================================
## 3. Semi-synthetic data generation
## =============================================================================

make_feature_covariance <- function(p) {
  groups <- rep(1:6, length.out = p)
  Sigma <- matrix(0.15, p, p)
  diag(Sigma) <- 1
  for (g in unique(groups)) {
    idx <- which(groups == g)
    Sigma[idx, idx] <- 0.45
    diag(Sigma)[idx] <- 1
  }
  Sigma
}

to_clinical_scale <- function(Z) {
  X <- copy(as.data.table(Z))
  setnames(X, names(X), feature_names)
  for (j in seq_along(feature_names)) {
    v <- feature_names[j]
    m <- feature_scale[feature == v, mean]
    s <- feature_scale[feature == v, sd]
    X[[v]] <- m + s * X[[v]]
    if (v %in% positive_features) X[[v]] <- pmax(X[[v]], 0.001)
    if (v %in% names(bounded_low)) X[[v]] <- pmax(X[[v]], bounded_low[[v]])
    if (v %in% names(bounded_high)) X[[v]] <- pmin(X[[v]], bounded_high[[v]])
  }
  X[]
}

SENS_MISSING_MECHANISM <- "MAR"
SENS_MISSING_INTENSITY <- 0.08
SENS_S3_VALIDATION_MISSING_MECHANISM <- "MNAR"
SENS_S3_VALIDATION_MISSING_INTENSITY <- 0.18

inject_missingness <- function(dt, vars, severity_z, mechanism = c("MCAR", "MAR", "MNAR"), intensity = 0.10) {
  mechanism <- match.arg(mechanism)
  dt <- copy(as.data.table(dt))
  intensity <- as.numeric(intensity)
  if (!is.finite(intensity) || intensity <= 0) return(dt[])
  intensity <- min(max(intensity, 0), 0.80)

  for (v in vars) {
    if (mechanism == "MCAR") {
      p_miss <- rep(intensity, nrow(dt))
    } else if (mechanism == "MAR") {
      ## In ICU EHRs, more severe patients are often measured more frequently.
      p_miss <- safe_plogis(logit(intensity) - 0.65 * severity_z)
    } else {
      ## MNAR sensitivity: missingness depends partly on latent severity and
      ## partly on the variable's own value before masking.
      z_feature <- as.numeric(scale(dt[[v]]))
      z_feature[!is.finite(z_feature)] <- 0
      p_miss <- safe_plogis(logit(intensity) + 0.35 * severity_z + 0.25 * z_feature)
    }
    miss <- stats::rbinom(nrow(dt), 1, p_miss) == 1
    dt[miss, (v) := NA_real_]
  }
  dt[]
}

simulate_latent_features <- function(
    n,
    severity_mean = 0,
    severity_sd = 1,
    true_subtype = NULL,
    subtype_delta = 0,
    external_shift = NULL,
    severity_signal_multiplier = 1,
    organ_signal_multiplier = 1,
    noise_sd = 0.65,
    incremental_group = NULL,
    incremental_delta = 0
) {
  Sigma <- make_feature_covariance(P_FEATURES)
  eps <- MASS::mvrnorm(n = n, mu = rep(0, P_FEATURES), Sigma = Sigma)

  severity_z <- as.numeric(stats::rnorm(n, mean = severity_mean, sd = severity_sd))

  ## The severity component in organ modules can be attenuated in Scenario 2
  ## so that the true discrete subtype structure is not overwhelmed by a
  ## single continuous severity axis.
  renal_factor  <- severity_signal_multiplier * 0.65 * severity_z + stats::rnorm(n, 0, 0.75)
  acid_factor   <- severity_signal_multiplier * 0.60 * severity_z + stats::rnorm(n, 0, 0.75)
  shock_factor  <- severity_signal_multiplier * 0.70 * severity_z + stats::rnorm(n, 0, 0.75)
  resp_factor   <- severity_signal_multiplier * 0.55 * severity_z + stats::rnorm(n, 0, 0.80)
  inflam_factor <- severity_signal_multiplier * 0.40 * severity_z + stats::rnorm(n, 0, 0.90)

  organ_load <- matrix(0, nrow = P_FEATURES, ncol = 5)
  rownames(organ_load) <- feature_names
  colnames(organ_load) <- c("renal", "acid", "shock", "resp", "inflam")

  organ_load[c("creatinine_max", "bun_max", "urine_output_24h", "potassium_max"), "renal"] <- c(0.8, 0.7, -0.8, 0.4)
  organ_load[c("lactate_max", "bicarbonate_min", "ph_min", "aniongap_max", "baseexcess_min"), "acid"] <- c(0.7, -0.7, -0.6, 0.6, -0.7)
  organ_load[c("mbp_min", "vasopressor_dose_max", "shock_index_max", "heart_rate_max"), "shock"] <- c(-0.7, 0.8, 0.8, 0.5)
  organ_load[c("pao2fio2_min", "spo2_min", "fio2_max", "resp_rate_max"), "resp"] <- c(-0.8, -0.4, 0.7, 0.4)
  organ_load[c("wbc_max", "platelets_min", "bilirubin_max", "inr_max", "ptt_max", "temperature_max"), "inflam"] <- c(0.5, -0.4, 0.4, 0.4, 0.3, 0.3)

  latent_mat <- cbind(renal_factor, acid_factor, shock_factor, resp_factor, inflam_factor)
  Z <- matrix(0, nrow = n, ncol = P_FEATURES)
  colnames(Z) <- feature_names

  for (j in seq_len(P_FEATURES)) {
    v <- feature_names[j]
    Z[, j] <- severity_signal_multiplier * severity_loading[[v]] * severity_z +
      organ_signal_multiplier * as.numeric(latent_mat %*% organ_load[v, ]) +
      noise_sd * eps[, j]
  }

  if (!is.null(true_subtype) && subtype_delta > 0) {
    G <- as.character(true_subtype)
    metabolic <- G == "metabolic"
    renal <- G == "renal_dysfunction"
    shock <- G == "shock"

    shift <- matrix(0, nrow = n, ncol = P_FEATURES)
    colnames(shift) <- feature_names

    ## Module-specific subtype shifts are intentionally strong in this positive
    ## control scenario. They create separable metabolic, renal, and shock-like
    ## groups while retaining within-subtype heterogeneity.
    shift[metabolic, c("lactate_max", "aniongap_max", "bicarbonate_min", "ph_min", "baseexcess_min",
                       "glucose_max", "chloride_max")] <-
      matrix(rep(c(1.25, 1.10, -1.15, -0.95, -1.10, 0.65, 0.45), sum(metabolic)), ncol = 7, byrow = TRUE)

    shift[renal, c("creatinine_max", "bun_max", "urine_output_24h", "potassium_max",
                   "ph_min", "bicarbonate_min", "aniongap_max")] <-
      matrix(rep(c(1.35, 1.35, -1.35, 0.75, -0.35, -0.45, 0.45), sum(renal)), ncol = 7, byrow = TRUE)

    shift[shock, c("mbp_min", "vasopressor_dose_max", "shock_index_max", "lactate_max", "gcs_min",
                   "pao2fio2_min", "fio2_max")] <-
      matrix(rep(c(-1.25, 1.35, 1.35, 0.85, -0.75, -0.55, 0.65), sum(shock)), ncol = 7, byrow = TRUE)

    Z <- Z + subtype_delta * shift
  }

  if (!is.null(incremental_group) && incremental_delta > 0) {
    ## Domain 2 positive control: a hidden state produces only a weak,
    ## coordinated multivariable pattern. The raw-feature model sees this noisy
    ## partial proxy, while Domain 2 additionally receives a noisy oracle K2
    ## label derived from the hidden state.
    G2 <- as.integer(incremental_group)
    high <- G2 == 1

    inc_shift <- matrix(0, nrow = n, ncol = P_FEATURES)
    colnames(inc_shift) <- feature_names

    inc_features <- c(
      "lactate_max", "aniongap_max", "bicarbonate_min", "ph_min",
      "mbp_min", "shock_index_max", "heart_rate_max", "resp_rate_max",
      "albumin_min", "platelets_min", "baseexcess_min", "fio2_max"
    )
    inc_values <- c(
      0.65, 0.55, -0.50, -0.40,
      -0.45, 0.50, 0.35, 0.35,
      -0.35, -0.30, -0.45, 0.35
    )

    inc_shift[high, inc_features] <-
      matrix(rep(inc_values, sum(high)), ncol = length(inc_features), byrow = TRUE)

    ## Low-risk latent patients are not simply the opposite phenotype; this
    ## asymmetry prevents the raw-feature linear model from perfectly replacing
    ## the discretized label.
    inc_shift[!high, c("lactate_max", "aniongap_max", "mbp_min", "albumin_min")] <-
      matrix(rep(c(-0.15, -0.10, 0.10, 0.10), sum(!high)), ncol = 4, byrow = TRUE)

    Z <- Z + incremental_delta * inc_shift
  }

  if (!is.null(external_shift)) {
    for (v in names(external_shift)) {
      if (v %in% colnames(Z)) Z[, v] <- Z[, v] + external_shift[[v]]
    }
    for (v in c("creatinine_max", "lactate_max", "bicarbonate_min", "mbp_min")) {
      Z[, v] <- Z[, v] + stats::rnorm(n, 0, 0.35)
    }
  }

  list(
    Z = Z,
    severity_z = severity_z,
    renal_factor = renal_factor,
    acid_factor = acid_factor,
    shock_factor = shock_factor
  )
}

generate_baseline_covariates <- function(severity_z) {
  n <- length(severity_z)
  age <- pmin(pmax(round(stats::rnorm(n, mean = 64 + 4 * severity_z, sd = 13)), 18), 95)
  sex <- stats::rbinom(n, 1, 0.55)
  sofa <- pmin(pmax(round(7 + 3.0 * severity_z + stats::rnorm(n, 0, 2.0)), 0), 24)
  p_stage3 <- safe_plogis(-1.2 + 1.25 * severity_z)
  p_stage2 <- safe_plogis(-0.5 + 0.75 * severity_z) * (1 - p_stage3)
  u <- stats::runif(n)
  aki_stage <- ifelse(u < p_stage3, 3, ifelse(u < p_stage3 + p_stage2, 2, 1))
  data.table(age = age, sex = sex, sofa = sofa, aki_stage = aki_stage)
}

generate_treatment <- function(severity_z, renal_factor, shock_factor,
                               strength = c("moderate", "strong")) {
  strength <- match.arg(strength)

  if (strength == "moderate") {
    gamma0 <- -3.2
    gamma_severity <- 1.35
    gamma_renal <- 0.65
    gamma_shock <- 0.65
  } else {
    ## Scenario 4 deliberately creates practical/near-structural non-positivity:
    ## treatment is almost never observed in low-severity / low-risk regions but
    ## becomes common among the sickest patients. This should drive the low-risk
    ## phenotype treated cell toward zero in some repeats.
    gamma0 <- S4_STRONG_GAMMA0
    gamma_severity <- S4_STRONG_GAMMA_SEVERITY
    gamma_renal <- S4_STRONG_GAMMA_RENAL
    gamma_shock <- S4_STRONG_GAMMA_SHOCK
  }

  p_t <- safe_plogis(
    gamma0 +
      gamma_severity * severity_z +
      gamma_renal * renal_factor +
      gamma_shock * shock_factor
  )
  stats::rbinom(length(severity_z), 1, p_t)
}

generate_outcomes <- function(severity_z, baseline_dt, true_subtype = NULL,
                              scenario, treatment = NULL,
                              domain2_latent_group = NULL,
                              domain2_latent_z = NULL) {
  n <- length(severity_z)
  subtype_effect <- rep(0, n)
  if (!is.null(true_subtype)) {
    ## Scenario 2 is intentionally a positive-control setting in which
    ## latent subtype membership retains residual prognostic signal beyond
    ## measured physiology. This makes Domain 2 capable of showing a positive
    ## incremental-value control.
    subtype_effect[true_subtype == "renal_dysfunction"] <- 0.65
    subtype_effect[true_subtype == "metabolic"] <- 0.45
    subtype_effect[true_subtype == "shock"] <- 0.95
  }

  lp_mort <- -2.35 +
    0.90 * severity_z +
    0.018 * (baseline_dt$age - 65) +
    0.06 * (baseline_dt$sofa - 7) +
    0.18 * (baseline_dt$aki_stage - 1)

  if (scenario == "S2_true_discrete_subtypes") lp_mort <- lp_mort + subtype_effect

  if (scenario == "S5_domain2_incremental_value_positive_control" &&
      !is.null(domain2_latent_group)) {
    ## Direct hidden-state outcome effect. This latent state is excluded from
    ## MICE, clustering, and all ordinary prediction models; the Scenario 5 K2
    ## label is a noisy proxy derived from the hidden latent state itself.
    lp_mort <- lp_mort +
      S5_LATENT_OUTCOME_LOGOR_MORT * as.integer(domain2_latent_group == 1)
    if (!is.null(domain2_latent_z)) {
      lp_mort <- lp_mort + 0.20 * pmax(domain2_latent_z, 0)
    }
  }

  if (!is.null(treatment)) {
    ## Weak average treatment effect. This is not intended to create a clean
    ## causal estimand; Scenario 4 focuses on non-positivity.
    lp_mort <- lp_mort - 0.20 * treatment
    if (!is.null(true_subtype)) lp_mort <- lp_mort + 0.15 * treatment * (true_subtype == "shock")
  }

  mortality_30d <- stats::rbinom(n, 1, safe_plogis(lp_mort))

  lp_make <- -1.10 +
    1.05 * severity_z +
    0.04 * (baseline_dt$sofa - 7) +
    0.45 * (baseline_dt$aki_stage - 1) +
    0.15 * subtype_effect

  if (scenario == "S5_domain2_incremental_value_positive_control" &&
      !is.null(domain2_latent_group)) {
    lp_make <- lp_make +
      S5_LATENT_OUTCOME_LOGOR_MAKE * as.integer(domain2_latent_group == 1)
    if (!is.null(domain2_latent_z)) {
      lp_make <- lp_make + 0.15 * pmax(domain2_latent_z, 0)
    }
  }

  make30 <- stats::rbinom(n, 1, safe_plogis(lp_make))

  data.table(mortality_30d = mortality_30d, make30 = make30)
}

generate_dataset <- function(scenario, cohort = c("derivation", "validation"),
                             n, repeat_id, seed, subtype_delta = S2_SUBTYPE_DELTA) {
  cohort <- match.arg(cohort)
  set.seed(seed)

  true_subtype <- NULL
  external_shift <- NULL
  severity_mean <- 0
  severity_sd <- 1
  domain2_latent_z <- rep(NA_real_, n)
  domain2_latent_group <- rep(0L, n)
  domain2_oracle_label <- factor(rep(NA_character_, n), levels = c("C1", "C2"))

  if (scenario == "S1_pure_severity_continuum") {
    true_subtype <- rep("continuum", n)
  }

  if (scenario == "S2_true_discrete_subtypes") {
    true_subtype <- sample(
      c("metabolic", "renal_dysfunction", "shock"),
      size = n, replace = TRUE, prob = c(0.35, 0.35, 0.30)
    )
  }

  if (scenario == "S3_transport_drift") {
    true_subtype <- rep("continuum", n)
    if (cohort == "validation") {
      severity_mean <- 0.55
      severity_sd <- 1.10
      external_shift <- c(
        creatinine_max = 0.45,
        lactate_max = 0.40,
        bicarbonate_min = -0.35,
        mbp_min = -0.35,
        urine_output_24h = -0.25,
        pao2fio2_min = -0.25
      )
    }
  }

  if (scenario == "S4_actionability_positivity_failure") {
    true_subtype <- rep("continuum", n)
  }

  if (scenario == "S5_domain2_incremental_value_positive_control") {
    ## Hidden latent state for the Domain 2 metric-sensitivity positive control.
    ## The noisy oracle label is generated once here from the hidden latent
    ## state and is then carried unchanged through all MICE completed datasets.
    domain2_latent_z <- stats::rnorm(n)
    threshold <- stats::qnorm(1 - S5_LATENT_PREVALENCE)
    domain2_latent_group <- as.integer(domain2_latent_z > threshold)
    true_subtype <- ifelse(domain2_latent_group == 1, "incremental_high", "incremental_low")
    domain2_oracle_label <- make_noisy_oracle_k2_label(
      domain2_latent_group,
      accuracy = S5_ORACLE_LABEL_ACCURACY,
      seed = seed + 909L
    )
  }

  latent <- simulate_latent_features(
    n = n,
    severity_mean = severity_mean,
    severity_sd = severity_sd,
    true_subtype = true_subtype,
    subtype_delta = ifelse(scenario == "S2_true_discrete_subtypes", subtype_delta, 0),
    external_shift = external_shift,
    severity_signal_multiplier = ifelse(
      scenario == "S2_true_discrete_subtypes",
      S2_SEVERITY_SIGNAL_MULTIPLIER,
      ifelse(
        scenario == "S5_domain2_incremental_value_positive_control",
        S5_SEVERITY_SIGNAL_MULTIPLIER,
        1
      )
    ),
    organ_signal_multiplier = ifelse(
      scenario == "S2_true_discrete_subtypes",
      S2_ORGAN_SIGNAL_MULTIPLIER,
      ifelse(
        scenario == "S5_domain2_incremental_value_positive_control",
        S5_ORGAN_SIGNAL_MULTIPLIER,
        1
      )
    ),
    noise_sd = ifelse(
      scenario == "S2_true_discrete_subtypes",
      S2_NOISE_SD,
      ifelse(
        scenario == "S5_domain2_incremental_value_positive_control",
        S5_NOISE_SD,
        0.65
      )
    ),
    incremental_group = if (scenario == "S5_domain2_incremental_value_positive_control") domain2_latent_group else NULL,
    incremental_delta = ifelse(
      scenario == "S5_domain2_incremental_value_positive_control",
      S5_FEATURE_DELTA,
      0
    )
  )

  features_dt <- to_clinical_scale(latent$Z)

  ## MICE is the imputation method. Here we create EHR-like missingness.
  ## Sensitivity scripts can override the mechanism/intensity through
  ## SENS_MISSING_* parameters without changing the main data-generating code.
  if (scenario == "S3_transport_drift" && cohort == "validation") {
    features_dt <- inject_missingness(
      features_dt, feature_names, latent$severity_z,
      mechanism = SENS_S3_VALIDATION_MISSING_MECHANISM,
      intensity = SENS_S3_VALIDATION_MISSING_INTENSITY
    )
  } else {
    features_dt <- inject_missingness(
      features_dt, feature_names, latent$severity_z,
      mechanism = SENS_MISSING_MECHANISM,
      intensity = SENS_MISSING_INTENSITY
    )
  }

  baseline_dt <- generate_baseline_covariates(latent$severity_z)

  treatment_strength <- ifelse(scenario == "S4_actionability_positivity_failure", "strong", "moderate")
  treatment <- generate_treatment(latent$severity_z, latent$renal_factor, latent$shock_factor, treatment_strength)

  outcomes_dt <- generate_outcomes(
    severity_z = latent$severity_z,
    baseline_dt = baseline_dt,
    true_subtype = true_subtype,
    scenario = scenario,
    treatment = if (scenario == "S4_actionability_positivity_failure") treatment else NULL,
    domain2_latent_group = if (scenario == "S5_domain2_incremental_value_positive_control") domain2_latent_group else NULL,
    domain2_latent_z = if (scenario == "S5_domain2_incremental_value_positive_control") domain2_latent_z else NULL
  )

  dt <- data.table(
    patient_id = sprintf("%s_rep%03d_%05d", ifelse(cohort == "derivation", "D", "V"), repeat_id, seq_len(n)),
    cohort = cohort,
    true_scenario = scenario,
    repeat_id = repeat_id,
    severity_z = latent$severity_z,
    true_subtype = true_subtype,
    domain2_latent_group = domain2_latent_group,
    domain2_latent_z = domain2_latent_z,
    domain2_oracle_label = domain2_oracle_label
  )

  dt <- cbind(dt, features_dt, baseline_dt, outcomes_dt)
  dt[, treatment := treatment]
  dt[, rrt_early := treatment]
  dt[]
}

## =============================================================================
## 4. MICE imputation
## =============================================================================

make_mice_imputations <- function(dt, scenario, repeat_id, cohort,
                                  m = MICE_M, maxit = MICE_MAXIT, seed = GLOBAL_SEED) {
  if (m != 5L) stop("This manuscript-aligned simulation requires MICE m=5.", call. = FALSE)
  dt <- copy(as.data.table(dt))

  ## Exclude all outcome, treatment, derived-label, and latent-truth variables.
  leak_vars <- c(
    "patient_id", "cohort", "true_scenario", "repeat_id",
    "severity_z", "true_subtype",
    "domain2_latent_group", "domain2_latent_z", "domain2_oracle_label",
    "mortality_30d", "make30",
    "treatment", "rrt_early",
    "cluster"
  )

  use_vars <- setdiff(intersect(c(feature_names, baseline_vars), names(dt)), leak_vars)
  imp_data <- as.data.frame(dt[, ..use_vars])

  ## Baseline covariates are fully observed by design; features have missingness.
  ## MICE sees only pre-outcome variables and does not see treatment/outcome/labels.
  ini <- tryCatch(mice::mice(imp_data, maxit = 0, printFlag = FALSE), error = function(e) NULL)
  if (is.null(ini)) {
    stop("MICE initialization failed for ", scenario, " repeat ", repeat_id, " cohort ", cohort, call. = FALSE)
  }

  meth <- ini$method
  pred <- ini$predictorMatrix

  ## Do not impute fully observed baseline variables unless missingness occurs.
  ## Keep defaults from mice, but prevent self-prediction is already handled.
  set.seed(seed)
  imp <- tryCatch(
    mice::mice(
      imp_data,
      m = m,
      maxit = maxit,
      method = meth,
      predictorMatrix = pred,
      seed = seed,
      printFlag = FALSE
    ),
    error = function(e) {
      stop("MICE failed for ", scenario, " repeat ", repeat_id, " cohort ", cohort,
           ": ", conditionMessage(e), call. = FALSE)
    }
  )

  out_list <- vector("list", m)
  for (ii in seq_len(m)) {
    completed <- as.data.table(mice::complete(imp, action = ii))
    out <- data.table::setalloccol(data.table::copy(dt))
    ## Replace only variables used by MICE. Outcomes, treatment, true labels, and
    ## severity remain untouched.
    for (v in use_vars) data.table::set(out, j = v, value = completed[[v]])
    data.table::set(out, j = "imputation_method", value = "mice")
    data.table::set(out, j = "imputation_id", value = ii)
    data.table::set(out, j = "is_primary_imputation", value = (ii == PRIMARY_IMPUTATION_ID))
    out_list[[ii]] <- out
  }

  out_list
}

## =============================================================================
## 5. Preprocessing, clustering, and centroid export
## =============================================================================

fit_kmeans_centroids <- function(dt, vars, k = 2,
                                 analysis_context = "domain",
                                 scenario = NA_character_,
                                 repeat_id = NA_integer_,
                                 imputation_id = NA_integer_) {
  dt <- copy(as.data.table(dt))
  assert_no_missing_features(dt, vars, where = paste(analysis_context, scenario, repeat_id, imputation_id))

  wb <- winsor_fit(dt, vars)
  dt_w <- winsor_apply(dt, wb)
  std <- standardize_train(dt_w, vars)

  km <- stats::kmeans(std$Xz, centers = k, nstart = 50, iter.max = 100)

  if (k == 2) {
    labs <- as.character(km$cluster)
    mean_sev <- tapply(dt$severity_z, labs, mean, na.rm = TRUE)
    high_old <- names(which.max(mean_sev))
    low_old <- setdiff(sort(unique(labs)), high_old)[1]
    label <- factor(ifelse(labs == high_old, "C1", "C2"), levels = c("C1", "C2"))
    centroids <- rbind(
      C1 = km$centers[high_old, ],
      C2 = km$centers[low_old, ]
    )
  } else {
    label <- factor(paste0("C", km$cluster))
    centroids <- km$centers
    rownames(centroids) <- paste0("C", seq_len(k))
  }

  centroids_dt <- as.data.table(centroids, keep.rownames = "cluster")
  centroids_dt[, `:=`(
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_method = "mice",
    imputation_id = imputation_id,
    analysis_context = analysis_context,
    cluster_k = k
  )]

  prep_dt <- safe_rbindlist(list(
    data.table(feature = names(std$center), parameter = "center", value = as.numeric(std$center)),
    data.table(feature = names(std$scale),  parameter = "scale",  value = as.numeric(std$scale)),
    wb[, .(feature, parameter = "winsor_low",  value = winsor_low)],
    wb[, .(feature, parameter = "winsor_high", value = winsor_high)]
  ), fill = TRUE)
  prep_dt[, `:=`(
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_method = "mice",
    imputation_id = imputation_id,
    analysis_context = analysis_context,
    cluster_k = k
  )]

  list(
    label = label,
    kmeans = km,
    centroids = centroids,
    center = std$center,
    scale = std$scale,
    winsor_bounds = wb,
    Xz = std$Xz,
    dt_winsorized = dt_w,
    centroids_dt = centroids_dt,
    preprocess_dt = prep_dt
  )
}

apply_frozen_centroids <- function(dt, vars, centroids, center, scale, winsor_bounds,
                                   where = "external") {
  dt <- copy(as.data.table(dt))
  assert_no_missing_features(dt, vars, where = where)
  dt_w <- winsor_apply(dt, winsor_bounds)
  Xz <- standardize_apply(dt_w, vars, center, scale)

  dmat <- sapply(seq_len(nrow(centroids)), function(i) {
    rowSums((sweep(Xz, 2, centroids[i, ], "-"))^2)
  })
  idx <- max.col(-dmat)
  labs <- rownames(centroids)[idx]
  list(label = factor(labs, levels = rownames(centroids)), Xz = Xz, dt_winsorized = dt_w)
}

## =============================================================================
## 6. Domain 1: discreteness diagnostics
## =============================================================================

cluster_index <- function(Xz, labels) {
  labs <- factor(labels)
  tot <- sum(rowSums((sweep(Xz, 2, colMeans(Xz), "-"))^2))
  within <- 0
  for (lv in levels(labs)) {
    Xi <- Xz[labs == lv, , drop = FALSE]
    if (nrow(Xi) > 0) within <- within + sum(rowSums((sweep(Xi, 2, colMeans(Xi), "-"))^2))
  }
  within / tot
}

run_gap_statistic <- function(Xz, max_k = K_GAP_MAX, B = GAP_BOOT, subsample_n = GAP_SUBSAMPLE_N) {
  Xuse <- Xz
  if (nrow(Xuse) > subsample_n) Xuse <- Xuse[sample.int(nrow(Xuse), subsample_n), , drop = FALSE]

  out <- tryCatch({
    gap <- cluster::clusGap(
      Xuse,
      FUNcluster = function(x, k) stats::kmeans(x, centers = k, nstart = 20, iter.max = 50),
      K.max = max_k,
      B = B,
      verbose = FALSE
    )
    tab <- as.data.table(gap$Tab)
    tab[, gap_k := seq_len(.N)]
    data.table::setcolorder(tab, "gap_k")
    list(
      table = tab,
      best_k_maxgap = tab[which.max(tab$gap), gap_k],
      best_k_firstse = cluster::maxSE(tab$gap, tab$SE.sim, method = "firstSEmax")
    )
  }, error = function(e) {
    list(table = data.table(), best_k_maxgap = NA_integer_, best_k_firstse = NA_integer_)
  })
  out
}

run_discriminant_axis <- function(Xz, label) {
  lab <- factor(label)
  if (nlevels(lab) != 2) return(rep(NA_real_, nrow(Xz)))
  mu1 <- colMeans(Xz[lab == levels(lab)[1], , drop = FALSE])
  mu2 <- colMeans(Xz[lab == levels(lab)[2], , drop = FALSE])
  w <- mu1 - mu2
  denom <- sqrt(sum(w^2))
  if (!is.finite(denom) || denom == 0) return(rep(NA_real_, nrow(Xz)))
  as.numeric(Xz %*% w / denom)
}

run_dip_tests <- function(Xz, label) {
  axis <- run_discriminant_axis(Xz, label)
  disc_p <- tryCatch(diptest::dip.test(axis)$p.value, error = function(e) NA_real_)

  pca <- tryCatch(stats::prcomp(Xz, center = FALSE, scale. = FALSE), error = function(e) NULL)
  pc_p <- rep(NA_real_, 5)
  if (!is.null(pca)) {
    for (i in seq_len(min(5, ncol(pca$x)))) {
      pc_p[i] <- tryCatch(diptest::dip.test(pca$x[, i])$p.value, error = function(e) NA_real_)
    }
  }

  data.table(
    dip_discriminant_p = disc_p,
    dip_pc1_p = pc_p[1], dip_pc2_p = pc_p[2], dip_pc3_p = pc_p[3],
    dip_pc4_p = pc_p[4], dip_pc5_p = pc_p[5]
  )
}

run_sigclust_like <- function(Xz, k = 2, B = SIGCLUST_NULL_B, subsample_n = SIGCLUST_SUBSAMPLE_N) {
  Xuse <- Xz
  if (nrow(Xuse) > subsample_n) Xuse <- Xuse[sample.int(nrow(Xuse), subsample_n), , drop = FALSE]

  obs_km <- tryCatch(stats::kmeans(Xuse, centers = k, nstart = 30, iter.max = 100), error = function(e) NULL)
  if (is.null(obs_km)) {
    return(data.table(
      sigclust_obs_index = NA_real_, sigclust_null_mean = NA_real_,
      sigclust_p = NA_real_, sigclust_effect_sd = NA_real_,
      sigclust_relative_effect_pct = NA_real_
    ))
  }

  obs_idx <- cluster_index(Xuse, obs_km$cluster)
  mu <- colMeans(Xuse)
  Sigma <- stats::cov(Xuse)
  diag(Sigma) <- diag(Sigma) + 1e-6

  null_idx <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    Xnull <- tryCatch(MASS::mvrnorm(nrow(Xuse), mu = mu, Sigma = Sigma), error = function(e) NULL)
    if (!is.null(Xnull)) {
      km_null <- tryCatch(stats::kmeans(Xnull, centers = k, nstart = 10, iter.max = 50), error = function(e) NULL)
      if (!is.null(km_null)) null_idx[b] <- cluster_index(Xnull, km_null$cluster)
    }
  }

  null_idx <- null_idx[is.finite(null_idx)]
  if (!length(null_idx)) {
    return(data.table(
      sigclust_obs_index = obs_idx, sigclust_null_mean = NA_real_,
      sigclust_p = NA_real_, sigclust_effect_sd = NA_real_,
      sigclust_relative_effect_pct = NA_real_
    ))
  }

  p_val <- mean(null_idx <= obs_idx)
  eff_sd <- (mean(null_idx) - obs_idx) / stats::sd(null_idx)
  rel_pct <- 100 * (mean(null_idx) - obs_idx) / mean(null_idx)

  data.table(
    sigclust_obs_index = obs_idx,
    sigclust_null_mean = mean(null_idx),
    sigclust_p = p_val,
    sigclust_effect_sd = eff_sd,
    sigclust_relative_effect_pct = rel_pct
  )
}

run_domain1 <- function(dt, scenario, repeat_id, imputation_id) {
  fit <- fit_kmeans_centroids(
    dt, feature_names, k = 2, analysis_context = "domain1_derivation",
    scenario = scenario, repeat_id = repeat_id, imputation_id = imputation_id
  )
  label <- fit$label
  Xz <- fit$Xz

  x <- copy(dt)
  x[, cluster := label]
  prev <- x[, .(
    n = .N,
    prevalence = .N / nrow(x),
    mortality_rate = mean(mortality_30d),
    make30_rate = mean(make30),
    mean_severity = mean(severity_z)
  ), by = cluster]

  c1_prev <- prev[cluster == "C1", prevalence]
  c1_mort <- prev[cluster == "C1", mortality_rate]
  c2_mort <- prev[cluster == "C2", mortality_rate]
  outcome_sep <- c1_mort - c2_mort

  gap <- run_gap_statistic(Xz)
  dip <- run_dip_tests(Xz, label)
  sig <- run_sigclust_like(Xz, k = 2)

  ari_k2 <- NA_real_
  ari_k3 <- NA_real_
  k3_n_clusters <- NA_integer_
  fit3 <- NULL

  if (scenario == "S2_true_discrete_subtypes") {
    ## The empirical EHR audit workflow keeps K=2 as the main analysis. However,
    ## this positive-control simulation is generated from three known latent
    ## subtypes, so K=3 recovery is explicitly reported as the "true-K"
    ## diagnostic.
    ari_k2 <- tryCatch(
      mclust::adjustedRandIndex(as.character(label), dt$true_subtype),
      error = function(e) NA_real_
    )

    fit3 <- tryCatch(
      fit_kmeans_centroids(
        dt, feature_names, k = 3, analysis_context = "domain1_derivation_k3_positive_control",
        scenario = scenario, repeat_id = repeat_id, imputation_id = imputation_id
      ),
      error = function(e) NULL
    )

    if (!is.null(fit3)) {
      ari_k3 <- tryCatch(
        mclust::adjustedRandIndex(as.character(fit3$label), dt$true_subtype),
        error = function(e) NA_real_
      )
      k3_n_clusters <- length(unique(fit3$label))
    }
  }

  summary <- data.table(
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_method = "mice",
    imputation_id = imputation_id,
    is_primary_imputation = (imputation_id == PRIMARY_IMPUTATION_ID),
    domain = "Domain1_Discreteness",
    cluster_k = 2,
    c1_prevalence = c1_prev,
    mortality_rate_C1 = c1_mort,
    mortality_rate_C2 = c2_mort,
    outcome_separation = outcome_sep,
    gap_best_k_maxgap = gap$best_k_maxgap,
    gap_best_k_firstse = gap$best_k_firstse,
    dip_discriminant_p = dip$dip_discriminant_p,
    dip_pc1_p = dip$dip_pc1_p,
    dip_pc2_p = dip$dip_pc2_p,
    dip_pc3_p = dip$dip_pc3_p,
    dip_pc4_p = dip$dip_pc4_p,
    dip_pc5_p = dip$dip_pc5_p,
    sigclust_obs_index = sig$sigclust_obs_index,
    sigclust_null_mean = sig$sigclust_null_mean,
    sigclust_p = sig$sigclust_p,
    sigclust_effect_sd = sig$sigclust_effect_sd,
    sigclust_relative_effect_pct = sig$sigclust_relative_effect_pct,
    ari_to_true_subtype = ari_k2,
    ari_k3_to_true_subtype = ari_k3,
    k3_n_clusters = k3_n_clusters
  )

  summary[, interpretation_flag := fifelse(
    scenario == "S5_domain2_incremental_value_positive_control",
    "not_domain1_positive_control_metric_sensitivity_only",
    fifelse(
      scenario == "S1_pure_severity_continuum" &
        mortality_rate_C1 > mortality_rate_C2 &
        (is.na(dip_discriminant_p) | dip_discriminant_p > 0.05) &
        (is.na(sigclust_effect_sd) | sigclust_effect_sd < 3),
      "stable_partition_not_discrete",
      fifelse(
        scenario == "S2_true_discrete_subtypes" &
          !is.na(ari_k3_to_true_subtype) & ari_k3_to_true_subtype > 0.5,
        "true_subtype_recovery_by_true_K",
        fifelse(
          scenario == "S2_true_discrete_subtypes" &
            !is.na(ari_to_true_subtype) & ari_to_true_subtype > 0.5,
          "true_subtype_recovery_by_K2",
          "mixed_or_indeterminate"
        )
      )
    )
  )]

  if (nrow(gap$table)) {
    gap$table[, `:=`(
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_method = "mice",
      imputation_id = imputation_id
    )]
  }

  centroids_out <- fit$centroids_dt
  preprocess_out <- fit$preprocess_dt
  if (!is.null(fit3)) {
    centroids_out <- safe_rbindlist(list(centroids_out, fit3$centroids_dt), fill = TRUE)
    preprocess_out <- safe_rbindlist(list(preprocess_out, fit3$preprocess_dt), fill = TRUE)
  }

  labels_dt <- data.table(
    patient_id = dt$patient_id,
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_id = imputation_id,
    is_primary_imputation = (imputation_id == PRIMARY_IMPUTATION_ID),
    cluster = as.character(label)
  )

  list(
    summary = summary,
    gap_table = gap$table,
    labels_dt = labels_dt,
    centroids_dt = centroids_out,
    preprocess_dt = preprocess_out
  )
}

## =============================================================================
## 7. Prediction metrics and Domain 2/3
## =============================================================================

calc_auc <- function(y, pred) {
  y <- as.integer(y)
  pred <- as.numeric(pred)
  ok <- is.finite(pred) & !is.na(y)
  if (length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(as.numeric(pROC::auc(pROC::roc(y[ok], pred[ok], quiet = TRUE))))
}

calc_brier <- function(y, pred) mean((as.numeric(y) - as.numeric(pred))^2, na.rm = TRUE)

calc_calibration_slope <- function(y, pred) {
  pred <- clamp_prob(pred)
  lp <- logit(pred)
  tryCatch({
    fit <- stats::glm(y ~ lp, family = stats::binomial())
    as.numeric(stats::coef(fit)[["lp"]])
  }, error = function(e) NA_real_)
}

calc_calibration_intercept <- function(y, pred) {
  pred <- clamp_prob(pred)
  lp <- logit(pred)
  tryCatch({
    fit <- stats::glm(y ~ 1, offset = lp, family = stats::binomial())
    as.numeric(stats::coef(fit)[1])
  }, error = function(e) NA_real_)
}

evaluate_predictions <- function(y, pred) {
  data.table(
    AUC = calc_auc(y, pred),
    Brier = calc_brier(y, pred),
    calibration_slope = calc_calibration_slope(y, pred),
    calibration_intercept = calc_calibration_intercept(y, pred)
  )
}

make_model_matrix <- function(dt, vars, include_intercept = FALSE) {
  d <- copy(as.data.table(dt))
  d[, sex := factor(sex)]
  d[, aki_stage := factor(aki_stage)]
  form <- stats::as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- stats::model.matrix(form, data = d)
  if (!include_intercept && "(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  mm
}

align_matrix_columns <- function(x_train, x_test) {
  all_cols <- union(colnames(x_train), colnames(x_test))
  add_missing <- function(x, cols) {
    miss <- setdiff(cols, colnames(x))
    if (length(miss)) {
      z <- matrix(0, nrow = nrow(x), ncol = length(miss))
      colnames(z) <- miss
      x <- cbind(x, z)
    }
    x[, cols, drop = FALSE]
  }
  list(train = add_missing(x_train, all_cols), test = add_missing(x_test, all_cols))
}

fit_predict_glm <- function(train_dt, test_dt, outcome, vars) {
  train <- copy(as.data.table(train_dt))
  test <- copy(as.data.table(test_dt))
  train[, sex := factor(sex)]
  test[, sex := factor(sex, levels = levels(train$sex))]
  train[, aki_stage := factor(aki_stage)]
  test[, aki_stage := factor(aki_stage, levels = levels(train$aki_stage))]

  form <- stats::as.formula(paste(outcome, "~", paste(vars, collapse = " + ")))
  fit <- tryCatch(stats::glm(form, data = train, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(test)))
  as.numeric(tryCatch(stats::predict(fit, newdata = test, type = "response"),
                      error = function(e) rep(NA_real_, nrow(test))))
}

fit_predict_glmnet <- function(train_dt, test_dt, outcome, vars, alpha = 0.5) {
  y_train <- as.integer(train_dt[[outcome]])
  if (length(unique(y_train)) < 2) return(rep(NA_real_, nrow(test_dt)))

  x_train <- make_model_matrix(train_dt, vars)
  x_test <- make_model_matrix(test_dt, vars)
  aligned <- align_matrix_columns(x_train, x_test)

  pred <- tryCatch({
    cvfit <- glmnet::cv.glmnet(
      x = aligned$train,
      y = y_train,
      family = "binomial",
      alpha = alpha,
      nfolds = GLMNET_NFOLD,
      type.measure = "deviance",
      standardize = TRUE
    )
    as.numeric(stats::predict(cvfit, newx = aligned$test, s = "lambda.min", type = "response"))
  }, error = function(e) rep(NA_real_, nrow(test_dt)))

  pred
}


paired_delong_delta_auc <- function(y, pred_reference, pred_augmented,
                                    conf_level = D2_RUBIN_CONF_LEVEL) {
  y <- as.integer(y)
  ok <- !is.na(y) & is.finite(pred_reference) & is.finite(pred_augmented)
  y <- y[ok]
  p0 <- as.numeric(pred_reference[ok])
  p1 <- as.numeric(pred_augmented[ok])

  if (length(unique(y)) < 2L) {
    return(data.table(
      delta_auc = NA_real_, within_var = NA_real_, within_se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, variance_status = "single_outcome_class"
    ))
  }

  roc0 <- tryCatch(pROC::roc(y, p0, quiet = TRUE, direction = "<"), error = function(e) NULL)
  roc1 <- tryCatch(pROC::roc(y, p1, quiet = TRUE, direction = "<"), error = function(e) NULL)
  if (is.null(roc0) || is.null(roc1)) {
    return(data.table(
      delta_auc = NA_real_, within_var = NA_real_, within_se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, variance_status = "roc_failure"
    ))
  }

  delta <- as.numeric(pROC::auc(roc1) - pROC::auc(roc0))
  tst <- tryCatch(
    pROC::roc.test(roc1, roc0, paired = TRUE, method = "delong",
                   conf.level = conf_level),
    error = function(e) NULL
  )

  if (is.null(tst) || is.null(tst$conf.int) || length(tst$conf.int) < 2L) {
    return(data.table(
      delta_auc = delta, within_var = NA_real_, within_se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, variance_status = "delong_ci_unavailable"
    ))
  }

  ci <- as.numeric(tst$conf.int[c(1, length(tst$conf.int))])
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  se <- abs(ci[2] - ci[1]) / (2 * zcrit)
  data.table(
    delta_auc = delta,
    within_var = se^2,
    within_se = se,
    ci_low = ci[1],
    ci_high = ci[2],
    variance_status = "paired_delong"
  )
}

rubin_pool_scalar <- function(theta, within_var, conf_level = D2_RUBIN_CONF_LEVEL) {
  ok <- is.finite(theta) & is.finite(within_var) & within_var >= 0
  theta <- as.numeric(theta[ok])
  within_var <- as.numeric(within_var[ok])
  m <- length(theta)

  if (m == 0L) {
    return(data.table(
      m_used = 0L, rubin_estimate = NA_real_, within_imp_var_W = NA_real_,
      between_imp_var_B = NA_real_, total_var_T = NA_real_,
      rubin_se = NA_real_, rubin_df = NA_real_,
      rubin_ci_low = NA_real_, rubin_ci_high = NA_real_,
      rubin_ci_covers_zero = NA
    ))
  }

  qbar <- mean(theta)
  W <- mean(within_var)
  B <- if (m > 1L) stats::var(theta) else 0
  Tvar <- W + (1 + 1 / m) * B
  se <- sqrt(max(Tvar, 0))

  if (m > 1L && is.finite(B) && B > 0 && is.finite(W)) {
    r <- ((1 + 1 / m) * B) / max(W, .Machine$double.eps)
    df <- (m - 1) * (1 + 1 / r)^2
  } else {
    df <- Inf
  }

  alpha <- 1 - conf_level
  crit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df = df) else stats::qnorm(1 - alpha / 2)
  lo <- qbar - crit * se
  hi <- qbar + crit * se

  data.table(
    m_used = m,
    rubin_estimate = qbar,
    within_imp_var_W = W,
    between_imp_var_B = B,
    total_var_T = Tvar,
    rubin_se = se,
    rubin_df = df,
    rubin_ci_low = lo,
    rubin_ci_high = hi,
    rubin_ci_covers_zero = lo <= 0 & hi >= 0
  )
}

run_prediction_models <- function(train_dt, test_dt, outcome, label_col = "cluster") {
  base_vars <- c("age", "sex", "sofa", "aki_stage")
  base_label_vars <- c(base_vars, label_col)
  raw_vars <- c(base_vars, feature_names)
  raw_label_vars <- c(base_vars, feature_names, label_col)

  pred_base       <- fit_predict_glm(train_dt, test_dt, outcome, base_vars)
  pred_base_label <- fit_predict_glm(train_dt, test_dt, outcome, base_label_vars)
  pred_raw        <- fit_predict_glmnet(train_dt, test_dt, outcome, raw_vars, alpha = 0.5)
  pred_raw_label  <- fit_predict_glmnet(train_dt, test_dt, outcome, raw_label_vars, alpha = 0.5)

  y <- test_dt[[outcome]]

  out <- safe_rbindlist(list(
    cbind(data.table(model = "Base"), evaluate_predictions(y, pred_base)),
    cbind(data.table(model = "Base+K2"), evaluate_predictions(y, pred_base_label)),
    cbind(data.table(model = "Raw-EN"), evaluate_predictions(y, pred_raw)),
    cbind(data.table(model = "Raw-EN+K2"), evaluate_predictions(y, pred_raw_label))
  ), fill = TRUE)

  decisive <- paired_delong_delta_auc(y, pred_raw, pred_raw_label)
  base_inc <- paired_delong_delta_auc(y, pred_base, pred_base_label)

  out[, `:=`(
    decisive_delta_auc_raw_label_minus_raw = decisive$delta_auc,
    decisive_delta_auc_within_var = decisive$within_var,
    decisive_delta_auc_within_se = decisive$within_se,
    decisive_delta_auc_imp_ci_low = decisive$ci_low,
    decisive_delta_auc_imp_ci_high = decisive$ci_high,
    decisive_variance_status = decisive$variance_status,
    delta_auc_base_label_minus_base = base_inc$delta_auc,
    base_delta_auc_within_var = base_inc$within_var
  )]
  out[]
}


make_noisy_oracle_k2_label <- function(domain2_latent_group,
                                       accuracy = S5_ORACLE_LABEL_ACCURACY,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- as.integer(domain2_latent_group)
  true_label <- ifelse(g == 1L, "C1", "C2")
  flip <- stats::rbinom(length(g), size = 1, prob = 1 - accuracy) == 1
  noisy_label <- true_label
  noisy_label[flip & true_label == "C1"] <- "C2"
  noisy_label[flip & true_label == "C2"] <- "C1"
  factor(noisy_label, levels = c("C1", "C2"))
}

run_domain2 <- function(dt, scenario, repeat_id, imputation_id, split_index) {
  train <- copy(dt[split_index])
  test <- copy(dt[-split_index])

  label_source <- "kmeans_frozen_centroids"
  centroids_dt <- data.table()
  preprocess_dt <- data.table()

  if (scenario == "S5_domain2_incremental_value_positive_control") {
    ## Scenario 5 is a metric-sensitivity positive control for Domain 2. The
    ## fixed noisy oracle label is derived from the hidden latent state, not
    ## from the observed 33-feature pattern. It is an indicator-sensitivity
    ## positive control and is deliberately unavailable to Raw-EN.
    label_source <- "fixed_noisy_oracle_hidden_state_label"
    train[, cluster := factor(domain2_oracle_label, levels = c("C1", "C2"))]
    test[, cluster := factor(domain2_oracle_label, levels = c("C1", "C2"))]

    centroids_dt <- data.table(
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_method = "mice",
      imputation_id = imputation_id,
      analysis_context = "domain2_metric_positive_control_noisy_oracle_label",
      cluster_k = 2,
      label_source = label_source,
      oracle_label_accuracy_target = S5_ORACLE_LABEL_ACCURACY
    )

    preprocess_dt <- data.table(
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_method = "mice",
      imputation_id = imputation_id,
      analysis_context = "domain2_metric_positive_control_noisy_oracle_label",
      feature = "domain2_hidden_state",
      parameter = "oracle_label_accuracy_target",
      value = S5_ORACLE_LABEL_ACCURACY,
      cluster_k = 2,
      label_source = label_source
    )
  } else {
    fit <- fit_kmeans_centroids(
      train, feature_names, k = 2, analysis_context = "domain2_train_split",
      scenario = scenario, repeat_id = repeat_id, imputation_id = imputation_id
    )
    train[, cluster := fit$label]
    applied <- apply_frozen_centroids(
      test, feature_names, fit$centroids, fit$center, fit$scale, fit$winsor_bounds,
      where = paste("domain2_test", scenario, repeat_id, imputation_id)
    )
    test[, cluster := applied$label]

    centroids_dt <- fit$centroids_dt
    preprocess_dt <- fit$preprocess_dt
  }

  outcomes <- c("mortality_30d", "make30")
  out <- safe_rbindlist(lapply(outcomes, function(y) {
    res <- run_prediction_models(train, test, outcome = y, label_col = "cluster")
    res[, outcome := y]
    res
  }), fill = TRUE)

  out[, `:=`(
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_method = "mice",
    imputation_id = imputation_id,
    is_primary_imputation = (imputation_id == PRIMARY_IMPUTATION_ID),
    domain = "Domain2_IncrementalValue",
    domain2_label_source = label_source,
    oracle_label_accuracy_target = ifelse(
      scenario == "S5_domain2_incremental_value_positive_control",
      S5_ORACLE_LABEL_ACCURACY,
      NA_real_
    ),
    oracle_label_fixed_across_imputations =
      scenario == "S5_domain2_incremental_value_positive_control"
  )]

  ## Do not classify Domain 2 from a single split/repeat point estimate.
  ## The decisive Delta AUC is interpreted only after Monte Carlo aggregation
  ## across repeats; otherwise ordinary AUC sampling noise can create false
  ## positive incremental-value flags.
  out[, interpretation_flag := "per_repeat_not_classified_mc_interval_required"]
  out[, interpretation_text := ifelse(
    scenario == "S5_domain2_incremental_value_positive_control",
    "Scenario 5 uses a noisy oracle K2 label to test whether the Domain 2 Delta AUC metric detects true incremental label value beyond raw features.",
    "Domain 2 incremental value is classified using the confidence interval for the across-repeat mean across repeats, not the per-repeat point estimate."
  )]

  list(summary = out, centroids_dt = centroids_dt, preprocess_dt = preprocess_dt)
}

run_domain3 <- function(deriv_dt, valid_dt, scenario, repeat_id, imputation_id) {
  deriv <- copy(deriv_dt)
  valid <- copy(valid_dt)

  fit <- fit_kmeans_centroids(
    deriv, feature_names, k = 2, analysis_context = "domain3_derivation",
    scenario = scenario, repeat_id = repeat_id, imputation_id = imputation_id
  )
  deriv[, cluster := fit$label]

  ## Strict frozen transport: validation uses its MICE-completed values, but uses
  ## derivation winsor bounds, derivation mean/SD, and derivation centroids.
  applied <- apply_frozen_centroids(
    valid, feature_names, fit$centroids, fit$center, fit$scale, fit$winsor_bounds,
    where = paste("domain3_validation", scenario, repeat_id, imputation_id)
  )
  valid[, cluster := applied$label]

  prev_deriv <- mean(deriv$cluster == "C1")
  prev_valid <- mean(valid$cluster == "C1")
  prev_drift <- prev_valid - prev_deriv

  outcomes <- c("mortality_30d", "make30")
  perf <- safe_rbindlist(lapply(outcomes, function(y) {
    res <- run_prediction_models(deriv, valid, outcome = y, label_col = "cluster")
    res[, outcome := y]
    res
  }), fill = TRUE)

  perf[, `:=`(
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_method = "mice",
    imputation_id = imputation_id,
    is_primary_imputation = (imputation_id == PRIMARY_IMPUTATION_ID),
    domain = "Domain3_Transportability",
    derivation_C1_prevalence = prev_deriv,
    validation_C1_prevalence = prev_valid,
    prevalence_drift = prev_drift
  )]

  auc_raw <- perf[model == "Raw-EN", .(raw_auc = AUC), by = outcome]
  auc_label <- perf[model == "Base+K2", .(label_auc = AUC), by = outcome]
  perf <- merge(perf, auc_raw, by = "outcome", all.x = TRUE)
  perf <- merge(perf, auc_label, by = "outcome", all.x = TRUE)
  perf[, raw_minus_label_auc := raw_auc - label_auc]
  perf[, interpretation_flag := fifelse(
    abs(prevalence_drift) > 0.10 | raw_minus_label_auc > 0.02,
    "transport_drift_or_raw_model_advantage",
    "apparently_transportable"
  )]

  list(summary = perf, centroids_dt = fit$centroids_dt, preprocess_dt = fit$preprocess_dt)
}

## =============================================================================
## 8. Domain 4: positivity and actionability diagnostics
## =============================================================================

effective_sample_size <- function(w) {
  w <- as.numeric(w)
  w <- w[is.finite(w)]
  if (!length(w)) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

fit_propensity <- function(dt) {
  vars <- c(
    "age", "sex", "sofa", "aki_stage",
    "creatinine_max", "bun_max", "lactate_max", "bicarbonate_min",
    "ph_min", "urine_output_24h", "mbp_min",
    "vasopressor_dose_max", "shock_index_max"
  )
  d <- copy(as.data.table(dt))
  d[, sex := factor(sex)]
  d[, aki_stage := factor(aki_stage)]
  form <- stats::as.formula(paste("treatment ~", paste(vars, collapse = " + ")))
  fit <- tryCatch(stats::glm(form, data = d, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, nrow(d)))
  ps <- tryCatch(stats::predict(fit, type = "response"), error = function(e) rep(NA_real_, nrow(d)))
  clamp_prob(ps, eps = 1e-4)
}

calc_positivity_summary <- function(dt, ps) {
  d <- copy(as.data.table(dt))
  d[, ps := ps]
  d[, ipw := ifelse(treatment == 1, 1 / ps, 1 / (1 - ps))]
  d[, ipw_trunc20 := pmin(ipw, 20)]

  count_by_cluster <- d[, .(
    n = .N,
    treated_n = sum(treatment == 1),
    untreated_n = sum(treatment == 0),
    treated_prop = mean(treatment == 1),
    min_cell_count = min(sum(treatment == 1), sum(treatment == 0))
  ), by = cluster]

  ps_treated <- d[treatment == 1, ps]
  ps_untreated <- d[treatment == 0, ps]

  overlap_low <- max(min(ps_treated, na.rm = TRUE), min(ps_untreated, na.rm = TRUE))
  overlap_high <- min(max(ps_treated, na.rm = TRUE), max(ps_untreated, na.rm = TRUE))
  prop_outside_overlap <- mean(d$ps < overlap_low | d$ps > overlap_high, na.rm = TRUE)

  weight_summary <- data.table(
    ps_min_treated = min(ps_treated, na.rm = TRUE),
    ps_max_treated = max(ps_treated, na.rm = TRUE),
    ps_min_untreated = min(ps_untreated, na.rm = TRUE),
    ps_max_untreated = max(ps_untreated, na.rm = TRUE),
    overlap_low = overlap_low,
    overlap_high = overlap_high,
    prop_outside_overlap = prop_outside_overlap,
    ipw_mean = mean(d$ipw, na.rm = TRUE),
    ipw_p95 = as.numeric(stats::quantile(d$ipw, 0.95, na.rm = TRUE)),
    ipw_p99 = as.numeric(stats::quantile(d$ipw, 0.99, na.rm = TRUE)),
    ipw_max = max(d$ipw, na.rm = TRUE),
    extreme_w_gt10 = mean(d$ipw > 10, na.rm = TRUE),
    extreme_w_gt20 = mean(d$ipw > 20, na.rm = TRUE),
    ess_ipw = effective_sample_size(d$ipw),
    ess_ipw_trunc20 = effective_sample_size(d$ipw_trunc20)
  )

  list(
    count_by_cluster = count_by_cluster,
    weight_summary = weight_summary,
    ps_data = d[, .(cluster, treatment, ps, ipw)]
  )
}

estimate_interaction_model <- function(dt) {
  d <- copy(as.data.table(dt))
  d[, treatment := as.integer(treatment)]
  d[, mortality_30d := as.integer(mortality_30d)]
  d[, sex := factor(sex)]
  d[, aki_stage := factor(aki_stage)]
  d[, cluster := factor(cluster, levels = c("C1", "C2"))]

  exposure_cells <- data.table::CJ(
    cluster = c("C1", "C2"),
    treatment = 0:1,
    unique = TRUE
  )
  observed_cells <- d[, .(
    cell_n = .N,
    event_n = sum(mortality_30d == 1, na.rm = TRUE),
    nonevent_n = sum(mortality_30d == 0, na.rm = TRUE)
  ), by = .(cluster = as.character(cluster), treatment)]
  exposure_cells <- merge(
    exposure_cells,
    observed_cells,
    by = c("cluster", "treatment"),
    all.x = TRUE
  )
  for (v in c("cell_n", "event_n", "nonevent_n")) {
    exposure_cells[is.na(get(v)), (v) := 0L]
  }

  get_cell <- function(cl, tr, col = "cell_n") {
    exposure_cells[cluster == cl & treatment == tr, get(col)][1]
  }

  n_c1_t0 <- get_cell("C1", 0)
  n_c1_t1 <- get_cell("C1", 1)
  n_c2_t0 <- get_cell("C2", 0)
  n_c2_t1 <- get_cell("C2", 1)
  min_exposure_cell_n <- min(exposure_cells$cell_n)
  empty_exposure_cell <- any(exposure_cells$cell_n == 0)
  empty_outcome_cell <- any(
    exposure_cells$cell_n > 0 &
      (exposure_cells$event_n == 0 | exposure_cells$nonevent_n == 0)
  )

  base_result <- data.table(
    interaction_estimate = NA_real_,
    interaction_se = NA_real_,
    interaction_p = NA_real_,
    interaction_or = NA_real_,
    interaction_ci_low = NA_real_,
    interaction_ci_high = NA_real_,
    interaction_model_status = NA_character_,
    interaction_identifiable = FALSE,
    glm_converged = FALSE,
    quasi_separation_flag = FALSE,
    empty_exposure_cell = empty_exposure_cell,
    empty_outcome_cell = empty_outcome_cell,
    min_exposure_cell_n = min_exposure_cell_n,
    n_C1_T0 = n_c1_t0,
    n_C1_T1 = n_c1_t1,
    n_C2_T0 = n_c2_t0,
    n_C2_T1 = n_c2_t1,
    interaction_warning = NA_character_
  )

  if (empty_exposure_cell) {
    base_result[, interaction_model_status := "nonidentifiable_empty_treatment_by_cluster_cell"]
    return(base_result)
  }

  warning_messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(
        mortality_30d ~ treatment * cluster + age + sex + sofa + aki_stage,
        data = d,
        family = stats::binomial()
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    base_result[, `:=`(
      interaction_model_status = "model_fit_error",
      interaction_warning = conditionMessage(fit)
    )]
    return(base_result)
  }

  sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  term <- if (!is.null(sm)) {
    grep("^treatment:cluster|^cluster.*:treatment", rownames(sm), value = TRUE)
  } else character()

  if (!length(term)) {
    base_result[, `:=`(
      interaction_model_status = "interaction_term_not_estimable",
      glm_converged = isTRUE(fit$converged),
      interaction_warning = paste(unique(warning_messages), collapse = " | ")
    )]
    return(base_result)
  }

  term <- term[1]
  est <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  p <- unname(sm[term, "Pr(>|z|)"])
  all_coef <- stats::coef(fit)
  fitted_p <- tryCatch(stats::fitted(fit), error = function(e) numeric())

  warn_text <- paste(unique(warning_messages), collapse = " | ")
  warning_sep <- grepl(
    "fitted probabilities numerically 0 or 1|did not converge|separation",
    warn_text,
    ignore.case = TRUE
  )
  extreme_coef <- any(is.finite(all_coef) & abs(all_coef) > 10)
  huge_se <- !is.finite(se) || se > 10
  extreme_fitted <- length(fitted_p) > 0 &&
    any(fitted_p < 1e-6 | fitted_p > 1 - 1e-6, na.rm = TRUE)
  quasi_sep <- empty_outcome_cell || warning_sep || extreme_coef || huge_se ||
    extreme_fitted || !isTRUE(fit$converged)

  status <- if (quasi_sep) {
    "quasi_separation_or_unstable_extrapolation"
  } else if (min_exposure_cell_n < 10) {
    "estimable_but_sparse_cell_unstable"
  } else {
    "estimable_no_major_separation_flag"
  }

  identifiable <- status == "estimable_no_major_separation_flag"

  data.table(
    interaction_estimate = est,
    interaction_se = se,
    interaction_p = p,
    interaction_or = if (is.finite(est)) exp(est) else NA_real_,
    interaction_ci_low = if (is.finite(est) && is.finite(se)) est - 1.96 * se else NA_real_,
    interaction_ci_high = if (is.finite(est) && is.finite(se)) est + 1.96 * se else NA_real_,
    interaction_model_status = status,
    interaction_identifiable = identifiable,
    glm_converged = isTRUE(fit$converged),
    quasi_separation_flag = quasi_sep,
    empty_exposure_cell = empty_exposure_cell,
    empty_outcome_cell = empty_outcome_cell,
    min_exposure_cell_n = min_exposure_cell_n,
    n_C1_T0 = n_c1_t0,
    n_C1_T1 = n_c1_t1,
    n_C2_T0 = n_c2_t0,
    n_C2_T1 = n_c2_t1,
    interaction_warning = if (nzchar(warn_text)) warn_text else NA_character_
  )
}


run_domain4 <- function(dt, scenario, repeat_id, imputation_id) {
  fit <- fit_kmeans_centroids(
    dt, feature_names, k = 2, analysis_context = "domain4_derivation",
    scenario = scenario, repeat_id = repeat_id, imputation_id = imputation_id
  )

  d <- copy(dt)
  d[, cluster := fit$label]

  ps <- fit_propensity(d)
  pos <- calc_positivity_summary(d, ps)
  interaction <- estimate_interaction_model(d)

  min_treated <- min(pos$count_by_cluster$treated_n, na.rm = TRUE)
  min_untreated <- min(pos$count_by_cluster$untreated_n, na.rm = TRUE)
  min_cell <- min(pos$count_by_cluster$min_cell_count, na.rm = TRUE)

  flag <- ifelse(
    min_cell < 30 ||
      pos$weight_summary$extreme_w_gt10 > 0.01 ||
      pos$weight_summary$ess_ipw < 0.25 * nrow(d) ||
      pos$weight_summary$prop_outside_overlap > 0.10 ||
      !isTRUE(interaction$interaction_identifiable[1]),
    "positivity_or_interaction_nonidentifiability_do_not_interpret_hte",
    "positivity_acceptable_for_exploratory_hte"
  )

  summary <- cbind(
    data.table(
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_method = "mice",
      imputation_id = imputation_id,
      is_primary_imputation = (imputation_id == PRIMARY_IMPUTATION_ID),
      domain = "Domain4_ActionabilityPositivity",
      min_treated_by_cluster = min_treated,
      min_untreated_by_cluster = min_untreated,
      min_cell_count_by_cluster = min_cell
    ),
    pos$weight_summary,
    interaction
  )

  if (scenario != "S4_actionability_positivity_failure") {
    flag <- "not_primary_actionability_scenario"
  }

  summary[, interpretation_flag := flag]
  summary[, interpretation_text := ifelse(
    flag == "not_primary_actionability_scenario",
    "Domain 4 is interpreted primarily in Scenario 4; this result is retained for completeness only.",
    ifelse(
      flag == "positivity_or_interaction_nonidentifiability_do_not_interpret_hte",
      "Treatment overlap or interaction identifiability is insufficient; the phenotype-by-treatment term is an instability diagnostic and should not be interpreted as causal HTE.",
      "Treatment overlap appears acceptable for exploratory interaction diagnostics; causal interpretation still requires assumptions."
    )
  )]

  counts <- copy(pos$count_by_cluster)
  counts[, `:=`(scenario = scenario, repeat_id = repeat_id, imputation_method = "mice", imputation_id = imputation_id)]

  ps_data <- copy(pos$ps_data)
  ps_data[, `:=`(scenario = scenario, repeat_id = repeat_id, imputation_method = "mice", imputation_id = imputation_id)]

  list(
    summary = summary,
    counts = counts,
    ps_data = ps_data,
    centroids_dt = fit$centroids_dt,
    preprocess_dt = fit$preprocess_dt
  )
}

compute_label_stability_ari <- function(label_list, scenario, repeat_id) {
  if (length(label_list) < 1L || is.null(label_list[[PRIMARY_IMPUTATION_ID]])) {
    return(list(by_imputation = data.table(), summary = data.table()))
  }

  primary <- copy(label_list[[PRIMARY_IMPUTATION_ID]][, .(
    patient_id, primary_cluster = cluster
  )])

  by_imp <- safe_rbindlist(lapply(seq_along(label_list), function(ii) {
    cur <- label_list[[ii]]
    if (is.null(cur) || !nrow(cur)) return(NULL)
    z <- merge(
      primary,
      cur[, .(patient_id, current_cluster = cluster)],
      by = "patient_id",
      all = FALSE
    )
    ari <- if (nrow(z)) {
      tryCatch(
        mclust::adjustedRandIndex(z$primary_cluster, z$current_cluster),
        error = function(e) NA_real_
      )
    } else NA_real_

    data.table(
      scenario = scenario,
      repeat_id = repeat_id,
      reference_imputation_id = PRIMARY_IMPUTATION_ID,
      imputation_id = ii,
      is_primary_imputation = (ii == PRIMARY_IMPUTATION_ID),
      matched_n = nrow(z),
      label_stability_ari = ari
    )
  }), fill = TRUE)

  nonprimary <- by_imp[imputation_id != PRIMARY_IMPUTATION_ID]
  summary <- data.table(
    scenario = scenario,
    repeat_id = repeat_id,
    primary_imputation_id = PRIMARY_IMPUTATION_ID,
    m = length(label_list),
    ari_mean_vs_primary = if (nrow(nonprimary)) mean(nonprimary$label_stability_ari, na.rm = TRUE) else NA_real_,
    ari_median_vs_primary = if (nrow(nonprimary)) stats::median(nonprimary$label_stability_ari, na.rm = TRUE) else NA_real_,
    ari_min_vs_primary = if (nrow(nonprimary)) min(nonprimary$label_stability_ari, na.rm = TRUE) else NA_real_,
    ari_prop_ge_0_90 = if (nrow(nonprimary)) mean(nonprimary$label_stability_ari >= 0.90, na.rm = TRUE) else NA_real_
  )

  list(by_imputation = by_imp, summary = summary)
}

## =============================================================================
## 9. One repeat: generate data, MICE m=5, analyze each imputation
## =============================================================================

run_one_repeat <- function(scenario, repeat_id, sensitivity_id = "main", params = list()) {
  apply_sensitivity_params(params)
  local_rep <- if ("local_rep_id" %in% names(params)) params$local_rep_id else repeat_id
  log_msg("Running ", sensitivity_id, " | ", scenario, " | repeat ", local_rep, "/", N_ACTIVE_REP)

  seed_offset <- if ("sens_setting_index" %in% names(params)) as.integer(params$sens_setting_index) * 10000000L else 0L
  seed_base <- GLOBAL_SEED + seed_offset + repeat_id * 1000 + match(scenario, SCENARIOS) * 100000

  deriv_raw <- generate_dataset(
    scenario = scenario, cohort = "derivation", n = N_TRAIN,
    repeat_id = repeat_id, seed = seed_base + 1
  )

  valid_raw <- generate_dataset(
    scenario = scenario, cohort = "validation", n = N_EXTERNAL,
    repeat_id = repeat_id, seed = seed_base + 2
  )

  deriv_imp_list <- make_mice_imputations(
    deriv_raw, scenario, repeat_id, cohort = "derivation",
    m = MICE_M, maxit = MICE_MAXIT, seed = seed_base + 101
  )

  valid_imp_list <- make_mice_imputations(
    valid_raw, scenario, repeat_id, cohort = "validation",
    m = MICE_M, maxit = MICE_MAXIT, seed = seed_base + 202
  )

  if (scenario == "S5_domain2_incremental_value_positive_control") {
    deriv_oracle_ref <- as.character(deriv_imp_list[[PRIMARY_IMPUTATION_ID]]$domain2_oracle_label)
    valid_oracle_ref <- as.character(valid_imp_list[[PRIMARY_IMPUTATION_ID]]$domain2_oracle_label)
    deriv_fixed <- all(vapply(deriv_imp_list, function(z) {
      identical(as.character(z$domain2_oracle_label), deriv_oracle_ref)
    }, logical(1)))
    valid_fixed <- all(vapply(valid_imp_list, function(z) {
      identical(as.character(z$domain2_oracle_label), valid_oracle_ref)
    }, logical(1)))
    if (!deriv_fixed || !valid_fixed) {
      stop("Scenario 5 oracle labels changed across MICE completed datasets.",
           call. = FALSE)
    }
  }

  set.seed(seed_base + 303)
  split_index <- sample.int(nrow(deriv_raw), size = floor(0.70 * nrow(deriv_raw)))

  out <- list(
    domain1 = list(), gap = list(), domain2 = list(), domain3 = list(),
    domain4 = list(), domain4_counts = list(), domain4_ps = list(),
    centroids = list(), preprocess = list(),
    label_stability = list(), label_stability_summary = list()
  )
  d1_label_list <- vector("list", MICE_M)

  for (ii in seq_len(MICE_M)) {
    deriv <- deriv_imp_list[[ii]]
    valid <- valid_imp_list[[ii]]

    d1 <- run_domain1(deriv, scenario, repeat_id, ii)
    d2 <- run_domain2(deriv, scenario, repeat_id, ii, split_index)
    d3 <- run_domain3(deriv, valid, scenario, repeat_id, ii)
    d4 <- run_domain4(deriv, scenario, repeat_id, ii)

    d1_label_list[[ii]] <- d1$labels_dt
    out$domain1[[ii]] <- d1$summary
    out$gap[[ii]] <- d1$gap_table
    out$domain2[[ii]] <- d2$summary
    out$domain3[[ii]] <- d3$summary
    out$domain4[[ii]] <- d4$summary
    out$domain4_counts[[ii]] <- d4$counts

    if (scenario == "S4_actionability_positivity_failure" &&
        repeat_id <= min(5, N_ACTIVE_REP)) {
      out$domain4_ps[[ii]] <- d4$ps_data
    }

    out$centroids[[length(out$centroids) + 1L]] <- d1$centroids_dt
    out$centroids[[length(out$centroids) + 1L]] <- d2$centroids_dt
    out$centroids[[length(out$centroids) + 1L]] <- d3$centroids_dt
    out$centroids[[length(out$centroids) + 1L]] <- d4$centroids_dt

    out$preprocess[[length(out$preprocess) + 1L]] <- d1$preprocess_dt
    out$preprocess[[length(out$preprocess) + 1L]] <- d2$preprocess_dt
    out$preprocess[[length(out$preprocess) + 1L]] <- d3$preprocess_dt
    out$preprocess[[length(out$preprocess) + 1L]] <- d4$preprocess_dt
  }

  stability <- compute_label_stability_ari(d1_label_list, scenario, repeat_id)
  out$label_stability[[1L]] <- stability$by_imputation
  out$label_stability_summary[[1L]] <- stability$summary

  res <- lapply(out, function(x) safe_rbindlist(x, fill = TRUE))
  attach_sensitivity_metadata(res, params)
}


## =============================================================================
## 10. Batch run
## =============================================================================

all_d1 <- list(); all_gap <- list(); all_d2 <- list(); all_d3 <- list()
all_d4 <- list(); all_d4_counts <- list(); all_d4_ps <- list()
all_centroids <- list(); all_preprocess <- list(); error_log <- list()
all_label_stability <- list(); all_label_stability_summary <- list()

counter <- 1L



tasks <- SENSITIVITY_GRID[, {
  data.table(
    scenario = rep(SCENARIOS, each = N_ACTIVE_REP),
    local_rep_id = rep(seq_len(N_ACTIVE_REP), times = length(SCENARIOS))
  )
}, by = names(SENSITIVITY_GRID)]
## repeat_id is made unique across sensitivity settings to prevent accidental
## pooling before the sensitivity metadata is re-attached to the master table.
tasks[, rep_id := as.integer((sens_setting_index - 1L) * 100000L + local_rep_id)]
data.table::fwrite(tasks, file.path(TAB_DIR, "sensitivity_task_manifest.csv"))
log_msg("Total parallel tasks: ", nrow(tasks))


res_list <- future.apply::future_lapply(
  seq_len(nrow(tasks)),
  function(i) {
    sc <- tasks$scenario[i]
    rp <- tasks$rep_id[i]
    params <- as.list(tasks[i])
    tryCatch(
      run_one_repeat(sc, rp, sensitivity_id = tasks$sensitivity_id[i], params = params),
      error = function(e) {


        message("ERROR in ", sc, " repeat ", rp, ": ", conditionMessage(e))
        list(
          error = data.table(
            scenario = sc,
            repeat_id = rp,
            error_message = conditionMessage(e)
          )
        )
      }
    )
  },
  future.seed = TRUE
)


counter <- 1L
for (res in res_list) {
  if (!is.null(res$error)) {
    error_log[[length(error_log) + 1L]] <- res$error
    next
  }
  if (!is.null(res)) {
    all_d1[[counter]] <- res$domain1
    all_gap[[counter]] <- res$gap
    all_d2[[counter]] <- res$domain2
    all_d3[[counter]] <- res$domain3
    all_d4[[counter]] <- res$domain4
    all_d4_counts[[counter]] <- res$domain4_counts
    all_d4_ps[[counter]] <- res$domain4_ps
    all_centroids[[counter]] <- res$centroids
    all_preprocess[[counter]] <- res$preprocess
    all_label_stability[[counter]] <- res$label_stability
    all_label_stability_summary[[counter]] <- res$label_stability_summary
    counter <- counter + 1L
  }
}

log_msg("Parallel execution completed.")

domain1_by_imp <- safe_rbindlist(all_d1, fill = TRUE)
gap_by_imp <- safe_rbindlist(all_gap, fill = TRUE)
domain2_by_imp <- safe_rbindlist(all_d2, fill = TRUE)
domain3_by_imp <- safe_rbindlist(all_d3, fill = TRUE)
domain4_by_imp <- safe_rbindlist(all_d4, fill = TRUE)
domain4_counts_by_imp <- safe_rbindlist(all_d4_counts, fill = TRUE)
domain4_ps_by_imp <- safe_rbindlist(all_d4_ps, fill = TRUE)
centroids_long <- safe_rbindlist(all_centroids, fill = TRUE)
preprocess_long <- safe_rbindlist(all_preprocess, fill = TRUE)
label_stability_by_imp <- safe_rbindlist(all_label_stability, fill = TRUE)
label_stability_summary <- safe_rbindlist(all_label_stability_summary, fill = TRUE)
error_dt <- safe_rbindlist(error_log, fill = TRUE)

domain1_primary <- domain1_by_imp[imputation_id == PRIMARY_IMPUTATION_ID]
domain2_primary <- domain2_by_imp[imputation_id == PRIMARY_IMPUTATION_ID]
domain3_primary <- domain3_by_imp[imputation_id == PRIMARY_IMPUTATION_ID]
domain4_primary <- domain4_by_imp[imputation_id == PRIMARY_IMPUTATION_ID]
primary_centroids_long <- centroids_long[imputation_id == PRIMARY_IMPUTATION_ID]
primary_preprocess_long <- preprocess_long[imputation_id == PRIMARY_IMPUTATION_ID]

## Check for duplicated names after major binding.
check_bound_table_names <- function(...) {
  objs <- list(...)
  obj_names <- names(objs)
  for (i in seq_along(objs)) {
    nm <- names(objs[[i]])
    dup <- unique(nm[duplicated(nm)])
    if (length(dup)) {
      stop("Duplicated column names in ", obj_names[i], ": ",
           paste(dup, collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

check_bound_table_names(
  domain1_by_imp = domain1_by_imp,
  gap_by_imp = gap_by_imp,
  domain2_by_imp = domain2_by_imp,
  domain3_by_imp = domain3_by_imp,
  domain4_by_imp = domain4_by_imp,
  domain4_counts_by_imp = domain4_counts_by_imp,
  label_stability_by_imp = label_stability_by_imp,
  label_stability_summary = label_stability_summary,
  centroids_long = centroids_long,
  preprocess_long = preprocess_long
)


## =============================================================================
## 11. Across-imputation summaries
## =============================================================================

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else as.numeric(stats::median(x))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

summarise_across_imputations <- function(dt, by_cols) {
  dt <- fix_duplicate_colnames(copy(as.data.table(dt)))
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  num_cols <- setdiff(
    num_cols,
    c(by_cols, "repeat_id", "imputation_id", "is_primary_imputation")
  )

  mean_dt <- dt[, lapply(.SD, safe_mean), by = by_cols, .SDcols = num_cols]
  median_dt <- dt[, lapply(.SD, safe_median), by = by_cols, .SDcols = num_cols]
  sd_dt <- dt[, lapply(.SD, safe_sd), by = by_cols, .SDcols = num_cols]
  primary_dt <- dt[imputation_id == PRIMARY_IMPUTATION_ID,
                   lapply(.SD, function(x) {
                     if (!length(x)) NA_real_ else as.numeric(x[1])
                   }),
                   by = by_cols, .SDcols = num_cols]

  data.table::setnames(median_dt, num_cols, paste0(num_cols, "_median"))
  data.table::setnames(sd_dt, num_cols, paste0(num_cols, "_sd_across_imputations"))
  data.table::setnames(primary_dt, num_cols, paste0(num_cols, "_primary_imp1"))

  out <- Reduce(
    function(x, y) merge(x, y, by = by_cols, all = TRUE),
    list(mean_dt, median_dt, sd_dt, primary_dt)
  )

  flag_cols <- intersect(
    c("interpretation_flag", "interpretation_text", "domain2_label_source",
      "interaction_model_status"),
    names(dt)
  )

  for (fc in flag_cols) {
    mode_tab <- dt[, .N, by = c(by_cols, fc)]
    mode_tab[, mode_prop := N / sum(N), by = by_cols]
    data.table::setorderv(
      mode_tab,
      c(by_cols, "N"),
      c(rep(1, length(by_cols)), -1)
    )
    mode_dt <- mode_tab[, .SD[1], by = by_cols]
    mode_dt[, N := NULL]
    data.table::setnames(
      mode_dt,
      c(fc, "mode_prop"),
      c(fc, paste0(fc, "_mode_proportion"))
    )
    out <- merge(out, mode_dt, by = by_cols, all.x = TRUE)
  }

  p_cols <- grep("_p$", num_cols, value = TRUE)
  if (length(p_cols)) {
    prop_dt <- dt[, lapply(.SD, function(x) {
      ok <- is.finite(x)
      if (!any(ok)) NA_real_ else mean(x[ok] < 0.05)
    }), by = by_cols, .SDcols = p_cols]
    data.table::setnames(prop_dt, p_cols, paste0(p_cols, "_prop_lt_0_05"))
    out <- merge(out, prop_dt, by = by_cols, all.x = TRUE)
  }

  logical_cols <- names(dt)[vapply(dt, is.logical, logical(1))]
  logical_cols <- setdiff(logical_cols, c(by_cols, "is_primary_imputation"))
  if (length(logical_cols)) {
    ## For logical diagnostics, retain the original unsuffixed column as the
    ## across-imputation proportion TRUE so downstream master-table code can
    ## use it directly. A suffixed audit copy is also retained explicitly.
    lprop <- dt[, lapply(.SD, function(x) {
      ok <- !is.na(x)
      if (!any(ok)) NA_real_ else mean(x[ok])
    }), by = by_cols, .SDcols = logical_cols]
    out <- merge(out, lprop, by = by_cols, all.x = TRUE)

    lprop_audit <- data.table::copy(lprop)
    data.table::setnames(
      lprop_audit, logical_cols, paste0(logical_cols, "_proportion_true")
    )
    out <- merge(out, lprop_audit, by = by_cols, all.x = TRUE)
  }

  out[, `:=`(
    imputation_method = "mice",
    imputation_id = 0L,
    primary_imputation_id = PRIMARY_IMPUTATION_ID,
    imputation_summary = paste0(
      "mean_median_proportion_across_", MICE_M,
      "_completed_datasets_with_imp1_primary_label"
    )
  )]
  out[]
}


domain1_pooled <- summarise_across_imputations(
  domain1_by_imp,
  by_cols = c("scenario", "repeat_id", "domain", "cluster_k")
)

domain2_pooled <- summarise_across_imputations(
  domain2_by_imp,
  by_cols = c("scenario", "repeat_id", "domain", "outcome", "model")
)

domain3_pooled <- summarise_across_imputations(
  domain3_by_imp,
  by_cols = c("scenario", "repeat_id", "domain", "outcome", "model")
)

domain4_pooled <- summarise_across_imputations(
  domain4_by_imp,
  by_cols = c("scenario", "repeat_id", "domain")
)

domain4_counts_pooled <- summarise_across_imputations(
  domain4_counts_by_imp,
  by_cols = c("scenario", "repeat_id", "cluster")
)

## Rubin-style pooling of the decisive paired Delta AUC across m=5.
domain2_delta_by_imp <- unique(domain2_by_imp[, .(
  scenario, repeat_id, outcome, imputation_id,
  is_primary_imputation,
  domain2_label_source,
  decisive_delta_auc_raw_label_minus_raw,
  decisive_delta_auc_within_var,
  decisive_delta_auc_within_se,
  decisive_delta_auc_imp_ci_low,
  decisive_delta_auc_imp_ci_high,
  decisive_variance_status
)])

domain2_rubin_by_repeat <- domain2_delta_by_imp[, {
  rp <- rubin_pool_scalar(
    theta = decisive_delta_auc_raw_label_minus_raw,
    within_var = decisive_delta_auc_within_var,
    conf_level = D2_RUBIN_CONF_LEVEL
  )
  rp[, `:=`(
    primary_imp1_delta_auc = decisive_delta_auc_raw_label_minus_raw[
      imputation_id == PRIMARY_IMPUTATION_ID
    ][1],
    across_imp_median_delta_auc = safe_median(
      decisive_delta_auc_raw_label_minus_raw
    ),
    across_imp_positive_proportion = mean(
      decisive_delta_auc_raw_label_minus_raw > 0,
      na.rm = TRUE
    ),
    variance_method = "paired_DeLong_within_imputation_plus_Rubin_pooling",
    domain2_label_source = domain2_label_source[1]
  )]
  rp
}, by = .(scenario, repeat_id, outcome)]

domain2_rubin_by_repeat[, interpretation_flag := fifelse(
  rubin_ci_low > 0,
  "repeat_level_rubin_ci_positive",
  fifelse(
    rubin_ci_high < 0,
    "repeat_level_rubin_ci_negative",
    "repeat_level_rubin_ci_includes_zero"
  )
)]

## Variability across imputations for key metrics.
mice_variability_summary <- safe_rbindlist(list(
  label_stability_by_imp[imputation_id != PRIMARY_IMPUTATION_ID, .(
    across_imp_mean = safe_mean(label_stability_ari),
    across_imp_median = safe_median(label_stability_ari),
    across_imp_sd = safe_sd(label_stability_ari),
    across_imp_proportion = mean(label_stability_ari >= 0.90, na.rm = TRUE),
    metric = "K2_label_ARI_vs_primary_imputation1"
  ), by = .(scenario, repeat_id)],
  domain1_by_imp[, .(
    across_imp_mean = safe_mean(ari_to_true_subtype),
    across_imp_median = safe_median(ari_to_true_subtype),
    across_imp_sd = safe_sd(ari_to_true_subtype),
    across_imp_proportion = mean(ari_to_true_subtype > 0.5, na.rm = TRUE),
    metric = "ari_to_true_subtype"
  ), by = .(scenario, repeat_id, domain)],
  domain1_by_imp[, .(
    across_imp_mean = safe_mean(ari_k3_to_true_subtype),
    across_imp_median = safe_median(ari_k3_to_true_subtype),
    across_imp_sd = safe_sd(ari_k3_to_true_subtype),
    across_imp_proportion = mean(ari_k3_to_true_subtype > 0.5, na.rm = TRUE),
    metric = "ari_k3_to_true_subtype"
  ), by = .(scenario, repeat_id, domain)],
  domain2_by_imp[, .(
    across_imp_mean = safe_mean(AUC),
    across_imp_median = safe_median(AUC),
    across_imp_sd = safe_sd(AUC),
    across_imp_proportion = NA_real_,
    metric = paste0(outcome, "_", model, "_AUC")
  ), by = .(scenario, repeat_id, domain, outcome, model)],
  domain3_by_imp[, .(
    across_imp_mean = safe_mean(AUC),
    across_imp_median = safe_median(AUC),
    across_imp_sd = safe_sd(AUC),
    across_imp_proportion = NA_real_,
    metric = paste0(outcome, "_external_", model, "_AUC")
  ), by = .(scenario, repeat_id, domain, outcome, model)],
  domain4_by_imp[, .(
    across_imp_mean = safe_mean(ipw_p99),
    across_imp_median = safe_median(ipw_p99),
    across_imp_sd = safe_sd(ipw_p99),
    across_imp_proportion = mean(
      interpretation_flag ==
        "positivity_or_interaction_nonidentifiability_do_not_interpret_hte",
      na.rm = TRUE
    ),
    metric = "ipw_p99"
  ), by = .(scenario, repeat_id, domain)]
), fill = TRUE)

## =============================================================================
## 12. Master summary table
## =============================================================================

make_master_rows <- function() {
  d1_rows <- safe_rbindlist(list(
    domain1_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain1_Discreteness",
      metric = "outcome_separation_C1_minus_C2",
      estimate = outcome_separation,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "Across-imputation mean; completed dataset 1 defines the primary label and median/proportion summaries are retained in the domain table."
    )],
    domain1_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain1_Discreteness",
      metric = "dip_discriminant_p",
      estimate = dip_discriminant_p,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "Across-imputation mean of the discriminant-axis dip p-value; rejection proportion is reported separately."
    )],
    domain1_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain1_Discreteness",
      metric = "sigclust_effect_sd",
      estimate = sigclust_effect_sd,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "Gaussian-null cluster-index sensitivity effect; interpreted with magnitude rather than p-value alone."
    )],
    domain1_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain1_Discreteness",
      metric = "ari_to_true_subtype",
      estimate = ari_to_true_subtype,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "K2 ARI to known truth; cross-imputation label stability is reported in a separate ARI table."
    )],
    domain1_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain1_Discreteness",
      metric = "ari_k3_to_true_subtype",
      estimate = ari_k3_to_true_subtype,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "K3 ARI is the true-K positive-control recovery diagnostic for Scenario 2."
    )]
  ), fill = TRUE)

  d2_delta <- domain2_rubin_by_repeat[, .(
    scenario, repeat_id,
    imputation_method = "mice",
    imputation_id = 0L,
    imputation_summary = "Rubin_style_pooling_across_5_completed_datasets",
    outcome,
    domain = "Domain2_IncrementalValue",
    metric = paste0(
      outcome,
      "_decisive_delta_auc_RawENplusK2_minus_RawEN"
    ),
    estimate = rubin_estimate,
    lower_ci = rubin_ci_low,
    upper_ci = rubin_ci_high,
    interpretation_flag,
    interpretation_text = "Paired Delta AUC pooled across five completed datasets using within-imputation paired DeLong variance and Rubin variance decomposition."
  )]

  d2_auc <- domain2_pooled[, .(
    scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
    domain = "Domain2_IncrementalValue",
    metric = paste0(outcome, "_", model, "_AUC"),
    estimate = AUC,
    lower_ci = NA_real_, upper_ci = NA_real_,
    interpretation_flag,
    interpretation_text = "AUC is summarized across all five completed datasets; the decisive Delta AUC uses Rubin-style pooling."
  )]

  ## Exactly one prevalence-drift record per scenario and repeat. The flag is
  ## derived from prevalence drift itself rather than model-specific flags.
  d3_prevalence <- domain3_pooled[, .(
    imputation_method = imputation_method[1],
    imputation_id = imputation_id[1],
    imputation_summary = imputation_summary[1],
    prevalence_drift = prevalence_drift[1]
  ), by = .(scenario, repeat_id)]
  d3_prevalence[, interpretation_flag := fifelse(
    abs(prevalence_drift) > 0.10,
    "material_prevalence_drift",
    "no_material_prevalence_drift"
  )]

  d3_rows <- safe_rbindlist(list(
    d3_prevalence[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain3_Transportability",
      metric = "C1_prevalence_drift_validation_minus_derivation",
      estimate = prevalence_drift,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "Exactly one frozen-label prevalence-drift estimate per simulation repeat."
    )],
    domain3_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain3_Transportability",
      metric = paste0(outcome, "_external_", model, "_AUC"),
      estimate = AUC,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "External AUC under frozen transport, summarized across five completed datasets."
    )],
    domain3_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain3_Transportability",
      metric = paste0(outcome, "_external_", model, "_calibration_slope"),
      estimate = calibration_slope,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag,
      interpretation_text = "Calibration slope under frozen transport."
    )]
  ), fill = TRUE)

  d4_rows <- safe_rbindlist(list(
    domain4_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain4_ActionabilityPositivity",
      metric = "min_cell_count_by_cluster",
      estimate = min_cell_count_by_cluster,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag, interpretation_text
    )],
    domain4_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain4_ActionabilityPositivity",
      metric = "ipw_p99",
      estimate = ipw_p99,
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag, interpretation_text
    )],
    domain4_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain4_ActionabilityPositivity",
      metric = "interaction_logOR",
      estimate = interaction_estimate,
      lower_ci = interaction_ci_low,
      upper_ci = interaction_ci_high,
      interpretation_flag, interpretation_text
    )],
    domain4_pooled[, .(
      scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
      domain = "Domain4_ActionabilityPositivity",
      metric = "interaction_identifiable",
      estimate = as.numeric(interaction_identifiable),
      lower_ci = NA_real_, upper_ci = NA_real_,
      interpretation_flag, interpretation_text
    )]
  ), fill = TRUE)

  safe_rbindlist(list(d1_rows, d2_auc, d2_delta, d3_rows, d4_rows), fill = TRUE)
}


master_summary <- make_master_rows()

sensitivity_cols_for_merge <- grep("^(sensitivity_id|sens_)", names(domain1_pooled), value = TRUE)
if (length(sensitivity_cols_for_merge)) {
  sensitivity_repeat_map <- unique(domain1_pooled[, c("scenario", "repeat_id", sensitivity_cols_for_merge), with = FALSE])
  master_summary <- merge(master_summary, sensitivity_repeat_map, by = c("scenario", "repeat_id"), all.x = TRUE)
}
aggregate_by_cols <- c(grep("^(sensitivity_id|sens_)", names(master_summary), value = TRUE), "scenario", "domain", "metric")

aggregate_summary <- master_summary[, {
  x <- estimate[is.finite(estimate)]
  n <- length(x)
  mu <- if (n) mean(x) else NA_real_
  med <- if (n) stats::median(x) else NA_real_
  sdx <- if (n > 1L) stats::sd(x) else NA_real_
  mcse <- if (n > 1L) sdx / sqrt(n) else NA_real_
  crit <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  mean_lo <- if (n > 1L) mu - crit * mcse else NA_real_
  mean_hi <- if (n > 1L) mu + crit * mcse else NA_real_
  emp <- if (n) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE) else c(NA_real_, NA_real_)

  ft <- table(interpretation_flag, useNA = "no")
  flag_mode <- if (length(ft)) names(sort(ft, decreasing = TRUE))[1] else NA_character_
  flag_mode_prop <- if (length(ft)) max(ft) / sum(ft) else NA_real_

  list(
    mean_estimate = mu,
    median_estimate = med,
    sd_across_repeats = sdx,
    mcse_mean = mcse,
    mean_ci_low = mean_lo,
    mean_ci_high = mean_hi,
    empirical_p025 = as.numeric(emp[1]),
    empirical_p975 = as.numeric(emp[2]),
    ## Backward-compatible aliases now refer to the CI for the mean, not the
    ## percentile spread across repetitions.
    mc_low = mean_lo,
    mc_high = mean_hi,
    nonmissing_n = n,
    flag_mode = flag_mode,
    flag_mode_proportion = flag_mode_prop
  )
}, by = aggregate_by_cols]

## Domain 2 flags are assigned from the 95% confidence interval for the
## across-repeat mean of the repeat-level Rubin-pooled Delta AUC. Empirical
## percentiles describe dispersion only and are not used for the flag.
aggregate_summary[, mean_ci_covers_zero := NA]
aggregate_summary[
  domain == "Domain2_IncrementalValue" & grepl("decisive_delta_auc", metric),
  mean_ci_covers_zero := mean_ci_low <= 0 & mean_ci_high >= 0
]
aggregate_summary[
  domain == "Domain2_IncrementalValue" & grepl("decisive_delta_auc", metric),
  flag_mode := fifelse(
    mean_ci_low > 0,
    "mean_ci_positive_incremental_value",
    fifelse(
      mean_ci_high < 0,
      "mean_ci_negative_delta",
      "mean_ci_includes_zero_no_reliable_incremental_value"
    )
  )
]

aggregate_summary[
  domain == "Domain2_IncrementalValue" &
    grepl("decisive_delta_auc", metric) &
    scenario == "S5_domain2_incremental_value_positive_control",
  expected_positive_control := TRUE
]
aggregate_summary[is.na(expected_positive_control), expected_positive_control := FALSE]

rubin_sensitivity_cols <- grep("^(sensitivity_id|sens_)", names(domain2_rubin_by_repeat), value = TRUE)
repeat_detection_by <- c(rubin_sensitivity_cols, "scenario", "outcome")
repeat_detection <- domain2_rubin_by_repeat[, .(
  repeat_rubin_ci_positive_proportion = mean(rubin_ci_low > 0, na.rm = TRUE),
  repeat_rubin_ci_includes_zero_proportion = mean(
    rubin_ci_low <= 0 & rubin_ci_high >= 0,
    na.rm = TRUE
  )
), by = repeat_detection_by]
repeat_detection[, metric := paste0(
  outcome, "_decisive_delta_auc_RawENplusK2_minus_RawEN"
)]
repeat_merge_by <- c(rubin_sensitivity_cols, "scenario", "metric")
aggregate_summary <- merge(
  aggregate_summary,
  repeat_detection[, c(repeat_merge_by,
    "repeat_rubin_ci_positive_proportion",
    "repeat_rubin_ci_includes_zero_proportion"), with = FALSE],
  by = repeat_merge_by,
  all.x = TRUE
)

domain2_delta_mean_ci_flags <- aggregate_summary[
  domain == "Domain2_IncrementalValue" & grepl("decisive_delta_auc", metric)
]

## =============================================================================
## 13. Save tables
## =============================================================================

data.table::fwrite(domain1_by_imp, file.path(TAB_DIR, "domain1_discreteness_by_imputation.csv"))
data.table::fwrite(domain1_primary, file.path(TAB_DIR, "domain1_primary_imputation1.csv"))
data.table::fwrite(domain1_pooled, file.path(TAB_DIR, "domain1_discreteness_summary.csv"))
data.table::fwrite(label_stability_by_imp, file.path(TAB_DIR, "mice_label_stability_ari_by_imputation.csv"))
data.table::fwrite(label_stability_summary, file.path(TAB_DIR, "mice_label_stability_ari_summary.csv"))
data.table::fwrite(gap_by_imp, file.path(TAB_DIR, "domain1_gap_tables_by_imputation.csv"))

data.table::fwrite(domain2_by_imp, file.path(TAB_DIR, "domain2_incremental_value_by_imputation.csv"))
data.table::fwrite(domain2_primary, file.path(TAB_DIR, "domain2_primary_imputation1.csv"))
data.table::fwrite(domain2_pooled, file.path(TAB_DIR, "domain2_incremental_value_summary.csv"))
data.table::fwrite(domain2_delta_by_imp, file.path(TAB_DIR, "domain2_delta_auc_by_imputation_within_variance.csv"))
data.table::fwrite(domain2_rubin_by_repeat, file.path(TAB_DIR, "domain2_delta_auc_rubin_by_repeat.csv"))
data.table::fwrite(domain2_delta_mean_ci_flags, file.path(TAB_DIR, "domain2_delta_auc_mean_ci_flags.csv"))

data.table::fwrite(domain3_by_imp, file.path(TAB_DIR, "domain3_transportability_by_imputation.csv"))
data.table::fwrite(domain3_primary, file.path(TAB_DIR, "domain3_primary_imputation1.csv"))
data.table::fwrite(domain3_pooled, file.path(TAB_DIR, "domain3_transportability_summary.csv"))
data.table::fwrite(
  unique(master_summary[
    domain == "Domain3_Transportability" &
      metric == "C1_prevalence_drift_validation_minus_derivation"
  ]),
  file.path(TAB_DIR, "domain3_prevalence_drift_unique_by_repeat.csv")
)

data.table::fwrite(domain4_by_imp, file.path(TAB_DIR, "domain4_actionability_positivity_by_imputation.csv"))
data.table::fwrite(domain4_primary, file.path(TAB_DIR, "domain4_primary_imputation1.csv"))
data.table::fwrite(domain4_pooled, file.path(TAB_DIR, "domain4_actionability_positivity_summary.csv"))
data.table::fwrite(
  domain4_pooled[, .(
    scenario, repeat_id, imputation_method, imputation_id, imputation_summary,
    interaction_model_status,
    interaction_model_status_mode_proportion,
    interaction_identifiable,
    interaction_identifiable_proportion_true,
    quasi_separation_flag,
    quasi_separation_flag_proportion_true,
    empty_exposure_cell,
    empty_exposure_cell_proportion_true,
    empty_outcome_cell,
    empty_outcome_cell_proportion_true,
    min_exposure_cell_n,
    min_exposure_cell_n_median,
    n_C1_T0, n_C1_T1, n_C2_T0, n_C2_T1,
    n_C1_T0_median, n_C1_T1_median, n_C2_T0_median, n_C2_T1_median
  )],
  file.path(TAB_DIR, "domain4_interaction_identifiability_summary.csv")
)
data.table::fwrite(domain4_counts_by_imp, file.path(TAB_DIR, "domain4_treatment_counts_by_cluster_by_imputation.csv"))
data.table::fwrite(domain4_counts_pooled, file.path(TAB_DIR, "domain4_treatment_counts_by_cluster.csv"))
if (nrow(domain4_ps_by_imp)) data.table::fwrite(domain4_ps_by_imp, file.path(TAB_DIR, "domain4_ps_plot_data_first_repeats.csv"))

data.table::fwrite(centroids_long, file.path(TAB_DIR, "frozen_centroids_long.csv"))
data.table::fwrite(preprocess_long, file.path(TAB_DIR, "frozen_preprocess_parameters_long.csv"))
data.table::fwrite(primary_centroids_long, file.path(TAB_DIR, "primary_imputation1_centroids.csv"))
data.table::fwrite(primary_preprocess_long, file.path(TAB_DIR, "primary_imputation1_preprocess_parameters.csv"))

data.table::fwrite(mice_variability_summary, file.path(TAB_DIR, "mice_imputation_variability_summary.csv"))
data.table::fwrite(master_summary, file.path(TAB_DIR, "simulation_four_domain_master_summary.csv"))
data.table::fwrite(aggregate_summary, file.path(TAB_DIR, "simulation_four_domain_aggregate_summary.csv"))

if (nrow(error_dt)) data.table::fwrite(error_dt, file.path(LOG_DIR, "simulation_error_log.csv"))

log_msg("Tables saved.")

if (isTRUE(CREATE_STANDARD_FIGURES)) {
  tryCatch({

## =============================================================================
## 14. Figures
## =============================================================================

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "bottom"
    )
}

save_pub_plot <- function(p, filename, width = 9, height = 6) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(filename, ".png")), plot = p, width = width, height = height, dpi = 600)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(filename, ".pdf")), plot = p, width = width, height = height)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(filename, ".tiff")), plot = p, width = width, height = height, dpi = 600, compression = "lzw")
}

scenario_labels <- c(
  S1_pure_severity_continuum = "S1: continuum",
  S2_true_discrete_subtypes = "S2: true subtypes",
  S3_transport_drift = "S3: transport drift",
  S4_actionability_positivity_failure = "S4: positivity failure",
  S5_domain2_incremental_value_positive_control = "S5: D2 positive"
)

## Figure S1: overview.
fig_s1_dt <- aggregate_summary[
  metric %in% c(
    "dip_discriminant_p",
    "sigclust_effect_sd",
    "mortality_30d_decisive_delta_auc_RawENplusK2_minus_RawEN",
    "C1_prevalence_drift_validation_minus_derivation",
    "min_cell_count_by_cluster",
    "ipw_p99"
  )
]
fig_s1_dt[, scenario_lab := scenario_labels[scenario]]

p_s1 <- ggplot(fig_s1_dt, aes(x = scenario_lab, y = mean_estimate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mc_low, ymax = mc_high), width = 0.15) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(title = "Figure S1. Four-domain MICE simulation overview",
       x = NULL, y = "Monte Carlo mean with 95% mean CI") +
  theme_pub()
save_pub_plot(p_s1, "Figure_S1_four_domain_overview_mice", width = 10, height = 7)

## Figure S2: Domain 1, S1 vs S2.
fig_s2 <- domain1_pooled[scenario %in% c("S1_pure_severity_continuum", "S2_true_discrete_subtypes")]
fig_s2_long <- melt(
  fig_s2,
  id.vars = c("scenario", "repeat_id"),
  measure.vars = c("dip_discriminant_p", "sigclust_effect_sd", "ari_to_true_subtype", "ari_k3_to_true_subtype"),
  variable.name = "metric", value.name = "estimate"
)
fig_s2_long[, scenario_lab := scenario_labels[scenario]]

p_s2 <- ggplot(fig_s2_long, aes(x = scenario_lab, y = estimate)) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  labs(title = "Figure S2. Discreteness diagnostics after MICE",
       x = NULL, y = "Estimate") +
  theme_pub()
save_pub_plot(p_s2, "Figure_S2_discreteness_diagnostics_mice", width = 12, height = 5)

## Figure S3: Domain 2.
fig_s3_auc <- domain2_pooled[outcome == "mortality_30d"]
fig_s3_auc[, scenario_lab := scenario_labels[scenario]]
fig_s3_auc[, model := factor(model, levels = c("Base", "Base+K2", "Raw-EN", "Raw-EN+K2"))]

p_s3a <- ggplot(fig_s3_auc, aes(x = model, y = AUC)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = function(x) {
    data.frame(y = mean(x, na.rm = TRUE),
               ymin = stats::quantile(x, 0.025, na.rm = TRUE),
               ymax = stats::quantile(x, 0.975, na.rm = TRUE))
  }, geom = "errorbar", width = 0.15) +
  facet_wrap(~ scenario_lab, ncol = 2) +
  labs(title = "Figure S3A. Internal AUC ladder after MICE",
       x = NULL, y = "AUC") +
  theme_pub()

fig_s3_delta <- unique(domain2_pooled[outcome == "mortality_30d", .(
  scenario, repeat_id, delta = decisive_delta_auc_raw_label_minus_raw
)])
fig_s3_delta[, scenario_lab := scenario_labels[scenario]]

p_s3b <- ggplot(fig_s3_delta, aes(x = scenario_lab, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(title = "Figure S3B. Decisive Delta AUC: Raw-EN+K2 minus Raw-EN",
       x = NULL, y = "Delta AUC") +
  theme_pub()

save_pub_plot(p_s3a / p_s3b + patchwork::plot_layout(heights = c(2, 1)),
              "Figure_S3_incremental_value_mice", width = 10, height = 8)

## Figure S4: Domain 3.
fig_s4_prev <- unique(domain3_pooled[, .(scenario, repeat_id, prevalence_drift)])
fig_s4_prev[, scenario_lab := scenario_labels[scenario]]

p_s4a <- ggplot(fig_s4_prev, aes(x = scenario_lab, y = prevalence_drift)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(title = "Figure S4A. Frozen-label prevalence drift after MICE",
       x = NULL, y = "Validation C1 prevalence − derivation C1 prevalence") +
  theme_pub()

fig_s4_perf <- domain3_pooled[outcome == "mortality_30d" & model %in% c("Base+K2", "Raw-EN")]
fig_s4_perf[, scenario_lab := scenario_labels[scenario]]
fig_s4_perf[, model := factor(model, levels = c("Base+K2", "Raw-EN"))]

p_s4b <- ggplot(fig_s4_perf, aes(x = model, y = AUC)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = function(x) {
    data.frame(y = mean(x, na.rm = TRUE),
               ymin = stats::quantile(x, 0.025, na.rm = TRUE),
               ymax = stats::quantile(x, 0.975, na.rm = TRUE))
  }, geom = "errorbar", width = 0.15) +
  facet_wrap(~ scenario_lab, ncol = 2) +
  labs(title = "Figure S4B. External AUC under frozen transport",
       x = NULL, y = "External AUC") +
  theme_pub()

p_s4c <- ggplot(fig_s4_perf, aes(x = model, y = calibration_slope)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = function(x) {
    data.frame(y = mean(x, na.rm = TRUE),
               ymin = stats::quantile(x, 0.025, na.rm = TRUE),
               ymax = stats::quantile(x, 0.975, na.rm = TRUE))
  }, geom = "errorbar", width = 0.15) +
  facet_wrap(~ scenario_lab, ncol = 2) +
  labs(title = "Figure S4C. External calibration slope",
       x = NULL, y = "Calibration slope") +
  theme_pub()

save_pub_plot(p_s4a / p_s4b / p_s4c, "Figure_S4_transport_drift_mice", width = 10, height = 10)

## Figure S5: Domain 4.
fig_s5_counts <- domain4_counts_pooled[scenario == "S4_actionability_positivity_failure"]
fig_s5_counts[, cluster := factor(cluster, levels = c("C1", "C2"))]

p_s5a <- ggplot(fig_s5_counts, aes(x = cluster, y = treated_prop)) +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(title = "Figure S5A. Treatment proportion by K2 phenotype",
       x = "K2 phenotype", y = "Treated proportion") +
  theme_pub()

p_s5b <- ggplot(fig_s5_counts, aes(x = cluster, y = treated_n)) +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(title = "Figure S5B. Treated count by K2 phenotype",
       x = "K2 phenotype", y = "Treated n") +
  theme_pub()

if (nrow(domain4_ps_by_imp)) {
  domain4_ps_by_imp[, treatment_lab := ifelse(treatment == 1, "Treated", "Untreated")]
  p_s5c <- ggplot(domain4_ps_by_imp[scenario == "S4_actionability_positivity_failure"],
                  aes(x = ps, linetype = treatment_lab)) +
    geom_density(linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ cluster, ncol = 2) +
    labs(title = "Figure S5C. Propensity-score overlap by phenotype",
         x = "Estimated propensity score", y = "Density", linetype = NULL) +
    theme_pub()
} else {
  p_s5c <- ggplot() + labs(title = "Figure S5C. No PS data available") + theme_pub()
}

p_s5d <- ggplot(domain4_pooled[scenario == "S4_actionability_positivity_failure"],
                aes(x = "", y = ipw_p99)) +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(title = "Figure S5D. IPW 99th percentile", x = NULL, y = "IPW p99") +
  theme_pub()

save_pub_plot((p_s5a | p_s5b) / (p_s5c | p_s5d),
              "Figure_S5_positivity_failure_mice", width = 11, height = 8)

log_msg("Figures saved.")

## Generic sensitivity plots that do not require all five primary scenarios.
tryCatch({
  sens_d2 <- domain2_delta_mean_ci_flags
  if (exists("sens_d2") && nrow(sens_d2)) {
    sens_d2[, setting_lab := if ("sensitivity_id" %in% names(sens_d2)) sensitivity_id else scenario]
    p_sens_d2 <- ggplot(
      sens_d2[grepl("mortality_30d|make30", metric)],
      aes(x = setting_lab, y = mean_estimate)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = mean_ci_low, ymax = mean_ci_high), width = 0.15) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(title = paste0("Sensitivity Domain 2 mean Delta AUC: ", SENSITIVITY_CATEGORY),
           x = NULL, y = "Mean Delta AUC with mean CI") +
      theme_pub()
    save_pub_plot(p_sens_d2, "Sensitivity_Domain2_delta_auc_mean_ci", width = 10, height = 6)
  }

  sens_d1 <- domain1_pooled
  if (exists("sens_d1") && nrow(sens_d1)) {
    sens_d1[, setting_lab := if ("sensitivity_id" %in% names(sens_d1)) sensitivity_id else scenario]
    d1_vars <- intersect(c("ari_k3_to_true_subtype", "sigclust_effect_sd", "dip_discriminant_p"), names(sens_d1))
    if (length(d1_vars)) {
      sens_d1_long <- data.table::melt(
        sens_d1,
        id.vars = intersect(c("setting_lab", "scenario", "repeat_id"), names(sens_d1)),
        measure.vars = d1_vars,
        variable.name = "metric",
        value.name = "estimate"
      )
      if (nrow(sens_d1_long[is.finite(estimate)])) {
        p_sens_d1 <- ggplot(sens_d1_long[is.finite(estimate)], aes(x = setting_lab, y = estimate)) +
          geom_boxplot(outlier.alpha = 0.25) +
          facet_wrap(~ metric, scales = "free_y") +
          labs(title = paste0("Sensitivity Domain 1 diagnostics: ", SENSITIVITY_CATEGORY),
               x = NULL, y = "Estimate") +
          theme_pub()
        save_pub_plot(p_sens_d1, "Sensitivity_Domain1_diagnostics", width = 10, height = 6)
      }
    }
  }

  log_msg("Sensitivity generic figures saved.")
}, error = function(e) {
  log_msg("Sensitivity generic figure generation skipped/failed: ", conditionMessage(e))
})



  }, error = function(e) {
    log_msg("Standard figure generation skipped/failed: ", conditionMessage(e))
    data.table::fwrite(
      data.table::data.table(
        time = as.character(Sys.time()),
        stage = "standard_figure_generation",
        message = conditionMessage(e)
      ),
      file.path(LOG_DIR, "figure_generation_warning.csv")
    )
  })
} else {
  log_msg("Standard figure generation disabled by CREATE_STANDARD_FIGURES = FALSE.")
}

## =============================================================================
## 15. Markdown report and R objects
## =============================================================================

get_mean_metric <- function(metric_name, scenario_name = NULL) {
  d <- aggregate_summary[metric == metric_name]
  if (!is.null(scenario_name)) d <- d[scenario == scenario_name]
  if (!nrow(d)) return("NA")
  sprintf("%.3f (%.3f to %.3f)", d$mean_estimate[1], d$mc_low[1], d$mc_high[1])
}

report_path <- file.path(OUT_DIR, "simulation_four_domain_report.md")

report_lines <- c(
  "# Semi-synthetic simulation with MICE for a four-domain EHR subphenotype validation framework",
  "",
  paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Main methodological updates",
  "",
  "- Primary imputation used MICE with m = 5 and maxit = 5.",
  "- Each imputed dataset was analyzed separately.",
  "- Completed dataset 1 defines the primary K2 label; all five completed datasets are re-analysed and label stability is quantified by ARI.",
  "- Domain 2 paired Delta AUC was Rubin-style pooled across the five completed datasets; other non-standard metrics were summarized by mean, median, and relevant proportions.",
  "- Formal 1st/99th percentile winsorization was performed after MICE and before standardization.",
  "- External validation used validation-cohort MICE but derivation-frozen winsor bounds, scaling parameters, and centroids. Validation cohorts were not re-clustered.",
  "- Scenario 5 was added strictly as a Domain 2 metric-sensitivity positive control, not as evidence of clustering validity: a hidden latent state was weakly expressed by the 33 raw features, directly affected outcomes, and was represented by a fixed noisy oracle K2 label unavailable to the raw-feature model.",
  "- Domain 2 interpretation flags were assigned from the confidence interval for the across-repeat mean of repeat-level Rubin-pooled Delta AUC; empirical percentiles were retained only as dispersion.",
  "",
  "## Selected aggregate results",
  "",
  paste0("- Scenario 1 discriminant-axis dip p-value: ", get_mean_metric("dip_discriminant_p", "S1_pure_severity_continuum")),
  paste0("- Scenario 2 K2 ARI to true subtype: ", get_mean_metric("ari_to_true_subtype", "S2_true_discrete_subtypes")),
  paste0("- Scenario 2 K3 true-K ARI to true subtype: ", get_mean_metric("ari_k3_to_true_subtype", "S2_true_discrete_subtypes")),
  paste0("- Scenario 3 C1 prevalence drift: ", get_mean_metric("C1_prevalence_drift_validation_minus_derivation", "S3_transport_drift")),
  paste0("- Scenario 4 minimum cell count by cluster: ", get_mean_metric("min_cell_count_by_cluster", "S4_actionability_positivity_failure")),
  paste0("- Scenario 4 IPW p99: ", get_mean_metric("ipw_p99", "S4_actionability_positivity_failure")),
  paste0("- Scenario 5 mortality decisive Delta AUC: ", get_mean_metric("mortality_30d_decisive_delta_auc_RawENplusK2_minus_RawEN", "S5_domain2_incremental_value_positive_control")),
  "",
  "## Manuscript-ready Methods paragraph",
  "",
  "We conducted a semi-synthetic simulation study to evaluate the proposed four-domain validation framework under known data-generating mechanisms. Each derivation dataset included 5,000 simulated patients and 33 EHR-like first-24-hour physiological features, and each validation dataset included 3,000 patients. Five scenarios were considered: a pure severity continuum, true discrete subtypes, external transport drift, severity-driven treatment assignment causing practical non-positivity, and a Domain 2 metric-level positive-control setting in which a hidden latent state was weakly expressed by the 33 raw features, directly affected outcomes, and was represented by a noisy oracle K2 label unavailable to the raw-feature model. In keeping with the empirical analyses, simulated EHR-like missingness was handled using multiple imputation by chained equations (MICE, m = 5, maxit = 5). The imputation model was restricted to pre-outcome physiological features and baseline covariates and did not include outcomes, treatment indicators, derived cluster labels, true latent subtype labels, latent severity, the hidden Domain 2 latent state, or the noisy oracle label except in the explicit Domain 2 positive-control comparison. Completed dataset 1 defined the primary K2 label; all five completed datasets were analyzed separately, cross-imputation label stability was quantified by ARI, paired Domain 2 Delta AUC was Rubin-style pooled, and other non-standard metrics were summarized using across-imputation means, medians, and relevant proportions. Formal 1st/99th percentile winsorization was applied after imputation and before derivation-set standardization. External validation used frozen derivation-set winsor bounds, means, standard deviations, and K-means centroids; external cohorts were not re-clustered. Domain 2 incremental value was classified using confidence interval for the across-repeat means across repeats rather than per-repeat point estimates.",
  "",
  "## Output files",
  "",
  paste0("- Tables: `", TAB_DIR, "`"),
  paste0("- Figures: `", FIG_DIR, "`"),
  paste0("- Logs: `", LOG_DIR, "`")
)

writeLines(report_lines, con = report_path)

saveRDS(
  list(
    domain1_by_imputation = domain1_by_imp,
    domain1_pooled = domain1_pooled,
    domain1_primary = domain1_primary,
    label_stability_by_imputation = label_stability_by_imp,
    label_stability_summary = label_stability_summary,
    domain2_by_imputation = domain2_by_imp,
    domain2_pooled = domain2_pooled,
    domain2_rubin_by_repeat = domain2_rubin_by_repeat,
    domain3_by_imputation = domain3_by_imp,
    domain3_pooled = domain3_pooled,
    domain4_by_imputation = domain4_by_imp,
    domain4_pooled = domain4_pooled,
    centroids_long = centroids_long,
    preprocess_long = preprocess_long,
    mice_variability_summary = mice_variability_summary,
    master_summary = master_summary,
    aggregate_summary = aggregate_summary,
    domain2_delta_mean_ci_flags = domain2_delta_mean_ci_flags,
    settings = list(
      PROJECT_DIR = PROJECT_DIR,
      SIM_ROOT = SIM_ROOT,
      N_TRAIN = N_TRAIN,
      N_EXTERNAL = N_EXTERNAL,
      N_ACTIVE_REP = N_ACTIVE_REP,
      IMPUTATION_METHOD = IMPUTATION_METHOD,
      MICE_M = MICE_M,
      MICE_MAXIT = MICE_MAXIT,
      PRIMARY_IMPUTATION_ID = PRIMARY_IMPUTATION_ID,
      D2_RUBIN_CONF_LEVEL = D2_RUBIN_CONF_LEVEL,
      WINSOR_PROBS = WINSOR_PROBS,
      S5_LATENT_PREVALENCE = S5_LATENT_PREVALENCE,
      S5_FEATURE_DELTA = S5_FEATURE_DELTA,
      S5_LATENT_OUTCOME_LOGOR_MORT = S5_LATENT_OUTCOME_LOGOR_MORT,
      S5_LATENT_OUTCOME_LOGOR_MAKE = S5_LATENT_OUTCOME_LOGOR_MAKE,
      S5_ORACLE_LABEL_ACCURACY = S5_ORACLE_LABEL_ACCURACY,
      S4_STRONG_GAMMA0 = S4_STRONG_GAMMA0,
      S4_STRONG_GAMMA_SEVERITY = S4_STRONG_GAMMA_SEVERITY,
      S4_STRONG_GAMMA_RENAL = S4_STRONG_GAMMA_RENAL,
      S4_STRONG_GAMMA_SHOCK = S4_STRONG_GAMMA_SHOCK,
      GLOBAL_SEED = GLOBAL_SEED
    )
  ),
  file = file.path(RAW_DIR, "simulation_four_domain_all_results_mice.rds")
)

sink(file.path(LOG_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

warnings_path <- file.path(LOG_DIR, "warnings_after_run.txt")
capture.output(warnings(), file = warnings_path)

log_msg("Markdown report saved: ", report_path)
log_msg("Session info saved.")
log_msg("Warnings saved: ", warnings_path)
log_msg("Simulation completed successfully.")
log_msg("Output folder: ", OUT_DIR)

cat("\n============================================================\n")
cat("MICE-based four-domain simulation completed.\n")
cat("Output folder:\n", OUT_DIR, "\n")
cat("Key files:\n")
cat("  ", file.path(TAB_DIR, "simulation_four_domain_master_summary.csv"), "\n")
cat("  ", file.path(TAB_DIR, "mice_imputation_variability_summary.csv"), "\n")
cat("  ", file.path(TAB_DIR, "mice_label_stability_ari_summary.csv"), "\n")
cat("  ", file.path(TAB_DIR, "domain2_delta_auc_rubin_by_repeat.csv"), "\n")
cat("  ", file.path(TAB_DIR, "domain2_delta_auc_mean_ci_flags.csv"), "\n")
cat("  ", file.path(TAB_DIR, "frozen_centroids_long.csv"), "\n")
cat("  ", file.path(TAB_DIR, "frozen_preprocess_parameters_long.csv"), "\n")
cat("  ", report_path, "\n")
cat("============================================================\n")
