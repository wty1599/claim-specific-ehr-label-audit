################################################################################
# 07_gmm_ward_paired_labels_transport_FINAL.R
#
# Combined, standalone redesign of sensitivity analyses 07 and 07b.
#
# Purpose
# -------
# Evaluate whether the behavior and external assignment transportability of
# EHR-derived K=2/K=3 labels are robust to two alternative clustering methods:
#
#   1. Gaussian mixture models (GMM), fitted strictly with mclust::Mclust(...)
#   2. Ward.D2 hierarchical clustering, fitted on a bounded subsample and
#      extended to all patients by nearest-centroid assignment
#
# Required design features
# ------------------------
# - Scenarios S1/S2/S3 are rerun from the same simulated datasets.
# - GMM and Ward labels are paired within the same repeat, imputation, cohort,
#   patient, and K.
# - K=2 and K=3 are both fitted.
# - Individual patient labels are saved.
# - Every multisession worker explicitly loads mclust.
# - Every GMM fit uses mclust::Mclust(...).
# - A failed/invalid GMM fit terminates that repeat.
# - There is no k-means fallback anywhere in this script.
# - MICE m=5, maxit=5; completed dataset 1 is the primary imputation.
# - Derivation preprocessing is frozen for external frozen assignment.
# - External reclustering is independently fitted after external-local
#   preprocessing.
# - Task-level checkpoints permit restart of failed/incomplete repeats.
#
# The script intentionally focuses on:
#   Domain 1: algorithm behavior, latent-truth recovery, label stability
#   Domain 3: external frozen-vs-reclustered assignment transportability
#
# Domain 2 and Domain 4 are not re-estimated here because they are addressed by
# the main simulation and dedicated sensitivities. This avoids making a
# clustering-algorithm sensitivity analysis carry unrelated estimands.
################################################################################

## Clear the workspace so this script is robust to being sourced after other
## scripts in the same R session (matches 03b/06). Without this, leftover
## objects from a prior run can trip future's global-export scan
## ("ls(envir = env): invalid 'envir' argument"). A fresh R session is still
## the cleanest way to run.
rm(list = ls())

options(stringsAsFactors = FALSE)
options(warn = 1)

## =============================================================================
## 0. User settings
## =============================================================================

PROJECT_DIR <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")

ANALYSIS_ID <- "sensitivity_07_gmm_ward_paired_labels_transport_FINAL"
SIM_ROOT <- file.path(PROJECT_DIR, ANALYSIS_ID)

USE_TEST_MODE <- FALSE

N_REP_FORMAL <- 100L          # aligned with the pre-registered plan (B = 100 primary)
N_REP_TEST <- 1L
N_ACTIVE_REP <- if (USE_TEST_MODE) N_REP_TEST else N_REP_FORMAL

SCENARIOS <- c(
  "S1_pure_severity_continuum",
  "S2_true_discrete_subtypes",
  "S3_transport_drift"
)

N_DERIVATION <- if (USE_TEST_MODE) 800L else 5000L
N_EXTERNAL <- if (USE_TEST_MODE) 500L else 3000L

MICE_M <- if (USE_TEST_MODE) 2L else 5L
MICE_MAXIT <- if (USE_TEST_MODE) 2L else 5L
PRIMARY_IMPUTATION_ID <- 1L

K_VALUES <- c(2L, 3L)

## Parallelism tuned for a 16-core / 32 GB workstation.
## GMM (mclust) and Ward are CPU-bound and memory-light here; cores, not RAM,
## are the binding constraint. Leave 2 logical cores for the master R session +
## OS/RStudio and cap at 12 to avoid oversubscription. On Windows,
## future::multisession is the correct plan (multicore/fork is unavailable).
N_WORKERS <- if (USE_TEST_MODE) {
  2L
} else {
  max(1L, min(12L, future::availableCores() - 2L))
}
GLOBAL_SEED <- 2026071207L

WINSOR_PROBS <- c(0.01, 0.99)

## GMM configuration.
## GMM is fitted on all available patients. Every fit is strict:
## an error or invalid return terminates the repeat.
GMM_MODEL_NAMES <- c("EII", "VII", "EEI", "VVI", "EEE")
GMM_VERBOSE <- FALSE

## Ward.D2 is fitted on a bounded subsample to control O(n^2) memory.
WARD_SUBSAMPLE_N <- if (USE_TEST_MODE) 300L else 2000L

## Operating-characteristic thresholds.
TRUE_K_ARI_THRESHOLD <- 0.80
ASSIGNMENT_ARI_THRESHOLD <- 0.80
ALGORITHM_AGREEMENT_ARI_THRESHOLD <- 0.80
MICE_STABILITY_ARI_THRESHOLD <- 0.90
PREVALENCE_DRIFT_THRESHOLD <- 0.10
MATCHED_AGREEMENT_THRESHOLD <- 0.80

## Output controls.
SAVE_ALL_IMPUTATION_LABELS_RDS <- TRUE
SAVE_PRIMARY_IMPUTATION_LABELS_CSV_GZ <- TRUE
SAVE_CLUSTER_COMPOSITION <- TRUE
RESUME_COMPLETED_TASKS <- TRUE

## A completed marker is written only when all tasks succeed.
ALLOW_PARTIAL_OUTPUTS <- TRUE

RUN_MODE <- if (USE_TEST_MODE) "test" else "formal"

OUT_DIR <- file.path(SIM_ROOT, paste0("output_", RUN_MODE))
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
CHECKPOINT_DIR <- file.path(OUT_DIR, "task_checkpoints")
LABEL_RDS_DIR <- file.path(OUT_DIR, "individual_labels_all_imputations_rds")
LABEL_PRIMARY_DIR <- file.path(
  OUT_DIR,
  "individual_labels_primary_imputation_csv_gz"
)

for (d in c(
  SIM_ROOT, OUT_DIR, TAB_DIR, FIG_DIR, LOG_DIR,
  CHECKPOINT_DIR, LABEL_RDS_DIR, LABEL_PRIMARY_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  LOG_DIR,
  paste0("sensitivity07_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

## A configuration signature prevents accidental reuse of stale checkpoints.
CONFIG_SIGNATURE <- paste(
  ANALYSIS_ID,
  RUN_MODE,
  paste(SCENARIOS, collapse = "|"),
  N_ACTIVE_REP,
  N_DERIVATION,
  N_EXTERNAL,
  MICE_M,
  MICE_MAXIT,
  paste(K_VALUES, collapse = "|"),
  paste(GMM_MODEL_NAMES, collapse = "|"),
  WARD_SUBSAMPLE_N,
  GLOBAL_SEED,
  sep = "__"
)

## =============================================================================
## 1. Package checks and logging
## =============================================================================

required_pkgs <- c(
  "data.table",
  "mice",
  "mclust",
  "cluster",
  "future",
  "future.apply",
  "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(mclust)
  library(cluster)
  library(future)
  library(future.apply)
  library(ggplot2)
})

options(future.globals.maxSize = 12 * 1024^3)
data.table::setDTthreads(1L)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  invisible(msg)
}

log_msg("Combined sensitivity 07/07b started.")
log_msg("Run mode: ", RUN_MODE)
log_msg("Scenarios: ", paste(SCENARIOS, collapse = ", "))
log_msg("Repeats per scenario: ", N_ACTIVE_REP)
log_msg("N derivation/external: ", N_DERIVATION, "/", N_EXTERNAL)
log_msg("MICE m/maxit: ", MICE_M, "/", MICE_MAXIT)
log_msg("K values: ", paste(K_VALUES, collapse = ", "))
log_msg("GMM modelNames: ", paste(GMM_MODEL_NAMES, collapse = ", "))
log_msg("Ward subsample N: ", WARD_SUBSAMPLE_N)
log_msg("Workers: ", N_WORKERS)
log_msg("Configuration signature: ", CONFIG_SIGNATURE)

## =============================================================================
## 2. General utilities
## =============================================================================

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_quantile <- function(x, p) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(
    x,
    probs = p,
    names = FALSE,
    na.rm = TRUE,
    type = 7
  ))
}

safe_int <- function(x, default = NA_integer_) {
  if (!length(x) || all(is.na(x))) return(default)
  as.integer(x[1L])
}

safe_num <- function(x, default = NA_real_) {
  if (!length(x) || all(is.na(x))) return(default)
  as.numeric(x[1L])
}

safe_rbindlist <- function(x) {
  x <- Filter(
    function(z) {
      !is.null(z) &&
        (is.data.frame(z) || data.table::is.data.table(z)) &&
        nrow(z) > 0L
    },
    x
  )
  if (!length(x)) return(data.table())
  x <- lapply(x, function(z) {
    z <- as.data.table(z)
    if (anyDuplicated(names(z))) {
      setnames(z, make.unique(names(z), sep = "_dup"))
    }
    z
  })
  rbindlist(x, fill = TRUE, use.names = TRUE)
}

summarise_numeric <- function(dt, group_cols, value_col, metric_name) {
  if (!nrow(dt) || !value_col %in% names(dt)) return(data.table())
  dt[, {
    z <- as.numeric(get(value_col))
    z <- z[is.finite(z)]
    n <- length(z)
    list(
      metric = metric_name,
      n = n,
      mean = if (n) mean(z) else NA_real_,
      sd = if (n > 1L) sd(z) else NA_real_,
      mcse = if (n > 1L) sd(z) / sqrt(n) else NA_real_,
      median = if (n) median(z) else NA_real_,
      empirical_p025 = if (n) safe_quantile(z, 0.025) else NA_real_,
      empirical_p975 = if (n) safe_quantile(z, 0.975) else NA_real_
    )
  }, by = group_cols]
}

summarise_binary <- function(dt, group_cols, value_col, metric_name) {
  if (!nrow(dt) || !value_col %in% names(dt)) return(data.table())
  dt[, {
    z <- as.numeric(get(value_col))
    z <- z[is.finite(z)]
    n <- length(z)
    p <- if (n) mean(z) else NA_real_
    ci <- if (n) {
      stats::binom.test(sum(z), n)$conf.int
    } else {
      c(NA_real_, NA_real_)
    }
    list(
      metric = metric_name,
      n = n,
      proportion = p,
      mcse = if (n) sqrt(p * (1 - p) / n) else NA_real_,
      exact_low = ci[1L],
      exact_high = ci[2L]
    )
  }, by = group_cols]
}

all_permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
  do.call(
    rbind,
    lapply(v, function(x) {
      rest <- all_permutations(v[v != x])
      cbind(x, rest)
    })
  )
}

best_label_match <- function(reference, candidate, k) {
  if (length(reference) != length(candidate) ||
      anyNA(reference) || anyNA(candidate)) {
    return(list(
      agreement = NA_real_,
      mapped = rep(NA_integer_, length(candidate)),
      map = rep(NA_integer_, k)
    ))
  }

  perms <- all_permutations(seq_len(k))
  agreement <- apply(perms, 1L, function(p) {
    mean(reference == p[candidate])
  })
  best <- which.max(agreement)

  list(
    agreement = agreement[best],
    mapped = as.integer(perms[best, candidate]),
    map = as.integer(perms[best, ])
  )
}

apply_label_map <- function(labels, map) {
  out <- rep(NA_integer_, length(labels))
  ok <- is.finite(labels) & labels >= 1L & labels <= length(map)
  out[ok] <- as.integer(map[labels[ok]])
  out
}

adjusted_rand <- function(a, b) {
  if (length(a) != length(b) || !length(a) || anyNA(a) || anyNA(b)) {
    return(NA_real_)
  }
  as.numeric(mclust::adjustedRandIndex(a, b))
}

normalized_mutual_information <- function(a, b) {
  if (length(a) != length(b) || !length(a) || anyNA(a) || anyNA(b)) {
    return(NA_real_)
  }

  tab <- table(a, b)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)

  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1L]
    c <- nz[i, 2L]
    pxy[r, c] * log(pxy[r, c] / (px[r] * py[c]))
  }, numeric(1)))

  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))

  if (!is.finite(hx) || !is.finite(hy) || hx <= 0 || hy <= 0) {
    return(NA_real_)
  }

  mi / sqrt(hx * hy)
}

cluster_purity <- function(labels, truth) {
  if (length(labels) != length(truth) || anyNA(labels) || anyNA(truth)) {
    return(NA_real_)
  }
  tab <- table(labels, truth)
  sum(apply(tab, 1L, max)) / sum(tab)
}

cluster_prevalence <- function(labels, k) {
  tab <- table(factor(labels, levels = seq_len(k)))
  as.numeric(tab) / sum(tab)
}

centroids_from_labels <- function(X, labels, k) {
  X <- as.matrix(X)
  out <- matrix(NA_real_, nrow = k, ncol = ncol(X))
  colnames(out) <- colnames(X)

  for (g in seq_len(k)) {
    idx <- which(labels == g)
    if (length(idx)) {
      out[g, ] <- colMeans(X[idx, , drop = FALSE])
    }
  }

  out
}

profile_correlation_after_match <- function(
  derivation_centers,
  external_centers,
  map_external_to_derivation
) {
  k <- nrow(derivation_centers)
  mapped_external <- matrix(
    NA_real_,
    nrow = k,
    ncol = ncol(external_centers)
  )

  for (old in seq_len(k)) {
    new <- map_external_to_derivation[old]
    mapped_external[new, ] <- external_centers[old, ]
  }

  cors <- vapply(seq_len(k), function(g) {
    suppressWarnings(stats::cor(
      derivation_centers[g, ],
      mapped_external[g, ],
      use = "pairwise.complete.obs"
    ))
  }, numeric(1))

  safe_mean(cors)
}

write_session_info <- function(path) {
  zz <- file(path, open = "wt", encoding = "UTF-8")
  ## Revert the sink BEFORE closing the connection. Closing an active sink
  ## target first leaves sink() pointing at an invalid connection, which throws
  ## "invalid connection" on exit.
  on.exit({
    sink()
    close(zz)
  }, add = TRUE)
  sink(zz)
  print(sessionInfo())
}

save_plot_both <- function(p, stem, width = 8, height = 5) {
  ggsave(
    file.path(FIG_DIR, paste0(stem, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
  ggsave(
    file.path(FIG_DIR, paste0(stem, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

theme_publication <- function() {
  theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

## =============================================================================
## 3. Data-generating mechanisms
## =============================================================================

feature_names <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min",
  "abs_lymphocytes_min", "lactate_max", "ph_min", "pco2_max",
  "calcium_min", "calcium_max", "ptt_max", "inr_max",
  "temperature_min", "temperature_max", "urine_output_24h_ml",
  "glucose_max", "aniongap_max", "potassium_min", "potassium_max",
  "hemoglobin_min", "sodium_min", "sodium_max", "wbc_max",
  "platelets_min", "bicarbonate_min", "chloride_min", "chloride_max",
  "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)

baseline_vars <- c("age", "sex", "sofa", "aki_stage")
imputation_vars <- c(feature_names, baseline_vars)

severity_loading <- c(
  bilirubin_total_max = 0.45,
  alt_max = 0.25,
  pao2fio2ratio_min = -0.65,
  abs_lymphocytes_min = -0.35,
  lactate_max = 0.85,
  ph_min = -0.75,
  pco2_max = 0.30,
  calcium_min = -0.20,
  calcium_max = 0.10,
  ptt_max = 0.40,
  inr_max = 0.45,
  temperature_min = -0.15,
  temperature_max = 0.20,
  urine_output_24h_ml = -0.95,
  glucose_max = 0.25,
  aniongap_max = 0.80,
  potassium_min = -0.10,
  potassium_max = 0.35,
  hemoglobin_min = -0.25,
  sodium_min = -0.10,
  sodium_max = 0.10,
  wbc_max = 0.35,
  platelets_min = -0.40,
  bicarbonate_min = -0.85,
  chloride_min = -0.10,
  chloride_max = 0.15,
  bun_max = 0.90,
  creatinine_max = 0.95,
  resp_rate_max = 0.50,
  gcs_min = -0.60,
  spo2_min = -0.35,
  heart_rate_max = 0.55,
  mbp_min = -0.70
)

stopifnot(
  length(feature_names) == 33L,
  identical(names(severity_loading), feature_names)
)

make_subtype_patterns <- function() {
  p <- length(feature_names)
  out <- matrix(0, nrow = 3L, ncol = p)
  colnames(out) <- feature_names
  rownames(out) <- paste0("Subtype", 1:3)

  out[1L, c(
    "lactate_max", "aniongap_max", "heart_rate_max",
    "resp_rate_max", "wbc_max", "temperature_max"
  )] <- 1
  out[1L, c(
    "mbp_min", "bicarbonate_min", "ph_min", "platelets_min"
  )] <- -1

  out[2L, c(
    "creatinine_max", "bun_max", "potassium_max",
    "inr_max", "ptt_max"
  )] <- 1
  out[2L, c(
    "urine_output_24h_ml", "calcium_min", "hemoglobin_min"
  )] <- -1

  out[3L, c(
    "bilirubin_total_max", "alt_max", "pco2_max", "chloride_max"
  )] <- 1
  out[3L, c(
    "pao2fio2ratio_min", "spo2_min", "gcs_min",
    "abs_lymphocytes_min"
  )] <- -1

  out <- out - matrix(
    colMeans(out),
    nrow = 3L,
    ncol = p,
    byrow = TRUE
  )

  row_norm <- sqrt(rowMeans(out^2))
  out / row_norm
}

subtype_patterns <- make_subtype_patterns()

set.seed(91001L)
factor_loadings <- matrix(
  rnorm(length(feature_names) * 5L, sd = 0.25),
  nrow = length(feature_names),
  ncol = 5L
)

scenario_parameters <- function(scenario) {
  switch(
    scenario,
    S1_pure_severity_continuum = list(
      delta = 0,
      severity_multiplier = 1.00,
      noise_sd = 0.75,
      derivation_prob = c(1/3, 1/3, 1/3),
      external_prob = c(1/3, 1/3, 1/3),
      external_severity_mean = 0,
      external_feature_drift = FALSE,
      external_missing_extra = 0
    ),
    S2_true_discrete_subtypes = list(
      delta = 3.20,
      severity_multiplier = 0.35,
      noise_sd = 0.45,
      derivation_prob = c(1/3, 1/3, 1/3),
      external_prob = c(1/3, 1/3, 1/3),
      external_severity_mean = 0,
      external_feature_drift = FALSE,
      external_missing_extra = 0
    ),
    S3_transport_drift = list(
      delta = 3.20,
      severity_multiplier = 0.35,
      noise_sd = 0.45,
      derivation_prob = c(1/3, 1/3, 1/3),
      external_prob = c(0.55, 0.30, 0.15),
      external_severity_mean = 0.35,
      external_feature_drift = TRUE,
      external_missing_extra = 0.05
    ),
    stop("Unknown scenario: ", scenario, call. = FALSE)
  )
}

generate_cohort <- function(
  n,
  patient_prefix,
  delta,
  subtype_prob,
  severity_mean,
  severity_multiplier,
  noise_sd,
  external_feature_drift,
  seed
) {
  set.seed(seed)

  latent_subtype <- sample.int(
    3L,
    n,
    replace = TRUE,
    prob = subtype_prob
  )

  severity_true <- rnorm(
    n,
    mean = severity_mean,
    sd = 1
  )

  latent_factors <- matrix(
    rnorm(n * ncol(factor_loadings)),
    nrow = n
  )

  correlated_noise <- latent_factors %*% t(factor_loadings)
  independent_noise <- matrix(
    rnorm(
      n * length(feature_names),
      sd = noise_sd
    ),
    nrow = n
  )

  X <- severity_multiplier *
    (severity_true %o% severity_loading) +
    correlated_noise +
    independent_noise +
    delta * subtype_patterns[
      latent_subtype,
      ,
      drop = FALSE
    ]

  if (isTRUE(external_feature_drift)) {
    shift <- setNames(rep(0, length(feature_names)), feature_names)
    shift[c(
      "lactate_max", "creatinine_max", "bun_max",
      "resp_rate_max", "heart_rate_max"
    )] <- c(0.35, 0.30, 0.30, 0.20, 0.20)

    shift[c(
      "mbp_min", "pao2fio2ratio_min", "urine_output_24h_ml"
    )] <- c(-0.30, -0.35, -0.25)

    X <- sweep(X, 2L, shift, FUN = "+")

    scale_mult <- setNames(rep(1, length(feature_names)), feature_names)
    scale_mult[c(
      "lactate_max", "creatinine_max",
      "pao2fio2ratio_min", "urine_output_24h_ml"
    )] <- c(1.20, 1.15, 1.20, 1.15)

    X <- sweep(X, 2L, scale_mult, FUN = "*")
  }

  age <- pmin(
    95,
    pmax(
      18,
      round(
        65 +
          11 * rnorm(n) +
          2.0 * severity_true
      )
    )
  )

  sex <- rbinom(
    n,
    1,
    plogis(-0.05 + 0.08 * severity_true)
  )

  sofa <- pmin(
    24,
    pmax(
      0,
      round(
        6 +
          2.8 * severity_true +
          rnorm(n, sd = 2)
      )
    )
  )

  aki_stage <- pmin(
    3L,
    pmax(
      0L,
      as.integer(cut(
        severity_true + rnorm(n, sd = 0.8),
        breaks = c(-Inf, -0.6, 0.2, 1.0, Inf),
        labels = FALSE
      )) - 1L
    )
  )

  lp_mortality <- -2.35 +
    0.70 * severity_true +
    0.018 * (age - 65) +
    0.08 * sex +
    0.065 * sofa +
    0.16 * aki_stage

  mortality_30d <- rbinom(
    n,
    1,
    plogis(lp_mortality)
  )

  dt <- as.data.table(X)
  setnames(dt, feature_names)

  dt[, patient_id := sprintf(
    "%s%06d",
    patient_prefix,
    seq_len(.N)
  )]
  dt[, latent_subtype := latent_subtype]
  dt[, severity_true := severity_true]
  dt[, age := age]
  dt[, sex := sex]
  dt[, sofa := sofa]
  dt[, aki_stage := aki_stage]
  dt[, mortality_30d := mortality_30d]

  setcolorder(
    dt,
    c(
      "patient_id",
      "latent_subtype",
      "severity_true",
      baseline_vars,
      feature_names,
      "mortality_30d"
    )
  )

  dt
}

inject_mar_missingness <- function(
  dt,
  seed,
  external_extra = 0
) {
  set.seed(seed)
  out <- copy(dt)

  p <- length(feature_names)
  base_rate <- seq(0.05, 0.28, length.out = p)
  base_rate <- base_rate[sample.int(p)]

  severity_z <- as.numeric(scale(out$severity_true))
  age_z <- as.numeric(scale(out$age))

  for (j in seq_along(feature_names)) {
    v <- feature_names[j]
    xz <- as.numeric(scale(out[[v]]))
    xz[!is.finite(xz)] <- 0

    lin <- qlogis(
      pmin(0.80, base_rate[j] + external_extra)
    ) +
      0.30 * severity_z +
      0.12 * age_z +
      0.10 * xz

    probability <- pmin(
      0.85,
      pmax(0.005, plogis(lin))
    )

    miss <- rbinom(
      nrow(out),
      1,
      probability
    ) == 1L

    set(
      out,
      i = which(miss),
      j = v,
      value = NA_real_
    )
  }

  for (v in c("sofa", "aki_stage")) {
    probability <- pmin(
      0.20,
      pmax(
        0.01,
        0.03 +
          external_extra / 2 +
          0.03 * plogis(severity_z)
      )
    )

    miss <- rbinom(
      nrow(out),
      1,
      probability
    ) == 1L

    set(
      out,
      i = which(miss),
      j = v,
      value = NA_real_
    )
  }

  out
}

mice_complete_sets <- function(dt, m, maxit, seed) {
  imp_cols <- intersect(imputation_vars, names(dt))
  imp_data <- as.data.frame(dt[, ..imp_cols])

  method <- rep("pmm", length(imp_cols))
  names(method) <- imp_cols

  no_missing <- vapply(
    imp_data,
    function(x) !anyNA(x),
    logical(1)
  )
  method[no_missing] <- ""

  predictor_matrix <- mice::make.predictorMatrix(imp_data)
  diag(predictor_matrix) <- 0

  fit <- mice::mice(
    imp_data,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = predictor_matrix,
    seed = seed,
    printFlag = FALSE
  )

  lapply(seq_len(m), function(i) {
    completed <- mice::complete(fit, action = i)
    out <- copy(dt)

    for (v in imp_cols) {
      set(
        out,
        j = v,
        value = completed[[v]]
      )
    }

    out[, imputation_id := i]
    out
  })
}

## =============================================================================
## 4. Preprocessing
## =============================================================================

fit_preprocess <- function(
  dt,
  features,
  probs = WINSOR_PROBS
) {
  bounds <- lapply(features, function(v) {
    as.numeric(stats::quantile(
      dt[[v]],
      probs = probs,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    ))
  })
  names(bounds) <- features

  Xw <- sapply(features, function(v) {
    b <- bounds[[v]]
    pmin(
      pmax(as.numeric(dt[[v]]), b[1L]),
      b[2L]
    )
  })

  Xw <- as.matrix(Xw)
  colnames(Xw) <- features

  center <- colMeans(Xw)
  scale_value <- apply(Xw, 2L, stats::sd)
  scale_value[
    !is.finite(scale_value) |
      scale_value < 1e-8
  ] <- 1

  list(
    bounds = bounds,
    center = center,
    scale = scale_value,
    features = features
  )
}

apply_preprocess <- function(dt, prep) {
  Xw <- sapply(prep$features, function(v) {
    b <- prep$bounds[[v]]
    pmin(
      pmax(as.numeric(dt[[v]]), b[1L]),
      b[2L]
    )
  })

  Xw <- as.matrix(Xw)
  colnames(Xw) <- prep$features

  Xz <- sweep(
    Xw,
    2L,
    prep$center,
    FUN = "-"
  )

  Xz <- sweep(
    Xz,
    2L,
    prep$scale,
    FUN = "/"
  )

  Xz[!is.finite(Xz)] <- 0
  Xz
}

preprocess_to_long <- function(
  prep,
  task_id,
  scenario,
  repeat_id,
  imputation_id,
  cohort_basis
) {
  rbindlist(lapply(prep$features, function(v) {
    data.table(
      task_id = task_id,
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_id = imputation_id,
      cohort_basis = cohort_basis,
      feature = v,
      winsor_low = prep$bounds[[v]][1L],
      winsor_high = prep$bounds[[v]][2L],
      center = prep$center[[v]],
      scale = prep$scale[[v]]
    )
  }))
}

## =============================================================================
## 5. Strict GMM and Ward fitting
## =============================================================================

extract_gmm_centers <- function(fit, k, feature_names) {
  mu <- fit$parameters$mean

  if (is.null(mu)) {
    stop("GMM parameters$mean is NULL.", call. = FALSE)
  }

  if (is.vector(mu)) {
    mu <- matrix(mu, ncol = 1L)
  }

  mu <- as.matrix(mu)

  if (nrow(mu) == length(feature_names) &&
      ncol(mu) == k) {
    centers <- t(mu)
  } else if (
    nrow(mu) == k &&
      ncol(mu) == length(feature_names)
  ) {
    centers <- mu
  } else {
    stop(
      "Unexpected GMM mean dimensions: ",
      paste(dim(mu), collapse = " x "),
      "; expected ",
      length(feature_names), " x ", k,
      " or ", k, " x ", length(feature_names),
      call. = FALSE
    )
  }

  colnames(centers) <- feature_names
  centers
}

strict_gmm_fit <- function(
  X,
  k,
  seed,
  stage_label
) {
  X <- as.matrix(X)

  if (nrow(X) <= k ||
      ncol(X) < 1L ||
      any(!is.finite(X))) {
    stop(
      "GMM invalid input at ", stage_label,
      call. = FALSE
    )
  }

  set.seed(seed)

  warning_messages <- character()

  fit <- tryCatch(
    withCallingHandlers(
      mclust::Mclust(
        data = X,
        G = k,
        modelNames = GMM_MODEL_NAMES,
        verbose = GMM_VERBOSE
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      stop(
        "GMM failed at ", stage_label,
        ": ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.null(fit) ||
      !inherits(fit, "Mclust") ||
      safe_int(fit$G) != k ||
      is.null(fit$classification) ||
      length(fit$classification) != nrow(X) ||
      anyNA(fit$classification) ||
      length(unique(fit$classification)) != k) {
    stop(
      "GMM returned an invalid fit at ", stage_label,
      ". No fallback is permitted.",
      call. = FALSE
    )
  }

  centers <- extract_gmm_centers(
    fit,
    k = k,
    feature_names = colnames(X)
  )

  posterior <- fit$z
  max_posterior <- if (!is.null(posterior)) {
    apply(as.matrix(posterior), 1L, max)
  } else {
    rep(NA_real_, nrow(X))
  }

  list(
    algorithm = "gmm",
    k = k,
    labels_raw = as.integer(fit$classification),
    centers = centers,
    fit_object = fit,
    model_name = as.character(fit$modelName)[1L],
    bic = safe_num(fit$bic),
    log_likelihood = safe_num(fit$loglik),
    entropy = {
      z <- as.matrix(fit$z)
      if (!is.null(z) && nrow(z)) {
        -sum(z * log(pmax(z, 1e-12))) /
          (nrow(z) * log(k))
      } else {
        NA_real_
      }
    },
    mean_max_posterior = safe_mean(max_posterior),
    min_max_posterior = suppressWarnings(min(
      max_posterior,
      na.rm = TRUE
    )),
    warnings = paste(
      unique(warning_messages),
      collapse = " | "
    )
  )
}

strict_gmm_predict <- function(
  gmm_fit,
  X_new,
  stage_label
) {
  X_new <- as.matrix(X_new)

  pred <- tryCatch(
    predict(
      gmm_fit$fit_object,
      newdata = X_new
    ),
    error = function(e) {
      stop(
        "Frozen GMM prediction failed at ",
        stage_label,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  classification <- pred$classification

  if (is.null(classification) ||
      length(classification) != nrow(X_new) ||
      anyNA(classification)) {
    stop(
      "Frozen GMM prediction returned invalid labels at ",
      stage_label,
      call. = FALSE
    )
  }

  as.integer(classification)
}

ward_fit <- function(
  X,
  k,
  seed,
  subsample_n,
  stage_label
) {
  X <- as.matrix(X)

  n_sub <- min(
    nrow(X),
    as.integer(subsample_n)
  )

  if (n_sub <= k + 5L) {
    stop(
      "Insufficient Ward subsample at ",
      stage_label,
      call. = FALSE
    )
  }

  set.seed(seed)
  idx <- sort(sample.int(
    nrow(X),
    n_sub,
    replace = FALSE
  ))

  X_sub <- X[idx, , drop = FALSE]

  hc <- tryCatch(
    stats::hclust(
      stats::dist(X_sub),
      method = "ward.D2"
    ),
    error = function(e) {
      stop(
        "Ward.D2 failed at ",
        stage_label,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  labels_sub <- as.integer(stats::cutree(hc, k = k))
  centers <- centroids_from_labels(
    X_sub,
    labels_sub,
    k
  )

  if (any(!is.finite(centers))) {
    stop(
      "Ward centroids are non-finite at ",
      stage_label,
      call. = FALSE
    )
  }

  labels_all <- assign_nearest_centroid(
    X,
    centers
  )

  silhouette_mean <- tryCatch(
    safe_mean(
      cluster::silhouette(
        labels_sub,
        stats::dist(X_sub)
      )[, "sil_width"]
    ),
    error = function(e) NA_real_
  )

  list(
    algorithm = "ward",
    k = k,
    labels_raw = labels_all,
    centers = centers,
    subsample_index = idx,
    subsample_labels = labels_sub,
    silhouette_mean = silhouette_mean
  )
}

assign_nearest_centroid <- function(X, centers) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)

  distance_matrix <- sapply(
    seq_len(nrow(centers)),
    function(g) {
      rowSums(
        (
          X -
            matrix(
              centers[g, ],
              nrow = nrow(X),
              ncol = ncol(X),
              byrow = TRUE
            )
        )^2
      )
    }
  )

  if (is.null(dim(distance_matrix))) {
    distance_matrix <- matrix(
      distance_matrix,
      ncol = 1L
    )
  }

  as.integer(max.col(
    -distance_matrix,
    ties.method = "first"
  ))
}

## =============================================================================
## 6. Label alignment, composition, and metrics
## =============================================================================

severity_canonical_map_k2 <- function(centers) {
  loading <- severity_loading[colnames(centers)]
  burden_score <- as.numeric(centers %*% loading)

  high_old <- which.max(burden_score)
  low_old <- setdiff(1:2, high_old)

  map <- integer(2L)
  map[high_old] <- 1L
  map[low_old] <- 2L
  map
}

outcome_separation_k2 <- function(
  dt,
  labels_canonical,
  outcome = "mortality_30d"
) {
  if (length(labels_canonical) != nrow(dt) ||
      anyNA(labels_canonical)) {
    return(NA_real_)
  }

  y <- as.numeric(dt[[outcome]])

  safe_mean(y[labels_canonical == 1L]) -
    safe_mean(y[labels_canonical == 2L])
}

cluster_composition_table <- function(
  labels,
  truth,
  task_id,
  scenario,
  repeat_id,
  imputation_id,
  cohort,
  label_source,
  algorithm,
  k
) {
  dt <- data.table(
    cluster = as.integer(labels),
    latent_subtype = as.integer(truth)
  )

  out <- dt[, .N, by = .(
    cluster,
    latent_subtype
  )]

  out[, cluster_n := sum(N), by = cluster]
  out[, subtype_n := sum(N), by = latent_subtype]
  out[, total_n := sum(N)]
  out[, proportion_within_cluster := N / cluster_n]
  out[, proportion_within_subtype := N / subtype_n]

  out[, `:=`(
    task_id = task_id,
    scenario = scenario,
    repeat_id = repeat_id,
    imputation_id = imputation_id,
    cohort = cohort,
    label_source = label_source,
    algorithm = algorithm,
    k = k
  )]

  setcolorder(
    out,
    c(
      "task_id", "scenario", "repeat_id", "imputation_id",
      "cohort", "label_source", "algorithm", "k",
      "cluster", "latent_subtype",
      "N", "cluster_n", "subtype_n", "total_n",
      "proportion_within_cluster",
      "proportion_within_subtype"
    )
  )

  out
}

transport_metrics <- function(
  frozen_labels,
  reclustered_labels,
  derivation_labels,
  X_external_local,
  derivation_centers,
  reclustered_centers,
  k
) {
  match <- best_label_match(
    frozen_labels,
    reclustered_labels,
    k
  )

  mapped_reclustered <- match$mapped

  prev_derivation <- cluster_prevalence(
    derivation_labels,
    k
  )

  prev_external <- cluster_prevalence(
    mapped_reclustered,
    k
  )

  data.table(
    k = k,
    frozen_vs_reclustered_ari = adjusted_rand(
      frozen_labels,
      reclustered_labels
    ),
    matched_agreement = match$agreement,
    mean_abs_prevalence_difference = mean(
      abs(prev_external - prev_derivation)
    ),
    max_abs_prevalence_difference = max(
      abs(prev_external - prev_derivation)
    ),
    centroid_profile_correlation =
      profile_correlation_after_match(
        derivation_centers,
        reclustered_centers,
        match$map
      )
  )
}

## =============================================================================
## 7. Task-level analysis
## =============================================================================

run_task <- function(task_row) {
  task_id <- as.character(task_row$task_id)
  scenario <- as.character(task_row$scenario)
  repeat_id <- as.integer(task_row$repeat_id)
  task_seed <- as.integer(task_row$task_seed)

  pars <- scenario_parameters(scenario)

  derivation_raw <- generate_cohort(
    n = N_DERIVATION,
    patient_prefix = "D",
    delta = pars$delta,
    subtype_prob = pars$derivation_prob,
    severity_mean = 0,
    severity_multiplier = pars$severity_multiplier,
    noise_sd = pars$noise_sd,
    external_feature_drift = FALSE,
    seed = task_seed + 1L
  )

  external_raw <- generate_cohort(
    n = N_EXTERNAL,
    patient_prefix = "E",
    delta = pars$delta,
    subtype_prob = pars$external_prob,
    severity_mean = pars$external_severity_mean,
    severity_multiplier = pars$severity_multiplier,
    noise_sd = pars$noise_sd,
    external_feature_drift = pars$external_feature_drift,
    seed = task_seed + 2L
  )

  derivation_missing <- inject_mar_missingness(
    derivation_raw,
    seed = task_seed + 10L,
    external_extra = 0
  )

  external_missing <- inject_mar_missingness(
    external_raw,
    seed = task_seed + 11L,
    external_extra = pars$external_missing_extra
  )

  derivation_sets <- mice_complete_sets(
    derivation_missing,
    m = MICE_M,
    maxit = MICE_MAXIT,
    seed = task_seed + 100L
  )

  external_sets <- mice_complete_sets(
    external_missing,
    m = MICE_M,
    maxit = MICE_MAXIT,
    seed = task_seed + 200L
  )

  labels_derivation <- vector("list", MICE_M)
  labels_external <- vector("list", MICE_M)

  domain1_rows <- list()
  transport_rows <- list()
  algorithm_agreement_rows <- list()
  composition_rows <- list()
  model_rows <- list()
  centroid_rows <- list()
  preprocess_rows <- list()

  ## Store raw labels for post-hoc cross-imputation and paired alignment.
  raw_label_store <- vector("list", MICE_M)

  for (imp in seq_len(MICE_M)) {
    derivation_dt <- derivation_sets[[imp]]
    external_dt <- external_sets[[imp]]

    prep_derivation <- fit_preprocess(
      derivation_dt,
      feature_names
    )

    prep_external <- fit_preprocess(
      external_dt,
      feature_names
    )

    X_derivation <- apply_preprocess(
      derivation_dt,
      prep_derivation
    )

    X_external_frozen_space <- apply_preprocess(
      external_dt,
      prep_derivation
    )

    X_external_local <- apply_preprocess(
      external_dt,
      prep_external
    )

    preprocess_rows[[length(preprocess_rows) + 1L]] <-
      preprocess_to_long(
        prep_derivation,
        task_id,
        scenario,
        repeat_id,
        imp,
        cohort_basis = "derivation"
      )

    preprocess_rows[[length(preprocess_rows) + 1L]] <-
      preprocess_to_long(
        prep_external,
        task_id,
        scenario,
        repeat_id,
        imp,
        cohort_basis = "external_local"
      )

    raw_label_store[[imp]] <- list(
      derivation = list(),
      external_frozen = list(),
      external_reclustered = list(),
      fit_objects = list()
    )

    for (k in K_VALUES) {
      ## -----------------------------------------------------------------------
      ## Strict GMM.
      ## Every call is mclust::Mclust(...); no fallback is defined.
      ## -----------------------------------------------------------------------

      gmm_derivation <- strict_gmm_fit(
        X_derivation,
        k = k,
        seed = task_seed + 100000L + 1000L * k + imp,
        stage_label = paste0(
          task_id,
          "/imp", imp,
          "/derivation/GMM/K", k
        )
      )

      gmm_external_frozen <- strict_gmm_predict(
        gmm_derivation,
        X_external_frozen_space,
        stage_label = paste0(
          task_id,
          "/imp", imp,
          "/external_frozen/GMM/K", k
        )
      )

      gmm_external_reclustered <- strict_gmm_fit(
        X_external_local,
        k = k,
        seed = task_seed + 200000L + 1000L * k + imp,
        stage_label = paste0(
          task_id,
          "/imp", imp,
          "/external_reclustered/GMM/K", k
        )
      )

      ## -----------------------------------------------------------------------
      ## Ward.D2.
      ## -----------------------------------------------------------------------

      ward_derivation <- ward_fit(
        X_derivation,
        k = k,
        seed = task_seed + 300000L + 1000L * k + imp,
        subsample_n = WARD_SUBSAMPLE_N,
        stage_label = paste0(
          task_id,
          "/imp", imp,
          "/derivation/Ward/K", k
        )
      )

      ward_external_frozen <- assign_nearest_centroid(
        X_external_frozen_space,
        ward_derivation$centers
      )

      ward_external_reclustered <- ward_fit(
        X_external_local,
        k = k,
        seed = task_seed + 400000L + 1000L * k + imp,
        subsample_n = WARD_SUBSAMPLE_N,
        stage_label = paste0(
          task_id,
          "/imp", imp,
          "/external_reclustered/Ward/K", k
        )
      )

      raw_label_store[[imp]]$derivation$gmm[[as.character(k)]] <-
        gmm_derivation$labels_raw

      raw_label_store[[imp]]$derivation$ward[[as.character(k)]] <-
        ward_derivation$labels_raw

      raw_label_store[[imp]]$external_frozen$gmm[[as.character(k)]] <-
        gmm_external_frozen

      raw_label_store[[imp]]$external_frozen$ward[[as.character(k)]] <-
        ward_external_frozen

      raw_label_store[[imp]]$external_reclustered$gmm[[as.character(k)]] <-
        gmm_external_reclustered$labels_raw

      raw_label_store[[imp]]$external_reclustered$ward[[as.character(k)]] <-
        ward_external_reclustered$labels_raw

      raw_label_store[[imp]]$fit_objects$gmm_derivation[[as.character(k)]] <-
        gmm_derivation

      raw_label_store[[imp]]$fit_objects$gmm_external[[as.character(k)]] <-
        gmm_external_reclustered

      raw_label_store[[imp]]$fit_objects$ward_derivation[[as.character(k)]] <-
        ward_derivation

      raw_label_store[[imp]]$fit_objects$ward_external[[as.character(k)]] <-
        ward_external_reclustered

      ## Canonical K2 labels are defined by prespecified severity burden.
      gmm_derivation_canonical <- if (k == 2L) {
        apply_label_map(
          gmm_derivation$labels_raw,
          severity_canonical_map_k2(
            gmm_derivation$centers
          )
        )
      } else {
        gmm_derivation$labels_raw
      }

      ward_derivation_canonical <- if (k == 2L) {
        apply_label_map(
          ward_derivation$labels_raw,
          severity_canonical_map_k2(
            ward_derivation$centers
          )
        )
      } else {
        ward_derivation$labels_raw
      }

      for (algorithm in c("gmm", "ward")) {
        fit_derivation <- if (algorithm == "gmm") {
          gmm_derivation
        } else {
          ward_derivation
        }

        labels_derivation_algorithm <- if (algorithm == "gmm") {
          gmm_derivation$labels_raw
        } else {
          ward_derivation$labels_raw
        }

        labels_external_frozen_algorithm <- if (algorithm == "gmm") {
          gmm_external_frozen
        } else {
          ward_external_frozen
        }

        fit_external <- if (algorithm == "gmm") {
          gmm_external_reclustered
        } else {
          ward_external_reclustered
        }

        labels_external_reclustered_algorithm <- fit_external$labels_raw

        labels_canonical <- if (algorithm == "gmm") {
          gmm_derivation_canonical
        } else {
          ward_derivation_canonical
        }

        domain1_rows[[length(domain1_rows) + 1L]] <- data.table(
          task_id = task_id,
          scenario = scenario,
          repeat_id = repeat_id,
          imputation_id = imp,
          algorithm = algorithm,
          k = k,
          ari_vs_latent_subtype = adjusted_rand(
            labels_derivation_algorithm,
            derivation_dt$latent_subtype
          ),
          nmi_vs_latent_subtype = normalized_mutual_information(
            labels_derivation_algorithm,
            derivation_dt$latent_subtype
          ),
          purity_vs_latent_subtype = cluster_purity(
            labels_derivation_algorithm,
            derivation_dt$latent_subtype
          ),
          mortality_30d_C1_minus_C2 = if (k == 2L) {
            outcome_separation_k2(
              derivation_dt,
              labels_canonical,
              outcome = "mortality_30d"
            )
          } else {
            NA_real_
          },
          cluster_1_prevalence = mean(
            labels_canonical == 1L
          ),
          cluster_min_prevalence = min(
            cluster_prevalence(
              labels_derivation_algorithm,
              k
            )
          ),
          ward_silhouette_mean = if (algorithm == "ward") {
            ward_derivation$silhouette_mean
          } else {
            NA_real_
          },
          gmm_model_name = if (algorithm == "gmm") {
            gmm_derivation$model_name
          } else {
            NA_character_
          },
          gmm_bic = if (algorithm == "gmm") {
            gmm_derivation$bic
          } else {
            NA_real_
          },
          gmm_entropy = if (algorithm == "gmm") {
            gmm_derivation$entropy
          } else {
            NA_real_
          },
          gmm_mean_max_posterior = if (algorithm == "gmm") {
            gmm_derivation$mean_max_posterior
          } else {
            NA_real_
          },
          gmm_min_max_posterior = if (algorithm == "gmm") {
            gmm_derivation$min_max_posterior
          } else {
            NA_real_
          },
          gmm_warnings = if (algorithm == "gmm") {
            gmm_derivation$warnings
          } else {
            NA_character_
          }
        )

        tm <- transport_metrics(
          frozen_labels = labels_external_frozen_algorithm,
          reclustered_labels = labels_external_reclustered_algorithm,
          derivation_labels = labels_derivation_algorithm,
          X_external_local = X_external_local,
          derivation_centers = fit_derivation$centers,
          reclustered_centers = fit_external$centers,
          k = k
        )

        tm[, `:=`(
          task_id = task_id,
          scenario = scenario,
          repeat_id = repeat_id,
          imputation_id = imp,
          algorithm = algorithm,
          frozen_ari_vs_external_truth = adjusted_rand(
            labels_external_frozen_algorithm,
            external_dt$latent_subtype
          ),
          reclustered_ari_vs_external_truth = adjusted_rand(
            labels_external_reclustered_algorithm,
            external_dt$latent_subtype
          ),
          frozen_nmi_vs_external_truth = normalized_mutual_information(
            labels_external_frozen_algorithm,
            external_dt$latent_subtype
          ),
          reclustered_nmi_vs_external_truth =
            normalized_mutual_information(
              labels_external_reclustered_algorithm,
              external_dt$latent_subtype
            )
        )]

        transport_rows[[length(transport_rows) + 1L]] <- tm

        if (isTRUE(SAVE_CLUSTER_COMPOSITION)) {
          composition_rows[[length(composition_rows) + 1L]] <-
            cluster_composition_table(
              labels_derivation_algorithm,
              derivation_dt$latent_subtype,
              task_id,
              scenario,
              repeat_id,
              imp,
              cohort = "derivation",
              label_source = "fitted",
              algorithm = algorithm,
              k = k
            )

          composition_rows[[length(composition_rows) + 1L]] <-
            cluster_composition_table(
              labels_external_frozen_algorithm,
              external_dt$latent_subtype,
              task_id,
              scenario,
              repeat_id,
              imp,
              cohort = "external",
              label_source = "frozen",
              algorithm = algorithm,
              k = k
            )

          composition_rows[[length(composition_rows) + 1L]] <-
            cluster_composition_table(
              labels_external_reclustered_algorithm,
              external_dt$latent_subtype,
              task_id,
              scenario,
              repeat_id,
              imp,
              cohort = "external",
              label_source = "reclustered",
              algorithm = algorithm,
              k = k
            )
        }

        centers_long <- rbindlist(lapply(seq_len(k), function(g) {
          data.table(
            task_id = task_id,
            scenario = scenario,
            repeat_id = repeat_id,
            imputation_id = imp,
            algorithm = algorithm,
            cohort_basis = "derivation",
            k = k,
            cluster = g,
            feature = colnames(fit_derivation$centers),
            center = as.numeric(
              fit_derivation$centers[g, ]
            )
          )
        }))

        centroid_rows[[length(centroid_rows) + 1L]] <- centers_long

        external_centers_long <- rbindlist(
          lapply(seq_len(k), function(g) {
            data.table(
              task_id = task_id,
              scenario = scenario,
              repeat_id = repeat_id,
              imputation_id = imp,
              algorithm = algorithm,
              cohort_basis = "external_local_reclustered",
              k = k,
              cluster = g,
              feature = colnames(fit_external$centers),
              center = as.numeric(
                fit_external$centers[g, ]
              )
            )
          })
        )

        centroid_rows[[length(centroid_rows) + 1L]] <-
          external_centers_long
      }

      ## Same-patient GMM-vs-Ward agreement.
      for (label_source in c(
        "derivation",
        "external_frozen",
        "external_reclustered"
      )) {
        gmm_labels <- raw_label_store[[imp]][[label_source]]$gmm[[as.character(k)]]

        ward_labels <- raw_label_store[[imp]][[label_source]]$ward[[as.character(k)]]

        match_ward_to_gmm <- best_label_match(
          gmm_labels,
          ward_labels,
          k
        )

        algorithm_agreement_rows[[length(algorithm_agreement_rows) + 1L]] <- data.table(
          task_id = task_id,
          scenario = scenario,
          repeat_id = repeat_id,
          imputation_id = imp,
          cohort_label_source = label_source,
          k = k,
          gmm_vs_ward_ari = adjusted_rand(
            gmm_labels,
            ward_labels
          ),
          gmm_vs_ward_matched_agreement =
            match_ward_to_gmm$agreement
        )
      }

      ## Model metadata.
      model_rows[[length(model_rows) + 1L]] <- data.table(
        task_id = task_id,
        scenario = scenario,
        repeat_id = repeat_id,
        imputation_id = imp,
        cohort_basis = "derivation",
        algorithm = "gmm",
        k = k,
        model_name = gmm_derivation$model_name,
        bic = gmm_derivation$bic,
        log_likelihood = gmm_derivation$log_likelihood,
        entropy = gmm_derivation$entropy,
        mean_max_posterior =
          gmm_derivation$mean_max_posterior,
        min_max_posterior =
          gmm_derivation$min_max_posterior,
        warnings = gmm_derivation$warnings
      )

      model_rows[[length(model_rows) + 1L]] <- data.table(
        task_id = task_id,
        scenario = scenario,
        repeat_id = repeat_id,
        imputation_id = imp,
        cohort_basis = "external_local_reclustered",
        algorithm = "gmm",
        k = k,
        model_name = gmm_external_reclustered$model_name,
        bic = gmm_external_reclustered$bic,
        log_likelihood =
          gmm_external_reclustered$log_likelihood,
        entropy = gmm_external_reclustered$entropy,
        mean_max_posterior =
          gmm_external_reclustered$mean_max_posterior,
        min_max_posterior =
          gmm_external_reclustered$min_max_posterior,
        warnings = gmm_external_reclustered$warnings
      )

      model_rows[[length(model_rows) + 1L]] <- data.table(
        task_id = task_id,
        scenario = scenario,
        repeat_id = repeat_id,
        imputation_id = imp,
        cohort_basis = "derivation",
        algorithm = "ward",
        k = k,
        model_name = "ward.D2",
        bic = NA_real_,
        log_likelihood = NA_real_,
        entropy = NA_real_,
        mean_max_posterior = NA_real_,
        min_max_posterior = NA_real_,
        warnings = NA_character_,
        ward_silhouette_mean =
          ward_derivation$silhouette_mean
      )

      model_rows[[length(model_rows) + 1L]] <- data.table(
        task_id = task_id,
        scenario = scenario,
        repeat_id = repeat_id,
        imputation_id = imp,
        cohort_basis = "external_local_reclustered",
        algorithm = "ward",
        k = k,
        model_name = "ward.D2",
        bic = NA_real_,
        log_likelihood = NA_real_,
        entropy = NA_real_,
        mean_max_posterior = NA_real_,
        min_max_posterior = NA_real_,
        warnings = NA_character_,
        ward_silhouette_mean =
          ward_external_reclustered$silhouette_mean
      )
    }
  }

  ## ---------------------------------------------------------------------------
  ## Individual-label tables.
  ## ---------------------------------------------------------------------------

  primary_raw <- raw_label_store[[PRIMARY_IMPUTATION_ID]]

  for (imp in seq_len(MICE_M)) {
    derivation_dt <- derivation_sets[[imp]]
    external_dt <- external_sets[[imp]]

    dlab <- data.table(
      task_id = task_id,
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_id = imp,
      cohort = "derivation",
      patient_id = derivation_dt$patient_id,
      latent_subtype = derivation_dt$latent_subtype,
      severity_true = derivation_dt$severity_true,
      mortality_30d = derivation_dt$mortality_30d
    )

    elab <- data.table(
      task_id = task_id,
      scenario = scenario,
      repeat_id = repeat_id,
      imputation_id = imp,
      cohort = "external",
      patient_id = external_dt$patient_id,
      latent_subtype = external_dt$latent_subtype,
      severity_true = external_dt$severity_true,
      mortality_30d = external_dt$mortality_30d
    )

    for (k in K_VALUES) {
      k_chr <- as.character(k)

      for (algorithm in c("gmm", "ward")) {
        deriv_raw <- raw_label_store[[imp]]$derivation[[algorithm]][[k_chr]]

        external_frozen_raw <- raw_label_store[[imp]]$external_frozen[[algorithm]][[k_chr]]

        external_reclustered_raw <- raw_label_store[[imp]]$external_reclustered[[algorithm]][[k_chr]]

        primary_deriv_raw <- primary_raw$derivation[[algorithm]][[k_chr]]

        primary_external_reclustered_raw <-
          primary_raw$external_reclustered[[algorithm]][[k_chr]]

        map_to_primary_derivation <- best_label_match(
          primary_deriv_raw,
          deriv_raw,
          k
        )$map

        map_recluster_to_frozen <- best_label_match(
          external_frozen_raw,
          external_reclustered_raw,
          k
        )$map

        map_external_recluster_to_primary <- best_label_match(
          primary_external_reclustered_raw,
          external_reclustered_raw,
          k
        )$map

        dlab[, (
          paste0(algorithm, "_k", k, "_raw")
        ) := deriv_raw]

        dlab[, (
          paste0(
            algorithm,
            "_k", k,
            "_aligned_to_primary_imputation"
          )
        ) := apply_label_map(
          deriv_raw,
          map_to_primary_derivation
        )]

        elab[, (
          paste0(
            algorithm,
            "_k", k,
            "_frozen_raw"
          )
        ) := external_frozen_raw]

        elab[, (
          paste0(
            algorithm,
            "_k", k,
            "_frozen_aligned_to_primary_imputation"
          )
        ) := apply_label_map(
          external_frozen_raw,
          map_to_primary_derivation
        )]

        elab[, (
          paste0(
            algorithm,
            "_k", k,
            "_reclustered_raw"
          )
        ) := external_reclustered_raw]

        elab[, (
          paste0(
            algorithm,
            "_k", k,
            "_reclustered_matched_to_frozen"
          )
        ) := apply_label_map(
          external_reclustered_raw,
          map_recluster_to_frozen
        )]

        elab[, (
          paste0(
            algorithm,
            "_k", k,
            "_reclustered_aligned_to_primary_imputation"
          )
        ) := apply_label_map(
          external_reclustered_raw,
          map_external_recluster_to_primary
        )]
      }

      ## Ward matched to GMM on the same patients.
      deriv_ward_to_gmm_map <- best_label_match(
        raw_label_store[[imp]]$derivation$gmm[[k_chr]],
        raw_label_store[[imp]]$derivation$ward[[k_chr]],
        k
      )$map

      external_recluster_ward_to_gmm_map <- best_label_match(
        raw_label_store[[imp]]$external_reclustered$gmm[[k_chr]],
        raw_label_store[[imp]]$external_reclustered$ward[[k_chr]],
        k
      )$map

      dlab[, (
        paste0(
          "ward_k", k,
          "_matched_to_gmm_same_dataset"
        )
      ) := apply_label_map(
        raw_label_store[[imp]]$derivation$ward[[k_chr]],
        deriv_ward_to_gmm_map
      )]

      elab[, (
        paste0(
          "ward_k", k,
          "_frozen_matched_to_gmm_same_dataset"
        )
      ) := apply_label_map(
        raw_label_store[[imp]]$external_frozen$ward[[k_chr]],
        deriv_ward_to_gmm_map
      )]

      elab[, (
        paste0(
          "ward_k", k,
          "_reclustered_matched_to_gmm_same_dataset"
        )
      ) := apply_label_map(
        raw_label_store[[imp]]$external_reclustered$ward[[k_chr]],
        external_recluster_ward_to_gmm_map
      )]
    }

    labels_derivation[[imp]] <- dlab
    labels_external[[imp]] <- elab
  }

  labels_all <- rbindlist(
    c(labels_derivation, labels_external),
    fill = TRUE,
    use.names = TRUE
  )

  label_rds_path <- file.path(
    LABEL_RDS_DIR,
    paste0(task_id, "_individual_labels_all_imputations.rds")
  )

  primary_label_csv_path <- file.path(
    LABEL_PRIMARY_DIR,
    paste0(
      task_id,
      "_individual_labels_primary_imputation.csv.gz"
    )
  )

  if (isTRUE(SAVE_ALL_IMPUTATION_LABELS_RDS)) {
    saveRDS(
      labels_all,
      label_rds_path,
      compress = "xz"
    )
  }

  if (isTRUE(SAVE_PRIMARY_IMPUTATION_LABELS_CSV_GZ)) {
    fwrite(
      labels_all[
        imputation_id == PRIMARY_IMPUTATION_ID
      ],
      primary_label_csv_path
    )
  }

  ## ---------------------------------------------------------------------------
  ## Cross-imputation stability.
  ## ---------------------------------------------------------------------------

  stability_rows <- list()

  for (algorithm in c("gmm", "ward")) {
    for (k in K_VALUES) {
      k_chr <- as.character(k)

      for (label_source in c(
        "derivation",
        "external_frozen",
        "external_reclustered"
      )) {
        reference <- raw_label_store[[PRIMARY_IMPUTATION_ID]][[label_source]][[algorithm]][[k_chr]]

        ## Exclude the trivial primary-vs-itself comparison (ARI = 1), which
        ## otherwise inflates mean stability and the >= threshold proportion.
        compare_ids <- setdiff(seq_len(MICE_M), PRIMARY_IMPUTATION_ID)

        ari_values <- if (length(compare_ids)) {
          vapply(
            compare_ids,
            function(imp) {
              candidate <- raw_label_store[[imp]][[label_source]][[algorithm]][[k_chr]]

              adjusted_rand(
                reference,
                candidate
              )
            },
            numeric(1)
          )
        } else {
          NA_real_
        }

        finite_ari <- ari_values[is.finite(ari_values)]

        stability_rows[[length(stability_rows) + 1L]] <-
          data.table(
            task_id = task_id,
            scenario = scenario,
            repeat_id = repeat_id,
            algorithm = algorithm,
            k = k,
            label_source = label_source,
            m = MICE_M,
            n_stability_comparisons = length(compare_ids),
            ari_mean_vs_primary = if (length(finite_ari)) {
              mean(finite_ari)
            } else {
              NA_real_
            },
            ari_median_vs_primary = if (length(finite_ari)) {
              median(finite_ari)
            } else {
              NA_real_
            },
            ari_min_vs_primary = if (length(finite_ari)) {
              min(finite_ari)
            } else {
              NA_real_
            },
            ari_prop_ge_threshold = if (length(finite_ari)) {
              mean(finite_ari >= MICE_STABILITY_ARI_THRESHOLD)
            } else {
              NA_real_
            }
          )
      }
    }
  }

  ## Label inventory row.
  label_inventory <- data.table(
    task_id = task_id,
    scenario = scenario,
    repeat_id = repeat_id,
    all_imputation_label_rds = if (
      isTRUE(SAVE_ALL_IMPUTATION_LABELS_RDS)
    ) {
      label_rds_path
    } else {
      NA_character_
    },
    primary_imputation_label_csv_gz = if (
      isTRUE(SAVE_PRIMARY_IMPUTATION_LABELS_CSV_GZ)
    ) {
      primary_label_csv_path
    } else {
      NA_character_
    },
    all_imputation_rows = nrow(labels_all),
    primary_imputation_rows = sum(
      labels_all$imputation_id ==
        PRIMARY_IMPUTATION_ID
    ),
    label_columns = paste(
      names(labels_all),
      collapse = "|"
    )
  )

  list(
    ok = TRUE,
    config_signature = CONFIG_SIGNATURE,
    task_id = task_id,
    scenario = scenario,
    repeat_id = repeat_id,
    domain1_by_imputation = safe_rbindlist(domain1_rows),
    transport_by_imputation = safe_rbindlist(transport_rows),
    algorithm_agreement_by_imputation =
      safe_rbindlist(algorithm_agreement_rows),
    mice_stability = safe_rbindlist(stability_rows),
    cluster_composition = safe_rbindlist(composition_rows),
    model_metadata = safe_rbindlist(model_rows),
    centroids_long = safe_rbindlist(centroid_rows),
    preprocess_long = safe_rbindlist(preprocess_rows),
    label_inventory = label_inventory
  )
}

## =============================================================================
## 8. Parallel task execution
## =============================================================================

tasks <- CJ(
  scenario = SCENARIOS,
  repeat_id = seq_len(N_ACTIVE_REP)
)

tasks[, scenario_index := match(
  scenario,
  SCENARIOS
)]

tasks[, task_id := sprintf(
  "%s_rep_%03d",
  c(
    S1_pure_severity_continuum = "S1",
    S2_true_discrete_subtypes = "S2",
    S3_transport_drift = "S3"
  )[scenario],
  repeat_id
)]

tasks[, task_seed := GLOBAL_SEED +
  scenario_index * 100000L +
  repeat_id]

fwrite(
  tasks,
  file.path(TAB_DIR, "sensitivity07_task_manifest.csv")
)

configuration <- data.table(
  parameter = c(
    "ANALYSIS_ID",
    "RUN_MODE",
    "N_REP_PER_SCENARIO",
    "N_DERIVATION",
    "N_EXTERNAL",
    "MICE_M",
    "MICE_MAXIT",
    "PRIMARY_IMPUTATION_ID",
    "K_VALUES",
    "GMM_MODEL_NAMES",
    "WARD_SUBSAMPLE_N",
    "N_WORKERS",
    "GLOBAL_SEED",
    "CONFIG_SIGNATURE"
  ),
  value = c(
    ANALYSIS_ID,
    RUN_MODE,
    N_ACTIVE_REP,
    N_DERIVATION,
    N_EXTERNAL,
    MICE_M,
    MICE_MAXIT,
    PRIMARY_IMPUTATION_ID,
    paste(K_VALUES, collapse = "|"),
    paste(GMM_MODEL_NAMES, collapse = "|"),
    WARD_SUBSAMPLE_N,
    N_WORKERS,
    GLOBAL_SEED,
    CONFIG_SIGNATURE
  )
)

fwrite(
  configuration,
  file.path(TAB_DIR, "sensitivity07_configuration.csv")
)

## Branch OUTSIDE plan() and pass a bare strategy. Passing an if/else compound
## expression as plan()'s first argument (with the default substitute = TRUE)
## makes future try to parse it as a strategy and fails with
## "ls(envir = env): invalid 'envir' argument".
if (N_WORKERS > 1L) {
  future::plan(future::multisession, workers = N_WORKERS)
} else {
  future::plan(future::sequential)
}

log_msg("Total paired repeat tasks: ", nrow(tasks))
log_msg(
  "Each task contains both algorithms, K=2/K=3, ",
  "derivation and external labels, and all MICE datasets."
)

worker_fun <- function(i) {
  ## REQUIRED: explicit package loading in every multisession worker.
  suppressPackageStartupMessages({
    library(data.table)
    library(mice)
    library(mclust)
    library(cluster)
  })

  ## Prevent thread oversubscription inside multisession workers.
  data.table::setDTthreads(1L)

  task_row <- tasks[i]

  checkpoint_path <- file.path(
    CHECKPOINT_DIR,
    paste0(task_row$task_id, ".rds")
  )

  if (isTRUE(RESUME_COMPLETED_TASKS) &&
      file.exists(checkpoint_path)) {
    old <- tryCatch(
      readRDS(checkpoint_path),
      error = function(e) NULL
    )

    old_label_files_exist <- FALSE

    if (is.list(old) &&
        isTRUE(old$ok) &&
        identical(
          old$config_signature,
          CONFIG_SIGNATURE
        ) &&
        is.data.frame(old$label_inventory) &&
        nrow(old$label_inventory) == 1L) {

      required_label_files <- c(
        old$label_inventory$all_imputation_label_rds,
        old$label_inventory$primary_imputation_label_csv_gz
      )

      required_label_files <- required_label_files[
        !is.na(required_label_files) &
          nzchar(required_label_files)
      ]

      old_label_files_exist <- (
        length(required_label_files) > 0L &&
          all(file.exists(required_label_files))
      )
    }

    if (isTRUE(old_label_files_exist)) {
      return(old)
    }
  }

  started_at <- Sys.time()

  result <- tryCatch(
    {
      out <- run_task(task_row)
      out$started_at <- format(
        started_at,
        "%Y-%m-%d %H:%M:%S"
      )
      out$finished_at <- timestamp()
      out$elapsed_minutes <- as.numeric(
        difftime(
          Sys.time(),
          started_at,
          units = "mins"
        )
      )
      out
    },
    error = function(e) {
      list(
        ok = FALSE,
        config_signature = CONFIG_SIGNATURE,
        task_id = as.character(task_row$task_id),
        scenario = as.character(task_row$scenario),
        repeat_id = as.integer(task_row$repeat_id),
        started_at = format(
          started_at,
          "%Y-%m-%d %H:%M:%S"
        ),
        finished_at = timestamp(),
        elapsed_minutes = as.numeric(
          difftime(
            Sys.time(),
            started_at,
            units = "mins"
          )
        ),
        error = conditionMessage(e)
      )
    }
  )

  saveRDS(
    result,
    checkpoint_path,
    compress = "xz"
  )

  result
}

results <- future.apply::future_lapply(
  seq_len(nrow(tasks)),
  worker_fun,
  future.seed = TRUE,
  future.scheduling = 1,
  future.packages = c(
    "data.table",
    "mice",
    "mclust",
    "cluster"
  )
)

future::plan(future::sequential)
invisible(gc(full = TRUE))

## =============================================================================
## 9. Completion QC and failure handling
## =============================================================================

error_log <- rbindlist(
  lapply(results, function(z) {
    if (is.list(z) && !isTRUE(z$ok)) {
      data.table(
        task_id = z$task_id,
        scenario = z$scenario,
        repeat_id = z$repeat_id,
        started_at = z$started_at,
        finished_at = z$finished_at,
        elapsed_minutes = z$elapsed_minutes,
        error = z$error
      )
    } else {
      NULL
    }
  }),
  fill = TRUE
)

## When every task succeeds, the rbindlist above yields a 0-row, 0-column
## table, which makes fwrite warn and later `by = scenario` fail with
## "object 'scenario' not found". Give it the expected empty schema.
if (!nrow(error_log)) {
  error_log <- data.table(
    task_id = character(0),
    scenario = character(0),
    repeat_id = integer(0),
    started_at = character(0),
    finished_at = character(0),
    elapsed_minutes = numeric(0),
    error = character(0)
  )
}

fwrite(
  error_log,
  file.path(
    LOG_DIR,
    "sensitivity07_repeat_failure_log.csv"
  )
)

success <- Filter(
  function(z) {
    is.list(z) &&
      isTRUE(z$ok) &&
      identical(
        z$config_signature,
        CONFIG_SIGNATURE
      )
  },
  results
)

completion_qc <- tasks[, .(
  expected_repeats = .N
), by = scenario]

observed_success <- rbindlist(
  lapply(success, function(z) {
    data.table(
      scenario = z$scenario,
      repeat_id = z$repeat_id
    )
  }),
  fill = TRUE
)

## Same empty-schema guard for the all-failure case.
if (!nrow(observed_success)) {
  observed_success <- data.table(
    scenario = character(0),
    repeat_id = integer(0)
  )
}

success_counts <- observed_success[, .(
  successful_repeats = uniqueN(repeat_id)
), by = scenario]

failure_counts <- error_log[, .(
  failed_repeats = uniqueN(repeat_id)
), by = scenario]

completion_qc <- merge(
  completion_qc,
  success_counts,
  by = "scenario",
  all.x = TRUE
)

completion_qc <- merge(
  completion_qc,
  failure_counts,
  by = "scenario",
  all.x = TRUE
)

completion_qc[is.na(successful_repeats), successful_repeats := 0L]
completion_qc[is.na(failed_repeats), failed_repeats := 0L]

completion_qc[, completion_proportion :=
  successful_repeats / expected_repeats]

completion_qc[, missing_repeats :=
  expected_repeats -
    successful_repeats -
    failed_repeats]

fwrite(
  completion_qc,
  file.path(
    TAB_DIR,
    "sensitivity07_run_completion_qc.csv"
  )
)

if (!length(success)) {
  stop(
    "All combined sensitivity 07/07b repeats failed. ",
    "See: ",
    file.path(
      LOG_DIR,
      "sensitivity07_repeat_failure_log.csv"
    ),
    call. = FALSE
  )
}

log_msg(
  "Successful paired repeats: ",
  length(success),
  " / ",
  nrow(tasks),
  "; failed repeats: ",
  nrow(error_log)
)

## =============================================================================
## 10. Aggregate successful repeats
## =============================================================================

domain1_imp <- safe_rbindlist(
  lapply(success, `[[`, "domain1_by_imputation")
)

transport_imp <- safe_rbindlist(
  lapply(success, `[[`, "transport_by_imputation")
)

algorithm_agreement_imp <- safe_rbindlist(
  lapply(
    success,
    `[[`,
    "algorithm_agreement_by_imputation"
  )
)

mice_stability <- safe_rbindlist(
  lapply(success, `[[`, "mice_stability")
)

cluster_composition <- safe_rbindlist(
  lapply(success, `[[`, "cluster_composition")
)

model_metadata <- safe_rbindlist(
  lapply(success, `[[`, "model_metadata")
)

centroids_long <- safe_rbindlist(
  lapply(success, `[[`, "centroids_long")
)

preprocess_long <- safe_rbindlist(
  lapply(success, `[[`, "preprocess_long")
)

label_inventory <- safe_rbindlist(
  lapply(success, `[[`, "label_inventory")
)

domain1_repeat <- domain1_imp[, .(
  ari_vs_latent_subtype =
    safe_mean(ari_vs_latent_subtype),
  nmi_vs_latent_subtype =
    safe_mean(nmi_vs_latent_subtype),
  purity_vs_latent_subtype =
    safe_mean(purity_vs_latent_subtype),
  mortality_30d_C1_minus_C2 =
    safe_mean(mortality_30d_C1_minus_C2),
  cluster_1_prevalence =
    safe_mean(cluster_1_prevalence),
  cluster_min_prevalence =
    safe_mean(cluster_min_prevalence),
  ward_silhouette_mean =
    safe_mean(ward_silhouette_mean),
  gmm_bic =
    safe_mean(gmm_bic),
  gmm_entropy =
    safe_mean(gmm_entropy),
  gmm_mean_max_posterior =
    safe_mean(gmm_mean_max_posterior),
  gmm_min_max_posterior =
    safe_mean(gmm_min_max_posterior)
), by = .(
  task_id,
  scenario,
  repeat_id,
  algorithm,
  k
)]

transport_repeat <- transport_imp[, .(
  frozen_vs_reclustered_ari =
    safe_mean(frozen_vs_reclustered_ari),
  matched_agreement =
    safe_mean(matched_agreement),
  mean_abs_prevalence_difference =
    safe_mean(mean_abs_prevalence_difference),
  max_abs_prevalence_difference =
    safe_mean(max_abs_prevalence_difference),
  centroid_profile_correlation =
    safe_mean(centroid_profile_correlation),
  frozen_ari_vs_external_truth =
    safe_mean(frozen_ari_vs_external_truth),
  reclustered_ari_vs_external_truth =
    safe_mean(reclustered_ari_vs_external_truth),
  frozen_nmi_vs_external_truth =
    safe_mean(frozen_nmi_vs_external_truth),
  reclustered_nmi_vs_external_truth =
    safe_mean(reclustered_nmi_vs_external_truth)
), by = .(
  task_id,
  scenario,
  repeat_id,
  algorithm,
  k
)]

algorithm_agreement_repeat <-
  algorithm_agreement_imp[, .(
    gmm_vs_ward_ari =
      safe_mean(gmm_vs_ward_ari),
    gmm_vs_ward_matched_agreement =
      safe_mean(gmm_vs_ward_matched_agreement)
  ), by = .(
    task_id,
    scenario,
    repeat_id,
    cohort_label_source,
    k
  )]

## Direct K=3 minus K=2 transportability contrast.
transport_wide <- dcast(
  transport_repeat,
  task_id + scenario + repeat_id + algorithm ~ k,
  value.var = c(
    "frozen_vs_reclustered_ari",
    "matched_agreement",
    "mean_abs_prevalence_difference",
    "max_abs_prevalence_difference",
    "centroid_profile_correlation"
  )
)

setnames(
  transport_wide,
  old = grep(
    "_2$",
    names(transport_wide),
    value = TRUE
  ),
  new = sub(
    "_2$",
    "_K2",
    grep(
      "_2$",
      names(transport_wide),
      value = TRUE
    )
  )
)

setnames(
  transport_wide,
  old = grep(
    "_3$",
    names(transport_wide),
    value = TRUE
  ),
  new = sub(
    "_3$",
    "_K3",
    grep(
      "_3$",
      names(transport_wide),
      value = TRUE
    )
  )
)

transport_k_comparison <- transport_wide[, .(
  task_id,
  scenario,
  repeat_id,
  algorithm,
  ari_K2 =
    frozen_vs_reclustered_ari_K2,
  ari_K3 =
    frozen_vs_reclustered_ari_K3,
  ari_K3_minus_K2 =
    frozen_vs_reclustered_ari_K3 -
    frozen_vs_reclustered_ari_K2,
  matched_agreement_K2 =
    matched_agreement_K2,
  matched_agreement_K3 =
    matched_agreement_K3,
  matched_agreement_K3_minus_K2 =
    matched_agreement_K3 -
    matched_agreement_K2,
  mean_abs_prevalence_difference_K2 =
    mean_abs_prevalence_difference_K2,
  mean_abs_prevalence_difference_K3 =
    mean_abs_prevalence_difference_K3,
  prevalence_error_K3_minus_K2 =
    mean_abs_prevalence_difference_K3 -
    mean_abs_prevalence_difference_K2,
  profile_correlation_K2 =
    centroid_profile_correlation_K2,
  profile_correlation_K3 =
    centroid_profile_correlation_K3,
  profile_correlation_K3_minus_K2 =
    centroid_profile_correlation_K3 -
    centroid_profile_correlation_K2
)]

## Same-repeat algorithm contrast: GMM minus Ward.
d1_algorithm_wide <- dcast(
  domain1_repeat,
  task_id + scenario + repeat_id + k ~ algorithm,
  value.var = c(
    "ari_vs_latent_subtype",
    "nmi_vs_latent_subtype",
    "purity_vs_latent_subtype",
    "mortality_30d_C1_minus_C2",
    "cluster_min_prevalence"
  )
)

algorithm_domain1_contrast <-
  d1_algorithm_wide[, .(
    task_id,
    scenario,
    repeat_id,
    k,
    ari_GMM_minus_Ward =
      ari_vs_latent_subtype_gmm -
      ari_vs_latent_subtype_ward,
    nmi_GMM_minus_Ward =
      nmi_vs_latent_subtype_gmm -
      nmi_vs_latent_subtype_ward,
    purity_GMM_minus_Ward =
      purity_vs_latent_subtype_gmm -
      purity_vs_latent_subtype_ward,
    outcome_separation_GMM_minus_Ward =
      mortality_30d_C1_minus_C2_gmm -
      mortality_30d_C1_minus_C2_ward,
    minimum_cluster_prevalence_GMM_minus_Ward =
      cluster_min_prevalence_gmm -
      cluster_min_prevalence_ward
  )]

transport_algorithm_wide <- dcast(
  transport_repeat,
  task_id + scenario + repeat_id + k ~ algorithm,
  value.var = c(
    "frozen_vs_reclustered_ari",
    "matched_agreement",
    "mean_abs_prevalence_difference",
    "centroid_profile_correlation"
  )
)

algorithm_transport_contrast <-
  transport_algorithm_wide[, .(
    task_id,
    scenario,
    repeat_id,
    k,
    transport_ARI_GMM_minus_Ward =
      frozen_vs_reclustered_ari_gmm -
      frozen_vs_reclustered_ari_ward,
    matched_agreement_GMM_minus_Ward =
      matched_agreement_gmm -
      matched_agreement_ward,
    prevalence_error_GMM_minus_Ward =
      mean_abs_prevalence_difference_gmm -
      mean_abs_prevalence_difference_ward,
    profile_correlation_GMM_minus_Ward =
      centroid_profile_correlation_gmm -
      centroid_profile_correlation_ward
  )]

## =============================================================================
## 11. Operating characteristics and Monte Carlo error
## =============================================================================

operating_characteristics <- safe_rbindlist(list(
  summarise_binary(
    domain1_repeat[
      scenario == "S2_true_discrete_subtypes" &
        k == 3L
    ][, success_flag :=
        ari_vs_latent_subtype >=
          TRUE_K_ARI_THRESHOLD],
    c("scenario", "algorithm", "k"),
    "success_flag",
    "S2_true_K3_recovery_ARI_ge_threshold"
  ),

  summarise_binary(
    domain1_repeat[
      scenario == "S1_pure_severity_continuum" &
        k == 3L
    ][, false_recovery_flag :=
        ari_vs_latent_subtype >=
          TRUE_K_ARI_THRESHOLD],
    c("scenario", "algorithm", "k"),
    "false_recovery_flag",
    "S1_false_true_K3_recovery_ARI_ge_threshold"
  ),

  summarise_binary(
    transport_repeat[
      ,
      transport_success :=
        frozen_vs_reclustered_ari >=
          ASSIGNMENT_ARI_THRESHOLD
    ],
    c("scenario", "algorithm", "k"),
    "transport_success",
    "assignment_transport_ARI_ge_threshold"
  ),

  summarise_binary(
    transport_repeat[
      ,
      material_prevalence_drift :=
        max_abs_prevalence_difference >
          PREVALENCE_DRIFT_THRESHOLD
    ],
    c("scenario", "algorithm", "k"),
    "material_prevalence_drift",
    "material_prevalence_drift"
  ),

  summarise_binary(
    algorithm_agreement_repeat[
      ,
      algorithm_agreement_success :=
        gmm_vs_ward_ari >=
          ALGORITHM_AGREEMENT_ARI_THRESHOLD
    ],
    c(
      "scenario",
      "cohort_label_source",
      "k"
    ),
    "algorithm_agreement_success",
    "GMM_vs_Ward_ARI_ge_threshold"
  ),

  summarise_binary(
    mice_stability[
      ,
      mice_stable :=
        ari_mean_vs_primary >=
          MICE_STABILITY_ARI_THRESHOLD
    ],
    c(
      "scenario",
      "algorithm",
      "label_source",
      "k"
    ),
    "mice_stable",
    "MICE_stability_ARI_ge_threshold"
  )
))

mcse_summary <- safe_rbindlist(list(
  summarise_numeric(
    domain1_repeat,
    c("scenario", "algorithm", "k"),
    "ari_vs_latent_subtype",
    "ARI_vs_latent_subtype"
  ),
  summarise_numeric(
    domain1_repeat,
    c("scenario", "algorithm", "k"),
    "nmi_vs_latent_subtype",
    "NMI_vs_latent_subtype"
  ),
  summarise_numeric(
    domain1_repeat,
    c("scenario", "algorithm", "k"),
    "purity_vs_latent_subtype",
    "Purity_vs_latent_subtype"
  ),
  summarise_numeric(
    domain1_repeat[k == 2L],
    c("scenario", "algorithm"),
    "mortality_30d_C1_minus_C2",
    "K2_mortality_separation"
  ),
  summarise_numeric(
    transport_repeat,
    c("scenario", "algorithm", "k"),
    "frozen_vs_reclustered_ari",
    "Frozen_vs_reclustered_ARI"
  ),
  summarise_numeric(
    transport_repeat,
    c("scenario", "algorithm", "k"),
    "max_abs_prevalence_difference",
    "Maximum_absolute_prevalence_difference"
  ),
  summarise_numeric(
    transport_repeat,
    c("scenario", "algorithm", "k"),
    "centroid_profile_correlation",
    "Centroid_profile_correlation"
  ),
  summarise_numeric(
    algorithm_agreement_repeat,
    c(
      "scenario",
      "cohort_label_source",
      "k"
    ),
    "gmm_vs_ward_ari",
    "GMM_vs_Ward_ARI"
  ),
  summarise_numeric(
    mice_stability,
    c(
      "scenario",
      "algorithm",
      "label_source",
      "k"
    ),
    "ari_mean_vs_primary",
    "MICE_label_stability_ARI"
  ),
  summarise_numeric(
    algorithm_domain1_contrast,
    c("scenario", "k"),
    "ari_GMM_minus_Ward",
    "Paired_ARI_GMM_minus_Ward"
  ),
  summarise_numeric(
    algorithm_transport_contrast,
    c("scenario", "k"),
    "transport_ARI_GMM_minus_Ward",
    "Paired_transport_ARI_GMM_minus_Ward"
  ),
  summarise_numeric(
    transport_k_comparison,
    c("scenario", "algorithm"),
    "ari_K3_minus_K2",
    "Transport_ARI_K3_minus_K2"
  )
))

## =============================================================================
## 12. ADEMP and label documentation
## =============================================================================

ademp <- data.table(
  component = c(
    "Aim",
    "Data-generating mechanisms",
    "Estimands",
    "Methods",
    "Performance measures",
    "Repetitions and Monte Carlo error",
    "Failure policy",
    "Individual-label preservation",
    "Primary interpretation rule"
  ),
  specification = c(
    paste(
      "Combine the original algorithm-behavior and transportability",
      "sensitivities by comparing GMM and Ward for K=2/K=3 on",
      "the same simulated patients."
    ),
    paste(
      "S1 pure severity continuum; S2 strong true K=3 structure;",
      "S3 true K=3 structure with external subtype-prevalence,",
      "severity, feature-distribution, and missingness drift."
    ),
    paste(
      "Latent-truth ARI/NMI/purity; K2 outcome separation;",
      "GMM-Ward paired agreement; MICE stability;",
      "frozen-vs-reclustered ARI and matched agreement;",
      "prevalence drift; centroid-profile correspondence;",
      "K3-minus-K2 and GMM-minus-Ward paired contrasts."
    ),
    paste(
      "MICE m=5; 1st/99th percentile winsorization;",
      "derivation-frozen preprocessing for external assignment;",
      "external-local preprocessing for independent reclustering;",
      "strict mclust::Mclust GMM and Ward.D2."
    ),
    paste(
      "Across-repeat mean, SD, empirical 2.5%-97.5% interval,",
      "MCSE, threshold-based operating characteristics,",
      "task-failure counts, and exact binomial intervals."
    ),
    paste0(
      N_ACTIVE_REP,
      " paired repeats per scenario in this run; algorithms and K",
      " are paired within repeat and imputation."
    ),
    paste(
      "Any derivation/external GMM failure or invalid GMM return",
      "terminates that repeat. No k-means fallback is permitted."
    ),
    paste(
      "Every successful task writes an all-imputation compressed RDS",
      "and a primary-imputation compressed CSV containing same-patient",
      "GMM/Ward K2/K3 labels, frozen labels, reclustered labels,",
      "cross-imputation alignment, and Ward-to-GMM matched labels."
    ),
    paste(
      "Stability, prognostic separation, or algorithm agreement alone",
      "does not establish discreteness. S1/S2/S3 must be interpreted",
      "jointly using truth recovery, imputation stability, paired",
      "algorithm agreement, and assignment transportability.",
      "Only S3 applies external drift; in S1/S2 the external cohort is an",
      "independent draw from the same data-generating mechanism, so",
      "frozen-vs-reclustered agreement there reflects same-distribution",
      "assignment reproducibility rather than transportability under drift."
    )
  )
)

label_dictionary <- data.table(
  field_pattern = c(
    "*_raw",
    "*_aligned_to_primary_imputation",
    "*_frozen_raw",
    "*_frozen_aligned_to_primary_imputation",
    "*_reclustered_raw",
    "*_reclustered_matched_to_frozen",
    "*_reclustered_aligned_to_primary_imputation",
    "ward_*_matched_to_gmm_same_dataset",
    "latent_subtype"
  ),
  definition = c(
    "Unmodified algorithm output label.",
    paste(
      "Label permuted to match the same algorithm/K in completed",
      "dataset 1 on the same patients."
    ),
    paste(
      "External label assigned by the derivation-fitted model and",
      "derivation preprocessing."
    ),
    paste(
      "Frozen external label transformed using the derivation",
      "cross-imputation permutation."
    ),
    "External label from independent external-local reclustering.",
    paste(
      "External reclustered label optimally permuted to the frozen",
      "external labels in the same imputation."
    ),
    paste(
      "External reclustered label permuted to completed dataset 1",
      "external reclustering for the same algorithm/K."
    ),
    paste(
      "Ward label optimally permuted to GMM on the exact same",
      "patients, cohort/source, imputation, and K."
    ),
    "Known latent class used only for simulation evaluation."
  )
)

## =============================================================================
## 13. Save tables
## =============================================================================

tables_to_write <- list(
  sensitivity07_domain1_by_imputation =
    domain1_imp,
  sensitivity07_domain1_by_repeat =
    domain1_repeat,
  sensitivity07_assignment_transportability_by_imputation =
    transport_imp,
  sensitivity07_assignment_transportability_by_repeat =
    transport_repeat,
  sensitivity07_k2_k3_transportability_comparison =
    transport_k_comparison,
  sensitivity07_algorithm_agreement_by_imputation =
    algorithm_agreement_imp,
  sensitivity07_algorithm_agreement_by_repeat =
    algorithm_agreement_repeat,
  sensitivity07_domain1_GMM_minus_Ward_paired_contrast =
    algorithm_domain1_contrast,
  sensitivity07_transport_GMM_minus_Ward_paired_contrast =
    algorithm_transport_contrast,
  sensitivity07_mice_label_stability =
    mice_stability,
  sensitivity07_cluster_composition =
    cluster_composition,
  sensitivity07_model_metadata =
    model_metadata,
  sensitivity07_centroids_long =
    centroids_long,
  sensitivity07_preprocess_parameters_long =
    preprocess_long,
  sensitivity07_individual_label_file_inventory =
    label_inventory,
  sensitivity07_individual_label_dictionary =
    label_dictionary,
  sensitivity07_operating_characteristics =
    operating_characteristics,
  sensitivity07_mcse_summary =
    mcse_summary,
  sensitivity07_ADEMP_table =
    ademp,
  sensitivity07_repeat_failure_log =
    error_log
)

for (nm in names(tables_to_write)) {
  fwrite(
    tables_to_write[[nm]],
    file.path(
      TAB_DIR,
      paste0(nm, ".csv")
    )
  )
}

## =============================================================================
## 14. Figures
## =============================================================================

if (nrow(domain1_repeat)) {
  d1_plot_dt <- domain1_repeat[, .(
    mean = safe_mean(ari_vs_latent_subtype),
    low = safe_quantile(
      ari_vs_latent_subtype,
      0.025
    ),
    high = safe_quantile(
      ari_vs_latent_subtype,
      0.975
    )
  ), by = .(
    scenario,
    algorithm,
    k
  )]

  p_d1 <- ggplot(
    d1_plot_dt,
    aes(
      x = factor(k),
      y = mean,
      shape = algorithm,
      group = algorithm
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.25),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = low,
        ymax = high
      ),
      width = 0.12,
      position = position_dodge(width = 0.25)
    ) +
    facet_wrap(~ scenario) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "K",
      y = "Adjusted Rand index versus latent truth",
      shape = "Algorithm",
      title = "Latent-structure recovery by GMM and Ward"
    ) +
    theme_publication()

  save_plot_both(
    p_d1,
    "Figure_S07_latent_truth_recovery",
    width = 10,
    height = 5
  )
}

if (nrow(transport_repeat)) {
  tr_plot_dt <- transport_repeat[, .(
    mean = safe_mean(
      frozen_vs_reclustered_ari
    ),
    low = safe_quantile(
      frozen_vs_reclustered_ari,
      0.025
    ),
    high = safe_quantile(
      frozen_vs_reclustered_ari,
      0.975
    )
  ), by = .(
    scenario,
    algorithm,
    k
  )]

  p_tr <- ggplot(
    tr_plot_dt,
    aes(
      x = factor(k),
      y = mean,
      shape = algorithm,
      group = algorithm
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.25),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = low,
        ymax = high
      ),
      width = 0.12,
      position = position_dodge(width = 0.25)
    ) +
    facet_wrap(~ scenario) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "K",
      y = "Frozen-versus-reclustered ARI",
      shape = "Algorithm",
      title = "External assignment transportability"
    ) +
    theme_publication()

  save_plot_both(
    p_tr,
    "Figure_S07_assignment_transportability",
    width = 10,
    height = 5
  )
}

if (nrow(algorithm_agreement_repeat)) {
  ag_plot_dt <- algorithm_agreement_repeat[, .(
    mean = safe_mean(gmm_vs_ward_ari),
    low = safe_quantile(
      gmm_vs_ward_ari,
      0.025
    ),
    high = safe_quantile(
      gmm_vs_ward_ari,
      0.975
    )
  ), by = .(
    scenario,
    cohort_label_source,
    k
  )]

  p_ag <- ggplot(
    ag_plot_dt,
    aes(
      x = factor(k),
      y = mean
    )
  ) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = low,
        ymax = high
      ),
      width = 0.12
    ) +
    facet_grid(
      cohort_label_source ~ scenario
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "K",
      y = "GMM-versus-Ward ARI",
      title = "Same-patient algorithm agreement"
    ) +
    theme_publication()

  save_plot_both(
    p_ag,
    "Figure_S07_same_patient_algorithm_agreement",
    width = 11,
    height = 7
  )
}

if (nrow(mice_stability)) {
  mice_plot_dt <- mice_stability[, .(
    mean = safe_mean(
      ari_mean_vs_primary
    ),
    low = safe_quantile(
      ari_mean_vs_primary,
      0.025
    ),
    high = safe_quantile(
      ari_mean_vs_primary,
      0.975
    )
  ), by = .(
    scenario,
    algorithm,
    label_source,
    k
  )]

  p_mice <- ggplot(
    mice_plot_dt,
    aes(
      x = factor(k),
      y = mean,
      shape = algorithm,
      group = algorithm
    )
  ) +
    geom_point(
      position = position_dodge(width = 0.25),
      size = 2
    ) +
    geom_errorbar(
      aes(
        ymin = low,
        ymax = high
      ),
      width = 0.12,
      position = position_dodge(width = 0.25)
    ) +
    facet_grid(
      label_source ~ scenario
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "K",
      y = "ARI versus completed dataset 1",
      shape = "Algorithm",
      title = "Cross-imputation label stability"
    ) +
    theme_publication()

  save_plot_both(
    p_mice,
    "Figure_S07_mice_label_stability",
    width = 11,
    height = 7
  )
}

if (nrow(transport_k_comparison)) {
  k_plot_dt <- transport_k_comparison[, .(
    mean = safe_mean(ari_K3_minus_K2),
    low = safe_quantile(
      ari_K3_minus_K2,
      0.025
    ),
    high = safe_quantile(
      ari_K3_minus_K2,
      0.975
    )
  ), by = .(
    scenario,
    algorithm
  )]

  p_k <- ggplot(
    k_plot_dt,
    aes(
      x = algorithm,
      y = mean
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = low,
        ymax = high
      ),
      width = 0.12
    ) +
    facet_wrap(~ scenario) +
    labs(
      x = NULL,
      y = "Assignment ARI: K3 minus K2",
      title = "K=3 versus K=2 external transportability"
    ) +
    theme_publication()

  save_plot_both(
    p_k,
    "Figure_S07_K3_minus_K2_transportability",
    width = 9,
    height = 5
  )
}

## =============================================================================
## 15. Completion markers
## =============================================================================

all_tasks_successful <- (
  length(success) == nrow(tasks) &&
    nrow(error_log) == 0L
)

completion_summary <- data.table(
  analysis_id = ANALYSIS_ID,
  run_mode = RUN_MODE,
  expected_tasks = nrow(tasks),
  successful_tasks = length(success),
  failed_tasks = nrow(error_log),
  all_tasks_successful = all_tasks_successful,
  completed_at = timestamp(),
  config_signature = CONFIG_SIGNATURE
)

fwrite(
  completion_summary,
  file.path(
    LOG_DIR,
    "sensitivity07_run_completion_summary.csv"
  )
)

completed_marker <- file.path(
  LOG_DIR,
  "run_completed.ok"
)

partial_marker <- file.path(
  LOG_DIR,
  "run_partial_with_failures.txt"
)

if (file.exists(completed_marker)) {
  file.remove(completed_marker)
}

if (file.exists(partial_marker)) {
  file.remove(partial_marker)
}

if (isTRUE(all_tasks_successful)) {
  writeLines(
    paste0(
      "completed_at=", timestamp(), "\n",
      "run_mode=", RUN_MODE, "\n",
      "expected_tasks=", nrow(tasks), "\n",
      "successful_tasks=", length(success), "\n",
      "failed_tasks=0\n",
      "config_signature=", CONFIG_SIGNATURE, "\n"
    ),
    completed_marker
  )
} else if (isTRUE(ALLOW_PARTIAL_OUTPUTS)) {
  writeLines(
    paste0(
      "partial_at=", timestamp(), "\n",
      "run_mode=", RUN_MODE, "\n",
      "expected_tasks=", nrow(tasks), "\n",
      "successful_tasks=", length(success), "\n",
      "failed_tasks=", nrow(error_log), "\n",
      "Rerun the same script to retry failed repeats.\n",
      "config_signature=", CONFIG_SIGNATURE, "\n"
    ),
    partial_marker
  )
}

write_session_info(
  file.path(
    LOG_DIR,
    "sessionInfo.txt"
  )
)

capture.output(
  warnings(),
  file = file.path(
    LOG_DIR,
    "warnings_after_run.txt"
  )
)

log_msg("Combined sensitivity 07/07b finished.")
log_msg("All tasks successful: ", all_tasks_successful)
log_msg("Tables: ", TAB_DIR)
log_msg("Figures: ", FIG_DIR)
log_msg("Individual labels RDS: ", LABEL_RDS_DIR)
log_msg("Primary-imputation labels CSV.GZ: ", LABEL_PRIMARY_DIR)

cat("\n============================================================\n")
cat("Combined sensitivity 07/07b finished.\n")
cat("Run mode:", RUN_MODE, "\n")
cat("Expected tasks:", nrow(tasks), "\n")
cat("Successful tasks:", length(success), "\n")
cat("Failed tasks:", nrow(error_log), "\n")
cat("All tasks successful:", all_tasks_successful, "\n")
cat("Output folder:\n", OUT_DIR, "\n")
cat("============================================================\n")
