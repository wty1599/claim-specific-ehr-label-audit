#!/usr/bin/env Rscript

# WP6: Distribution fidelity of the five highest-loading variables.
#
# This is a formal post hoc descriptive audit. It compares the locked empirical MIMIC
# primary-imputation matrix with an exactly regenerated formal S1 derivation
# matrix (repeat 1, primary imputation). It does not change the D1 verdict and
# does not treat point masses as evidence of a biological subtype.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(ggplot2)
})

required_pkgs <- c("data.table", "mice", "ggplot2", "digest")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       call. = FALSE)
}

RUN_MODE <- "formal"
RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid())
GLOBAL_AUDIT_SEED <- 20260718L

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
AUDIT_ROOT <- file.path(
  PROJECT_ROOT,
  "domain1_renal_identity_and_method_closure_20260717"
)
SCRIPT_PATH <- file.path(
  AUDIT_ROOT, "scripts", "29e_domain1_top5_distribution_fidelity.R"
)
TABLE_DIR <- file.path(AUDIT_ROOT, "tables")
FIGURE_DIR <- file.path(AUDIT_ROOT, "figures")
LOG_DIR <- file.path(AUDIT_ROOT, "logs")
PROVENANCE_DIR <- file.path(AUDIT_ROOT, "provenance")
CHECKPOINT_DIR <- file.path(AUDIT_ROOT, "checkpoints")

for (d in c(TABLE_DIR, FIGURE_DIR, LOG_DIR, PROVENANCE_DIR, CHECKPOINT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

INPUT <- list(
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  mice_primary = file.path(PROJECT_ROOT, "output", "model", "mice_primary.rds"),
  locked_mimic_matrix = file.path(
    PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"
  ),
  feature_sets = file.path(PROJECT_ROOT, "output", "model", "feature_sets.rds"),
  locked_simulation_script = file.path(
    PROJECT_ROOT, "analysis_archive", "simulations", "locked_results_20260712",
    "05_method_source_archive_20260712", "A_locked_result_sources",
    "simulation_four_domain_framework_mice_parallel_v9_extended_validation_NO_PAC_FORMAL.R"
  ),
  locked_simulation_results = file.path(
    PROJECT_ROOT, "analysis_archive", "simulations", "main_extended_validation_no_pac_formal",
    "output", "raw_results", "simulation_four_domain_all_results_mice.rds"
  )
)

OUTPUT <- list(
  pointmass = file.path(TABLE_DIR, "Table_D1_real_vs_S1_pointmass.csv"),
  distance = file.path(TABLE_DIR, "Table_D1_real_vs_S1_distribution_distance.csv"),
  figure_source = file.path(
    TABLE_DIR, "29e_Figure_D1_top5_distribution_fidelity_source_data.csv"
  ),
  figure_pdf = file.path(FIGURE_DIR, "Figure_D1_top5_distribution_fidelity.pdf"),
  figure_png = file.path(FIGURE_DIR, "Figure_D1_top5_distribution_fidelity.png"),
  figure_tiff = file.path(FIGURE_DIR, "Figure_D1_top5_distribution_fidelity.tiff"),
  log = file.path(LOG_DIR, "29e_distribution_fidelity_log.md"),
  qc = file.path(LOG_DIR, "29e_distribution_fidelity_QC.csv"),
  session = file.path(LOG_DIR, "29e_sessionInfo.txt"),
  input_checksums = file.path(PROVENANCE_DIR, "29e_WP6_input_checksums.csv"),
  run_metadata = file.path(PROVENANCE_DIR, "29e_WP6_run_metadata.csv"),
  output_manifest = file.path(PROVENANCE_DIR, "29e_WP6_output_sha256_manifest.csv"),
  completion = file.path(AUDIT_ROOT, "29e_run_completed.ok"),
  in_progress = file.path(CHECKPOINT_DIR, "29e_run_in_progress.txt")
)

final_targets <- unlist(OUTPUT[setdiff(names(OUTPUT), "in_progress")], use.names = FALSE)
already_present <- final_targets[file.exists(final_targets)]
if (length(already_present)) {
  stop(
    "WP6 outputs already exist; refusing to overwrite:\n",
    paste(already_present, collapse = "\n"), call. = FALSE
  )
}

missing_inputs <- unlist(INPUT, use.names = FALSE)
missing_inputs <- missing_inputs[!file.exists(missing_inputs)]
if (length(missing_inputs)) {
  stop("Missing authoritative inputs:\n", paste(missing_inputs, collapse = "\n"),
       call. = FALSE)
}

writeLines(c(
  "WP6 RUN IN PROGRESS",
  paste0("run_id=", RUN_ID),
  paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script=", SCRIPT_PATH)
), OUTPUT$in_progress, useBytes = TRUE)

sha256_file <- function(path) {
  unname(tools::md5sum(path))
}

sha256_windows <- function(path) {
  digest::digest(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
}

script_hash_start <- sha256_windows(SCRIPT_PATH)

input_hash_start <- rbindlist(lapply(names(INPUT), function(nm) {
  data.table(
    input_name = nm,
    path = normalizePath(INPUT[[nm]], winslash = "/", mustWork = TRUE),
    sha256 = sha256_windows(INPUT[[nm]]),
    bytes = file.info(INPUT[[nm]])$size,
    modified = format(file.info(INPUT[[nm]])$mtime, "%Y-%m-%d %H:%M:%S %Z")
  )
}))

TOP5 <- c(
  "bun_max", "creatinine_max", "aniongap_max",
  "bicarbonate_min", "lactate_max"
)

VARIABLE_LABELS <- c(
  bun_max = "Maximum BUN",
  creatinine_max = "Maximum creatinine",
  aniongap_max = "Maximum anion gap",
  bicarbonate_min = "Minimum bicarbonate",
  lactate_max = "Maximum lactate"
)

DATASET_LABELS <- c(
  MIMIC = "MIMIC-IV empirical",
  S1 = "S1 continuous-severity"
)

COLORS <- c(
  "MIMIC-IV empirical" = "#2E2E2E",
  "S1 continuous-severity" = "#738A96"
)

LINE_TYPES <- c(
  "MIMIC-IV empirical" = "solid",
  "S1 continuous-severity" = "22"
)

message("WP6 palette: ", paste(names(COLORS), COLORS, sep = "=", collapse = "; "))
message("Figure size: 183 x 210 mm; PNG/TIFF 600 dpi")

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

safe_prop_at <- function(x, value, tol = 1e-12) {
  x <- x[is.finite(x)]
  if (!length(x) || !is.finite(value)) return(c(n = NA_real_, prop = NA_real_))
  n <- sum(abs(x - value) <= tol * pmax(1, abs(value)))
  c(n = n, prop = n / length(x))
}

distribution_moments <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) {
    return(c(skewness = NA_real_, excess_kurtosis = NA_real_))
  }
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(c(skewness = NA_real_, excess_kurtosis = NA_real_))
  }
  z <- (x - mean(x)) / s
  c(skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3)
}

common_quantile_wasserstein <- function(x, y, grid_n = 10001L) {
  probs <- seq(0, 1, length.out = grid_n)
  qx <- stats::quantile(x, probs = probs, na.rm = TRUE, type = 7, names = FALSE)
  qy <- stats::quantile(y, probs = probs, na.rm = TRUE, type = 7, names = FALSE)
  mean(abs(qx - qy))
}

extract_assignment_name <- function(expr) {
  if (!is.call(expr)) return(NA_character_)
  op <- as.character(expr[[1L]])[1L]
  if (!op %in% c("<-", "=")) return(NA_character_)
  lhs <- expr[[2L]]
  if (!is.symbol(lhs)) return(NA_character_)
  as.character(lhs)
}

load_locked_simulation_engine <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  engine <- new.env(parent = globalenv())
  keep <- c(
    "GLOBAL_SEED", "N_TRAIN", "N_EXTERNAL", "P_FEATURES", "N_REP",
    "USE_TEST_MODE", "MICE_M", "MICE_MAXIT", "PRIMARY_IMPUTATION_ID",
    "S2_SUBTYPE_DELTA", "S2_SEVERITY_SIGNAL_MULTIPLIER",
    "S2_ORGAN_SIGNAL_MULTIPLIER", "S2_NOISE_SD",
    "S5_LATENT_PREVALENCE", "S5_FEATURE_DELTA",
    "S5_LATENT_OUTCOME_LOGOR_MORT", "S5_LATENT_OUTCOME_LOGOR_MAKE",
    "S5_SEVERITY_SIGNAL_MULTIPLIER", "S5_ORGAN_SIGNAL_MULTIPLIER",
    "S5_NOISE_SD", "S5_ORACLE_LABEL_ACCURACY",
    "S4_STRONG_GAMMA0", "S4_STRONG_GAMMA_SEVERITY",
    "S4_STRONG_GAMMA_RENAL", "S4_STRONG_GAMMA_SHOCK",
    "WINSOR_PROBS", "SCENARIOS", "fix_duplicate_colnames", "safe_rbindlist", "feature_names",
    "baseline_vars", "severity_loading", "feature_scale", "positive_features",
    "bounded_low", "bounded_high", "safe_plogis", "clamp_prob", "logit",
    "winsor_fit", "winsor_apply", "standardize_train", "standardize_apply",
    "assert_no_missing_features", "make_feature_covariance",
    "to_clinical_scale", "inject_missingness", "simulate_latent_features",
    "generate_baseline_covariates", "generate_treatment", "generate_outcomes",
    "generate_dataset", "make_mice_imputations"
  )
  found <- character()
  for (expr in exprs) {
    nm <- extract_assignment_name(expr)
    if (!is.na(nm) && nm %in% keep) {
      eval(expr, envir = engine)
      found <- c(found, nm)
    }
  }
  absent <- setdiff(keep, found)
  if (length(absent)) {
    stop("Locked simulation source is missing required definitions: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  engine
}

qc_rows <- list()
add_qc <- function(check, pass, observed, expected, detail = "") {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table(
    check = check,
    pass = isTRUE(pass),
    observed = as.character(observed),
    expected = as.character(expected),
    detail = as.character(detail)
  )
}

# -------------------------------------------------------------------------
# MIMIC: recover original values, locked winsor-before-MICE bounds, and the
# exact primary-imputation standardized matrix.
# -------------------------------------------------------------------------
feature_sets <- readRDS(INPUT$feature_sets)
assert("primary" %in% names(feature_sets), "feature_sets.rds has no primary set")
assert(identical(length(feature_sets$primary), 33L), "Primary feature count is not 33")
assert(all(TOP5 %in% feature_sets$primary), "Top-five variables are not all primary features")
primary_vars <- as.character(feature_sets$primary)

mimic_original <- fread(INPUT$final_full, select = c("stay_id", TOP5))
mimic_original[, (TOP5) := lapply(.SD, as.numeric), .SDcols = TOP5]

mids <- readRDS(INPUT$mice_primary)
assert(inherits(mids, "mids"), "mice_primary.rds is not a mids object")
assert(mids$m >= 5L, "MIMIC MICE object contains fewer than five imputations")
assert(mids$iteration == 10L, "MIMIC MICE object maxit is not 10")
assert(nrow(mids$data) == 20049L, "MIMIC MICE denominator is not 20,049")
assert(all(TOP5 %in% names(mids$data)), "Top-five variables missing from MIMIC MICE object")

mimic_completed <- as.data.table(mice::complete(mids, action = 1L))
mimic_locked <- as.data.table(readRDS(INPUT$locked_mimic_matrix))
assert(nrow(mimic_locked) == 20049L, "Locked MIMIC matrix denominator is not 20,049")
assert(all(c("stay_id", TOP5) %in% names(mimic_locked)),
       "Locked MIMIC matrix is missing required columns")
assert(!anyDuplicated(mimic_original$stay_id), "final_full contains duplicate stay_id values")
assert(!anyDuplicated(mimic_locked$stay_id), "Locked MIMIC matrix contains duplicate stay_id values")
mimic_match <- match(as.character(mimic_locked$stay_id),
                     as.character(mimic_original$stay_id))
assert(!anyNA(mimic_match), "Some locked MIMIC stay_id values are absent from final_full")
mimic_original <- mimic_original[mimic_match]
assert(identical(as.character(mimic_original$stay_id), as.character(mimic_locked$stay_id)),
       "MIMIC final_full and locked matrix row order differ")

mimic_bounds <- rbindlist(lapply(TOP5, function(v) {
  q <- stats::quantile(mimic_original[[v]], probs = c(0.01, 0.99),
                       na.rm = TRUE, type = 7, names = FALSE)
  data.table(variable = v, winsor_low = q[1], winsor_high = q[2])
}))

mimic_alignment <- rbindlist(lapply(TOP5, function(v) {
  b <- mimic_bounds[variable == v]
  raw_w <- pmin(pmax(mimic_original[[v]], b$winsor_low), b$winsor_high)
  obs <- which(!is.na(mids$data[[v]]))
  data.table(
    variable = v,
    observed_n = length(obs),
    max_abs_diff = if (length(obs)) {
      max(abs(raw_w[obs] - mids$data[[v]][obs]), na.rm = TRUE)
    } else NA_real_
  )
}))

mimic_z_reconstructed <- scale(as.matrix(mimic_completed[, ..primary_vars]))
mimic_z_locked <- as.matrix(mimic_locked[, ..primary_vars])
mimic_matrix_max_abs <- max(abs(mimic_z_reconstructed - mimic_z_locked))

add_qc(
  "mimic_denominator", nrow(mimic_locked) == 20049L,
  nrow(mimic_locked), 20049L,
  "Locked empirical matrix"
)
add_qc(
  "mimic_mice_settings", mids$m >= 5L && mids$iteration == 10L,
  paste0("m=", mids$m, "; maxit=", mids$iteration, "; primary=1"),
  "m>=5; maxit=10; primary=1",
  "Primary imputation is fixed at completed dataset 1"
)
add_qc(
  "mimic_pre_mice_winsor_alignment",
  all(is.finite(mimic_alignment$max_abs_diff)) &&
    max(mimic_alignment$max_abs_diff) <= 1e-10,
  format(max(mimic_alignment$max_abs_diff), scientific = TRUE), "<=1e-10",
  "Observed values in mids$data reproduce 1st/99th-percentile winsorized final_full"
)
add_qc(
  "mimic_locked_matrix_reproduction", mimic_matrix_max_abs <= 1e-12,
  format(mimic_matrix_max_abs, scientific = TRUE), "<=1e-12",
  "scale(complete(mice_primary, 1)) versus locked X_primary33_std_mice"
)

# -------------------------------------------------------------------------
# S1: load definitions from the locked formal source, regenerate formal
# derivation repeat 1 and MICE primary imputation, then verify every formal
# preprocessing parameter against the archived main output.
# -------------------------------------------------------------------------
sim <- load_locked_simulation_engine(INPUT$locked_simulation_script)

add_qc("simulation_formal_mode", identical(sim$USE_TEST_MODE, FALSE),
       sim$USE_TEST_MODE, FALSE, "Read directly from locked source")
add_qc("simulation_repeat_count", identical(as.integer(sim$N_REP), 100L),
       sim$N_REP, 100L, "Read directly from locked source")
add_qc("simulation_train_n", identical(as.integer(sim$N_TRAIN), 5000L),
       sim$N_TRAIN, 5000L, "Read directly from locked source")
add_qc("simulation_mice_settings",
       sim$MICE_M == 5L && sim$MICE_MAXIT == 5L && sim$PRIMARY_IMPUTATION_ID == 1L,
       paste0("m=", sim$MICE_M, "; maxit=", sim$MICE_MAXIT,
              "; primary=", sim$PRIMARY_IMPUTATION_ID),
       "m=5; maxit=5; primary=1", "Read directly from locked source")
add_qc("simulation_winsor_settings", identical(as.numeric(sim$WINSOR_PROBS), c(0.01, 0.99)),
       paste(sim$WINSOR_PROBS, collapse = ","), "0.01,0.99",
       "Read directly from locked source")

s1_name <- "S1_pure_severity_continuum"
s1_repeat <- 1L
s1_seed_base <- sim$GLOBAL_SEED + s1_repeat * 1000L +
  match(s1_name, sim$SCENARIOS) * 1e5

# future.apply(future.seed = TRUE) executes the formal main harness under the
# L'Ecuyer-CMRG RNG kind. Reproducing only set.seed() under Mersenne-Twister
# would regenerate a different, albeit valid, S1 dataset.
original_rng_kind <- RNGkind()
RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = original_rng_kind[2L],
        sample.kind = original_rng_kind[3L])
s1_raw <- sim$generate_dataset(
  scenario = s1_name,
  cohort = "derivation",
  n = sim$N_TRAIN,
  repeat_id = s1_repeat,
  seed = s1_seed_base + 1L
)
s1_imp_list <- sim$make_mice_imputations(
  s1_raw,
  scenario = s1_name,
  repeat_id = s1_repeat,
  cohort = "derivation",
  m = sim$MICE_M,
  maxit = sim$MICE_MAXIT,
  seed = s1_seed_base + 101L
)
do.call(RNGkind, as.list(original_rng_kind))
s1_completed <- as.data.table(s1_imp_list[[sim$PRIMARY_IMPUTATION_ID]])
sim$assert_no_missing_features(s1_completed, sim$feature_names, "WP6 regenerated S1")

s1_bounds <- sim$winsor_fit(s1_completed, sim$feature_names)
s1_winsor <- sim$winsor_apply(s1_completed, s1_bounds)
s1_standardized <- sim$standardize_train(s1_winsor, sim$feature_names)

main_results <- readRDS(INPUT$locked_simulation_results)
assert("preprocess_long" %in% names(main_results),
       "Locked simulation result has no preprocess_long table")
archived_preprocess <- as.data.table(main_results$preprocess_long)[
  scenario == s1_name &
    repeat_id == s1_repeat &
    imputation_id == sim$PRIMARY_IMPUTATION_ID &
    analysis_context == "domain1_derivation" &
    cluster_k == 2L,
  .(feature, parameter, archived_value = value)
]
assert(nrow(archived_preprocess) == 4L * sim$P_FEATURES,
       "Archived S1 repeat-1 preprocessing table is incomplete")

generated_preprocess <- rbindlist(list(
  data.table(feature = names(s1_standardized$center), parameter = "center",
             generated_value = as.numeric(s1_standardized$center)),
  data.table(feature = names(s1_standardized$scale), parameter = "scale",
             generated_value = as.numeric(s1_standardized$scale)),
  s1_bounds[, .(feature, parameter = "winsor_low", generated_value = winsor_low)],
  s1_bounds[, .(feature, parameter = "winsor_high", generated_value = winsor_high)]
))

s1_prep_check <- merge(
  archived_preprocess, generated_preprocess,
  by = c("feature", "parameter"), all = TRUE
)
s1_prep_check[, abs_diff := abs(archived_value - generated_value)]
s1_preprocess_max_abs <- max(s1_prep_check$abs_diff, na.rm = TRUE)

add_qc("s1_regenerated_denominator", nrow(s1_completed) == 5000L,
       nrow(s1_completed), 5000L,
       "Formal S1 derivation repeat 1, primary imputation")
add_qc(
  "s1_archived_preprocess_reproduction",
  nrow(s1_prep_check) == 4L * sim$P_FEATURES &&
    all(is.finite(s1_prep_check$abs_diff)) &&
    s1_preprocess_max_abs <= 1e-10,
  format(s1_preprocess_max_abs, scientific = TRUE), "<=1e-10",
  "All 33 features x center/scale/winsor-low/winsor-high"
)

# -------------------------------------------------------------------------
# Point-mass, tail, and distribution-shape summaries.
# -------------------------------------------------------------------------
mimic_z_top5 <- as.data.table(mimic_z_locked[, TOP5, drop = FALSE])
s1_z_top5 <- as.data.table(s1_standardized$Xz[, TOP5, drop = FALSE])

dataset_objects <- list(
  MIMIC = list(
    raw = mimic_original,
    processed_raw = mimic_completed,
    standardized = mimic_z_top5,
    bounds = mimic_bounds,
    n_total = nrow(mimic_completed),
    pipeline = "1/99 winsor before MICE; primary imputation 1; z standardization"
  ),
  S1 = list(
    raw = s1_raw,
    processed_raw = s1_winsor,
    standardized = s1_z_top5,
    bounds = s1_bounds[feature %in% TOP5,
                       .(variable = feature, winsor_low, winsor_high)],
    n_total = nrow(s1_completed),
    pipeline = "MICE; primary imputation 1; 1/99 winsor; z standardization"
  )
)

for (dataset_key in names(dataset_objects)) {
  if ("feature" %in% names(dataset_objects[[dataset_key]]$bounds)) {
    setnames(dataset_objects[[dataset_key]]$bounds, "feature", "variable")
  }
}

summarize_one <- function(dataset_key, variable) {
  obj <- dataset_objects[[dataset_key]]
  target_variable <- variable
  raw_x <- as.numeric(obj$raw[[variable]])
  raw_x <- raw_x[is.finite(raw_x)]
  proc_x <- as.numeric(obj$processed_raw[[variable]])
  proc_x <- proc_x[is.finite(proc_x)]
  z_x <- as.numeric(obj$standardized[[variable]])
  z_x <- z_x[is.finite(z_x)]
  if (anyDuplicated(obj$bounds[["variable"]])) {
    stop("Duplicate winsor-bound rows for dataset ", dataset_key, call. = FALSE)
  }
  bound_idx <- match(target_variable, obj$bounds[["variable"]])
  if (is.na(bound_idx)) {
    stop("No winsor bounds for ", dataset_key, ": ", target_variable,
         call. = FALSE)
  }
  b <- obj$bounds[bound_idx]
  lo <- b$winsor_low[1L]
  hi <- b$winsor_high[1L]
  raw_min <- min(raw_x)
  raw_max <- max(raw_x)
  raw_min_mass <- safe_prop_at(raw_x, raw_min)
  raw_max_mass <- safe_prop_at(raw_x, raw_max)
  low_mass <- safe_prop_at(proc_x, lo)
  high_mass <- safe_prop_at(proc_x, hi)
  zero_mass <- safe_prop_at(proc_x, 0)
  positive <- proc_x[proc_x > 0]
  smallest_positive <- if (length(positive)) min(positive) else NA_real_
  smallest_positive_mass <- safe_prop_at(proc_x, smallest_positive)
  mom <- distribution_moments(z_x)
  q_raw <- stats::quantile(proc_x, probs = c(0.50, 0.75, 0.90, 0.95, 0.99),
                           names = FALSE, type = 7)
  q_z <- stats::quantile(z_x, probs = c(0.01, 0.50, 0.75, 0.90, 0.95, 0.99),
                         names = FALSE, type = 7)
  data.table(
    dataset = DATASET_LABELS[[dataset_key]],
    dataset_key = dataset_key,
    variable = variable,
    variable_label = VARIABLE_LABELS[[variable]],
    n_total = obj$n_total,
    n_raw_observed = length(raw_x),
    pipeline = obj$pipeline,
    raw_min = raw_min,
    raw_min_mass_n = as.integer(raw_min_mass[["n"]]),
    raw_min_mass_prop = raw_min_mass[["prop"]],
    raw_max = raw_max,
    raw_max_mass_n = as.integer(raw_max_mass[["n"]]),
    raw_max_mass_prop = raw_max_mass[["prop"]],
    winsor_low = lo,
    winsor_low_mass_n = as.integer(low_mass[["n"]]),
    winsor_low_mass_prop = low_mass[["prop"]],
    winsor_high = hi,
    winsor_high_mass_n = as.integer(high_mass[["n"]]),
    winsor_high_mass_prop = high_mass[["prop"]],
    processed_min = min(proc_x),
    processed_max = max(proc_x),
    n_unique_raw_observed = uniqueN(raw_x),
    n_unique_processed = uniqueN(proc_x),
    n_unique_standardized = uniqueN(z_x),
    zero_mass_n = as.integer(zero_mass[["n"]]),
    zero_mass_prop = zero_mass[["prop"]],
    smallest_positive_processed_value = smallest_positive,
    smallest_positive_value_mass_n = as.integer(smallest_positive_mass[["n"]]),
    smallest_positive_value_mass_prop = smallest_positive_mass[["prop"]],
    skewness_standardized = mom[["skewness"]],
    excess_kurtosis_standardized = mom[["excess_kurtosis"]],
    processed_p50 = q_raw[1],
    processed_p75 = q_raw[2],
    processed_p90 = q_raw[3],
    processed_p95 = q_raw[4],
    processed_p99 = q_raw[5],
    standardized_p01 = q_z[1],
    standardized_p50 = q_z[2],
    standardized_p75 = q_z[3],
    standardized_p90 = q_z[4],
    standardized_p95 = q_z[5],
    standardized_p99 = q_z[6]
  )
}

# data.table's column scoping is deliberately avoided in the call wrapper.
pointmass_rows <- list()
for (dataset_key in names(dataset_objects)) {
  for (v in TOP5) {
    pointmass_rows[[length(pointmass_rows) + 1L]] <- summarize_one(dataset_key, v)
  }
}
pointmass <- rbindlist(pointmass_rows, use.names = TRUE, fill = TRUE)

expected_bounds <- rbindlist(lapply(names(dataset_objects), function(dataset_key) {
  dataset_key_value <- dataset_key
  b <- copy(dataset_objects[[dataset_key]]$bounds)
  b[variable %in% TOP5, .(
    dataset_key = dataset_key_value,
    variable,
    expected_winsor_low = winsor_low,
    expected_winsor_high = winsor_high
  )]
}))
bound_audit <- merge(
  pointmass[, .(dataset_key, variable, winsor_low, winsor_high,
                processed_min, processed_max,
                winsor_low_mass_n, winsor_high_mass_n)],
  expected_bounds,
  by = c("dataset_key", "variable"),
  all = TRUE
)
bound_tolerance <- 1e-10
add_qc(
  "pointmass_variable_specific_bounds_match",
  nrow(bound_audit) == 10L &&
    all(is.finite(bound_audit$expected_winsor_low)) &&
    all(is.finite(bound_audit$expected_winsor_high)) &&
    max(abs(bound_audit$winsor_low - bound_audit$expected_winsor_low)) <= bound_tolerance &&
    max(abs(bound_audit$winsor_high - bound_audit$expected_winsor_high)) <= bound_tolerance,
  paste0("rows=", nrow(bound_audit), "; max_abs_diff=",
         format(max(c(abs(bound_audit$winsor_low - bound_audit$expected_winsor_low),
                      abs(bound_audit$winsor_high - bound_audit$expected_winsor_high)),
                    na.rm = TRUE), scientific = TRUE)),
  "10 rows; max_abs_diff<=1e-10",
  "Each dataset-variable row must use its own fitted winsor bounds"
)
add_qc(
  "processed_values_within_variable_specific_bounds",
  nrow(bound_audit) == 10L &&
    all(bound_audit$processed_min >= bound_audit$winsor_low - bound_tolerance) &&
    all(bound_audit$processed_max <= bound_audit$winsor_high + bound_tolerance),
  paste0("violations=", sum(
    bound_audit$processed_min < bound_audit$winsor_low - bound_tolerance |
      bound_audit$processed_max > bound_audit$winsor_high + bound_tolerance,
    na.rm = TRUE
  )),
  "violations=0",
  "Post-winsor values must be bounded by the matched variable-specific limits"
)
add_qc(
  "winsor_boundary_masses_observed",
  nrow(bound_audit) == 10L &&
    all(bound_audit$winsor_low_mass_n >= 1L) &&
    all(bound_audit$winsor_high_mass_n >= 1L),
  paste0("positive_rows=", sum(
    bound_audit$winsor_low_mass_n >= 1L & bound_audit$winsor_high_mass_n >= 1L
  ), "/", nrow(bound_audit)),
  "positive_rows=10/10",
  "Both fitted winsor boundaries must be represented after clipping"
)

distance <- rbindlist(lapply(TOP5, function(v) {
  x <- as.numeric(mimic_z_top5[[v]])
  y <- as.numeric(s1_z_top5[[v]])
  support <- sort(unique(c(x, y)))
  ks_d <- max(abs(stats::ecdf(x)(support) - stats::ecdf(y)(support)))
  probs <- c(0.01, 0.50, 0.75, 0.90, 0.95, 0.99)
  qx <- stats::quantile(x, probs = probs, names = FALSE, type = 7)
  qy <- stats::quantile(y, probs = probs, names = FALSE, type = 7)
  data.table(
    variable = v,
    variable_label = VARIABLE_LABELS[[v]],
    mimic_n = length(x),
    s1_n = length(y),
    comparison_scale = "within-dataset z-standardized formal processed values",
    ecdf_sup_distance = ks_d,
    ks_statistic = ks_d,
    ks_p_value_interpretation = "NOT_INTERPRETED_DUE_TO_TIES_AND_LARGE_N",
    wasserstein_1_quantile_approx = common_quantile_wasserstein(x, y),
    wasserstein_grid_n = 10001L,
    mimic_p01 = qx[1], s1_p01 = qy[1], p01_difference_mimic_minus_s1 = qx[1] - qy[1],
    mimic_p50 = qx[2], s1_p50 = qy[2], p50_difference_mimic_minus_s1 = qx[2] - qy[2],
    mimic_p75 = qx[3], s1_p75 = qy[3], p75_difference_mimic_minus_s1 = qx[3] - qy[3],
    mimic_p90 = qx[4], s1_p90 = qy[4], p90_difference_mimic_minus_s1 = qx[4] - qy[4],
    mimic_p95 = qx[5], s1_p95 = qy[5], p95_difference_mimic_minus_s1 = qx[5] - qy[5],
    mimic_p99 = qx[6], s1_p99 = qy[6], p99_difference_mimic_minus_s1 = qx[6] - qy[6],
    interpretation_boundary = paste(
      "Prespecified formal S1 repeat-1 descriptive marginal-shape comparison only;",
      "no pass/fail threshold; point masses do not establish a subtype"
    )
  )
}))

add_qc("pointmass_table_complete", nrow(pointmass) == 10L &&
         all(is.finite(pointmass$skewness_standardized)) &&
         all(is.finite(pointmass$excess_kurtosis_standardized)),
       paste0(nrow(pointmass), " rows"), "10 rows",
       "Five variables x two datasets")
add_qc("distribution_distance_complete", nrow(distance) == 5L &&
         all(is.finite(distance$ecdf_sup_distance)) &&
         all(is.finite(distance$wasserstein_1_quantile_approx)),
       paste0(nrow(distance), " rows"), "5 rows",
       "Five variables")

# -------------------------------------------------------------------------
# Figure source data and publication export.
# Density is normalized within dataset-variable to emphasize shape. ECDF and
# density use exactly the locked/regenerated standardized values.
# -------------------------------------------------------------------------
standardized_long <- rbindlist(list(
  melt(cbind(data.table(row_id = seq_len(nrow(mimic_z_top5))), mimic_z_top5),
       id.vars = "row_id", variable.name = "variable", value.name = "z_value")[,
         `:=`(dataset = DATASET_LABELS[["MIMIC"]])],
  melt(cbind(data.table(row_id = seq_len(nrow(s1_z_top5))), s1_z_top5),
       id.vars = "row_id", variable.name = "variable", value.name = "z_value")[,
         `:=`(dataset = DATASET_LABELS[["S1"]])]
))
standardized_long[, variable_label := factor(
  VARIABLE_LABELS[as.character(variable)],
  levels = unname(VARIABLE_LABELS[TOP5])
)]

global_x <- stats::quantile(
  standardized_long$z_value,
  probs = c(0.001, 0.999), na.rm = TRUE, names = FALSE, type = 7
)
global_x <- c(floor(global_x[1] * 2) / 2, ceiling(global_x[2] * 2) / 2)
grid <- seq(global_x[1], global_x[2], length.out = 512L)

curve_data <- rbindlist(lapply(split(standardized_long,
                                     by = c("variable", "dataset"), keep.by = TRUE),
                               function(d) {
  x <- d$z_value[is.finite(d$z_value)]
  den <- stats::density(x, from = global_x[1], to = global_x[2], n = 512L,
                        adjust = 1, na.rm = TRUE)
  den_y <- den$y / max(den$y)
  e <- stats::ecdf(x)
  rbindlist(list(
    data.table(variable = d$variable[1], dataset = d$dataset[1],
               view = "Relative density", x = den$x, y = den_y),
    data.table(variable = d$variable[1], dataset = d$dataset[1],
               view = "Empirical CDF", x = grid, y = e(grid))
  ))
}))
curve_data[, variable_label := factor(
  VARIABLE_LABELS[as.character(variable)],
  levels = unname(VARIABLE_LABELS[TOP5])
)]
curve_data[, view := factor(view, levels = c("Relative density", "Empirical CDF"))]

theme_pub <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "#333333"),
    strip.background = element_rect(fill = "#F1F3F4", colour = "#B8B8B8", linewidth = 0.35),
    strip.text.x = element_text(size = 9.5, face = "bold", margin = margin(4, 3, 4, 3)),
    strip.text.y = element_text(size = 9.2, face = "bold", angle = 0,
                                margin = margin(3, 5, 3, 5)),
    panel.spacing.x = grid::unit(7, "pt"),
    panel.spacing.y = grid::unit(5, "pt"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.key.width = grid::unit(16, "pt"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", hjust = 0,
                                 margin = margin(b = 7)),
    plot.caption = element_text(size = 7.8, colour = "grey35", hjust = 0,
                                lineheight = 1.12, margin = margin(t = 7)),
    plot.margin = margin(8, 10, 8, 8)
  )

p <- ggplot(curve_data, aes(x = x, y = y, colour = dataset,
                            linetype = dataset, group = dataset)) +
  geom_line(linewidth = 0.65, lineend = "round") +
  facet_grid(rows = vars(variable_label), cols = vars(view), switch = "y") +
  scale_colour_manual(values = COLORS, breaks = names(COLORS)) +
  scale_linetype_manual(values = LINE_TYPES, breaks = names(LINE_TYPES)) +
  scale_x_continuous(limits = global_x, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Distribution fidelity of the five highest-loading variables",
    subtitle = paste(
      "Locked MIMIC-IV primary-imputation matrix versus prespecified formal S1",
      "derivation repeat 1, primary imputation"
    ),
    x = "Within-dataset standardized value",
    y = "Relative density or cumulative probability",
    caption = paste0(
      "Density curves are normalized within each dataset-variable pair.\n",
      "Separate standardization evaluates marginal shape, not clinical scale or joint structure.\n",
      "Point masses do not establish a subtype."
    )
  ) +
  theme_pub

width_in <- 183 / 25.4
height_in <- 210 / 25.4

grDevices::cairo_pdf(OUTPUT$figure_pdf, width = width_in, height = height_in,
                     family = "sans")
print(p)
grDevices::dev.off()

grDevices::png(OUTPUT$figure_png, width = width_in, height = height_in,
               units = "in", res = 600, type = "cairo-png", bg = "white")
print(p)
grDevices::dev.off()

grDevices::tiff(OUTPUT$figure_tiff, width = width_in, height = height_in,
                units = "in", res = 600, compression = "lzw", bg = "white")
print(p)
grDevices::dev.off()

fwrite(pointmass, OUTPUT$pointmass)
fwrite(distance, OUTPUT$distance)
fwrite(curve_data, OUTPUT$figure_source)

qc <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
if (any(!qc$pass)) {
  fwrite(qc, OUTPUT$qc)
  stop("WP6 hard QC failed; see: ", OUTPUT$qc, call. = FALSE)
}

input_hash_end <- rbindlist(lapply(names(INPUT), function(nm) {
  data.table(input_name = nm, sha256_end = sha256_windows(INPUT[[nm]]))
}))
input_hashes <- merge(input_hash_start, input_hash_end, by = "input_name", all = TRUE)
input_hashes[, unchanged_during_run := sha256 == sha256_end]
assert(all(input_hashes$unchanged_during_run), "An authoritative input changed during WP6")

script_hash_end <- sha256_windows(SCRIPT_PATH)
assert(identical(script_hash_start, script_hash_end),
       "The WP6 script changed during execution")

add_qc("input_hashes_unchanged", all(input_hashes$unchanged_during_run),
       paste0(sum(input_hashes$unchanged_during_run), "/", nrow(input_hashes)),
       paste0(nrow(input_hashes), "/", nrow(input_hashes)),
       "SHA-256 before versus after run")
qc <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
fwrite(qc, OUTPUT$qc)
fwrite(input_hashes, OUTPUT$input_checksums)

run_metadata <- data.table(
  run_id = RUN_ID,
  run_mode = RUN_MODE,
  started_from_locked_inputs = TRUE,
  mimic_n = nrow(mimic_locked),
  s1_n = nrow(s1_completed),
  s1_scenario = s1_name,
  s1_repeat_id = s1_repeat,
  s1_primary_imputation_id = sim$PRIMARY_IMPUTATION_ID,
  s1_seed_base = s1_seed_base,
  s1_generation_seed = s1_seed_base + 1L,
  s1_mice_seed = s1_seed_base + 101L,
  s1_rng_kind = "L'Ecuyer-CMRG (formal future.apply worker semantics)",
  empirical_pipeline = dataset_objects$MIMIC$pipeline,
  s1_pipeline = dataset_objects$S1$pipeline,
  pointmass_is_subtype_evidence = FALSE,
  automatic_fidelity_pass_fail_threshold = FALSE,
  s1_scope = "prespecified formal repeat 1 descriptive marginal-shape comparison",
  distance_scope = paste(
    "separately standardized marginal shape, ties, truncation, and tails;",
    "not absolute clinical location/scale or joint dependence"
  ),
  script_sha256_start = script_hash_start,
  script_sha256_end = script_hash_end
)
fwrite(run_metadata, OUTPUT$run_metadata)

log_lines <- c(
  "# WP6 distribution-fidelity audit",
  "",
  paste0("- Run ID: `", RUN_ID, "`"),
  paste0("- Run mode: `", RUN_MODE, "`"),
  paste0("- Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Scope: five highest-loading variables only.",
  "- Empirical source: locked MIMIC-IV primary-imputation standardized matrix.",
  "- S1 source: exact regeneration from the locked formal main-simulation source, prespecified derivation repeat 1, primary imputation.",
  "",
  "## Pipeline identity",
  "",
  paste0("- MIMIC: ", dataset_objects$MIMIC$pipeline, "."),
  paste0("- S1: ", dataset_objects$S1$pipeline, "."),
  paste0("- Locked MIMIC matrix maximum absolute reproduction difference: ",
         format(mimic_matrix_max_abs, scientific = TRUE), "."),
  paste0("- Locked S1 preprocessing maximum absolute reproduction difference: ",
         format(s1_preprocess_max_abs, scientific = TRUE), "."),
  "",
  "## Definitions",
  "",
  "- Raw min/max point mass uses finite observed values before imputation.",
  "- Winsor point mass uses the formal post-winsor, primary-imputation values.",
  "- Smallest positive processed value is the smallest strictly positive post-winsor value; it is not an assay detection limit.",
  "- Skewness and excess kurtosis are calculated on the formal within-dataset standardized values.",
  "- ECDF distance is the two-sample Kolmogorov-Smirnov D statistic; its P value is not interpreted because of ties and large sample sizes.",
  "- Wasserstein-1 distance is approximated by the mean absolute difference between 10,001 matched empirical quantiles on the standardized scale.",
  "",
  "## Interpretation boundary",
  "",
  "This analysis has no prespecified fidelity pass/fail threshold. It is a descriptive marginal-shape comparison using prespecified formal S1 repeat 1. Separately standardized distances assess marginal shape, ties, truncation, and tail behavior, not absolute clinical location or scale, joint dependence, subtype status, or biological mechanism. Point masses, ties, or tail differences do not establish a biological subtype and do not by themselves change the empirical Domain 1 verdict.",
  "",
  "## Outputs",
  "",
  paste0("- `", basename(OUTPUT$pointmass), "`"),
  paste0("- `", basename(OUTPUT$distance), "`"),
  paste0("- `", basename(OUTPUT$figure_pdf), "`"),
  paste0("- `", basename(OUTPUT$figure_png), "`"),
  paste0("- `", basename(OUTPUT$figure_tiff), "`"),
  paste0("- `", basename(OUTPUT$figure_source), "`"),
  "",
  "## QC",
  "",
  paste0("All ", nrow(qc), " hard QC checks passed."),
  ""
)
writeLines(log_lines, OUTPUT$log, useBytes = TRUE)
writeLines(capture.output(sessionInfo()), OUTPUT$session, useBytes = TRUE)

manifest_targets <- unlist(OUTPUT[c(
  "pointmass", "distance", "figure_source", "figure_pdf", "figure_png",
  "figure_tiff", "log", "qc", "session", "input_checksums", "run_metadata"
)], use.names = FALSE)

output_manifest <- rbindlist(lapply(manifest_targets, function(path) {
  data.table(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = sha256_windows(path)
  )
}))
fwrite(output_manifest, OUTPUT$output_manifest)
manifest_sha <- sha256_windows(OUTPUT$output_manifest)

writeLines(c(
  "RUN COMPLETED",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", sha256_windows(SCRIPT_PATH)),
  paste0("output_manifest_sha256=", manifest_sha),
  "interpretation_boundary=prespecified formal S1 repeat-1 descriptive marginal-shape audit; point masses do not establish a subtype"
), OUTPUT$completion, useBytes = TRUE)

if (file.exists(OUTPUT$in_progress)) unlink(OUTPUT$in_progress)

message("WP6 completed successfully.")
message("Point-mass table: ", OUTPUT$pointmass)
message("Distance table: ", OUTPUT$distance)
message("Figure: ", OUTPUT$figure_pdf)
message("Completion marker: ", OUTPUT$completion)
