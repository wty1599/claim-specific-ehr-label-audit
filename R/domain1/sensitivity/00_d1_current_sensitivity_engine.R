#!/usr/bin/env Rscript

# Shared governed engine for the current-rule D1 portions of SA-01 and SA-05.
# This script writes only inside d1_current_rule_sensitivity_20260803 and never
# overwrites the archived July sensitivity results.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(diptest)
  library(Matrix)
  library(mclust)
  library(digest)
  library(future)
  library(future.apply)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv("MIMIC_IV_ROOT", unset = Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())),
  winslash = "/", mustWork = TRUE
)
PACKAGE_ROOT <- file.path(PROJECT_ROOT, "d1_current_rule_sensitivity_20260803")
ANALYSIS <- toupper(Sys.getenv("D1_SENS_ANALYSIS", unset = "SA01"))
RUN_MODE <- tolower(Sys.getenv("D1_SENS_RUN_MODE", unset = "test"))
if (!ANALYSIS %in% c("SA01", "SA05")) {
  stop("D1_SENS_ANALYSIS must be SA01 or SA05.", call. = FALSE)
}
if (!RUN_MODE %in% c("test", "smoke", "formal")) {
  stop("D1_SENS_RUN_MODE must be test, smoke, or formal.", call. = FALSE)
}
IS_FORMAL <- identical(RUN_MODE, "formal")

public_repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
SCRIPT_PATH <- file.path(public_repo, "R/domain1/sensitivity/00_d1_current_sensitivity_engine.R")
WRAPPER_PATH <- file.path(
  public_repo, "R/domain1/sensitivity",
  if (ANALYSIS == "SA01") "01_sa01_sample_size_current_d1.R" else
    "05_sa05_missingness_current_d1.R"
)
UTILS_PATH <- file.path(
  public_repo, "R/domain1/empirical",
  "00_domain1_revision_utils_v4_lloyd.R"
)
DECISION_PATH <- file.path(
  public_repo, "R/domain1/empirical",
  "01_audit_discreteness_v2.R"
)
LOADER_PATH <- file.path(
  public_repo, "R/domain1/empirical/dependencies",
  "00_domain1_locked_engine_v3.R"
)
LOCKED_ENGINE_PATH <- file.path(
  public_repo, "R/simulation/core",
  "01_S7_discrete_outcome_null_formal.R"
)
REGISTRY_PATH <- file.path(PACKAGE_ROOT, "config", "G1_G7_CANONICAL_REGISTRY.csv")
CROSSWALK_PATH <- file.path(PACKAGE_ROOT, "config", "G3_G4_G6_SUBDESIGN_CROSSWALK.csv")

required_files <- c(
  SCRIPT_PATH, WRAPPER_PATH, UTILS_PATH, DECISION_PATH,
  LOADER_PATH, LOCKED_ENGINE_PATH, REGISTRY_PATH, CROSSWALK_PATH
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing governed input(s): ", paste(missing_files, collapse = " | "),
       call. = FALSE)
}

OUT_ROOT <- file.path(PACKAGE_ROOT, paste0("output_", tolower(ANALYSIS), "_", RUN_MODE))
TABLE_DIR <- file.path(OUT_ROOT, "tables")
QC_DIR <- file.path(OUT_ROOT, "qc")
LOG_DIR <- file.path(OUT_ROOT, "logs")
PROV_DIR <- file.path(OUT_ROOT, "provenance")
CHECKPOINT_DIR <- file.path(OUT_ROOT, "checkpoints")
for (d in c(TABLE_DIR, QC_DIR, LOG_DIR, PROV_DIR, CHECKPOINT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(LOG_DIR, paste0(ANALYSIS, "_", RUN_MODE, ".log"))
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
SCRIPT_HASH <- sha256_file(SCRIPT_PATH)
WRAPPER_HASH <- sha256_file(WRAPPER_PATH)
SPEC_HASH <- NA_character_
AUTHORITY_HASH <- digest::digest(
  vapply(c(UTILS_PATH, DECISION_PATH, LOADER_PATH, LOCKED_ENGINE_PATH),
         sha256_file, character(1)),
  algo = "sha256"
)
RUN_SIGNATURE <- digest::digest(list(
  analysis = ANALYSIS, run_mode = RUN_MODE, script = SCRIPT_HASH,
  wrapper = WRAPPER_HASH, spec = SPEC_HASH, authority = AUTHORITY_HASH
), algo = "sha256")

# Governed constants required by the authority utilities and loader.
D1_DIP_FDR_ALPHA <- 0.05
D1_EVALUATION_N <- 1000L
D1_TRAIN_FRACTION <- 0.80
D1_KMEANS_NSTART <- 25L
D1_KMEANS_ITERMAX <- 100L
D1_MAHALANOBIS_RIDGE_MULTIPLIER <- 1e-4
D1_MICE_M <- 5L
D1_MICE_MAXIT <- 5L
D1_PRIMARY_IMPUTATION_ID <- 1L

source(UTILS_PATH, encoding = "UTF-8")
source(LOADER_PATH, encoding = "UTF-8")
source(DECISION_PATH, encoding = "UTF-8")

engine <- load_locked_main_s7_engine(LOCKED_ENGINE_PATH)
stopifnot(length(engine$feature_names) == 33L)

GLOBAL_SEED <- 202608030L
GATE_QUANTILE <- 0.95
QUANTILE_TYPE <- 8L
B_GATE <- switch(RUN_MODE, test = 4L, smoke = 20L, formal = 200L)
B_EVAL <- if (ANALYSIS == "SA01") {
  switch(RUN_MODE, test = 2L, smoke = 5L, formal = 100L)
} else {
  switch(RUN_MODE, test = 1L, smoke = 3L, formal = 50L)
}
B_Q95_BOOT <- switch(RUN_MODE, test = 50L, smoke = 200L, formal = 2000L)
detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores)) detected_cores <- 2L
N_WORKERS <- as.integer(Sys.getenv(
  "D1_SENS_WORKERS",
  unset = as.character(if (RUN_MODE == "test") 1L else
    max(1L, min(6L, detected_cores - 2L)))
))
if (!is.finite(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L

if (ANALYSIS == "SA01") {
  conditions <- data.table(
    condition_order = 1:4,
    condition_id = c("n1000_600", "n2500_1500", "n5000_3000", "n10000_6000"),
    n_train = c(1000L, 2500L, 5000L, 10000L),
    n_eval = c(600L, 1500L, 3000L, 6000L),
    missingness_mechanism = "LOCKED_MAR",
    target_missingness = 0.08
  )
} else {
  conditions <- CJ(
    missingness_mechanism = c("MCAR", "MAR", "MNAR"),
    target_missingness = c(0, 0.15, 0.30), unique = TRUE
  )
  conditions[, `:=`(
    condition_order = seq_len(.N),
    condition_id = sprintf(
      "%s_miss%02d", missingness_mechanism,
      as.integer(round(100 * target_missingness))
    ),
    n_train = 4000L,
    n_eval = 1000L
  )]
  setcolorder(conditions, c(
    "condition_order", "condition_id", "n_train", "n_eval",
    "missingness_mechanism", "target_missingness"
  ))
}

controls <- data.table(
  scenario_order = 1:2,
  reader_id = c("G1", "G7"),
  engine_scenario = c(
    "S1_pure_severity_continuum",
    "S7_discrete_structure_outcome_null"
  ),
  dgm_family = c("S1_continuum", "S7_like_K2"),
  truth = c("negative", "positive"),
  expected_state = c("NO_DISCRETE_EVIDENCE", "DISCRETE_EVIDENCE")
)

# At zero missingness the MCAR/MAR/MNAR labels describe the same scientific
# condition. Use one canonical condition index for every stochastic step so
# that initialisation and gate calibration cannot create artificial mechanism
# differences when no values are missing.
canonical_condition_index <- function(condition) {
  if (ANALYSIS == "SA05" && isTRUE(condition$target_missingness == 0)) {
    return(conditions[
      missingness_mechanism == "MCAR" & target_missingness == 0,
      condition_order
    ][1L])
  }
  condition$condition_order
}

write_csv_atomic(conditions, file.path(PROV_DIR, "conditions.csv"))
write_csv_atomic(controls, file.path(PROV_DIR, "controls.csv"))
file.copy(REGISTRY_PATH, file.path(PROV_DIR, basename(REGISTRY_PATH)), overwrite = TRUE)
file.copy(CROSSWALK_PATH, file.path(PROV_DIR, basename(CROSSWALK_PATH)), overwrite = TRUE)

run_contract <- data.table(
  analysis = ANALYSIS,
  run_mode = RUN_MODE,
  B_gate = B_GATE,
  B_evaluation = B_EVAL,
  B_q95_bootstrap = B_Q95_BOOT,
  gate_quantile = GATE_QUANTILE,
  quantile_type = QUANTILE_TYPE,
  mice_m = D1_MICE_M,
  mice_maxit = D1_MICE_MAXIT,
  primary_imputation = D1_PRIMARY_IMPUTATION_ID,
  kmeans_algorithm = "Lloyd",
  kmeans_nstart = D1_KMEANS_NSTART,
  kmeans_itermax = D1_KMEANS_ITERMAX,
  shape_contract = "heldout_fisher_discriminant_bh6",
  reference_scope = "processed_space_fixed_k2_n_specific",
  workers = N_WORKERS,
  global_seed = GLOBAL_SEED,
  script_sha256 = SCRIPT_HASH,
  wrapper_sha256 = WRAPPER_HASH,
  spec_sha256 = SPEC_HASH,
  authority_bundle_sha256 = AUTHORITY_HASH,
  run_signature = RUN_SIGNATURE,
  r_version = R.version.string
)
write_csv_atomic(run_contract, file.path(PROV_DIR, "run_contract.csv"))

randomized_pit <- function(x, seed) {
  set.seed(seed)
  n <- length(x)
  r_min <- rank(x, ties.method = "min")
  r_max <- rank(x, ties.method = "max")
  u <- (r_min - 1 + runif(n) * (r_max - r_min + 1)) / n
  pmin(pmax(u, 1 / (2 * n)), 1 - 1 / (2 * n))
}

fit_gaussian_copula <- function(X, seed) {
  X <- as.matrix(X)
  U <- vapply(seq_len(ncol(X)), function(j) {
    randomized_pit(X[, j], seed + 1009L * j)
  }, numeric(nrow(X)))
  colnames(U) <- colnames(X)
  Z <- qnorm(U)
  R_raw <- cor(Z)
  eig <- eigen(R_raw, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 1e-10) {
    near <- Matrix::nearPD(R_raw, corr = TRUE, keepDiag = TRUE)
    R_used <- as.matrix(near$mat)
    nearpd_used <- TRUE
    adjustment <- norm(R_used - R_raw, type = "F")
  } else {
    R_used <- R_raw
    nearpd_used <- FALSE
    adjustment <- 0
  }
  list(
    R = R_used,
    margins = setNames(lapply(seq_len(ncol(X)), function(j) sort(X[, j])),
                      colnames(X)),
    nearpd_used = nearpd_used,
    nearpd_frobenius = adjustment,
    min_eigen_raw = min(eig)
  )
}

simulate_gaussian_copula <- function(fit, n, seed) {
  set.seed(seed)
  p <- ncol(fit$R)
  Z <- matrix(rnorm(n * p), nrow = n) %*% chol(fit$R)
  U <- pnorm(Z)
  X <- vapply(seq_len(p), function(j) {
    margin <- fit$margins[[j]]
    idx <- pmax(1L, pmin(length(margin), ceiling(U[, j] * length(margin))))
    margin[idx]
  }, numeric(n))
  colnames(X) <- names(fit$margins)
  X
}

stdz <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

set.seed(91002L)
missing_weights <- sample(seq(0.60, 1.40, length.out = length(engine$feature_names)))
missing_weights <- missing_weights / mean(missing_weights)
names(missing_weights) <- engine$feature_names

calibrate_intercept <- function(score, target) {
  if (target <= 0) return(-Inf)
  uniroot(function(a) mean(plogis(a + score)) - target,
          c(-40, 40), tol = 1e-12)$root
}

target_by_variable <- function(target) {
  if (target == 0) {
    return(setNames(rep(0, length(engine$feature_names)), engine$feature_names))
  }
  p <- pmin(target * missing_weights, 0.60)
  p <- p * target / mean(p)
  p <- pmin(p, 0.65)
  p <- p * target / mean(p)
  setNames(p, engine$feature_names)
}

missing_scores <- function(dt, variable, mechanism, reader_id) {
  mar <- 0.55 * stdz(dt$sofa) + 0.25 * stdz(dt$age) +
    0.20 * stdz(dt$aki_stage) + 0.15 * (dt$sex - mean(dt$sex))
  own <- stdz(dt[[variable]])
  subtype_term <- if (reader_id == "G7") {
    0.20 * stdz(as.integer(as.factor(dt$true_subtype)))
  } else 0
  mnar <- 0.80 * own + 0.25 * stdz(dt$severity_z) + subtype_term
  selected <- switch(
    mechanism,
    MCAR = rep(0, nrow(dt)), MAR = mar, MNAR = mnar,
    stop("Unknown missingness mechanism: ", mechanism, call. = FALSE)
  )
  list(selected = selected, mar = mar, own = own)
}

inject_condition_missingness <- function(dt, reader_id, mechanism, target,
                                         mask_seed) {
  out <- copy(as.data.table(dt))
  targets <- target_by_variable(target)
  masks <- matrix(FALSE, nrow(out), length(engine$feature_names))
  var_rows <- vector("list", length(engine$feature_names))
  for (j in seq_along(engine$feature_names)) {
    v <- engine$feature_names[j]
    tv <- targets[[v]]
    scores <- missing_scores(dt, v, mechanism, reader_id)
    if (tv == 0) {
      prob <- rep(0, nrow(dt))
      mask <- rep(FALSE, nrow(dt))
    } else {
      alpha <- calibrate_intercept(scores$selected, tv)
      prob <- plogis(alpha + scores$selected)
      set.seed(mask_seed + 1000L * j)
      mask <- runif(nrow(dt)) < prob
    }
    masks[, j] <- mask
    if (any(mask)) set(out, which(mask), v, NA_real_)
    var_rows[[j]] <- data.table(
      variable = v,
      target_variable_missingness = tv,
      realized_variable_missingness = mean(mask),
      missing_cells = sum(mask),
      eligible_cells = length(mask),
      selected_signal_correlation = suppressWarnings(
        cor(as.numeric(mask), scores$selected)
      )
    )
  }
  list(
    data = out,
    overall = data.table(
      realized_overall_missingness = mean(masks),
      missing_cells = sum(masks),
      eligible_cells = length(masks)
    ),
    variable = rbindlist(var_rows)
  )
}

add_baseline <- function(dt, cohort, scenario, repeat_id) {
  dt <- copy(as.data.table(dt))
  baseline <- as.data.table(engine$generate_baseline_covariates(dt$severity_z))
  out <- cbind(dt, baseline)
  out[, `:=`(
    patient_id = sprintf("%s_%s_%04d_%06d", scenario, cohort,
                         repeat_id, seq_len(.N)),
    cohort = cohort,
    true_scenario = scenario,
    repeat_id = repeat_id
  )]
  out
}

complete_imputations <- function(dt, scenario, repeat_id, cohort, seed) {
  use_vars <- intersect(c(engine$feature_names, engine$baseline_vars), names(dt))
  if (!anyNA(dt[, ..use_vars])) {
    return(list(
      imputations = replicate(D1_MICE_M, copy(dt[, ..use_vars]), simplify = FALSE),
      warning_count = 0L,
      warning_text = ""
    ))
  }
  warnings <- character()
  imps <- withCallingHandlers(
    engine$make_mice_imputations(
      dt, scenario, repeat_id, cohort,
      m = D1_MICE_M, maxit = D1_MICE_MAXIT, seed = seed
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    imputations = imps,
    warning_count = length(warnings),
    warning_text = paste(unique(warnings), collapse = " | ")
  )
}

make_custom_missingness_pair <- function(control, condition, repeat_id,
                                         complete_seed, mask_seed, mice_seed) {
  complete <- simulate_complete_d1_pair(
    engine = engine,
    dgm_family = control$dgm_family,
    n_train = condition$n_train,
    n_evaluation = condition$n_eval,
    delta = if (control$reader_id == "G7") 3.2 else 0,
    seed = complete_seed
  )
  train0 <- add_baseline(
    complete$train_dt, "derivation", control$engine_scenario, repeat_id
  )
  eval0 <- add_baseline(
    complete$evaluation_dt, "same_distribution_evaluation",
    control$engine_scenario, repeat_id
  )
  trm <- inject_condition_missingness(
    train0, control$reader_id, condition$missingness_mechanism,
    condition$target_missingness, mask_seed + 100000L
  )
  evm <- inject_condition_missingness(
    eval0, control$reader_id, condition$missingness_mechanism,
    condition$target_missingness, mask_seed + 200000L
  )
  tri <- complete_imputations(
    trm$data, control$engine_scenario, repeat_id, "derivation",
    mice_seed + 1000L
  )
  evi <- complete_imputations(
    evm$data, control$engine_scenario, repeat_id,
    "same_distribution_evaluation", mice_seed + 2000L
  )
  train <- as.data.table(tri$imputations[[D1_PRIMARY_IMPUTATION_ID]])
  evaluation <- as.data.table(evi$imputations[[D1_PRIMARY_IMPUTATION_ID]])
  prep <- preprocess_train_evaluation(engine, train, evaluation)
  list(
    X_train = prep$X_train,
    X_evaluation = prep$X_evaluation,
    evaluation_true_labels = eval0$true_subtype,
    realized_train = trm$overall$realized_overall_missingness,
    realized_evaluation = evm$overall$realized_overall_missingness,
    realized_overall = (trm$overall$missing_cells + evm$overall$missing_cells) /
      (trm$overall$eligible_cells + evm$overall$eligible_cells),
    mice_warning_count = tri$warning_count + evi$warning_count,
    mice_warning_text = paste(
      c(tri$warning_text, evi$warning_text)[
        nzchar(c(tri$warning_text, evi$warning_text))
      ], collapse = " | "
    )
  )
}

seed_components <- function(condition, control, repeat_id, stage) {
  scenario_i <- control$scenario_order
  if (ANALYSIS == "SA01") {
    complete_group <- condition$condition_order
    mask_group <- condition$condition_order
  } else {
    # Complete data are paired across all mechanism/intensity settings.
    complete_group <- 1L
    mechanism_i <- match(condition$missingness_mechanism, c("MCAR", "MAR", "MNAR"))
    mask_group <- mechanism_i
  }
  stage_i <- match(stage, c("reference", "evaluation"))
  complete_seed <- GLOBAL_SEED + scenario_i * 10000000L +
    complete_group * 100000L + repeat_id * 100L + stage_i
  mask_seed <- GLOBAL_SEED + scenario_i * 10000000L +
    mask_group * 1000000L + repeat_id * 100L + stage_i
  target_i <- if (condition$target_missingness == 0) 0L else
    as.integer(round(condition$target_missingness * 100))
  mechanism_i <- if (condition$target_missingness == 0) 0L else
    match(condition$missingness_mechanism, c("MCAR", "MAR", "MNAR"), nomatch = 0L)
  mice_seed <- GLOBAL_SEED + scenario_i * 10000000L +
    mechanism_i * 1000000L + target_i * 10000L + repeat_id * 100L + stage_i
  list(complete = complete_seed, mask = mask_seed, mice = mice_seed)
}

make_source_pair <- function(condition, control, repeat_id, stage) {
  seeds <- seed_components(condition, control, repeat_id, stage)
  if (ANALYSIS == "SA01") {
    pair <- simulate_mice_d1_pair(
      engine = engine,
      scenario = control$engine_scenario,
      repeat_id = repeat_id,
      seed = seeds$complete,
      n_train = condition$n_train,
      n_evaluation = condition$n_eval
    )
    return(list(
      X_train = pair$X_train,
      X_evaluation = pair$X_evaluation,
      evaluation_true_labels = pair$evaluation_dt$true_subtype,
      realized_train = NA_real_,
      realized_evaluation = NA_real_,
      realized_overall = NA_real_,
      mice_warning_count = pair$mice_warning_count,
      mice_warning_text = pair$mice_warning_text
    ))
  }
  make_custom_missingness_pair(
    control, condition, repeat_id,
    seeds$complete, seeds$mask, seeds$mice
  )
}

checkpoint_path <- function(stage, id) {
  d <- file.path(CHECKPOINT_DIR, stage)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, paste0(id, "_", substr(RUN_SIGNATURE, 1, 12), ".rds"))
}

cached_compute <- function(path, signature, fun) {
  if (file.exists(path)) {
    x <- readRDS(path)
    if (identical(attr(x, "signature"), signature)) return(x)
  }
  x <- fun()
  attr(x, "signature") <- signature
  saveRDS_atomic(x, path)
  x
}

diagnose_pair <- function(X_train, X_evaluation, true_labels, seed) {
  domain1_external_diagnostics(
    X_train, X_evaluation,
    evaluation_labels = true_labels,
    seed = seed
  )
}

log_msg("Starting ", ANALYSIS, " current-rule D1; mode=", RUN_MODE,
        "; conditions=", nrow(conditions), "; B_gate=", B_GATE,
        "; B_eval=", B_EVAL, "; workers=", N_WORKERS)

# Fit one condition- and control-matched empirical-margin N1 reference before
# any evaluation repeat is interpreted.
options(future.globals.maxSize = 4 * 1024^3)
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
if (N_WORKERS > 1L) future::plan(multisession, workers = N_WORKERS)

reference_tasks <- CJ(
  condition_index = seq_len(nrow(conditions)),
  scenario_index = seq_len(nrow(controls)), unique = TRUE
)

run_reference_task <- function(i) {
    ci <- reference_tasks$condition_index[i]
    si <- reference_tasks$scenario_index[i]
    condition <- conditions[ci]
    control <- controls[si]
    id <- paste(condition$condition_id, control$reader_id, sep = "__")
    sig <- digest::digest(list(RUN_SIGNATURE, "reference", condition, control))
    path <- checkpoint_path("reference", id)
    ref <- cached_compute(path, sig, function() {
      pair <- make_source_pair(condition, control, 0L, "reference")
      scientific_ci <- canonical_condition_index(condition)
      fit <- fit_gaussian_copula(
        pair$X_train,
        seed = GLOBAL_SEED + scientific_ci * 10000L + si * 1000L + 77L
      )
      list(
        fit = fit,
        reference = data.table(
          condition_id = condition$condition_id,
          reader_id = control$reader_id,
          n_train = condition$n_train,
          n_eval = condition$n_eval,
          missingness_mechanism = condition$missingness_mechanism,
          target_missingness = condition$target_missingness,
          realized_missingness = pair$realized_overall,
          mice_warning_count = pair$mice_warning_count,
          copula_nearpd_used = fit$nearpd_used,
          copula_nearpd_frobenius = fit$nearpd_frobenius,
          copula_min_eigen_raw = fit$min_eigen_raw
        )
      )
    })
    list(id = id, reference = ref$reference, fit = ref$fit)
}
reference_list <- future.apply::future_lapply(
  seq_len(nrow(reference_tasks)), run_reference_task,
  future.seed = NULL, future.scheduling = 1
)
reference_rows <- setNames(lapply(reference_list, `[[`, "reference"),
                           vapply(reference_list, `[[`, character(1), "id"))
reference_fits <- setNames(lapply(reference_list, `[[`, "fit"),
                           vapply(reference_list, `[[`, character(1), "id"))
reference_table <- rbindlist(reference_rows, fill = TRUE)
write_csv_atomic(reference_table, file.path(TABLE_DIR, "D1_reference_fit_audit.csv"))

gate_tasks <- CJ(
  condition_index = seq_len(nrow(conditions)),
  scenario_index = seq_len(nrow(controls)),
  gate_repeat = seq_len(B_GATE), unique = TRUE
)
gate_tasks[, task_id := sprintf(
  "gate_%s_%s_%04d",
  conditions$condition_id[condition_index],
  controls$reader_id[scenario_index], gate_repeat
)]

run_gate_task <- function(i) {
  task <- gate_tasks[i]
  condition <- conditions[task$condition_index]
  control <- controls[task$scenario_index]
  key <- paste(condition$condition_id, control$reader_id, sep = "__")
  sig <- digest::digest(list(RUN_SIGNATURE, "gate", task))
  path <- checkpoint_path("gate", task$task_id)
  cached_compute(path, sig, function() {
    scientific_ci <- canonical_condition_index(condition)
    seed0 <- GLOBAL_SEED + scientific_ci * 1000000L +
      task$scenario_index * 100000L + task$gate_repeat * 10L
    Xtr <- simulate_gaussian_copula(
      reference_fits[[key]], condition$n_train, seed0 + 1L
    )
    Xev <- simulate_gaussian_copula(
      reference_fits[[key]], condition$n_eval, seed0 + 2L
    )
    diag <- diagnose_pair(Xtr, Xev, NULL, seed0 + 3L)
    cbind(
      task, condition,
      control[, .(scenario_order, reader_id, truth)],
      diag
    )
  })
}

gate_list <- future.apply::future_lapply(
  seq_len(nrow(gate_tasks)), run_gate_task,
  future.seed = NULL, future.scheduling = 2
)
gate_rows <- rbindlist(gate_list, fill = TRUE)

gate_summary <- gate_rows[, .(
  B_gate = .N,
  q95 = quantile(
    mahalanobis_centroid_delta, GATE_QUANTILE,
    type = QUANTILE_TYPE, names = FALSE
  ),
  delta_mean = mean(mahalanobis_centroid_delta),
  delta_sd = sd(mahalanobis_centroid_delta),
  delta_mcse = sd(mahalanobis_centroid_delta) / sqrt(.N),
  fisher_shape_alert_rate = mean(dip_discriminant_p_fdr < D1_DIP_FDR_ALPHA)
), by = .(
  condition_order, condition_id, n_train, n_eval,
  missingness_mechanism, target_missingness, scenario_order, reader_id, truth
)]

boot_rows <- rbindlist(lapply(seq_len(nrow(gate_summary)), function(i) {
  g <- gate_summary[i]
  x <- gate_rows[
    condition_id == g$condition_id & reader_id == g$reader_id,
    mahalanobis_centroid_delta
  ]
  scientific_ci <- canonical_condition_index(conditions[
    condition_id == g$condition_id
  ])
  set.seed(
    GLOBAL_SEED + 80000000L + scientific_ci * 100L + g$scenario_order
  )
  q <- replicate(B_Q95_BOOT, quantile(
    sample(x, length(x), replace = TRUE), GATE_QUANTILE,
    type = QUANTILE_TYPE, names = FALSE
  ))
  data.table(
    condition_id = g$condition_id,
    reader_id = g$reader_id,
    q95_bootstrap_B = B_Q95_BOOT,
    q95_bootstrap_se = sd(q),
    q95_bootstrap_low = quantile(q, 0.025, names = FALSE),
    q95_bootstrap_high = quantile(q, 0.975, names = FALSE)
  )
}))
gate_summary <- merge(
  gate_summary, boot_rows,
  by = c("condition_id", "reader_id"), all.x = TRUE, sort = FALSE
)
write_csv_atomic(gate_rows, file.path(TABLE_DIR, "D1_N1_gate_by_repeat.csv"))
write_csv_atomic(gate_summary, file.path(TABLE_DIR, "D1_N1_gate_summary.csv"))

eval_tasks <- CJ(
  condition_index = seq_len(nrow(conditions)),
  scenario_index = seq_len(nrow(controls)),
  repeat_id = seq_len(B_EVAL), unique = TRUE
)
eval_tasks[, task_id := sprintf(
  "eval_%s_%s_%04d",
  conditions$condition_id[condition_index],
  controls$reader_id[scenario_index], repeat_id
)]

run_eval_task <- function(i) {
  task <- eval_tasks[i]
  condition <- conditions[task$condition_index]
  control <- controls[task$scenario_index]
  sig <- digest::digest(list(RUN_SIGNATURE, "evaluation", task))
  path <- checkpoint_path("evaluation", task$task_id)
  cached_compute(path, sig, function() {
    warnings <- character()
    tryCatch(withCallingHandlers({
      pair <- make_source_pair(condition, control, task$repeat_id, "evaluation")
      scientific_ci <- canonical_condition_index(condition)
      seed0 <- GLOBAL_SEED + scientific_ci * 1000000L +
        task$scenario_index * 100000L + task$repeat_id * 10L + 7L
      diag <- diagnose_pair(
        pair$X_train, pair$X_evaluation,
        if (control$reader_id == "G7") pair$evaluation_true_labels else NULL,
        seed0
      )
      gate <- gate_summary[
        condition_id == condition$condition_id &
          reader_id == control$reader_id, q95
      ]
      if (length(gate) != 1L || !is.finite(gate)) {
        stop("Condition-specific q95 lookup failed.", call. = FALSE)
      }
      verdict <- audit_discreteness_v2(list(
        separation = diag$mahalanobis_centroid_delta,
        gate = gate,
        shape_reject = diag$dip_discriminant_p_fdr < D1_DIP_FDR_ALPHA,
        shape_source = diag$shape_source,
        reference_scope = diag$reference_scope,
        partition_source = diag$partition_source,
        kmeans_algorithm = diag$kmeans_algorithm,
        kmeans_nstart = diag$kmeans_nstart,
        kmeans_itermax = diag$kmeans_itermax
      ))
      if (verdict$state == "NOT_EVALUATED") {
        stop("Governed wrapper returned NOT_EVALUATED: ",
             verdict$failure_code, call. = FALSE)
      }
      cbind(
        task, condition,
        control[, .(
          scenario_order, reader_id, engine_scenario, truth, expected_state
        )],
        diag,
        data.table(
          separation_gate_q95 = gate,
          separation_alert = verdict$separation_alert,
          shape_alert = verdict$shape_alert,
          state = verdict$state,
          strict_correct = verdict$state == control$expected_state,
          alert_positive = verdict$state == "DISCRETE_EVIDENCE",
          realized_train_missingness = pair$realized_train,
          realized_evaluation_missingness = pair$realized_evaluation,
          realized_overall_missingness = pair$realized_overall,
          mice_warning_count = pair$mice_warning_count,
          mice_warning_text = pair$mice_warning_text,
          task_warning_count = length(warnings),
          task_warning_text = paste(unique(warnings), collapse = " | "),
          status = "completed",
          error = NA_character_,
          decision_rule_version = verdict$decision_rule_version
        )
      )
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      cbind(
        task, condition,
        control[, .(
          scenario_order, reader_id, engine_scenario, truth, expected_state
        )],
        data.table(
          state = "FAILED", strict_correct = FALSE,
          alert_positive = NA, status = "failed",
          error = conditionMessage(e),
          task_warning_count = length(warnings),
          task_warning_text = paste(unique(warnings), collapse = " | ")
        )
      )
    })
  })
}

if (N_WORKERS > 1L) future::plan(multisession, workers = N_WORKERS)
eval_list <- future.apply::future_lapply(
  seq_len(nrow(eval_tasks)), run_eval_task,
  future.seed = NULL, future.scheduling = 2
)
future::plan(sequential)
eval_rows <- rbindlist(eval_list, fill = TRUE)
write_csv_atomic(eval_rows, file.path(TABLE_DIR, "D1_current_rule_by_repeat.csv"))

wilson_cols <- function(x, n, prefix) {
  ci <- wilson_interval(x, n)
  setNames(data.table(x = x, n = n, rate = if (n) x / n else NA_real_,
                      low = ci[1L], high = ci[2L]),
           paste0(prefix, c("_events", "_n", "_rate", "_low", "_high")))
}

operating_rows <- eval_rows[, {
  valid <- status == "completed"
  z <- .SD[valid]
  strict_n <- sum(z$strict_correct)
  discrete_n <- sum(z$state == "DISCRETE_EVIDENCE")
  no_discrete_n <- sum(z$state == "NO_DISCRETE_EVIDENCE")
  inconclusive_n <- sum(z$state == "INCONCLUSIVE")
  cbind(
    data.table(
      expected_n = .N,
      valid_n = nrow(z),
      failed_n = sum(!valid),
      expected_state = unique(expected_state)
    ),
    wilson_cols(strict_n, nrow(z), "strict_correct"),
    wilson_cols(discrete_n, nrow(z), "discrete"),
    wilson_cols(no_discrete_n, nrow(z), "no_discrete"),
    wilson_cols(inconclusive_n, nrow(z), "inconclusive"),
    wilson_cols(strict_n, .N, "failure_equals_error_correct")
  )
}, by = .(
  condition_order, condition_id, n_train, n_eval,
  missingness_mechanism, target_missingness, scenario_order, reader_id, truth
)]
setorder(operating_rows, condition_order, scenario_order)
write_csv_atomic(
  operating_rows,
  file.path(TABLE_DIR, "D1_current_rule_operating_characteristics.csv")
)

continuous_summary <- eval_rows[status == "completed", .(
  n = .N,
  separation_mean = mean(mahalanobis_centroid_delta),
  separation_sd = sd(mahalanobis_centroid_delta),
  separation_mcse = sd(mahalanobis_centroid_delta) / sqrt(.N),
  separation_p025 = quantile(mahalanobis_centroid_delta, 0.025),
  separation_p975 = quantile(mahalanobis_centroid_delta, 0.975),
  separation_margin_mean = mean(
    mahalanobis_centroid_delta - separation_gate_q95
  ),
  fisher_p_fdr_mean = mean(dip_discriminant_p_fdr),
  fisher_p_fdr_median = median(dip_discriminant_p_fdr),
  fisher_p_fdr_p025 = quantile(dip_discriminant_p_fdr, 0.025),
  fisher_p_fdr_p975 = quantile(dip_discriminant_p_fdr, 0.975),
  true_label_ari_mean_descriptive = if (all(is.na(refit_vs_input_label_ari_descriptive)))
    NA_real_ else mean(refit_vs_input_label_ari_descriptive, na.rm = TRUE),
  true_label_ari_p025_descriptive = if (all(is.na(refit_vs_input_label_ari_descriptive)))
    NA_real_ else quantile(refit_vs_input_label_ari_descriptive, 0.025, na.rm = TRUE),
  true_label_ari_p975_descriptive = if (all(is.na(refit_vs_input_label_ari_descriptive)))
    NA_real_ else quantile(refit_vs_input_label_ari_descriptive, 0.975, na.rm = TRUE)
), by = .(
  condition_order, condition_id, n_train, n_eval,
  missingness_mechanism, target_missingness, scenario_order, reader_id, truth
)]
write_csv_atomic(
  continuous_summary,
  file.path(TABLE_DIR, "D1_current_rule_continuous_summary.csv")
)

failure_table <- eval_rows[status != "completed", .N, by = .(
  condition_id, reader_id, error
)]
if (!nrow(failure_table)) {
  failure_table <- data.table(
    condition_id = character(), reader_id = character(),
    error = character(), N = integer()
  )
}
write_csv_atomic(failure_table, file.path(QC_DIR, "failure_reasons.csv"))

if (ANALYSIS == "SA05") {
  missingness_by_repeat <- eval_rows[status == "completed", .(
    condition_id, reader_id, repeat_id, missingness_mechanism,
    target_missingness, realized_train_missingness,
    realized_evaluation_missingness, realized_overall_missingness
  )]
  write_csv_atomic(
    missingness_by_repeat,
    file.path(TABLE_DIR, "SA05_realized_missingness_by_repeat.csv")
  )
  missingness_summary <- missingness_by_repeat[, .(
    n = .N,
    realized_mean = mean(realized_overall_missingness),
    realized_sd = sd(realized_overall_missingness),
    realized_mcse = sd(realized_overall_missingness) / sqrt(.N),
    realized_p025 = quantile(realized_overall_missingness, 0.025),
    realized_p975 = quantile(realized_overall_missingness, 0.975),
    mean_absolute_target_error = mean(abs(
      realized_overall_missingness - target_missingness
    ))
  ), by = .(reader_id, missingness_mechanism, target_missingness)]
  write_csv_atomic(
    missingness_summary,
    file.path(TABLE_DIR, "SA05_realized_missingness_summary.csv")
  )
  zero_keys <- c("reader_id", "repeat_id")
  zero_metrics <- c(
    "realized_overall_missingness", "mahalanobis_centroid_delta",
    "dip_discriminant_p_fdr", "separation_gate_q95", "separation_alert",
    "shape_alert", "state", "strict_correct"
  )
  zero_long <- eval_rows[target_missingness == 0 & status == "completed"]
  zero_checks <- lapply(zero_metrics, function(metric) {
    z <- dcast(
      zero_long,
      reader_id + repeat_id ~ missingness_mechanism,
      value.var = metric
    )
    z[, metric := metric]
    z[, exact_metric_equivalence :=
        (is.na(MCAR) & is.na(MAR) & is.na(MNAR)) |
        (!is.na(MCAR) & !is.na(MAR) & !is.na(MNAR) &
           as.character(MCAR) == as.character(MAR) &
           as.character(MCAR) == as.character(MNAR))]
    z
  })
  zero_metric_qc <- rbindlist(zero_checks, fill = TRUE)
  zero_wide <- zero_metric_qc[, .(
    exact_zero_equivalence = all(exact_metric_equivalence),
    metrics_checked = paste(metric, collapse = ";")
  ), by = zero_keys]
  write_csv_atomic(
    zero_metric_qc,
    file.path(QC_DIR, "SA05_zero_missingness_equivalence.csv")
  )
  monotonic <- dcast(
    missingness_by_repeat,
    reader_id + missingness_mechanism + repeat_id ~ target_missingness,
    value.var = "realized_overall_missingness"
  )
  setnames(monotonic, c("0", "0.15", "0.3"), c("miss0", "miss15", "miss30"))
  monotonic[, monotonic_pass := miss0 <= miss15 & miss15 <= miss30]
  write_csv_atomic(
    monotonic,
    file.path(QC_DIR, "SA05_missingness_monotonicity.csv")
  )
}

qc <- rbindlist(list(
  data.table(
    item = "gate_row_count",
    pass = nrow(gate_rows) == nrow(conditions) * nrow(controls) * B_GATE,
    observed = nrow(gate_rows),
    expected = nrow(conditions) * nrow(controls) * B_GATE
  ),
  data.table(
    item = "evaluation_row_count",
    pass = nrow(eval_rows) == nrow(conditions) * nrow(controls) * B_EVAL,
    observed = nrow(eval_rows),
    expected = nrow(conditions) * nrow(controls) * B_EVAL
  ),
  data.table(
    item = "finite_q95",
    pass = all(is.finite(gate_summary$q95)),
    observed = sum(is.finite(gate_summary$q95)),
    expected = nrow(gate_summary)
  ),
  data.table(
    item = "no_failed_evaluation_repeats",
    pass = all(eval_rows$status == "completed"),
    observed = sum(eval_rows$status == "completed"),
    expected = nrow(eval_rows)
  ),
  data.table(
    item = "no_NOT_EVALUATED_states",
    pass = !any(eval_rows$state == "NOT_EVALUATED", na.rm = TRUE),
    observed = sum(eval_rows$state == "NOT_EVALUATED", na.rm = TRUE),
    expected = 0L
  ),
  data.table(
    item = "truth_ari_excluded_from_verdict",
    pass = all(eval_rows$oracle_ari_used_in_verdict[eval_rows$status == "completed"] == FALSE),
    observed = sum(eval_rows$oracle_ari_used_in_verdict[eval_rows$status == "completed"]),
    expected = 0L
  )
), fill = TRUE)

if (ANALYSIS == "SA05") {
  qc <- rbind(qc,
    data.table(
      item = "zero_missingness_exact_equivalence",
      pass = all(zero_wide$exact_zero_equivalence),
      observed = sum(zero_wide$exact_zero_equivalence),
      expected = nrow(zero_wide)
    ),
    data.table(
      item = "missingness_intensity_monotonicity",
      pass = all(monotonic$monotonic_pass),
      observed = sum(monotonic$monotonic_pass),
      expected = nrow(monotonic)
    ),
    data.table(
      item = "realized_missingness_fidelity",
      pass = all(missingness_summary$mean_absolute_target_error <= 0.01),
      observed = max(missingness_summary$mean_absolute_target_error),
      expected = 0.01
    ), fill = TRUE
  )
}
write_csv_atomic(qc, file.path(QC_DIR, "critical_qc.csv"))

if (RUN_MODE == "smoke" && all(qc$pass)) {
  writeLines(c(
    "SMOKE_READY_FOR_FORMAL",
    paste0("analysis=", ANALYSIS),
    paste0("script_sha256=", SCRIPT_HASH),
    paste0("spec_sha256=", SPEC_HASH),
    paste0("run_signature=", RUN_SIGNATURE),
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  ), file.path(PROV_DIR, "SMOKE_READY_FOR_FORMAL.txt"))
}

writeLines(capture.output(sessionInfo()), file.path(PROV_DIR, "sessionInfo.txt"))
writeLines(c(
  paste0("analysis=", ANALYSIS),
  paste0("run_mode=", RUN_MODE),
  paste0("critical_qc_pass=", all(qc$pass)),
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "Archived pre-Lloyd SA-01/SA-05 outputs were not overwritten."
), file.path(PROV_DIR, "completion_status.txt"))

log_msg("Completed ", ANALYSIS, " ", RUN_MODE,
        "; critical_qc_pass=", all(qc$pass),
        "; output=", OUT_ROOT)

if (!all(qc$pass)) {
  stop("Critical QC failed; output is not eligible for formal interpretation.",
       call. = FALSE)
}
