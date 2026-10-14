public_repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
#!/usr/bin/env Rscript

# WP7 Lloyd-aligned direct negative/positive controls for a Gaussian-copula
# empirical-margin candidate null. This script runs Domain 1 only. It never
# replaces locked results and never uses test/smoke output as a formal result.

options(stringsAsFactors = FALSE, warn = 1)

PACKAGE_STARTUP_WARNINGS <- character()
withCallingHandlers(
  suppressPackageStartupMessages({
    library(data.table)
    library(mice)
    library(diptest)
    library(Matrix)
    library(mclust)
    library(ggplot2)
    library(patchwork)
    library(digest)
  }),
  warning = function(w) {
    PACKAGE_STARTUP_WARNINGS <<- c(
      PACKAGE_STARTUP_WARNINGS, conditionMessage(w)
    )
    invokeRestart("muffleWarning")
  }
)
PACKAGE_STARTUP_WARNINGS <- unique(PACKAGE_STARTUP_WARNINGS)
PACKAGE_STARTUP_WARNING_CLASS <- if (!length(PACKAGE_STARTUP_WARNINGS)) {
  "none"
} else if (all(grepl(
  "^package .+ was built under R version [0-9.]+$",
  PACKAGE_STARTUP_WARNINGS
))) {
  "package_built_under_newer_R_patch_release"
} else {
  stop(
    "Unclassified package-startup warning(s): ",
    paste(PACKAGE_STARTUP_WARNINGS, collapse = " | "), call. = FALSE
  )
}

required_pkgs <- c(
  "data.table", "mice", "diptest", "Matrix", "mclust", "ggplot2",
  "patchwork", "digest", "scales"
)
missing_pkgs <- required_pkgs[!vapply(
  required_pkgs, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_pkgs)) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Run contract
# -------------------------------------------------------------------------

RUN_MODE <- tolower(Sys.getenv("D1_WP7_RUN_MODE", unset = "test"))
if (!RUN_MODE %in% c("test", "smoke", "formal")) {
  stop("D1_WP7_RUN_MODE must be test, smoke, or formal.", call. = FALSE)
}
IS_FORMAL <- identical(RUN_MODE, "formal")
RUN_ID <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid())

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT)) stop("Set EHR_AUDIT_WORK_ROOT to the authorized workspace.")
AUDIT_ROOT <- file.path(PROJECT_ROOT, "domain1_lloyd_alignment_20260720")
SCRIPT_PATH <- file.path(public_repo, "R/domain1/empirical/02_copula_controls_formal.R")




V3_ROOT <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations",
  "Domain1_revision_v3_observable_heldout_20260717"
)
V4_UTILS <- file.path(public_repo, "R/domain1/empirical/00_domain1_revision_utils_v4_lloyd.R")
STRICT_WRAPPER <- file.path(public_repo, "R/domain1/empirical/01_audit_discreteness_v2.R")
V3_LOADER <- file.path(public_repo, "R/domain1/empirical/dependencies/00_domain1_locked_engine_v3.R")
V3_CONFIG <- file.path(public_repo, "R/domain1/empirical/dependencies/00_domain1_revision_config.R")
LOCKED_S7_ENGINE <- file.path(public_repo, "R/simulation/core/01_S7_discrete_outcome_null_formal.R")
LOCKED_MAIN_ENGINE <- file.path(public_repo, "R/simulation/core/00_main_S1_S5_formal.R")
LOCKED_MAIN_RESULTS <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations", "main_extended_validation_no_pac_formal",
  "output", "raw_results", "simulation_four_domain_all_results_mice.rds"
)
LOCKED_S7_RESULTS <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations", "S7_audit_20260712", "extracted",
  "simulation_v9_scenario_S7_discrete_structure_outcome_null", "output",
  "raw_results", "simulation_four_domain_all_results_mice.rds"
)
LOCKED_EMPIRICAL_MATRIX <- file.path(
  PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"
)
LOCKED_EMPIRICAL_LABELS <- file.path(
  PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"
)
LOCKED_EMPIRICAL_MAHALANOBIS <- file.path(
  PROJECT_ROOT, "output", "qc_domain1_mice",
  "posthoc_locked_k2_pooled_within_mahalanobis_20260717.csv"
)



TABLE_DIR <- file.path(AUDIT_ROOT, "tables")
FIGURE_DIR <- file.path(AUDIT_ROOT, "figures")
LOG_DIR <- file.path(AUDIT_ROOT, "logs")
PROVENANCE_DIR <- file.path(AUDIT_ROOT, "provenance")
CHECKPOINT_ROOT <- file.path(AUDIT_ROOT, "checkpoints")

preflight_sha256 <- function(path) {
  digest::digest(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
}

if (IS_FORMAL) {
  OUT_TABLE_DIR <- TABLE_DIR
  OUT_FIGURE_DIR <- FIGURE_DIR
  OUT_LOG_DIR <- LOG_DIR
  OUT_PROVENANCE_DIR <- PROVENANCE_DIR
  TASK_CHECKPOINT_DIR <- file.path(CHECKPOINT_ROOT, "wp7", "29f_lloyd_formal_tasks")
} else {
  NONFORMAL_ROOT <- file.path(
    AUDIT_ROOT, "nonformal_runs", RUN_MODE, RUN_ID
  )
  OUT_TABLE_DIR <- file.path(NONFORMAL_ROOT, "tables")
  OUT_FIGURE_DIR <- file.path(NONFORMAL_ROOT, "figures")
  OUT_LOG_DIR <- file.path(NONFORMAL_ROOT, "logs")
  OUT_PROVENANCE_DIR <- file.path(NONFORMAL_ROOT, "provenance")
  TASK_CHECKPOINT_DIR <- file.path(NONFORMAL_ROOT, "tasks")
}

for (d in c(
  OUT_TABLE_DIR, OUT_FIGURE_DIR, OUT_LOG_DIR, OUT_PROVENANCE_DIR,
  TASK_CHECKPOINT_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

B_CALIBRATION <- if (IS_FORMAL) 200L else if (RUN_MODE == "smoke") 8L else 2L
B_EVALUATION <- if (IS_FORMAL) 200L else if (RUN_MODE == "smoke") 8L else 2L
SOURCE_N <- 5000L
EVALUATION_N <- 1000L
S7_DELTA <- 3.2
KMEANS_NSTART <- 25L
KMEANS_ITERMAX <- 100L
RIDGE_MULTIPLIER <- 1e-4
DIP_ALPHA <- 0.05
NULL_QUANTILE <- 0.95
MICE_M <- 5L
MICE_MAXIT <- 5L
PRIMARY_IMPUTATION_ID <- 1L
PIT_BASE_SEED <- 314159265L
TASK_BASE_SEED <- 2026071721L
RULE_VERSION <- "wp7_gaussian_copula_direct_controls_lloyd_v2_20260722"
CHECKPOINT_SCHEMA <- "29f_wp7_lloyd_schema_20260722_1"

detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores)) detected_cores <- 2L
default_workers <- if (IS_FORMAL) {
  max(1L, min(8L, detected_cores - 2L))
} else if (RUN_MODE == "smoke") {
  2L
} else {
  1L
}
N_WORKERS <- as.integer(Sys.getenv(
  "D1_WP7_WORKERS", unset = as.character(default_workers)
))
if (!is.finite(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L

# Names used by the sourced v3 utilities.
D1_KMEANS_NSTART <- KMEANS_NSTART
D1_KMEANS_ITERMAX <- KMEANS_ITERMAX
D1_MAHALANOBIS_RIDGE_MULTIPLIER <- RIDGE_MULTIPLIER
D1_DIP_FDR_ALPHA <- DIP_ALPHA
D1_MICE_M <- MICE_M
D1_MICE_MAXIT <- MICE_MAXIT
D1_PRIMARY_IMPUTATION_ID <- PRIMARY_IMPUTATION_ID

for (path in c(
  SCRIPT_PATH, V3_CONFIG,
  V4_UTILS, STRICT_WRAPPER, V3_LOADER,
  LOCKED_MAIN_ENGINE,
  LOCKED_MAIN_RESULTS, LOCKED_S7_ENGINE, LOCKED_S7_RESULTS,
  LOCKED_EMPIRICAL_MATRIX, LOCKED_EMPIRICAL_LABELS,
  LOCKED_EMPIRICAL_MAHALANOBIS
)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}
INPUT <- list(
  locked_empirical_matrix = LOCKED_EMPIRICAL_MATRIX,
  locked_empirical_labels = LOCKED_EMPIRICAL_LABELS,
  locked_empirical_mahalanobis = LOCKED_EMPIRICAL_MAHALANOBIS,
  locked_main_engine = LOCKED_MAIN_ENGINE,
  locked_main_results_reference = LOCKED_MAIN_RESULTS,
  locked_s7_engine = LOCKED_S7_ENGINE,
  locked_s7_results_reference = LOCKED_S7_RESULTS,
  v3_config_reference = V3_CONFIG,
  v4_lloyd_utils = V4_UTILS,
  strict_decision_wrapper = STRICT_WRAPPER,
  v3_loader = V3_LOADER
)

INPUT_ROLE <- c(
  locked_empirical_matrix = "computational_input",
  locked_empirical_labels = "provenance_reference_not_used_for_inference",
  locked_empirical_mahalanobis = "locked_authority_reference_not_used_for_inference",
  locked_main_engine = "computational_engine_definition",
  locked_main_results_reference = "provenance_reference_not_used_for_inference",
  locked_s7_engine = "computational_engine_definition",
  locked_s7_results_reference = "provenance_reference_not_used_for_inference",
  v3_config_reference = "provenance_reference_not_used_for_inference",
  v4_lloyd_utils = "computational_utility",
  strict_decision_wrapper = "governed_decision_contract_reference",
  v3_loader = "computational_engine_loader"
)
stopifnot(identical(names(INPUT), names(INPUT_ROLE)))

OUTPUT <- list(
  negative = file.path(OUT_TABLE_DIR, "Table_D1_copula_negative_control_LLOYD.csv"),
  positive = file.path(OUT_TABLE_DIR, "Table_D1_copula_positive_control_LLOYD.csv"),
  multimodality = file.path(
    OUT_TABLE_DIR, "Table_D1_copula_null_multimodality_LLOYD.csv"
  ),
  by_repeat = file.path(OUT_TABLE_DIR, "29f_copula_controls_by_repeat_LLOYD.csv"),
  marginal = file.path(
    OUT_TABLE_DIR, "Table_D1_copula_empirical_marginal_dip_descriptive.csv"
  ),
  source_marginal = file.path(
    OUT_TABLE_DIR, "29f_copula_control_source_marginal_dip_by_repeat.csv"
  ),
  figure_source = file.path(
    OUT_TABLE_DIR, "29f_Figure_D1_copula_controls_source_data.csv"
  ),
  figure_pdf = file.path(OUT_FIGURE_DIR, "Figure_D1_copula_controls.pdf"),
  figure_png = file.path(OUT_FIGURE_DIR, "Figure_D1_copula_controls.png"),
  figure_tiff = file.path(OUT_FIGURE_DIR, "Figure_D1_copula_controls.tiff"),
  decision = file.path(OUT_LOG_DIR, "D1_COPULA_GO_NO_GO.txt"),
  log = file.path(OUT_LOG_DIR, "29f_copula_log_LLOYD.md"),
  qc = file.path(OUT_LOG_DIR, "29f_copula_QC_LLOYD.csv"),
  failures = file.path(OUT_LOG_DIR, "29f_failure_log_LLOYD.csv"),
  top_level_failure = file.path(
    OUT_LOG_DIR, "29f_top_level_failure.txt"
  ),
  session = file.path(OUT_LOG_DIR, "29f_sessionInfo.txt"),
  task_manifest = file.path(OUT_PROVENANCE_DIR, "29f_task_manifest.csv"),
  seed_registry = file.path(OUT_PROVENANCE_DIR, "29f_global_seed_registry_LLOYD.csv"),
  run_contract = file.path(OUT_PROVENANCE_DIR, "29f_WP7_run_contract.csv"),
  checkpoint_inventory = file.path(
    OUT_PROVENANCE_DIR, "29f_WP7_checkpoint_inventory.csv"
  ),
  input_checksums = file.path(
    OUT_PROVENANCE_DIR, "29f_input_checksums_LLOYD.csv"
  ),
  run_metadata = file.path(OUT_PROVENANCE_DIR, "29f_WP7_run_metadata.csv"),
  output_manifest = file.path(
    OUT_PROVENANCE_DIR, "29f_output_sha256_manifest_LLOYD.csv"
  ),
  in_progress = file.path(TASK_CHECKPOINT_DIR, "29f_run_in_progress.txt"),
  partial = file.path(OUT_LOG_DIR, "29f_run_partial_or_qc_failure.txt"),
  completion = if (IS_FORMAL) {
    file.path(AUDIT_ROOT, "29f_run_completed_LLOYD.ok")
  } else {
    file.path(dirname(OUT_TABLE_DIR), "29f_nonformal_run_completed.ok")
  }
)

if (IS_FORMAL) {
  protected <- unique(unlist(OUTPUT, use.names = FALSE))
  present <- protected[file.exists(protected)]
  if (length(present)) {
    stop(
      "Formal WP7 outputs already exist; refusing to overwrite:\n",
      paste(present, collapse = "\n"), call. = FALSE
    )
  }
}

message(
  "WP7 run contract: mode=", RUN_MODE,
  "; calibration/evaluation B per control=", B_CALIBRATION, "/", B_EVALUATION,
  "; source/evaluation n=", SOURCE_N, "/", EVALUATION_N,
  "; algorithm=Lloyd; nstart=", KMEANS_NSTART,
  "; iter.max=", KMEANS_ITERMAX,
  "; workers=", N_WORKERS
)
message(
  "WP7 palette: source=#2E2E2E; copula-null=#728995; ",
  "positive-control accent=#C77A30; reference=grey70"
)
message("Figure size: 183 x 178 mm; PNG/TIFF 600 dpi")

# -------------------------------------------------------------------------
# Integrity helpers
# -------------------------------------------------------------------------

sha256_file <- function(path) {
  digest::digest(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
}

md5_file <- function(path) {
  digest::digest(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    algo = "md5", file = TRUE, serialize = FALSE
  )
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (anyDuplicated(names(x))) {
    stop("Duplicate output columns: ",
         paste(unique(names(x)[duplicated(names(x))]), collapse = ", "),
         call. = FALSE)
  }
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(as.data.table(x), tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize: ", path, call. = FALSE)
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize: ", path, call. = FALSE)
  invisible(path)
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize: ", path, call. = FALSE)
  invisible(path)
}

wilson <- function(x, n, conf.level = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0 || x < 0 || x > n) {
    return(c(low = NA_real_, high = NA_real_))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, centre - half), high = min(1, centre + half))
}

newcombe_difference <- function(x1, n1, x0, n0, conf.level = 0.95) {
  p1 <- x1 / n1
  p0 <- x0 / n0
  ci1 <- wilson(x1, n1, conf.level)
  ci0 <- wilson(x0, n0, conf.level)
  difference <- p1 - p0
  lower <- difference - sqrt((p1 - ci1["low"])^2 +
                               (ci0["high"] - p0)^2)
  upper <- difference + sqrt((ci1["high"] - p1)^2 +
                               (p0 - ci0["low"])^2)
  c(estimate = difference, low = max(-1, lower), high = min(1, upper))
}

wp7_control_state <- function(shape_positive, delta_positive) {
  out <- rep("WP7_CONTROL_ONLY_DISCORDANT", length(shape_positive))
  out[shape_positive & delta_positive] <- "WP7_CONTROL_ONLY_BOTH_POSITIVE"
  out[!shape_positive & !delta_positive] <- "WP7_CONTROL_ONLY_BOTH_NEGATIVE"
  out
}

numeric_summary <- function(x, prefix) {
  x <- x[is.finite(x)]
  n <- length(x)
  vals <- c(
    n = n,
    mean = if (n) mean(x) else NA_real_,
    sd = if (n > 1L) sd(x) else NA_real_,
    mcse = if (n > 1L) sd(x) / sqrt(n) else NA_real_,
    median = if (n) median(x) else NA_real_,
    p025 = if (n) quantile(x, 0.025, type = 8, names = FALSE) else NA_real_,
    p975 = if (n) quantile(x, 0.975, type = 8, names = FALSE) else NA_real_
  )
  as.data.table(as.list(setNames(vals, paste0(prefix, "_", names(vals)))))
}

qc_rows <- list()
add_qc <- function(check, pass, observed, expected, detail = "") {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table(
    check = check, pass = isTRUE(pass), observed = as.character(observed),
    expected = as.character(expected), detail = as.character(detail)
  )
}

script_hash_start <- sha256_file(SCRIPT_PATH)
input_hash_start <- rbindlist(lapply(names(INPUT), function(nm) {
  path <- INPUT[[nm]]
  data.table(
    input_name = nm,
    input_role = unname(INPUT_ROLE[[nm]]),
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(path), md5 = md5_file(path), bytes = file.info(path)$size,
    modified = format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %Z")
  )
}))

EXPECTED_SHA256 <- c(
  locked_empirical_matrix = "9896ad4906207ca910e54eb187507037106bff212299586141908f387cfe6279",
  locked_empirical_labels = "159576d0b2796e9786db0eb4c8747af987be5e53af200419e3e37567c1d875f9",
  locked_empirical_mahalanobis = "28c33cc58306a2662bb1dbe17edf49d7c5707b897af1eff01599a3d88148e72b",
  locked_main_engine = "dc1cd7f4c5beb2cf0f22a95a7c582fc2770f600491fec0ab5740d64ba6d20488",
  locked_s7_engine = "b55cf80d2ac4fdb166cc1b7c448629837d2f6364199d079848ded8dd9a0fc54a",
  v4_lloyd_utils = "a58e730f70a16cb95b9b3fd9d241e121270066e2f36e61826ac76b11b4014281",
  strict_decision_wrapper = "99a7b6e54c889185d0b0e91fd499a7677fcfdce7cd74c4c8ff016864de456c3c",
  v3_loader = "e897194e5ca6cf051a9b7456c8749ee9c99077c07f91075f89850249d7f51c91"
)
observed_expected <- input_hash_start[input_name %in% names(EXPECTED_SHA256)]
observed_expected[, expected_sha256 := unname(EXPECTED_SHA256[input_name])]
if (nrow(observed_expected) != length(EXPECTED_SHA256) ||
    any(tolower(observed_expected$sha256) !=
        tolower(observed_expected$expected_sha256))) {
  mismatch <- observed_expected[tolower(sha256) != tolower(expected_sha256)]
  stop(
    "Frozen authority hash mismatch before task execution: ",
    paste(mismatch$input_name, collapse = ", "), call. = FALSE
  )
}

OLD_ERROR_OPTION <- getOption("error")
RUN_STAGE <- "load_hashed_compute_tools"
top_level_failure_handler <- function() {
  failure_lines <- c(
    "WP7 TOP-LEVEL FAILURE",
    "run_status=RUN_FAILED",
    "decision=NOT_EVALUATED",
    paste0("run_id=", RUN_ID),
    paste0("run_mode=", RUN_MODE),
    paste0("stage=", RUN_STAGE),
    paste0("error=", trimws(geterrmessage())),
    paste0("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "No scientific go/no-go decision is valid."
  )
  for (path in c(OUTPUT$top_level_failure, OUTPUT$partial)) {
    try({
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      cat(paste(failure_lines, collapse = "\n"), "\n",
          file = path, append = file.exists(path))
    }, silent = TRUE)
  }
  invisible(NULL)
}
top_level_abort <- function() {
  try(top_level_failure_handler(), silent = TRUE)
  options(error = NULL)
  quit(save = "no", status = 1L, runLast = FALSE)
}
options(error = top_level_abort)

# Hashes are frozen before executable analysis helpers are loaded.
source(V4_UTILS, encoding = "UTF-8")
source(V3_LOADER, encoding = "UTF-8")
contract_probe_X <- rbind(
  cbind(seq(-2, -1, length.out = 50), seq(-1, -2, length.out = 50)),
  cbind(seq(1, 2, length.out = 50), seq(2, 1, length.out = 50))
)
contract_probe <- safe_kmeans(
  contract_probe_X, centers = 2L, seed = 20260722L,
  nstart = 25L, iter.max = 100L
)
contract_probe_audit <- kmeans_audit_record(contract_probe, "probe")
if (inherits(contract_probe, "d1_kmeans_error") ||
    !identical(contract_probe_audit$probe_algorithm_used, "Lloyd") ||
    !identical(contract_probe_audit$probe_fallback_used, FALSE) ||
    contract_probe_audit$probe_warning_count != 0L ||
    contract_probe_audit$probe_nstart != 25L ||
    contract_probe_audit$probe_itermax != 100L ||
    !identical(contract_probe_audit$probe_failed, FALSE)) {
  stop("Loaded clustering utility does not satisfy the frozen Lloyd contract.",
       call. = FALSE)
}
rm(contract_probe_X, contract_probe, contract_probe_audit)
RUN_STAGE <- "initialization_complete"

atomic_write_lines(c(
  "WP7 RUN IN PROGRESS",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", script_hash_start),
  paste0("prerun_spec_sha256=", NA_character_)
), OUTPUT$in_progress)

input_contract_hash <- digest::digest(
  input_hash_start[order(input_name), .(input_name, sha256, md5)],
  algo = "sha256", serialize = TRUE
)
run_contract <- data.table(
  run_id = RUN_ID,
  run_mode = RUN_MODE,
  reportable_formal_result = IS_FORMAL,
  output_root = normalizePath(
    if (IS_FORMAL) AUDIT_ROOT else NONFORMAL_ROOT,
    winslash = "/", mustWork = TRUE
  ),
  checkpoint_root = normalizePath(
    TASK_CHECKPOINT_DIR, winslash = "/", mustWork = TRUE
  ),
  calibration_repetitions_per_control = B_CALIBRATION,
  evaluation_repetitions_per_control = B_EVALUATION,
  source_n = SOURCE_N,
  evaluation_n = EVALUATION_N,
  kmeans_algorithm = "Lloyd",
  kmeans_nstart = KMEANS_NSTART,
  kmeans_itermax = KMEANS_ITERMAX,
  kmeans_fallback_allowed = FALSE,
  mice_m = MICE_M,
  mice_maxit = MICE_MAXIT,
  primary_imputation_id = PRIMARY_IMPUTATION_ID,
  s7_delta = S7_DELTA,
  s7_prevalence = "0.50/0.50",
  ridge_multiplier = RIDGE_MULTIPLIER,
  dip_alpha = DIP_ALPHA,
  null_quantile = NULL_QUANTILE,
  task_base_seed = TASK_BASE_SEED,
  pit_base_seed = PIT_BASE_SEED,
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA,
  script_sha256 = script_hash_start,
  frozen_spec_sha256 = NA_character_,
  utils_sha256 = sha256_file(V4_UTILS),
  wrapper_sha256 = sha256_file(STRICT_WRAPPER),
  input_contract_sha256 = input_contract_hash,
  verdict_scope = "WP7_CONTROL_ONLY_NOT_WP8_FORMAL_EMPIRICAL_GATE",
  package_startup_warning_class = PACKAGE_STARTUP_WARNING_CLASS,
  package_startup_warning_count = length(PACKAGE_STARTUP_WARNINGS),
  package_startup_warning_text = paste(PACKAGE_STARTUP_WARNINGS, collapse = " | ")
)
atomic_fwrite(run_contract, OUTPUT$run_contract)
run_contract_sha256 <- sha256_file(OUTPUT$run_contract)

# -------------------------------------------------------------------------
# Copula and held-out diagnostics
# -------------------------------------------------------------------------

randomized_pit <- function(x, seed) {
  set.seed(seed)
  n <- length(x)
  r_min <- rank(x, ties.method = "min")
  r_max <- rank(x, ties.method = "max")
  u <- (r_min - 1 + runif(n) * (r_max - r_min + 1)) / n
  pmin(pmax(u, 1 / (2 * n)), 1 - 1 / (2 * n))
}

fit_gaussian_copula <- function(X, pit_seed) {
  X <- as.matrix(X)
  p <- ncol(X)
  U <- vapply(seq_len(p), function(j) {
    randomized_pit(X[, j], pit_seed + j * 1009L)
  }, numeric(nrow(X)))
  colnames(U) <- colnames(X)
  Z <- qnorm(U)
  R_raw <- cor(Z)
  eig_raw <- eigen(R_raw, symmetric = TRUE, only.values = TRUE)$values
  nearpd_used <- min(eig_raw) <= 1e-10
  if (nearpd_used) {
    near <- Matrix::nearPD(R_raw, corr = TRUE, keepDiag = TRUE)
    R_used <- as.matrix(near$mat)
    adjustment <- norm(R_used - R_raw, type = "F")
  } else {
    R_used <- R_raw
    adjustment <- 0
  }
  eig_used <- eigen(R_used, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(R_used)) || min(eig_used) <= 0) {
    stop("Copula correlation matrix is not positive definite.", call. = FALSE)
  }
  margins <- lapply(seq_len(p), function(j) sort(X[, j]))
  names(margins) <- colnames(X)
  list(
    R_raw = R_raw, R_used = R_used, sorted_margins = margins,
    pit_seed = pit_seed, min_eigen_raw = min(eig_raw),
    min_eigen_used = min(eig_used), nearpd_used = nearpd_used,
    nearpd_frobenius = adjustment,
    pit_method = "randomized_Rueschendorf_distribution_transform"
  )
}

simulate_gaussian_copula <- function(fit, n, seed) {
  set.seed(seed)
  p <- ncol(fit$R_used)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(fit$R_used)
  U <- pnorm(Z)
  X <- vapply(seq_len(p), function(j) {
    margin <- fit$sorted_margins[[j]]
    idx <- pmax(1L, pmin(length(margin), ceiling(U[, j] * length(margin))))
    margin[idx]
  }, numeric(n))
  colnames(X) <- names(fit$sorted_margins)
  X
}

max_point_mass <- function(x) {
  tab <- table(x, useNA = "no")
  if (!length(tab)) return(NA_real_)
  max(tab) / length(x)
}

copula_fidelity <- function(source_X, null_X) {
  source_X <- as.matrix(source_X)
  null_X <- as.matrix(null_X)
  upper <- upper.tri(cor(source_X, method = "spearman"), diag = FALSE)
  rs <- cor(source_X, method = "spearman")
  rn <- cor(null_X, method = "spearman")
  data.table(
    mean_abs_column_mean_difference = mean(abs(
      colMeans(source_X) - colMeans(null_X)
    )),
    mean_abs_column_sd_difference = mean(abs(
      apply(source_X, 2L, sd) - apply(null_X, 2L, sd)
    )),
    spearman_upper_triangle_rmse = sqrt(mean((rs[upper] - rn[upper])^2)),
    mean_abs_max_point_mass_difference = mean(abs(
      vapply(seq_len(ncol(source_X)), function(j) max_point_mass(source_X[, j]), numeric(1)) -
        vapply(seq_len(ncol(null_X)), function(j) max_point_mass(null_X[, j]), numeric(1))
    ))
  )
}

marginal_dip_table <- function(X, control, repeat_id, role) {
  X <- as.matrix(X)
  rows <- rbindlist(lapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    test <- diptest::dip.test(x)
    data.table(
      control = control, repeat_id = repeat_id, role = role,
      feature = colnames(X)[j], n = length(x), unique_n = uniqueN(x),
      max_point_mass = max_point_mass(x), dip_statistic = unname(test$statistic),
      dip_p_raw = test$p.value
    )
  }))
  rows[, dip_p_bh := p.adjust(dip_p_raw, method = "BH")]
  rows[, bh_reject_0_05 := dip_p_bh < DIP_ALPHA]
  rows[, used_in_go_no_go := FALSE]
  rows
}

heldout_diagnostics <- function(X_train, X_evaluation,
                                evaluation_true_labels = NULL, seed) {
  X_train <- as.matrix(X_train)
  X_evaluation <- as.matrix(X_evaluation)
  if (ncol(X_train) != ncol(X_evaluation) ||
      any(!is.finite(X_train)) || any(!is.finite(X_evaluation))) {
    stop("Invalid train/evaluation matrices.", call. = FALSE)
  }
  fit <- safe_kmeans(
    X_train, centers = 2L, seed = seed, nstart = KMEANS_NSTART,
    iter.max = KMEANS_ITERMAX
  )
  if (inherits(fit, "d1_kmeans_error")) {
    stop("K-means failed: ", fit$error, call. = FALSE)
  }
  audit <- kmeans_audit_record(fit, "kmeans")
  centres <- fit$centers
  evaluation_labels <- nearest_centroid_labels(X_evaluation, centres)
  separation <- pooled_within_mahalanobis(
    X_train, fit$cluster, centres, ridge_multiplier = RIDGE_MULTIPLIER
  )
  if (!is.finite(separation$delta) || !is.finite(separation$ridge) ||
      any(!is.finite(evaluation_labels))) {
    stop("Non-finite separation or assignment output.", call. = FALSE)
  }
  direction <- separation$fisher_direction
  if (sum(direction^2) <= 1e-12) direction[1L] <- 1
  discriminant <- as.numeric(X_evaluation %*% direction)

  pc_fit <- prcomp(X_train, center = TRUE, scale. = FALSE)
  pc_eval <- predict(pc_fit, newdata = X_evaluation)
  if (ncol(pc_eval) < 5L || any(!is.finite(discriminant)) ||
      any(!is.finite(pc_eval[, seq_len(5L), drop = FALSE]))) {
    stop("Non-finite held-out discriminant/PCA projection.", call. = FALSE)
  }
  p_raw <- c(
    discriminant = diptest::dip.test(discriminant)$p.value,
    vapply(seq_len(5L), function(j) {
      diptest::dip.test(pc_eval[, j])$p.value
    }, numeric(1))
  )
  names(p_raw)[-1L] <- paste0("pc", seq_len(5L))
  p_bh <- p.adjust(p_raw, method = "BH")
  if (length(p_raw) != 6L || length(p_bh) != 6L ||
      any(!is.finite(p_raw)) || any(!is.finite(p_bh))) {
    stop("Non-finite held-out dip P value.", call. = FALSE)
  }
  any_axis_bh_reject <- any(p_bh < DIP_ALPHA)
  if (length(any_axis_bh_reject) != 1L || is.na(any_axis_bh_reject)) {
    stop("Invalid held-out any-axis rejection event.", call. = FALSE)
  }

  true_ari <- if (is.null(evaluation_true_labels)) {
    NA_real_
  } else {
    mclust::adjustedRandIndex(
      as.integer(as.factor(evaluation_true_labels)), evaluation_labels
    )
  }

  out <- data.table(
    delta = separation$delta,
    mahalanobis_ridge = separation$ridge,
    holdout_discriminant_dip_p_raw = unname(p_raw["discriminant"]),
    holdout_discriminant_dip_p_bh = unname(p_bh["discriminant"]),
    holdout_pc1_dip_p_raw = unname(p_raw["pc1"]),
    holdout_pc1_dip_p_bh = unname(p_bh["pc1"]),
    holdout_pc2_dip_p_raw = unname(p_raw["pc2"]),
    holdout_pc2_dip_p_bh = unname(p_bh["pc2"]),
    holdout_pc3_dip_p_raw = unname(p_raw["pc3"]),
    holdout_pc3_dip_p_bh = unname(p_bh["pc3"]),
    holdout_pc4_dip_p_raw = unname(p_raw["pc4"]),
    holdout_pc4_dip_p_bh = unname(p_bh["pc4"]),
    holdout_pc5_dip_p_raw = unname(p_raw["pc5"]),
    holdout_pc5_dip_p_bh = unname(p_bh["pc5"]),
    holdout_any_axis_bh_reject = any_axis_bh_reject,
    holdout_assigned_cluster1_proportion = mean(evaluation_labels == 1L),
    holdout_true_label_ari_descriptive = true_ari,
    dip_source = "independent evaluation observations",
    pca_source = "training rotation frozen to independent evaluation",
    discriminant_source = "training Fisher direction frozen to independent evaluation",
    reference_scope =
      "processed_space_fixed_k2_empirical_margin_gaussian_copula",
    shape_source = "heldout_fisher_discriminant_bh6",
    partition_source = paste(
      "Lloyd K2 fitted on training data; nearest-centroid assignment on",
      "held-out evaluation"
    ),
    oracle_ari_used_in_verdict = FALSE
  )
  cbind(out, audit)
}

prefix_columns <- function(dt, prefix) {
  dt <- copy(as.data.table(dt))
  setnames(dt, names(dt), paste0(prefix, "_", names(dt)))
  dt
}

checkpoint_signature <- function(task) {
  digest::digest(list(
    rule_version = RULE_VERSION,
    checkpoint_schema = CHECKPOINT_SCHEMA,
    script_sha256 = script_hash_start,
    input_sha256 = input_hash_start[, .(input_name, sha256)],
    task = task,
    run_mode = RUN_MODE,
    b_calibration = B_CALIBRATION,
    b_evaluation = B_EVALUATION,
    source_n = SOURCE_N,
    evaluation_n = EVALUATION_N,
    s7_delta = S7_DELTA,
    algorithm = "Lloyd",
    nstart = KMEANS_NSTART,
    itermax = KMEANS_ITERMAX,
    mice_m = MICE_M,
    mice_maxit = MICE_MAXIT,
    primary_imputation_id = PRIMARY_IMPUTATION_ID,
    ridge = RIDGE_MULTIPLIER,
    dip_alpha = DIP_ALPHA,
    pit_base_seed = PIT_BASE_SEED
  ), algo = "sha256")
}

read_checkpoint <- function(path, signature) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(
    readRDS(path),
    error = function(e) {
      stop("Existing checkpoint is unreadable: ", path, " | ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (!identical(attr(x, "checkpoint_signature"), signature)) {
    stop("Existing checkpoint signature mismatch; refusing silent reuse: ",
         path, call. = FALSE)
  }
  if (!is.list(x) || is.null(x$summary) ||
      !all(as.character(x$summary$status) == "completed")) {
    stop("Existing checkpoint is incomplete/failed; refusing silent reuse: ",
         path, call. = FALSE)
  }
  attr(x, "checkpoint_action") <- "reused_valid"
  x
}

write_checkpoint <- function(x, path, signature) {
  attr(x, "checkpoint_signature") <- signature
  save_rds_atomic(x, path)
  x
}

extract_assignment_name <- function(expr) {
  if (!is.call(expr)) return(NA_character_)
  op <- as.character(expr[[1L]])[1L]
  if (!op %in% c("<-", "=")) return(NA_character_)
  lhs <- expr[[2L]]
  if (!is.symbol(lhs)) return(NA_character_)
  as.character(lhs)
}

load_locked_main_engine <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  env$PROJECT_DIR <- dirname(path)
  env$log_msg <- function(...) invisible(NULL)
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
    "WINSOR_PROBS", "SCENARIOS", "fix_duplicate_colnames", "safe_rbindlist",
    "feature_names", "baseline_vars", "severity_loading", "feature_scale",
    "positive_features", "bounded_low", "bounded_high", "safe_plogis",
    "clamp_prob", "logit", "winsor_fit", "winsor_apply",
    "standardize_train", "standardize_apply", "assert_no_missing_features",
    "make_feature_covariance", "to_clinical_scale", "inject_missingness",
    "simulate_latent_features", "generate_baseline_covariates",
    "generate_treatment", "generate_outcomes", "generate_dataset",
    "make_mice_imputations"
  )
  found <- character()
  for (expr in exprs) {
    nm <- extract_assignment_name(expr)
    if (!is.na(nm) && nm %in% keep) {
      eval(expr, envir = env)
      found <- c(found, nm)
    }
  }
  absent <- setdiff(keep, found)
  if (length(absent)) {
    stop(
      "Locked main engine is missing definitions: ",
      paste(absent, collapse = ", "), call. = FALSE
    )
  }
  if (length(env$feature_names) != 33L) {
    stop("Locked main engine does not define 33 features.", call. = FALSE)
  }
  env
}

# -------------------------------------------------------------------------
# Inputs and task manifest
# -------------------------------------------------------------------------

engine_s1 <- load_locked_main_engine(LOCKED_MAIN_ENGINE)
engine_s7 <- load_locked_main_s7_engine(LOCKED_S7_ENGINE)
if (length(engine_s1$feature_names) != 33L ||
    length(engine_s7$feature_names) != 33L ||
    !identical(engine_s1$feature_names, engine_s7$feature_names)) {
  stop("Locked S1/S7 engine feature definitions disagree.", call. = FALSE)
}
if (!isTRUE(all.equal(as.numeric(engine_s7$S7_CLUSTER_DELTA), S7_DELTA,
                      tolerance = 0))) {
  stop("Locked S7 engine delta does not equal the frozen WP7 value 3.2.",
       call. = FALSE)
}
if (!isTRUE(all.equal(
  as.numeric(engine_s7$S7_CLUSTER_PREVALENCE), c(0.50, 0.50), tolerance = 0
))) {
  stop("Locked S7 engine prevalence does not equal the frozen 0.50/0.50.",
       call. = FALSE)
}

empirical <- as.data.table(readRDS(LOCKED_EMPIRICAL_MATRIX))
if (!"stay_id" %in% names(empirical) || anyDuplicated(empirical$stay_id)) {
  stop("Locked empirical matrix integrity failure.", call. = FALSE)
}
empirical_features <- setdiff(names(empirical), "stay_id")
if (length(empirical_features) != 33L || anyDuplicated(empirical_features)) {
  stop("Locked empirical matrix must contain 33 unique feature columns.",
       call. = FALSE)
}
# The empirical MIMIC feature names need not equal the semi-synthetic engine
# names. This matrix is used only for descriptive marginal dip tests; each
# control-specific copula is fitted to its own S1/S7 processed source matrix.
X_empirical <- as.matrix(empirical[, ..empirical_features])
storage.mode(X_empirical) <- "double"
if (any(!is.finite(X_empirical))) {
  stop("Locked empirical matrix contains non-finite values.", call. = FALSE)
}

empirical_marginal <- marginal_dip_table(
  X_empirical, control = "MIMIC_empirical_descriptive",
  repeat_id = 0L, role = "empirical_processed_margin"
)

tasks <- rbindlist(list(
  CJ(
    control = c("S1_negative_control", "S7_positive_control_delta_3p2"),
    pool = "calibration_pool", repeat_id = seq_len(B_CALIBRATION),
    sorted = TRUE
  ),
  CJ(
    control = c("S1_negative_control", "S7_positive_control_delta_3p2"),
    pool = "evaluation_pool", repeat_id = seq_len(B_EVALUATION),
    sorted = TRUE
  )
))
setorder(tasks, control, pool, repeat_id)
tasks[, task_id := sprintf("%s_%s_rep%03d", control, pool, repeat_id)]
tasks[, seed := as.integer(TASK_BASE_SEED + .I * 10007L)]
tasks[, dgm_family := fifelse(
  control == "S1_negative_control", "S1_continuum", "S7_like_K2"
)]
tasks[, injected_delta := fifelse(
  control == "S1_negative_control", 0, S7_DELTA
)]
tasks[, `:=`(
  source_n = SOURCE_N,
  evaluation_n = EVALUATION_N,
  pit_seed = as.integer(
    PIT_BASE_SEED + repeat_id * 7919L +
      fifelse(control == "S1_negative_control", 0L, 5000000L) +
      fifelse(pool == "calibration_pool", 0L, 10000000L)
  ),
  source_kmeans_seed = fifelse(
    pool == "evaluation_pool", as.integer(seed + 301L), NA_integer_
  ),
  calibration_null_train_seed = fifelse(
    pool == "calibration_pool", as.integer(seed + 401L), NA_integer_
  ),
  calibration_null_evaluation_seed = fifelse(
    pool == "calibration_pool", as.integer(seed + 402L), NA_integer_
  ),
  calibration_null_kmeans_seed = fifelse(
    pool == "calibration_pool", as.integer(seed + 403L), NA_integer_
  ),
  audit_null_train_seed = fifelse(
    pool == "evaluation_pool", as.integer(seed + 501L), NA_integer_
  ),
  audit_null_evaluation_seed = fifelse(
    pool == "evaluation_pool", as.integer(seed + 502L), NA_integer_
  ),
  audit_null_kmeans_seed = fifelse(
    pool == "evaluation_pool", as.integer(seed + 503L), NA_integer_
  ),
  mice_m = MICE_M,
  mice_maxit = MICE_MAXIT,
  primary_imputation_id = PRIMARY_IMPUTATION_ID,
  kmeans_algorithm = "Lloyd",
  kmeans_nstart = KMEANS_NSTART,
  kmeans_itermax = KMEANS_ITERMAX,
  kmeans_fallback_allowed = FALSE,
  run_mode = RUN_MODE,
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA
)]

seed_registry <- rbindlist(lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i]
  common <- data.table(
    module = c("dgm_train", "dgm_evaluation", "mice_train", "mice_evaluation"),
    feature_index = NA_integer_,
    seed = as.integer(c(
      task$seed + 1L, task$seed + 2L,
      task$seed + 101L, task$seed + 202L
    ))
  )
  phase_specific <- if (task$pool == "calibration_pool") {
    data.table(
      module = c(
        "calibration_null_train", "calibration_null_evaluation",
        "calibration_null_kmeans"
      ),
      feature_index = NA_integer_,
      seed = as.integer(c(
        task$calibration_null_train_seed,
        task$calibration_null_evaluation_seed,
        task$calibration_null_kmeans_seed
      ))
    )
  } else {
    data.table(
      module = c(
        "source_kmeans", "audit_null_train", "audit_null_evaluation",
        "audit_null_kmeans"
      ),
      feature_index = NA_integer_,
      seed = as.integer(c(
        task$source_kmeans_seed, task$audit_null_train_seed,
        task$audit_null_evaluation_seed, task$audit_null_kmeans_seed
      ))
    )
  }
  pit <- data.table(
    module = "pit_feature",
    feature_index = seq_len(33L),
    seed = as.integer(task$pit_seed + seq_len(33L) * 1009L)
  )
  registry <- rbindlist(list(common, phase_specific, pit), use.names = TRUE)
  registry[, `:=`(
    task_id = task$task_id,
    control = task$control,
    pool = task$pool,
    repeat_id = task$repeat_id
  )]
  setcolorder(
    registry,
    c("task_id", "control", "pool", "repeat_id", "module",
      "feature_index", "seed")
  )
  registry
}))
expected_seed_registry_rows <-
  2L * B_CALIBRATION * (4L + 3L + 33L) +
  2L * B_EVALUATION * (4L + 4L + 33L)
invalid_seed_n <- seed_registry[
  !is.finite(seed) | seed <= 0L | seed >= .Machine$integer.max, .N
]
duplicate_seed_values <- seed_registry[, .N, by = seed][N > 1L]
cross_module_seed_values <- seed_registry[
  , .(module_n = uniqueN(module)), by = seed
][module_n > 1L]
atomic_fwrite(seed_registry, OUTPUT$seed_registry)
if (invalid_seed_n > 0L || nrow(duplicate_seed_values) > 0L ||
    nrow(cross_module_seed_values) > 0L) {
  atomic_write_lines(c(
    "WP7 GLOBAL SEED QC FAILED",
    "run_status=RUN_FAILED",
    "decision=NOT_EVALUATED",
    paste0("invalid_seed_rows=", invalid_seed_n),
    paste0("duplicated_seed_values=", nrow(duplicate_seed_values)),
    paste0("cross_module_seed_values=", nrow(cross_module_seed_values)),
    "Formal task execution is prohibited."
  ), OUTPUT$partial)
  stop("WP7 global seed registry failed hard QC.", call. = FALSE)
}
atomic_fwrite(tasks, OUTPUT$task_manifest)

run_task <- function(i) {
  task <- tasks[i]
  signature <- checkpoint_signature(task)
  checkpoint <- file.path(TASK_CHECKPOINT_DIR, paste0(task$task_id, ".rds"))
  cached <- read_checkpoint(checkpoint, signature)
  if (!is.null(cached)) return(cached)

  warnings <- character()
  out <- tryCatch(
    withCallingHandlers({
      source_engine <- if (task$control == "S1_negative_control") {
        engine_s1
      } else {
        engine_s7
      }
      source_scenario <- if (task$control == "S1_negative_control") {
        "S1_pure_severity_continuum"
      } else {
        "S7_discrete_structure_outcome_null"
      }
      pair <- simulate_mice_d1_pair(
        engine = source_engine,
        scenario = source_scenario,
        repeat_id = task$repeat_id,
        seed = task$seed,
        n_train = SOURCE_N,
        n_evaluation = EVALUATION_N
      )

      copula <- fit_gaussian_copula(pair$X_train, pit_seed = task$pit_seed)
      common <- data.table(
        copula_nearpd_used = copula$nearpd_used,
        copula_nearpd_frobenius = copula$nearpd_frobenius,
        copula_min_eigen_raw = copula$min_eigen_raw,
        copula_min_eigen_used = copula$min_eigen_used,
        pit_method = copula$pit_method,
        marginal_dip_used_in_go_no_go = FALSE,
        source_pipeline = paste0(
          "formal MAR+MICE m=", MICE_M, "; maxit=", MICE_MAXIT,
          "; primary imputation=", PRIMARY_IMPUTATION_ID,
          "; training preprocessing frozen to evaluation"
        ),
        outcomes_generated_but_excluded = TRUE,
        source_mice_warning_count = pair$mice_warning_count,
        source_mice_warning_text = pair$mice_warning_text,
        task_warning_count = length(warnings),
        task_warning_text = paste(unique(warnings), collapse = " | "),
        status = "completed", error = NA_character_
      )

      if (task$pool == "calibration_pool") {
        calibration_null_train <- simulate_gaussian_copula(
          copula, SOURCE_N, seed = task$calibration_null_train_seed
        )
        calibration_null_evaluation <- simulate_gaussian_copula(
          copula, EVALUATION_N,
          seed = task$calibration_null_evaluation_seed
        )
        calibration_null_diag <- heldout_diagnostics(
          calibration_null_train, calibration_null_evaluation,
          evaluation_true_labels = NULL,
          seed = task$calibration_null_kmeans_seed
        )
        calibration_fidelity <- prefix_columns(
          copula_fidelity(pair$X_train, calibration_null_train),
          "calibration_null_fidelity"
        )
        row <- cbind(
          task,
          prefix_columns(calibration_null_diag, "calibration_null"),
          calibration_fidelity,
          common
        )
        margins <- NULL
      } else {
        true_eval <- if (task$dgm_family == "S7_like_K2") {
          pair$evaluation_dt$true_subtype
        } else {
          NULL
        }
        source_diag <- heldout_diagnostics(
          pair$X_train, pair$X_evaluation,
          evaluation_true_labels = true_eval,
          seed = task$source_kmeans_seed
        )
        audit_null_train <- simulate_gaussian_copula(
          copula, SOURCE_N, seed = task$audit_null_train_seed
        )
        audit_null_evaluation <- simulate_gaussian_copula(
          copula, EVALUATION_N, seed = task$audit_null_evaluation_seed
        )
        audit_null_diag <- heldout_diagnostics(
          audit_null_train, audit_null_evaluation,
          evaluation_true_labels = NULL,
          seed = task$audit_null_kmeans_seed
        )
        audit_fidelity <- prefix_columns(
          copula_fidelity(pair$X_train, audit_null_train),
          "audit_null_fidelity"
        )
        margins <- marginal_dip_table(
          pair$X_train, control = task$control,
          repeat_id = task$repeat_id,
          role = "evaluation_pool_source_training_margin_descriptive"
        )
        row <- cbind(
          task,
          prefix_columns(source_diag, "source"),
          prefix_columns(audit_null_diag, "audit_null"),
          audit_fidelity,
          data.table(
            source_marginal_bh_reject_n = sum(margins$bh_reject_0_05),
            source_marginal_any_bh_reject = any(margins$bh_reject_0_05)
          ),
          common
        )
      }
      list(summary = row, margins = margins)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      list(
        summary = cbind(task, data.table(
          task_warning_count = length(warnings),
          task_warning_text = paste(unique(warnings), collapse = " | "),
          status = "failed", error = conditionMessage(e)
        )),
        margins = NULL
      )
    }
  )
  if (!is.null(out$summary)) {
    out$summary[, `:=`(
      task_warning_count = length(warnings),
      task_warning_text = paste(unique(warnings), collapse = " | ")
    )]
  }
  attr(out, "checkpoint_action") <- "computed_new"
  write_checkpoint(out, checkpoint, signature)
}

# -------------------------------------------------------------------------
# Parallel formal execution
# -------------------------------------------------------------------------

RUN_STAGE <- "task_execution"
ids <- seq_len(nrow(tasks))
if (N_WORKERS > 1L) {
  cl <- parallel::makeCluster(N_WORKERS)
  rows <- tryCatch({
    parallel::clusterCall(cl, function(pkgs) {
      for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      }
      data.table::setDTthreads(1L)
      NULL
    }, c("data.table", "mice", "diptest", "Matrix", "mclust", "digest"))
    parallel::clusterExport(
      cl,
      varlist = c(
        "tasks", "engine_s1", "engine_s7", "TASK_CHECKPOINT_DIR", "RULE_VERSION",
        "CHECKPOINT_SCHEMA", "script_hash_start", "input_hash_start",
        "RUN_MODE", "B_CALIBRATION", "B_EVALUATION",
        "SOURCE_N", "EVALUATION_N", "S7_DELTA", "KMEANS_NSTART",
        "KMEANS_ITERMAX",
        "RIDGE_MULTIPLIER", "DIP_ALPHA", "PIT_BASE_SEED",
        "MICE_M", "MICE_MAXIT", "PRIMARY_IMPUTATION_ID",
        "D1_KMEANS_NSTART", "D1_KMEANS_ITERMAX",
        "D1_MAHALANOBIS_RIDGE_MULTIPLIER",
        "D1_DIP_FDR_ALPHA", "D1_MICE_M", "D1_MICE_MAXIT",
        "D1_PRIMARY_IMPUTATION_ID", "simulate_mice_d1_pair",
        "run_locked_mice_with_warning_audit", "preprocess_train_evaluation",
        "safe_kmeans", "kmeans_audit_record",
        "nearest_centroid_labels", "pooled_within_mahalanobis",
        "randomized_pit", "fit_gaussian_copula", "simulate_gaussian_copula",
        "max_point_mass", "copula_fidelity", "marginal_dip_table",
        "heldout_diagnostics", "prefix_columns", "checkpoint_signature",
        "read_checkpoint", "write_checkpoint", "save_rds_atomic", "run_task"
      ),
      envir = environment()
    )
    parallel::parLapply(cl, ids, run_task)
  }, finally = parallel::stopCluster(cl))
} else {
  rows <- lapply(ids, run_task)
}

results <- rbindlist(lapply(rows, `[[`, "summary"), fill = TRUE)
source_margins <- rbindlist(lapply(rows, `[[`, "margins"), fill = TRUE)
checkpoint_inventory <- rbindlist(lapply(seq_along(rows), function(i) {
  checkpoint_path <- file.path(
    TASK_CHECKPOINT_DIR, paste0(tasks$task_id[i], ".rds")
  )
  data.table(
    task_id = tasks$task_id[i],
    control = tasks$control[i],
    pool = tasks$pool[i],
    repeat_id = tasks$repeat_id[i],
    checkpoint_action = attr(rows[[i]], "checkpoint_action") %||% "unknown",
    checkpoint_path = normalizePath(
      checkpoint_path, winslash = "/", mustWork = TRUE
    ),
    checkpoint_sha256 = sha256_file(checkpoint_path),
    checkpoint_signature = attr(rows[[i]], "checkpoint_signature") %||% NA_character_,
    status = as.character(rows[[i]]$summary$status[1L])
  )
}))
atomic_fwrite(checkpoint_inventory, OUTPUT$checkpoint_inventory)
failures <- results[status != "completed"]
atomic_fwrite(failures, OUTPUT$failures)

if (nrow(results) != nrow(tasks) || nrow(failures)) {
  atomic_write_lines(c(
    "WP7 CONTROL RUN FAILED",
    "run_status=RUN_FAILED",
    "decision=NOT_EVALUATED",
    paste0("run_id=", RUN_ID),
    paste0("expected_tasks=", nrow(tasks)),
    paste0("observed_tasks=", nrow(results)),
    paste0("failed_tasks=", nrow(failures)),
    "No scientific go/no-go decision is valid."
  ), OUTPUT$partial)
  stop("WP7 contains missing or failed tasks; no scientific outputs created.",
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Prespecified control-only thresholds and decisions
# -------------------------------------------------------------------------

calibration <- results[pool == "calibration_pool"]
evaluation <- results[pool == "evaluation_pool"]
calibration_thresholds <- calibration[, .(
  calibration_null_delta_q95 = quantile(
    calibration_null_delta, NULL_QUANTILE, type = 8, names = FALSE
  ),
  calibration_null_delta_median = median(calibration_null_delta)
), by = control]
results <- merge(
  results, calibration_thresholds, by = "control", all.x = TRUE, sort = FALSE
)
results[, `:=`(
  source_shape_positive = NA,
  source_delta_positive = NA,
  audit_null_shape_positive = NA,
  audit_null_delta_positive = NA,
  source_wp7_control_state = NA_character_,
  audit_null_wp7_control_state = NA_character_,
  source_delta_covered_by_calibration_q95 = NA,
  source_delta_minus_calibration_q95 = NA_real_
)]
results[pool == "evaluation_pool", `:=`(
  source_shape_positive = source_holdout_discriminant_dip_p_bh < DIP_ALPHA,
  source_delta_positive = source_delta >= calibration_null_delta_q95,
  audit_null_shape_positive =
    audit_null_holdout_discriminant_dip_p_bh < DIP_ALPHA,
  audit_null_delta_positive = audit_null_delta >= calibration_null_delta_q95,
  source_delta_covered_by_calibration_q95 =
    source_delta < calibration_null_delta_q95,
  source_delta_minus_calibration_q95 =
    source_delta - calibration_null_delta_q95
)]
results[pool == "evaluation_pool", source_wp7_control_state :=
          wp7_control_state(source_shape_positive, source_delta_positive)]
results[pool == "evaluation_pool", audit_null_wp7_control_state :=
          wp7_control_state(audit_null_shape_positive, audit_null_delta_positive)]
results[, verdict_scope := paste(
  "WP7_CONTROL_ONLY components: held-out discriminant BH dip plus a",
  "control-specific q95 frozen from an independent calibration pool;",
  "state names are not the formal v3 empirical verdict"
)]

calibration <- results[pool == "calibration_pool"]
evaluation <- results[pool == "evaluation_pool"]
neg_cal <- calibration[control == "S1_negative_control"]
pos_cal <- calibration[control == "S7_positive_control_delta_3p2"]
neg <- evaluation[control == "S1_negative_control"]
pos <- evaluation[control == "S7_positive_control_delta_3p2"]

neg_delta_fp_n <- sum(neg$source_delta_positive)
neg_delta_fp_ci <- wilson(neg_delta_fp_n, nrow(neg))
neg_dip_fp_n <- sum(neg$source_shape_positive)
neg_dip_fp_ci <- wilson(neg_dip_fp_n, nrow(neg))
neg_joint_fp_n <- sum(
  neg$source_wp7_control_state == "WP7_CONTROL_ONLY_BOTH_POSITIVE"
)
neg_joint_fp_ci <- wilson(neg_joint_fp_n, nrow(neg))
neg_audit_delta_fp_n <- sum(neg$audit_null_delta_positive)
neg_audit_delta_fp_ci <- wilson(neg_audit_delta_fp_n, nrow(neg))
neg_audit_dip_fp_n <- sum(neg$audit_null_shape_positive)
neg_audit_dip_fp_ci <- wilson(neg_audit_dip_fp_n, nrow(neg))
neg_audit_joint_fp_n <- sum(
  neg$audit_null_wp7_control_state == "WP7_CONTROL_ONLY_BOTH_POSITIVE"
)
neg_audit_joint_fp_ci <- wilson(neg_audit_joint_fp_n, nrow(neg))

negative_table <- cbind(
  data.table(
    control = "S1_negative_control",
    truth = "continuous severity; no discrete subtype",
    n_calibration_rep = nrow(neg_cal),
    n_evaluation_rep = nrow(neg),
    control_only_threshold = paste0(
      "type-8 q", NULL_QUANTILE,
      " from an independent S1 calibration pool"
    ),
    calibration_null_delta_q95 = unique(neg$calibration_null_delta_q95),
    source_delta_false_positive_n = neg_delta_fp_n,
    source_delta_false_positive_rate = neg_delta_fp_n / nrow(neg),
    source_delta_false_positive_wilson_low = neg_delta_fp_ci["low"],
    source_delta_false_positive_wilson_high = neg_delta_fp_ci["high"],
    source_dip_false_positive_n = neg_dip_fp_n,
    source_dip_false_positive_rate = neg_dip_fp_n / nrow(neg),
    source_dip_false_positive_wilson_low = neg_dip_fp_ci["low"],
    source_dip_false_positive_wilson_high = neg_dip_fp_ci["high"],
    source_joint_alert_false_positive_n = neg_joint_fp_n,
    source_joint_alert_false_positive_rate = neg_joint_fp_n / nrow(neg),
    source_joint_alert_wilson_low = neg_joint_fp_ci["low"],
    source_joint_alert_wilson_high = neg_joint_fp_ci["high"],
    audit_null_delta_false_positive_n = neg_audit_delta_fp_n,
    audit_null_delta_false_positive_rate = neg_audit_delta_fp_n / nrow(neg),
    audit_null_delta_wilson_low = neg_audit_delta_fp_ci["low"],
    audit_null_delta_wilson_high = neg_audit_delta_fp_ci["high"],
    audit_null_dip_false_positive_n = neg_audit_dip_fp_n,
    audit_null_dip_false_positive_rate = neg_audit_dip_fp_n / nrow(neg),
    audit_null_dip_wilson_low = neg_audit_dip_fp_ci["low"],
    audit_null_dip_wilson_high = neg_audit_dip_fp_ci["high"],
    audit_null_joint_alert_false_positive_n = neg_audit_joint_fp_n,
    audit_null_joint_alert_false_positive_rate =
      neg_audit_joint_fp_n / nrow(neg),
    audit_null_joint_alert_wilson_low = neg_audit_joint_fp_ci["low"],
    audit_null_joint_alert_wilson_high = neg_audit_joint_fp_ci["high"]
  ),
  numeric_summary(neg$source_delta, "source_delta"),
  numeric_summary(neg_cal$calibration_null_delta, "calibration_null_delta"),
  numeric_summary(neg$audit_null_delta, "audit_null_delta"),
  numeric_summary(
    neg$source_delta_minus_calibration_q95,
    "source_minus_calibration_q95"
  )
)
negative_table[, interpretation_boundary := paste(
  "Direct S1 control only; the independent calibration q95 is not the formal",
  "empirical Domain 1 gate and does not assign an empirical Domain 1 state"
)]

pos_covered_n <- sum(pos$source_delta_covered_by_calibration_q95)
pos_covered_ci <- wilson(pos_covered_n, nrow(pos))
pos_detect_n <- sum(pos$source_delta_positive)
pos_detect_ci <- wilson(pos_detect_n, nrow(pos))
pos_joint_n <- sum(
  pos$source_wp7_control_state == "WP7_CONTROL_ONLY_BOTH_POSITIVE"
)
pos_joint_ci <- wilson(pos_joint_n, nrow(pos))
pos_null_q95 <- unique(pos$calibration_null_delta_q95)
pos_source_median <- median(pos$source_delta)
if (length(pos_null_q95) != 1L || length(pos_source_median) != 1L) {
  stop("S7 global threshold summary is not unique.", call. = FALSE)
}

positive_table <- cbind(
  data.table(
    control = "S7_positive_control_delta_3p2",
    truth = "two discrete groups; injected delta=3.2",
    n_calibration_rep = nrow(pos_cal),
    n_evaluation_rep = nrow(pos),
    source_true_label_ari_mean = mean(
      pos$source_holdout_true_label_ari_descriptive, na.rm = TRUE
    ),
    source_true_label_ari_p025 = quantile(
      pos$source_holdout_true_label_ari_descriptive, 0.025,
      type = 8, names = FALSE, na.rm = TRUE
    ),
    source_true_label_ari_p975 = quantile(
      pos$source_holdout_true_label_ari_descriptive, 0.975,
      type = 8, names = FALSE, na.rm = TRUE
    ),
    null_delta_q95 = pos_null_q95,
    source_delta_median_global = pos_source_median,
    null_q95_catches_or_covers_source_median =
      pos_null_q95 >= pos_source_median,
    source_covered_n = pos_covered_n,
    source_covered_rate = pos_covered_n / nrow(pos),
    source_covered_wilson_low = pos_covered_ci["low"],
    source_covered_wilson_high = pos_covered_ci["high"],
    source_delta_detection_n = pos_detect_n,
    source_delta_detection_rate = pos_detect_n / nrow(pos),
    source_delta_detection_wilson_low = pos_detect_ci["low"],
    source_delta_detection_wilson_high = pos_detect_ci["high"],
    source_joint_alert_n_control_only = pos_joint_n,
    source_joint_alert_rate_control_only = pos_joint_n / nrow(pos),
    source_joint_alert_wilson_low = pos_joint_ci["low"],
    source_joint_alert_wilson_high = pos_joint_ci["high"]
  ),
  numeric_summary(pos$source_delta, "source_delta"),
  numeric_summary(pos_cal$calibration_null_delta, "calibration_null_delta"),
  numeric_summary(pos$audit_null_delta, "audit_null_delta"),
  numeric_summary(
    pos$source_delta_minus_calibration_q95,
    "source_minus_calibration_q95"
  )
)
positive_table[, interpretation_boundary := paste(
  "Absolute S7 coverage is reported; the frozen go/no-go uses the Newcombe",
  "interval for S7 Delta detection minus S1 Delta false-positive rate"
)]

make_axis_long <- function(dt, role) {
  prefix <- switch(
    role,
    source = "source",
    calibration_null = "calibration_null",
    audit_null = "audit_null"
  )
  axes <- c("discriminant", paste0("pc", 1:5))
  rows <- rbindlist(lapply(axes, function(axis) {
    pcol <- paste0(prefix, "_holdout_", axis, "_dip_p_bh")
    data.table(
      control = dt$control,
      repeat_id = dt$repeat_id,
      dataset_role = role,
      axis = axis,
      p_bh = dt[[pcol]],
      rejected = dt[[pcol]] < DIP_ALPHA
    )
  }))
  any_col <- paste0(prefix, "_holdout_any_axis_bh_reject")
  rbind(
    rows,
    data.table(
      control = dt$control, repeat_id = dt$repeat_id,
      dataset_role = role, axis = "any_axis_BH",
      p_bh = NA_real_, rejected = dt[[any_col]]
    )
  )
}

axis_long <- rbindlist(list(
  make_axis_long(evaluation, "source"),
  make_axis_long(calibration, "calibration_null"),
  make_axis_long(evaluation, "audit_null")
), use.names = TRUE, fill = TRUE)
multimodality_table <- axis_long[, {
  x <- sum(rejected)
  ci <- wilson(x, .N)
  .(
    n_rep = .N, reject_n = x, reject_rate = x / .N,
    wilson_low = ci["low"], wilson_high = ci["high"],
    mean_bh_p = if (all(is.na(p_bh))) NA_real_ else mean(p_bh, na.rm = TRUE),
    median_bh_p = if (all(is.na(p_bh))) NA_real_ else median(p_bh, na.rm = TRUE)
  )
}, by = .(control, dataset_role, axis)]
multimodality_table[, null_multimodality_status := fifelse(
  dataset_role == "audit_null" & axis == "any_axis_BH" &
    wilson_low > DIP_ALPHA,
  "SYSTEMATIC_ABOVE_NOMINAL",
  fifelse(
    dataset_role == "audit_null" & axis == "any_axis_BH" &
      wilson_high <= DIP_ALPHA,
    "NO_SYSTEMATIC_INCREASE_DETECTED",
    fifelse(
      dataset_role == "audit_null" & axis == "any_axis_BH",
      "INDETERMINATE", "DESCRIPTIVE"
    )
  )
)]
multimodality_table[, go_no_go_role := fifelse(
  dataset_role == "audit_null" & axis == "any_axis_BH",
  "PRIMARY_NULL_MULTIMODALITY_EVENT", "DESCRIPTIVE"
)]

discrimination_difference <- newcombe_difference(
  pos_detect_n, nrow(pos), neg_delta_fp_n, nrow(neg)
)
positive_table[, `:=`(
  s1_delta_false_positive_rate = neg_delta_fp_n / nrow(neg),
  s7_detection_minus_s1_delta_fpr = discrimination_difference["estimate"],
  discrimination_newcombe_low = discrimination_difference["low"],
  discrimination_newcombe_high = discrimination_difference["high"],
  discrimination_status = fifelse(
    discrimination_difference["low"] > 0,
    "RELATIVE_DISCRIMINATION_DEMONSTRATED",
    fifelse(
      discrimination_difference["high"] <= 0,
      "RELATIVE_DISCRIMINATION_FAILED", "INDETERMINATE"
    )
  )
)]

null_shape_status <- multimodality_table[
  dataset_role == "audit_null" & axis == "any_axis_BH",
  setNames(null_multimodality_status, control)
]
null_shape_invalid <- any(null_shape_status == "SYSTEMATIC_ABOVE_NOMINAL")
null_shape_clear <- length(null_shape_status) == 2L && all(
  null_shape_status == "NO_SYSTEMATIC_INCREASE_DETECTED"
)
discrimination_clear <- discrimination_difference["low"] > 0
discrimination_failed <- discrimination_difference["high"] <= 0
s7_absorbed_descriptive <- isTRUE(
  positive_table$null_q95_catches_or_covers_source_median
)

decision <- if (discrimination_failed || null_shape_invalid) {
  "N1_CONTROL_INVALID"
} else if (discrimination_clear && null_shape_clear) {
  "N1_CONTROL_CONTINUE"
} else {
  "N1_CONTROL_INDETERMINATE"
}

# -------------------------------------------------------------------------
# Hard QC before final outputs
# -------------------------------------------------------------------------

pool_counts <- results[, .N, by = .(control, pool)]
pool_counts[, expected_n := fifelse(
  pool == "calibration_pool", B_CALIBRATION, B_EVALUATION
)]
diagnostic_columns <- function(prefix) {
  axes <- c("discriminant", paste0("pc", 1:5))
  c(
    paste0(prefix, "_delta"),
    paste0(prefix, "_holdout_", axes, "_dip_p_raw"),
    paste0(prefix, "_holdout_", axes, "_dip_p_bh")
  )
}
diagnostics_complete <- function(dt, prefix) {
  cols <- diagnostic_columns(prefix)
  if (!all(cols %in% names(dt))) return(FALSE)
  numeric_ok <- all(vapply(cols, function(nm) {
    all(is.finite(dt[[nm]]))
  }, logical(1)))
  any_col <- paste0(prefix, "_holdout_any_axis_bh_reject")
  logical_ok <- any_col %in% names(dt) &&
    all(!is.na(dt[[any_col]])) && all(dt[[any_col]] %in% c(TRUE, FALSE))
  numeric_ok && logical_ok
}

strict_lloyd_audit_complete <- function(dt, prefix) {
  cols <- paste0(prefix, c(
    "_kmeans_algorithm_used", "_kmeans_fallback_used",
    "_kmeans_warning_count", "_kmeans_nstart", "_kmeans_itermax",
    "_kmeans_failed"
  ))
  if (!all(cols %in% names(dt)) || !nrow(dt)) return(FALSE)
  all(dt[[cols[1L]]] == "Lloyd") &&
    all(dt[[cols[2L]]] %in% FALSE) &&
    all(dt[[cols[3L]]] == 0L) &&
    all(dt[[cols[4L]]] == 25L) &&
    all(dt[[cols[5L]]] == 100L) &&
    all(dt[[cols[6L]]] %in% FALSE)
}
seed_collision_n <- tasks[, {
  used <- na.omit(c(
    seed + 1L, seed + 2L, seed + 101L, seed + 202L,
    source_kmeans_seed,
    calibration_null_train_seed, calibration_null_evaluation_seed,
    calibration_null_kmeans_seed,
    audit_null_train_seed, audit_null_evaluation_seed, audit_null_kmeans_seed
  ))
  .(collision = anyDuplicated(used) > 0L)
}, by = task_id][collision == TRUE, .N]
if (!length(seed_collision_n)) seed_collision_n <- 0L

add_qc(
  "task_count",
  nrow(results) == 2L * (B_CALIBRATION + B_EVALUATION),
  nrow(results), 2L * (B_CALIBRATION + B_EVALUATION)
)
add_qc(
  "checkpoint_inventory_complete",
  nrow(checkpoint_inventory) == nrow(tasks) &&
    all(checkpoint_inventory$status == "completed") &&
    all(checkpoint_inventory$checkpoint_action %in%
          c("computed_new", "reused_valid")) &&
    all(nzchar(checkpoint_inventory$checkpoint_sha256)),
  paste0("rows=", nrow(checkpoint_inventory),
         ";completed=", sum(checkpoint_inventory$status == "completed")),
  paste0("rows/completed=", nrow(tasks), "/", nrow(tasks))
)
if (!IS_FORMAL) {
  nonformal_prefix <- paste0(
    normalizePath(NONFORMAL_ROOT, winslash = "/", mustWork = TRUE), "/"
  )
  nonformal_paths <- normalizePath(
    unique(unlist(OUTPUT, use.names = FALSE)),
    winslash = "/", mustWork = FALSE
  )
  add_qc(
    "nonformal_output_path_isolation",
    all(startsWith(nonformal_paths, nonformal_prefix)) &&
      !file.exists(file.path(AUDIT_ROOT, "29f_run_completed_LLOYD.ok")),
    paste0("paths=", length(nonformal_paths),
           ";formal_marker_exists=",
           file.exists(file.path(AUDIT_ROOT, "29f_run_completed_LLOYD.ok"))),
    "all test/smoke outputs under the active nonformal RUN_ID; no formal marker"
  )
}
add_qc(
  "per_control_pool_count",
  all(pool_counts$N == pool_counts$expected_n) && nrow(pool_counts) == 4L,
  paste(pool_counts[, paste0(control, "/", pool, "=", N)], collapse = ";"),
  paste0("two controls x calibration=", B_CALIBRATION,
         " and evaluation=", B_EVALUATION)
)
add_qc("all_tasks_completed", all(results$status == "completed"),
       sum(results$status == "completed"), nrow(results))
add_qc(
  "strict_lloyd_contract_all_cluster_fits",
  strict_lloyd_audit_complete(calibration, "calibration_null") &&
    strict_lloyd_audit_complete(evaluation, "source") &&
    strict_lloyd_audit_complete(evaluation, "audit_null"),
  paste0(
    "calibration=", strict_lloyd_audit_complete(calibration, "calibration_null"),
    ";source=", strict_lloyd_audit_complete(evaluation, "source"),
    ";audit_null=", strict_lloyd_audit_complete(evaluation, "audit_null")
  ),
  "all fits Lloyd; nstart=25; iter.max=100; no warnings/fallback/failure"
)
add_qc("feature_count", ncol(X_empirical) == 33L,
       ncol(X_empirical), 33L)
add_qc("empirical_n", nrow(X_empirical) == 20049L,
       nrow(X_empirical), 20049L)
add_qc("source_evaluation_sizes",
       all(results$source_n == SOURCE_N) && all(results$evaluation_n == EVALUATION_N),
       paste(unique(results$source_n), unique(results$evaluation_n), sep = "/"),
       paste(SOURCE_N, EVALUATION_N, sep = "/"))
add_qc(
  "all_six_axis_diagnostics_finite",
  diagnostics_complete(evaluation, "source") &&
    diagnostics_complete(calibration, "calibration_null") &&
    diagnostics_complete(evaluation, "audit_null"),
  "source/evaluation; calibration-null/calibration; audit-null/evaluation",
  "six raw P, six BH P, delta, and any-axis event complete"
)
add_qc(
  "independent_calibration_thresholds",
  nrow(calibration_thresholds) == 2L &&
    all(is.finite(calibration_thresholds$calibration_null_delta_q95)) &&
    all(evaluation[, .N, by = control]$N == B_EVALUATION),
  paste(signif(calibration_thresholds$calibration_null_delta_q95, 6),
        collapse = ";"),
  "two finite q95 gates from calibration-only repeats"
)
add_qc("randomized_pit_only",
       all(results$pit_method == "randomized_Rueschendorf_distribution_transform"),
       paste(unique(results$pit_method), collapse = ";"),
       "randomized_Rueschendorf_distribution_transform")
add_qc("pit_seeds_unique_and_saved", uniqueN(results$pit_seed) == nrow(results),
       uniqueN(results$pit_seed), nrow(results))
add_qc(
  "global_seed_registry_row_count",
  nrow(seed_registry) == expected_seed_registry_rows,
  nrow(seed_registry), expected_seed_registry_rows,
  "calibration: 4 common + 3 null + 33 PIT; evaluation: 4 common + 4 analysis + 33 PIT"
)
add_qc(
  "global_seed_values_valid",
  invalid_seed_n == 0L,
  invalid_seed_n, 0L,
  "all actual seeds finite, positive, and below the R integer maximum"
)
add_qc(
  "global_seed_values_unique",
  nrow(duplicate_seed_values) == 0L &&
    uniqueN(seed_registry$seed) == nrow(seed_registry),
  nrow(duplicate_seed_values), 0L,
  "zero reused seed values across all tasks, controls, pools, features, and modules"
)
add_qc(
  "cross_module_seed_intersections_zero",
  nrow(cross_module_seed_values) == 0L,
  nrow(cross_module_seed_values), 0L,
  "zero exact seed intersections between named computational modules"
)
add_qc(
  "disjoint_dgm_mice_kmeans_null_seed_offsets",
  seed_collision_n == 0L,
  seed_collision_n, 0L,
  "within-task DGM +1/+2, MICE +101/+202, k-means +301, null +401:503"
)
add_qc("formal_mice_configuration",
       all(results$mice_m == 5L) && all(results$mice_maxit == 5L) &&
         all(results$primary_imputation_id == 1L),
       paste(unique(results$mice_m), unique(results$mice_maxit),
             unique(results$primary_imputation_id), sep = "/"),
       "5/5/1")
add_qc("marginal_dip_not_used_in_go_no_go",
       all(results$marginal_dip_used_in_go_no_go %in% FALSE) &&
         all(empirical_marginal$used_in_go_no_go %in% FALSE),
       "FALSE", "FALSE")
add_qc("source_marginal_rows",
       nrow(source_margins) == 2L * B_EVALUATION * 33L,
       nrow(source_margins), 2L * B_EVALUATION * 33L)
add_qc("valid_control_only_verdicts",
       all(evaluation$source_wp7_control_state %in% c(
         "WP7_CONTROL_ONLY_BOTH_POSITIVE",
         "WP7_CONTROL_ONLY_BOTH_NEGATIVE",
         "WP7_CONTROL_ONLY_DISCORDANT"
       )) && all(evaluation$audit_null_wp7_control_state %in% c(
         "WP7_CONTROL_ONLY_BOTH_POSITIVE",
         "WP7_CONTROL_ONLY_BOTH_NEGATIVE",
         "WP7_CONTROL_ONLY_DISCORDANT"
       )) && all(is.na(calibration$source_wp7_control_state)) &&
         all(is.na(calibration$audit_null_wp7_control_state)),
       "checked", "evaluation-only three WP7 states; calibration states NA")
add_qc("formal_gate_not_claimed",
       all(grepl("WP7_CONTROL_ONLY", results$verdict_scope)),
       "WP7_CONTROL_ONLY", "not WP8 formal gate")
add_qc("s7_delta_fixed", all(pos$injected_delta == 3.2),
       paste(unique(pos$injected_delta), collapse = ","), "3.2")
add_qc("s1_delta_fixed", all(neg$injected_delta == 0),
       paste(unique(neg$injected_delta), collapse = ","), "0")
add_qc("two_audit_null_shape_statuses",
       length(null_shape_status) == 2L,
       length(null_shape_status), 2L)
add_qc("newcombe_interval_finite",
       all(is.finite(discrimination_difference)),
       paste(signif(discrimination_difference, 5), collapse = "/"),
       "estimate/low/high finite")
add_qc("technical_failures_zero", nrow(failures) == 0L,
       nrow(failures), 0L)

qc <- rbindlist(qc_rows)
if (any(!qc$pass)) {
  atomic_fwrite(qc, OUTPUT$qc)
  atomic_write_lines(c(
    "WP7 HARD QC FAILED",
    "run_status=RUN_FAILED",
    "decision=NOT_EVALUATED",
    paste0("run_id=", RUN_ID),
    paste0("failed_checks=", paste(qc[!pass, check], collapse = "; ")),
    "No scientific go/no-go decision is valid."
  ), OUTPUT$partial)
  stop("WP7 hard QC failed; see ", OUTPUT$qc, call. = FALSE)
}

# -------------------------------------------------------------------------
# Publication figure: quantitative grid
# -------------------------------------------------------------------------

CONTROL_LABELS <- c(
  S1_negative_control = "S1 negative control",
  S7_positive_control_delta_3p2 = "S7-like positive control (delta = 3.2)"
)
ROLE_LABELS <- c(
  source = "Source",
  calibration_null = "Calibration null",
  audit_null = "Audit null"
)
COLORS <- c(
  Source = "#2E2E2E", `Calibration null` = "#728995",
  `Audit null` = "#B2BEC4"
)
ACCENT <- "#C77A30"

delta_long <- rbind(
  evaluation[, .(control, repeat_id, role = "Source", delta = source_delta)],
  calibration[, .(
    control, repeat_id, role = "Calibration null",
    delta = calibration_null_delta
  )],
  evaluation[, .(
    control, repeat_id, role = "Audit null", delta = audit_null_delta
  )]
)
delta_long[, control_label := factor(
  CONTROL_LABELS[control], levels = unname(CONTROL_LABELS)
)]
delta_long[, role := factor(role, levels = names(COLORS))]

difference_long <- evaluation[, .(
  control, repeat_id,
  separation_margin = source_delta_minus_calibration_q95
)]
difference_long[, control_label := factor(
  CONTROL_LABELS[control], levels = unname(CONTROL_LABELS)
)]

null_multi_plot <- multimodality_table[
  dataset_role == "audit_null" & axis %in% c(
    "discriminant", paste0("pc", 1:5), "any_axis_BH"
  )
]
null_multi_plot[, axis_label := factor(
  axis,
  levels = c("discriminant", paste0("pc", 1:5), "any_axis_BH"),
  labels = c("Discriminant", paste0("PC", 1:5), "Any axis")
)]
null_multi_plot[, control_label := factor(
  CONTROL_LABELS[control], levels = unname(CONTROL_LABELS)
)]

theme_pub <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "#333333"),
    strip.background = element_rect(fill = "#F1F3F4", colour = "#B8B8B8",
                                    linewidth = 0.35),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", hjust = 0),
    plot.caption = element_text(size = 7.7, colour = "grey35", hjust = 0,
                                lineheight = 1.12),
    plot.tag = element_text(size = 11, face = "bold"),
    panel.spacing = grid::unit(6, "pt"),
    plot.margin = margin(7, 8, 7, 12)
  )

p_a <- ggplot(delta_long, aes(x = role, y = delta, fill = role)) +
  geom_violin(width = 0.82, linewidth = 0.35, colour = "grey35", trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35,
               fill = "white", colour = "#333333") +
  facet_wrap(vars(control_label), nrow = 1) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Source separation versus fitted copula-null separation",
    x = NULL, y = "Mahalanobis centroid distance"
  ) +
  theme_pub + theme(legend.position = "none")

p_b <- ggplot(difference_long, aes(x = control_label, y = separation_margin,
                                   fill = control_label)) +
  geom_hline(yintercept = 0, linetype = "22", colour = "grey70",
             linewidth = 0.45) +
  geom_violin(width = 0.68, trim = TRUE, linewidth = 0.35,
              colour = "grey35") +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35,
               fill = "white", colour = "#333333") +
  scale_fill_manual(values = c("#777777", ACCENT)) +
  labs(
    title = "Source separation above the frozen null gate",
    x = NULL, y = "Delta above q95 gate"
  ) +
  scale_x_discrete(labels = c(
    "S1 negative control" = "S1 negative",
    "S7-like positive control (delta = 3.2)" =
      "S7-like positive\n(delta = 3.2)"
  )) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(7, 16, 7, 12)
  )

p_c <- ggplot(null_multi_plot,
              aes(x = axis_label, y = reject_rate, group = control_label,
                  colour = control_label)) +
  geom_hline(yintercept = DIP_ALPHA, linetype = "22", colour = "grey70",
             linewidth = 0.45) +
  geom_errorbar(aes(ymin = wilson_low, ymax = wilson_high),
                width = 0, linewidth = 0.45, position = position_dodge(0.35)) +
  geom_point(size = 1.9, position = position_dodge(0.35)) +
  scale_colour_manual(
    values = c("#4A4A4A", ACCENT),
    labels = c("S1 negative", "S7-like positive (delta = 3.2)"),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Audit-null multimodality",
    x = NULL, y = "BH rejection (95% Wilson CI)",
    colour = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box.margin = margin(t = -3, r = 0, b = 0, l = 0),
    plot.margin = margin(7, 8, 7, 16)
  )

figure <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.05, 1), guides = "keep") +
  plot_annotation(
    title = "Direct controls for the Gaussian-copula candidate null",
    subtitle = paste0(
      "Formal MAR+MICE (m=5, primary imputation 1); randomized PIT; K=2\n",
      "n=", SOURCE_N, " training + ", EVALUATION_N,
      " independent evaluation observations; calibration/evaluation repeats=",
      B_CALIBRATION, "/", B_EVALUATION, " per control"
    ),
    caption = paste0(
      "Delta is calculated in training data; dip tests use held-out projections only. ",
      "The 5% line is the nominal BH level.\n",
      "Independent calibration pools estimate frozen q95 gates; evaluation-pool audit nulls assess null shape.\n",
      "WP7 control states are direct-control diagnostics, not the WP8 formal empirical gate.\n",
      "Marginal dip tests are descriptive and do not enter the go/no-go decision."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        family = "sans", size = 12, face = "bold", hjust = 0,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        family = "sans", size = 8.5, colour = "grey35", hjust = 0,
        lineheight = 1.08, margin = margin(b = 5)
      ),
      plot.caption = element_text(
        family = "sans", size = 7.7, colour = "grey35", hjust = 0,
        lineheight = 1.12, margin = margin(t = 5)
      ),
      plot.margin = margin(7, 10, 7, 15)
    )
  )

figure_source <- rbindlist(list(
  delta_long[, .(
    panel = "A", control = as.character(control), repeat_id,
    category = as.character(role), estimate = delta,
    ci_low = NA_real_, ci_high = NA_real_
  )],
  difference_long[, .(
    panel = "B", control = as.character(control), repeat_id,
    category = "source_delta_minus_frozen_calibration_q95",
    estimate = separation_margin,
    ci_low = NA_real_, ci_high = NA_real_
  )],
  null_multi_plot[, .(
    panel = "C", control = as.character(control), repeat_id = NA_integer_,
    category = as.character(axis_label), estimate = reject_rate,
    ci_low = wilson_low, ci_high = wilson_high
  )]
), fill = TRUE)

# -------------------------------------------------------------------------
# Final writes, figure export, provenance, completion
# -------------------------------------------------------------------------

RUN_STAGE <- "pre_output_integrity_recheck"
input_hash_end <- rbindlist(lapply(names(INPUT), function(nm) {
  data.table(
    input_name = nm,
    sha256_end = sha256_file(INPUT[[nm]]),
    md5_end = md5_file(INPUT[[nm]])
  )
}))
input_hashes <- merge(input_hash_start, input_hash_end, by = "input_name",
                      all = TRUE)
input_hashes[, unchanged_during_run :=
               sha256 == sha256_end & md5 == md5_end]
if (!all(input_hashes$unchanged_during_run)) {
  stop("An authoritative WP7 input changed during the run.", call. = FALSE)
}
script_hash_end <- sha256_file(SCRIPT_PATH)
if (!identical(script_hash_start, script_hash_end)) {
  stop("The WP7 script changed during execution.", call. = FALSE)
}

RUN_STAGE <- "scientific_output_write"
atomic_fwrite(negative_table, OUTPUT$negative)
atomic_fwrite(positive_table, OUTPUT$positive)
atomic_fwrite(multimodality_table, OUTPUT$multimodality)
atomic_fwrite(results, OUTPUT$by_repeat)
atomic_fwrite(empirical_marginal, OUTPUT$marginal)
atomic_fwrite(source_margins, OUTPUT$source_marginal)
atomic_fwrite(figure_source, OUTPUT$figure_source)
atomic_fwrite(qc, OUTPUT$qc)

width_in <- 183 / 25.4
height_in <- 178 / 25.4
grDevices::cairo_pdf(
  OUTPUT$figure_pdf, width = width_in, height = height_in, family = "sans"
)
print(figure)
grDevices::dev.off()
grDevices::png(
  OUTPUT$figure_png, width = width_in, height = height_in,
  units = "in", res = 600, type = "cairo-png", bg = "white"
)
print(figure)
grDevices::dev.off()
grDevices::tiff(
  OUTPUT$figure_tiff, width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", bg = "white"
)
print(figure)
grDevices::dev.off()

decision_lines <- c(
  paste0("decision=", decision),
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("calibration_repetitions_per_control=", B_CALIBRATION),
  paste0("evaluation_repetitions_per_control=", B_EVALUATION),
  paste0("S1_source_delta_false_positive_rate=", signif(
    negative_table$source_delta_false_positive_rate, 8
  )),
  paste0("S1_source_delta_false_positive_wilson95=",
         signif(negative_table$source_delta_false_positive_wilson_low, 8), ",",
         signif(negative_table$source_delta_false_positive_wilson_high, 8)),
  paste0("S1_source_dip_false_positive_rate=", signif(
    negative_table$source_dip_false_positive_rate, 8
  )),
  paste0("S1_source_joint_alert_false_positive_rate=", signif(
    negative_table$source_joint_alert_false_positive_rate, 8
  )),
  paste0("S7_source_delta_detection_rate=", signif(
    positive_table$source_delta_detection_rate, 8
  )),
  paste0("S7_detection_minus_S1_delta_FPR=",
         signif(discrimination_difference["estimate"], 8)),
  paste0("discrimination_newcombe95=",
         signif(discrimination_difference["low"], 8), ",",
         signif(discrimination_difference["high"], 8)),
  paste0("S7_source_delta_median=", signif(pos_source_median, 8)),
  paste0("S7_calibration_null_delta_q95=", signif(pos_null_q95, 8)),
  paste0("S7_null_q95_catches_or_covers_source_median_descriptive=",
         s7_absorbed_descriptive),
  paste0("audit_null_multimodality_status=",
         paste(null_shape_status, collapse = ";")),
  "marginal_dip_role=DESCRIPTIVE_ONLY_NOT_IN_GO_NO_GO",
  paste0("pit_base_seed=", PIT_BASE_SEED),
  "verdict_scope=WP7_CONTROL_ONLY_NOT_WP8_FORMAL_EMPIRICAL_GATE",
  paste0(
    "interpretation=",
    if (decision == "N1_CONTROL_CONTINUE") {
      paste(
        "Direct controls demonstrated relative discrimination and both",
        "audit-null families ruled out a systematic increase above the",
        "nominal multimodality rate;",
        "this permits later calibration but does not validate N1 as the final",
        "empirical null and does not assign the downstream empirical D1 state."
      )
    } else if (decision == "N1_CONTROL_INDETERMINATE") {
      paste(
        "At least one frozen direct-control component was indeterminate;",
        "N1 cannot proceed as validated or be declared invalid without",
        "reporting that uncertainty, and no threshold revision is allowed."
      )
    } else {
      paste(
        "At least one frozen direct-control failure invalidated N1 as a",
        "candidate null; no parameter tuning or threshold revision is allowed."
      )
    }
  )
)
atomic_write_lines(decision_lines, OUTPUT$decision)
atomic_fwrite(input_hashes, OUTPUT$input_checksums)

run_metadata <- data.table(
  run_id = RUN_ID, run_mode = RUN_MODE,
  reportable_formal_result = IS_FORMAL,
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA,
  calibration_repetitions_per_control = B_CALIBRATION,
  evaluation_repetitions_per_control = B_EVALUATION,
  source_n = SOURCE_N, evaluation_n = EVALUATION_N,
  s7_injected_delta = S7_DELTA,
  mice_m = MICE_M,
  mice_maxit = MICE_MAXIT,
  primary_imputation_id = PRIMARY_IMPUTATION_ID,
  kmeans_algorithm = "Lloyd",
  kmeans_nstart = KMEANS_NSTART,
  kmeans_itermax = KMEANS_ITERMAX,
  kmeans_fallback_allowed = FALSE,
  ridge_multiplier = RIDGE_MULTIPLIER,
  dip_alpha = DIP_ALPHA,
  null_quantile = NULL_QUANTILE,
  pit_method = "randomized_Rueschendorf_distribution_transform",
  pit_base_seed = PIT_BASE_SEED,
  task_base_seed = TASK_BASE_SEED,
  seed_registry_rows = nrow(seed_registry),
  expected_seed_registry_rows = expected_seed_registry_rows,
  invalid_seed_rows = invalid_seed_n,
  duplicated_seed_values = nrow(duplicate_seed_values),
  cross_module_seed_values = nrow(cross_module_seed_values),
  global_seed_registry_pass =
    nrow(seed_registry) == expected_seed_registry_rows &&
    invalid_seed_n == 0L && nrow(duplicate_seed_values) == 0L &&
    nrow(cross_module_seed_values) == 0L,
  workers = N_WORKERS,
  package_startup_warning_class = PACKAGE_STARTUP_WARNING_CLASS,
  package_startup_warning_count = length(PACKAGE_STARTUP_WARNINGS),
  package_startup_warning_text = paste(PACKAGE_STARTUP_WARNINGS, collapse = " | "),
  marginal_dip_used_in_go_no_go = FALSE,
  mutually_independent_calibration_and_evaluation_pools = TRUE,
  wp7_control_only_gate_used = TRUE,
  wp8_formal_empirical_gate_claimed = FALSE,
  reference_scope =
    "processed_space_fixed_k2_empirical_margin_gaussian_copula",
  nonformal_outputs_are_reportable = FALSE,
  go_no_go_decision = decision,
  script_sha256_start = script_hash_start,
  script_sha256_end = script_hash_end,
  prerun_spec_sha256 = NA_character_
)
atomic_fwrite(run_metadata, OUTPUT$run_metadata)

log_lines <- c(
  "# WP7 Gaussian-copula direct-control audit",
  "",
  paste0("- Run ID: `", RUN_ID, "`"),
  paste0("- Run mode: `", RUN_MODE, "`"),
  paste0("- Decision: `", decision, "`"),
  paste0("- Calibration/evaluation repetitions per control: ",
         B_CALIBRATION, "/", B_EVALUATION),
  paste0("- Formal MICE: m=", MICE_M, ", maxit=", MICE_MAXIT,
         ", primary imputation=", PRIMARY_IMPUTATION_ID, "."),
  paste0("- Training/evaluation n: ", SOURCE_N, "/", EVALUATION_N),
  paste0("- K-means: Lloyd; nstart=", KMEANS_NSTART,
         "; iter.max=", KMEANS_ITERMAX, "; fallback=FALSE."),
  paste0("- Workers: ", N_WORKERS),
  paste0("- Package-startup warning class/count: ",
         PACKAGE_STARTUP_WARNING_CLASS, "/", length(PACKAGE_STARTUP_WARNINGS), "."),
  "",
  "## Frozen decision components",
  "",
  paste0(
    "- S1 source Delta false-positive rate: ",
    sprintf("%.4f", negative_table$source_delta_false_positive_rate),
    " (95% Wilson ",
    sprintf("%.4f", negative_table$source_delta_false_positive_wilson_low), " to ",
    sprintf("%.4f", negative_table$source_delta_false_positive_wilson_high), ")."
  ),
  paste0(
    "- S7 source Delta detection rate: ",
    sprintf("%.4f", positive_table$source_delta_detection_rate),
    "; detection minus S1 Delta FPR: ",
    sprintf("%.4f", discrimination_difference["estimate"]),
    " (Newcombe 95% ", sprintf("%.4f", discrimination_difference["low"]),
    " to ", sprintf("%.4f", discrimination_difference["high"]), ")."
  ),
  paste0(
    "- S7 source median delta: ", sprintf("%.4f", pos_source_median),
    "; S7-fitted copula-null q95: ", sprintf("%.4f", pos_null_q95),
    "; descriptively catches/covers the source median: `",
    s7_absorbed_descriptive, "`."
  ),
  paste0(
    "- Audit-null held-out multimodality statuses: `",
    paste(null_shape_status, collapse = "; "), "`."
  ),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "Marginal dip tests are descriptive only. The three WP7 control states use",
    "control-specific q95 thresholds frozen from independent calibration pools and avoid",
    "the formal v3 verdict names. They are not the WP8",
    "formal empirical gate, do not replace locked results, and do not assign",
    "the downstream empirical Domain 1 state."
  ),
  "",
  "## Reproducibility",
  "",
  paste0("- Randomized PIT base seed: `", PIT_BASE_SEED, "`."),
  paste0(
    "- Global actual-seed registry: ", nrow(seed_registry), "/",
    expected_seed_registry_rows,
    " rows; invalid=", invalid_seed_n,
    "; duplicated values=", nrow(duplicate_seed_values),
    "; cross-module intersections=", nrow(cross_module_seed_values), "."
  ),
  "- Every task-level seed and all 33 feature-level PIT seeds are retained in the global seed registry.",
  paste0("- Hard QC: ", nrow(qc), "/", nrow(qc), " checks passed."),
  "- Technical failures: 0.",
  ""
)
atomic_write_lines(log_lines, OUTPUT$log)
atomic_write_lines(capture.output(sessionInfo()), OUTPUT$session)

manifest_targets <- unlist(OUTPUT[c(
  "negative", "positive", "multimodality", "by_repeat", "marginal",
  "source_marginal", "figure_source", "figure_pdf", "figure_png",
  "figure_tiff", "decision", "log", "qc", "failures", "session",
  "task_manifest", "seed_registry", "run_contract", "checkpoint_inventory",
  "input_checksums", "run_metadata"
)], use.names = FALSE)
output_manifest <- rbindlist(lapply(manifest_targets, function(path) {
  if (!file.exists(path)) stop("Missing manifest target: ", path, call. = FALSE)
  data.table(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = sha256_file(path)
  )
}))
atomic_fwrite(output_manifest, OUTPUT$output_manifest)
manifest_sha <- sha256_file(OUTPUT$output_manifest)
fallback_cols <- grep("_kmeans_fallback_used$", names(results), value = TRUE)
cluster_warning_cols <- grep("_kmeans_warning_count$", names(results), value = TRUE)
fallback_n <- sum(vapply(fallback_cols, function(nm) {
  sum(results[[nm]] %in% TRUE, na.rm = TRUE)
}, numeric(1)))
cluster_warning_n <- sum(vapply(cluster_warning_cols, function(nm) {
  sum(results[[nm]], na.rm = TRUE)
}, numeric(1)))
checkpoint_inventory_sha256 <- sha256_file(OUTPUT$checkpoint_inventory)

atomic_write_lines(c(
  if (IS_FORMAL) "RUN COMPLETED" else "NONFORMAL TEST/SMOKE RUN COMPLETED",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("decision=", decision),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", script_hash_end),
  paste0("prerun_spec_sha256=", NA_character_),
  paste0("output_manifest_sha256=", manifest_sha),
  paste0("run_contract_sha256=", run_contract_sha256),
  paste0("checkpoint_inventory_sha256=", checkpoint_inventory_sha256),
  paste0("expected_task_n=", nrow(tasks)),
  paste0("completed_task_n=", sum(results$status == "completed")),
  paste0("failed_task_n=", nrow(failures)),
  paste0("fallback_n=", fallback_n),
  paste0("cluster_warning_n=", cluster_warning_n),
  paste0("package_startup_warning_class=", PACKAGE_STARTUP_WARNING_CLASS),
  paste0("package_startup_warning_n=", length(PACKAGE_STARTUP_WARNINGS)),
  paste0("reportable_formal_result=", IS_FORMAL),
  "verdict_scope=WP7_CONTROL_ONLY_NOT_WP8_FORMAL_EMPIRICAL_GATE"
), OUTPUT$completion)

if (file.exists(OUTPUT$in_progress)) unlink(OUTPUT$in_progress)
RUN_STAGE <- "completed"
options(error = OLD_ERROR_OPTION)

message("WP7 completed successfully: ", decision)
message("Negative control: ", OUTPUT$negative)
message("Positive control: ", OUTPUT$positive)
message("Null multimodality: ", OUTPUT$multimodality)
message("Figure: ", OUTPUT$figure_pdf)
message("Decision: ", OUTPUT$decision)
