public_repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
#!/usr/bin/env Rscript

# WP8: sample-size-specific empirical Gaussian-copula (N1) gates, S7-like K2
# three-state operating surface, null-induced asymptote c0, and empirical
# delta inversion.
# This is an independent complete-data, same-distribution analysis. It does
# not modify WP7 or any locked empirical/simulation output.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(diptest)
  library(Matrix)
  library(mclust)
  library(ggplot2)
  library(patchwork)
  library(digest)
  library(future)
  library(future.apply)
})

required_pkgs <- c(
  "data.table", "diptest", "Matrix", "mclust", "ggplot2", "patchwork", "digest",
  "future", "future.apply", "ps"
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

RUN_MODE <- tolower(Sys.getenv("D1_WP8_RUN_MODE", unset = "test"))
if (!RUN_MODE %in% c("test", "smoke", "formal")) {
  stop("D1_WP8_RUN_MODE must be test, smoke, or formal.", call. = FALSE)
}
IS_FORMAL <- identical(RUN_MODE, "formal")
FORCE_TOP_LEVEL_TEST_FAILURE <- identical(
  Sys.getenv("D1_WP8_TEST_FORCE_TOP_LEVEL_ERROR", unset = "NO"), "YES"
)
if (IS_FORMAL && FORCE_TOP_LEVEL_TEST_FAILURE) {
  stop("The intentional fail-fast test hook is prohibited in formal mode.",
       call. = FALSE)
}
RUN_ID_ENV <- Sys.getenv("D1_WP8_RUN_ID", unset = "")
RUN_ID <- if (nzchar(RUN_ID_ENV)) RUN_ID_ENV else paste0(
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid()
)
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$", RUN_ID)) {
  stop(
    "D1_WP8_RUN_ID must contain only letters, digits, dot, underscore, or ",
    "hyphen (1-80 characters, starting with a letter or digit).",
    call. = FALSE
  )
}

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT)) stop("Set EHR_AUDIT_WORK_ROOT to the authorized workspace.")
AUDIT_ROOT <- file.path(PROJECT_ROOT, "domain1_lloyd_alignment_20260720")
SCRIPT_PATH <- file.path(public_repo, "R/domain1/empirical/03_gate_detection_surface_formal.R")



TEST_EVIDENCE_RUN_ID <- "20260722_test_contract3"
SMOKE_EVIDENCE_RUN_ID <- "20260722_smoke_contract2"
TEST_EVIDENCE_ROOT <- file.path(
  AUDIT_ROOT, "nonformal_runs", "29g", "test", TEST_EVIDENCE_RUN_ID
)
SMOKE_EVIDENCE_ROOT <- file.path(
  AUDIT_ROOT, "nonformal_runs", "29g", "smoke", SMOKE_EVIDENCE_RUN_ID
)












V3_ROOT <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations",
  "Domain1_revision_v3_observable_heldout_20260717"
)
V4_UTILS <- file.path(public_repo, "R/domain1/empirical/00_domain1_revision_utils_v4_lloyd.R")
STRICT_WRAPPER <- file.path(public_repo, "R/domain1/empirical/01_audit_discreteness_v2.R")
V3_LOADER <- file.path(public_repo, "R/domain1/empirical/dependencies/00_domain1_locked_engine_v3.R")


LOCKED_ENGINE <- file.path(public_repo, "R/simulation/core/01_S7_discrete_outcome_null_formal.R")
EMPIRICAL_DELTA_SOURCE <- file.path(
  PROJECT_ROOT, "output", "qc_domain1_mice",
  "posthoc_locked_k2_pooled_within_mahalanobis_20260717.csv"
)
LOCKED_EMPIRICAL_MATRIX <- file.path(
  PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"
)

TABLE_DIR <- file.path(AUDIT_ROOT, "tables")
FIGURE_DIR <- file.path(AUDIT_ROOT, "figures")
LOG_DIR <- file.path(AUDIT_ROOT, "logs")
PROVENANCE_DIR <- file.path(AUDIT_ROOT, "provenance")
CHECKPOINT_ROOT <- file.path(AUDIT_ROOT, "checkpoints", "wp8")
NONFORMAL_ROOT <- file.path(AUDIT_ROOT, "nonformal_runs", "29g")
FORMAL_LOCK_DIR <- file.path(CHECKPOINT_ROOT, "29g_lloyd_formal_exclusive.lock")
FORMAL_LOCK_OWNER_FILE <- file.path(FORMAL_LOCK_DIR, "owner.csv")

if (IS_FORMAL) {
  RUN_ROOT <- file.path(CHECKPOINT_ROOT, "29g_lloyd_formal")
  OUT_TABLE_DIR <- TABLE_DIR
  OUT_FIGURE_DIR <- FIGURE_DIR
  OUT_LOG_DIR <- LOG_DIR
  OUT_PROVENANCE_DIR <- PROVENANCE_DIR
} else {
  RUN_ROOT <- file.path(
    NONFORMAL_ROOT, RUN_MODE, RUN_ID
  )
  OUT_TABLE_DIR <- file.path(RUN_ROOT, "tables")
  OUT_FIGURE_DIR <- file.path(RUN_ROOT, "figures")
  OUT_LOG_DIR <- file.path(RUN_ROOT, "logs")
  OUT_PROVENANCE_DIR <- file.path(RUN_ROOT, "provenance")
}
GATE_CHECKPOINT_DIR <- file.path(RUN_ROOT, "gate_tasks")
SCAN_CHECKPOINT_DIR <- file.path(RUN_ROOT, "scan_tasks")
if (IS_FORMAL && dir.exists(RUN_ROOT)) {
  existing_formal_state <- list.files(
    RUN_ROOT, all.files = TRUE, no.. = TRUE, recursive = TRUE,
    full.names = TRUE
  )
  if (length(existing_formal_state)) {
    stop(
      "Formal 29g run root is not empty. Preserve and audit the existing state; ",
      "do not resume or overwrite it:\n",
      paste(existing_formal_state, collapse = "\n"), call. = FALSE
    )
  }
}
if (!IS_FORMAL && dir.exists(RUN_ROOT)) {
  nonformal_reuse_blockers <- file.path(RUN_ROOT, c(
    "FAILED_DO_NOT_USE.txt",
    "29g_top_level_failure.txt",
    "29g_top_level_failure_LLOYD.txt",
    "run_partial_or_qc_failure.txt",
    "run_partial_or_qc_failure_LLOYD.txt",
    "run_completed.ok",
    "run_completed_LLOYD.ok",
    "29g_nonformal_run_completed.ok",
    "29g_nonformal_run_completed_LLOYD.ok",
    "29g_run_in_progress.txt",
    "29g_run_in_progress_LLOYD.txt"
  ))
  present_reuse_blockers <- nonformal_reuse_blockers[
    file.exists(nonformal_reuse_blockers)
  ]
  if (length(present_reuse_blockers)) {
    stop(
      "Refusing to reuse a failed or completed nonformal WP8 run directory. ",
      "Choose a new D1_WP8_RUN_ID. Blocking marker(s):\n",
      paste(present_reuse_blockers, collapse = "\n"), call. = FALSE
    )
  }
}
for (d in c(
  RUN_ROOT, OUT_TABLE_DIR, OUT_FIGURE_DIR, OUT_LOG_DIR,
  OUT_PROVENANCE_DIR, GATE_CHECKPOINT_DIR, SCAN_CHECKPOINT_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

N_GRID <- if (RUN_MODE == "test") {
  c(2000L, 19049L, 20049L)
} else {
  c(2000L, 5000L, 19049L, 20049L, 50000L)
}
DELTA_GRID <- if (RUN_MODE == "test") {
  c(0, 2.5, 3.2)
} else {
  c(0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.20)
}
PI_GRID <- c(0.20, 0.50)
B_GATE <- if (IS_FORMAL) 200L else 2L
B_SCAN <- if (IS_FORMAL) 100L else 2L
B_BOOT <- if (IS_FORMAL) 2000L else 100L
EVALUATION_N <- 1000L
KMEANS_NSTART <- 25L
KMEANS_ITERMAX <- 100L
RIDGE_MULTIPLIER <- 1e-4
DIP_ALPHA <- 0.05
GATE_QUANTILE <- 0.95
QUANTILE_TYPE <- 8L
EMPIRICAL_EXPECTED_N <- 20049L
EMPIRICAL_EXPECTED_CLUSTER_1_N <- 3992L
EMPIRICAL_EXPECTED_CLUSTER_2_N <- 16057L
EMPIRICAL_EXPECTED_DELTA <- 3.52107344267228
EMPIRICAL_EXPECTED_P <- 33L
EMPIRICAL_EXPECTED_MATRIX_MD5 <- "15a6f073e7e6d6c2d611b0b2d04b8f99"
REQUESTED_PROBABILITY_N <- 19049L
UNIMODALITY_CONTEXT_DELTA <- 2.98
PIT_BASE_SEED <- 314159265L
FUTURE_GATE_SEED <- 930000001L
FUTURE_SCAN_SEED <- 930000002L
RULE_VERSION <- "wp8_empirical_N1_n_specific_gate_three_state_surface_lloyd_v2_20260722"
CHECKPOINT_SCHEMA <- "29g_wp8_lloyd_schema_20260722_1"
TASK_BASE_GATE <- 100000000L
TASK_BASE_SCAN <- 500000000L
Q95_BOOT_SEED_BASE <- 900000000L
C0_BOOT_SEED <- 910000001L
INVERSION_BOOT_SEED <- 920000001L
FORMAL_LOCK_STALE_MINUTES <- 10
FUTURE_SCHEDULING <- if (IS_FORMAL) 4 else 1

detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores)) detected_cores <- 2L
default_workers <- if (IS_FORMAL) 4L else if (RUN_MODE == "smoke") 2L else 1L
N_WORKERS <- as.integer(Sys.getenv(
  "D1_WP8_WORKERS", unset = as.character(default_workers)
))
if (!is.finite(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L
if (IS_FORMAL && N_WORKERS <= 1L) {
  stop("Formal WP8 requires future::multisession with at least 2 workers.",
       call. = FALSE)
}

# Globals required by the sourced v3 diagnostic utilities.
D1_KMEANS_NSTART <- KMEANS_NSTART
D1_KMEANS_ITERMAX <- KMEANS_ITERMAX
D1_MAHALANOBIS_RIDGE_MULTIPLIER <- RIDGE_MULTIPLIER
D1_DIP_FDR_ALPHA <- DIP_ALPHA

for (path in c(
  SCRIPT_PATH, V4_UTILS, STRICT_WRAPPER, V3_LOADER,
  LOCKED_ENGINE, EMPIRICAL_DELTA_SOURCE,
  LOCKED_EMPIRICAL_MATRIX
)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}
active_script_text <- readLines(SCRIPT_PATH, warn = FALSE, encoding = "UTF-8")
if (any(grepl(
  "max\\s*\\([^)]*q95", active_script_text, ignore.case = TRUE
))) {
  stop("Forbidden cross-n q95 aggregation pattern detected in active WP8 script.",
       call. = FALSE)
}
if (any(grepl(
  "source\\s*\\(\\s*V3_OLD_(CALIBRATION|SCAN)_REFERENCE",
  active_script_text
))) {
  stop("Historical buggy gate/scan script must never be sourced by WP8.",
       call. = FALSE)
}

INPUT <- list(
  locked_s7_engine = LOCKED_ENGINE,
  v4_lloyd_utils = V4_UTILS,
  strict_decision_wrapper = STRICT_WRAPPER,
  v3_definition_loader = V3_LOADER,
  empirical_delta_source = EMPIRICAL_DELTA_SOURCE,
  locked_empirical_matrix = LOCKED_EMPIRICAL_MATRIX
)
INPUT_ROLE <- c(
  locked_s7_engine = "computational_engine_definition",
  v4_lloyd_utils = "strict_Lloyd_computational_utility",
  strict_decision_wrapper = "governed_three_state_decision_contract",
  v3_definition_loader = "computational_engine_loader",
  empirical_delta_source = "locked_empirical_inversion_input_filtered_to_compatible_row",
  locked_empirical_matrix = "computational_input_for_fixed_empirical_N1_fit"
)
stopifnot(identical(names(INPUT), names(INPUT_ROLE)))

EXPECTED_FIXED_SHA256 <- c(
  empirical_delta_source =
    "28c33cc58306a2662bb1dbe17edf49d7c5707b897af1eff01599a3d88148e72b",
  locked_empirical_matrix =
    "9896ad4906207ca910e54eb187507037106bff212299586141908f387cfe6279"
)
if (!all(names(EXPECTED_FIXED_SHA256) %in% names(INPUT))) {
  stop("Expected fixed-hash map does not match the active input registry.",
       call. = FALSE)
}

OUTPUT <- list(
  gate = file.path(OUT_TABLE_DIR, "Table_D1_formal_gate_by_n_LLOYD.csv"),
  c0 = file.path(OUT_TABLE_DIR, "Table_D1_c0_formal_LLOYD.csv"),
  surface = file.path(OUT_TABLE_DIR, "Table_D1_three_state_surface_LLOYD.csv"),
  inversion = file.path(
    OUT_TABLE_DIR, "Table_D1_empirical_delta_inversion_LLOYD.csv"
  ),
  requested = file.path(
    OUT_TABLE_DIR, "29g_requested_probabilities_LLOYD.csv"
  ),
  gate_by_repeat = file.path(
    OUT_TABLE_DIR, "29g_formal_gate_by_repeat_LLOYD.csv"
  ),
  scan_by_repeat = file.path(
    OUT_TABLE_DIR, "29g_detection_surface_by_repeat_LLOYD.csv"
  ),
  inversion_mapping = file.path(
    OUT_TABLE_DIR, "29g_empirical_delta_mapping_source_LLOYD.csv"
  ),
  n1_fit_summary = file.path(
    OUT_TABLE_DIR, "29g_n1_fit_summary_LLOYD.csv"
  ),
  feature_space_audit = file.path(
    OUT_TABLE_DIR, "29g_feature_space_role_audit_LLOYD.csv"
  ),
  figure_gate_source = file.path(
    OUT_TABLE_DIR, "29g_Figure_D1_gate_by_n_source_data_LLOYD.csv"
  ),
  figure_surface_source = file.path(
    OUT_TABLE_DIR, "29g_Figure_D1_three_state_surface_source_data_LLOYD.csv"
  ),
  figure_inversion_source = file.path(
    OUT_TABLE_DIR, "29g_Figure_D1_delta_inversion_source_data_LLOYD.csv"
  ),
  gate_pdf = file.path(OUT_FIGURE_DIR, "Figure_D1_gate_by_n_LLOYD.pdf"),
  gate_png = file.path(OUT_FIGURE_DIR, "Figure_D1_gate_by_n_LLOYD.png"),
  gate_tiff = file.path(OUT_FIGURE_DIR, "Figure_D1_gate_by_n_LLOYD.tiff"),
  surface_pdf = file.path(
    OUT_FIGURE_DIR, "Figure_D1_three_state_surface_LLOYD.pdf"
  ),
  surface_png = file.path(
    OUT_FIGURE_DIR, "Figure_D1_three_state_surface_LLOYD.png"
  ),
  surface_tiff = file.path(
    OUT_FIGURE_DIR, "Figure_D1_three_state_surface_LLOYD.tiff"
  ),
  inversion_pdf = file.path(
    OUT_FIGURE_DIR, "Figure_D1_delta_inversion_LLOYD.pdf"
  ),
  inversion_png = file.path(
    OUT_FIGURE_DIR, "Figure_D1_delta_inversion_LLOYD.png"
  ),
  inversion_tiff = file.path(
    OUT_FIGURE_DIR, "Figure_D1_delta_inversion_LLOYD.tiff"
  ),
  log = file.path(OUT_LOG_DIR, "29g_formal_gate_detection_surface_log_LLOYD.md"),
  qc = file.path(OUT_LOG_DIR, "29g_WP8_QC_LLOYD.csv"),
  failures = file.path(OUT_LOG_DIR, "29g_failure_log_LLOYD.csv"),
  session = file.path(OUT_LOG_DIR, "29g_sessionInfo_LLOYD.txt"),
  top_level_failure = file.path(RUN_ROOT, "29g_top_level_failure_LLOYD.txt"),
  partial = file.path(RUN_ROOT, "run_partial_or_qc_failure_LLOYD.txt"),
  in_progress = file.path(RUN_ROOT, "29g_run_in_progress_LLOYD.txt"),
  run_completion = file.path(RUN_ROOT, "run_completed_LLOYD.ok"),
  completion = if (IS_FORMAL) {
    file.path(AUDIT_ROOT, "29g_run_completed_LLOYD.ok")
  } else {
    file.path(RUN_ROOT, "29g_nonformal_run_completed_LLOYD.ok")
  },
  planned_manifest = file.path(
    RUN_ROOT, "task_manifest_planned.csv"
  ),
  task_manifest = file.path(
    OUT_PROVENANCE_DIR, "29g_task_manifest_LLOYD.csv"
  ),
  exact_task_manifest = file.path(RUN_ROOT, "task_manifest.csv"),
  seed_registry = file.path(
    OUT_PROVENANCE_DIR, "29g_global_seed_registry_LLOYD.csv"
  ),
  input_checksums = file.path(
    OUT_PROVENANCE_DIR, "29g_input_checksums_LLOYD.csv"
  ),
  run_metadata = file.path(
    OUT_PROVENANCE_DIR, "29g_run_metadata_LLOYD.csv"
  ),
  n1_fit_rds = file.path(
    OUT_PROVENANCE_DIR, "29g_empirical_gaussian_copula_N1_fit_LLOYD.rds"
  ),
  checkpoint_inventory = file.path(
    OUT_PROVENANCE_DIR, "29g_checkpoint_inventory_LLOYD.csv"
  ),
  output_manifest = file.path(
    OUT_PROVENANCE_DIR, "29g_output_sha256_manifest_LLOYD.csv"
  ),
  runtime_qc = file.path(RUN_ROOT, "29g_WP8_QC_runtime_LLOYD.csv"),
  runtime_failures = file.path(RUN_ROOT, "29g_failure_log_runtime_LLOYD.csv"),
  runtime_seed_registry = file.path(RUN_ROOT, "29g_global_seed_registry_runtime_LLOYD.csv")
)

if (IS_FORMAL) {
  final_targets <- unique(unlist(OUTPUT[c(
    "gate", "c0", "surface", "inversion", "requested",
    "gate_by_repeat", "scan_by_repeat", "inversion_mapping", "n1_fit_summary",
    "feature_space_audit",
    "figure_gate_source", "figure_surface_source",
    "figure_inversion_source", "gate_pdf", "gate_png", "gate_tiff",
    "surface_pdf", "surface_png", "surface_tiff", "inversion_pdf",
    "inversion_png", "inversion_tiff", "log", "qc", "failures",
    "session", "task_manifest", "seed_registry", "input_checksums",
    "run_metadata", "n1_fit_rds", "checkpoint_inventory",
    "output_manifest", "completion"
  )], use.names = FALSE))
  present <- final_targets[file.exists(final_targets)]
  if (length(present)) {
    stop(
      "Formal WP8 final outputs already exist; refusing to overwrite:\n",
      paste(present, collapse = "\n"), call. = FALSE
    )
  }
}

message(
  "WP8 run contract: mode=", RUN_MODE,
  "; n=", paste(N_GRID, collapse = "/"),
  "; delta settings=", length(DELTA_GRID),
  "; pi=", paste(PI_GRID, collapse = "/"),
  "; B gate/scan=", B_GATE, "/", B_SCAN,
  "; n_eval=", EVALUATION_N,
  "; nstart=", KMEANS_NSTART,
  "; future multisession workers=", N_WORKERS
)
message(
  "WP8 palette: no-evidence=#6F7782; inconclusive=#C77A30; ",
  "discrete=#3A5A8C; empirical=#1F3A5F; reference=grey70"
)
message(
  "Figure sizes (mm): gate 183x112; surface 183x150; inversion 183x112"
)

# -------------------------------------------------------------------------
# Integrity and statistical helpers
# -------------------------------------------------------------------------

sha256_file <- function(path) {
  digest::digest(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    algo = "sha256", file = TRUE, serialize = FALSE
  )
}

verify_sha256_manifest <- function(path) {
  manifest <- fread(path)
  required <- c("path", "sha256")
  if (!all(required %in% names(manifest)) || !nrow(manifest)) {
    stop("Invalid SHA-256 manifest: ", path, call. = FALSE)
  }
  manifest[, exists_now := file.exists(path)]
  manifest[, observed_sha256 := vapply(
    path, function(p) if (file.exists(p)) sha256_file(p) else NA_character_,
    character(1)
  )]
  manifest[, hash_match := exists_now & sha256 == observed_sha256]
  manifest
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.table(x)
  if (anyDuplicated(names(x))) {
    stop("Duplicate output columns: ",
         paste(unique(names(x)[duplicated(names(x))]), collapse = ", "),
         call. = FALSE)
  }
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(x, tmp)
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

FORMAL_LOCK_ACQUIRED <- FALSE
FORMAL_LOCK_EVER_ACQUIRED <- FALSE
FORMAL_LOCK_TOKEN <- NA_character_

formal_owner_is_active <- function(owner) {
  if (!nrow(owner) || !all(c(
    "host", "pid", "process_create_time_unix"
  ) %in% names(owner))) {
    return(NA)
  }
  current_host <- unname(Sys.info()[["nodename"]])
  if (!identical(as.character(owner$host[[1L]]), current_host)) {
    return(NA)
  }
  pid <- suppressWarnings(as.integer(owner$pid[[1L]]))
  expected_created <- suppressWarnings(as.numeric(
    owner$process_create_time_unix[[1L]]
  ))
  if (!is.finite(pid) || pid <= 0L || !is.finite(expected_created)) {
    return(NA)
  }
  handle_result <- tryCatch(
    list(handle = ps::ps_handle(pid), error = NULL),
    error = function(e) list(handle = NULL, error = e)
  )
  if (!is.null(handle_result$error)) {
    if (inherits(handle_result$error, "no_such_process")) return(FALSE)
    return(NA)
  }
  handle <- handle_result$handle
  running <- tryCatch(
    ps::ps_is_running(handle),
    error = function(e) {
      if (inherits(e, "no_such_process")) FALSE else NA
    }
  )
  if (is.na(running)) return(NA)
  if (!isTRUE(running)) {
    return(FALSE)
  }
  created_result <- tryCatch(
    list(value = as.numeric(ps::ps_create_time(handle)), error = NULL),
    error = function(e) list(value = NA_real_, error = e)
  )
  if (!is.null(created_result$error)) {
    if (inherits(created_result$error, "no_such_process")) return(FALSE)
    return(NA)
  }
  observed_created <- created_result$value
  if (!is.finite(observed_created)) return(NA)
  abs(observed_created - expected_created) < 1
}

release_formal_lock <- function() {
  if (!IS_FORMAL || !isTRUE(FORMAL_LOCK_ACQUIRED) ||
      !dir.exists(FORMAL_LOCK_DIR)) {
    return(invisible(FALSE))
  }
  owner <- tryCatch(fread(FORMAL_LOCK_OWNER_FILE), error = function(e) NULL)
  token_matches <- is.data.frame(owner) && nrow(owner) == 1L &&
    "token" %in% names(owner) &&
    identical(as.character(owner$token[[1L]]), FORMAL_LOCK_TOKEN)
  if (!token_matches) {
    warning("Formal lock ownership changed; lock was not removed.")
    return(invisible(FALSE))
  }
  expected_parent <- normalizePath(
    CHECKPOINT_ROOT, winslash = "/", mustWork = TRUE
  )
  observed_parent <- normalizePath(
    dirname(FORMAL_LOCK_DIR), winslash = "/", mustWork = TRUE
  )
  if (!identical(expected_parent, observed_parent) ||
      basename(FORMAL_LOCK_DIR) != "29g_formal_exclusive.lock") {
    stop("Refusing to remove an unexpected formal lock path.", call. = FALSE)
  }
  unlink(FORMAL_LOCK_DIR, recursive = TRUE, force = TRUE)
  FORMAL_LOCK_ACQUIRED <<- FALSE
  invisible(!dir.exists(FORMAL_LOCK_DIR))
}

acquire_formal_lock <- function() {
  if (!IS_FORMAL) return(invisible(FALSE))
  dir.create(CHECKPOINT_ROOT, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(FORMAL_LOCK_DIR)) {
    owner <- tryCatch(fread(FORMAL_LOCK_OWNER_FILE), error = function(e) NULL)
    active <- if (is.data.frame(owner)) {
      formal_owner_is_active(owner)
    } else {
      NA
    }
    recover_requested <- identical(
      toupper(Sys.getenv("D1_WP8_RECOVER_STALE_LOCK", unset = "NO")),
      "YES"
    )
    lock_age_minutes <- as.numeric(difftime(
      Sys.time(), file.info(FORMAL_LOCK_DIR)$mtime, units = "mins"
    ))
    same_host <- is.data.frame(owner) && nrow(owner) == 1L &&
      "host" %in% names(owner) &&
      identical(
        as.character(owner$host[[1L]]), unname(Sys.info()[["nodename"]])
      )
    can_recover <- recover_requested && same_host &&
      identical(active, FALSE) && is.finite(lock_age_minutes) &&
      lock_age_minutes >= FORMAL_LOCK_STALE_MINUTES &&
      !file.exists(OUTPUT$completion)
    if (!can_recover) {
      stop(
        "Formal WP8 lock already exists. Active/unknown locks are never ",
        "overwritten. For a confirmed stale same-host lock older than ",
        FORMAL_LOCK_STALE_MINUTES, " minutes, set ",
        "D1_WP8_RECOVER_STALE_LOCK=YES.", call. = FALSE
      )
    }
    stale_path <- paste0(
      FORMAL_LOCK_DIR, ".stale_", format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_pid", Sys.getpid()
    )
    if (!file.rename(FORMAL_LOCK_DIR, stale_path)) {
      stop("Could not archive confirmed stale formal lock.", call. = FALSE)
    }
  }
  if (!dir.create(FORMAL_LOCK_DIR, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not acquire exclusive formal WP8 lock.", call. = FALSE)
  }
  handle <- ps::ps_handle(Sys.getpid())
  process_create_time <- tryCatch(
    as.numeric(ps::ps_create_time(handle)), error = function(e) NA_real_
  )
  if (!is.finite(process_create_time)) {
    unlink(FORMAL_LOCK_DIR, recursive = TRUE, force = TRUE)
    stop(
      "Exclusive formal lock was created but process identity could not be recorded.",
      call. = FALSE
    )
  }
  FORMAL_LOCK_TOKEN <<- digest::digest(
    list(Sys.info()[["nodename"]], Sys.getpid(), RUN_ID, Sys.time()),
    algo = "sha256", serialize = TRUE
  )
  owner <- data.table(
    token = FORMAL_LOCK_TOKEN,
    host = unname(Sys.info()[["nodename"]]),
    pid = Sys.getpid(),
    process_create_time_unix = process_create_time,
    run_id = RUN_ID,
    acquired = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stale_recovery_minimum_minutes = FORMAL_LOCK_STALE_MINUTES
  )
  FORMAL_LOCK_ACQUIRED <<- TRUE
  owner_write_ok <- tryCatch({
    atomic_fwrite(owner, FORMAL_LOCK_OWNER_FILE)
    TRUE
  }, error = function(e) FALSE)
  if (!owner_write_ok) {
    unlink(FORMAL_LOCK_DIR, recursive = TRUE, force = TRUE)
    FORMAL_LOCK_ACQUIRED <<- FALSE
    stop("Exclusive formal lock was created but owner metadata could not be written.",
         call. = FALSE)
  }
  FORMAL_LOCK_EVER_ACQUIRED <<- TRUE
  invisible(TRUE)
}

wp8_checkpoint_signature <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

wp8_read_checkpoint <- function(path, signature) {
  if (!file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.list(obj) && identical(obj$signature, signature) &&
      identical(obj$checkpoint_schema, CHECKPOINT_SCHEMA) &&
      identical(obj$rule_version, RULE_VERSION) &&
      identical(obj$run_id, RUN_ID) &&
      is.data.frame(obj$row)) {
    row <- as.data.table(obj$row)
    if (nrow(row) == 1L && "status" %in% names(row) &&
        identical(as.character(row$status[[1L]]), "completed")) {
      return(row)
    }
  }
  NULL
}

wp8_write_checkpoint <- function(row, path, signature) {
  if (file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    old_sig <- if (is.list(old)) old$signature else NA_character_
    old_status <- if (is.list(old) && is.data.frame(old$row) &&
                      "status" %in% names(old$row)) {
      as.character(old$row$status[[1L]])
    } else {
      NA_character_
    }
    if (!identical(old_sig, signature) ||
        !identical(old_status, "completed")) {
      archive <- paste0(
        path, ".superseded_or_failed_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        "_pid", Sys.getpid()
      )
      if (!file.rename(path, archive)) {
        stop("Could not archive an invalid or superseded checkpoint: ", path,
             call. = FALSE)
      }
    }
  }
  save_rds_atomic(list(
    signature = signature,
    checkpoint_schema = CHECKPOINT_SCHEMA,
    rule_version = RULE_VERSION,
    run_id = RUN_ID,
    row = row
  ), path)
}

wilson <- function(x, n, conf.level = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0 || x < 0 || x > n) {
    return(unname(c(NA_real_, NA_real_)))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  unname(c(max(0, centre - half), min(1, centre + half)))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) {
    as.numeric(quantile(x, p, type = QUANTILE_TYPE, names = FALSE))
  } else {
    NA_real_
  }
}

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
    stop("Empirical N1 correlation matrix is not positive definite.",
         call. = FALSE)
  }
  margins <- lapply(seq_len(p), function(j) sort(X[, j]))
  names(margins) <- colnames(X)
  list(
    R_raw = R_raw,
    R_used = R_used,
    sorted_margins = margins,
    pit_seed = pit_seed,
    min_eigen_raw = min(eig_raw),
    min_eigen_used = min(eig_used),
    nearpd_used = nearpd_used,
    nearpd_frobenius = adjustment,
    pit_method = "randomized_Rueschendorf_distribution_transform",
    source_n = nrow(X),
    source_p = ncol(X)
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

invert_monotone <- function(delta, fitted_mean, target) {
  keep <- is.finite(delta) & is.finite(fitted_mean)
  delta <- delta[keep]
  fitted_mean <- fitted_mean[keep]
  ord <- order(delta)
  delta <- delta[ord]
  fitted_mean <- fitted_mean[ord]
  if (!length(delta) || target < min(fitted_mean) || target > max(fitted_mean)) {
    return(NA_real_)
  }
  inverse_grid <- data.table(
    fitted_mean = fitted_mean,
    delta = delta
  )[, .(
    delta_midpoint = (min(delta) + max(delta)) / 2
  ), by = fitted_mean]
  setorder(inverse_grid, fitted_mean)
  if (nrow(inverse_grid) == 1L) {
    return(if (isTRUE(all.equal(target, inverse_grid$fitted_mean))) {
      inverse_grid$delta_midpoint
    } else {
      NA_real_
    })
  }
  as.numeric(approx(
    x = inverse_grid$fitted_mean,
    y = inverse_grid$delta_midpoint,
    xout = target,
    ties = "ordered", rule = 1
  )$y)
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

script_hash_start <- sha256_file(SCRIPT_PATH)
input_hash_start <- rbindlist(lapply(names(INPUT), function(nm) {
  path <- INPUT[[nm]]
  data.table(
    input_name = nm,
    input_role = unname(INPUT_ROLE[[nm]]),
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(path),
    bytes = file.info(path)$size,
    modified = format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %Z")
  )
}))
input_hash_start[, expected_fixed_sha256 :=
  unname(EXPECTED_FIXED_SHA256[input_name])]
input_hash_start[, fixed_hash_match :=
  is.na(expected_fixed_sha256) | sha256 == expected_fixed_sha256]
fixed_hash_rows <- input_hash_start[
  input_name %in% names(EXPECTED_FIXED_SHA256)
]
if (nrow(fixed_hash_rows) != length(EXPECTED_FIXED_SHA256) ||
    any(!fixed_hash_rows$fixed_hash_match)) {
  mismatch <- fixed_hash_rows[!fixed_hash_match, .(
    input_name, observed = sha256, expected = expected_fixed_sha256
  )]
  stop(
    "At least one fixed authoritative input has the wrong pre-run SHA-256:\n",
    paste(capture.output(print(mismatch)), collapse = "\n"), call. = FALSE
  )
}
input_chain_sha256 <- digest::digest(
  input_hash_start[order(input_name), .(
    input_name, input_role, path, sha256, bytes
  )],
  algo = "sha256", serialize = TRUE
)

OLD_ERROR_OPTION <- getOption("error")
RUN_STAGE <- "initialization"
top_level_failure_handler <- function() {
  failure_lines <- c(
    "WP8 TOP-LEVEL FAILURE",
    paste0("run_id=", RUN_ID),
    paste0("run_mode=", RUN_MODE),
    paste0("stage=", RUN_STAGE),
    paste0("error=", trimws(geterrmessage())),
    paste0("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "No formal WP8 result is lockable."
  )
  if (!IS_FORMAL || isTRUE(FORMAL_LOCK_EVER_ACQUIRED)) {
    for (path in c(OUTPUT$top_level_failure, OUTPUT$partial)) {
      try(atomic_write_lines(failure_lines, path), silent = TRUE)
    }
  }
  try(future::plan(future::sequential), silent = TRUE)
  try(release_formal_lock(), silent = TRUE)
  # A top-level error handler that returns lets Rscript continue evaluating
  # later expressions. Terminate explicitly so a failed run cannot publish
  # tail-end QC, figures, or a false success message.
  options(error = NULL)
  quit(save = "no", status = 1L, runLast = FALSE)
}
options(error = top_level_failure_handler)
acquire_formal_lock()
atomic_write_lines(c(
  "WP8 RUN IN PROGRESS",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", script_hash_start),
  paste0("frozen_spec_sha256=", NA_character_)
), OUTPUT$in_progress)
if (FORCE_TOP_LEVEL_TEST_FAILURE) {
  RUN_STAGE <- "intentional_nonformal_fail_fast_test"
  stop("Intentional nonformal WP8 top-level failure for handler validation.",
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Source the locked observable engine only after hashing inputs
# -------------------------------------------------------------------------

RUN_STAGE <- "source_locked_engine"
source(V4_UTILS, encoding = "UTF-8")
source(STRICT_WRAPPER, encoding = "UTF-8")
source(V3_LOADER, encoding = "UTF-8")
engine <- load_locked_main_s7_engine(LOCKED_ENGINE)
if (!exists("S7_CLUSTER_PREVALENCE", envir = engine, inherits = FALSE) ||
    !exists("S7_CLUSTER_DELTA", envir = engine, inherits = FALSE)) {
  stop("Locked engine lacks mutable S7 prevalence/delta constants.",
       call. = FALSE)
}
engine_feature_names <- as.character(engine$feature_names)
if (length(engine_feature_names) != 33L ||
    anyDuplicated(engine_feature_names) || anyNA(engine_feature_names) ||
    any(!nzchar(engine_feature_names))) {
  stop("Locked engine does not expose exactly 33 features.", call. = FALSE)
}

empirical_all <- fread(EMPIRICAL_DELTA_SOURCE)
required_empirical_fields <- c(
  "compatible_with_full_cohort_matrix", "interpretation_status", "n",
  "p", "cluster_1_n", "cluster_2_n", "pooled_within_mahalanobis",
  "matrix_md5"
)
if (!all(required_empirical_fields %in% names(empirical_all))) {
  stop("Empirical Mahalanobis source lacks compatibility fields.",
       call. = FALSE)
}
empirical <- empirical_all[
  compatible_with_full_cohort_matrix %in% TRUE &
    interpretation_status == "valid locked full-cohort descriptive calculation"
]
if (nrow(empirical) != 1L) {
  stop(
    "Empirical Mahalanobis source must yield exactly one compatible valid row; got ",
    nrow(empirical), call. = FALSE
  )
}
empirical_n <- as.integer(empirical$n)
empirical_p <- as.integer(empirical$p)
empirical_cluster_1_n <- as.integer(empirical$cluster_1_n)
empirical_cluster_2_n <- as.integer(empirical$cluster_2_n)
empirical_pi <- min(empirical_cluster_1_n, empirical_cluster_2_n) / empirical_n
empirical_delta <- as.numeric(empirical$pooled_within_mahalanobis)

if (empirical_n != EMPIRICAL_EXPECTED_N || empirical_p != EMPIRICAL_EXPECTED_P ||
    empirical_cluster_1_n != EMPIRICAL_EXPECTED_CLUSTER_1_N ||
    empirical_cluster_2_n != EMPIRICAL_EXPECTED_CLUSTER_2_N ||
    empirical_cluster_1_n + empirical_cluster_2_n != empirical_n ||
    abs(empirical_delta - EMPIRICAL_EXPECTED_DELTA) > 1e-10) {
  stop("Locked empirical Mahalanobis source differs from the frozen anchor.",
       call. = FALSE)
}
if (!empirical_n %in% N_GRID || !0.20 %in% PI_GRID) {
  stop("Empirical inversion setting is absent from the frozen grid.",
       call. = FALSE)
}
if (!REQUESTED_PROBABILITY_N %in% N_GRID) {
  stop("Requested-probability n is absent from the frozen grid.", call. = FALSE)
}

empirical_matrix_dt <- as.data.table(readRDS(LOCKED_EMPIRICAL_MATRIX))
empirical_matrix_md5_observed <- digest::digest(
  LOCKED_EMPIRICAL_MATRIX, algo = "md5", file = TRUE, serialize = FALSE
)
empirical_matrix_md5_source <- tolower(as.character(empirical$matrix_md5))
if (!identical(
  tolower(empirical_matrix_md5_observed), EMPIRICAL_EXPECTED_MATRIX_MD5
) || !identical(
  empirical_matrix_md5_source, EMPIRICAL_EXPECTED_MATRIX_MD5
)) {
  stop("Locked empirical matrix MD5 identity check failed.", call. = FALSE)
}
if (!"stay_id" %in% names(empirical_matrix_dt)) {
  stop("Locked empirical matrix lacks stay_id.", call. = FALSE)
}
if (nrow(empirical_matrix_dt) != EMPIRICAL_EXPECTED_N ||
    anyDuplicated(empirical_matrix_dt$stay_id)) {
  stop("Locked empirical matrix row integrity failure.", call. = FALSE)
}
empirical_feature_names <- setdiff(names(empirical_matrix_dt), "stay_id")
if (length(empirical_feature_names) != EMPIRICAL_EXPECTED_P ||
    anyDuplicated(empirical_feature_names) || anyNA(empirical_feature_names) ||
    any(!nzchar(empirical_feature_names))) {
  stop("Locked empirical matrix must expose exactly 33 unique named features.",
       call. = FALSE)
}
X_empirical <- as.matrix(empirical_matrix_dt[, ..empirical_feature_names])
storage.mode(X_empirical) <- "double"
if (any(!is.finite(X_empirical))) {
  stop("Locked empirical matrix contains non-finite values.", call. = FALSE)
}

# The empirical N1 gate and the locked S7-like alternative have separate,
# prespecified roles. Exact feature-name identity is not assumed: N1 preserves
# the locked empirical 33-feature geometry, whereas S7 supplies a generator-
# specific 33-dimensional detection surface. Their names and exact overlap are
# recorded so this distinction cannot be hidden or misreported.
feature_name_overlap <- intersect(empirical_feature_names, engine_feature_names)
feature_space_audit <- rbindlist(list(
  data.table(
    feature_space = "empirical_N1_gate",
    feature_index = seq_along(empirical_feature_names),
    feature_name = empirical_feature_names,
    exact_name_in_other_space = empirical_feature_names %in% engine_feature_names,
    scientific_role = "fit fixed empirical Gaussian-copula N1 and calibrate corresponding-n q95"
  ),
  data.table(
    feature_space = "locked_S7_like_alternative",
    feature_index = seq_along(engine_feature_names),
    feature_name = engine_feature_names,
    exact_name_in_other_space = engine_feature_names %in% empirical_feature_names,
    scientific_role = "generate the locked semi-synthetic K2 detection surface and delta inversion scale"
  )
), use.names = TRUE)
feature_space_audit[, `:=`(
  features_in_space = EMPIRICAL_EXPECTED_P,
  exact_name_overlap_n = length(feature_name_overlap),
  exact_feature_name_identity_required = FALSE,
  cross_space_claim = "same dimension and observable diagnostic pipeline; no feature-name identity claim"
)]

RUN_STAGE <- "fit_fixed_empirical_N1"
n1_fit <- fit_gaussian_copula(X_empirical, pit_seed = PIT_BASE_SEED)
n1_fit_summary <- data.table(
  null_id = "N1_empirical_gaussian_copula",
  source_matrix = normalizePath(
    LOCKED_EMPIRICAL_MATRIX, winslash = "/", mustWork = TRUE
  ),
  source_n = n1_fit$source_n,
  source_p = n1_fit$source_p,
  pit_method = n1_fit$pit_method,
  pit_base_seed = PIT_BASE_SEED,
  min_eigen_raw = n1_fit$min_eigen_raw,
  min_eigen_used = n1_fit$min_eigen_used,
  nearpd_used = n1_fit$nearpd_used,
  nearpd_frobenius = n1_fit$nearpd_frobenius,
  preprocessing = paste(
    "no new imputation/winsorization/standardization; null fitted to the",
    "locked processed MICE-imputation-1 standardized matrix"
  ),
  fit_scope = "fixed once before all n-specific calibration replicates",
  empirical_feature_names = paste(empirical_feature_names, collapse = ";"),
  s7_engine_feature_names = paste(engine_feature_names, collapse = ";"),
  exact_feature_name_overlap_n = length(feature_name_overlap),
  exact_feature_name_identity_required = FALSE,
  feature_space_boundary = paste(
    "N1 uses the locked empirical primary33 geometry; S7 uses the locked",
    "semi-synthetic 33-dimensional generator; delta inversion is generator-specific"
  )
)

# -------------------------------------------------------------------------
# Task and seed manifests
# -------------------------------------------------------------------------

RUN_STAGE <- "build_manifests_and_seed_registry"
gate_tasks <- CJ(
  n_training_generated = as.integer(N_GRID),
  repeat_id = seq_len(B_GATE),
  sorted = TRUE
)
gate_tasks[, task_index := .I]
gate_tasks[, task_id := sprintf(
  "gate_N1_n%05d_rep%03d", n_training_generated, repeat_id
)]
gate_tasks[, task_seed := as.integer(TASK_BASE_GATE + task_index * 101L)]
gate_tasks[, `:=`(
  train_generation_seed = task_seed + 1L,
  evaluation_generation_seed = task_seed + 2L,
  primary_kmeans_seed = task_seed + 19L
)]
gate_tasks[, `:=`(
  task_type = "gate",
  dgm_family = "N1_empirical_gaussian_copula",
  delta = 0,
  pi_minor = NA_real_,
  evaluation_n = EVALUATION_N,
  gate_threshold = NA_real_,
  gate_source = "fixed_empirical_N1_corresponding_n_q95_after_gate_stage"
)]

scan_tasks <- CJ(
  n_training_generated = as.integer(N_GRID),
  delta = as.numeric(DELTA_GRID),
  pi_minor = as.numeric(PI_GRID),
  repeat_id = seq_len(B_SCAN),
  sorted = TRUE
)
scan_tasks[, task_index := .I]
scan_tasks[, task_id := sprintf(
  "scan_S7K2_n%05d_delta%04.2f_pi%03.2f_rep%03d",
  n_training_generated, delta, pi_minor, repeat_id
)]
scan_tasks[, task_seed := as.integer(TASK_BASE_SCAN + task_index * 101L)]
scan_tasks[, `:=`(
  train_generation_seed = task_seed + 1L,
  evaluation_generation_seed = task_seed + 2L,
  primary_kmeans_seed = task_seed + 19L
)]
scan_tasks[, `:=`(
  task_type = "scan",
  dgm_family = "S7_like_K2",
  evaluation_n = EVALUATION_N,
  gate_threshold = NA_real_,
  gate_source = "corresponding_n_fixed_empirical_N1_q95"
)]

planned_manifest <- rbindlist(list(gate_tasks, scan_tasks), fill = TRUE)
planned_manifest[, `:=`(
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA,
  run_mode = RUN_MODE
)]
setcolorder(planned_manifest, c(
  "task_id", "task_type", "dgm_family", "n_training_generated",
  "evaluation_n", "delta", "pi_minor", "repeat_id", "task_index",
  "task_seed", "train_generation_seed", "evaluation_generation_seed",
  "primary_kmeans_seed", "gate_source", "gate_threshold",
  "rule_version", "checkpoint_schema", "run_mode"
))
atomic_fwrite(planned_manifest, OUTPUT$planned_manifest)

task_seed_registry <- rbindlist(list(
  planned_manifest[, .(
    task_id, task_type, module = "train_generation",
    seed = train_generation_seed
  )],
  planned_manifest[, .(
    task_id, task_type, module = "evaluation_generation",
    seed = evaluation_generation_seed
  )],
  planned_manifest[, .(
    task_id, task_type, module = "primary_kmeans",
    seed = primary_kmeans_seed
  )]
), fill = TRUE)
postprocess_seed_registry <- rbindlist(list(
  data.table(
    task_id = sprintf("N1_PIT_feature_%02d", seq_along(empirical_feature_names)),
    task_type = "N1_fit",
    module = "randomized_PIT",
    feature = empirical_feature_names,
    seed = as.integer(
      PIT_BASE_SEED + seq_along(empirical_feature_names) * 1009L
    )
  ),
  data.table(
    task_id = sprintf("q95_boot_n%05d", N_GRID),
    task_type = "postprocess",
    module = "q95_bootstrap",
    seed = as.integer(Q95_BOOT_SEED_BASE + seq_along(N_GRID) * 1009L)
  ),
  data.table(
    task_id = c("c0_bootstrap", "empirical_delta_inversion_bootstrap"),
      task_type = "postprocess",
      module = c("c0_bootstrap", "inversion_bootstrap"),
      seed = c(C0_BOOT_SEED, INVERSION_BOOT_SEED)
  ),
  data.table(
    task_id = c("future_gate_scheduler", "future_scan_scheduler"),
    task_type = "scheduler",
    module = c("future_gate_scheduler", "future_scan_scheduler"),
    seed = c(FUTURE_GATE_SEED, FUTURE_SCAN_SEED)
  )
), fill = TRUE)
seed_registry <- rbindlist(
  list(task_seed_registry, postprocess_seed_registry), fill = TRUE
)
seed_registry[, seed := as.integer(seed)]

invalid_seed_n <- seed_registry[
  !is.finite(seed) | seed <= 0L | seed > .Machine$integer.max, .N
]
duplicate_seed_values <- seed_registry[, .N, by = seed][N > 1L]
cross_module_seed_values <- seed_registry[, uniqueN(module), by = seed][V1 > 1L]
expected_seed_rows <- nrow(planned_manifest) * 3L +
  length(empirical_feature_names) +
  length(N_GRID) + 4L
seed_qc_pass <- nrow(seed_registry) == expected_seed_rows &&
  invalid_seed_n == 0L && nrow(duplicate_seed_values) == 0L &&
  nrow(cross_module_seed_values) == 0L
if (!seed_qc_pass) {
  atomic_fwrite(seed_registry, OUTPUT$runtime_seed_registry)
  stop("Global WP8 actual-seed registry failed hard QC.", call. = FALSE)
}

expected_gate_tasks <- length(N_GRID) * B_GATE
expected_scan_tasks <- length(N_GRID) * length(DELTA_GRID) *
  length(PI_GRID) * B_SCAN
if (IS_FORMAL && (expected_gate_tasks != 1000L ||
                  expected_scan_tasks != 11000L ||
                  nrow(planned_manifest) != 12000L)) {
  stop("Formal task-count contract mismatch.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Task functions
# -------------------------------------------------------------------------

analysis_signature <- digest::digest(list(
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA,
  script_sha256 = script_hash_start,
  frozen_spec_sha256 = NA_character_,
  complete_input_chain_sha256 = input_chain_sha256,
  empirical_delta_source_sha256 = sha256_file(EMPIRICAL_DELTA_SOURCE),
  empirical_matrix_sha256 = sha256_file(LOCKED_EMPIRICAL_MATRIX),
  empirical_N1_fit_sha256 = digest::digest(
    n1_fit, algo = "sha256", serialize = TRUE
  ),
  pit_base_seed = PIT_BASE_SEED,
  locked_engine_sha256 = sha256_file(LOCKED_ENGINE),
  v4_lloyd_utils_sha256 = sha256_file(V4_UTILS),
  strict_decision_wrapper_sha256 = sha256_file(STRICT_WRAPPER),
  v3_loader_sha256 = sha256_file(V3_LOADER),
  n_grid = N_GRID,
  delta_grid = DELTA_GRID,
  pi_grid = PI_GRID,
  evaluation_n = EVALUATION_N,
  kmeans_nstart = KMEANS_NSTART,
  ridge_multiplier = RIDGE_MULTIPLIER,
  dip_alpha = DIP_ALPHA
), algo = "sha256", serialize = TRUE)

expected_task_signature <- function(task, task_type) {
  if (identical(task_type, "gate")) {
    task_fields <- list(
      task_id = as.character(task$task_id[[1L]]),
      n_training_generated = as.integer(task$n_training_generated[[1L]]),
      repeat_id = as.integer(task$repeat_id[[1L]]),
      task_seed = as.integer(task$task_seed[[1L]]),
      evaluation_n = as.integer(task$evaluation_n[[1L]])
    )
  } else if (identical(task_type, "scan")) {
    task_fields <- list(
      task_id = as.character(task$task_id[[1L]]),
      n_training_generated = as.integer(task$n_training_generated[[1L]]),
      delta = as.numeric(task$delta[[1L]]),
      pi_minor = as.numeric(task$pi_minor[[1L]]),
      repeat_id = as.integer(task$repeat_id[[1L]]),
      task_seed = as.integer(task$task_seed[[1L]]),
      evaluation_n = as.integer(task$evaluation_n[[1L]]),
      gate_threshold = as.numeric(task$gate_threshold[[1L]])
    )
  } else {
    stop("Unknown task type for checkpoint signature: ", task_type,
         call. = FALSE)
  }
  wp8_checkpoint_signature(list(
    analysis_signature = analysis_signature,
    task = task_fields,
    nstart = KMEANS_NSTART,
    ridge = RIDGE_MULTIPLIER,
    dip_alpha = DIP_ALPHA
  ))
}

checkpoint_path_for <- function(task_id, task_type) {
  dir <- if (identical(task_type, "gate")) {
    GATE_CHECKPOINT_DIR
  } else if (identical(task_type, "scan")) {
    SCAN_CHECKPOINT_DIR
  } else {
    stop("Unknown task type for checkpoint path: ", task_type,
         call. = FALSE)
  }
  file.path(dir, paste0(task_id, ".rds"))
}

assert_strict_lloyd_diagnostic <- function(diag) {
  required <- c(
    "kmeans_primary_algorithm_used", "kmeans_primary_fallback_used",
    "kmeans_primary_warning_count", "kmeans_primary_nstart",
    "kmeans_primary_itermax", "kmeans_primary_failed",
    "kmeans_algorithm", "kmeans_nstart", "kmeans_itermax",
    "shape_source", "reference_scope", "partition_source",
    "oracle_ari_used_in_verdict"
  )
  missing <- setdiff(required, names(diag))
  if (length(missing)) {
    stop("Strict Lloyd diagnostic fields are missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  ok <-
    identical(as.character(diag$kmeans_primary_algorithm_used[[1L]]), "Lloyd") &&
    identical(as.logical(diag$kmeans_primary_fallback_used[[1L]]), FALSE) &&
    identical(as.integer(diag$kmeans_primary_warning_count[[1L]]), 0L) &&
    identical(as.integer(diag$kmeans_primary_nstart[[1L]]), 25L) &&
    identical(as.integer(diag$kmeans_primary_itermax[[1L]]), 100L) &&
    identical(as.logical(diag$kmeans_primary_failed[[1L]]), FALSE) &&
    identical(as.character(diag$kmeans_algorithm[[1L]]), "Lloyd") &&
    identical(as.integer(diag$kmeans_nstart[[1L]]), 25L) &&
    identical(as.integer(diag$kmeans_itermax[[1L]]), 100L) &&
    identical(as.character(diag$shape_source[[1L]]),
              D1_REQUIRED_SHAPE_SOURCE) &&
    identical(as.character(diag$reference_scope[[1L]]),
              D1_REQUIRED_REFERENCE_SCOPE) &&
    identical(as.character(diag$partition_source[[1L]]),
              D1_REQUIRED_PARTITION_SOURCE) &&
    identical(as.logical(diag$oracle_ari_used_in_verdict[[1L]]), FALSE)
  if (!ok) {
    stop("Diagnostic violates the strict Lloyd/held-out contract.",
         call. = FALSE)
  }
  invisible(TRUE)
}

run_gate_task <- function(i) {
  data.table::setDTthreads(1L)
  task <- gate_tasks[i]
  signature <- expected_task_signature(task, "gate")
  checkpoint <- checkpoint_path_for(task$task_id, "gate")
  cached <- wp8_read_checkpoint(checkpoint, signature)
  if (!is.null(cached)) return(cached)

  warning_text <- character()
  row <- tryCatch(withCallingHandlers({
    X_train <- simulate_gaussian_copula(
      n1_fit, n = task$n_training_generated,
      seed = task$train_generation_seed
    )
    X_evaluation <- simulate_gaussian_copula(
      n1_fit, n = EVALUATION_N,
      seed = task$evaluation_generation_seed
    )
    diag <- domain1_external_diagnostics(
      X_train, X_evaluation,
      evaluation_labels = NULL,
      seed = task$task_seed + 17L
    )
    if (!is.finite(diag$mahalanobis_centroid_delta) ||
        !is.finite(diag$dip_discriminant_p_fdr)) {
      stop("N1 gate task produced a non-finite required D1 metric.",
           call. = FALSE)
    }
    assert_strict_lloyd_diagnostic(diag)
    cbind(task, diag, data.table(
      truth_structure = "fixed_empirical_gaussian_copula_candidate_null",
      null_id = "N1_empirical_gaussian_copula",
      analysis_signature = analysis_signature,
      status = "completed",
      error = NA_character_
    ))
  }, warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cbind(task, data.table(
      truth_structure = "fixed_empirical_gaussian_copula_candidate_null",
      null_id = "N1_empirical_gaussian_copula",
      analysis_signature = analysis_signature,
      status = "failed",
      error = conditionMessage(e)
    ))
  })
  row[, `:=`(
    task_warning_count = length(warning_text),
    task_warning_text = paste(unique(warning_text), collapse = " | ")
  )]
  if (identical(as.character(row$status[[1L]]), "completed") &&
      length(warning_text) > 0L) {
    row[, `:=`(
      status = "failed",
      error = paste0(
        "Task emitted warning(s) under the strict formal contract: ",
        paste(unique(warning_text), collapse = " | ")
      )
    )]
  }
  wp8_write_checkpoint(row, checkpoint, signature)
  row
}

run_scan_task <- function(i) {
  data.table::setDTthreads(1L)
  task <- scan_tasks[i]
  if (!is.finite(task$gate_threshold)) {
    stop("Scan task gate was not resolved before execution.", call. = FALSE)
  }
  signature <- expected_task_signature(task, "scan")
  checkpoint <- checkpoint_path_for(task$task_id, "scan")
  cached <- wp8_read_checkpoint(checkpoint, signature)
  if (!is.null(cached)) return(cached)

  warning_text <- character()
  row <- tryCatch(withCallingHandlers({
    old_prevalence <- engine$S7_CLUSTER_PREVALENCE
    old_delta <- engine$S7_CLUSTER_DELTA
    on.exit({
      assign("S7_CLUSTER_PREVALENCE", old_prevalence, envir = engine)
      assign("S7_CLUSTER_DELTA", old_delta, envir = engine)
    }, add = TRUE)
    assign(
      "S7_CLUSTER_PREVALENCE",
      c(task$pi_minor, 1 - task$pi_minor),
      envir = engine
    )
    assign("S7_CLUSTER_DELTA", task$delta, envir = engine)
    pair <- simulate_complete_d1_pair(
      engine = engine,
      dgm_family = "S7_like_K2",
      n_train = task$n_training_generated,
      n_evaluation = EVALUATION_N,
      delta = task$delta,
      seed = task$task_seed
    )
    diag <- domain1_external_diagnostics(
      pair$X_train, pair$X_evaluation,
      evaluation_labels = pair$evaluation_dt$true_subtype,
      seed = task$task_seed + 17L
    )
    if (!is.finite(diag$mahalanobis_centroid_delta) ||
        !is.finite(diag$dip_discriminant_p_fdr)) {
      stop("S7 scan task produced a non-finite required D1 metric.",
           call. = FALSE)
    }
    assert_strict_lloyd_diagnostic(diag)
    governed <- audit_discreteness_v2(list(
      separation = as.numeric(diag$mahalanobis_centroid_delta[[1L]]),
      gate = as.numeric(task$gate_threshold[[1L]]),
      shape_reject = as.logical(
        diag$dip_discriminant_p_fdr[[1L]] < DIP_ALPHA
      ),
      shape_source = as.character(diag$shape_source[[1L]]),
      reference_scope = as.character(diag$reference_scope[[1L]]),
      partition_source = as.character(diag$partition_source[[1L]]),
      kmeans_algorithm = as.character(diag$kmeans_algorithm[[1L]]),
      kmeans_nstart = as.integer(diag$kmeans_nstart[[1L]]),
      kmeans_itermax = as.integer(diag$kmeans_itermax[[1L]])
    ))
    if (identical(governed$state, "NOT_EVALUATED")) {
      stop("Governed Domain 1 wrapper rejected the diagnostic: ",
           governed$failure_code, ": ", governed$failure_detail,
           call. = FALSE)
    }
    verdict <- unname(c(
      NO_DISCRETE_EVIDENCE = "no_discrete_evidence",
      INCONCLUSIVE = "inconclusive",
      DISCRETE_EVIDENCE = "discrete_evidence"
    )[[governed$state]])
    train_props <- prop.table(table(pair$train_dt$true_subtype))
    eval_props <- prop.table(table(pair$evaluation_dt$true_subtype))
    cbind(task, diag, data.table(
      truth_structure = if (task$delta > 0) {
        "two_component_discrete_generator"
      } else {
        "zero_shift_latent_labels_no_observed_separation"
      },
      realized_train_minor_prevalence = min(as.numeric(train_props)),
      realized_evaluation_minor_prevalence = min(as.numeric(eval_props)),
      shape_signal = diag$dip_discriminant_p_fdr < DIP_ALPHA,
      separation_signal =
        diag$mahalanobis_centroid_delta >= task$gate_threshold,
      revised_verdict = verdict,
      decision_rule_version = governed$decision_rule_version,
      analysis_signature = analysis_signature,
      status = "completed",
      error = NA_character_
    ))
  }, warning = function(w) {
    warning_text <<- c(warning_text, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cbind(task, data.table(
      truth_structure = if (task$delta > 0) {
        "two_component_discrete_generator"
      } else {
        "zero_shift_latent_labels_no_observed_separation"
      },
      oracle_ari_used_in_verdict = FALSE,
      analysis_signature = analysis_signature,
      status = "failed",
      error = conditionMessage(e)
    ))
  })
  row[, `:=`(
    task_warning_count = length(warning_text),
    task_warning_text = paste(unique(warning_text), collapse = " | ")
  )]
  if (identical(as.character(row$status[[1L]]), "completed") &&
      length(warning_text) > 0L) {
    row[, `:=`(
      status = "failed",
      error = paste0(
        "Task emitted warning(s) under the strict formal contract: ",
        paste(unique(warning_text), collapse = " | ")
      )
    )]
  }
  wp8_write_checkpoint(row, checkpoint, signature)
  row
}

run_multisession <- function(ids, fun, label, scheduler_seed) {
  message(
    "Starting ", label, ": ", length(ids),
    " tasks with future::multisession workers=", N_WORKERS
  )
  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 2 * 1024^3)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)
  future::plan(future::multisession, workers = N_WORKERS)
  on.exit(future::plan(future::sequential), add = TRUE)
  future.apply::future_lapply(
    ids, fun,
    future.seed = scheduler_seed,
    future.scheduling = FUTURE_SCHEDULING,
    future.packages = c("data.table", "diptest", "Matrix", "mclust")
  )
}

# -------------------------------------------------------------------------
# Gate calibration
# -------------------------------------------------------------------------

RUN_STAGE <- "gate_multisession"
gate_rows <- run_multisession(
  seq_len(nrow(gate_tasks)), run_gate_task,
  "sample-size-specific empirical N1 gate", FUTURE_GATE_SEED
)
future::plan(future::sequential)
gate_results <- rbindlist(gate_rows, fill = TRUE)
if (nrow(gate_results) != nrow(gate_tasks) ||
    uniqueN(gate_results$task_id) != nrow(gate_tasks)) {
  stop("Gate task count or task-ID uniqueness mismatch.", call. = FALSE)
}

accepted_task_status <- "completed"
gate_failures <- gate_results[!status %in% accepted_task_status]
if (nrow(gate_failures)) {
  atomic_fwrite(gate_failures, OUTPUT$runtime_failures)
  stop("At least one gate task failed; formal gate cannot be released.",
       call. = FALSE)
}
if (any(!is.finite(gate_results$mahalanobis_centroid_delta))) {
  stop("Gate results contain non-finite separation values.", call. = FALSE)
}

RUN_STAGE <- "gate_summary_and_c0"
gate_summary_rows <- lapply(seq_along(N_GRID), function(j) {
  n_value <- N_GRID[j]
  gate_n <- gate_results[n_training_generated == n_value]
  x <- gate_n$mahalanobis_centroid_delta
  seed <- as.integer(Q95_BOOT_SEED_BASE + j * 1009L)
  set.seed(seed)
  boot_q95 <- replicate(B_BOOT, safe_quantile(
    sample(x, length(x), replace = TRUE), GATE_QUANTILE
  ))
  data.table(
    n_training_generated = n_value,
    n_planned = B_GATE,
    n_completed = length(x),
    n_failed = B_GATE - length(x),
    separation_mean = mean(x),
    separation_sd = sd(x),
    separation_mcse = sd(x) / sqrt(length(x)),
    separation_p025 = safe_quantile(x, 0.025),
    separation_p975 = safe_quantile(x, 0.975),
    q95 = safe_quantile(x, GATE_QUANTILE),
    q95_bootstrap_B = B_BOOT,
    q95_bootstrap_seed = seed,
    q95_bootstrap_se = sd(boot_q95),
    q95_bootstrap_ci_low = safe_quantile(boot_q95, 0.025),
    q95_bootstrap_ci_high = safe_quantile(boot_q95, 0.975),
    quantile_type = QUANTILE_TYPE,
    dip_rejection_count = sum(gate_n$dip_discriminant_p_fdr < DIP_ALPHA),
    dip_rejection_rate = mean(gate_n$dip_discriminant_p_fdr < DIP_ALPHA),
    fallback_count = sum(gate_n$kmeans_primary_fallback_used %in% TRUE),
    warning_task_count = sum(
      gate_n$task_warning_count > 0L |
        gate_n$kmeans_primary_warning_count > 0L
    ),
    warning_total_count = sum(
      gate_n$task_warning_count + gate_n$kmeans_primary_warning_count
    ),
    gate_role = "sample_size_specific_fixed_empirical_N1_q95",
    gate_null = "fixed_empirical_gaussian_copula_N1"
  )
})
gate_table <- rbindlist(gate_summary_rows)
if (nrow(gate_table) != length(N_GRID) ||
    any(gate_table$n_completed != B_GATE) || any(gate_table$n_failed != 0L)) {
  stop("Gate summary does not contain complete B_gate at every n.",
       call. = FALSE)
}

gate_means <- gate_table[, .(
  n_training_generated,
  inv_n = 1 / n_training_generated,
  inv_sqrt_n = 1 / sqrt(n_training_generated),
  mean_delta = separation_mean,
  mean_mcse = separation_mcse
)]
if (any(!is.finite(gate_means$mean_mcse) | gate_means$mean_mcse <= 0)) {
  stop("c0 primary inverse-MC-variance weights are non-finite or non-positive.",
       call. = FALSE)
}
gate_means[, primary_weight := 1 / mean_mcse^2]
c0_fit_primary <- lm(
  mean_delta ~ inv_n, data = gate_means, weights = primary_weight
)
c0_fit_sensitivity <- lm(
  mean_delta ~ inv_sqrt_n, data = gate_means
)
c0_estimate <- unname(coef(c0_fit_primary)[1L])
c0_slope <- unname(coef(c0_fit_primary)[2L])
c0_sensitivity_estimate <- unname(coef(c0_fit_sensitivity)[1L])
c0_sensitivity_slope <- unname(coef(c0_fit_sensitivity)[2L])

gate_split <- split(
  gate_results$mahalanobis_centroid_delta,
  gate_results$n_training_generated
)
set.seed(C0_BOOT_SEED)
c0_boot <- replicate(B_BOOT, {
  boot_means <- vapply(as.character(N_GRID), function(nm) {
    x <- gate_split[[nm]]
    mean(sample(x, length(x), replace = TRUE))
  }, numeric(1))
  primary <- lm(
    boot_means ~ gate_means$inv_n,
    weights = gate_means$primary_weight
  )
  sensitivity <- lm(boot_means ~ gate_means$inv_sqrt_n)
  c(
    primary_c0 = unname(coef(primary)[1L]),
    sensitivity_c0 = unname(coef(sensitivity)[1L])
  )
})
c0_boot <- t(c0_boot)
primary_c0_boot <- c0_boot[, "primary_c0"]
sensitivity_c0_boot <- c0_boot[, "sensitivity_c0"]
c0_table <- rbindlist(list(
  data.table(
    model_id = "primary_inverse_n_inverse_mc_variance",
    is_primary = TRUE,
    model = "E(Delta_hat|n)=c0+c1/n",
    fit_unit = paste0(
      length(N_GRID), "_setting_level_means_inverse_MC_variance_weighted"
    ),
    weighting = "inverse_original_setting_level_MCSE_squared",
    c0 = c0_estimate,
    c1 = c0_slope,
    model_r_squared = summary(c0_fit_primary)$r.squared,
    c0_bootstrap_se = safe_sd(primary_c0_boot),
    c0_bootstrap_ci_low = safe_quantile(primary_c0_boot, 0.025),
    c0_bootstrap_ci_high = safe_quantile(primary_c0_boot, 0.975)
  ),
  data.table(
    model_id = "sensitivity_inverse_sqrt_n_equal_weight",
    is_primary = FALSE,
    model = "E(Delta_hat|n)=c0+c1/sqrt(n)",
    fit_unit = paste0(length(N_GRID), "_setting_level_means_equal_weight"),
    weighting = "equal_setting_weight",
    c0 = c0_sensitivity_estimate,
    c1 = c0_sensitivity_slope,
    model_r_squared = summary(c0_fit_sensitivity)$r.squared,
    c0_bootstrap_se = safe_sd(sensitivity_c0_boot),
    c0_bootstrap_ci_low = safe_quantile(sensitivity_c0_boot, 0.025),
    c0_bootstrap_ci_high = safe_quantile(sensitivity_c0_boot, 0.975)
  )
), fill = TRUE)
c0_table[, `:=`(
  n_settings = length(N_GRID),
  B_gate_per_n = B_GATE,
  c0_bootstrap_B = B_BOOT,
  c0_bootstrap_seed = C0_BOOT_SEED,
  interpretation = paste(
    "N1_pipeline_specific_asymptotic_partition_separation;",
    "not_a_gate_and_not_a_universal_constant"
  )
)]
gate_table[, `:=`(
  c0_primary = c0_estimate,
  c0_primary_bootstrap_ci_low = c0_table[
    is_primary == TRUE, c0_bootstrap_ci_low
  ],
  c0_primary_bootstrap_ci_high = c0_table[
    is_primary == TRUE, c0_bootstrap_ci_high
  ],
  c0_sensitivity = c0_sensitivity_estimate
)]

# Resolve each scan task to its corresponding-n gate.
scan_tasks <- merge(
  scan_tasks,
  gate_table[, .(
    n_training_generated,
    gate_threshold_resolved = q95
  )],
  by = "n_training_generated", all.x = TRUE, sort = FALSE
)
scan_tasks[, gate_threshold := gate_threshold_resolved]
scan_tasks[, gate_threshold_resolved := NULL]
setorder(scan_tasks, task_index)
if (any(!is.finite(scan_tasks$gate_threshold))) {
  stop("At least one scan task lacks a corresponding-n q95.", call. = FALSE)
}

resolved_manifest <- copy(planned_manifest)
resolved_manifest[task_type == "scan", gate_threshold :=
  scan_tasks$gate_threshold[match(task_id, scan_tasks$task_id)]]
resolved_manifest[task_type == "gate", gate_threshold :=
  gate_table$q95[match(n_training_generated, gate_table$n_training_generated)]]
resolved_manifest[, manifest_stage := "resolved_before_scan"]
atomic_fwrite(resolved_manifest, OUTPUT$exact_task_manifest)

# -------------------------------------------------------------------------
# Detection surface
# -------------------------------------------------------------------------

RUN_STAGE <- "scan_multisession"
scan_rows <- run_multisession(
  seq_len(nrow(scan_tasks)), run_scan_task,
  "S7-like K2 detection surface", FUTURE_SCAN_SEED
)
future::plan(future::sequential)
scan_results <- rbindlist(scan_rows, fill = TRUE)
if (nrow(scan_results) != nrow(scan_tasks) ||
    uniqueN(scan_results$task_id) != nrow(scan_tasks)) {
  stop("Scan task count or task-ID uniqueness mismatch.", call. = FALSE)
}

scan_failures <- scan_results[!status %in% accepted_task_status]
failure_table <- rbindlist(list(gate_failures, scan_failures), fill = TRUE)
if (!nrow(failure_table)) {
  failure_table <- data.table(
    task_id = character(), task_type = character(), status = character(),
    error = character()
  )
}
if (nrow(scan_failures)) {
  atomic_fwrite(failure_table, OUTPUT$runtime_failures)
  stop("At least one scan task failed; formal surface is partial.",
       call. = FALSE)
}

valid_verdicts <- c(
  "no_discrete_evidence", "inconclusive", "discrete_evidence"
)
if (any(!scan_results$revised_verdict %in% valid_verdicts) ||
    any(!is.finite(scan_results$mahalanobis_centroid_delta)) ||
    any(!is.finite(scan_results$dip_discriminant_p_fdr))) {
  stop("Completed scan results contain invalid verdicts or metrics.",
       call. = FALSE)
}

surface_rows <- scan_results[, {
  n_completed <- .N
  counts <- table(factor(revised_verdict, levels = valid_verdicts))
  make_state <- function(state) {
    x <- as.integer(counts[[state]])
    rate <- x / n_completed
    ci <- wilson(x, n_completed)
    c(
      count = x,
      rate = rate,
      mcse = sqrt(rate * (1 - rate) / n_completed),
      ci_low = ci[[1L]],
      ci_high = ci[[2L]]
    )
  }
  no_ev <- make_state("no_discrete_evidence")
  inc <- make_state("inconclusive")
  disc <- make_state("discrete_evidence")
  list(
    n_planned = B_SCAN,
    n_completed = n_completed,
    n_failed = B_SCAN - n_completed,
    gate_threshold = unique(gate_threshold),
    separation_mean = mean(mahalanobis_centroid_delta),
    separation_sd = sd(mahalanobis_centroid_delta),
    separation_mcse = sd(mahalanobis_centroid_delta) / sqrt(n_completed),
    separation_p025 = safe_quantile(mahalanobis_centroid_delta, 0.025),
    separation_p975 = safe_quantile(mahalanobis_centroid_delta, 0.975),
    dip_rejection_count = sum(dip_discriminant_p_fdr < DIP_ALPHA),
    dip_rejection_rate = mean(dip_discriminant_p_fdr < DIP_ALPHA),
    fallback_count = sum(kmeans_primary_fallback_used %in% TRUE),
    warning_task_count = sum(
      task_warning_count > 0L | kmeans_primary_warning_count > 0L
    ),
    warning_total_count = sum(
      task_warning_count + kmeans_primary_warning_count
    ),
    no_evidence_count = no_ev["count"],
    no_evidence_rate = no_ev["rate"],
    no_evidence_mcse = no_ev["mcse"],
    no_evidence_wilson_low = no_ev["ci_low"],
    no_evidence_wilson_high = no_ev["ci_high"],
    inconclusive_count = inc["count"],
    inconclusive_rate = inc["rate"],
    inconclusive_mcse = inc["mcse"],
    inconclusive_wilson_low = inc["ci_low"],
    inconclusive_wilson_high = inc["ci_high"],
    discrete_count = disc["count"],
    discrete_rate = disc["rate"],
    discrete_mcse = disc["mcse"],
    discrete_wilson_low = disc["ci_low"],
    discrete_wilson_high = disc["ci_high"],
    realized_train_minor_prevalence_mean =
      mean(realized_train_minor_prevalence),
    realized_evaluation_minor_prevalence_mean =
      mean(realized_evaluation_minor_prevalence)
  )
}, by = .(n_training_generated, delta, pi_minor)]
setorder(surface_rows, pi_minor, n_training_generated, delta)

if (nrow(surface_rows) != length(N_GRID) * length(DELTA_GRID) *
    length(PI_GRID) || any(surface_rows$n_completed != B_SCAN) ||
    any(surface_rows$n_failed != 0L)) {
  stop("Three-state surface setting completeness failed.", call. = FALSE)
}
if (any(abs(
  surface_rows$no_evidence_rate + surface_rows$inconclusive_rate +
    surface_rows$discrete_rate - 1
) > 1e-12)) {
  stop("Three-state probabilities do not sum to one.", call. = FALSE)
}

requested_table <- rbindlist(list(
  surface_rows[
    n_training_generated == REQUESTED_PROBABILITY_N &
      pi_minor == 0.20 & delta == 2.5,
    .(
      request = paste0(
        "P(no_discrete_evidence|delta=2.5,pi=0.2,n=",
        REQUESTED_PROBABILITY_N, ")"
      ),
      n_training_generated, delta, pi_minor,
      state = "no_discrete_evidence", n_planned, n_completed, n_failed,
      count = no_evidence_count, probability = no_evidence_rate,
      mcse = no_evidence_mcse,
      wilson_low = no_evidence_wilson_low,
      wilson_high = no_evidence_wilson_high
    )
  ],
  surface_rows[
    n_training_generated == REQUESTED_PROBABILITY_N &
      pi_minor == 0.20 & delta == 2.5,
    .(
      request = paste0(
        "P(inconclusive|delta=2.5,pi=0.2,n=",
        REQUESTED_PROBABILITY_N, ")"
      ),
      n_training_generated, delta, pi_minor,
      state = "inconclusive", n_planned, n_completed, n_failed,
      count = inconclusive_count, probability = inconclusive_rate,
      mcse = inconclusive_mcse,
      wilson_low = inconclusive_wilson_low,
      wilson_high = inconclusive_wilson_high
    )
  ],
  surface_rows[
    n_training_generated == REQUESTED_PROBABILITY_N &
      pi_minor == 0.20 & delta == 3.2,
    .(
      request = paste0(
        "P(discrete_evidence|delta=3.2,pi=0.2,n=",
        REQUESTED_PROBABILITY_N, ")"
      ),
      n_training_generated, delta, pi_minor,
      state = "discrete_evidence", n_planned, n_completed, n_failed,
      count = discrete_count, probability = discrete_rate,
      mcse = discrete_mcse,
      wilson_low = discrete_wilson_low,
      wilson_high = discrete_wilson_high
    )
  ]
), fill = TRUE)
if (nrow(requested_table) != 3L) {
  stop("Requested-probability extraction failed.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Empirical delta inversion
# -------------------------------------------------------------------------

RUN_STAGE <- "empirical_delta_inversion"
mapping <- scan_results[
  n_training_generated == empirical_n & pi_minor == 0.20,
  .(
    n = .N,
    observed_delta_mean = mean(mahalanobis_centroid_delta),
    observed_delta_sd = sd(mahalanobis_centroid_delta),
    observed_delta_mcse = sd(mahalanobis_centroid_delta) / sqrt(.N),
    observed_delta_p025 = safe_quantile(mahalanobis_centroid_delta, 0.025),
    observed_delta_p975 = safe_quantile(mahalanobis_centroid_delta, 0.975)
  ),
  by = delta
]
setorder(mapping, delta)
raw_monotonicity_violations <- sum(diff(mapping$observed_delta_mean) < 0)
iso <- isoreg(mapping$delta, mapping$observed_delta_mean)
mapping[, observed_delta_isotonic := as.numeric(iso$yf)]
delta_inversion <- invert_monotone(
  mapping$delta, mapping$observed_delta_isotonic, empirical_delta
)
inversion_status <- if (is.finite(delta_inversion)) {
  "point_estimate_bracketed"
} else {
  "point_estimate_not_bracketed_no_extrapolation"
}

mapping_split <- split(
  scan_results[
    n_training_generated == empirical_n & pi_minor == 0.20,
    mahalanobis_centroid_delta
  ],
  scan_results[
    n_training_generated == empirical_n & pi_minor == 0.20,
    as.character(delta)
  ]
)
set.seed(INVERSION_BOOT_SEED)
inversion_boot <- replicate(B_BOOT, {
  boot_means <- vapply(as.character(mapping$delta), function(nm) {
    x <- mapping_split[[nm]]
    mean(sample(x, length(x), replace = TRUE))
  }, numeric(1))
  iso_b <- isoreg(mapping$delta, boot_means)
  invert_monotone(mapping$delta, as.numeric(iso_b$yf), empirical_delta)
})
valid_inversion_boot <- inversion_boot[is.finite(inversion_boot)]

relation_to_critical <- if (!is.finite(delta_inversion)) {
  "not_estimable"
} else if (delta_inversion > UNIMODALITY_CONTEXT_DELTA) {
  "point_estimate_above_contextual_critical_value"
} else if (delta_inversion < UNIMODALITY_CONTEXT_DELTA) {
  "point_estimate_below_contextual_critical_value"
} else {
  "point_estimate_equal_to_contextual_critical_value"
}

inversion_table <- data.table(
  empirical_source = normalizePath(
    EMPIRICAL_DELTA_SOURCE, winslash = "/", mustWork = TRUE
  ),
  empirical_n = empirical_n,
  empirical_cluster_1_n = empirical_cluster_1_n,
  empirical_cluster_2_n = empirical_cluster_2_n,
  empirical_minor_prevalence = empirical_pi,
  simulated_pi_used = 0.20,
  pi_approximation_absolute_difference = abs(empirical_pi - 0.20),
  empirical_pooled_within_mahalanobis = empirical_delta,
  inversion_mapping =
    "isotonic_setting_means_then_piecewise_linear_inverse",
  isotonic_plateau_inverse_policy = "midpoint_of_injected_delta_values_on_plateau",
  delta_inversion = delta_inversion,
  inversion_status = inversion_status,
  raw_mean_monotonicity_violations = raw_monotonicity_violations,
  bootstrap_B = B_BOOT,
  bootstrap_seed = INVERSION_BOOT_SEED,
  bootstrap_valid_n = length(valid_inversion_boot),
  bootstrap_not_bracketed_n = B_BOOT - length(valid_inversion_boot),
  bootstrap_valid_fraction = length(valid_inversion_boot) / B_BOOT,
  bootstrap_interval_available = is.finite(delta_inversion) &&
    length(valid_inversion_boot) > 1L,
  bootstrap_interval_scope = if (is.finite(delta_inversion) &&
                                   length(valid_inversion_boot) > 1L) {
    "conditional_on_successful_bootstrap_bracketing"
  } else if (is.finite(delta_inversion)) {
    "not_reported_fewer_than_two_successfully_bracketed_bootstrap_replicates"
  } else {
    "not_reported_because_point_estimate_not_bracketed"
  },
  bootstrap_se = if (is.finite(delta_inversion) &&
                     length(valid_inversion_boot) > 1L) {
    sd(valid_inversion_boot)
  } else {
    NA_real_
  },
  bootstrap_ci_low = if (is.finite(delta_inversion) &&
                         length(valid_inversion_boot) > 1L) {
    safe_quantile(valid_inversion_boot, 0.025)
  } else {
    NA_real_
  },
  bootstrap_ci_high = if (is.finite(delta_inversion) &&
                          length(valid_inversion_boot) > 1L) {
    safe_quantile(valid_inversion_boot, 0.975)
  } else {
    NA_real_
  },
  contextual_unimodality_critical_delta = UNIMODALITY_CONTEXT_DELTA,
  relation_to_contextual_critical_delta = relation_to_critical,
  reporting_boundary = paste(
    "DGM-scale calibration only; not a biological effect estimate and not",
    "an empirical subtype verdict"
  )
)

# Every planned task must have one readable, signature-matched final checkpoint.
RUN_STAGE <- "checkpoint_inventory"
checkpoint_inventory <- rbindlist(lapply(seq_len(nrow(resolved_manifest)), function(i) {
  task <- resolved_manifest[i]
  task_type <- as.character(task$task_type)
  path <- checkpoint_path_for(as.character(task$task_id), task_type)
  expected_signature <- expected_task_signature(task, task_type)
  obj <- if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(e) e)
  } else {
    NULL
  }
  readable <- is.list(obj) && !inherits(obj, "error") &&
    is.data.frame(obj$row) && nrow(obj$row) == 1L
  observed_signature <- if (readable) as.character(obj$signature) else NA_character_
  observed_schema <- if (readable) as.character(obj$checkpoint_schema) else NA_character_
  observed_rule_version <- if (readable) as.character(obj$rule_version) else NA_character_
  observed_run_id <- if (readable) as.character(obj$run_id) else NA_character_
  row_status <- if (readable && "status" %in% names(obj$row)) {
    as.character(obj$row$status[[1L]])
  } else {
    NA_character_
  }
  row_task_id <- if (readable && "task_id" %in% names(obj$row)) {
    as.character(obj$row$task_id[[1L]])
  } else {
    NA_character_
  }
  data.table(
    task_id = as.character(task$task_id),
    task_type = task_type,
    checkpoint_path = normalizePath(
      path, winslash = "/", mustWork = FALSE
    ),
    exists = file.exists(path),
    readable = readable,
    expected_signature = expected_signature,
    observed_signature = observed_signature,
    signature_match = readable && identical(
      observed_signature, expected_signature
    ),
    checkpoint_schema = observed_schema,
    checkpoint_schema_match = readable && identical(
      observed_schema, CHECKPOINT_SCHEMA
    ),
    rule_version = observed_rule_version,
    rule_version_match = readable && identical(
      observed_rule_version, RULE_VERSION
    ),
    run_id = observed_run_id,
    run_id_match = readable && identical(observed_run_id, RUN_ID),
    row_task_id = row_task_id,
    row_task_id_match = readable && identical(
      row_task_id, as.character(task$task_id)
    ),
    status = row_status,
    accepted_status = readable && row_status %in% accepted_task_status,
    bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
    sha256 = if (file.exists(path)) sha256_file(path) else NA_character_
  )
}))
checkpoint_inventory[, checkpoint_valid :=
  exists & readable & signature_match & checkpoint_schema_match &
    rule_version_match & run_id_match & row_task_id_match & accepted_status]

# -------------------------------------------------------------------------
# Hard QC before final scientific output
# -------------------------------------------------------------------------

RUN_STAGE <- "hard_qc"
resolved_threshold_check <- merge(
  scan_results[, .(
    task_id, n_training_generated, observed_gate = gate_threshold
  )],
  gate_table[, .(
    n_training_generated, expected_gate = q95
  )],
  by = "n_training_generated", all.x = TRUE
)

add_qc("fixed authoritative pre-run SHA-256 map verified",
       all(fixed_hash_rows$fixed_hash_match),
       sum(fixed_hash_rows$fixed_hash_match), nrow(fixed_hash_rows))
add_qc("locked S7 engine exposes 33 unique named features",
       length(engine_feature_names) == EMPIRICAL_EXPECTED_P &&
         !anyDuplicated(engine_feature_names) &&
         all(nzchar(engine_feature_names)),
       paste(length(engine_feature_names), uniqueN(engine_feature_names), sep = "/"),
       paste(EMPIRICAL_EXPECTED_P, EMPIRICAL_EXPECTED_P, sep = "/"))
add_qc("locked empirical N1 matrix exposes 33 unique named features",
       length(empirical_feature_names) == EMPIRICAL_EXPECTED_P &&
         !anyDuplicated(empirical_feature_names) &&
         all(nzchar(empirical_feature_names)),
       paste(length(empirical_feature_names), uniqueN(empirical_feature_names),
             sep = "/"),
       paste(EMPIRICAL_EXPECTED_P, EMPIRICAL_EXPECTED_P, sep = "/"))
add_qc("feature-space roles are explicit and exact-name identity is not claimed",
       nrow(feature_space_audit) == 2L * EMPIRICAL_EXPECTED_P &&
         all(feature_space_audit$exact_feature_name_identity_required %in% FALSE) &&
         uniqueN(feature_space_audit$feature_space) == 2L,
       paste0("rows=", nrow(feature_space_audit),
              "; exact_overlap=", length(feature_name_overlap)),
       paste0("rows=", 2L * EMPIRICAL_EXPECTED_P,
              "; two role-separated feature spaces"))
add_qc("N1 PIT registry uses empirical feature names only",
       setequal(
         seed_registry[module == "randomized_PIT", feature],
         empirical_feature_names
       ) &&
         !any(seed_registry[module == "randomized_PIT", feature] %in%
                setdiff(engine_feature_names, empirical_feature_names)),
       paste(seed_registry[module == "randomized_PIT", feature], collapse = ";"),
       paste(empirical_feature_names, collapse = ";"))
add_qc("empirical source has exactly one compatible valid row",
       nrow(empirical) == 1L, nrow(empirical), 1)
add_qc("fixed empirical N1 source dimensions locked",
       n1_fit$source_n == EMPIRICAL_EXPECTED_N &&
         n1_fit$source_p == EMPIRICAL_EXPECTED_P,
       paste(n1_fit$source_n, n1_fit$source_p, sep = "x"),
       paste(EMPIRICAL_EXPECTED_N, EMPIRICAL_EXPECTED_P, sep = "x"))
add_qc("fixed empirical N1 PIT method recorded",
       identical(
         n1_fit$pit_method,
         "randomized_Rueschendorf_distribution_transform"
       ), n1_fit$pit_method,
       "randomized_Rueschendorf_distribution_transform")
add_qc("global actual-seed registry rows", nrow(seed_registry) == expected_seed_rows,
       nrow(seed_registry), expected_seed_rows)
add_qc("global actual-seed values valid", invalid_seed_n == 0L,
       invalid_seed_n, 0)
add_qc("global actual-seed values unique", nrow(duplicate_seed_values) == 0L,
       nrow(duplicate_seed_values), 0)
add_qc("cross-module seed intersections absent",
       nrow(cross_module_seed_values) == 0L,
       nrow(cross_module_seed_values), 0)
add_qc("gate task count complete", nrow(gate_results) == expected_gate_tasks,
       nrow(gate_results), expected_gate_tasks)
add_qc("gate task failures absent", nrow(gate_failures) == 0L,
       nrow(gate_failures), 0)
add_qc("gate tasks use strict Lloyd without fallback or warnings",
       all(gate_results$kmeans_primary_algorithm_used == "Lloyd") &&
         all(gate_results$kmeans_primary_nstart == 25L) &&
         all(gate_results$kmeans_primary_itermax == 100L) &&
         all(gate_results$kmeans_primary_fallback_used %in% FALSE) &&
         all(gate_results$kmeans_primary_warning_count == 0L) &&
         all(gate_results$task_warning_count == 0L),
       paste0(
         "fallback=", sum(gate_results$kmeans_primary_fallback_used %in% TRUE),
         ";warnings=", sum(
           gate_results$kmeans_primary_warning_count +
             gate_results$task_warning_count
         )
       ), "fallback=0;warnings=0")
add_qc("gate source is fixed empirical N1",
       all(gate_results$null_id == "N1_empirical_gaussian_copula"),
       paste(unique(gate_results$null_id), collapse = ";"),
       "N1_empirical_gaussian_copula")
add_qc("sample-size-specific gate rows complete",
       nrow(gate_table) == length(N_GRID), nrow(gate_table), length(N_GRID))
add_qc("each n uses B_gate completed replicates",
       all(gate_table$n_completed == B_GATE),
       paste(gate_table$n_completed, collapse = ";"), B_GATE)
gate_uncertainty_ok <- all(is.finite(gate_table$q95)) &&
  all(is.finite(gate_table$q95_bootstrap_se)) &&
  all(is.finite(gate_table$q95_bootstrap_ci_low)) &&
  all(is.finite(gate_table$q95_bootstrap_ci_high)) &&
  all(gate_table$q95_bootstrap_ci_low <= gate_table$q95_bootstrap_ci_high)
add_qc("gate q95 estimates and bootstrap uncertainty finite and ordered",
       gate_uncertainty_ok, gate_uncertainty_ok, TRUE)
add_qc("scan task count complete", nrow(scan_results) == expected_scan_tasks,
       nrow(scan_results), expected_scan_tasks)
add_qc("scan task failures absent", nrow(scan_failures) == 0L,
       nrow(scan_failures), 0)
add_qc("scan tasks use strict Lloyd without fallback or warnings",
       all(scan_results$kmeans_primary_algorithm_used == "Lloyd") &&
         all(scan_results$kmeans_primary_nstart == 25L) &&
         all(scan_results$kmeans_primary_itermax == 100L) &&
         all(scan_results$kmeans_primary_fallback_used %in% FALSE) &&
         all(scan_results$kmeans_primary_warning_count == 0L) &&
         all(scan_results$task_warning_count == 0L),
       paste0(
         "fallback=", sum(scan_results$kmeans_primary_fallback_used %in% TRUE),
         ";warnings=", sum(
           scan_results$kmeans_primary_warning_count +
             scan_results$task_warning_count
         )
       ), "fallback=0;warnings=0")
add_qc("checkpoint inventory count complete",
       nrow(checkpoint_inventory) == nrow(resolved_manifest),
       nrow(checkpoint_inventory), nrow(resolved_manifest))
add_qc("all final checkpoints readable and signature matched",
       all(checkpoint_inventory$checkpoint_valid),
       sum(checkpoint_inventory$checkpoint_valid), nrow(checkpoint_inventory))
add_qc("every scan task uses corresponding-n q95",
       all(abs(
         resolved_threshold_check$observed_gate -
           resolved_threshold_check$expected_gate
       ) < 1e-12),
       max(abs(
         resolved_threshold_check$observed_gate -
           resolved_threshold_check$expected_gate
       )), 0)
add_qc("no recovery ARI enters verdict",
       all(scan_results$oracle_ari_used_in_verdict %in% FALSE),
       sum(scan_results$oracle_ari_used_in_verdict %in% TRUE), 0)
add_qc("three-state surface complete",
       nrow(surface_rows) == length(N_GRID) * length(DELTA_GRID) *
         length(PI_GRID),
       nrow(surface_rows),
       length(N_GRID) * length(DELTA_GRID) * length(PI_GRID))
add_qc("three-state probabilities sum to one",
       all(abs(
         surface_rows$no_evidence_rate + surface_rows$inconclusive_rate +
           surface_rows$discrete_rate - 1
       ) < 1e-12),
       max(abs(
         surface_rows$no_evidence_rate + surface_rows$inconclusive_rate +
           surface_rows$discrete_rate - 1
       )), 0)
state_count_sum <- surface_rows$no_evidence_count +
  surface_rows$inconclusive_count + surface_rows$discrete_count
add_qc("three-state counts sum to completed n",
       all(state_count_sum == surface_rows$n_completed),
       max(abs(state_count_sum - surface_rows$n_completed)), 0)
wilson_long_qc <- rbindlist(list(
  surface_rows[, .(
    rate = no_evidence_rate,
    low = no_evidence_wilson_low,
    high = no_evidence_wilson_high
  )],
  surface_rows[, .(
    rate = inconclusive_rate,
    low = inconclusive_wilson_low,
    high = inconclusive_wilson_high
  )],
  surface_rows[, .(
    rate = discrete_rate,
    low = discrete_wilson_low,
    high = discrete_wilson_high
  )]
))
wilson_surface_ok <- all(is.finite(unlist(wilson_long_qc))) &&
  all(wilson_long_qc$low >= 0 & wilson_long_qc$high <= 1) &&
  all(wilson_long_qc$low <= wilson_long_qc$rate) &&
  all(wilson_long_qc$rate <= wilson_long_qc$high)
add_qc("Wilson intervals are finite, legal, ordered, and contain rates",
       wilson_surface_ok, wilson_surface_ok, TRUE)
wilson_zero <- wilson(0L, 100L)
wilson_full <- wilson(100L, 100L)
wilson_boundary_ok <- length(wilson_zero) == 2L &&
  length(wilson_full) == 2L &&
  all(is.finite(c(wilson_zero, wilson_full))) &&
  abs(wilson_zero[[1L]]) < 1e-15 &&
  wilson_zero[[2L]] > 0 && wilson_zero[[2L]] < 1 &&
  wilson_full[[1L]] > 0 && wilson_full[[1L]] < 1 &&
  abs(wilson_full[[2L]] - 1) < 1e-15
add_qc("Wilson 0/100 and 100/100 boundary behavior valid",
       wilson_boundary_ok,
       paste(signif(c(wilson_zero, wilson_full), 8), collapse = ";"),
       "0;<1;>0;1")
add_qc("requested probabilities extracted", nrow(requested_table) == 3L,
       nrow(requested_table), 3)
add_qc("requested probabilities use n=19049",
       all(requested_table$n_training_generated == REQUESTED_PROBABILITY_N),
       paste(unique(requested_table$n_training_generated), collapse = ";"),
       REQUESTED_PROBABILITY_N)
expected_c0_models <- c(
  "primary_inverse_n_inverse_mc_variance",
  "sensitivity_inverse_sqrt_n_equal_weight"
)
add_qc("c0 table contains exactly two prespecified models",
       nrow(c0_table) == 2L &&
         setequal(c0_table$model_id, expected_c0_models),
       paste(c0_table$model_id, collapse = ";"),
       paste(expected_c0_models, collapse = ";"))
c0_numeric_ok <- all(is.finite(c0_table$c0)) &&
  all(is.finite(c0_table$c1)) &&
  all(is.finite(c0_table$model_r_squared)) &&
  all(is.finite(c0_table$c0_bootstrap_se)) &&
  all(is.finite(c0_table$c0_bootstrap_ci_low)) &&
  all(is.finite(c0_table$c0_bootstrap_ci_high)) &&
  all(c0_table$c0_bootstrap_ci_low <= c0_table$c0_bootstrap_ci_high)
add_qc("both c0 model estimates and uncertainty are finite and ordered",
       c0_numeric_ok, c0_numeric_ok, TRUE)
add_qc("inversion mapping contains every delta setting exactly once",
       nrow(mapping) == length(DELTA_GRID) &&
         setequal(mapping$delta, DELTA_GRID),
       paste(mapping$delta, collapse = ";"),
       paste(DELTA_GRID, collapse = ";"))
add_qc("inversion mapping uses B_scan replicates per delta",
       all(mapping$n == B_SCAN), paste(mapping$n, collapse = ";"), B_SCAN)
add_qc("isotonic inversion mapping is nondecreasing",
       all(diff(mapping$observed_delta_isotonic) >= -1e-12),
       min(diff(mapping$observed_delta_isotonic)), ">=-1e-12")
add_qc("inversion bootstrap accounting complete",
       inversion_table$bootstrap_valid_n +
         inversion_table$bootstrap_not_bracketed_n == B_BOOT,
       inversion_table$bootstrap_valid_n +
         inversion_table$bootstrap_not_bracketed_n, B_BOOT)
inversion_ci_logic_ok <- if (is.finite(delta_inversion)) {
  identical(inversion_status, "point_estimate_bracketed") &&
    if (length(valid_inversion_boot) > 1L) {
      isTRUE(inversion_table$bootstrap_interval_available) &&
        identical(
          inversion_table$bootstrap_interval_scope,
          "conditional_on_successful_bootstrap_bracketing"
        ) &&
      all(is.finite(c(
        inversion_table$bootstrap_se,
        inversion_table$bootstrap_ci_low,
        inversion_table$bootstrap_ci_high
      ))) &&
        inversion_table$bootstrap_ci_low <= inversion_table$bootstrap_ci_high
    } else {
      !isTRUE(inversion_table$bootstrap_interval_available) &&
        identical(
          inversion_table$bootstrap_interval_scope,
          "not_reported_fewer_than_two_successfully_bracketed_bootstrap_replicates"
        ) &&
      all(is.na(c(
        inversion_table$bootstrap_se,
        inversion_table$bootstrap_ci_low,
        inversion_table$bootstrap_ci_high
      )))
    }
} else {
  identical(
    inversion_status,
    "point_estimate_not_bracketed_no_extrapolation"
  ) && !isTRUE(inversion_table$bootstrap_interval_available) &&
    identical(
      inversion_table$bootstrap_interval_scope,
      "not_reported_because_point_estimate_not_bracketed"
    ) && all(is.na(c(
    inversion_table$bootstrap_se,
    inversion_table$bootstrap_ci_low,
    inversion_table$bootstrap_ci_high
  )))
}
add_qc("point inversion status and conditional interval logic consistent",
       inversion_ci_logic_ok, inversion_ci_logic_ok, TRUE)
add_qc("empirical source n locked", empirical_n == EMPIRICAL_EXPECTED_N,
       empirical_n, EMPIRICAL_EXPECTED_N)
add_qc("empirical source cluster 1 n locked",
       empirical_cluster_1_n == EMPIRICAL_EXPECTED_CLUSTER_1_N,
       empirical_cluster_1_n, EMPIRICAL_EXPECTED_CLUSTER_1_N)
add_qc("empirical source cluster 2 n locked",
       empirical_cluster_2_n == EMPIRICAL_EXPECTED_CLUSTER_2_N,
       empirical_cluster_2_n, EMPIRICAL_EXPECTED_CLUSTER_2_N)
add_qc("empirical cluster counts sum to n",
       empirical_cluster_1_n + empirical_cluster_2_n == empirical_n,
       empirical_cluster_1_n + empirical_cluster_2_n, empirical_n)
add_qc("empirical CSV matrix MD5 matches locked RDS",
       identical(
         empirical_matrix_md5_source,
         tolower(empirical_matrix_md5_observed)
       ) && identical(
         empirical_matrix_md5_source, EMPIRICAL_EXPECTED_MATRIX_MD5
       ), empirical_matrix_md5_source, EMPIRICAL_EXPECTED_MATRIX_MD5)
add_qc("empirical source delta locked",
       abs(empirical_delta - EMPIRICAL_EXPECTED_DELTA) <= 1e-10,
       signif(empirical_delta, 15), signif(EMPIRICAL_EXPECTED_DELTA, 15))
add_qc("future multisession used", N_WORKERS >= 2L || !IS_FORMAL,
       paste0("workers=", N_WORKERS), "formal workers >=2")
if (IS_FORMAL) {
  add_qc("formal gate task contract", expected_gate_tasks == 1000L,
         expected_gate_tasks, 1000)
  add_qc("formal scan task contract", expected_scan_tasks == 11000L,
         expected_scan_tasks, 11000)
}

qc <- rbindlist(qc_rows, fill = TRUE)
if (!all(qc$pass)) {
  atomic_fwrite(qc, OUTPUT$runtime_qc)
  atomic_fwrite(failure_table, OUTPUT$runtime_failures)
  stop("WP8 hard QC failed; no formal completion marker is allowed.",
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Publication figures (R only)
# -------------------------------------------------------------------------

RUN_STAGE <- "figure_construction"
COL_NO <- "#6F7782"
COL_INC <- "#C77A30"
COL_DISC <- "#3A5A8C"
COL_EMP <- "#1F3A5F"
COL_POINT <- "#2E2E2E"
COL_REF <- "grey70"

theme_pub <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "#333333"),
    strip.background = element_rect(
      fill = "#F1F3F4", colour = "#B8B8B8", linewidth = 0.35
    ),
    strip.text = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "grey35", hjust = 0),
    plot.caption = element_text(size = 7.7, colour = "grey35", hjust = 0,
                                lineheight = 1.12),
    panel.spacing = grid::unit(6, "pt"),
    plot.margin = margin(7, 10, 7, 10)
  )

gate_figure_source <- copy(gate_table)
p_gate <- ggplot(gate_figure_source, aes(x = n_training_generated, y = q95)) +
  geom_ribbon(
    aes(ymin = q95_bootstrap_ci_low, ymax = q95_bootstrap_ci_high),
    fill = COL_DISC, alpha = 0.14
  ) +
  geom_line(colour = COL_DISC, linewidth = 0.65) +
  geom_point(colour = COL_DISC, size = 2) +
  geom_line(aes(y = separation_mean), colour = "grey55", linewidth = 0.55,
            linetype = "22") +
  geom_point(aes(y = separation_mean), colour = COL_POINT, size = 1.7,
             shape = 21, fill = "white", stroke = 0.45) +
  geom_hline(yintercept = c0_estimate, colour = COL_EMP, linetype = "33",
             linewidth = 0.55) +
  scale_x_log10(
    breaks = N_GRID,
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    guide = guide_axis(n.dodge = 2)
  ) +
  labs(
    title = "Sample-size-specific empirical Gaussian-copula gates",
    subtitle = paste0(
      "Fixed N1 candidate null; B=", B_GATE,
      " per n; q95 uncertainty from ", B_BOOT, " bootstrap resamples"
    ),
    x = "Generated training sample size (log scale)",
    y = "Pooled-within Mahalanobis separation",
    caption = paste0(
      "Solid blue: corresponding-n q95 (bootstrap 95% interval); ",
      "dashed grey: mean.\n",
      "Dotted dark line: primary c0 asymptote from ",
      "E(Delta_hat|n)=c0+c1/n."
    )
  ) +
  theme_pub

surface_long <- rbindlist(list(
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "No discrete evidence", probability = no_evidence_rate,
    mcse = no_evidence_mcse,
    wilson_low = no_evidence_wilson_low,
    wilson_high = no_evidence_wilson_high
  )],
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "Inconclusive", probability = inconclusive_rate,
    mcse = inconclusive_mcse,
    wilson_low = inconclusive_wilson_low,
    wilson_high = inconclusive_wilson_high
  )],
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "Discrete evidence", probability = discrete_rate,
    mcse = discrete_mcse,
    wilson_low = discrete_wilson_low,
    wilson_high = discrete_wilson_high
  )]
))
surface_long[, state := factor(
  state,
  levels = c("No discrete evidence", "Inconclusive", "Discrete evidence")
)]
surface_long[, n_label := factor(
  n_training_generated,
  levels = N_GRID,
  labels = paste0("n=", format(N_GRID, big.mark = ","))
)]
surface_long[, pi_label := factor(
  pi_minor,
  levels = PI_GRID,
  labels = c("Minor prevalence 0.20", "Balanced prevalence 0.50")
)]
surface_expected_for_figure <- rbindlist(list(
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "No discrete evidence", probability = no_evidence_rate,
    mcse = no_evidence_mcse,
    wilson_low = no_evidence_wilson_low,
    wilson_high = no_evidence_wilson_high
  )],
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "Inconclusive", probability = inconclusive_rate,
    mcse = inconclusive_mcse,
    wilson_low = inconclusive_wilson_low,
    wilson_high = inconclusive_wilson_high
  )],
  surface_rows[, .(
    n_training_generated, delta, pi_minor,
    state = "Discrete evidence", probability = discrete_rate,
    mcse = discrete_mcse,
    wilson_low = discrete_wilson_low,
    wilson_high = discrete_wilson_high
  )]
))
surface_observed_for_figure <- surface_long[, .(
  n_training_generated, delta, pi_minor,
  state = as.character(state), probability, mcse, wilson_low, wilson_high
)]
setorder(
  surface_expected_for_figure,
  n_training_generated, delta, pi_minor, state
)
setorder(
  surface_observed_for_figure,
  n_training_generated, delta, pi_minor, state
)
surface_figure_consistent <-
  identical(
    surface_observed_for_figure[, .(
      n_training_generated, delta, pi_minor, state
    )],
    surface_expected_for_figure[, .(
      n_training_generated, delta, pi_minor, state
    )]
  ) &&
  all(abs(
    as.matrix(surface_observed_for_figure[, .(
      probability, mcse, wilson_low, wilson_high
    )]) -
      as.matrix(surface_expected_for_figure[, .(
        probability, mcse, wilson_low, wilson_high
      )])
  ) < 1e-12)
add_qc(
  "figure source exactly matches the three-state surface",
  surface_figure_consistent, surface_figure_consistent, TRUE
)
qc <- rbindlist(qc_rows, fill = TRUE)
if (!all(qc$pass)) {
  atomic_fwrite(qc, OUTPUT$runtime_qc)
  stop("Figure-source consistency QC failed.", call. = FALSE)
}

p_surface <- ggplot(
  surface_long,
  aes(x = delta, y = probability, colour = state, linetype = state)
) +
  geom_hline(yintercept = c(0, 0.5, 1), colour = "grey90", linewidth = 0.3) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.35) +
  facet_grid(pi_label ~ n_label) +
  scale_colour_manual(values = c(
    "No discrete evidence" = COL_NO,
    "Inconclusive" = COL_INC,
    "Discrete evidence" = COL_DISC
  )) +
  scale_linetype_manual(values = c(
    "No discrete evidence" = "22",
    "Inconclusive" = "solid",
    "Discrete evidence" = "solid"
  )) +
  scale_x_continuous(breaks = c(0, 1, 2, 3.2)) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    title = "Three-state Domain 1 operating surface",
    subtitle = paste0(
      "S7-like K2 generator; B=", B_SCAN,
      " per setting; each n uses its fixed-N1 corresponding-n q95"
    ),
    x = "Injected separation (delta)", y = "Verdict probability",
    colour = NULL, linetype = NULL,
    caption = paste(
      "Dip evidence is evaluated on independent held-out observations.",
      "Recovery ARI does not enter the verdict."
    )
  ) +
  theme_pub +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 7.8),
    strip.text = element_text(size = 7.8)
  )

inversion_figure_source <- copy(mapping)
p_inversion <- ggplot(
  inversion_figure_source,
  aes(x = delta, y = observed_delta_mean)
) +
  geom_ribbon(
    aes(ymin = observed_delta_p025, ymax = observed_delta_p975),
    fill = COL_DISC, alpha = 0.12
  ) +
  geom_line(aes(y = observed_delta_isotonic), colour = COL_DISC,
            linewidth = 0.75) +
  geom_point(colour = COL_POINT, size = 1.8) +
  geom_hline(yintercept = empirical_delta, colour = COL_EMP,
             linetype = "22", linewidth = 0.6) +
  geom_vline(xintercept = UNIMODALITY_CONTEXT_DELTA, colour = COL_REF,
             linetype = "33", linewidth = 0.5) +
  {if (is.finite(delta_inversion)) {
    geom_vline(xintercept = delta_inversion, colour = COL_INC,
               linetype = "solid", linewidth = 0.6)
  }} +
  annotate(
    "text", x = min(mapping$delta), y = empirical_delta,
    label = paste0("Empirical Delta_hat = ", sprintf("%.3f", empirical_delta)),
    hjust = 0, vjust = -0.6, size = 3, colour = COL_EMP
  ) +
  labs(
    title = "Empirical separation mapped to the injected-delta scale",
    subtitle = paste0(
      "n=", format(empirical_n, big.mark = ","),
      "; simulated pi=0.20; isotonic mean mapping with bootstrap uncertainty"
    ),
    x = "Injected separation (delta)",
    y = "Observed pooled-within Mahalanobis separation",
    caption = paste0(
      "Blue: isotonic setting means; band: empirical 2.5th-97.5th ",
      "percentiles. Amber: inverted empirical estimate.\n",
      "Grey: contextual critical delta 2.98. This is DGM-scale calibration, ",
      "not a biological effect estimate."
    )
  ) +
  theme_pub

save_figure <- function(plot, pdf_path, png_path, tiff_path,
                        width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(
    pdf_path, width = width_in, height = height_in, family = "sans"
  )
  print(plot)
  grDevices::dev.off()
  grDevices::png(
    png_path, width = width_in, height = height_in,
    units = "in", res = 600, type = "cairo-png", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  grDevices::tiff(
    tiff_path, width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

# -------------------------------------------------------------------------
# Final writes and provenance
# -------------------------------------------------------------------------

RUN_STAGE <- "pre_output_integrity_recheck"
input_hash_end <- rbindlist(lapply(names(INPUT), function(nm) {
  data.table(input_name = nm, sha256_end = sha256_file(INPUT[[nm]]))
}))
input_hashes <- merge(input_hash_start, input_hash_end, by = "input_name",
                      all = TRUE)
input_hashes[, unchanged_during_run := sha256 == sha256_end]
if (!all(input_hashes$unchanged_during_run)) {
  stop("An authoritative WP8 input changed during execution.", call. = FALSE)
}
script_hash_end <- sha256_file(SCRIPT_PATH)
if (!identical(script_hash_start, script_hash_end)) {
  stop("The WP8 script changed during execution.", call. = FALSE)
}
add_qc("all authoritative inputs unchanged during run",
       all(input_hashes$unchanged_during_run),
       sum(input_hashes$unchanged_during_run), nrow(input_hashes))
add_qc("script unchanged during run",
       identical(script_hash_start, script_hash_end),
       script_hash_end, script_hash_start)

qc <- rbindlist(qc_rows, fill = TRUE)
if (!all(qc$pass)) {
  atomic_fwrite(qc, OUTPUT$runtime_qc)
  stop("Final integrity QC failed.", call. = FALSE)
}

RUN_STAGE <- "scientific_output_write"
if (IS_FORMAL) {
  present_now <- final_targets[file.exists(final_targets)]
  if (length(present_now)) {
    stop(
      "A formal final target appeared during computation; refusing overwrite:\n",
      paste(present_now, collapse = "\n"), call. = FALSE
    )
  }
}
atomic_fwrite(gate_table, OUTPUT$gate)
atomic_fwrite(c0_table, OUTPUT$c0)
atomic_fwrite(surface_rows, OUTPUT$surface)
atomic_fwrite(inversion_table, OUTPUT$inversion)
atomic_fwrite(requested_table, OUTPUT$requested)
atomic_fwrite(gate_results, OUTPUT$gate_by_repeat)
atomic_fwrite(scan_results, OUTPUT$scan_by_repeat)
atomic_fwrite(mapping, OUTPUT$inversion_mapping)
atomic_fwrite(n1_fit_summary, OUTPUT$n1_fit_summary)
atomic_fwrite(feature_space_audit, OUTPUT$feature_space_audit)
atomic_fwrite(gate_figure_source, OUTPUT$figure_gate_source)
atomic_fwrite(surface_long, OUTPUT$figure_surface_source)
atomic_fwrite(inversion_figure_source, OUTPUT$figure_inversion_source)
atomic_fwrite(qc, OUTPUT$qc)
atomic_fwrite(failure_table, OUTPUT$failures)
atomic_fwrite(resolved_manifest, OUTPUT$task_manifest)
atomic_fwrite(seed_registry, OUTPUT$seed_registry)
atomic_fwrite(input_hashes, OUTPUT$input_checksums)
atomic_fwrite(checkpoint_inventory, OUTPUT$checkpoint_inventory)
save_rds_atomic(n1_fit, OUTPUT$n1_fit_rds)

save_figure(
  p_gate, OUTPUT$gate_pdf, OUTPUT$gate_png, OUTPUT$gate_tiff, 183, 112
)
save_figure(
  p_surface, OUTPUT$surface_pdf, OUTPUT$surface_png, OUTPUT$surface_tiff,
  183, 150
)
save_figure(
  p_inversion, OUTPUT$inversion_pdf, OUTPUT$inversion_png,
  OUTPUT$inversion_tiff, 183, 112
)

run_metadata <- data.table(
  run_id = RUN_ID,
  run_mode = RUN_MODE,
  reportable_formal_result = IS_FORMAL,
  intentional_fail_fast_test_hook_active = FORCE_TOP_LEVEL_TEST_FAILURE,
  rule_version = RULE_VERSION,
  checkpoint_schema = CHECKPOINT_SCHEMA,
  n_grid = paste(N_GRID, collapse = ";"),
  delta_grid = paste(DELTA_GRID, collapse = ";"),
  pi_grid = paste(PI_GRID, collapse = ";"),
  evaluation_n = EVALUATION_N,
  B_gate = B_GATE,
  B_scan = B_SCAN,
  B_bootstrap = B_BOOT,
  kmeans_nstart = KMEANS_NSTART,
  ridge_multiplier = RIDGE_MULTIPLIER,
  dip_alpha = DIP_ALPHA,
  gate_quantile = GATE_QUANTILE,
  quantile_type = QUANTILE_TYPE,
        workers = N_WORKERS,
        future_scheduling = FUTURE_SCHEDULING,
        parallel_backend = "future::multisession",
  future_version = as.character(packageVersion("future")),
  future_apply_version = as.character(packageVersion("future.apply")),
  seed_registry_rows = nrow(seed_registry),
  expected_seed_registry_rows = expected_seed_rows,
  invalid_seed_rows = invalid_seed_n,
  duplicated_seed_values = nrow(duplicate_seed_values),
  cross_module_seed_values = nrow(cross_module_seed_values),
  empirical_n = empirical_n,
  empirical_minor_prevalence = empirical_pi,
  empirical_delta = empirical_delta,
  requested_probability_n = REQUESTED_PROBABILITY_N,
  gate_null = "fixed_empirical_gaussian_copula_N1",
  n1_source_n = n1_fit$source_n,
  n1_source_p = n1_fit$source_p,
  empirical_feature_space_n = length(empirical_feature_names),
  s7_engine_feature_space_n = length(engine_feature_names),
  exact_feature_name_overlap_n = length(feature_name_overlap),
  exact_feature_name_identity_required = FALSE,
  empirical_feature_space_role =
    "fixed empirical N1 gate calibration from locked primary33 matrix",
  s7_feature_space_role =
    "locked semi-synthetic alternative for detection surface and generator-specific inversion",
  n1_pit_method = n1_fit$pit_method,
  n1_nearpd_used = n1_fit$nearpd_used,
  c0_primary_model = c0_table[is_primary == TRUE, model],
  c0_primary = c0_estimate,
  c0_sensitivity_model = c0_table[is_primary == FALSE, model],
  c0_sensitivity = c0_sensitivity_estimate,
  delta_inversion = delta_inversion,
  inversion_status = inversion_status,
  raw_mapping_monotonicity_violations = raw_monotonicity_violations,
  gate_from_corresponding_n_only = TRUE,
  pooled_or_worst_case_gate_used = FALSE,
  gate_fallback_task_n = sum(
    gate_results$kmeans_primary_fallback_used %in% TRUE
  ),
  scan_fallback_task_n = sum(
    scan_results$kmeans_primary_fallback_used %in% TRUE
  ),
  gate_warning_task_n = sum(
    gate_results$task_warning_count > 0L |
      gate_results$kmeans_primary_warning_count > 0L
  ),
  scan_warning_task_n = sum(
    scan_results$task_warning_count > 0L |
      scan_results$kmeans_primary_warning_count > 0L
  ),
  checkpoint_valid_n = sum(checkpoint_inventory$checkpoint_valid),
  checkpoint_planned_n = nrow(checkpoint_inventory),
  complete_input_chain_sha256 = input_chain_sha256,
        formal_exclusive_lock_used = isTRUE(FORMAL_LOCK_EVER_ACQUIRED),
  formal_lock_stale_recovery_minimum_minutes = FORMAL_LOCK_STALE_MINUTES,
  recovery_ari_used_in_verdict = FALSE,
  wp7_outputs_modified = FALSE,
  calibration_scope =
    "pipeline_specific_fixed_K2_processed_space_Lloyd_reference",
  uncertainty_propagation_scope =
    "does_not_propagate_original_cleaning_MICE_K_selection_or_full_pipeline_uncertainty",
  empirical_d1_verdict_changed_by_wp8 = FALSE,
  script_35_authorized_by_wp8_completion = FALSE,
  script_sha256_start = script_hash_start,
  script_sha256_end = script_hash_end,
  frozen_spec_sha256 = NA_character_
)
atomic_fwrite(run_metadata, OUTPUT$run_metadata)
atomic_write_lines(capture.output(sessionInfo()), OUTPUT$session)

log_lines <- c(
  paste0(
    "# WP8 ", if (IS_FORMAL) "formal" else "nonformal validation",
    " gate, detection surface, and delta inversion"
  ),
  "",
  paste0("- Run ID: `", RUN_ID, "`."),
  paste0("- Run mode: `", RUN_MODE, "`."),
  paste0("- Parallel backend: `future::multisession`, workers=", N_WORKERS, "."),
  paste0("- Gate tasks: ", nrow(gate_results), "/", expected_gate_tasks, "."),
  paste0("- Scan tasks: ", nrow(scan_results), "/", expected_scan_tasks, "."),
  paste0("- Failures: ", nrow(failure_table), "."),
  paste0("- Sample-size-specific q95 values: ", paste(
    paste0(gate_table$n_training_generated, "=", sprintf("%.5f", gate_table$q95)),
    collapse = "; "
  ), "."),
  paste0(
    "- Primary N1 c0: ", sprintf("%.5f", c0_estimate),
    " (bootstrap 95% ",
    sprintf("%.5f", c0_table[is_primary == TRUE, c0_bootstrap_ci_low]),
    " to ",
    sprintf("%.5f", c0_table[is_primary == TRUE, c0_bootstrap_ci_high]),
    "); sensitivity c0=", sprintf("%.5f", c0_sensitivity_estimate), "."
  ),
  paste0(
    "- Empirical delta inversion: ",
    if (is.finite(delta_inversion)) sprintf("%.5f", delta_inversion) else "not bracketed",
    "; bootstrap valid=", length(valid_inversion_boot), "/", B_BOOT, "."
  ),
  "",
  "## Requested operating points",
  "",
  paste0(
    "- ", requested_table$request, ": ",
    sprintf("%.4f", requested_table$probability), " (Wilson 95% ",
    sprintf("%.4f", requested_table$wilson_low), " to ",
    sprintf("%.4f", requested_table$wilson_high), ")."
  ),
  "",
  "## Interpretation boundary",
  "",
  paste(
    "WP8 calibrates a pipeline- and sample-size-specific diagnostic surface.",
    "The inversion is a DGM-scale calibration, not a biological effect estimate.",
    "WP8 does not by itself change the locked empirical D1 verdict."
  ),
  paste0(
    "- Feature-space boundary: empirical N1 p=", length(empirical_feature_names),
    "; S7 engine p=", length(engine_feature_names),
    "; exact-name overlap=", length(feature_name_overlap),
    "; exact name identity is neither required nor claimed."
  ),
  "",
  "## Reproducibility",
  "",
  paste0("- Actual-seed registry: ", nrow(seed_registry), "/",
         expected_seed_rows, "; duplicate values=0; cross-module intersections=0."),
  paste0("- Hard QC: ", sum(qc$pass), "/", nrow(qc), " passed."),
  paste0(
    "- K-means fallbacks: gate=",
    sum(gate_results$kmeans_primary_fallback_used %in% TRUE),
    "; scan=", sum(scan_results$kmeans_primary_fallback_used %in% TRUE), "."
  ),
  paste0(
    "- Warning-bearing tasks: gate=",
    sum(
      gate_results$task_warning_count > 0L |
        gate_results$kmeans_primary_warning_count > 0L
    ),
    "; scan=", sum(
      scan_results$task_warning_count > 0L |
        scan_results$kmeans_primary_warning_count > 0L
    ), "."
  ),
  "- Every task has an atomic, signature-validated checkpoint.",
  ""
)
atomic_write_lines(log_lines, OUTPUT$log)

manifest_targets <- unlist(OUTPUT[c(
  "gate", "c0", "surface", "inversion", "requested",
  "gate_by_repeat", "scan_by_repeat", "inversion_mapping", "n1_fit_summary",
  "feature_space_audit",
  "figure_gate_source", "figure_surface_source", "figure_inversion_source",
  "gate_pdf", "gate_png", "gate_tiff", "surface_pdf", "surface_png",
  "surface_tiff", "inversion_pdf", "inversion_png", "inversion_tiff",
  "log", "qc", "failures", "session", "task_manifest", "seed_registry",
  "input_checksums", "run_metadata", "n1_fit_rds", "checkpoint_inventory"
)], use.names = FALSE)
output_manifest <- rbindlist(lapply(manifest_targets, function(path) {
  if (!file.exists(path)) stop("Missing manifest target: ", path, call. = FALSE)
  data.table(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = file.info(path)$size,
    sha256 = sha256_file(path)
  )
}))
if (any(!is.finite(output_manifest$bytes) | output_manifest$bytes <= 0)) {
  stop("At least one WP8 output is empty; output manifest was not released.",
       call. = FALSE)
}
atomic_fwrite(output_manifest, OUTPUT$output_manifest)

manifest_verify <- output_manifest[, .(
  path, expected_sha256 = sha256,
  observed_sha256 = vapply(path, sha256_file, character(1))
)]
if (!all(manifest_verify$expected_sha256 == manifest_verify$observed_sha256)) {
  stop("WP8 output-manifest verification failed.", call. = FALSE)
}
completion_lines <- c(
  if (IS_FORMAL) "RUN COMPLETED" else "NONFORMAL RUN COMPLETED",
  paste0("run_id=", RUN_ID),
  paste0("run_mode=", RUN_MODE),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("script_sha256=", script_hash_end),
  paste0("frozen_spec_sha256=", NA_character_),
  paste0("output_manifest_sha256=", sha256_file(OUTPUT$output_manifest)),
  paste0("gate_tasks_completed=", nrow(gate_results)),
  paste0("scan_tasks_completed=", nrow(scan_results)),
  paste0("expected_task_n=", expected_gate_tasks + expected_scan_tasks),
  paste0("completed_task_n=", nrow(gate_results) + nrow(scan_results)),
  "failed_task_n=0",
  "fallback_n=0",
  "cluster_warning_n=0",
  paste0("reportable_formal_result=", if (IS_FORMAL) "TRUE" else "FALSE"),
  "gate_null=FIXED_EMPIRICAL_GAUSSIAN_COPULA_N1",
  "gate_role=SAMPLE_SIZE_SPECIFIC_CORRESPONDING_N_Q95",
  paste0("empirical_feature_space_n=", length(empirical_feature_names)),
  paste0("s7_engine_feature_space_n=", length(engine_feature_names)),
  paste0("exact_feature_name_overlap_n=", length(feature_name_overlap)),
  "exact_feature_name_identity_claim=FALSE",
  paste0("requested_probability_n=", REQUESTED_PROBABILITY_N),
  paste0("empirical_inversion_n=", empirical_n),
  "calibration_scope=PIPELINE_SPECIFIC_FIXED_K2_PROCESSED_SPACE_LLOYD_REFERENCE",
  "uncertainty_propagation_scope=DOES_NOT_PROPAGATE_ORIGINAL_CLEANING_MICE_K_SELECTION_OR_FULL_PIPELINE_UNCERTAINTY",
  "empirical_verdict_changed_by_wp8=FALSE",
  "script_35_authorized_by_wp8_completion=FALSE"
)
atomic_write_lines(completion_lines, OUTPUT$run_completion)
lock_release_ok <- release_formal_lock()
if (IS_FORMAL && !isTRUE(lock_release_ok)) {
  stop("Formal output is complete but exclusive lock release failed; ",
       "the audit-root completion marker was not created.", call. = FALSE)
}
atomic_write_lines(completion_lines, OUTPUT$completion)

if (file.exists(OUTPUT$in_progress)) unlink(OUTPUT$in_progress)
if (file.exists(OUTPUT$partial)) unlink(OUTPUT$partial)
if (file.exists(OUTPUT$top_level_failure)) unlink(OUTPUT$top_level_failure)
RUN_STAGE <- "completed"
options(error = OLD_ERROR_OPTION)
future::plan(future::sequential)

message("WP8 completed successfully in mode: ", RUN_MODE)
message("Gate table: ", OUTPUT$gate)
message("Three-state surface: ", OUTPUT$surface)
message("Empirical inversion: ", OUTPUT$inversion)
message("Completion marker: ", OUTPUT$completion)
