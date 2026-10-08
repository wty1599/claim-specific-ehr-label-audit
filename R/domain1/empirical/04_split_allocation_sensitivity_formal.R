public_repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
#!/usr/bin/env Rscript

# WP9: Lloyd-aligned split-allocation sensitivity for the governed Domain 1
# three-state diagnostic. This script never selects a split post hoc.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(diptest)
  library(Matrix)
  library(mclust)
  library(digest)
  library(future)
  library(future.apply)
  library(ps)
})

REQUIRED_PACKAGES <- c(
  "data.table", "diptest", "Matrix", "mclust", "digest", "future",
  "future.apply", "ps"
)
missing_packages <- REQUIRED_PACKAGES[!vapply(
  REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

# -------------------------------------------------------------------------
# Frozen paths and run mode
# -------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
if (!nzchar(PROJECT_ROOT)) stop("Set EHR_AUDIT_WORK_ROOT to the authorized workspace.")
AUDIT_ROOT <- file.path(PROJECT_ROOT, "domain1_lloyd_alignment_20260720")
SCRIPT_PATH <- file.path(public_repo, "R/domain1/empirical/04_split_allocation_sensitivity_formal.R")

V4_UTILS <- file.path(public_repo, "R/domain1/empirical/00_domain1_revision_utils_v4_lloyd.R")
STRICT_WRAPPER <- file.path(public_repo, "R/domain1/empirical/01_audit_discreteness_v2.R")
V3_ROOT <- file.path(
  PROJECT_ROOT, "analysis_archive", "simulations",
  "Domain1_revision_v3_observable_heldout_20260717"
)
V3_LOADER <- file.path(public_repo, "R/domain1/empirical/dependencies/00_domain1_locked_engine_v3.R")
LOCKED_S7_ENGINE <- file.path(public_repo, "R/simulation/core/01_S7_discrete_outcome_null_formal.R")
LOCKED_MATRIX <- file.path(
  PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"
)
LOCKED_LABELS <- file.path(
  PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"
)


LOCKED_N1_FIT <- file.path(
  AUDIT_ROOT, "provenance",
  "29g_empirical_gaussian_copula_N1_fit_LLOYD.rds"
)
WP7_SEED_REGISTRY <- file.path(
  AUDIT_ROOT, "provenance", "29f_global_seed_registry_LLOYD.csv"
)
WP8_SEED_REGISTRY <- file.path(
  AUDIT_ROOT, "provenance", "29g_global_seed_registry_LLOYD.csv"
)







RUN_MODE <- tolower(Sys.getenv("D1_WP9_RUN_MODE", unset = "test"))
if (!RUN_MODE %in% c("test", "smoke", "formal")) {
  stop("D1_WP9_RUN_MODE must be test, smoke, or formal.", call. = FALSE)
}
IS_FORMAL <- identical(RUN_MODE, "formal")
RUN_ID_ENV <- Sys.getenv("D1_WP9_RUN_ID", unset = "")
RUN_ID <- if (nzchar(RUN_ID_ENV)) RUN_ID_ENV else paste0(
  format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid()
)
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$", RUN_ID)) {
  stop("D1_WP9_RUN_ID has an invalid format.", call. = FALSE)
}

NONFORMAL_ROOT <- file.path(AUDIT_ROOT, "nonformal_runs", "35")
FORMAL_RUN_ROOT <- file.path(
  AUDIT_ROOT, "checkpoints", "wp9", "35_lloyd_formal"
)
RUN_ROOT <- if (IS_FORMAL) FORMAL_RUN_ROOT else file.path(
  NONFORMAL_ROOT, RUN_MODE, RUN_ID
)
OUT_TABLE_DIR <- if (IS_FORMAL) file.path(AUDIT_ROOT, "tables") else
  file.path(RUN_ROOT, "tables")
OUT_LOG_DIR <- if (IS_FORMAL) file.path(AUDIT_ROOT, "logs") else
  file.path(RUN_ROOT, "logs")
OUT_PROVENANCE_DIR <- if (IS_FORMAL) file.path(AUDIT_ROOT, "provenance") else
  file.path(RUN_ROOT, "provenance")
CHECKPOINT_DIR <- file.path(RUN_ROOT, "tasks")
FAILURE_MARKER <- file.path(RUN_ROOT, "FAILED_DO_NOT_USE.txt")
LOCK_DIR <- file.path(
  AUDIT_ROOT, "checkpoints", "wp9", "35_lloyd_formal_exclusive.lock"
)
LOCK_OWNER <- file.path(LOCK_DIR, "owner.csv")

# -------------------------------------------------------------------------
# Frozen constants
# -------------------------------------------------------------------------

N_TOTAL <- 20049L
SPLITS <- data.table(
  split_order = 1:3,
  split_id = c(
    "train19049_eval1000", "train15049_eval5000",
    "train10024_eval10025"
  ),
  split_label = c(
    "Train 19,049 / evaluate 1,000",
    "Train 15,049 / evaluate 5,000",
    "Train 10,024 / evaluate 10,025"
  ),
  n_train = c(19049L, 15049L, 10024L),
  n_eval = c(1000L, 5000L, 10025L)
)
stopifnot(all(SPLITS$n_train + SPLITS$n_eval == N_TOTAL))

B_GATE <- switch(RUN_MODE, test = 2L, smoke = 20L, formal = 200L)
B_EVALUATION <- switch(RUN_MODE, test = 2L, smoke = 10L, formal = 100L)
B_BOOT <- switch(RUN_MODE, test = 100L, smoke = 200L, formal = 2000L)
SEED_NAMESPACE_START <- switch(
  RUN_MODE, test = 1100000000L, smoke = 1150000000L,
  formal = 1200000000L
)
EXPECTED_TASKS <- B_GATE + 3L * B_EVALUATION
EXPECTED_REPEAT_ROWS <- 3L * EXPECTED_TASKS
EXPECTED_SEEDS <- 5L * B_GATE + 14L * B_EVALUATION + 3L
if (!identical(
  c(EXPECTED_TASKS, EXPECTED_REPEAT_ROWS, EXPECTED_SEEDS),
  switch(
    RUN_MODE,
    test = c(8L, 24L, 41L),
    smoke = c(50L, 150L, 243L),
    formal = c(500L, 1500L, 2403L)
  )
)) stop("Mode-specific task/row/seed contract mismatch.", call. = FALSE)

KMEANS_NSTART <- 25L
KMEANS_ITERMAX <- 100L
RIDGE_MULTIPLIER <- 1e-4
DIP_ALPHA <- 0.05
GATE_QUANTILE <- 0.95
QUANTILE_TYPE <- 8L
S7_DELTA <- 3.2
S7_MINOR_PREVALENCE <- 0.20
RULE_VERSION <- "WP9_SPLIT_ALLOCATION_LLOYD_V2_20260722"
CHECKPOINT_SCHEMA <- "WP9_SPLIT_ALLOCATION_CHECKPOINT_V2_20260722"
SPEC_SHA256 <- "d0ca2811f6e3aa279cfcb413c71a6ce74fc4716f998ef4b1e50bf40f65fd8910"

D1_KMEANS_NSTART <- KMEANS_NSTART
D1_MAHALANOBIS_RIDGE_MULTIPLIER <- RIDGE_MULTIPLIER
D1_DIP_FDR_ALPHA <- DIP_ALPHA

# -------------------------------------------------------------------------
# Strict local file and hash utilities
# -------------------------------------------------------------------------

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

strict_tmp_path <- function(path) {
  token <- substr(digest::digest(
    list(path, Sys.getpid(), format(Sys.time(), "%OS6")),
    algo = "sha256", serialize = TRUE
  ), 1L, 16L)
  paste0(path, ".tmp_", Sys.getpid(), "_", token)
}

publish_with_writer <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Refusing to overwrite existing target: ", path, call. = FALSE)
  }
  tmp <- strict_tmp_path(path)
  if (file.exists(tmp) || dir.exists(tmp)) {
    stop("Temporary target collision: ", tmp, call. = FALSE)
  }
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)
  writer(tmp)
  if (!file.exists(tmp)) stop("Writer did not create temporary file.", call. = FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Target appeared before atomic publication: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Atomic no-overwrite publication failed: ", path, call. = FALSE)
  }
  invisible(path)
}

publish_csv <- function(x, path) publish_with_writer(
  path, function(tmp) data.table::fwrite(x, tmp, na = "")
)
publish_text <- function(x, path) publish_with_writer(
  path, function(tmp) writeLines(x, tmp, useBytes = TRUE)
)
publish_rds <- function(x, path) publish_with_writer(
  path, function(tmp) saveRDS(x, tmp, version = 3)
)

safe_q <- function(x, p) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(
    x, p, type = QUANTILE_TYPE, names = FALSE
  ))
}

numeric_summary <- function(x, prefix) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  out <- list(
    n = n,
    mean = if (n) mean(x) else NA_real_,
    sd = if (n > 1L) stats::sd(x) else NA_real_,
    mcse = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_,
    p025 = if (n) safe_q(x, 0.025) else NA_real_,
    p975 = if (n) safe_q(x, 0.975) else NA_real_
  )
  names(out) <- paste0(prefix, "_", names(out))
  as.data.table(out)
}

rate_summary <- function(x, n, prefix) {
  ci <- wilson_interval(x, n)
  out <- list(
    count = as.integer(x), n = as.integer(n),
    rate = if (n > 0L) x / n else NA_real_,
    ci_low = unname(ci["low"]), ci_high = unname(ci["high"])
  )
  names(out) <- paste0(prefix, "_", names(out))
  as.data.table(out)
}

read_key_values <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  keep <- grepl("^[A-Za-z0-9_.-]+=", lines)
  parts <- strsplit(lines[keep], "=", fixed = TRUE)
  keys <- vapply(parts, `[[`, character(1), 1L)
  vals <- vapply(parts, function(z) paste(z[-1L], collapse = "="), character(1))
  stats::setNames(vals, keys)
}

STATE <- new.env(parent = emptyenv())
STATE$stage <- "initialization"
STATE$lock_acquired <- FALSE
STATE$lock_ever_acquired <- FALSE
STATE$lock_token <- NA_character_
STATE$log_enabled <- FALSE
STATE$task_ids <- character()
STATE$completed_task_ids <- character()

log_msg <- function(...) {
  msg <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    paste0(..., collapse = "")
  )
  cat(msg, "\n")
  if (isTRUE(STATE$log_enabled)) {
    cat(msg, "\n", file = OUTPUT$run_log, append = TRUE)
  }
  invisible(msg)
}

write_failure_marker <- function(error_text) {
  if (file.exists(FAILURE_MARKER)) return(invisible(FALSE))
  dir.create(dirname(FAILURE_MARKER), recursive = TRUE, showWarnings = FALSE)
  task_manifest_available <- length(STATE$task_ids) > 0L
  checkpoint_completed <- if (task_manifest_available) {
    STATE$task_ids[vapply(STATE$task_ids, function(id) {
      p <- file.path(CHECKPOINT_DIR, paste0(id, ".rds"))
      if (!file.exists(p)) return(FALSE)
      obj <- tryCatch(readRDS(p), error = function(e) NULL)
      is.list(obj) && identical(obj$task_id, id) &&
        identical(obj$status, "completed") &&
        identical(obj$result_row_count, 3L)
    }, logical(1))]
  } else character()
  completed <- if (task_manifest_available) {
    intersect(
      STATE$task_ids,
      unique(c(STATE$completed_task_ids, checkpoint_completed))
    )
  } else character()
  uncompleted <- if (task_manifest_available) {
    setdiff(STATE$task_ids, completed)
  } else character()
  lines <- c(
    "WP9 SPLIT-ALLOCATION FAILURE - DO NOT USE",
    paste0("run_mode=", RUN_MODE),
    paste0("run_id=", RUN_ID),
    paste0("stage=", STATE$stage),
    paste0("error=", gsub("[\r\n]+", " | ", error_text)),
    paste0("task_manifest_available=", task_manifest_available),
    paste0("completed_task_n=", length(completed)),
    paste0("completed_task_ids=", if (length(completed)) {
      paste(completed, collapse = ";")
    } else if (task_manifest_available) "<none>" else "<not_constructed>"),
    paste0("uncompleted_task_n=", if (task_manifest_available) {
      length(uncompleted)
    } else NA_integer_),
    paste0("uncompleted_task_ids=", if (length(uncompleted)) {
      paste(uncompleted, collapse = ";")
    } else if (task_manifest_available) "<none>" else "<not_constructed>"),
    paste0("time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "reportable_scientific_result=FALSE"
  )
  publish_text(lines, FAILURE_MARKER)
  invisible(TRUE)
}

acquire_formal_lock <- function(script_hash) {
  if (!IS_FORMAL) return(invisible(NULL))
  token <- tryCatch(
    system2(
      "powershell.exe",
      c(
        "-NoProfile", "-NonInteractive", "-Command",
        shQuote("[guid]::NewGuid().ToString('N')")
      ),
      stdout = TRUE, stderr = FALSE
    ),
    error = function(e) character()
  )
  token <- trimws(paste(token, collapse = ""))
  if (length(token) != 1L || !grepl("^[0-9a-fA-F]{32}$", token)) {
    stop("Could not obtain a random OS owner token.", call. = FALSE)
  }
  dir.create(dirname(LOCK_DIR), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(LOCK_DIR) || file.exists(LOCK_DIR)) {
    stop("Formal lock exists; state is not auto-reclaimed: ", LOCK_DIR,
         call. = FALSE)
  }
  if (!dir.create(LOCK_DIR, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not acquire formal exclusive lock.", call. = FALSE)
  }
  acquisition_complete <- FALSE
  on.exit({
    if (!acquisition_complete && dir.exists(LOCK_DIR)) {
      expected <- normalizePath(
        file.path(AUDIT_ROOT, "checkpoints", "wp9",
                  "35_lloyd_formal_exclusive.lock"),
        winslash = "/", mustWork = TRUE
      )
      observed <- normalizePath(LOCK_DIR, winslash = "/", mustWork = TRUE)
      if (!identical(observed, expected) ||
          basename(observed) != "35_lloyd_formal_exclusive.lock") {
        stop("Partial-lock rollback path guard failed.", call. = FALSE)
      }
      contents <- list.files(LOCK_DIR, all.files = TRUE, no.. = TRUE)
      if (length(setdiff(contents, "owner.csv"))) {
        stop("Unexpected content in partial formal lock; refusing rollback.",
             call. = FALSE)
      }
      if (file.exists(LOCK_OWNER)) {
        partial_owner <- data.table::fread(LOCK_OWNER)
        if (nrow(partial_owner) != 1L ||
            !identical(as.character(partial_owner$owner_token[[1L]]), token)) {
          stop("Partial formal-lock token mismatch; refusing rollback.",
               call. = FALSE)
        }
      }
      unlink(LOCK_DIR, recursive = TRUE, force = TRUE)
      if (dir.exists(LOCK_DIR) || file.exists(LOCK_DIR)) {
        stop("Partial formal-lock rollback failed.", call. = FALSE)
      }
    }
  }, add = TRUE)
  owner <- data.table(
    host = Sys.info()[["nodename"]],
    pid = as.integer(Sys.getpid()),
    process_create_time = as.numeric(ps::ps_create_time(ps::ps_handle())),
    owner_token = token,
    run_id = RUN_ID,
    script_sha256 = script_hash,
    spec_sha256 = SPEC_SHA256,
    acquired = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  publish_csv(owner, LOCK_OWNER)
  STATE$lock_acquired <- TRUE
  STATE$lock_ever_acquired <- TRUE
  STATE$lock_token <- token
  acquisition_complete <- TRUE
  invisible(owner)
}

release_formal_lock <- function() {
  if (!IS_FORMAL || !isTRUE(STATE$lock_acquired)) return(invisible(NULL))
  expected <- normalizePath(
    file.path(AUDIT_ROOT, "checkpoints", "wp9",
              "35_lloyd_formal_exclusive.lock"),
    winslash = "/", mustWork = TRUE
  )
  observed <- normalizePath(LOCK_DIR, winslash = "/", mustWork = TRUE)
  if (!identical(observed, expected) || basename(observed) !=
      "35_lloyd_formal_exclusive.lock") {
    stop("Formal lock path failed the exact deletion guard.", call. = FALSE)
  }
  owner <- data.table::fread(LOCK_OWNER)
  if (nrow(owner) != 1L || !identical(
    as.character(owner$owner_token[[1L]]), STATE$lock_token
  )) stop("Formal lock owner token changed; refusing release.", call. = FALSE)
  unlink(LOCK_DIR, recursive = TRUE, force = TRUE)
  if (dir.exists(LOCK_DIR) || file.exists(LOCK_DIR)) {
    stop("Formal lock release failed.", call. = FALSE)
  }
  STATE$lock_acquired <- FALSE
  invisible(NULL)
}

# -------------------------------------------------------------------------
# Output contract and stage evidence
# -------------------------------------------------------------------------

OUTPUT <- list(
  gate = file.path(OUT_TABLE_DIR, "35_WP9_gate_by_split_LLOYD.csv"),
  repeat_level = file.path(OUT_TABLE_DIR, "35_WP9_repeat_level_LLOYD.csv"),
  empirical = file.path(OUT_TABLE_DIR, "35_WP9_empirical_summary_LLOYD.csv"),
  controls = file.path(OUT_TABLE_DIR, "35_WP9_control_summary_LLOYD.csv"),
  states = file.path(OUT_TABLE_DIR, "35_WP9_state_distribution_LLOYD.csv"),
  q95_contrasts = file.path(OUT_TABLE_DIR, "35_WP9_q95_contrasts_LLOYD.csv"),
  transitions = file.path(OUT_TABLE_DIR, "35_WP9_state_transitions_LLOYD.csv"),
  failures = file.path(OUT_LOG_DIR, "35_WP9_failure_log_LLOYD.csv"),
  qc = file.path(OUT_LOG_DIR, "35_WP9_QC_LLOYD.csv"),
  run_log = file.path(OUT_LOG_DIR, "35_WP9_run_log_LLOYD.txt"),
  input_manifest = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_input_sha256_manifest_LLOYD.csv"
  ),
  output_manifest = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_output_sha256_manifest_LLOYD.csv"
  ),
  run_contract = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_run_contract_LLOYD.csv"
  ),
  task_manifest = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_task_manifest_LLOYD.csv"
  ),
  seed_registry = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_global_seed_registry_LLOYD.csv"
  ),
  checkpoint_inventory = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_checkpoint_inventory_LLOYD.csv"
  ),
  session = file.path(
    OUT_PROVENANCE_DIR, "35_WP9_sessionInfo_LLOYD.txt"
  ),
  completion = if (IS_FORMAL) {
    file.path(AUDIT_ROOT, "35_WP9_run_completed_LLOYD.ok")
  } else {
    file.path(RUN_ROOT, "35_WP9_nonformal_run_completed_LLOYD.ok")
  }
)

MANIFEST_KEYS <- setdiff(names(OUTPUT), c("output_manifest", "completion"))
if (length(OUTPUT) != 18L || length(MANIFEST_KEYS) != 16L) {
  stop("Mandatory artifact-count contract mismatch.", call. = FALSE)
}

expected_manifest_paths_for_root <- function(root) {
  c(
    gate = file.path(root, "tables", "35_WP9_gate_by_split_LLOYD.csv"),
    repeat_level = file.path(root, "tables", "35_WP9_repeat_level_LLOYD.csv"),
    empirical = file.path(root, "tables", "35_WP9_empirical_summary_LLOYD.csv"),
    controls = file.path(root, "tables", "35_WP9_control_summary_LLOYD.csv"),
    states = file.path(root, "tables", "35_WP9_state_distribution_LLOYD.csv"),
    q95_contrasts = file.path(root, "tables", "35_WP9_q95_contrasts_LLOYD.csv"),
    transitions = file.path(root, "tables", "35_WP9_state_transitions_LLOYD.csv"),
    failures = file.path(root, "logs", "35_WP9_failure_log_LLOYD.csv"),
    qc = file.path(root, "logs", "35_WP9_QC_LLOYD.csv"),
    run_log = file.path(root, "logs", "35_WP9_run_log_LLOYD.txt"),
    input_manifest = file.path(
      root, "provenance", "35_WP9_input_sha256_manifest_LLOYD.csv"
    ),
    run_contract = file.path(
      root, "provenance", "35_WP9_run_contract_LLOYD.csv"
    ),
    task_manifest = file.path(
      root, "provenance", "35_WP9_task_manifest_LLOYD.csv"
    ),
    seed_registry = file.path(
      root, "provenance", "35_WP9_global_seed_registry_LLOYD.csv"
    ),
    checkpoint_inventory = file.path(
      root, "provenance", "35_WP9_checkpoint_inventory_LLOYD.csv"
    ),
    session = file.path(root, "provenance", "35_WP9_sessionInfo_LLOYD.txt")
  )
}

verify_output_manifest_file <- function(manifest_path, expected_root) {
  x <- fread(manifest_path)
  required <- c("artifact", "path", "bytes", "sha256")
  expected_paths <- expected_manifest_paths_for_root(expected_root)
  expected_artifacts <- names(expected_paths)
  if (!identical(names(x), required) || nrow(x) != 16L ||
      anyDuplicated(x$artifact) || anyDuplicated(x$path) ||
      !identical(as.character(x$artifact), as.character(expected_artifacts))) {
    stop("Output manifest structure/artifact contract failed: ", manifest_path,
         call. = FALSE)
  }
  root <- paste0(normalizePath(
    expected_root, winslash = "/", mustWork = TRUE
  ), "/")
  observed_paths <- normalizePath(x$path, winslash = "/", mustWork = FALSE)
  expected_paths <- normalizePath(
    expected_paths, winslash = "/", mustWork = FALSE
  )
  forbidden_path <- grepl("/(tasks|35_lloyd_formal_exclusive\\.lock)/",
                          observed_paths) |
    grepl("\\.tmp_", observed_paths) |
    grepl("run_completed_LLOYD\\.ok$", observed_paths)
  if (!identical(observed_paths, unname(expected_paths)) ||
      any(!startsWith(observed_paths, root)) || any(forbidden_path) ||
      any(!file.exists(observed_paths))) {
    stop("Output manifest path-scope/existence contract failed: ", manifest_path,
         call. = FALSE)
  }
  observed_hash <- vapply(observed_paths, sha256_file, character(1))
  observed_bytes <- as.numeric(file.info(observed_paths)$size)
  if (!identical(as.character(x$path), observed_paths) ||
      !identical(as.character(x$sha256), unname(observed_hash)) ||
      !isTRUE(all.equal(as.numeric(x$bytes), observed_bytes, tolerance = 0))) {
    stop("Output manifest file hash/size verification failed: ", manifest_path,
         call. = FALSE)
  }
  invisible(x)
}

# -------------------------------------------------------------------------
# Input loading and deterministic task/seed manifests
# -------------------------------------------------------------------------

EXPECTED_HASH <- c(
  v4_lloyd_utils = sha256_file(V4_UTILS),
  strict_wrapper = sha256_file(STRICT_WRAPPER),
  locked_s7_engine = sha256_file(LOCKED_S7_ENGINE),
  v3_loader = sha256_file(V3_LOADER),
  locked_matrix = "9896ad4906207ca910e54eb187507037106bff212299586141908f387cfe6279",
  locked_labels = "159576d0b2796e9786db0eb4c8747af987be5e53af200419e3e37567c1d875f9",
  locked_n1_fit = "a9e6441869fe17dd160380804a1070408ef4036172ac8a316ec81d32f07add6a",
  wp7_seed_registry = "b59177ddc1457dedfd25f7c1a39dac2a4e2d2ee573be080eb15b092793f0174d",
  wp8_seed_registry = "b8137491fc7ad6cd85caceacb208939afad9408bd1cccf1739599debd1c3b511"
)
INPUT_PATH <- c(
  v4_lloyd_utils = V4_UTILS,
  strict_wrapper = STRICT_WRAPPER,
  locked_s7_engine = LOCKED_S7_ENGINE,
  v3_loader = V3_LOADER,
  locked_matrix = LOCKED_MATRIX,
  locked_labels = LOCKED_LABELS,
  locked_n1_fit = LOCKED_N1_FIT,
  wp7_seed_registry = WP7_SEED_REGISTRY,
  wp8_seed_registry = WP8_SEED_REGISTRY
)

build_input_manifest <- function() {
  if (!identical(names(INPUT_PATH), names(EXPECTED_HASH))) {
    stop("Input path/hash maps are not identical.", call. = FALSE)
  }
  missing <- INPUT_PATH[!file.exists(INPUT_PATH)]
  if (length(missing)) {
    stop("Missing governed input(s):\n", paste(missing, collapse = "\n"),
         call. = FALSE)
  }
  dt <- rbindlist(lapply(names(INPUT_PATH), function(nm) {
    p <- INPUT_PATH[[nm]]
    data.table(
      input_name = nm,
      path = normalizePath(p, winslash = "/", mustWork = TRUE),
      sha256 = sha256_file(p),
      expected_sha256 = EXPECTED_HASH[[nm]],
      bytes = as.numeric(file.info(p)$size)
    )
  }))
  dt[, hash_match := sha256 == expected_sha256]
  if (any(!dt$hash_match)) {
    stop("At least one governed input hash changed.", call. = FALSE)
  }
  dt
}

build_tasks <- function() {
  rbindlist(list(
    data.table(stage = "N1_gate_calibration", repeat_id = seq_len(B_GATE)),
    data.table(stage = "N1_independent_control", repeat_id = seq_len(B_EVALUATION)),
    data.table(stage = "S7_positive_control", repeat_id = seq_len(B_EVALUATION)),
    data.table(stage = "locked_empirical", repeat_id = seq_len(B_EVALUATION))
  ))[, `:=`(
    task_order = .I,
    task_id = sprintf("%s_rep%03d", stage, repeat_id),
    run_mode = RUN_MODE,
    run_id = RUN_ID
  )]
}

build_seed_registry <- function(tasks) {
  rows <- list()
  add <- function(module, stage, repeat_id, split_id, role) {
    rows[[length(rows) + 1L]] <<- data.table(
      module = module, stage = stage, repeat_id = as.integer(repeat_id),
      split_id = split_id, seed_role = role
    )
  }
  for (i in seq_len(nrow(tasks))) {
    stage <- tasks$stage[[i]]
    rep <- tasks$repeat_id[[i]]
    if (stage != "locked_empirical") {
      add("source_generation", stage, rep, NA_character_, "generation")
    }
    add("nested_split", stage, rep, NA_character_, "split")
    for (sid in SPLITS$split_id) {
      add("lloyd_k2", stage, rep, sid, "actual_lloyd_seed")
    }
  }
  add("paired_q95_bootstrap", "postprocess", NA_integer_, NA_character_,
      "paired_q95_bootstrap")
  add("future_scheduler", "scheduler", NA_integer_, NA_character_,
      "future_gate")
  add("future_scheduler", "scheduler", NA_integer_, NA_character_,
      "future_evaluation")
  registry <- rbindlist(rows, fill = TRUE)
  registry[, seed := as.integer(SEED_NAMESPACE_START + seq_len(.N) - 1L)]
  registry[, seed_registry_order := .I]
  registry
}

task_seed <- function(registry, stage_value, repeat_value, role_value,
                      split_value = NA_character_) {
  split_match <- if (is.na(split_value)) {
    is.na(registry[["split_id"]])
  } else {
    !is.na(registry[["split_id"]]) &
      registry[["split_id"]] == split_value
  }
  seed_rows <- registry[["stage"]] == stage_value &
    registry[["repeat_id"]] == repeat_value &
    registry[["seed_role"]] == role_value & split_match
  x <- registry[seed_rows, seed]
  if (length(x) != 1L) stop("Seed lookup did not return exactly one seed.",
                            call. = FALSE)
  as.integer(x)
}

module_seed <- function(registry, role) {
  seed_rows <- registry[["seed_role"]] == role
  x <- registry[seed_rows, seed]
  if (length(x) != 1L) stop("Module seed lookup failed: ", role, call. = FALSE)
  as.integer(x)
}

# -------------------------------------------------------------------------
# Scientific generators and paired task execution
# -------------------------------------------------------------------------

simulate_gaussian_copula <- function(fit, n, seed) {
  set.seed(seed)
  p <- ncol(fit$R_used)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*%
    chol(fit$R_used)
  u <- stats::pnorm(z)
  x <- vapply(seq_len(p), function(j) {
    margin <- fit$sorted_margins[[j]]
    idx <- pmax(1L, pmin(length(margin), ceiling(u[, j] * length(margin))))
    margin[idx]
  }, numeric(n))
  colnames(x) <- names(fit$sorted_margins)
  x
}

make_s7_raw <- function(engine, seed) {
  old_delta <- engine$S7_CLUSTER_DELTA
  old_prev <- engine$S7_CLUSTER_PREVALENCE
  on.exit({
    assign("S7_CLUSTER_DELTA", old_delta, envir = engine)
    assign("S7_CLUSTER_PREVALENCE", old_prev, envir = engine)
  }, add = TRUE)
  assign("S7_CLUSTER_DELTA", S7_DELTA, envir = engine)
  assign(
    "S7_CLUSTER_PREVALENCE",
    c(S7_MINOR_PREVALENCE, 1 - S7_MINOR_PREVALENCE),
    envir = engine
  )
  set.seed(seed)
  latent <- engine$simulate_latent_features_s7(N_TOTAL)
  list(
    features = as.data.table(engine$to_clinical_scale(latent$Z)),
    labels = latent$true_subtype
  )
}

make_nested_splits <- function(seed) {
  set.seed(seed)
  perm <- sample.int(N_TOTAL, N_TOTAL, replace = FALSE)
  out <- lapply(seq_len(nrow(SPLITS)), function(j) {
    ne <- SPLITS$n_eval[[j]]
    list(
      meta = SPLITS[j],
      eval = perm[seq_len(ne)],
      train = perm[seq.int(ne + 1L, N_TOTAL)]
    )
  })
  if (!all(out[[1L]]$eval %in% out[[2L]]$eval) ||
      !all(out[[2L]]$eval %in% out[[3L]]$eval) ||
      !all(out[[3L]]$train %in% out[[2L]]$train) ||
      !all(out[[2L]]$train %in% out[[1L]]$train)) {
    stop("Nested evaluation/reverse-training contract failed.", call. = FALSE)
  }
  for (x in out) {
    if (length(x$train) != x$meta$n_train ||
        length(x$eval) != x$meta$n_eval ||
        length(intersect(x$train, x$eval)) != 0L ||
        length(unique(c(x$train, x$eval))) != N_TOTAL) {
      stop("Split complement or coverage contract failed.", call. = FALSE)
    }
  }
  out
}

checkpoint_path <- function(task) file.path(
  CHECKPOINT_DIR, paste0(task$task_id[[1L]], ".rds")
)

checkpoint_signature <- function(task, seeds, script_hash, input_chain_hash) {
  digest::digest(list(
    checkpoint_schema = CHECKPOINT_SCHEMA,
    rule_version = RULE_VERSION,
    script_sha256 = script_hash,
    spec_sha256 = SPEC_SHA256,
    input_chain_sha256 = input_chain_hash,
    run_mode = RUN_MODE, run_id = RUN_ID,
    task_id = task$task_id[[1L]], stage = task$stage[[1L]],
    repeat_id = task$repeat_id[[1L]], seeds = seeds,
    split_ids = SPLITS$split_id
  ), algo = "sha256", serialize = TRUE)
}

read_checkpoint <- function(path, signature, task, script_hash,
                            input_hash_vector, expected_seeds) {
  if (!file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  expected_names <- c(
    "signature", "schema_version", "rule_version", "script_sha256",
    "frozen_spec_sha256", "run_mode", "run_id", "task_id", "stage",
    "repeat_id", "governed_input_hashes", "actual_seeds_used",
    "expected_split_ids", "status", "result_row_count", "result"
  )
  valid <- is.list(obj) && identical(names(obj), expected_names) &&
    identical(obj$signature, signature) &&
    identical(obj$schema_version, CHECKPOINT_SCHEMA) &&
    identical(obj$rule_version, RULE_VERSION) &&
    identical(obj$script_sha256, script_hash) &&
    identical(obj$frozen_spec_sha256, SPEC_SHA256) &&
    identical(obj$run_mode, RUN_MODE) && identical(obj$run_id, RUN_ID) &&
    identical(obj$task_id, task$task_id[[1L]]) &&
    identical(obj$stage, task$stage[[1L]]) &&
    identical(obj$repeat_id, as.integer(task$repeat_id[[1L]])) &&
    identical(obj$governed_input_hashes, input_hash_vector) &&
    is.data.frame(obj$actual_seeds_used) &&
    identical(
      as.data.frame(obj$actual_seeds_used), as.data.frame(expected_seeds)
    ) &&
    identical(obj$expected_split_ids, SPLITS$split_id) &&
    identical(obj$status, "completed") && is.data.frame(obj$result) &&
    identical(obj$result_row_count, 3L) && nrow(obj$result) == 3L &&
    setequal(as.character(obj$result$split_id), SPLITS$split_id)
  if (!valid) {
    stop("Existing checkpoint failed immutable validation: ", path,
         call. = FALSE)
  }
  as.data.table(obj$result)
}

run_one_task <- function(task, registry, script_hash, input_hash_vector,
                         input_chain_hash, n1_fit, engine, X_empirical,
                         labels_empirical) {
  data.table::setDTthreads(1L)
  stage <- task$stage[[1L]]
  rep_id <- task$repeat_id[[1L]]
  seed_rows <- registry[["stage"]] == stage &
    registry[["repeat_id"]] == rep_id
  seeds <- registry[seed_rows, .(
    module, split_id, seed_role, seed
  )]
  signature <- checkpoint_signature(
    task, seeds, script_hash, input_chain_hash
  )
  cp <- checkpoint_path(task)
  if (file.exists(cp)) {
    if (IS_FORMAL) stop("Formal checkpoint root was not empty.", call. = FALSE)
    cached <- read_checkpoint(
      cp, signature, task, script_hash, input_hash_vector, seeds
    )
    return(list(status = "completed", result = cached, error = NA_character_))
  }

  result <- tryCatch(withCallingHandlers({
    split_seed <- task_seed(registry, stage, rep_id, "split")
    split_list <- make_nested_splits(split_seed)
    generation_seed <- if (stage != "locked_empirical") {
      task_seed(registry, stage, rep_id, "generation")
    } else NA_integer_

    source <- if (stage %in% c(
      "N1_gate_calibration", "N1_independent_control"
    )) {
      simulate_gaussian_copula(n1_fit, N_TOTAL, generation_seed)
    } else if (stage == "S7_positive_control") {
      make_s7_raw(engine, generation_seed)
    } else if (stage == "locked_empirical") {
      X_empirical
    } else stop("Unknown task stage.", call. = FALSE)

    rows <- lapply(split_list, function(sp) {
      sid <- sp$meta$split_id[[1L]]
      actual_lloyd_seed <- task_seed(
        registry, stage, rep_id, "actual_lloyd_seed", sid
      )
      utility_seed_argument <- as.integer(actual_lloyd_seed - 2L)
      if (utility_seed_argument + 2L != actual_lloyd_seed) {
        stop("Frozen-v4 seed adapter equality failed.", call. = FALSE)
      }

      if (stage == "S7_positive_control") {
        prep <- preprocess_train_evaluation(
          engine, source$features[sp$train], source$features[sp$eval]
        )
        x_train <- prep$X_train
        x_eval <- prep$X_evaluation
        eval_labels <- source$labels[sp$eval]
        train_minor <- mean(source$labels[sp$train] == "S7_cluster_A")
        eval_minor <- mean(source$labels[sp$eval] == "S7_cluster_A")
      } else {
        x_train <- source[sp$train, , drop = FALSE]
        x_eval <- source[sp$eval, , drop = FALSE]
        eval_labels <- if (stage == "locked_empirical") {
          labels_empirical[sp$eval]
        } else NULL
        train_minor <- NA_real_
        eval_minor <- NA_real_
      }

      diag <- domain1_external_diagnostics(
        X_train = x_train, X_evaluation = x_eval,
        evaluation_labels = eval_labels, seed = utility_seed_argument
      )
      if (nrow(diag) != 1L || diag$kmeans_primary_algorithm_used[[1L]] !=
          "Lloyd" || isTRUE(diag$kmeans_primary_fallback_used[[1L]]) ||
          diag$kmeans_primary_warning_count[[1L]] != 0L ||
          diag$kmeans_primary_nstart[[1L]] != KMEANS_NSTART ||
          diag$kmeans_primary_itermax[[1L]] != KMEANS_ITERMAX) {
        stop("Strict Lloyd diagnostic metadata failed.", call. = FALSE)
      }
      cbind(data.table(
        stage = stage, repeat_id = as.integer(rep_id), task_id = task$task_id[[1L]],
        split_order = sp$meta$split_order[[1L]], split_id = sid,
        split_label = sp$meta$split_label[[1L]],
        n_train = sp$meta$n_train[[1L]], n_eval = sp$meta$n_eval[[1L]],
        generation_seed = generation_seed, split_seed = split_seed,
        registered_actual_lloyd_seed = actual_lloyd_seed,
        utility_seed_argument = utility_seed_argument,
        seed_adapter_match = utility_seed_argument + 2L == actual_lloyd_seed,
        shape_alert = diag$dip_discriminant_p_fdr[[1L]] < DIP_ALPHA,
        s7_train_minor_prevalence = train_minor,
        s7_eval_minor_prevalence = eval_minor,
        status = "completed", task_warning_count = 0L,
        run_mode = RUN_MODE, run_id = RUN_ID,
        rule_version = RULE_VERSION
      ), diag)
    })
    rbindlist(rows, fill = TRUE)
  }, warning = function(w) {
    stop("Task emitted warning under strict contract: ", conditionMessage(w),
         call. = FALSE)
  }), error = function(e) e)

  if (inherits(result, "error")) {
    return(list(
      status = "failed", result = NULL,
      error = conditionMessage(result), task_id = task$task_id[[1L]],
      stage = stage, repeat_id = rep_id
    ))
  }
  envelope <- list(
    signature = signature,
    schema_version = CHECKPOINT_SCHEMA,
    rule_version = RULE_VERSION,
    script_sha256 = script_hash,
    frozen_spec_sha256 = SPEC_SHA256,
    run_mode = RUN_MODE, run_id = RUN_ID,
    task_id = task$task_id[[1L]], stage = stage,
    repeat_id = as.integer(rep_id), governed_input_hashes = input_hash_vector,
    actual_seeds_used = seeds, expected_split_ids = SPLITS$split_id,
    status = "completed", result_row_count = 3L, result = result
  )
  publish_rds(envelope, cp)
  list(status = "completed", result = result, error = NA_character_)
}

run_task_block <- function(tasks, scheduler_seed, registry, script_hash,
                           input_hash_vector, input_chain_hash, n1_fit, engine,
                           X_empirical, labels_empirical, workers) {
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  out <- future.apply::future_lapply(
    seq_len(nrow(tasks)),
    function(i) run_one_task(
      tasks[i], registry, script_hash, input_hash_vector, input_chain_hash,
      n1_fit, engine, X_empirical, labels_empirical
    ),
    future.seed = scheduler_seed,
    future.scheduling = if (RUN_MODE == "test") 1 else 4,
    future.packages = c(
      "data.table", "diptest", "Matrix", "mclust", "digest"
    )
  )
  future::plan(future::sequential)
  out
}

# -------------------------------------------------------------------------
# Summaries
# -------------------------------------------------------------------------

summarise_source_split <- function(dt) {
  rbindlist(lapply(split(dt, by = c("stage", "split_id"), keep.by = TRUE),
    function(z) {
      base <- unique(z[, .(
        stage, split_order, split_id, split_label, n_train, n_eval
      )])
      if (nrow(base) != 1L) stop("Summary grouping metadata mismatch.")
      cbind(
        base,
        rate_summary(sum(z$separation_alert), nrow(z), "separation_alert"),
        rate_summary(sum(z$shape_alert), nrow(z), "shape_alert"),
        numeric_summary(z$mahalanobis_centroid_delta, "separation"),
        numeric_summary(z$dip_discriminant_p_fdr, "fisher_p_fdr"),
        numeric_summary(z$refit_vs_input_label_ari_descriptive, "ari"),
        numeric_summary(z$s7_train_minor_prevalence, "s7_train_minor_prev"),
        numeric_summary(z$s7_eval_minor_prevalence, "s7_eval_minor_prev")
      )
    }
  ), fill = TRUE)
}

make_state_distribution <- function(dt) {
  out <- dt[, .(count = .N), by = .(
    stage, split_id, state
  )]
  grid <- CJ(
    stage = unique(dt$stage), split_id = SPLITS$split_id,
    state = c("DISCRETE_EVIDENCE", "NO_DISCRETE_EVIDENCE", "INCONCLUSIVE"),
    unique = TRUE
  )
  out <- merge(grid, out, by = c("stage", "split_id", "state"), all.x = TRUE)
  out[is.na(count), count := 0L]
  out <- merge(out, SPLITS, by = "split_id", all.x = TRUE)
  out[, n := B_EVALUATION]
  out[, rate := count / n]
  cis <- t(vapply(seq_len(nrow(out)), function(i) {
    wilson_interval(out$count[[i]], out$n[[i]])
  }, numeric(2)))
  out[, `:=`(ci_low = cis[, 1L], ci_high = cis[, 2L])]
  setorder(out, stage, split_order, state)
  out[, .(
    stage, split_order, split_id, split_label, n_train, n_eval,
    state, count, n, rate, ci_low, ci_high
  )]
}

make_transitions <- function(dt) {
  pairs <- data.table(
    earlier = SPLITS$split_id[c(1L, 2L, 1L)],
    later = SPLITS$split_id[c(2L, 3L, 3L)]
  )
  rows <- list()
  for (st in unique(dt$stage)) for (i in seq_len(nrow(pairs))) {
    a <- dt[stage == st & split_id == pairs$earlier[[i]], .(
      repeat_id, row_state = state
    )]
    b <- dt[stage == st & split_id == pairs$later[[i]], .(
      repeat_id, column_state = state
    )]
    x <- merge(a, b, by = "repeat_id")
    tab <- x[, .(count = .N), by = .(row_state, column_state)]
    tab[, `:=`(
      stage = st, earlier_split = pairs$earlier[[i]],
      later_split = pairs$later[[i]], n_matched = nrow(x)
    )]
    rows[[length(rows) + 1L]] <- tab
  }
  rbindlist(rows, fill = TRUE)
}

checkpoint_inventory <- function(tasks, script_hash, input_hash_vector,
                                 input_chain_hash, registry) {
  rbindlist(lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i]
    stage <- task$stage[[1L]]
    rep_id <- task$repeat_id[[1L]]
    seed_rows <- registry[["stage"]] == stage &
      registry[["repeat_id"]] == rep_id
    seeds <- registry[seed_rows, .(
      module, split_id, seed_role, seed
    )]
    sig <- checkpoint_signature(task, seeds, script_hash, input_chain_hash)
    p <- checkpoint_path(task)
    obj <- tryCatch(readRDS(p), error = function(e) NULL)
    result <- if (is.list(obj) && is.data.frame(obj$result)) {
      as.data.table(obj$result)
    } else NULL
    expected_names <- c(
      "signature", "schema_version", "rule_version", "script_sha256",
      "frozen_spec_sha256", "run_mode", "run_id", "task_id", "stage",
      "repeat_id", "governed_input_hashes", "actual_seeds_used",
      "expected_split_ids", "status", "result_row_count", "result"
    )
    data.table(
      path = normalizePath(p, winslash = "/", mustWork = FALSE),
      bytes = if (file.exists(p)) as.numeric(file.info(p)$size) else NA_real_,
      sha256 = if (file.exists(p)) sha256_file(p) else NA_character_,
      readable = !is.null(obj), task_id = task$task_id[[1L]], stage = stage,
      repeat_id = as.integer(rep_id),
      schema_match = is.list(obj) && identical(names(obj), expected_names) &&
        identical(obj$signature, sig) &&
        identical(obj$schema_version, CHECKPOINT_SCHEMA) &&
        identical(obj$task_id, task$task_id[[1L]]) &&
        identical(obj$stage, stage) &&
        identical(obj$repeat_id, as.integer(rep_id)) &&
        is.data.frame(obj$actual_seeds_used) &&
        identical(
          as.data.frame(obj$actual_seeds_used), as.data.frame(seeds)
        ) &&
        identical(obj$expected_split_ids, SPLITS$split_id) &&
        identical(obj$result_row_count, 3L),
      rule_match = is.list(obj) && identical(obj$rule_version, RULE_VERSION),
      run_id_match = is.list(obj) && identical(obj$run_id, RUN_ID),
      run_mode_match = is.list(obj) && identical(obj$run_mode, RUN_MODE),
      script_hash_match = is.list(obj) && identical(obj$script_sha256, script_hash),
      spec_hash_match = is.list(obj) &&
        identical(obj$frozen_spec_sha256, SPEC_SHA256),
      input_hash_match = is.list(obj) &&
        identical(obj$governed_input_hashes, input_hash_vector),
      status = if (is.list(obj)) obj$status %||% NA_character_ else NA_character_,
      result_nrow = if (is.null(result)) NA_integer_ else nrow(result),
      split_id_match = !is.null(result) && nrow(result) == 3L &&
        !anyDuplicated(result$split_id) &&
        setequal(as.character(result$split_id), SPLITS$split_id)
    )
  }), fill = TRUE)
}

qc_row <- function(check, pass, observed, expected) data.table(
  check = check, pass = isTRUE(pass), observed = as.character(observed),
  expected = as.character(expected)
)

# -------------------------------------------------------------------------
# Main governed run
# -------------------------------------------------------------------------

wp9_main <- function() {
  STATE$stage <- "preflight"
  script_hash <- sha256_file(SCRIPT_PATH)
  prior_evidence <- NULL

  if (IS_FORMAL) {
    if (dir.exists(RUN_ROOT) && length(list.files(
      RUN_ROOT, all.files = TRUE, no.. = TRUE, recursive = TRUE
    ))) stop("Formal run root is not empty.", call. = FALSE)
    collisions <- unlist(OUTPUT[c(MANIFEST_KEYS, "output_manifest", "completion")])
    collisions <- collisions[file.exists(collisions) | dir.exists(collisions)]
    if (length(collisions)) {
      stop("Formal target collision(s):\n", paste(collisions, collapse = "\n"),
           call. = FALSE)
    }
  } else if (dir.exists(RUN_ROOT)) {
    blockers <- c(FAILURE_MARKER, OUTPUT$completion)
    if (any(file.exists(blockers))) {
      stop("Refusing reuse of completed or failed nonformal run root.",
           call. = FALSE)
    }
  }

  # The formal exclusive lock is the first formal filesystem mutation.
  # A competing process that cannot acquire it must not touch RUN_ROOT.
  acquire_formal_lock(script_hash)

  for (d in c(RUN_ROOT, OUT_TABLE_DIR, OUT_LOG_DIR, OUT_PROVENANCE_DIR,
              CHECKPOINT_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(OUTPUT$run_log)) stop("Run log target already exists.")
  STATE$log_enabled <- TRUE
  log_msg("WP9 start; mode=", RUN_MODE, "; run_id=", RUN_ID)

  STATE$stage <- "input_verification"
  input_manifest <- build_input_manifest()
  input_hash_vector <- stats::setNames(
    input_manifest$sha256, input_manifest$input_name
  )
  input_chain_hash <- digest::digest(
    input_hash_vector, algo = "sha256", serialize = TRUE
  )

  source(V4_UTILS, local = .GlobalEnv)
  source(STRICT_WRAPPER, local = .GlobalEnv)
  source(V3_LOADER, local = .GlobalEnv)
  engine <- load_locked_main_s7_engine(LOCKED_S7_ENGINE)
  n1_fit <- readRDS(LOCKED_N1_FIT)

  x_dt <- as.data.table(readRDS(LOCKED_MATRIX))
  y_dt <- as.data.table(readRDS(LOCKED_LABELS))
  if (!all(c("stay_id") %in% names(x_dt)) ||
      !all(c("stay_id", "cluster_k2") %in% names(y_dt)) ||
      nrow(x_dt) != N_TOTAL || nrow(y_dt) != N_TOTAL ||
      anyNA(x_dt$stay_id) || anyNA(y_dt$stay_id) ||
      anyDuplicated(x_dt$stay_id) || anyDuplicated(y_dt$stay_id) ||
      !setequal(x_dt$stay_id, y_dt$stay_id)) {
    stop("Locked empirical stay_id integrity contract failed.", call. = FALSE)
  }
  y_dt <- y_dt[match(x_dt$stay_id, y_dt$stay_id)]
  if (!identical(x_dt$stay_id, y_dt$stay_id)) {
    stop("Locked empirical post-match row order failed.", call. = FALSE)
  }
  feature_names <- setdiff(names(x_dt), "stay_id")
  X_empirical <- as.matrix(x_dt[, ..feature_names])
  labels_empirical <- as.integer(y_dt$cluster_k2)
  if (ncol(X_empirical) != 33L || !is.numeric(X_empirical) ||
      any(!is.finite(X_empirical)) ||
      !identical(as.integer(table(labels_empirical)), c(3992L, 16057L))) {
    stop("Locked empirical feature/label contract failed.", call. = FALSE)
  }
  if (length(engine$feature_names) != 33L ||
      length(intersect(feature_names, engine$feature_names)) != 23L) {
    stop("Empirical/S7 feature-space boundary contract failed.", call. = FALSE)
  }

  STATE$stage <- "manifest_construction"
  tasks <- build_tasks()
  STATE$task_ids <- as.character(tasks$task_id)
  registry <- build_seed_registry(tasks)
  old7 <- fread(WP7_SEED_REGISTRY)$seed
  old8 <- fread(WP8_SEED_REGISTRY)$seed
  test_seed <- smoke_seed <- integer()
  seed_pass <- nrow(registry) == EXPECTED_SEEDS &&
    all(is.finite(registry$seed)) && all(registry$seed > 0L) &&
    all(registry$seed <= .Machine$integer.max) &&
    !anyDuplicated(registry$seed) &&
    nrow(registry[, uniqueN(module), by = seed][V1 > 1L]) == 0L &&
    length(intersect(registry$seed, old7)) == 0L &&
    length(intersect(registry$seed, old8)) == 0L &&
    length(intersect(registry$seed, test_seed)) == 0L &&
    length(intersect(registry$seed, smoke_seed)) == 0L
  if (!seed_pass) stop("Global actual-seed registry hard QC failed.")

  workers_default <- max(1L, min(6L, parallel::detectCores(logical = TRUE) - 2L))
  workers <- as.integer(Sys.getenv(
    "D1_WP9_WORKERS", unset = as.character(workers_default)
  ))
  if (!is.finite(workers) || workers < 1L) stop("Invalid worker count.")
  if (IS_FORMAL && workers < 2L) stop("Formal WP9 requires >=2 workers.")

  run_contract <- data.table(
    run_mode = RUN_MODE, run_id = RUN_ID, n_total = N_TOTAL,
    B_gate = B_GATE, B_evaluation = B_EVALUATION, B_boot = B_BOOT,
    expected_tasks = EXPECTED_TASKS, expected_repeat_rows = EXPECTED_REPEAT_ROWS,
    expected_seeds = EXPECTED_SEEDS, workers = workers,
    kmeans_algorithm = "Lloyd", kmeans_nstart = KMEANS_NSTART,
    kmeans_itermax = KMEANS_ITERMAX, fallback_allowed = FALSE,
    ridge_multiplier = RIDGE_MULTIPLIER, dip_alpha = DIP_ALPHA,
    q95_probability = GATE_QUANTILE, quantile_type = QUANTILE_TYPE,
    s7_delta = S7_DELTA, s7_minor_prevalence = S7_MINOR_PREVALENCE,
    script_sha256 = script_hash, spec_sha256 = SPEC_SHA256,
    input_chain_sha256 = input_chain_hash
  )

  STATE$stage <- "gate_tasks"
  gate_tasks <- tasks[stage == "N1_gate_calibration"]
  gate_out <- run_task_block(
    gate_tasks, module_seed(registry, "future_gate"), registry, script_hash,
    input_hash_vector, input_chain_hash, n1_fit, engine, X_empirical,
    labels_empirical, workers
  )
  STATE$completed_task_ids <- unique(c(
    STATE$completed_task_ids,
    gate_tasks$task_id[vapply(
      gate_out, function(x) identical(x$status, "completed"), logical(1)
    )]
  ))
  gate_fail <- rbindlist(lapply(gate_out, function(x) if (x$status == "failed") {
    data.table(task_id = x$task_id, stage = x$stage, repeat_id = x$repeat_id,
               error = x$error)
  } else NULL), fill = TRUE)
  if (nrow(gate_fail)) {
    publish_csv(gate_fail, OUTPUT$failures)
    stop("Gate task failure(s); evaluation was not run.", call. = FALSE)
  }
  gate_rows <- rbindlist(lapply(gate_out, `[[`, "result"), fill = TRUE)
  if (nrow(gate_rows) != 3L * B_GATE) stop("Gate row-count mismatch.")

  STATE$stage <- "paired_q95_bootstrap"
  boot_seed <- module_seed(registry, "paired_q95_bootstrap")
  gate_by_split <- stats::setNames(lapply(SPLITS$split_id, function(sid) {
    x <- copy(gate_rows[split_id == sid])
    setorder(x, repeat_id)
    if (!identical(x$repeat_id, seq_len(B_GATE))) {
      stop("Gate repeats are not complete/aligned by repeat_id: ", sid,
           call. = FALSE)
    }
    x
  }), SPLITS$split_id)
  set.seed(boot_seed)
  boot_index <- replicate(
    B_BOOT, sample.int(B_GATE, B_GATE, replace = TRUE), simplify = FALSE
  )
  q95_boot <- rbindlist(lapply(seq_len(B_BOOT), function(b) {
    idx <- boot_index[[b]]
    rbindlist(lapply(SPLITS$split_id, function(sid) data.table(
       bootstrap_id = b, split_id = sid,
       q95 = safe_q(
         gate_by_split[[sid]]$mahalanobis_centroid_delta[idx],
         GATE_QUANTILE
       )
    )))
  }))
  gate_table <- gate_rows[, .(
    n_gate = .N,
    separation_mean = mean(mahalanobis_centroid_delta),
    separation_sd = sd(mahalanobis_centroid_delta),
    separation_mcse = sd(mahalanobis_centroid_delta) / sqrt(.N),
    separation_p025 = safe_q(mahalanobis_centroid_delta, 0.025),
    separation_p975 = safe_q(mahalanobis_centroid_delta, 0.975),
    q95 = safe_q(mahalanobis_centroid_delta, GATE_QUANTILE),
    fisher_axis_p_fdr_mean = mean(dip_discriminant_p_fdr),
    fisher_axis_p_fdr_sd = sd(dip_discriminant_p_fdr),
    fisher_axis_p_fdr_mcse = sd(dip_discriminant_p_fdr) / sqrt(.N),
    fisher_axis_p_fdr_p025 = safe_q(dip_discriminant_p_fdr, 0.025),
    fisher_axis_p_fdr_p975 = safe_q(dip_discriminant_p_fdr, 0.975)
  ), by = .(split_order, split_id, split_label, n_train, n_eval)]
  boot_summary <- q95_boot[, .(
    q95_bootstrap_B = .N,
    q95_bootstrap_seed = boot_seed,
    q95_bootstrap_se = sd(q95),
    q95_bootstrap_ci_low = safe_q(q95, 0.025),
    q95_bootstrap_ci_high = safe_q(q95, 0.975)
  ), by = split_id]
  gate_table <- merge(gate_table, boot_summary, by = "split_id")
  setorder(gate_table, split_order)

  contrast_pairs <- data.table(
    contrast_order = 1:3,
    earlier = SPLITS$split_id[c(1L, 2L, 1L)],
    later = SPLITS$split_id[c(2L, 3L, 3L)]
  )
  contrast_rows <- rbindlist(lapply(seq_len(nrow(contrast_pairs)), function(i) {
    early <- contrast_pairs$earlier[[i]]
    late <- contrast_pairs$later[[i]]
    point <- gate_table[split_id == late, q95] -
      gate_table[split_id == early, q95]
    wide <- dcast(
      q95_boot[split_id %in% c(early, late)], bootstrap_id ~ split_id,
      value.var = "q95"
    )
    delta <- wide[[late]] - wide[[early]]
    data.table(
      contrast_order = i, earlier_split = early, later_split = late,
      contrast = paste0(late, " minus ", early),
      q95_difference = point, bootstrap_B = length(delta),
      bootstrap_seed = boot_seed, bootstrap_se = sd(delta),
      ci_low = safe_q(delta, 0.025), ci_high = safe_q(delta, 0.975)
    )
  }))

  STATE$stage <- "combined_evaluation_tasks"
  eval_tasks <- tasks[stage != "N1_gate_calibration"]
  eval_out <- run_task_block(
    eval_tasks, module_seed(registry, "future_evaluation"), registry,
    script_hash, input_hash_vector, input_chain_hash, n1_fit, engine,
    X_empirical, labels_empirical, workers
  )
  STATE$completed_task_ids <- unique(c(
    STATE$completed_task_ids,
    eval_tasks$task_id[vapply(
      eval_out, function(x) identical(x$status, "completed"), logical(1)
    )]
  ))
  eval_fail <- rbindlist(lapply(eval_out, function(x) if (x$status == "failed") {
    data.table(task_id = x$task_id, stage = x$stage, repeat_id = x$repeat_id,
               error = x$error)
  } else NULL), fill = TRUE)
  if (nrow(eval_fail)) {
    publish_csv(eval_fail, OUTPUT$failures)
    stop("Evaluation task failure(s).", call. = FALSE)
  }
  eval_rows <- rbindlist(lapply(eval_out, `[[`, "result"), fill = TRUE)
  if (nrow(eval_rows) != 9L * B_EVALUATION) {
    stop("Evaluation row-count mismatch.")
  }

  STATE$stage <- "three_state_decision"
  eval_rows <- merge(
    eval_rows,
    gate_table[, .(split_id, separation_gate = q95)],
    by = "split_id", all.x = TRUE
  )
  decisions <- lapply(seq_len(nrow(eval_rows)), function(i) {
    x <- eval_rows[i]
    audit_discreteness_v2(list(
      separation = x$mahalanobis_centroid_delta[[1L]],
      gate = x$separation_gate[[1L]],
      shape_reject = x$shape_alert[[1L]],
      shape_source = x$shape_source[[1L]],
      reference_scope = x$reference_scope[[1L]],
      partition_source = x$partition_source[[1L]],
      kmeans_algorithm = x$kmeans_algorithm[[1L]],
      kmeans_nstart = x$kmeans_nstart[[1L]],
      kmeans_itermax = x$kmeans_itermax[[1L]]
    ))
  })
  if (any(vapply(decisions, function(x) x$state == "NOT_EVALUATED", logical(1)))) {
    stop("Strict wrapper returned NOT_EVALUATED.")
  }
  eval_rows[, `:=`(
    separation_alert = vapply(decisions, `[[`, logical(1), "separation_alert"),
    wrapper_shape_alert = vapply(decisions, `[[`, logical(1), "shape_alert"),
    state = vapply(decisions, `[[`, character(1), "state"),
    decision_rule_version = vapply(
      decisions, `[[`, character(1), "decision_rule_version"
    )
  )]
  if (!identical(eval_rows$shape_alert, eval_rows$wrapper_shape_alert)) {
    stop("Shape-alert handoff mismatch.")
  }
  repeat_rows <- rbindlist(list(gate_rows, eval_rows), fill = TRUE)

  STATE$stage <- "summary_and_qc"
  empirical_summary <- summarise_source_split(
    eval_rows[stage == "locked_empirical"]
  )
  control_summary <- summarise_source_split(
    eval_rows[stage %in% c("N1_independent_control", "S7_positive_control")]
  )
  state_distribution <- make_state_distribution(eval_rows)
  transitions <- make_transitions(eval_rows)
  failure_log <- data.table(
    task_id = character(), stage = character(), repeat_id = integer(),
    error = character()
  )
  cp_inventory <- checkpoint_inventory(
    tasks, script_hash, input_hash_vector, input_chain_hash, registry
  )
  expected_stage_tasks <- data.table(
    stage = c(
      "N1_gate_calibration", "N1_independent_control",
      "S7_positive_control", "locked_empirical"
    ),
    expected_n = c(B_GATE, B_EVALUATION, B_EVALUATION, B_EVALUATION)
  )
  observed_stage_tasks <- tasks[, .(observed_n = .N), by = stage]
  stage_task_audit <- merge(
    expected_stage_tasks, observed_stage_tasks, by = "stage", all = TRUE
  )
  stage_task_pass <- nrow(stage_task_audit) == 4L &&
    all(stage_task_audit$observed_n == stage_task_audit$expected_n)
  diagnostic_columns <- c(
    "mahalanobis_centroid_delta", "dip_discriminant_p_fdr"
  )

  qc <- rbindlist(list(
    qc_row("task_count", nrow(tasks) == EXPECTED_TASKS, nrow(tasks), EXPECTED_TASKS),
    qc_row("unique_stage_repeat", uniqueN(tasks, by = c("stage", "repeat_id")) ==
             EXPECTED_TASKS, uniqueN(tasks, by = c("stage", "repeat_id")), EXPECTED_TASKS),
    qc_row("stage_task_counts", stage_task_pass,
           paste(stage_task_audit$stage, stage_task_audit$observed_n,
                 collapse = ";"),
           paste(expected_stage_tasks$stage, expected_stage_tasks$expected_n,
                 collapse = ";")),
    qc_row("repeat_row_count", nrow(repeat_rows) ==
             EXPECTED_REPEAT_ROWS, nrow(repeat_rows),
           EXPECTED_REPEAT_ROWS),
    qc_row("three_rows_and_splits_per_task",
      nrow(repeat_rows[, .(
        n_row = .N, n_split = uniqueN(split_id)
      ), by = .(stage, repeat_id)][n_row != 3L | n_split != 3L]) == 0L,
      "checked", "3 rows / 3 splits per task"),
    qc_row("checkpoint_count", nrow(cp_inventory) == EXPECTED_TASKS,
           nrow(cp_inventory), EXPECTED_TASKS),
    qc_row("checkpoint_all_valid", all(cp_inventory[, readable &
      schema_match & rule_match & run_id_match & run_mode_match &
      script_hash_match & spec_hash_match & input_hash_match &
      status == "completed" & result_nrow == 3L & split_id_match]),
      sum(cp_inventory$readable), EXPECTED_TASKS),
    qc_row("seed_registry_rows", nrow(registry) == EXPECTED_SEEDS,
           nrow(registry), EXPECTED_SEEDS),
    qc_row("seed_registry_contract", seed_pass, seed_pass, TRUE),
    qc_row("seed_namespace_exact_consecutive_range",
      identical(
        registry$seed,
        as.integer(SEED_NAMESPACE_START + seq_len(EXPECTED_SEEDS) - 1L)
      ),
      paste0(min(registry$seed), "-", max(registry$seed)),
      paste0(SEED_NAMESPACE_START, "-",
             SEED_NAMESPACE_START + EXPECTED_SEEDS - 1L)),
    qc_row("failed_tasks", nrow(failure_log) == 0L, nrow(failure_log), 0L),
    qc_row("fallback_count", sum(eval_rows$kmeans_primary_fallback_used %in% TRUE) +
             sum(gate_rows$kmeans_primary_fallback_used %in% TRUE) == 0L,
           sum(eval_rows$kmeans_primary_fallback_used %in% TRUE) +
             sum(gate_rows$kmeans_primary_fallback_used %in% TRUE), 0L),
    qc_row("warning_count", sum(eval_rows$task_warning_count) +
             sum(gate_rows$task_warning_count) +
             sum(eval_rows$kmeans_primary_warning_count) +
             sum(gate_rows$kmeans_primary_warning_count) == 0L,
           sum(eval_rows$task_warning_count) + sum(gate_rows$task_warning_count) +
             sum(eval_rows$kmeans_primary_warning_count) +
             sum(gate_rows$kmeans_primary_warning_count), 0L),
    qc_row("strict_lloyd_metadata", all(c(gate_rows$kmeans_algorithm,
      eval_rows$kmeans_algorithm) == "Lloyd") &&
      all(c(gate_rows$kmeans_nstart, eval_rows$kmeans_nstart) == 25L) &&
      all(c(gate_rows$kmeans_itermax, eval_rows$kmeans_itermax) == 100L),
      "checked", "Lloyd/25/100"),
    qc_row("diagnostic_completeness",
      all(diagnostic_columns %in% names(repeat_rows)) &&
        all(vapply(diagnostic_columns, function(nm) {
        all(is.finite(repeat_rows[[nm]]))
      }, logical(1))) &&
        all(repeat_rows$dip_discriminant_p_fdr >= 0 &
              repeat_rows$dip_discriminant_p_fdr <= 1),
      "checked", "finite separation and adjusted Fisher-axis P in [0,1]"),
    qc_row("seed_adapter_rows", all(c(gate_rows$seed_adapter_match,
      eval_rows$seed_adapter_match)),
      sum(c(gate_rows$seed_adapter_match, eval_rows$seed_adapter_match)),
      EXPECTED_REPEAT_ROWS),
    qc_row("q95_finite_ordered", nrow(gate_table) == 3L &&
      all(is.finite(unlist(gate_table[, .(
        q95, q95_bootstrap_se, q95_bootstrap_ci_low,
        q95_bootstrap_ci_high, fisher_axis_p_fdr_mean,
        fisher_axis_p_fdr_sd, fisher_axis_p_fdr_mcse,
        fisher_axis_p_fdr_p025, fisher_axis_p_fdr_p975
      )]))) && all(gate_table$q95_bootstrap_ci_low <=
                     gate_table$q95_bootstrap_ci_high) &&
        all(gate_table$fisher_axis_p_fdr_p025 <=
              gate_table$fisher_axis_p_fdr_p975), "checked", TRUE),
    qc_row("q95_contrasts", nrow(contrast_rows) == 3L &&
      all(is.finite(unlist(contrast_rows[, .(
        q95_difference, bootstrap_se, ci_low, ci_high
      )]))) && all(contrast_rows$ci_low <= contrast_rows$ci_high),
      nrow(contrast_rows), 3L),
    qc_row("state_cell_totals", all(state_distribution[, sum(count),
      by = .(stage, split_id)]$V1 == B_EVALUATION), "checked", B_EVALUATION),
    qc_row("state_rate_totals", all(abs(state_distribution[, sum(rate),
      by = .(stage, split_id)]$V1 - 1) < 1e-12), "checked", 1),
    qc_row("wilson_bounds", all(state_distribution$ci_low >= 0 &
      state_distribution$ci_high <= 1), "checked", "[0,1]"),
    qc_row("transition_totals", all(transitions$count >= 0) &&
      all(transitions[, sum(count), by = .(
        stage, earlier_split, later_split
      )]$V1 == B_EVALUATION), "checked", B_EVALUATION),
    qc_row("ari_contract", all(is.na(eval_rows[
      stage == "N1_independent_control",
      refit_vs_input_label_ari_descriptive
    ])) && all(is.na(gate_rows$refit_vs_input_label_ari_descriptive)) &&
      all(is.finite(eval_rows[
      stage %in% c("S7_positive_control", "locked_empirical"),
      refit_vs_input_label_ari_descriptive
    ])) && all(eval_rows[
      stage %in% c("S7_positive_control", "locked_empirical"),
      refit_vs_input_label_ari_descriptive
    ] >= -1 & eval_rows[
      stage %in% c("S7_positive_control", "locked_empirical"),
      refit_vs_input_label_ari_descriptive
    ] <= 1), "checked", TRUE)
  ))
  if (any(!qc$pass)) stop("WP9 hard QC failed before output publication.")

  STATE$stage <- "output_publication"
  start_hash <- build_input_manifest()
  if (!identical(input_manifest$sha256, start_hash$sha256)) {
    stop("Governed input hash changed during run.")
  }
  publish_csv(gate_table, OUTPUT$gate)
  publish_csv(repeat_rows, OUTPUT$repeat_level)
  publish_csv(empirical_summary, OUTPUT$empirical)
  publish_csv(control_summary, OUTPUT$controls)
  publish_csv(state_distribution, OUTPUT$states)
  publish_csv(contrast_rows, OUTPUT$q95_contrasts)
  publish_csv(transitions, OUTPUT$transitions)
  publish_csv(failure_log, OUTPUT$failures)
  publish_csv(qc, OUTPUT$qc)
  publish_csv(input_manifest, OUTPUT$input_manifest)
  publish_csv(run_contract, OUTPUT$run_contract)
  publish_csv(tasks, OUTPUT$task_manifest)
  publish_csv(registry, OUTPUT$seed_registry)
  publish_csv(cp_inventory, OUTPUT$checkpoint_inventory)
  publish_text(capture.output(sessionInfo()), OUTPUT$session)

  log_msg("All scientific outputs and QC published; closing run log.")
  STATE$log_enabled <- FALSE

  STATE$stage <- "output_manifest"
  manifest_paths <- unlist(OUTPUT[MANIFEST_KEYS], use.names = TRUE)
  if (length(manifest_paths) != 16L || any(!file.exists(manifest_paths))) {
    stop("Explicit 16-file output-manifest contract failed.")
  }
  output_manifest <- data.table(
    artifact = names(manifest_paths),
    path = normalizePath(manifest_paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(manifest_paths)$size),
    sha256 = vapply(manifest_paths, sha256_file, character(1))
  )
  publish_csv(output_manifest, OUTPUT$output_manifest)
  if (nrow(output_manifest) != 16L ||
      any(output_manifest$artifact %in% c("output_manifest", "completion"))) {
    stop("Output manifest self/completion exclusion failed.")
  }
  verify_output_manifest_file(
    OUTPUT$output_manifest,
    if (IS_FORMAL) AUDIT_ROOT else RUN_ROOT
  )

  STATE$stage <- "final_input_reverification"
  end_input_manifest <- build_input_manifest()
  if (!identical(input_manifest[, .(input_name, path, sha256, expected_sha256, bytes)],
                 end_input_manifest[, .(input_name, path, sha256,
                                         expected_sha256, bytes)])) {
    stop("Governed input hashes/metadata changed before completion.")
  }

  completion_lines <- c(
    if (IS_FORMAL) "RUN COMPLETED" else "NONFORMAL RUN COMPLETED",
    paste0("run_mode=", RUN_MODE), paste0("run_id=", RUN_ID),
    paste0("script_sha256=", script_hash),
    paste0("spec_sha256=", SPEC_SHA256),
    paste0("input_manifest_sha256=", sha256_file(OUTPUT$input_manifest)),
    paste0("output_manifest_sha256=", sha256_file(OUTPUT$output_manifest)),
    paste0("checkpoint_inventory_sha256=", sha256_file(OUTPUT$checkpoint_inventory)),
    paste0("task_manifest_sha256=", sha256_file(OUTPUT$task_manifest)),
    paste0("seed_registry_sha256=", sha256_file(OUTPUT$seed_registry)),
    paste0("qc_sha256=", sha256_file(OUTPUT$qc)),
    paste0("task_n=", nrow(tasks)),
    paste0("repeat_split_row_n=", nrow(repeat_rows)),
    paste0("seed_n=", nrow(registry)),
    "failed_n=0", "fallback_n=0", "warning_n=0",
    "kmeans_algorithm=Lloyd", "kmeans_nstart=25", "kmeans_itermax=100",
    paste0("reportable_formal_result=", if (IS_FORMAL) "TRUE" else "FALSE"),
    "downstream_authorization=FALSE",
    paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )

  STATE$stage <- "lock_release_and_completion"
  release_formal_lock()
  publish_text(completion_lines, OUTPUT$completion)
  invisible(list(gate = gate_table, evaluation = eval_rows, qc = qc))
}

tryCatch(
  withCallingHandlers(
    wp9_main(),
    warning = function(w) {
      stop("Top-level warning under strict fail-closed contract: ",
           conditionMessage(w), call. = FALSE)
    }
  ),
  error = function(e) {
    try(future::plan(future::sequential), silent = TRUE)
    STATE$log_enabled <- FALSE
    # A formal process that never acquired the exclusive lock must leave the
    # active formal run and its output root completely untouched.
    marker_error <- NULL
    if (!IS_FORMAL || isTRUE(STATE$lock_ever_acquired)) {
      tryCatch(
        write_failure_marker(conditionMessage(e)),
        error = function(m) marker_error <<- conditionMessage(m)
      )
    }
    if (isTRUE(STATE$lock_acquired)) {
      try(release_formal_lock(), silent = TRUE)
    }
    if (!is.null(marker_error)) {
      stop("Original failure: ", conditionMessage(e),
           " | failure-marker publication also failed: ", marker_error,
           call. = FALSE)
    }
    stop(e)
  }
)
