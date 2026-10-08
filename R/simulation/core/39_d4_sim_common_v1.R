options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

SIM_PROJECT_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
SIM_ROOT <- file.path(SIM_PROJECT_ROOT, "R", "simulation")
SIM_SCRIPT_DIR <- file.path(SIM_ROOT, "core")
SIM_SPEC <- file.path(
  SIM_ROOT, "spec", "D4_SIM_PRERUN_SPEC_FROZEN_20260801_v2.md"
)
SIM_EXPECTED_HASHES <- file.path(
  SIM_ROOT, "spec", "D4_SIM_PUBLIC_SOURCE_HASHES_v2.csv"
)
SIM_EMPIRICAL_COMMON <- file.path(
  SIM_PROJECT_ROOT, "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
)
SIM_SUPPORT_ENGINE <- file.path(
  SIM_SCRIPT_DIR, "39_d4_support_engine_frozen_v1.R"
)

SIM_RUN_MODE <- Sys.getenv("D4_SIM_RUN_MODE", unset = "smoke")
if (!SIM_RUN_MODE %in% c("smoke", "formal")) {
  stop("D4_SIM_RUN_MODE must be smoke or formal.", call. = FALSE)
}

SIM_OUTPUT_ROOT <- Sys.getenv(
  "D4_SIM_OUTPUT_ROOT",
  unset = if (SIM_RUN_MODE == "formal") {
    file.path(SIM_PROJECT_ROOT, "outputs", "simulation", "D4_G5_G6_formal")
  } else {
    file.path(SIM_PROJECT_ROOT, "outputs", "simulation", paste0(
      "output_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S")
    ))
  }
)
SIM_TABLE_DIR <- file.path(SIM_OUTPUT_ROOT, "tables")
SIM_FIGURE_DIR <- file.path(SIM_OUTPUT_ROOT, "figures")
SIM_LOG_DIR <- file.path(SIM_OUTPUT_ROOT, "logs")
SIM_PROVENANCE_DIR <- file.path(SIM_OUTPUT_ROOT, "provenance")
SIM_CHECKPOINT_DIR <- file.path(SIM_OUTPUT_ROOT, "checkpoints")

SIM_BASE_SEED <- 20260801L
SIM_N <- if (SIM_RUN_MODE == "formal") 20049L else 2000L
SIM_REPS <- if (SIM_RUN_MODE == "formal") {
  c(G4L = 300L, G6L = 300L, D4LQ = 100L)
} else {
  c(G4L = 5L, G6L = 5L, D4LQ = 3L)
}
SIM_SCENARIOS <- c(
  G4L = "G4L_longitudinal_support_failure",
  G6L = "G6L_longitudinal_adequate_support",
  D4LQ = "D4LQ_ledger_stress_control"
)
SIM_EXPECTED_STATE <- c(
  G4L_longitudinal_support_failure = "SUPPORT_INADEQUATE",
  G6L_longitudinal_adequate_support =
    "SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT",
  D4LQ_ledger_stress_control = NA_character_
)

sim_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      ..., "\n", sep = "")
}

sim_require_files <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss)) {
    stop("Missing required files:\n", paste(miss, collapse = "\n"),
         call. = FALSE)
  }
}

sim_init_dirs <- function() {
  if (SIM_RUN_MODE == "formal" && file.exists(SIM_OUTPUT_ROOT)) {
    stop("Refusing to overwrite formal output: ", SIM_OUTPUT_ROOT,
         call. = FALSE)
  }
  for (d in c(
    SIM_OUTPUT_ROOT, SIM_TABLE_DIR, SIM_FIGURE_DIR, SIM_LOG_DIR,
    SIM_PROVENANCE_DIR, SIM_CHECKPOINT_DIR
  )) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sim_atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(as.data.table(x), tmp, na = "")
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path)
  invisible(path)
}

sim_atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path)
  invisible(path)
}

sim_sha256 <- function(path) {
  toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

sim_seed <- function(scenario, repeat_id) {
  offset <- match(scenario, unname(SIM_SCENARIOS)) * 100000L
  SIM_BASE_SEED + offset + as.integer(repeat_id)
}

sim_clamp <- function(x, low, high) pmin(pmax(x, low), high)

sim_wilson <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, ctr - half), high = min(1, ctr + half))
}

sim_safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

sim_safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) sd(x) else NA_real_
}

sim_q <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x)) as.numeric(quantile(x, p, names = FALSE)) else NA_real_
}

sim_source_empirical_engine <- function() {
  sim_require_files(c(SIM_SPEC, SIM_EMPIRICAL_COMMON, SIM_SUPPORT_ENGINE))
  source(SIM_SUPPORT_ENGINE, local = .GlobalEnv)
  required <- c("d4_fit_support", "d4_apply_state", "D4_GATES",
                "D4_SUPPORT_THRESHOLDS", "D4_PS_COVARIATES")
  absent <- required[!vapply(required, exists, logical(1), inherits = TRUE)]
  if (length(absent)) stop("Empirical D4 engine missing: ", paste(absent, collapse = ", "))
  invisible(TRUE)
}
