## Backward-compatible path aliases for historical analysis scripts.
##
## Preferred configuration:
##   EHR_AUDIT_REPO_ROOT            public repository root
##   EHR_AUDIT_RESTRICTED_DATA_ROOT credentialed/patient-level inputs
##   EHR_AUDIT_OUTPUT_ROOT          writable analysis output directory

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  this_file <- if (length(arg)) sub("^--file=", "", arg[[1]]) else NA_character_
}

configured_root <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = "")
if (nzchar(configured_root)) {
  repo_root <- normalizePath(configured_root, winslash = "/", mustWork = TRUE)
} else if (!is.na(this_file)) {
  repo_root <- normalizePath(
    file.path(dirname(this_file), "..", ".."),
    winslash = "/", mustWork = TRUE
  )
} else {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "common", "00_repo_paths.R"), encoding = "UTF-8")

PROJ_ROOT <- EHR_AUDIT_REPO_ROOT
DIR_DATA <- EHR_AUDIT_RESTRICTED_DATA_ROOT
DIR_R <- file.path(PROJ_ROOT, "R")
DIR_OUTPUT <- EHR_AUDIT_OUTPUT_ROOT
DIR_QC <- file.path(DIR_OUTPUT, "qc")
DIR_FIG <- file.path(DIR_OUTPUT, "figures")
DIR_MODEL <- file.path(DIR_OUTPUT, "model")
PATH_FINAL_FULL <- file.path(DIR_DATA, "final_full.csv")

for (d in c(DIR_OUTPUT, DIR_QC, DIR_FIG, DIR_MODEL)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("EHR audit repository root: ", PROJ_ROOT)
