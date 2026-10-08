## Public repository path helpers.
## Set EHR_AUDIT_REPO_ROOT to run scripts from outside the repository root.

ehr_audit_find_repo_root <- function(start = getwd()) {
  configured <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = "")
  if (nzchar(configured)) {
    return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  }

  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "README.md")) &&
        dir.exists(file.path(current, "R")) &&
        dir.exists(file.path(current, "sql"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }

  stop(
    "Repository root was not found. Run from the repository or set ",
    "EHR_AUDIT_REPO_ROOT.",
    call. = FALSE
  )
}

EHR_AUDIT_REPO_ROOT <- ehr_audit_find_repo_root()
EHR_AUDIT_RESTRICTED_DATA_ROOT <- Sys.getenv(
  "EHR_AUDIT_RESTRICTED_DATA_ROOT",
  unset = file.path(EHR_AUDIT_REPO_ROOT, "data_restricted")
)
EHR_AUDIT_OUTPUT_ROOT <- Sys.getenv(
  "EHR_AUDIT_OUTPUT_ROOT",
  unset = file.path(EHR_AUDIT_REPO_ROOT, "outputs")
)

ehr_audit_repo_path <- function(...) file.path(EHR_AUDIT_REPO_ROOT, ...)
ehr_audit_data_path <- function(...) file.path(EHR_AUDIT_RESTRICTED_DATA_ROOT, ...)
ehr_audit_output_path <- function(...) file.path(EHR_AUDIT_OUTPUT_ROOT, ...)
