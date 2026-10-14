ehr_audit_assert_external_output <- function(output_root, repo_root) {
  if (!nzchar(output_root) || !nzchar(repo_root) ||
      !dir.exists(output_root) || !dir.exists(repo_root)) {
    stop("Set existing EHR_AUDIT_OUTPUT_ROOT and EHR_AUDIT_REPO_ROOT directories.")
  }
  out <- tolower(normalizePath(output_root, winslash = "/", mustWork = TRUE))
  repo <- tolower(normalizePath(repo_root, winslash = "/", mustWork = TRUE))
  if (identical(out, repo) || startsWith(out, paste0(repo, "/"))) {
    stop("Restricted analysis output must be outside the Git repository.")
  }
  invisible(out)
}
