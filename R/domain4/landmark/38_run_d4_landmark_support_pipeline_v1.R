args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else Sys.getenv("D4_RUN_MODE", "smoke")
if (!mode %in% c("smoke", "formal")) {
  stop("Usage: Rscript 38_run_d4_landmark_support_pipeline_v1.R smoke|formal")
}

repo_root <- normalizePath(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
root <- file.path(repo_root, "R", "domain4")
script_dir <- file.path(root, "landmark")
output_root <- if (mode == "formal") {
  file.path(root, "output_formal_v1")
} else {
  file.path(root, paste0("output_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
if (dir.exists(output_root) && length(list.files(output_root, all.files = TRUE)) > 0L) {
  stop("Output directory is not empty: ", output_root, call. = FALSE)
}
Sys.setenv(D4_RUN_MODE = mode, D4_OUTPUT_ROOT = output_root)

scripts <- c(
  "38a_d4_landmark_source_audit_v1.R",
  "38b_d4_build_24h_opportunity_riskset_v1.R",
  "38c_d4_fit_k2_propensity_support_v1.R",
  "38d_d4_apply_locked_gates_v1.R",
  "38e_d4_landmark_sensitivity_v1.R",
  "38f_d4_release_and_manuscript_tables_v1.R",
  "38g_d4_independent_qc_v1.R"
)
cat("D4 pipeline mode:", mode, "\n")
cat("D4 output:", output_root, "\n")
for (s in scripts) {
  cat("\n--- running", s, "---\n")
  source(file.path(script_dir, s), encoding = "UTF-8")
}
cat("\nD4 pipeline completed:", output_root, "\n")
