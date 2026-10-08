#!/usr/bin/env Rscript
# Publication numbering wrapper; the underlying plotting scripts are unchanged.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
root <- dirname(normalizePath(sub("^--file=", "", arg), winslash = "/"))
requested <- commandArgs(TRUE)
ids <- if (length(requested)) as.integer(requested) else seq_len(13L)
stopifnot(length(ids) > 0L, all(!is.na(ids)), all(ids %in% seq_len(13L)))
legacy <- c(NA_integer_, 8L, 6L, 7L, 12L, 1L, 2L, 3L, 4L, 5L,
            9L, 10L, 11L)
stems <- c(
  "S1_MICE_diagnostics_v5_submission_main",
  "S2_profile_harmonization_v4_submission",
  "S3_D1_lineage_ecdf_v5",
  "S4_domain2_label_lineage_v5_submission_main",
  "S5_surrogate_application_sensitivity_v4_submission",
  "S6_D4_support_v5_submission_main",
  "S7_D1_split_sensitivity_v5_submission_main",
  "S8_oracle_accuracy_dose_response_v4_submission",
  "S9_separation_reproducibility_v4_submission",
  "S10_missingness_robustness_v4_submission",
  "S11_algorithm_robustness_v4_submission",
  "S12_D3_reconstructed_generator_validation_v4_submission"
)
exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")
run <- function(script, args = character()) {
  status <- system2(exe, c("--vanilla", shQuote(script), args), wait = TRUE)
  if (status != 0L) stop("Plot failed: ", script)
}
if (1L %in% ids) {
  run(file.path(root, "..", "main", "scripts", "SF_RiskSet.R"))
}
old <- unique(legacy[ids[ids != 1L]])
if (length(old)) {
  run(file.path(root, "run_all_supp_figures.R"), as.character(old))
}
out <- file.path(root, "Figures_current")
for (i in ids) {
  for (ext in c("pdf", "png", "svg")) {
    if (i == 1L) {
      source <- file.path(root, "..", "main", if (ext == "png") "qa" else "figures",
                          paste0("SF_RiskSet.", ext))
    } else {
      source <- file.path(root, "figures", "supplement",
                          paste0("Supplementary_Figure_", stems[legacy[i]], ".", ext))
    }
    if (i == 1L && ext == "svg") next
    if (!file.exists(source)) stop("Expected figure output is missing: ", source)
    target_dir <- file.path(out, toupper(ext))
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    target <- file.path(target_dir, paste0("Supplementary_Figure_S", i, ".", ext))
    stopifnot(file.copy(source, target, overwrite = TRUE))
    stopifnot(identical(unname(tools::md5sum(source)),
                        unname(tools::md5sum(target))))
  }
}
message("Current supplementary numbering rendered for S", paste(ids, collapse = ", S"))
