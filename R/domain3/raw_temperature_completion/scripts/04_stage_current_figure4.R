# Byte-preserving source staging only; no metrics or plots are recomputed.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(sub("^--file=", "", arg)), "../../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/raw_temperature_completion/00_raw_paths.R"))
source(file.path(repo, "R/common/02_figure_paths.R"))
destination <- ehr_audit_figure_paths(repo, require_sources = character())$figure4_source
files <- c("L0_L1_L2_fidelity.csv", "Figure4_panelB_public_funnel_points.csv",
  "Figure4_panelB_funnel_summary.csv", "Figure4_pooled_L1_prevalence.csv",
  "heldout_calibration_curves.csv", "heldout_risk_distribution_screened.csv",
  "heldout_performance_point.csv", "hospital_holdout_split_summary.csv",
  "Figure4_panelC_L1_transport.csv", "Figure4_panelC_L1_calibration.csv", "Figure4_reference_ranges.csv")
source_files <- file.path(CORRECTION_FIGURE_SOURCE, files)
raw_rows <- grepl("^Figure4_panelC_", files)
source_files[raw_rows] <- file.path(ROOT, "figure_source", files[raw_rows])
stopifnot(all(file.exists(source_files)))
qc <- jsonlite::fromJSON(file.path(ROOT, "provenance/completion_qc.json"))
stopifnot(qc$status == "PASS", qc$bootstrap_failures == 0L)
if (!dir.exists(destination)) dir.create(destination, recursive = TRUE, showWarnings = FALSE)
ehr_audit_assert_external_output(destination, repo)
targets <- file.path(destination, files)
for (i in seq_along(files)) {
  if (file.exists(targets[i])) {
    stopifnot(unname(tools::md5sum(targets[i])) == unname(tools::md5sum(source_files[i])))
  } else stopifnot(file.copy(source_files[i], targets[i], overwrite = FALSE))
}
stopifnot(identical(unname(tools::md5sum(source_files)), unname(tools::md5sum(targets))))
write.csv(data.frame(file = files, origin = ifelse(raw_rows, "raw_completion", "source_correction"),
  md5 = unname(tools::md5sum(targets))), file.path(ROOT, "provenance/current_figure_source_manifest.csv"),
  row.names = FALSE)
message("Current Figure 4 sources staged without changing their bytes.")
