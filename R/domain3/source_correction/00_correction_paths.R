# Shared publication-path adapter; no clinical input is read here.
source(file.path(repo, "R/common/00_repo_paths.R"))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(EHR_AUDIT_OUTPUT_ROOT, repo)
ehr_audit_assert_external_output(EHR_AUDIT_RESTRICTED_DATA_ROOT, repo)
base <- EHR_AUDIT_RESTRICTED_DATA_ROOT
project <- workspace <- EHR_AUDIT_OUTPUT_ROOT
root <- run_root <- file.path(project, "domain3_source_correction")
private <- private_dir <- private_out <- file.path(root, "private_local")
aggregate <- aggregate_dir <- aggregate_out <- file.path(root, "aggregate_outputs")
qa <- file.path(root, "qa")
source_root <- file.path(project, "domain3_L1")
recovery_root <- file.path(project, "domain3_parameter_recovery")
old_root <- file.path(project, "domain3_deployable_generator_v1_20260728")
source(file.path(repo, "R/common/02_figure_paths.R"))
correction_figure_source <- Sys.getenv("EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT", unset = "")
if (!nzchar(correction_figure_source)) {
  correction_figure_source <- file.path(root, "figure4/source_data")
  dir.create(correction_figure_source, recursive = TRUE, showWarnings = FALSE)
}
figure_source <- ehr_audit_figure_paths(repo, figure4_source = correction_figure_source,
                                      require_sources = character())$figure4_source
for (directory in c(private, aggregate, qa, figure_source)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
spec_root <- file.path(repo, "R/domain3/source_correction/spec")
