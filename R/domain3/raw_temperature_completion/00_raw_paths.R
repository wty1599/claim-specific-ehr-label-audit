# Public dependency and restricted-artifact mapping; no clinical rows are read.
source(file.path(repo, "R/common/00_repo_paths.R"))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(EHR_AUDIT_RESTRICTED_DATA_ROOT, repo)
ehr_audit_assert_external_output(EHR_AUDIT_OUTPUT_ROOT, repo)
RAW_CODE_ROOT <- file.path(repo, "R/domain3/raw_temperature_completion")
DATA_ROOT <- normalizePath(EHR_AUDIT_RESTRICTED_DATA_ROOT, winslash = "/", mustWork = TRUE)
ROOT <- file.path(EHR_AUDIT_OUTPUT_ROOT, "domain3_raw_temperature_completion")
OLD_ROOT <- file.path(EHR_AUDIT_OUTPUT_ROOT, "domain3_deployable_generator_v1_20260728")
AUDIT_ROOT <- file.path(EHR_AUDIT_OUTPUT_ROOT, "domain3_source_correction")
CORRECTION_FIGURE_SOURCE <- Sys.getenv("EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT",
  unset = file.path(AUDIT_ROOT, "figure4/source_data"))
ehr_audit_assert_external_output(CORRECTION_FIGURE_SOURCE, repo)
for (folder in c("private_local", "aggregate_outputs", "figure_source", "provenance", "logs")) {
  dir.create(file.path(ROOT, folder), recursive = TRUE, showWarnings = FALSE)
}
PATHS <- c(
  specification = file.path(RAW_CODE_ROOT, "spec/RAW_MODEL_COMPLETION_SPEC.md"),
  original_script = file.path(repo, "R/domain3/generator/04_analyze_d3_generator_transport_v1.R"),
  original_config = file.path(repo, "R/domain3/generator/00_d3_generator_config_v1.R"),
  original_utils = file.path(repo, "R/domain3/generator/00_d3_generator_utils_v1.R"),
  original_session = file.path(OLD_ROOT, "logs/d3_generator_v1_transport_sessionInfo.txt"),
  development = file.path(DATA_ROOT, "final_full.csv"),
  baseline = file.path(DATA_ROOT, "baseline_covars.csv"),
  external = file.path(DATA_ROOT, "eicu_external.csv"),
  generator = file.path(OLD_ROOT, "model/saaki_k2_deployable_surrogate_generator_v1.rds"),
  external_labels = file.path(OLD_ROOT, "labels/eicu_d3_generator_v1_labels.csv"),
  label_freeze = file.path(OLD_ROOT, "qc/eicu_d3_generator_v1_label_fingerprint.csv"),
  archived_predictions = file.path(OLD_ROOT, "private_local/d3_generator_v1_external_predictions_local.csv"),
  archived_performance = file.path(OLD_ROOT, "tables/d3_generator_v1_model_performance.csv"),
  corrected_temperature = file.path(AUDIT_ROOT, "private_local/eicu_temperature_corrected_v1.csv"),
  corrected_base_k2_predictions = file.path(AUDIT_ROOT, "private_local/corrected_external_predictions_local.csv"),
  figure_transport = file.path(CORRECTION_FIGURE_SOURCE, "Figure4_panelC_L1_transport.csv"),
  figure_calibration = file.path(CORRECTION_FIGURE_SOURCE, "Figure4_panelC_L1_calibration.csv")
)
