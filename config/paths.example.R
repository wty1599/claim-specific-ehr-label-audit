# Copy to config/local_paths.R and adapt locally.
# Correction scripts require DATA_ROOT and OUTPUT_ROOT to exist outside the
# public code directory. EHR_AUDIT_WORK_ROOT is for historical intermediates.
# Main figures: EHR_AUDIT_MAIN_FIGURE_SOURCE_ROOT (existing external input).
# Current Figure 4: EHR_AUDIT_FIGURE4_SOURCE_ROOT (existing external input).
# Earlier correction stage: EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT.
# Raw recovery/evaluation executable: EHR_AUDIT_MODEL_RSCRIPT (R 4.6.0).
# Graphics runtime is separate (R 4.5.3). Explicit Figure 4 source arguments
# override environment/default selection in both total and single entry points.

PROJECT_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = "."),
  winslash = "/",
  mustWork = FALSE
)

DATA_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = file.path(PROJECT_ROOT, "data_restricted")),
  winslash = "/",
  mustWork = FALSE
)

OUTPUT_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_OUTPUT_ROOT", unset = file.path(PROJECT_ROOT, "outputs")),
  winslash = "/",
  mustWork = FALSE
)
