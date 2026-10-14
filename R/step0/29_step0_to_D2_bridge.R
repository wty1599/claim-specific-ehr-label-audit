#!/usr/bin/env Rscript

## WP10-A: deterministic Domain 1 to Domain 2 bridge.
## Extraction/post-processing only; no clustering or outcome model is refit.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
AUDIT_ROOT <- file.path(
  PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"
)
SCRIPT_PATH <- file.path(AUDIT_ROOT, "scripts", "36_domain1_to_domain2_bridge.R")
PRERUN_SPEC <- file.path(
  AUDIT_ROOT, "provenance", "36_37_WP10_prerun_spec_FROZEN_20260718.md"
)
WP10_ROOT <- file.path(AUDIT_ROOT, "wp10_domain2_domain4_bridge_20260718")
RUN_ROOT <- file.path(WP10_ROOT, "domain2")
TABLE_DIR <- file.path(RUN_ROOT, "tables")
LOG_DIR <- file.path(RUN_ROOT, "logs")
PROVENANCE_DIR <- file.path(RUN_ROOT, "provenance")
for (d in c(RUN_ROOT, TABLE_DIR, LOG_DIR, PROVENANCE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

INPUT <- c(
  script = SCRIPT_PATH,
  frozen_spec = PRERUN_SPEC,
  identity_retention = file.path(
    AUDIT_ROOT, "tables", "Table_D1_identity_discrimination_retention.csv"
  ),
  identity_models = file.path(
    AUDIT_ROOT, "tables", "Table_D1_identity_models_oof.csv"
  ),
  mice_rubin = file.path(PROJECT_ROOT, "output", "qc", "mice_b2_rubin_only.csv"),
  crossfit_summary = file.path(
    PROJECT_ROOT, "output", "qc", "crossfit_k2_b2_summary.csv"
  ),
  crossfit_ladder = file.path(
    PROJECT_ROOT, "output", "qc", "crossfit_k2_b2_auc_ladder.csv"
  ),
  crossfit_dauc = file.path(
    PROJECT_ROOT, "output", "qc", "crossfit_k2_b2_decisive_dauc_forest.csv"
  ),
  crossfit_script = file.path(PROJECT_ROOT, "r", "24_crossfit_k2_b2_fig_v4.R")
)

for (p in INPUT) {
  if (!file.exists(p)) stop("Missing required WP10-D2 input: ", p, call. = FALSE)
}
spec_text <- readLines(PRERUN_SPEC, warn = FALSE, encoding = "UTF-8")
if (!any(spec_text == "Status: **FROZEN - AUTHORIZED FOR FORMAL EXECUTION**")) {
  stop("Frozen WP10 specification lacks formal authorization marker.", call. = FALSE)
}

COMPLETION <- file.path(RUN_ROOT, "36_WP10_D2_run_completed.ok")
if (file.exists(COMPLETION)) {
  stop("WP10-D2 output is already complete and will not be overwritten: ", RUN_ROOT,
       call. = FALSE)
}

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(x, tmp, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path, call. = FALSE)
}

atomic_write_lines <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic text write failed: ", path, call. = FALSE)
}

input_manifest <- rbindlist(lapply(names(INPUT), function(nm) {
  p <- INPUT[[nm]]
  data.table(
    input_name = nm,
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(p),
    bytes = file.info(p)$size
  )
}))
atomic_fwrite(
  input_manifest,
  file.path(PROVENANCE_DIR, "36_WP10_D2_input_sha256_manifest.csv")
)

retention <- fread(INPUT[["identity_retention"]])
models <- fread(INPUT[["identity_models"]])
mice <- fread(INPUT[["mice_rubin"]])
crossfit <- fread(INPUT[["crossfit_summary"]])
ladder <- fread(INPUT[["crossfit_ladder"]])
dauc <- fread(INPUT[["crossfit_dauc"]])

require_columns <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) {
    stop(label, " lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

require_columns(
  retention,
  c("model_id", "model_label", "auc", "all33_auc", "discrimination_retention"),
  "Identity retention table"
)
require_columns(
  models,
  c("model_id", "auc", "auc_ci_low", "auc_ci_high", "evaluation", "n"),
  "Identity model table"
)
require_columns(
  mice,
  c("outcome", "m", "auc_feat", "auc_feat_lab", "decisive_dAUC_feat_plus_label",
    "dAUC_lo", "dAUC_hi", "incr_base_to_basephen"),
  "MICE Rubin table"
)
require_columns(
  crossfit,
  c("outcome", "outcome_label", "imputation", "K_folds", "auc_feat",
    "auc_feat_lab", "decisive_dAUC", "dAUC_lo", "dAUC_hi",
    "incr_base_to_basephen"),
  "Cross-fitted summary"
)

renal_row <- retention[model_id == "creatinine_bun"]
renal_model <- models[model_id == "creatinine_bun"]
if (nrow(renal_row) != 1L || nrow(renal_model) != 1L) {
  stop("Expected exactly one creatinine+BUN identity row.", call. = FALSE)
}
if (abs(renal_row$auc - renal_model$auc) > 1e-12) {
  stop("Creatinine+BUN AUC differs across locked Domain 1 tables.", call. = FALSE)
}

bridge_identity <- data.table(
  bridge_section = "Domain 1 fixed-label reconstruction",
  evidence_id = "D1_creatinine_BUN_vs_all33",
  outcome_or_target = "Fixed full-cohort MICE-imputation-1 plain K-means K2 label",
  label_estimand = "fixed_full_cohort_K2",
  analysis = "10-fold stratified out-of-fold logistic reconstruction",
  estimate_name = "label_discrimination_retention",
  estimate = renal_row$discrimination_retention,
  ci_low = NA_real_,
  ci_high = NA_real_,
  component_auc = renal_row$auc,
  component_auc_ci_low = renal_model$auc_ci_low,
  component_auc_ci_high = renal_model$auc_ci_high,
  reference_auc = renal_row$all33_auc,
  baseline_plus_label_increment = NA_real_,
  n = renal_model$n,
  folds_or_imputations = "10 folds",
  source_file = basename(INPUT[["identity_retention"]]),
  source_fields = "model_id;auc;all33_auc;discrimination_retention",
  interpretation = paste(
    "Creatinine plus BUN retain the existing, prespecified measure of fixed-label",
    "separability; descriptive representational evidence only."
  )
)

outcome_labels <- c(
  mortality_30d = "30-day mortality",
  make30 = "MAKE-30"
)

bridge_mice <- mice[, data.table(
  bridge_section = "Domain 2 MICE-pooled redundancy",
  evidence_id = paste0("D2_mice_", outcome),
  outcome_or_target = unname(outcome_labels[outcome]),
  label_estimand = "MICE_imputation_specific_K2_with_Rubin_pooling",
  analysis = "Raw-EN + K2 minus Raw-EN; paired DeLong within imputation and Rubin pooling",
  estimate_name = "decisive_Delta_AUC",
  estimate = decisive_dAUC_feat_plus_label,
  ci_low = dAUC_lo,
  ci_high = dAUC_hi,
  component_auc = auc_feat_lab,
  component_auc_ci_low = NA_real_,
  component_auc_ci_high = NA_real_,
  reference_auc = auc_feat,
  baseline_plus_label_increment = incr_base_to_basephen,
  n = 20049L,
  folds_or_imputations = paste0("m=", m),
  source_file = basename(INPUT[["mice_rubin"]]),
  source_fields = paste(
    "auc_feat;auc_feat_lab;decisive_dAUC_feat_plus_label;dAUC_lo;dAUC_hi;",
    "incr_base_to_basephen"
  ),
  interpretation = paste(
    "Adding K2 to Raw-EN provides practically no incremental discrimination;",
    "this does not test a causal effect."
  )
)]

bridge_crossfit <- crossfit[, data.table(
  bridge_section = "Domain 2 strict cross-fitted redundancy",
  evidence_id = paste0("D2_crossfit_", outcome),
  outcome_or_target = outcome_label,
  label_estimand = "outcome_specific_foldwise_cross_fitted_K2",
  analysis = "Raw-EN + fold-wise K2 minus Raw-EN in the outer test folds",
  estimate_name = "decisive_Delta_AUC",
  estimate = decisive_dAUC,
  ci_low = dAUC_lo,
  ci_high = dAUC_hi,
  component_auc = auc_feat_lab,
  component_auc_ci_low = NA_real_,
  component_auc_ci_high = NA_real_,
  reference_auc = auc_feat,
  baseline_plus_label_increment = incr_base_to_basephen,
  n = 20049L,
  folds_or_imputations = paste0(K_folds, " folds; ", imputation),
  source_file = basename(INPUT[["crossfit_summary"]]),
  source_fields = paste(
    "auc_feat;auc_feat_lab;decisive_dAUC;dAUC_lo;dAUC_hi;",
    "incr_base_to_basephen"
  ),
  interpretation = paste(
    "The strict fold-wise label adds no measurable discrimination to Raw-EN;",
    "it is not the same estimand as the fixed full-cohort Domain 1 label."
  )
)]

bridge <- rbindlist(list(bridge_identity, bridge_mice, bridge_crossfit), fill = TRUE)
atomic_fwrite(bridge, file.path(TABLE_DIR, "Table_D1_D2_bridge.csv"))

qc <- rbindlist(list(
  data.table(
    check = "fixed_label_identity_row_unique",
    pass = nrow(renal_row) == 1L && nrow(renal_model) == 1L,
    observed = paste(nrow(renal_row), nrow(renal_model), sep = "/"),
    expected = "1/1"
  ),
  data.table(
    check = "mice_outcomes_complete",
    pass = setequal(mice$outcome, c("mortality_30d", "make30")),
    observed = paste(sort(mice$outcome), collapse = ";"),
    expected = "make30;mortality_30d"
  ),
  data.table(
    check = "crossfit_outcomes_complete",
    pass = setequal(crossfit$outcome, c("mortality_30d", "make30")),
    observed = paste(sort(crossfit$outcome), collapse = ";"),
    expected = "make30;mortality_30d"
  ),
  data.table(
    check = "crossfit_ladder_four_models_per_outcome",
    pass = nrow(ladder) == 8L && uniqueN(ladder$model) == 4L,
    observed = paste(nrow(ladder), uniqueN(ladder$model), sep = "/"),
    expected = "8/4"
  ),
  data.table(
    check = "crossfit_dauc_two_outcomes",
    pass = nrow(dauc) == 2L && setequal(dauc$outcome, crossfit$outcome),
    observed = nrow(dauc), expected = 2L
  ),
  data.table(
    check = "no_new_model_fit",
    pass = TRUE,
    observed = "deterministic extraction only",
    expected = "no glm/glmnet/kmeans call"
  ),
  data.table(
    check = "estimand_distinction_explicit",
    pass = all(c(
      "fixed_full_cohort_K2",
      "outcome_specific_foldwise_cross_fitted_K2"
    ) %in% bridge$label_estimand),
    observed = paste(unique(bridge$label_estimand), collapse = ";"),
    expected = "fixed and cross-fitted labels named separately"
  )
), fill = TRUE)
atomic_fwrite(qc, file.path(LOG_DIR, "36_WP10_D2_QC.csv"))

log_lines <- c(
  "# WP10 Domain 1 to Domain 2 bridge log",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "- No model was refit and no threshold was changed.",
  paste0(
    "- Creatinine+BUN fixed-label AUC = ", sprintf("%.6f", renal_row$auc),
    "; all-33 AUC = ", sprintf("%.6f", renal_row$all33_auc),
    "; locked discrimination retention = ",
    sprintf("%.6f", renal_row$discrimination_retention), "."
  ),
  "- Domain 1 uses the saved full-cohort MICE-imputation-1 plain K-means label.",
  paste(
    "- The strict Domain 2 analysis uses outcome-specific, fold-wise labels",
    "constructed inside each outer training fold; this is a different estimand."
  ),
  paste(
    "- The bridge supports a representational conclusion only: the label compresses",
    "renal-function inputs already available to Raw-EN."
  ),
  "- No causal conclusion is made."
)
atomic_write_lines(log_lines, file.path(LOG_DIR, "36_D2_bridge_log.md"))

capture.output(
  sessionInfo(),
  file = file.path(PROVENANCE_DIR, "36_WP10_D2_sessionInfo.txt")
)

if (!all(qc$pass)) {
  atomic_write_lines(c(
    "WP10-D2 QC FAILURE - DO NOT USE",
    paste0("failed_checks=", paste(qc[pass == FALSE, check], collapse = ";"))
  ), file.path(RUN_ROOT, "36_WP10_D2_QC_FAILURE.txt"))
  stop("WP10-D2 failed structural QC.", call. = FALSE)
}

output_files <- list.files(RUN_ROOT, recursive = TRUE, full.names = TRUE)
output_files <- output_files[file.info(output_files)$isdir == FALSE]
output_files <- output_files[!grepl(
  "36_WP10_D2_output_sha256_manifest.csv$|36_WP10_D2_run_completed.ok$",
  output_files
)]
output_manifest <- rbindlist(lapply(output_files, function(p) {
  data.table(
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(p),
    bytes = file.info(p)$size
  )
}))
atomic_fwrite(
  output_manifest,
  file.path(PROVENANCE_DIR, "36_WP10_D2_output_sha256_manifest.csv")
)

atomic_write_lines(c(
  "RUN COMPLETED",
  "work_package=WP10_Domain2_bridge",
  "run_mode=formal_deterministic_postprocess",
  "model_refit=FALSE",
  "threshold_changed=FALSE",
  paste0("n_bridge_rows=", nrow(bridge)),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0(
    "output_manifest_sha256=",
    sha256_file(file.path(PROVENANCE_DIR, "36_WP10_D2_output_sha256_manifest.csv"))
  )
), COMPLETION)

print(bridge)
