## ============================================================
## 25_domain2_fix_missing_outputs_v2.R
## Harmonize Domain 2 outputs after MICE-primary re-run
##
## Purpose:
## 1) Use actual output filenames from 24_crossfit_k2_b2_fig_v4.R.
## 2) Create canonical filenames expected by downstream Figure/Table scripts.
## 3) Avoid unnecessary heavy re-runs when existing results are already present.
## 4) Create sensitivity_b2_redundancy.csv from cross-fitted B2 summary if
##    22_sensitivity_aggregate_mice_safe.R did not explicitly write it.
## ============================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

suppressPackageStartupMessages({
  library(data.table)
})

QC_DIR  <- file.path(DIR_OUTPUT, "qc")
FIG_DIR <- file.path(DIR_OUTPUT, "fig")
dir.create(QC_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n========== Domain 2 harmonization v2 ==========" , "\n")

copy_file_safe <- function(from, to, overwrite = TRUE) {
  if (length(from) == 0 || is.na(from) || !file.exists(from)) return(FALSE)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = overwrite)
  cat("Copied:\n  ", from, "\n  -> ", to, "\n", sep = "")
  isTRUE(ok)
}

read_dt_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(fread(path), error = function(e) NULL)
}

## ------------------------------------------------------------
## 1. Inventory all Domain 2-related files
## ------------------------------------------------------------
d2_all <- list.files(
  DIR_OUTPUT,
  pattern = "mice_b2|crossfit|figure4|posctrl|redund|sensitivity|incremental|dauc|auc_ladder",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

d2_inventory <- data.table(
  file = basename(d2_all),
  path = d2_all,
  size_MB = round(file.info(d2_all)$size / 1024^2, 4),
  modified_time = as.character(file.info(d2_all)$mtime)
)

d2_inventory <- d2_inventory[order(-as.POSIXct(modified_time))]
fwrite(d2_inventory, file.path(QC_DIR, "domain2_all_related_outputs_inventory_v2.csv"))
cat("Inventory written: ", file.path(QC_DIR, "domain2_all_related_outputs_inventory_v2.csv"), "\n", sep = "")

## ------------------------------------------------------------
## 2. Canonicalize cross-fit data filenames
## Actual 24 script writes:
##   crossfit_k2_b2_figdata_auc_long.csv
##   crossfit_k2_b2_figdata_decisive_dauc.csv
## Downstream scripts may expect:
##   crossfit_k2_b2_auc_ladder.csv
##   crossfit_k2_b2_decisive_dauc_forest.csv
## ------------------------------------------------------------
actual_auc_long <- file.path(QC_DIR, "crossfit_k2_b2_figdata_auc_long.csv")
actual_dauc     <- file.path(QC_DIR, "crossfit_k2_b2_figdata_decisive_dauc.csv")
canon_auc       <- file.path(QC_DIR, "crossfit_k2_b2_auc_ladder.csv")
canon_dauc      <- file.path(QC_DIR, "crossfit_k2_b2_decisive_dauc_forest.csv")

if (!file.exists(canon_auc) && file.exists(actual_auc_long)) {
  copy_file_safe(actual_auc_long, canon_auc)
}

if (!file.exists(canon_dauc) && file.exists(actual_dauc)) {
  copy_file_safe(actual_dauc, canon_dauc)
}

## If still missing, rerun 24 only once; it uses checkpoints and should be fast.
if (!file.exists(canon_auc) || !file.exists(canon_dauc)) {
  script24 <- file.path(file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "r"), "24_crossfit_k2_b2_fig_v4.R")
  if (file.exists(script24)) {
    cat("\nCanonical cross-fit files still incomplete. Rerunning 24 script using checkpoints...\n")
    source(script24)
    if (!file.exists(canon_auc) && file.exists(actual_auc_long)) copy_file_safe(actual_auc_long, canon_auc)
    if (!file.exists(canon_dauc) && file.exists(actual_dauc)) copy_file_safe(actual_dauc, canon_dauc)
  } else {
    warning("24_crossfit_k2_b2_fig_v4.R not found: ", script24)
  }
}

## ------------------------------------------------------------
## 3. Ensure Figure 4 publication PNG exists
## Prefer rerunning publication Figure 4 script. If unavailable, convert PDF
## via magick when available. As a last resort, copy existing figure4_incremental.png
## but record this in the audit note.
## ------------------------------------------------------------
fig4_png <- file.path(FIG_DIR, "figure4_incremental_publication.png")
fig4_pdf <- file.path(FIG_DIR, "figure4_incremental_publication.pdf")
fig4_png_source_note <- NA_character_

if (!file.exists(fig4_png)) {
  script_fig4 <- file.path(file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "r"), "fig4_incremental_publication_v4_tight_aligned.R")
  if (file.exists(script_fig4)) {
    cat("\nRerunning Figure 4 publication script...\n")
    source(script_fig4)
  }
}

if (!file.exists(fig4_png) && file.exists(fig4_pdf) && requireNamespace("magick", quietly = TRUE)) {
  cat("\nConverting Figure 4 publication PDF to PNG via magick...\n")
  img <- magick::image_read_pdf(fig4_pdf, density = 300)
  magick::image_write(img[1], path = fig4_png, format = "png")
  fig4_png_source_note <- "converted_from_figure4_incremental_publication.pdf"
}

if (!file.exists(fig4_png)) {
  fallback_png <- file.path(FIG_DIR, "figure4_incremental.png")
  if (file.exists(fallback_png)) {
    copy_file_safe(fallback_png, fig4_png)
    fig4_png_source_note <- "copied_from_figure4_incremental.png_fallback_check_visual_consistency"
  }
}

## ------------------------------------------------------------
## 4. Create sensitivity_b2_redundancy.csv when absent
## Source priority:
##   A) existing 22 output if already present;
##   B) crossfit_k2_b2_summary.csv, which directly contains auc_feat,
##      auc_feat_lab, and decisive_dAUC for mortality_30d and make30.
## ------------------------------------------------------------
redun_file <- file.path(QC_DIR, "sensitivity_b2_redundancy.csv")

if (!file.exists(redun_file)) {
  cf_sum_file <- file.path(QC_DIR, "crossfit_k2_b2_summary.csv")
  cf <- read_dt_safe(cf_sum_file)
  lab_file <- file.path(DIR_MODEL, "labels_primary_mice.rds")
  n_mimic <- NA_integer_
  if (file.exists(lab_file)) {
    lab <- as.data.table(readRDS(lab_file))
    n_mimic <- nrow(lab)
  }
  
  if (!is.null(cf) && all(c("outcome", "auc_feat", "auc_feat_lab", "decisive_dAUC") %in% names(cf))) {
    get1 <- function(outcome, col) {
      val <- cf[outcome == !!outcome, get(col)]
      if (length(val) == 0) return(NA_real_)
      suppressWarnings(as.numeric(val[1]))
    }
    get1c <- function(outcome, col) {
      val <- cf[outcome == !!outcome, get(col)]
      if (length(val) == 0) return(NA_real_)
      suppressWarnings(as.numeric(val[1]))
    }
    redun <- data.table(
      clustering = "kmeans_k2_mice_crossfit",
      label_source = "strict fold-wise MICE K2 labels; summary from crossfit_k2_b2_summary.csv",
      label_column = "cluster_k2",
      n = n_mimic,
      mort_auc_feat = get1("mortality_30d", "auc_feat"),
      mort_auc_feat_lab = get1("mortality_30d", "auc_feat_lab"),
      mort_dAUC = get1("mortality_30d", "decisive_dAUC"),
      mort_dAUC_lo = if ("dAUC_lo" %in% names(cf)) get1("mortality_30d", "dAUC_lo") else NA_real_,
      mort_dAUC_hi = if ("dAUC_hi" %in% names(cf)) get1("mortality_30d", "dAUC_hi") else NA_real_,
      make_auc_feat = get1("make30", "auc_feat"),
      make_auc_feat_lab = get1("make30", "auc_feat_lab"),
      make_dAUC = get1("make30", "decisive_dAUC"),
      make_dAUC_lo = if ("dAUC_lo" %in% names(cf)) get1("make30", "dAUC_lo") else NA_real_,
      make_dAUC_hi = if ("dAUC_hi" %in% names(cf)) get1("make30", "dAUC_hi") else NA_real_,
      interpretation = "Adding the K2 subphenotype label to the 33-feature model produced no material AUC increment."
    )
    fwrite(redun, redun_file)
    cat("\nCreated sensitivity_b2_redundancy.csv from crossfit_k2_b2_summary.csv:\n")
    print(redun)
  } else {
    warning("Cannot create sensitivity_b2_redundancy.csv: crossfit_k2_b2_summary.csv missing or lacks required columns.")
  }
}

## ------------------------------------------------------------
## 5. Optional compact Domain 2 interpretation table
## ------------------------------------------------------------
domain2_interpretation_file <- file.path(QC_DIR, "domain2_decisive_interpretation.csv")

mice_b2 <- read_dt_safe(file.path(QC_DIR, "mice_b2_rubin_only.csv"))
cf_sum  <- read_dt_safe(file.path(QC_DIR, "crossfit_k2_b2_summary.csv"))
redun   <- read_dt_safe(redun_file)

interpretation <- data.table(
  result_component = c(
    "Baseline + phenotype incremental value",
    "Raw 33-feature redundancy control",
    "Strict cross-fitted MICE sensitivity",
    "Positive-control sanity check"
  ),
  primary_file = c(
    "mice_b2_rubin_only.csv",
    "sensitivity_b2_redundancy.csv",
    "crossfit_k2_b2_summary.csv",
    "posctrl_incremental_power.csv"
  ),
  conclusion = c(
    "K2 labels improve prediction beyond sparse clinical baseline.",
    "Adding K2 labels to a 33-feature model provides no material additional discrimination.",
    "Redundancy persists under strict fold-wise MICE and fold-wise label assignment.",
    "Positive controls verify that the pipeline can detect incremental structure when it exists."
  )
)

## Add compact numeric fields when available
if (!is.null(redun) && nrow(redun) > 0) {
  interpretation[result_component == "Raw 33-feature redundancy control", 
                 numeric_summary := paste0(
                   "mortality dAUC=", signif(redun$mort_dAUC[1], 3),
                   "; MAKE dAUC=", signif(redun$make_dAUC[1], 3)
                 )]
}

if (!is.null(cf_sum) && nrow(cf_sum) > 0) {
  mrow <- cf_sum[outcome == "mortality_30d"]
  mkrow <- cf_sum[outcome == "make30"]
  if (nrow(mrow) > 0 && nrow(mkrow) > 0) {
    interpretation[result_component == "Strict cross-fitted MICE sensitivity",
                   numeric_summary := paste0(
                     "mortality base+phen AUC=", sprintf("%.4f", mrow$auc_base_phen[1]),
                     ", raw33 AUC=", sprintf("%.4f", mrow$auc_feat[1]),
                     "; MAKE base+phen AUC=", sprintf("%.4f", mkrow$auc_base_phen[1]),
                     ", raw33 AUC=", sprintf("%.4f", mkrow$auc_feat[1])
                   )]
  }
}

fwrite(interpretation, domain2_interpretation_file)
cat("\nWrote: ", domain2_interpretation_file, "\n", sep = "")

## ------------------------------------------------------------
## 6. Final canonical audit
## ------------------------------------------------------------
domain2_files <- data.table(
  module = c(
    "MICE B2 Rubin primary",
    "MICE B2 per-imputation",
    "MICE B2 bootstrap variance",
    "Cross-fit K2/B2 summary",
    "Cross-fit AUC ladder data",
    "Cross-fit decisive dAUC data",
    "Positive control power",
    "Figure 4 publication PNG",
    "Figure 4 publication PDF",
    "Figure 4 panel A data",
    "Figure 4 panel B data",
    "Figure 4 panel C data",
    "Sensitivity aggregate B2 redundancy",
    "Sensitivity provenance audit",
    "Domain 2 interpretation table"
  ),
  path = c(
    file.path(QC_DIR, "mice_b2_rubin_only.csv"),
    file.path(QC_DIR, "mice_b2_per_imputation.csv"),
    file.path(QC_DIR, "mice_b2_bootvar_byimp.csv"),
    file.path(QC_DIR, "crossfit_k2_b2_summary.csv"),
    file.path(QC_DIR, "crossfit_k2_b2_auc_ladder.csv"),
    file.path(QC_DIR, "crossfit_k2_b2_decisive_dauc_forest.csv"),
    file.path(QC_DIR, "posctrl_incremental_power.csv"),
    file.path(FIG_DIR, "figure4_incremental_publication.png"),
    file.path(FIG_DIR, "figure4_incremental_publication.pdf"),
    file.path(QC_DIR, "figure4_panelA_auc_ladder.csv"),
    file.path(QC_DIR, "figure4_panelB_decisive_dauc.csv"),
    file.path(QC_DIR, "figure4_panelC_posctrl.csv"),
    file.path(QC_DIR, "sensitivity_b2_redundancy.csv"),
    file.path(QC_DIR, "sensitivity_mice_provenance_audit.csv"),
    domain2_interpretation_file
  )
)

domain2_files[, exists := file.exists(path)]
domain2_files[, size_MB := ifelse(exists, round(file.info(path)$size / 1024^2, 4), NA_real_)]
domain2_files[, modified_time := ifelse(exists, as.character(file.info(path)$mtime), NA_character_)]
domain2_files[, note := NA_character_]
domain2_files[module == "Figure 4 publication PNG", note := fig4_png_source_note]

cat("\n========== Final Domain 2 audit v2 ==========" , "\n")
print(domain2_files)

fwrite(domain2_files, file.path(QC_DIR, "domain2_required_outputs_audit_after_fix_v2.csv"))

missing <- domain2_files[exists == FALSE]
if (nrow(missing) > 0) {
  cat("\nStill missing canonical outputs:\n")
  print(missing)
} else {
  cat("\nAll canonical Domain 2 outputs are present.\n")
}

cat("\nFinished Domain 2 harmonization v2.\n")
