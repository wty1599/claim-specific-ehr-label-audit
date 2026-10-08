## ============================================================
## 25_domain2_finalize_canonical_outputs_v3.R
## Finalize canonical Domain 2 output filenames after running
## 25_domain2_fix_missing_outputs.R / 24_crossfit_k2_b2_fig_v4.R.
##
## Purpose:
##   1) Create canonical crossfit_k2_b2_auc_ladder.csv from the
##      actual 24-script output crossfit_k2_b2_figdata_auc_long.csv.
##   2) Create figure4_incremental_publication.png from the publication PDF
##      when possible; if conversion packages are unavailable, keep PDF as
##      the authoritative submission figure and record the PNG missingness.
##   3) Produce a final Domain 2 audit and file map.
## ============================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

suppressPackageStartupMessages({
  library(data.table)
})

QC_DIR  <- file.path(DIR_OUTPUT, "qc")
FIG_DIR <- file.path(DIR_OUTPUT, "fig")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n========== Domain 2 canonical finalization v3 ==========\n")

copy_if_exists <- function(from, to, label) {
  if (file.exists(from)) {
    ok <- file.copy(from, to, overwrite = TRUE)
    cat(sprintf("%s: %s -> %s | copied=%s\n", label, from, to, ok))
    return(ok)
  }
  cat(sprintf("%s source missing: %s\n", label, from))
  FALSE
}

## ------------------------------------------------------------
## 1. Canonical cross-fit AUC ladder CSV
## ------------------------------------------------------------
canon_auc <- file.path(QC_DIR, "crossfit_k2_b2_auc_ladder.csv")
actual_auc <- file.path(QC_DIR, "crossfit_k2_b2_figdata_auc_long.csv")

if (!file.exists(canon_auc)) {
  copy_if_exists(actual_auc, canon_auc, "Create canonical cross-fit AUC ladder")
}

## Validate structure lightly
if (file.exists(canon_auc)) {
  auc_dt <- fread(canon_auc)
  cat("\nCanonical AUC ladder file exists. Dimensions:\n")
  print(dim(auc_dt))
  print(names(auc_dt))
}

## ------------------------------------------------------------
## 2. Canonical cross-fit decisive dAUC CSV
## ------------------------------------------------------------
canon_dauc <- file.path(QC_DIR, "crossfit_k2_b2_decisive_dauc_forest.csv")
actual_dauc <- file.path(QC_DIR, "crossfit_k2_b2_figdata_decisive_dauc.csv")

if (!file.exists(canon_dauc)) {
  copy_if_exists(actual_dauc, canon_dauc, "Create canonical cross-fit decisive dAUC")
}

## ------------------------------------------------------------
## 3. Publication Figure 4 PNG
## ------------------------------------------------------------
fig4_pdf <- file.path(FIG_DIR, "figure4_incremental_publication.pdf")
fig4_png <- file.path(FIG_DIR, "figure4_incremental_publication.png")
fig4_png_note <- "not_attempted"

## If PNG is absent, first rerun the Figure 4 publication script.
if (!file.exists(fig4_png)) {
  fig4_script <- file.path(file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "r"), "fig4_incremental_publication_v4_tight_aligned.R")
  if (file.exists(fig4_script)) {
    cat("\nRerunning Figure 4 publication script to request PNG/PDF output...\n")
    tryCatch(
      source(fig4_script),
      error = function(e) warning("Figure 4 script failed: ", conditionMessage(e))
    )
  } else {
    warning("Figure 4 script not found: ", fig4_script)
  }
}

## If still absent, convert PDF -> PNG if possible.
if (!file.exists(fig4_png) && file.exists(fig4_pdf)) {
  if (requireNamespace("magick", quietly = TRUE)) {
    cat("\nConverting Figure 4 PDF to PNG using magick...\n")
    tryCatch({
      img <- magick::image_read_pdf(fig4_pdf, density = 300)
      magick::image_write(img[1], path = fig4_png, format = "png")
      fig4_png_note <- "created_from_pdf_by_magick"
    }, error = function(e) {
      warning("magick PDF conversion failed: ", conditionMessage(e))
      fig4_png_note <<- "magick_conversion_failed"
    })
  } else if (requireNamespace("pdftools", quietly = TRUE) && requireNamespace("png", quietly = TRUE)) {
    cat("\nConverting Figure 4 PDF to PNG using pdftools + png...\n")
    tryCatch({
      bitmap <- pdftools::pdf_render_page(fig4_pdf, page = 1, dpi = 300)
      arr <- aperm(bitmap, c(2, 3, 1)) / 255
      png::writePNG(arr, target = fig4_png)
      fig4_png_note <- "created_from_pdf_by_pdftools_png"
    }, error = function(e) {
      warning("pdftools/png PDF conversion failed: ", conditionMessage(e))
      fig4_png_note <<- "pdftools_conversion_failed"
    })
  } else {
    fig4_png_note <- "png_missing_pdf_available_no_converter_package"
    warning(
      "figure4_incremental_publication.pdf exists, but PNG is missing and neither magick nor pdftools+png is installed. ",
      "PDF remains the authoritative publication figure. Install magick or pdftools+png only if PNG is required."
    )
  }
}

## If PDF is absent but PNG exists, note it. If both absent, rerun fig4 script manually.
if (file.exists(fig4_png)) {
  cat("\nFigure 4 PNG exists: ", fig4_png, "\n", sep = "")
} else {
  cat("\nFigure 4 PNG still missing. PDF exists = ", file.exists(fig4_pdf), "\n", sep = "")
}

## ------------------------------------------------------------
## 4. Final file map
## ------------------------------------------------------------
file_map <- data.table(
  canonical_output = c(
    "crossfit_k2_b2_auc_ladder.csv",
    "crossfit_k2_b2_decisive_dauc_forest.csv",
    "figure4_incremental_publication.png",
    "figure4_incremental_publication.pdf"
  ),
  canonical_path = c(
    canon_auc,
    canon_dauc,
    fig4_png,
    fig4_pdf
  ),
  source_used = c(
    actual_auc,
    actual_dauc,
    ifelse(file.exists(fig4_png), fig4_pdf, NA_character_),
    fig4_pdf
  ),
  note = c(
    "Copied from 24-script output crossfit_k2_b2_figdata_auc_long.csv if needed.",
    "Copied from 24-script output crossfit_k2_b2_figdata_decisive_dauc.csv if needed.",
    fig4_png_note,
    "Publication PDF is acceptable as authoritative figure output."
  )
)
file_map[, exists := file.exists(canonical_path)]
file_map[, size_MB := ifelse(exists, round(file.info(canonical_path)$size / 1024^2, 4), NA_real_)]
file_map[, modified_time := ifelse(exists, as.character(file.info(canonical_path)$mtime), NA_character_)]

fwrite(file_map, file.path(QC_DIR, "domain2_canonical_file_map_v3.csv"))
cat("\n========== Canonical file map ==========" , "\n")
print(file_map)

## ------------------------------------------------------------
## 5. Final Domain 2 audit
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
    "Sensitivity provenance audit"
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
    file.path(QC_DIR, "sensitivity_mice_provenance_audit.csv")
  )
)

domain2_files[, exists := file.exists(path)]
domain2_files[, size_MB := ifelse(exists, round(file.info(path)$size / 1024^2, 4), NA_real_)]
domain2_files[, modified_time := ifelse(exists, as.character(file.info(path)$mtime), NA_character_)]

## Treat PNG as optional when PDF exists. Keep strict_exists for literal file presence.
domain2_files[, required_for_submission := TRUE]
domain2_files[module == "Figure 4 publication PNG", required_for_submission := FALSE]
domain2_files[, submission_ready := fifelse(module == "Figure 4 publication PNG", TRUE, exists)]

fwrite(domain2_files, file.path(QC_DIR, "domain2_required_outputs_audit_after_fix_v3.csv"))
cat("\n========== Final Domain 2 audit v3 ==========" , "\n")
print(domain2_files)

missing_required <- domain2_files[required_for_submission == TRUE & exists == FALSE]
if (nrow(missing_required) > 0) {
  cat("\nStill missing required outputs:\n")
  print(missing_required)
} else {
  cat("\nAll required Domain 2 outputs are present. Figure 4 PNG is optional if PDF exists.\n")
}

cat("\nDone. Key audit files:\n")
cat("\n - ", file.path(QC_DIR, "domain2_canonical_file_map_v3.csv"), sep = "")
cat("\n - ", file.path(QC_DIR, "domain2_required_outputs_audit_after_fix_v3.csv"), "\n", sep = "")
