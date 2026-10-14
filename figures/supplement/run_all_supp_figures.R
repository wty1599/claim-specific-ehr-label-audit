#!/usr/bin/env Rscript
# Plotting-only entry point. No fitting, imputation, simulation, or resampling.
options(stringsAsFactors = FALSE)
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
root <- dirname(normalizePath(sub("^--file=", "", args), winslash = "/"))
Sys.setenv(LC_ALL = "C")
required <- c("ggplot2","patchwork","data.table","jsonlite","svglite","digest")
if (!all(vapply(required, requireNamespace, logical(1), quietly=TRUE))) stop("Missing R plotting dependency.")
scripts <- c("SF1_mice_diagnostics_v4_countpanels.R","SF2_phenotype_harmonization.R",
 "SF3_D1_lineage_fisher_ecdf_v5.R","SF4_domain2_label_lineage_v4_countpanels.R",
 "SF5_external_validation_sensitivity.R","SF6_domain4_gate_overlap_v4_countpanels.R",
 "SF7_D1_split_scenario_ids_v4_countpanels.R","SF8_oracle_accuracy_dose_response.R",
 "SF9_subtype_separation_reproducibility.R","SF10_missingness_scenario_ids.R",
 "SF11_algorithm_scenario_ids.R","SF12_D3_scenario_ids.R")
stems <- c("S1_MICE_diagnostics_v5_submission_main", "S2_profile_harmonization_v4_submission",
 "S3_D1_lineage_ecdf_v5","S4_domain2_label_lineage_v5_submission_main",
 "S5_surrogate_application_sensitivity_v4_submission","S6_D4_support_v5_submission_main",
 "S7_D1_split_sensitivity_v5_submission_main","S8_oracle_accuracy_dose_response_v4_submission",
 "S9_separation_reproducibility_v4_submission","S10_missingness_robustness_v4_submission",
 "S11_algorithm_robustness_v4_submission","S12_D3_reconstructed_generator_validation_v4_submission")
stems <- paste0("Supplementary_Figure_", stems)
data_files <- list.files(file.path(root,"data"), recursive=TRUE, full.names=TRUE)
hash <- function(x) vapply(x, digest::digest, character(1),algo="sha256",file=TRUE)
before <- hash(data_files)
dir.create(file.path(root,"logs"),showWarnings=FALSE)
dir.create(file.path(root,"qa"),showWarnings=FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")
only <- commandArgs(TRUE)
ids <- if(length(only)) as.integer(only) else 1:12
stopifnot(all(ids %in% 1:12))
for(i in ids) {
  message("Rendering Supplementary Figure S",i)
  script <- file.path(root,"R","figures","supplement",scripts[i])
  log <- file.path(root,"logs",paste0("S",i,".log"))
  status <- system2(rscript,c("--vanilla",shQuote(script)),stdout=log,stderr=log)
  if(status != 0) stop("S",i," failed; inspect ",log)
}
after <- hash(data_files)
stopifnot(identical(before,after))
data.table::fwrite(data.table::data.table(path=sub(paste0(root,"/"),"",data_files,fixed=TRUE),
  before=before,after=after,unchanged=before==after),file.path(root,"qa","locked_input_hashes.csv"))
manifest <- list()
for(i in 1:12) for(ext in c("pdf","svg","png")) {
  src <- file.path(root,"figures","supplement",paste0(stems[i],".",ext))
  if(!file.exists(src)) next
  folder <- c(pdf="PDF_vector",svg="SVG_vector",png="PNG_600dpi")[[ext]]
  dir.create(file.path(root,"Figures",folder),recursive=TRUE,showWarnings=FALSE)
  dst <- file.path(root,"Figures",folder,paste0("Supplementary_Figure_S",i,".",ext))
  stopifnot(file.copy(src,dst,overwrite=TRUE))
  snap <- readRDS(file.path(root,"qa","snapshots",paste0(stems[i],".rds")))
  manifest[[length(manifest)+1L]] <- data.table::data.table(figure=paste0("S",i),format=ext,
    width_mm=snap$width_mm,height_mm=snap$height_mm,dpi=ifelse(ext=="png",600,NA),
    source_script=paste0("R/figures/supplement/",scripts[i]),
    file=sub(paste0(root,"/"),"",dst,fixed=TRUE),sha256=hash(dst))
}
data.table::fwrite(data.table::rbindlist(manifest),file.path(root,"supp_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()),file.path(root,"sessionInfo.txt"))
message("Plot-only run completed; locked input files are unchanged.")
