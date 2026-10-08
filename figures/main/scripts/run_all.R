#!/usr/bin/env Rscript
# Plotting-only entry point; restricted sources and outputs remain external.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg), winslash = "/"))
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(script_dir, "../../.."), winslash = "/"))
source(file.path(repo, "R/common/02_figure_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript run_all.R [authorized_figure4_aggregate_dir]")
figure_paths <- ehr_audit_figure_paths(repo,
  figure4_source = if (length(args)) args[[1]] else NULL, create_outputs = TRUE)
scripts <- file.path(repo, c("figures/main/scripts/Figure1.R",
  "figures/main/scripts/Figure2.R", "figures/main/scripts/Figure3.R",
  "figures/D3/scripts/Figure4.R", "figures/main/scripts/Figure5.R"))
stopifnot(all(file.exists(scripts)))
before <- ehr_audit_figure_snapshot(figure_paths)
run_id <- paste(format(Sys.time(), "%Y%m%dT%H%M%OS6"), Sys.getpid(), sep = "-")
Sys.setenv(EHR_AUDIT_REPO_ROOT = repo, EHR_AUDIT_FIGURE_INPUT_RUN_ID = run_id)
executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (script in scripts) {
  message("Rendering ", basename(script))
  script_args <- if (identical(basename(script), "Figure4.R"))
    shQuote(figure_paths$figure4_source) else character()
  status <- system2(executable, c("--vanilla", shQuote(script), script_args), wait = TRUE)
  if (status != 0L) stop("Figure redraw failed: ", script)
}
manifests <- c(file.path(figure_paths$main_qa, paste0("Figure", c(1, 2, 3, 5), "_input_fingerprints.csv")),
               file.path(figure_paths$figure4_qa, "Figure4_input_fingerprints.csv"))
consumed <- ehr_audit_verify_figure_inputs(figure_paths, before, manifests, run_id)
write.csv(consumed, file.path(figure_paths$main_qa, "main_figure_consumed_inputs.csv"), row.names = FALSE)
message("Five main figures rendered; actual consumed inputs matched the unchanged external sources.")
