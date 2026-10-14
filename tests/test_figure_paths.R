args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
repo <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
audit <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/02_figure_paths.R"))
ehr_audit_assert_external_output(audit, repo)
fixture <- tempfile("figure_paths_", tmpdir = audit)
dir.create(fixture)
main <- file.path(fixture, "main_sources")
raw <- file.path(fixture, "raw_sources")
unused <- file.path(fixture, "unused_sources")
output <- file.path(fixture, "generated")
for (path in c(main, raw, unused, output)) dir.create(path)
Sys.setenv(EHR_AUDIT_REPO_ROOT = repo, EHR_AUDIT_OUTPUT_ROOT = output,
  EHR_AUDIT_MAIN_FIGURE_SOURCE_ROOT = main, EHR_AUDIT_FIGURE4_SOURCE_ROOT = unused,
  EHR_AUDIT_FIGURE_INPUT_RUN_ID = "artificial-test")
checks <- character()
check <- function(name, condition) {
  stopifnot(isTRUE(condition)); checks <<- c(checks, name)
}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
paths <- ehr_audit_figure_paths(repo, figure4_source = raw, create_outputs = TRUE)
check("explicit Figure4 source overrides environment", identical(paths$figure4_source,
  normalizePath(raw, winslash = "/")))
check("environment fallback is shared", identical(ehr_audit_figure_paths(repo)$figure4_source,
  normalizePath(unused, winslash = "/")))
check("repository source rejected", fails(ehr_audit_figure_paths(repo, figure4_source = repo)))
check("overlapping source/output rejected", fails(ehr_audit_figure_paths(repo, figure4_source = output)))
Sys.setenv(EHR_AUDIT_OUTPUT_ROOT = repo)
check("repository output rejected before creation", fails(ehr_audit_figure_paths(repo, figure4_source = raw,
  create_outputs = TRUE)))
Sys.setenv(EHR_AUDIT_OUTPUT_ROOT = output)
sample <- data.frame(a = c(1, 2), b = c("x", "y"))
write.csv(sample, file.path(raw, "fixture.csv"), row.names = FALSE)
write.csv(sample, file.path(main, "fixture.csv"), row.names = FALSE)
write.csv(sample, file.path(unused, "fixture.csv"), row.names = FALSE)
manifest <- file.path(paths$figure4_qa, "fixture_inputs.csv")
reader <- ehr_audit_figure_reader(raw, manifest, repo)
before <- ehr_audit_figure_snapshot(paths)
check("reader preserves CSV values and options", identical(reader(file.path(raw, "fixture.csv"),
  check.names = FALSE, stringsAsFactors = FALSE), read.csv(file.path(raw, "fixture.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)))
check("consumed files covered by selected snapshot", nrow(ehr_audit_verify_figure_inputs(paths, before,
  manifest, "artificial-test")) == 1L)
check("unselected input rejected", fails(reader(file.path(unused, "fixture.csv"))))
check("stale consumption record rejected", fails(ehr_audit_verify_figure_inputs(paths, before,
  manifest, "different-run")))
rogue <- ehr_audit_figure_reader(unused, file.path(paths$figure4_qa, "rogue.csv"), repo)
invisible(rogue(file.path(unused, "fixture.csv")))
check("fingerprint cannot certify an unselected consumed file", fails(ehr_audit_verify_figure_inputs(paths,
  before, file.path(paths$figure4_qa, "rogue.csv"), "artificial-test")))
write.csv(data.frame(a = 3, b = "changed"), file.path(raw, "fixture.csv"), row.names = FALSE)
check("changed source detected", fails(ehr_audit_verify_figure_inputs(paths, before, manifest, "artificial-test")))
write.csv(sample, file.path(raw, "fixture.csv"), row.names = FALSE)
write.csv(sample, file.path(raw, "added.csv"), row.names = FALSE)
check("changed source membership detected", fails(ehr_audit_verify_figure_inputs(paths, before, manifest,
  "artificial-test")))

# Execute entry setup only, stopping before any plotting or clinical CSV read.
setup <- function(relative, anchor, trailing) {
  file <- file.path(repo, relative)
  e <- parse(file, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) {
    if (trailingOnly) trailing else paste0("--file=", file)
  }
  for (x in e) {
    if (anchor(x)) break
    eval(x, env)
  }
  env
}
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name(name))
renderer <- setup("figures/D3/scripts/Figure4.R", function(x) is_assignment(x, "fidelity"), raw)
entry <- setup("figures/main/scripts/run_all.R", function(x) is.call(x) && identical(x[[1]], as.name("for")), raw)
check("entry and Figure4 consume the identical argument root", identical(renderer$src,
  entry$figure_paths$figure4_source))
common <- setup("figures/main/scripts/common.R", function(x) is.call(x) &&
  identical(x[[1]], as.name("suppressPackageStartupMessages")), character())
check("common sources and outputs are external", identical(common$DATA, normalizePath(main, winslash = "/")) &&
  startsWith(common$OUT, normalizePath(output, winslash = "/")) &&
  startsWith(common$QA, normalizePath(output, winslash = "/")))
Sys.setenv(EHR_AUDIT_FIGURE4_SOURCE_ROOT = repo)
main_only <- ehr_audit_figure_paths(repo, require_sources = "main", create_outputs = FALSE)
check("unused Figure4 environment variable does not affect main-only setup",
      identical(main_only$main_source, normalizePath(main, winslash = "/")))
Sys.unsetenv("EHR_AUDIT_FIGURE4_SOURCE_ROOT")
for (relative in c("figures/main/scripts/Figure2.R", "figures/main/scripts/Figure3.R")) {
  definition <- Filter(function(x) is_assignment(x, "src"), as.list(parse(file.path(repo, relative))))[[1]]
  env <- new.env(parent = globalenv()); env$DATA <- main
  env$figure_reader <- function(path, ...) read.csv(path, ...)
  eval(definition, env)
  dir.create(file.path(main, "fig23"), showWarnings = FALSE)
  write.csv(sample, file.path(main, "fig23/fixture.csv"), row.names = FALSE)
  check(paste(basename(relative), "reader options preserved"), identical(env$src("fixture.csv"),
    read.csv(file.path(main, "fig23/fixture.csv"), check.names = FALSE,
             stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")))
}
write.csv(data.frame(check = checks, passed = TRUE), file.path(audit, "figure_path_tests.csv"), row.names = FALSE)
cat(length(checks), "figure path/consumption checks passed; artificial files only, no rendering.\n")
