# Shared external paths and consumed-input records for main-figure rendering.
source(file.path(repo, "R/common/01_assert_external_output.R"))

ehr_audit_figure_paths <- function(repo, figure4_source = NULL,
                                   require_sources = c("main", "figure4"),
                                   create_outputs = FALSE) {
  repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
  output <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT", unset = "")
  ehr_audit_assert_external_output(output, repo)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  select_source <- function(variable, default, required, explicit = NULL) {
    if (!required) return(default)
    configured <- if (is.null(explicit)) Sys.getenv(variable, unset = "") else explicit
    path <- if (nzchar(configured)) configured else default
    if (required || nzchar(configured) || dir.exists(path)) {
      ehr_audit_assert_external_output(path, repo)
      path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    }
    path
  }
  main_root <- file.path(output, "main_figures")
  d3_root <- file.path(output, "domain3_source_correction/figure4")
  paths <- list(repo = repo,
    main_source = select_source("EHR_AUDIT_MAIN_FIGURE_SOURCE_ROOT",
      file.path(main_root, "source_data"), "main" %in% require_sources),
    figure4_source = select_source("EHR_AUDIT_FIGURE4_SOURCE_ROOT",
      file.path(output, "domain3_raw_temperature_completion/figure4/source_data"),
      "figure4" %in% require_sources, figure4_source),
    main_output = file.path(main_root, "figures"), main_qa = file.path(main_root, "qa"),
    figure4_output = file.path(d3_root, "figures"), figure4_qa = file.path(d3_root, "qa"))
  # Reject overlapping source/output trees before creating any output directory.
  for (input in unlist(paths[c("main_source", "figure4_source")])) {
    for (destination in unlist(paths[c("main_output", "main_qa", "figure4_output", "figure4_qa")])) {
      a <- tolower(input); b <- tolower(destination)
      if (identical(a, b) || startsWith(a, paste0(b, "/")) || startsWith(b, paste0(a, "/"))) {
        stop("Figure source and generated-output directories must not overlap.")
      }
    }
  }
  if (create_outputs) {
    directories <- c(if ("main" %in% require_sources) unlist(paths[c("main_output", "main_qa")]),
                     if ("figure4" %in% require_sources) unlist(paths[c("figure4_output", "figure4_qa")]))
    for (directory in directories) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      ehr_audit_assert_external_output(directory, repo)
    }
  }
  paths
}

ehr_audit_figure_reader <- function(source_root, manifest_path, repo) {
  ehr_audit_assert_external_output(source_root, repo)
  ehr_audit_assert_external_output(dirname(manifest_path), repo)
  source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
  records <- data.frame(input_path = character(), md5 = character(), run_id = character())
  utils::write.csv(records, manifest_path, row.names = FALSE)
  function(path, ...) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    if (!startsWith(tolower(path), paste0(tolower(source_root), "/")) ||
        isTRUE(file.info(path)$isdir)) stop("Figure input must be a file within the selected source root.")
    before <- unname(tools::md5sum(path))
    result <- utils::read.csv(path, ...)
    after <- unname(tools::md5sum(path))
    if (is.na(before) || !identical(before, after)) stop("Figure source changed while being read.")
    records <<- records[records$input_path != path, , drop = FALSE]
    records <<- rbind(records, data.frame(input_path = path, md5 = before,
      run_id = Sys.getenv("EHR_AUDIT_FIGURE_INPUT_RUN_ID", unset = "standalone")))
    utils::write.csv(records, manifest_path, row.names = FALSE)
    result
  }
}

ehr_audit_figure_snapshot <- function(paths) {
  inputs <- unique(c(list.files(paths$main_source, recursive = TRUE, full.names = TRUE),
                     list.files(paths$figure4_source, recursive = TRUE, full.names = TRUE),
                     file.path(paths$repo, "figures/main/figure_spec.json")))
  inputs <- inputs[!is.na(file.info(inputs)$isdir) & !file.info(inputs)$isdir]
  inputs <- sort(normalizePath(inputs, winslash = "/", mustWork = TRUE))
  tools::md5sum(inputs)
}

ehr_audit_verify_figure_inputs <- function(paths, before, manifests, run_id) {
  if (!identical(before, ehr_audit_figure_snapshot(paths))) {
    stop("Selected figure sources or their membership changed during rendering.")
  }
  consumed <- lapply(manifests, function(manifest) {
    if (!file.exists(manifest)) stop("Renderer did not record its consumed inputs: ", manifest)
    records <- utils::read.csv(manifest, stringsAsFactors = FALSE)
    if (!all(c("input_path", "md5", "run_id") %in% names(records)) ||
        !nrow(records) || !all(records$run_id == run_id)) {
      stop("Missing or stale consumed-input records: ", manifest)
    }
    expected <- unname(before[records$input_path])
    if (anyNA(expected) || !identical(expected, records$md5)) {
      stop("Renderer consumed inputs not covered by the selected-source fingerprints.")
    }
    records
  })
  do.call(rbind, consumed)
}
