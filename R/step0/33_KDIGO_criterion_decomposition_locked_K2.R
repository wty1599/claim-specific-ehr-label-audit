## KDIGO 0-24 h criterion-source decomposition for locked K=2 labels.
## Descriptive post-processing only: no cohort rebuild, imputation, or clustering.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
ANALYSIS_ROOT <- file.path(
  PROJECT_ROOT,
  "kdigo_criterion_decomposition_20260725"
)

INPUT <- list(
  components = file.path(
    ANALYSIS_ROOT,
    "data_extract",
    "kdigo_component_sources_0_24h.csv"
  ),
  labels = file.path(
    PROJECT_ROOT,
    "output",
    "model",
    "labels_primary_mice.rds"
  ),
  authority_sql = file.path(
    PROJECT_ROOT,
    "sql",
    "build_baseline_covars_v2.sql"
  ),
  export_sql = file.path(
    ANALYSIS_ROOT,
    "scripts",
    "51a_export_kdigo_component_sources_0_24h.sql"
  ),
  postprocess_r = file.path(
    ANALYSIS_ROOT,
    "scripts",
    "51b_kdigo_criterion_decomposition_locked_k2.R"
  )
)

DIRS <- list(
  tables = file.path(ANALYSIS_ROOT, "tables"),
  logs = file.path(ANALYSIS_ROOT, "logs"),
  provenance = file.path(ANALYSIS_ROOT, "provenance")
)
lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE)

missing_input <- names(INPUT)[!file.exists(unlist(INPUT))]
if (length(missing_input)) {
  stop(
    "Missing input(s): ",
    paste(missing_input, collapse = ", "),
    ". Run 51a_export_kdigo_component_sources_0_24h.sql first."
  )
}

EXPECTED <- list(
  n = 20049L,
  c1 = 3992L,
  c2 = 16057L,
  stage_counts = c(`0` = 9056L, `1` = 4108L, `2` = 6030L, `3` = 855L)
)

write_csv_atomic <- function(x, path) {
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path)
  }
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (!file.rename(tmp, path)) {
    stop("Could not write: ", path)
  }
}

write_lines_atomic <- function(x, path) {
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path)
  }
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) {
    stop("Could not write: ", path)
  }
}

wilson_ci <- function(x, n, conf_level = 0.95) {
  if (n == 0L) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(1 - (1 - conf_level) / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) /
    denominator
  c(low = center - half, high = center + half)
}

summarize_pattern <- function(dt, pattern_col, pattern_levels) {
  counts <- dt[, .(n = .N), by = c(pattern_col, "cluster")]
  setnames(counts, pattern_col, "criterion_pattern")

  grid <- CJ(
    criterion_pattern = pattern_levels,
    cluster = c("C1 higher-risk", "C2 lower-risk"),
    unique = TRUE
  )
  counts <- merge(
    grid,
    counts,
    by = c("criterion_pattern", "cluster"),
    all.x = TRUE,
    sort = FALSE
  )
  counts[is.na(n), n := 0L]
  counts[, pattern_total := sum(n), by = criterion_pattern]
  counts[, cluster_total := sum(n), by = cluster]
  counts[, proportion_cluster_given_pattern := n / pattern_total]
  counts[, proportion_pattern_given_cluster := n / cluster_total]

  ci <- t(mapply(
    wilson_ci,
    x = counts$n,
    n = counts$pattern_total
  ))
  counts[, `:=`(
    cluster_given_pattern_wilson_low = ci[, "low"],
    cluster_given_pattern_wilson_high = ci[, "high"]
  )]

  counts[, criterion_pattern := factor(
    criterion_pattern,
    levels = pattern_levels
  )]
  setorder(counts, criterion_pattern, cluster)
  counts[, criterion_pattern := as.character(criterion_pattern)]
  counts[]
}

components <- fread(INPUT$components, showProgress = FALSE)
labels <- as.data.table(readRDS(INPUT$labels))

required_components <- c(
  "stay_id",
  "stage_creat_0_24h",
  "stage_uo_0_24h",
  "rrt_0_24h",
  "aki_stage_0_24h_existing",
  "aki_stage_0_24h_reconstructed"
)
required_labels <- c("stay_id", "cluster_k2")

if (!all(required_components %in% names(components))) {
  stop(
    "Component export is missing columns: ",
    paste(setdiff(required_components, names(components)), collapse = ", ")
  )
}
if (!all(required_labels %in% names(labels))) {
  stop(
    "Locked labels are missing columns: ",
    paste(setdiff(required_labels, names(labels)), collapse = ", ")
  )
}
if (anyDuplicated(components$stay_id)) {
  stop("Component export contains duplicate stay_id values.")
}
if (anyDuplicated(labels$stay_id)) {
  stop("Locked labels contain duplicate stay_id values.")
}

matched <- components[match(labels$stay_id, components$stay_id)]
if (anyNA(matched$stay_id)) {
  stop("At least one locked stay_id is absent from the component export.")
}
matched[, cluster_k2 := as.integer(labels$cluster_k2)]

integer_fields <- setdiff(required_components, "stay_id")
matched[, (integer_fields) := lapply(.SD, as.integer), .SDcols = integer_fields]

if (anyNA(matched[, ..integer_fields])) {
  stop("Missing component stage values remain after matching.")
}
if (!all(matched$stage_creat_0_24h %in% 0:3)) {
  stop("Invalid creatinine component stage.")
}
if (!all(matched$stage_uo_0_24h %in% 0:3)) {
  stop("Invalid urine-output component stage.")
}
if (!all(matched$rrt_0_24h %in% 0:1)) {
  stop("Invalid RRT flag.")
}
if (!all(matched$cluster_k2 %in% 1:2)) {
  stop("Invalid locked K=2 label.")
}

matched[, reconstructed_in_r := pmax(
  stage_creat_0_24h,
  stage_uo_0_24h,
  fifelse(rrt_0_24h == 1L, 3L, 0L)
)]

observed_stage_counts <- table(factor(
  matched$aki_stage_0_24h_existing,
  levels = 0:3
))

qc <- data.table(
  check = c(
    "locked_rows_matched",
    "locked_C1_count",
    "locked_C2_count",
    "sql_reconstructed_equals_existing",
    "R_reconstructed_equals_existing",
    "stage_0_count",
    "stage_1_count",
    "stage_2_count",
    "stage_3_count",
    "component_values_complete"
  ),
  observed = c(
    nrow(matched),
    sum(matched$cluster_k2 == 1L),
    sum(matched$cluster_k2 == 2L),
    sum(
      matched$aki_stage_0_24h_reconstructed ==
        matched$aki_stage_0_24h_existing
    ),
    sum(matched$reconstructed_in_r == matched$aki_stage_0_24h_existing),
    as.integer(observed_stage_counts[1L]),
    as.integer(observed_stage_counts[2L]),
    as.integer(observed_stage_counts[3L]),
    as.integer(observed_stage_counts[4L]),
    sum(complete.cases(matched[, ..integer_fields]))
  ),
  required = c(
    EXPECTED$n,
    EXPECTED$c1,
    EXPECTED$c2,
    EXPECTED$n,
    EXPECTED$n,
    unname(EXPECTED$stage_counts),
    EXPECTED$n
  )
)
qc[, pass := observed == required]

if (!all(qc$pass)) {
  failure_path <- file.path(
    DIRS$logs,
    "51b_QC_FAILURE_no_formal_tables_generated.csv"
  )
  fwrite(qc, failure_path, bom = TRUE)
  stop(
    "KDIGO component-source QC failed. No formal decomposition tables were ",
    "generated. See: ",
    failure_path
  )
}

matched[, `:=`(
  cluster = fifelse(
    cluster_k2 == 1L,
    "C1 higher-risk",
    "C2 lower-risk"
  ),
  creat_positive = stage_creat_0_24h > 0L,
  uo_positive = stage_uo_0_24h > 0L,
  rrt_positive = rrt_0_24h == 1L
)]

matched[, criterion4 := fcase(
  creat_positive & uo_positive, "Creatinine + urine output",
  creat_positive, "Creatinine only",
  uo_positive, "Urine output only",
  default = "Neither creatinine nor urine output"
)]

matched[, criterion8 := fcase(
  creat_positive & uo_positive & rrt_positive,
  "Creatinine + urine output + RRT",
  creat_positive & uo_positive,
  "Creatinine + urine output",
  creat_positive & rrt_positive,
  "Creatinine + RRT",
  uo_positive & rrt_positive,
  "Urine output + RRT",
  creat_positive,
  "Creatinine only",
  uo_positive,
  "Urine output only",
  rrt_positive,
  "RRT only",
  default = "No KDIGO trigger"
)]

matched[, source_at_combined_max := mapply(
  FUN = function(creat_stage, uo_stage, rrt_flag, combined_stage) {
    if (combined_stage == 0L) return("No KDIGO trigger")
    rrt_stage <- if (rrt_flag == 1L) 3L else 0L
    source <- character()
    if (creat_stage == combined_stage) source <- c(source, "Creatinine")
    if (uo_stage == combined_stage) source <- c(source, "Urine output")
    if (rrt_stage == combined_stage) source <- c(source, "RRT")
    paste(source, collapse = " + ")
  },
  stage_creat_0_24h,
  stage_uo_0_24h,
  rrt_0_24h,
  aki_stage_0_24h_existing,
  USE.NAMES = FALSE
)]

criterion4_levels <- c(
  "Neither creatinine nor urine output",
  "Creatinine only",
  "Urine output only",
  "Creatinine + urine output"
)

criterion8_levels <- c(
  "No KDIGO trigger",
  "Creatinine only",
  "Urine output only",
  "Creatinine + urine output",
  "RRT only",
  "Creatinine + RRT",
  "Urine output + RRT",
  "Creatinine + urine output + RRT"
)

table4 <- summarize_pattern(
  matched,
  "criterion4",
  criterion4_levels
)
table8 <- summarize_pattern(matched, "criterion8", criterion8_levels)

component_long <- melt(
  matched,
  id.vars = c("stay_id", "cluster"),
  measure.vars = c("stage_creat_0_24h", "stage_uo_0_24h"),
  variable.name = "component",
  value.name = "component_stage"
)
component_long[, component := fifelse(
  component == "stage_creat_0_24h",
  "Creatinine",
  "Urine output"
)]
component_stage <- component_long[, .(n = .N), by = .(
  component,
  component_stage,
  cluster
)]
component_stage[, stage_total := sum(n), by = .(component, component_stage)]
component_stage[, cluster_total := sum(n), by = .(component, cluster)]
component_stage[, `:=`(
  proportion_cluster_given_stage = n / stage_total,
  proportion_stage_given_cluster = n / cluster_total
)]
setorder(component_stage, component, component_stage, cluster)

source_at_max <- summarize_pattern(
  matched,
  "source_at_combined_max",
  unique(matched$source_at_combined_max)
)

combined_stage_source <- matched[, .(n = .N), by = .(
  combined_stage = aki_stage_0_24h_existing,
  source_at_combined_max,
  cluster
)]
combined_stage_source[, source_total := sum(n), by = .(
  combined_stage,
  source_at_combined_max
)]
combined_stage_source[, combined_stage_total := sum(n), by = combined_stage]
combined_stage_source[, cluster_total := sum(n), by = cluster]
combined_stage_source[, `:=`(
  proportion_cluster_given_stage_source = n / source_total,
  proportion_source_given_combined_stage = n / combined_stage_total,
  proportion_stage_source_given_cluster = n / cluster_total
)]
stage_source_ci <- t(mapply(
  wilson_ci,
  x = combined_stage_source$n,
  n = combined_stage_source$source_total
))
combined_stage_source[, `:=`(
  cluster_given_stage_source_wilson_low = stage_source_ci[, "low"],
  cluster_given_stage_source_wilson_high = stage_source_ci[, "high"]
)]
setorder(
  combined_stage_source,
  combined_stage,
  source_at_combined_max,
  cluster
)

rrt_table <- matched[, .(n = .N), by = .(rrt_0_24h, cluster)]
rrt_table[, rrt_total := sum(n), by = rrt_0_24h]
rrt_table[, cluster_total := sum(n), by = cluster]
rrt_table[, `:=`(
  proportion_cluster_given_rrt = n / rrt_total,
  proportion_rrt_given_cluster = n / cluster_total
)]
setorder(rrt_table, rrt_0_24h, cluster)

output_paths <- c(
  file.path(DIRS$tables, "51b_Table_KDIGO_criterion4_by_locked_K2.csv"),
  file.path(DIRS$tables, "51b_Table_KDIGO_criterion8_by_locked_K2.csv"),
  file.path(DIRS$tables, "51b_Table_KDIGO_component_stage_by_locked_K2.csv"),
  file.path(DIRS$tables, "51b_Table_KDIGO_combined_stage_source_by_locked_K2.csv"),
  file.path(DIRS$tables, "51b_Table_KDIGO_source_at_max_by_locked_K2.csv"),
  file.path(DIRS$tables, "51b_Table_KDIGO_RRT_by_locked_K2.csv"),
  file.path(DIRS$logs, "51b_KDIGO_component_source_QC.csv"),
  file.path(DIRS$provenance, "51b_input_manifest.csv"),
  file.path(DIRS$logs, "51b_KDIGO_criterion_decomposition_report.md")
)
if (any(file.exists(output_paths))) {
  stop(
    "Refusing to overwrite existing outputs:\n",
    paste(output_paths[file.exists(output_paths)], collapse = "\n")
  )
}

write_csv_atomic(table4, output_paths[1L])
write_csv_atomic(table8, output_paths[2L])
write_csv_atomic(component_stage, output_paths[3L])
write_csv_atomic(combined_stage_source, output_paths[4L])
write_csv_atomic(source_at_max, output_paths[5L])
write_csv_atomic(rrt_table, output_paths[6L])
write_csv_atomic(qc, output_paths[7L])

input_info <- file.info(unlist(INPUT))
manifest <- data.table(
  input_name = names(INPUT),
  path = unname(unlist(INPUT)),
  size_bytes = input_info$size,
  modified_time = format(input_info$mtime, "%Y-%m-%d %H:%M:%S"),
  md5 = unname(tools::md5sum(unlist(INPUT)))
)
write_csv_atomic(manifest, output_paths[8L])

report <- c(
  "# KDIGO 0-24 h criterion-source decomposition",
  "",
  "Status: completed descriptive post-processing.",
  "",
  "## Analysis boundary",
  "",
  "- Existing component stages were re-exported; raw KDIGO was not recomputed.",
  "- The existing combined 0-24 h stage was not changed.",
  "- Locked MICE-imputation-1 plain-k-means labels were not changed.",
  "- No cohort, MICE, clustering, outcome, or causal model was rerun.",
  "",
  "## QC",
  "",
  "All prespecified alignment and stage-reconstruction checks passed.",
  "",
  "## Outputs",
  "",
  paste0("- ", basename(output_paths[1:8])),
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string)
)
write_lines_atomic(report, output_paths[9L])

cat("KDIGO criterion-source decomposition completed.\n")
cat("All QC checks passed.\n")
print(qc)
print(table4)
