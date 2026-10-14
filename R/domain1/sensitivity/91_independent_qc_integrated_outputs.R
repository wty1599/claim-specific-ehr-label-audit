#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

root <- normalizePath(
  Sys.getenv("MIMIC_IV_ROOT", unset = Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())),
  winslash = "/", mustWork = TRUE
)
pkg <- file.path(root, "d1_current_rule_sensitivity_20260803")
out <- file.path(pkg, "integrated_reader_outputs_20260804")

files <- list(
  registry = file.path(out, "Table_G1_G7_canonical_reader_registry.csv"),
  domain_map = file.path(out, "Table_G1_G7_by_domain_reader_map.csv"),
  sa01 = file.path(out, "Table_D1_SA01_current_rule_reader_summary.csv"),
  sa05 = file.path(out, "Table_D1_SA05_current_rule_reader_summary.csv"),
  controls = file.path(out, "Table_current_formal_control_results_cross_domain.csv"),
  manifest = file.path(out, "CURRENT_RESULT_SOURCE_MANIFEST.csv"),
  zero_qc = file.path(pkg, "output_sa05_formal", "qc",
                      "SA05_zero_missingness_equivalence.csv"),
  sa01_qc = file.path(pkg, "output_sa01_formal", "qc", "critical_qc.csv"),
  sa05_qc = file.path(pkg, "output_sa05_formal", "qc", "critical_qc.csv")
)
missing <- names(files)[!vapply(files, file.exists, logical(1))]
if (length(missing)) {
  stop("Missing QC input(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

registry <- fread(files$registry)
domain_map <- fread(files$domain_map)
sa01 <- fread(files$sa01)
sa05 <- fread(files$sa05)
controls <- fread(files$controls)
manifest <- fread(files$manifest)
zero_qc <- fread(files$zero_qc)
sa01_qc <- fread(files$sa01_qc)
sa05_qc <- fread(files$sa05_qc)

qc <- list()
add_qc <- function(item, pass, observed, expected) {
  qc[[length(qc) + 1L]] <<- data.table(
    item = item, pass = isTRUE(pass),
    observed = as.character(observed), expected = as.character(expected)
  )
}

add_qc("canonical_registry_exact_G1_G7",
       identical(registry$canonical_id, paste0("G", 1:7)),
       paste(registry$canonical_id, collapse = ";"), paste0("G", 1:7) |> paste(collapse = ";"))
add_qc("domain_map_row_count", nrow(domain_map) == 15L, nrow(domain_map), 15L)
add_qc("domain_map_source_ids_no_absolute_paths",
       all(!grepl("[:\\\\/]", domain_map$current_source_id)),
       sum(grepl("[:\\\\/]", domain_map$current_source_id)), 0L)
add_qc("G2_excluded_from_current_D1_verdict",
       domain_map[domain == "D1 Discreteness" & reader_scenario == "G2",
                  current_formal_use] == "Excluded from current fixed-K2 verdict",
       domain_map[domain == "D1 Discreteness" & reader_scenario == "G2",
                  current_formal_use],
       "Excluded from current fixed-K2 verdict")
add_qc("G7_is_current_D1_positive_control",
       domain_map[domain == "D1 Discreteness" & reader_scenario == "G7",
                  grepl("Included", current_formal_use)] %in% TRUE,
       domain_map[domain == "D1 Discreteness" & reader_scenario == "G7",
                  current_formal_use], "Included")

add_qc("SA01_row_count", nrow(sa01) == 8L, nrow(sa01), 8L)
add_qc("SA01_all_valid_no_failures",
       sum(sa01$valid_n) == 800L && sum(sa01$failed_n) == 0L,
       paste(sum(sa01$valid_n), sum(sa01$failed_n), sep = "/"), "800/0")
add_qc("SA01_G7_discrete_400_of_400",
       sa01[reader_scenario == "G7", sum(discrete_n)] == 400L,
       sa01[reader_scenario == "G7", sum(discrete_n)], 400L)
add_qc("SA01_G1_false_discrete_zero",
       sa01[reader_scenario == "G1", sum(discrete_n)] == 0L,
       sa01[reader_scenario == "G1", sum(discrete_n)], 0L)
add_qc("SA01_G1_inconclusive_preserved",
       sa01[reader_scenario == "G1", sum(inconclusive_n)] == 153L,
       sa01[reader_scenario == "G1", sum(inconclusive_n)], 153L)

add_qc("SA05_row_count", nrow(sa05) == 18L, nrow(sa05), 18L)
add_qc("SA05_all_valid_no_failures",
       sum(sa05$valid_n) == 900L && sum(sa05$failed_n) == 0L,
       paste(sum(sa05$valid_n), sum(sa05$failed_n), sep = "/"), "900/0")
add_qc("SA05_G7_discrete_450_of_450",
       sa05[reader_scenario == "G7", sum(discrete_n)] == 450L,
       sa05[reader_scenario == "G7", sum(discrete_n)], 450L)
add_qc("SA05_G1_false_discrete_zero",
       sa05[reader_scenario == "G1", sum(discrete_n)] == 0L,
       sa05[reader_scenario == "G1", sum(discrete_n)], 0L)
add_qc("SA05_G1_inconclusive_preserved",
       sa05[reader_scenario == "G1", sum(inconclusive_n)] == 126L,
       sa05[reader_scenario == "G1", sum(inconclusive_n)], 126L)
add_qc("SA05_zero_eight_metrics_exact",
       nrow(zero_qc) == 800L && all(zero_qc$exact_metric_equivalence),
       paste(nrow(zero_qc), sum(zero_qc$exact_metric_equivalence), sep = "/"),
       "800/800")
add_qc("SA01_critical_qc_all_pass", all(sa01_qc$pass),
       sum(sa01_qc$pass), nrow(sa01_qc))
add_qc("SA05_critical_qc_all_pass", all(sa05_qc$pass),
       sum(sa05_qc$pass), nrow(sa05_qc))

d3_rows <- controls[domain == "D3"]
d4_rows <- controls[domain == "D4"]
add_qc("D3_current_controls_3x100_correct",
       nrow(d3_rows) == 3L && sum(d3_rows$n) == 300L &&
         sum(d3_rows$strict_expected_state) == 300L && sum(d3_rows$failed) == 0L,
       paste(nrow(d3_rows), sum(d3_rows$n), sum(d3_rows$strict_expected_state),
             sum(d3_rows$failed), sep = "/"), "3/300/300/0")
add_qc("D4_current_controls_2x300_correct",
       nrow(d4_rows) == 2L && sum(d4_rows$n) == 600L &&
         sum(d4_rows$strict_expected_state) == 600L && sum(d4_rows$failed) == 0L,
       paste(nrow(d4_rows), sum(d4_rows$n), sum(d4_rows$strict_expected_state),
             sum(d4_rows$failed), sep = "/"), "2/600/600/0")

hash_now <- vapply(manifest$absolute_path, function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}, character(1))
add_qc("source_manifest_hashes_current",
       all(hash_now == manifest$sha256), sum(hash_now == manifest$sha256), nrow(manifest))
add_qc("source_manifest_excludes_invalid_archive",
       !any(grepl("archive_invalid", manifest$absolute_path, ignore.case = TRUE)),
       sum(grepl("archive_invalid", manifest$absolute_path, ignore.case = TRUE)), 0L)

qc_dt <- rbindlist(qc)
fwrite(qc_dt, file.path(out, "INDEPENDENT_INTEGRATION_QC.csv"))

report <- c(
  "# Independent integration QC",
  "",
  paste0("Overall pass: ", all(qc_dt$pass)),
  paste0("Checks passed: ", sum(qc_dt$pass), "/", nrow(qc_dt)),
  "",
  "The current D1 results preserve INCONCLUSIVE as a separate state. G1 had",
  "no false DISCRETE_EVIDENCE alert in SA-01 or SA-05, but INCONCLUSIVE",
  "repeats are not reclassified as correct negatives. G7 is the aligned",
  "fixed-K2 positive control. G2 remains a legacy three-component recovery",
  "diagnostic and is excluded from the current D1 verdict.",
  "",
  "The SA-05 outputs archived under",
  "archive_invalid_zero_missingness_equivalence_20260804 are provenance only",
  "and are absent from the source manifest.",
  "",
  "Runtime note: R 4.6.0 was used. Several installed packages report that",
  "they were built under R 4.6.1; no package-load or computation failure was",
  "observed. This patch-level build warning should be retained in provenance."
)
writeLines(report, file.path(out, "INDEPENDENT_INTEGRATION_QC_REPORT.md"))

if (!all(qc_dt$pass)) {
  stop("Independent integration QC failed.", call. = FALSE)
}
cat("Independent integration QC passed: ", nrow(qc_dt), " checks.\n", sep = "")
