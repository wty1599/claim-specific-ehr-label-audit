#!/usr/bin/env Rscript

# Read-only post-processing audit for the locked 29f, 29g, and 35 outputs.
# This script does not refit models or modify any source result.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

ROOT <- file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "domain1_lloyd_alignment_20260720")
SCRIPT_ARG <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_PATH <- if (length(SCRIPT_ARG) == 1L) {
  normalizePath(sub("^--file=", "", SCRIPT_ARG), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(ROOT, "scripts", "36_domain1_cross_stage_comparison_audit_lloyd.R"), winslash = "/", mustWork = TRUE)
}
DIR_TABLE <- file.path(ROOT, "tables")
DIR_PROV <- file.path(ROOT, "provenance")
DIR_LOG <- file.path(ROOT, "logs")
DIR_MANUSCRIPT <- file.path(ROOT, "manuscript_update")
DIR_TMP <- file.path(ROOT, "checkpoints", "36_cross_stage_postprocess_tmp")

dir.create(DIR_MANUSCRIPT, recursive = TRUE, showWarnings = FALSE)

OUT <- c(
  cross_stage = file.path(DIR_TABLE, "36_D1_cross_stage_comparison_summary_LLOYD.csv"),
  gate = file.path(DIR_TABLE, "36_D1_gate_split_figure_source_LLOYD.csv"),
  states = file.path(DIR_TABLE, "36_D1_state_split_figure_source_LLOYD.csv"),
  empirical = file.path(DIR_TABLE, "36_D1_empirical_split_figure_source_LLOYD.csv"),
  transitions = file.path(DIR_TABLE, "36_D1_transition_stability_figure_source_LLOYD.csv"),
  runtimes = file.path(DIR_PROV, "36_D1_module_runtime_provenance_LLOYD.csv"),
  inputs = file.path(DIR_PROV, "36_D1_comparison_input_manifest_LLOYD.csv"),
  qc = file.path(DIR_PROV, "36_D1_comparison_QC_LLOYD.csv"),
  replacement_csv = file.path(DIR_MANUSCRIPT, "D1_LLOYD_final_replacement_map_20260722.csv"),
  replacement_md = file.path(DIR_MANUSCRIPT, "D1_LLOYD_final_replacement_map_20260722.md"),
  manifest = file.path(DIR_PROV, "36_D1_output_sha256_manifest_LLOYD.csv"),
  completion = file.path(ROOT, "36_D1_cross_stage_postprocess_completed_LLOYD.ok")
)

if (any(file.exists(OUT))) {
  stop("One or more 36 outputs already exist. No files were overwritten.")
}

if (dir.exists(DIR_TMP)) {
  stop("Temporary directory already exists: ", DIR_TMP)
}
dir.create(DIR_TMP, recursive = TRUE, showWarnings = FALSE)
on.exit({
  if (dir.exists(DIR_TMP)) unlink(DIR_TMP, recursive = TRUE, force = TRUE)
}, add = TRUE)

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

parse_kv <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- x[grepl("=", x, fixed = TRUE)]
  keys <- sub("=.*$", "", x)
  vals <- sub("^[^=]*=", "", x)
  setNames(vals, keys)
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  fread(path, na.strings = c("NA", ""))
}

verify_manifest <- function(path) {
  m <- read_csv(path)
  if (!all(c("path", "sha256") %in% names(m))) {
    stop("Manifest lacks path/sha256 columns: ", path)
  }
  m[, exists := file.exists(path)]
  m[, observed_sha256 := ifelse(exists, vapply(path, sha256_file, character(1)), NA_character_)]
  m[, hash_match := exists & observed_sha256 == sha256]
  if (!all(m$hash_match)) {
    print(m[!hash_match])
    stop("Locked output manifest verification failed: ", path)
  }
  m
}

wilson <- function(x, n, z = qnorm(0.975)) {
  if (n <= 0) return(c(low = NA_real_, high = NA_real_))
  p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, ctr - half), high = min(1, ctr + half))
}

qc_rows <- list()
add_qc <- function(check, pass, detail) {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table(
    check = check,
    pass = isTRUE(pass),
    detail = as.character(detail)
  )
  if (!isTRUE(pass)) stop("QC failed [", check, "]: ", detail)
}

completion_paths <- c(
  `29f` = file.path(ROOT, "29f_run_completed_LLOYD.ok"),
  `29g` = file.path(ROOT, "29g_run_completed_LLOYD.ok"),
  `35` = file.path(ROOT, "35_WP9_run_completed_LLOYD.ok")
)
manifest_paths <- c(
  `29f` = file.path(DIR_PROV, "29f_output_sha256_manifest_LLOYD.csv"),
  `29g` = file.path(DIR_PROV, "29g_output_sha256_manifest_LLOYD.csv"),
  `35` = file.path(DIR_PROV, "35_WP9_output_sha256_manifest_LLOYD.csv")
)

for (nm in names(completion_paths)) {
  add_qc(paste0(nm, "_completion_exists"), file.exists(completion_paths[[nm]]), completion_paths[[nm]])
  kv <- parse_kv(completion_paths[[nm]])
  add_qc(paste0(nm, "_formal_mode"), identical(unname(kv[["run_mode"]]), "formal"), kv[["run_mode"]])
  add_qc(
    paste0(nm, "_manifest_hash_in_completion"),
    identical(unname(kv[["output_manifest_sha256"]]), unname(sha256_file(manifest_paths[[nm]]))),
    sha256_file(manifest_paths[[nm]])
  )
  verified <- verify_manifest(manifest_paths[[nm]])
  add_qc(paste0(nm, "_all_locked_hashes_match"), all(verified$hash_match), nrow(verified))
}

audit_markers <- c(
  `29f` = file.path(DIR_PROV, "29f_formal_independent_audit_20260722_142013_pid31436_v2", "audit_completed.ok"),
  `29g` = file.path(DIR_PROV, "29g_post_finalization_audit_20260722_194720_pid5904", "audit_completed.ok"),
  `35` = file.path(DIR_PROV, "35_WP9_formal_independent_audit_20260722_2255_r453", "audit_completed.ok")
)
audit_required <- c(
  `29f` = "downstream_authorization=GO",
  `29g` = "formal_29g_finalization=GO",
  `35` = "downstream_authorization=GO"
)
for (nm in names(audit_markers)) {
  add_qc(paste0(nm, "_audit_marker_exists"), file.exists(audit_markers[[nm]]), audit_markers[[nm]])
  marker_text <- paste(readLines(audit_markers[[nm]], warn = FALSE), collapse = "\n")
  add_qc(
    paste0(nm, "_independent_audit_go"),
    grepl(audit_required[[nm]], marker_text, fixed = TRUE),
    audit_required[[nm]]
  )
}

neg <- read_csv(file.path(DIR_TABLE, "Table_D1_copula_negative_control_LLOYD.csv"))
pos <- read_csv(file.path(DIR_TABLE, "Table_D1_copula_positive_control_LLOYD.csv"))
gate_n <- read_csv(file.path(DIR_TABLE, "Table_D1_formal_gate_by_n_LLOYD.csv"))
inv <- read_csv(file.path(DIR_TABLE, "Table_D1_empirical_delta_inversion_LLOYD.csv"))
gate_split <- read_csv(file.path(DIR_TABLE, "35_WP9_gate_by_split_LLOYD.csv"))
state_split <- read_csv(file.path(DIR_TABLE, "35_WP9_state_distribution_LLOYD.csv"))
emp_split <- read_csv(file.path(DIR_TABLE, "35_WP9_empirical_summary_LLOYD.csv"))
control_split <- read_csv(file.path(DIR_TABLE, "35_WP9_control_summary_LLOYD.csv"))
q95_contrast <- read_csv(file.path(DIR_TABLE, "35_WP9_q95_contrasts_LLOYD.csv"))
transitions <- read_csv(file.path(DIR_TABLE, "35_WP9_state_transitions_LLOYD.csv"))
failure_35 <- read_csv(file.path(DIR_LOG, "35_WP9_failure_log_LLOYD.csv"))

add_qc("29f_source_delta_fpr_locked", isTRUE(all.equal(neg$source_delta_false_positive_rate, 0.495)), neg$source_delta_false_positive_rate)
add_qc("29f_source_joint_fpr_zero", neg$source_joint_alert_false_positive_rate == 0, neg$source_joint_alert_false_positive_rate)
add_qc("29f_audit_null_delta_fpr", isTRUE(all.equal(neg$audit_null_delta_false_positive_rate, 0.065)), neg$audit_null_delta_false_positive_rate)
add_qc("29f_audit_null_joint_fpr_zero", neg$audit_null_joint_alert_false_positive_rate == 0, neg$audit_null_joint_alert_false_positive_rate)
add_qc("29f_s7_detection_complete", pos$source_joint_alert_rate_control_only == 1, pos$source_joint_alert_rate_control_only)
add_qc("29f_s7_true_label_ari_one", pos$source_true_label_ari_mean == 1, pos$source_true_label_ari_mean)

add_qc("29g_gate_has_n19049", sum(gate_n$n_training_generated == 19049L) == 1L, paste(gate_n$n_training_generated, collapse = ","))
add_qc("29g_empirical_single_row", nrow(inv) == 1L, nrow(inv))
add_qc("29g_empirical_boundary_not_verdict", grepl("not an empirical subtype verdict", inv$reporting_boundary, fixed = TRUE), inv$reporting_boundary)

add_qc("35_three_splits", nrow(gate_split) == 3L && uniqueN(gate_split$split_id) == 3L, nrow(gate_split))
add_qc("35_gate_rep_200_each", all(gate_split$n_gate == 200L), paste(gate_split$n_gate, collapse = ","))
add_qc("35_state_rows_complete", nrow(state_split) == 27L, nrow(state_split))
add_qc("35_state_denominator_100", all(state_split$n == 100L), paste(unique(state_split$n), collapse = ","))
add_qc("35_failure_log_empty", nrow(failure_35) == 0L, nrow(failure_35))

state_wide <- dcast(state_split, stage + split_order + split_id + split_label + n_train + n_eval ~ state, value.var = "count")
for (col in c("DISCRETE_EVIDENCE", "INCONCLUSIVE", "NO_DISCRETE_EVIDENCE")) {
  if (!col %in% names(state_wide)) state_wide[, (col) := 0L]
}
add_qc(
  "35_empirical_all_inconclusive",
  all(state_wide[stage == "locked_empirical"]$INCONCLUSIVE == 100L),
  paste(state_wide[stage == "locked_empirical"]$INCONCLUSIVE, collapse = ",")
)
add_qc(
  "35_s7_all_discrete",
  all(state_wide[stage == "S7_positive_control"]$DISCRETE_EVIDENCE == 100L),
  paste(state_wide[stage == "S7_positive_control"]$DISCRETE_EVIDENCE, collapse = ",")
)
add_qc(
  "35_n1_never_discrete",
  all(state_wide[stage == "N1_independent_control"]$DISCRETE_EVIDENCE == 0L),
  paste(state_wide[stage == "N1_independent_control"]$DISCRETE_EVIDENCE, collapse = ",")
)

gate_19049 <- gate_n[n_training_generated == 19049L]
split_19049 <- gate_split[n_train == 19049L]
add_qc("29g_35_n19049_gate_ci_overlap", max(gate_19049$q95_bootstrap_ci_low, split_19049$q95_bootstrap_ci_low) <= min(gate_19049$q95_bootstrap_ci_high, split_19049$q95_bootstrap_ci_high), paste(gate_19049$q95, split_19049$q95))

cross_stage <- rbindlist(list(
  data.table(
    module = "29f", evidence = "S1 source Delta-only false-positive rate",
    estimate = neg$source_delta_false_positive_rate,
    ci_low = neg$source_delta_false_positive_wilson_low,
    ci_high = neg$source_delta_false_positive_wilson_high,
    numerator = neg$source_delta_false_positive_n, denominator = neg$n_evaluation_rep,
    interpretation = "Single-statistic failure remains visible; it is not the joint D1 verdict.",
    source_file = "Table_D1_copula_negative_control_LLOYD.csv",
    source_field = "source_delta_false_positive_rate"
  ),
  data.table(
    module = "29f", evidence = "S1 source joint-alert false-positive rate",
    estimate = neg$source_joint_alert_false_positive_rate,
    ci_low = neg$source_joint_alert_wilson_low,
    ci_high = neg$source_joint_alert_wilson_high,
    numerator = neg$source_joint_alert_false_positive_n, denominator = neg$n_evaluation_rep,
    interpretation = "The joint separation-and-shape rule did not flag the source S1 control.",
    source_file = "Table_D1_copula_negative_control_LLOYD.csv",
    source_field = "source_joint_alert_false_positive_rate"
  ),
  data.table(
    module = "29f", evidence = "S7 joint-alert detection rate",
    estimate = pos$source_joint_alert_rate_control_only,
    ci_low = pos$source_joint_alert_wilson_low,
    ci_high = pos$source_joint_alert_wilson_high,
    numerator = pos$source_joint_alert_n_control_only, denominator = pos$n_evaluation_rep,
    interpretation = "Strong injected discrete structure was detected under the joint rule.",
    source_file = "Table_D1_copula_positive_control_LLOYD.csv",
    source_field = "source_joint_alert_rate_control_only"
  ),
  data.table(
    module = "29g", evidence = "Formal N1 gate at generated n=19,049",
    estimate = gate_19049$q95,
    ci_low = gate_19049$q95_bootstrap_ci_low,
    ci_high = gate_19049$q95_bootstrap_ci_high,
    numerator = NA_integer_, denominator = gate_19049$n_completed,
    interpretation = "Full-size N1 reference threshold; not an empirical subtype verdict.",
    source_file = "Table_D1_formal_gate_by_n_LLOYD.csv",
    source_field = "q95"
  ),
  data.table(
    module = "29g", evidence = "Locked full-data pooled-within Mahalanobis separation",
    estimate = inv$empirical_pooled_within_mahalanobis,
    ci_low = NA_real_, ci_high = NA_real_, numerator = NA_integer_, denominator = NA_integer_,
    interpretation = "Descriptive full-data separation used for DGM calibration only.",
    source_file = "Table_D1_empirical_delta_inversion_LLOYD.csv",
    source_field = "empirical_pooled_within_mahalanobis"
  ),
  data.table(
    module = "29g", evidence = "DGM-scale inverted injected Delta",
    estimate = inv$delta_inversion,
    ci_low = inv$bootstrap_ci_low,
    ci_high = inv$bootstrap_ci_high,
    numerator = NA_integer_, denominator = inv$bootstrap_valid_n,
    interpretation = "Calibration on the simulation DGM scale; not a biological effect estimate.",
    source_file = "Table_D1_empirical_delta_inversion_LLOYD.csv",
    source_field = "delta_inversion"
  )
), fill = TRUE)

gate_source <- copy(gate_split)
gate_source[, `:=`(
  module = "35",
  x_label = sprintf("Train %s\nEvaluate %s", format(n_train, big.mark = ","), format(n_eval, big.mark = ",")),
  cross_stage_n19049_q95 = ifelse(n_train == 19049L, gate_19049$q95, NA_real_),
  cross_stage_n19049_q95_ci_low = ifelse(n_train == 19049L, gate_19049$q95_bootstrap_ci_low, NA_real_),
  cross_stage_n19049_q95_ci_high = ifelse(n_train == 19049L, gate_19049$q95_bootstrap_ci_high, NA_real_)
)]

stage_labels <- c(
  N1_independent_control = "N1 continuous null",
  S7_positive_control = "S7 discrete control",
  locked_empirical = "Locked empirical partition"
)
state_labels <- c(
  DISCRETE_EVIDENCE = "Discrete evidence",
  INCONCLUSIVE = "Inconclusive",
  NO_DISCRETE_EVIDENCE = "No discrete evidence"
)
state_source <- copy(state_split)
state_source[, stage_label := unname(stage_labels[stage])]
state_source[, state_label := unname(state_labels[state])]
state_source[, split_short := sprintf("%s / %s", format(n_train, big.mark = ","), format(n_eval, big.mark = ","))]

emp_source <- merge(
  emp_split,
  gate_split[, .(split_id, q95, q95_bootstrap_ci_low, q95_bootstrap_ci_high)],
  by = "split_id",
  all.x = TRUE,
  sort = FALSE
)
emp_source[, split_short := sprintf("%s / %s", format(n_train, big.mark = ","), format(n_eval, big.mark = ","))]
emp_source[, final_state := "Inconclusive"]

transition_source <- transitions[, .(
  count = sum(count),
  stable_count = sum(count[row_state == column_state]),
  n_matched = unique(n_matched)
), by = .(stage, earlier_split, later_split)]
transition_source[, stable_rate := stable_count / n_matched]
transition_source[, stage_label := unname(stage_labels[stage])]
transition_source[, transition_label := paste(
  sub("train", "", earlier_split),
  "to",
  sub("train", "", later_split)
)]

runtime_files <- c(
  `29f` = file.path(DIR_LOG, "29f_sessionInfo.txt"),
  `29g` = file.path(DIR_LOG, "29g_sessionInfo_LLOYD.txt"),
  `35` = file.path(DIR_PROV, "35_WP9_sessionInfo_LLOYD.txt")
)
runtimes <- rbindlist(lapply(names(runtime_files), function(nm) {
  first <- readLines(runtime_files[[nm]], warn = FALSE, n = 1L)
  data.table(module = nm, r_version_line = first, session_file = runtime_files[[nm]], session_sha256 = sha256_file(runtime_files[[nm]]))
}))
add_qc("runtime_versions_disclosed", uniqueN(runtimes$r_version_line) == 2L, paste(runtimes$r_version_line, collapse = " | "))

replacement <- data.table(
  item = c(
    "Marginal-unimodality qualification",
    "Gap statistic",
    "SigClust-only interpretation",
    "Empirical D1 verdict",
    "Single-statistic control performance",
    "Joint-rule control performance",
    "Split-allocation sensitivity",
    "DGM inversion",
    "Runtime reproducibility"
  ),
  deprecated_or_unsafe = c(
    "The initial marginal-unimodality precondition was unmet because 33/33 marginal tests rejected.",
    "Gap statistic was used as the decisive D1 evidence.",
    "SigClust significance alone established or excluded latent discreteness.",
    "The empirical K2 partition was D1-negative or non-discrete.",
    "The separation statistic had acceptable null behavior.",
    "Control performance was perfect without qualification.",
    "The D1 verdict was independent of train/evaluation allocation because the threshold was constant.",
    "The inverted Delta was a biological separation effect or subtype verdict.",
    "All D1 modules used the same R patch release."
  ),
  required_replacement = c(
    "The reference uses empirical margins by construction; marginal non-unimodality (33/33) is descriptive and is not a qualification criterion. The reference was qualified by S1 negative and S7 positive controls.",
    "Gap was retired because it addresses partition-number selection rather than the prespecified null-scale separation estimand used here.",
    "D1 combined pooled-within Mahalanobis separation against an empirical-margin Gaussian-copula N1 reference with Fisher-axis shape evidence.",
    "Across all three examined train/evaluation allocations, the locked empirical partition triggered the separation gate but not the shape gate and was classified as inconclusive in 100/100 repeats.",
    "The source S1 Delta-only false-positive rate was 99/200 (49.5%; Wilson 42.6%-56.4%), so separation alone was not accepted as a D1 verdict.",
    "The joint rule flagged 0/200 source S1 repeats and 200/200 S7 repeats; these controls qualify the reference but do not validate biological subtype status.",
    "The split-specific q95 increased from 2.7169 to 2.7353 as the evaluation set increased; the empirical three-state conclusion remained inconclusive in every examined allocation.",
    "The DGM-scale inversion was 2.517 (bootstrap interval 1.310-2.517) and is reported only as simulation calibration, not as a biological effect estimate.",
    "Modules 29f and 29g ran under R 4.5.1; module 35 and the post-processing/figure scripts ran under R 4.5.3. Module-specific session files are archived."
  ),
  manuscript_location = c(
    "Methods and D1 reference paragraph",
    "Methods and deprecation map",
    "Methods and D1 Results",
    "Abstract/Results/Discussion",
    "D1 Results and Supplement",
    "D1 Results and Supplement",
    "D1 sensitivity Results",
    "Supplement only",
    "Reproducibility statement"
  )
)

qc <- rbindlist(qc_rows)
add_qc("all_qc_before_write", all(qc$pass), nrow(qc))
qc <- rbindlist(qc_rows)

tmp_path <- function(final_path) file.path(DIR_TMP, basename(final_path))
fwrite(cross_stage, tmp_path(OUT[["cross_stage"]]))
fwrite(gate_source, tmp_path(OUT[["gate"]]))
fwrite(state_source, tmp_path(OUT[["states"]]))
fwrite(emp_source, tmp_path(OUT[["empirical"]]))
fwrite(transition_source, tmp_path(OUT[["transitions"]]))
fwrite(runtimes, tmp_path(OUT[["runtimes"]]))
fwrite(qc, tmp_path(OUT[["qc"]]))
fwrite(replacement, tmp_path(OUT[["replacement_csv"]]))

input_files <- unique(c(
  completion_paths, manifest_paths, audit_markers, runtime_files,
  file.path(DIR_TABLE, c(
    "Table_D1_copula_negative_control_LLOYD.csv",
    "Table_D1_copula_positive_control_LLOYD.csv",
    "Table_D1_formal_gate_by_n_LLOYD.csv",
    "Table_D1_empirical_delta_inversion_LLOYD.csv",
    "35_WP9_gate_by_split_LLOYD.csv",
    "35_WP9_state_distribution_LLOYD.csv",
    "35_WP9_empirical_summary_LLOYD.csv",
    "35_WP9_control_summary_LLOYD.csv",
    "35_WP9_q95_contrasts_LLOYD.csv",
    "35_WP9_state_transitions_LLOYD.csv"
  )),
  file.path(DIR_LOG, "35_WP9_failure_log_LLOYD.csv")
))
inputs <- data.table(
  path = input_files,
  bytes = as.numeric(file.info(input_files)$size),
  sha256 = vapply(input_files, sha256_file, character(1))
)
fwrite(inputs, tmp_path(OUT[["inputs"]]))

md <- c(
  "# D1 Lloyd-aligned final replacement map",
  "",
  "This file is a post-processing artifact. It does not replace or modify the locked 29f, 29g, or 35 outputs.",
  "",
  "## Locked interpretation",
  "",
  "- The 29f source S1 separation-only false-positive rate was 99/200 (49.5%); this failure remains visible.",
  "- The prespecified joint rule flagged 0/200 source S1 controls and 200/200 S7 controls.",
  "- Split-specific N1 q95 values were 2.7169, 2.7267, and 2.7353.",
  "- The locked empirical partition was inconclusive in 100/100 repeats under each of the three split allocations.",
  "- The full-data pooled-within Mahalanobis separation was 3.5211. Its DGM-scale inversion (2.517; bootstrap interval 1.310-2.517) is calibration only.",
  "",
  "## Replacement list",
  ""
)
for (i in seq_len(nrow(replacement))) {
  md <- c(
    md,
    paste0("### ", i, ". ", replacement$item[i]),
    "",
    paste0("**Deprecated/unsafe:** ", replacement$deprecated_or_unsafe[i]),
    "",
    paste0("**Use instead:** ", replacement$required_replacement[i]),
    "",
    paste0("**Location:** ", replacement$manuscript_location[i]),
    ""
  )
}
writeLines(md, tmp_path(OUT[["replacement_md"]]), useBytes = TRUE)

move_names <- setdiff(names(OUT), c("manifest", "completion"))
for (nm in move_names) {
  src <- tmp_path(OUT[[nm]])
  if (!file.exists(src)) stop("Temporary output missing: ", src)
  if (!file.rename(src, OUT[[nm]])) stop("Could not finalize output: ", OUT[[nm]])
}

manifest_targets <- unname(OUT[move_names])
out_manifest <- data.table(
  artifact = move_names,
  path = manifest_targets,
  bytes = as.numeric(file.info(manifest_targets)$size),
  sha256 = vapply(manifest_targets, sha256_file, character(1))
)
fwrite(out_manifest, OUT[["manifest"]])

completion <- c(
  "POST-PROCESSING COMPLETED",
  "run_mode=read_only_postprocess",
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("r_version=", R.version.string),
  paste0("script_sha256=", sha256_file(SCRIPT_PATH)),
  paste0("input_manifest_sha256=", sha256_file(OUT[["inputs"]])),
  paste0("output_manifest_sha256=", sha256_file(OUT[["manifest"]])),
  paste0("qc_passed=", sum(qc$pass), "/", nrow(qc)),
  "locked_source_outputs_modified=NO",
  "scientific_models_refit=NO"
)
writeLines(completion, OUT[["completion"]], useBytes = TRUE)

cat("36 cross-stage comparison audit completed.\n")
cat("QC:", sum(qc$pass), "/", nrow(qc), "passed.\n")
cat("Runtimes disclosed:", paste(unique(runtimes$r_version_line), collapse = " | "), "\n")
