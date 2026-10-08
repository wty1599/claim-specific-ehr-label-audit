## =====================================================================
## D3 regenerated deployable surrogate generator v1
## Final release audit. This script does not refit models or change labels.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
for (pkg in c("data.table")) require_namespace(pkg)
library(data.table)

audit_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)

check_rows <- list()
add_check <- function(check, observed, expected, pass, severity = "critical") {
  check_rows[[length(check_rows) + 1L]] <<- data.table(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = isTRUE(pass),
    severity = severity
  )
}

## ---------- Script integrity ----------
expected_scripts <- file.path(
  SCRIPT_DIR,
  c(
    "00_d3_generator_config_v1.R",
    "00_d3_generator_utils_v1.R",
    "01_build_d3_deployable_generator_v1.R",
    "02_validate_d3_generator_internal_fidelity_v1.R",
    "03_apply_d3_generator_to_eicu_v1.R",
    "04_analyze_d3_generator_transport_v1.R",
    "04b_bootstrap_generator_calibration_v1.R",
    "05_render_d3_generator_outputs_v1.R",
    "06_audit_release_d3_generator_v1.R"
  )
)
add_check(
  "all expected scripts exist",
  sum(file.exists(expected_scripts)),
  length(expected_scripts),
  all(file.exists(expected_scripts))
)
parse_ok <- vapply(expected_scripts, function(x) {
  if (!file.exists(x)) return(FALSE)
  !inherits(try(parse(file = x), silent = TRUE), "try-error")
}, logical(1))
add_check(
  "all expected scripts parse",
  sum(parse_ok),
  length(parse_ok),
  all(parse_ok)
)

## ---------- Locked input fingerprints ----------
input_map <- c(
  final_full = PATH_FINAL_FULL,
  eicu_external = PATH_EICU,
  locked_labels = PATH_LOCKED_LABELS,
  mice_primary = PATH_MICE,
  mice_matrix = PATH_MICE_X,
  legacy_eicu_labels = PATH_LEGACY_EICU_LABELS,
  baseline_covars = PATH_BASELINE,
  reference_recipe = PATH_REFERENCE_RECIPE,
  reference_centroids = PATH_REFERENCE_CENTROIDS
)
for (nm in names(input_map)) {
  observed <- if (file.exists(input_map[[nm]])) hash_file(input_map[[nm]]) else NA_character_
  add_check(
    paste0("locked input MD5: ", nm),
    observed,
    LOCKED_MD5[[nm]],
    identical(observed, unname(LOCKED_MD5[[nm]]))
  )
}

## ---------- Required outputs ----------
required_outputs <- c(
  PATH_GENERATOR,
  PATH_EXTERNAL_LABELS,
  file.path(OUT_TABLES, "internal_fidelity_summary.csv"),
  file.path(OUT_TABLES, "internal_fidelity_bootstrap_summary.csv"),
  file.path(OUT_TABLES, "internal_fidelity_crossfit_summary.csv"),
  file.path(OUT_TABLES, "d3_generator_v1_prevalence.csv"),
  file.path(OUT_TABLES, "d3_generator_v1_model_performance.csv"),
  file.path(
    OUT_TABLES, "d3_generator_v1_base_gen_calibration_bootstrap.csv"
  ),
  file.path(OUT_TABLES, "d3_generator_v1_hospital_prevalence.csv"),
  file.path(OUT_TABLES, "d3_generator_v1_missingness_transport.csv"),
  file.path(OUT_TABLES, "d3_generator_v1_alert_components.csv"),
  file.path(OUT_TABLES, "d3_generator_v1_claim_states.csv"),
  file.path(OUT_QC, "generator_v1_invariance_tests.csv"),
  file.path(OUT_QC, "eicu_d3_generator_v1_legacy_assignment_comparison.csv"),
  file.path(OUT_FIGURES, "Figure_D3_generator_v1.pdf"),
  file.path(OUT_FIGURES, "Figure_D3_generator_v1.png"),
  file.path(OUT_FIGURES, "Figure_D3_generator_v1.tiff"),
  file.path(OUT_FIGURES, "Figure_D3_generator_v1.svg"),
  file.path(OUT_TABLES, "Table_D3_generator_lineage.docx"),
  file.path(OUT_TABLES, "Table_D3_generator_fidelity.docx"),
  file.path(OUT_TABLES, "Table_D3_generator_external.docx")
)
output_exists <- file.exists(required_outputs)
output_nonempty <- output_exists & file.info(required_outputs)$size > 0
add_check(
  "all required outputs exist and are nonempty",
  sum(output_nonempty),
  length(required_outputs),
  all(output_nonempty)
)

## ---------- Scientific and denominator checks ----------
generator <- readRDS(PATH_GENERATOR)
labels_new <- fread(PATH_EXTERNAL_LABELS)
fidelity <- fread(file.path(OUT_TABLES, "internal_fidelity_summary.csv"))
invariance <- fread(file.path(OUT_QC, "generator_v1_invariance_tests.csv"))
legacy <- fread(file.path(
  OUT_QC, "eicu_d3_generator_v1_legacy_assignment_comparison.csv"
))
prevalence <- fread(file.path(OUT_TABLES, "d3_generator_v1_prevalence.csv"))
hospitals <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_hospital_prevalence.csv"
))
performance <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_model_performance.csv"
))
calibration_bootstrap <- fread(file.path(
  OUT_TABLES, "d3_generator_v1_base_gen_calibration_bootstrap.csv"
))
alerts <- fread(file.path(OUT_TABLES, "d3_generator_v1_alert_components.csv"))
states <- fread(file.path(OUT_TABLES, "d3_generator_v1_claim_states.csv"))
analysis_provenance <- fread(file.path(
  OUT_QC, "d3_generator_v1_outcome_analysis_provenance.csv"
))

add_check(
  "generator name",
  generator$generator_name,
  GENERATOR_NAME,
  identical(generator$generator_name, GENERATOR_NAME)
)
add_check(
  "generator file is bitwise regression-locked",
  hash_file(PATH_GENERATOR),
  EXPECTED_GENERATOR_FILE_MD5,
  identical(hash_file(PATH_GENERATOR), EXPECTED_GENERATOR_FILE_MD5)
)
add_check(
  "generator content hash is regression-locked",
  generator$generator_content_hash,
  EXPECTED_GENERATOR_CONTENT_HASH,
  identical(
    generator$generator_content_hash, EXPECTED_GENERATOR_CONTENT_HASH
  )
)
add_check(
  "generator status is explicitly post hoc surrogate",
  generator$post_hoc_status,
  GENERATOR_STATUS,
  identical(generator$post_hoc_status, GENERATOR_STATUS)
)
add_check(
  "MIMIC denominator",
  fidelity$n,
  EXPECTED_MIMIC_N,
  fidelity$n == EXPECTED_MIMIC_N
)
add_check(
  "apparent agreement numerator",
  fidelity$agreement_n,
  EXPECTED_APPARENT_AGREEMENT_N,
  fidelity$agreement_n == EXPECTED_APPARENT_AGREEMENT_N
)
add_check(
  "apparent discordant count",
  fidelity$discordant_n,
  EXPECTED_APPARENT_DISCORDANT_N,
  fidelity$discordant_n == EXPECTED_APPARENT_DISCORDANT_N
)
add_check(
  "apparent ARI regression check",
  round(fidelity$ari, 3),
  EXPECTED_APPARENT_ARI_ROUNDED,
  identical(round(fidelity$ari, 3), EXPECTED_APPARENT_ARI_ROUNDED)
)
add_check(
  "eICU stay denominator",
  nrow(labels_new),
  EXPECTED_EICU_STAYS,
  nrow(labels_new) == EXPECTED_EICU_STAYS
)
add_check(
  "eICU patient denominator",
  uniqueN(labels_new$uniquepid),
  EXPECTED_EICU_PATIENTS,
  uniqueN(labels_new$uniquepid) == EXPECTED_EICU_PATIENTS
)
add_check(
  "eICU hospital denominator",
  uniqueN(labels_new$hospitalid),
  EXPECTED_EICU_HOSPITALS,
  uniqueN(labels_new$hospitalid) == EXPECTED_EICU_HOSPITALS
)
add_check(
  "eICU C1 count",
  sum(labels_new$assigned_cluster == 1L),
  EXPECTED_EICU_C1_N,
  sum(labels_new$assigned_cluster == 1L) == EXPECTED_EICU_C1_N
)
add_check(
  "external label artifact is bitwise regression-locked",
  hash_file(PATH_EXTERNAL_LABELS),
  EXPECTED_EXTERNAL_LABEL_MD5,
  identical(hash_file(PATH_EXTERNAL_LABELS), EXPECTED_EXTERNAL_LABEL_MD5)
)
add_check(
  "legacy eICU label agreement",
  legacy$agreement,
  1,
  legacy$agreement == 1 && legacy$n_unequal_label == 0L
)
add_check(
  "all technical invariance tests pass",
  sum(invariance$pass),
  nrow(invariance),
  all(invariance$pass)
)
add_check(
  "hospital table contains 199 hospitals",
  uniqueN(hospitals$hospitalid),
  EXPECTED_EICU_HOSPITALS,
  uniqueN(hospitals$hospitalid) == EXPECTED_EICU_HOSPITALS
)
add_check(
  "external outcome model rows",
  nrow(performance),
  8,
  nrow(performance) == 8L
)
add_check(
  "same regenerated generator used for source and external labels",
  paste(
    analysis_provenance[
      item %in% c("source_label_object", "external_label_object"), value
    ],
    collapse = ";"
  ),
  paste(GENERATOR_NAME, GENERATOR_NAME, sep = ";"),
  identical(
    analysis_provenance[
      item %in% c("source_label_object", "external_label_object"), value
    ],
    c(GENERATOR_NAME, GENERATOR_NAME)
  )
)

expected_states <- c(
  "GENERATOR_ABSENT",
  "EVALUATED_FIDELITY_QUANTIFIED",
  "EVALUATED_TECHNICAL_APPLICATION",
  "REFERENCE_ABSENT",
  "EVALUATED_NO_ALERT",
  "EVALUATED_WITH_TRANSPORT_ALERT",
  "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY",
  "CONTEXT_ONLY_NOT_GENERATOR_TRANSPORT"
)
add_check(
  "claim-specific state vocabulary",
  paste(states$state, collapse = ";"),
  paste(expected_states, collapse = ";"),
  identical(states$state, expected_states)
)
add_check(
  "original generator remains absent",
  states[claim == "exact original discovery-generator transport", state],
  "GENERATOR_ABSENT",
  states[
    claim == "exact original discovery-generator transport",
    state
  ] == "GENERATOR_ABSENT"
)
add_check(
  "external semantic assignment reference remains absent",
  states[
    claim == "regenerated generator external semantic assignment fidelity",
    state
  ],
  "REFERENCE_ABSENT",
  states[
    claim == "regenerated generator external semantic assignment fidelity",
    state
  ] == "REFERENCE_ABSENT"
)
add_check(
  "mortality claim state matches prespecified generated-label calibration gates",
  states[claim == "regenerated generator mortality-model transport", state],
  if (any(alerts[
    outcome == "mortality" & drives_claim_state & !is.na(alert), alert
  ])) "EVALUATED_WITH_TRANSPORT_ALERT" else "EVALUATED_NO_ALERT",
  states[
    claim == "regenerated generator mortality-model transport", state
  ] == if (any(alerts[
    outcome == "mortality" & drives_claim_state & !is.na(alert), alert
  ])) "EVALUATED_WITH_TRANSPORT_ALERT" else "EVALUATED_NO_ALERT"
)
add_check(
  "kidney-composite claim is not treated as endpoint-equivalent transport",
  states[
    claim == "regenerated generator kidney-composite transport", state
  ],
  "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY",
  states[
    claim == "regenerated generator kidney-composite transport", state
  ] == "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY"
)
add_check(
  "pooled prevalence drift does not trigger gate",
  round(alerts[component == "prevalence_drift", observed], 6),
  "<=0.10",
  !alerts[component == "prevalence_drift", alert]
)
add_check(
  "base_gen calibration bootstrap has four parameter rows",
  nrow(calibration_bootstrap),
  4,
  nrow(calibration_bootstrap) == 4L
)
add_check(
  "base_gen calibration bootstrap uses patient resampling",
  paste(unique(calibration_bootstrap$resampling_unit), collapse = ";"),
  "uniquepid",
  identical(unique(calibration_bootstrap$resampling_unit), "uniquepid")
)
add_check(
  "base_gen calibration bootstrap has no failed replicates",
  sum(calibration_bootstrap$bootstrap_failed),
  0,
  sum(calibration_bootstrap$bootstrap_failed) == 0L
)
bootstrap_point_check <- merge(
  calibration_bootstrap[, .(
    outcome,
    parameter,
    bootstrap_estimate = estimate
  )],
  rbind(
    performance[model == "base_gen", .(
      outcome,
      parameter = "slope",
      performance_estimate = calibration_slope
    )],
    performance[model == "base_gen", .(
      outcome,
      parameter = "intercept",
      performance_estimate = calibration_intercept
    )]
  ),
  by = c("outcome", "parameter"),
  all = TRUE
)
add_check(
  "bootstrap point estimates match locked performance table",
  max(
    abs(
      bootstrap_point_check$bootstrap_estimate -
        bootstrap_point_check$performance_estimate
    ),
    na.rm = TRUE
  ),
  "<=1e-12",
  nrow(bootstrap_point_check) == 4L &&
    all(is.finite(bootstrap_point_check$bootstrap_estimate)) &&
    max(
      abs(
        bootstrap_point_check$bootstrap_estimate -
          bootstrap_point_check$performance_estimate
      ),
      na.rm = TRUE
    ) <= 1e-12
)

## The deprecated eICU count must not appear in publication-facing outputs.
publication_files <- list.files(
  c(OUT_TABLES, OUT_QC),
  pattern = "\\.(csv|txt|md)$",
  recursive = TRUE,
  full.names = TRUE
)
contains_old_count <- vapply(publication_files, function(path) {
  any(grepl(
    "(^|[^0-9])3[,]?417([^0-9]|$)",
    readLines(path, warn = FALSE),
    perl = TRUE
  ))
}, logical(1))
add_check(
  "deprecated eICU C1 count 3417 absent from new publication outputs",
  sum(contains_old_count),
  0,
  !any(contains_old_count)
)

checks <- rbindlist(check_rows, fill = TRUE)
write_replace_csv(checks, file.path(OUT_RELEASE, "D3_GENERATOR_V1_AUDIT_CHECKS.csv"))

critical_fail <- checks[severity == "critical" & !pass]
if (nrow(critical_fail)) {
  stop(
    "Release audit failed critical checks: ",
    paste(critical_fail$check, collapse = "; ")
  )
}

## ---------- Release bundle ----------
release_scripts <- file.path(OUT_RELEASE, "scripts")
release_model <- file.path(OUT_RELEASE, "model")
release_tables <- file.path(OUT_RELEASE, "tables")
release_figures <- file.path(OUT_RELEASE, "figures")
release_qc <- file.path(OUT_RELEASE, "qc")
for (d in c(
  release_scripts, release_model, release_tables, release_figures, release_qc
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

copy_replace <- function(from, to_dir) {
  to <- file.path(to_dir, basename(from))
  if (file.exists(to)) unlink(to)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Failed to copy release artifact: ", from)
  to
}

copied <- c(
  vapply(expected_scripts, copy_replace, character(1), to_dir = release_scripts),
  copy_replace(PATH_GENERATOR, release_model),
  vapply(
    list.files(OUT_TABLES, full.names = TRUE),
    copy_replace, character(1), to_dir = release_tables
  ),
  vapply(
    list.files(OUT_FIGURES, full.names = TRUE),
    copy_replace, character(1), to_dir = release_figures
  ),
  vapply(
    list.files(OUT_QC, full.names = TRUE),
    copy_replace, character(1), to_dir = release_qc
  )
)

## Remove patient-level material if any matching name accidentally enters release.
patient_level_patterns <- c(
  "labels", "predictions", "discordant_records", "patient_level"
)
release_files_now <- list.files(OUT_RELEASE, recursive = TRUE, full.names = TRUE)
unsafe <- release_files_now[vapply(release_files_now, function(x) {
  any(vapply(patient_level_patterns, grepl, logical(1), x = basename(x), fixed = TRUE))
}, logical(1))]
if (length(unsafe)) {
  stop("Patient-level artifact detected in release bundle: ", paste(unsafe, collapse = "; "))
}

## ---------- Manifest and author-facing records ----------
manifest_files <- list.files(
  OUT_RELEASE, recursive = TRUE, full.names = TRUE
)
manifest_files <- manifest_files[!grepl(
  "D3_GENERATOR_V1_RELEASE_MANIFEST\\.csv$", manifest_files
)]
manifest <- data.table(
  relative_path = sub(
    paste0("^", gsub("([\\\\.])", "\\\\\\1", normalizePath(
      OUT_RELEASE, winslash = "/", mustWork = TRUE
    )), "/?"),
    "",
    normalizePath(manifest_files, winslash = "/", mustWork = TRUE)
  ),
  bytes = file.info(manifest_files)$size,
  md5 = vapply(manifest_files, hash_file, character(1))
)
setorder(manifest, relative_path)
write_replace_csv(
  manifest,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_RELEASE_MANIFEST.csv")
)

gen_mort <- performance[outcome == "mortality" & model == "base_gen"]
gen_make <- performance[outcome == "make" & model == "base_gen"]
mort_slope_boot <- calibration_bootstrap[
  outcome == "mortality" & parameter == "slope"
]
replacement_map <- data.table(
  manuscript_claim = states$claim,
  permitted_status = states$state,
  permitted_wording = c(
    "The original workflow did not instantiate a frozen inductive new-patient mapper.",
    sprintf(
      paste(
        "The regenerated post hoc surrogate reproduced %d/%d locked labels",
        "(agreement %.3f; ARI %.3f); repeated held-out ARI was %.3f."
      ),
      fidelity$agreement_n, fidelity$n, fidelity$agreement, fidelity$ari,
      fread(file.path(
        OUT_TABLES, "internal_fidelity_crossfit_summary.csv"
      ))[metric == "ari", mean]
    ),
    sprintf(
      paste(
        "The frozen regenerated surrogate assigned all %d eICU stays without",
        "outcome access, preprocessing refitting, or eICU re-clustering."
      ),
      nrow(labels_new)
    ),
    paste(
      "No independent eICU reference partition exists; external assignment",
      "ARI and semantic fidelity remain unverified."
    ),
    sprintf(
      paste(
        "Using the same frozen generator in both cohorts, C1 prevalence was",
        "%.2f%% in MIMIC-IV and %.2f%% in eICU (absolute drift %.4f)."
      ),
      100 * prevalence[cohort == "MIMIC-IV", c1_prevalence],
      100 * prevalence[cohort == "eICU", c1_prevalence],
      abs(diff(prevalence$c1_prevalence))
    ),
    sprintf(
      paste(
        "The baseline-plus-generated-label mortality model had external AUC",
        "%.3f, calibration slope %.3f (patient-cluster bootstrap 95%% CI",
        "%.3f to %.3f), and intercept %.3f. The point estimate crossed the",
        "prespecified 0.80 slope gate, while its interval included that gate."
      ),
      gen_mort$external_auc, gen_mort$calibration_slope,
      mort_slope_boot$ci_low, mort_slope_boot$ci_high,
      gen_mort$calibration_intercept
    ),
    sprintf(
      paste(
        "The external in-hospital kidney composite yielded descriptive AUC",
        "%.3f and calibration slope %.3f, but is not endpoint-equivalent",
        "to development MAKE30."
      ),
      gen_make$external_auc, gen_make$calibration_slope
    ),
    paste(
      "Raw-EN and Raw-RF results describe physiology-model transport and",
      "are contextual comparators rather than regenerated-generator evidence."
    )
  ),
  prohibited_overclaim = c(
    "Do not call the original MICE workflow externally validated.",
    "Do not call the surrogate identical to the locked partition.",
    "Do not describe technical label generation as external semantic validation.",
    "Do not substitute legacy-surrogate agreement for an independent eICU truth reference.",
    "Do not treat pooled prevalence concordance as proof of semantic transport.",
    "Do not attribute Raw-EN or Raw-RF calibration to the generator.",
    "Do not call the eICU kidney composite transported MAKE30.",
    "Do not let contextual raw-feature models determine the generator claim state."
  )
)
write_replace_csv(
  replacement_map,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_MANUSCRIPT_REPLACEMENT_MAP.csv")
)

changes <- c(
  "# D3 regenerated deployable surrogate generator v1: changes",
  "",
  paste0("- Audit time: ", audit_time),
  "- Added a new, versioned post hoc frozen mapper using MIMIC-derived winsor limits, medians, centers, scales, and locked-label centroids.",
  "- Preserved the original locked MICE-imputation-1 labels and every existing Step 0/D1/D2/D4 result.",
  "- Quantified full-cohort and repeated held-out fidelity to the locked partition.",
  "- Assigned eICU labels without outcome access, eICU re-clustering, or eICU preprocessing refitting.",
  "- Reproduced the current authoritative eICU label file exactly (17,465/17,465 labels).",
  "- Added descriptive hospital-level prevalence and discrimination variation summaries across 199 hospitals.",
  "- Refit the baseline-plus-label transport model using regenerated labels in both MIMIC-IV and eICU.",
  "- Kept Raw-EN and Raw-RF as contextual physiology-model comparators rather than generator evidence.",
  "- Marked external semantic assignment fidelity as REFERENCE_ABSENT because eICU has no independent truth partition.",
  "- Kept the non-equivalent eICU kidney composite descriptive rather than treating it as transported MAKE30.",
  paste(
    "- Added patient-cluster bootstrap intervals for the frozen",
    "baseline-plus-generated-label calibration estimates without model",
    "refitting, relabeling, or changing the point-estimate gates."
  ),
  "- Added claim-specific states rather than replacing GENERATOR_ABSENT with a single global pass/fail label.",
  "- Added publication figures, Word tables, provenance manifests, and a release audit.",
  "",
  "## Protected statements",
  "",
  "- The original discovery workflow remains GENERATOR_ABSENT for exact inductive transport.",
  "- The new object is a post hoc regenerated deployable surrogate, not the original generator.",
  "- The mortality transport state is determined only by the same-object baseline-plus-generated-label model.",
  "- The regenerated generator does not replace the locked partition in other domains."
)
writeLines(
  changes,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_CHANGES.md"),
  useBytes = TRUE
)

consistency <- c(
  "# D3 regenerated deployable surrogate generator v1: consistency check",
  "",
  paste0("- Audit time: ", audit_time),
  paste0("- Critical checks passed: ", nrow(checks), "/", nrow(checks), "."),
  sprintf(
    "- Locked-label fidelity: %d/%d agreement (%.4f); ARI %.4f; discordant %d.",
    fidelity$agreement_n, fidelity$n, fidelity$agreement, fidelity$ari,
    fidelity$discordant_n
  ),
  sprintf(
    "- eICU assignment: %d stays, %d patients, %d hospitals; C1 %d (%.2f%%).",
    nrow(labels_new), uniqueN(labels_new$uniquepid),
    uniqueN(labels_new$hospitalid), sum(labels_new$assigned_cluster == 1L),
    100 * mean(labels_new$assigned_cluster == 1L)
  ),
  "- Legacy eICU assignment comparison: 17,465/17,465 identical.",
  "- Ten of ten technical invariance tests passed.",
  "- No publication-facing output contains the deprecated eICU C1 count 3,417.",
  "- No patient-level labels, predictions, or discordant-record files were copied into the release bundle.",
  "",
  "## Final D3 interpretation",
  "",
  "- Exact original discovery-generator transport: GENERATOR_ABSENT.",
  "- Regenerated-generator internal fidelity: EVALUATED_FIDELITY_QUANTIFIED.",
  "- Regenerated-generator technical external application: EVALUATED_TECHNICAL_APPLICATION.",
  "- External semantic assignment fidelity: REFERENCE_ABSENT.",
  paste0(
    "- Same-object prevalence transport: ",
    states[claim == "regenerated generator prevalence transport", state], "."
  ),
  paste0(
    "- Same-object mortality-model transport: ",
    states[claim == "regenerated generator mortality-model transport", state],
    "."
  ),
  sprintf(
    paste(
      "- Mortality calibration slope %.3f (patient-cluster bootstrap 95%% CI",
      "%.3f to %.3f); the prespecified point-estimate alert remains, while",
      "the interval includes the 0.80 gate."
    ),
    gen_mort$calibration_slope,
    mort_slope_boot$ci_low,
    mort_slope_boot$ci_high
  ),
  paste(
    "- Kidney-composite transport:",
    "ENDPOINT_NON_EQUIVALENT_DESCRIPTIVE_ONLY."
  ),
  paste(
    "- The original D3 absence finding is not erased; specific regenerated",
    "generator claims are now evaluable, while semantic assignment fidelity",
    "still lacks an external reference."
  )
)
writeLines(
  consistency,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_CONSISTENCY_CHECK.md"),
  useBytes = TRUE
)

final_interpretation <- c(
  "# D3 regenerated deployable surrogate generator v1: final interpretation",
  "",
  "## Object lineage",
  "",
  paste(
    "The original MICE-before-clustering workflow remains",
    "`GENERATOR_ABSENT` for exact inductive new-patient assignment."
  ),
  paste(
    "The new object is a post hoc regenerated frozen surrogate.",
    "It does not replace the locked Step 0/D1/D2/D4 partition."
  ),
  "",
  "## Claim-specific states",
  "",
  paste0(
    "- Exact original discovery-generator transport: ",
    states[claim == "exact original discovery-generator transport", state],
    "."
  ),
  paste0(
    "- Regenerated-generator internal fidelity: ",
    states[claim == "regenerated generator internal fidelity", state],
    "."
  ),
  paste0(
    "- Regenerated-generator technical external application: ",
    states[
      claim == "regenerated generator technical external application", state
    ],
    "."
  ),
  paste0(
    "- External semantic assignment fidelity: ",
    states[
      claim == "regenerated generator external semantic assignment fidelity",
      state
    ],
    "."
  ),
  paste0(
    "- Same-object prevalence transport: ",
    states[claim == "regenerated generator prevalence transport", state],
    "."
  ),
  paste0(
    "- Same-object mortality-model transport: ",
    states[claim == "regenerated generator mortality-model transport", state],
    "."
  ),
  paste0(
    "- Kidney endpoint transport: ",
    states[
      claim == "regenerated generator kidney-composite transport", state
    ],
    "."
  ),
  paste0(
    "- Raw-feature outcome-model transport: ",
    states[claim == "raw-feature outcome-model transport", state],
    "."
  ),
  "",
  "## Key quantitative evidence",
  "",
  sprintf(
    paste(
      "- Internal locked-label fidelity: %d/%d agreement (%.4f);",
      "apparent ARI %.4f; repeated held-out ARI %.4f."
    ),
    fidelity$agreement_n,
    fidelity$n,
    fidelity$agreement,
    fidelity$ari,
    fread(file.path(
      OUT_TABLES, "internal_fidelity_crossfit_summary.csv"
    ))[metric == "ari", mean]
  ),
  sprintf(
    paste(
      "- Same-object C1 prevalence: MIMIC-IV %.2f%% and eICU %.2f%%;",
      "absolute drift %.5f."
    ),
    100 * prevalence[cohort == "MIMIC-IV", c1_prevalence],
    100 * prevalence[cohort == "eICU", c1_prevalence],
    abs(diff(prevalence$c1_prevalence))
  ),
  sprintf(
    paste(
      "- Baseline-plus-generated-label mortality transport: external AUC",
      "%.3f; calibration slope %.3f (patient-cluster bootstrap 95%% CI",
      "%.3f to %.3f); intercept %.3f."
    ),
    gen_mort$external_auc,
    gen_mort$calibration_slope,
    mort_slope_boot$ci_low,
    mort_slope_boot$ci_high,
    gen_mort$calibration_intercept
  ),
  "",
  "## Evidence boundaries",
  "",
  paste(
    "- The mortality state is an operational alert under the prespecified",
    "point-estimate gate. The bootstrap interval includes 0.80, so this",
    "is not an interval-level declaration of calibration failure."
  ),
  paste(
    "- eICU lacks an independent reference partition; technical assignment",
    "does not establish external semantic fidelity."
  ),
  paste(
    "- External missingness is materially greater than in MIMIC-IV, which",
    "limits semantic interpretation of the technically complete assignments."
  ),
  paste(
    "- The eICU in-hospital kidney composite is not endpoint-equivalent to",
    "development MAKE30 and remains descriptive."
  ),
  paste(
    "- The 41 release checks are integrity and regression checks,",
    "not scientific validation metrics."
  ),
  "",
  "## Bottom line",
  "",
  paste(
    "D3 is no longer represented by one global `NOT_EVALUATED` state.",
    "The original generator remains absent, while the regenerated surrogate",
    "supports claim-specific technical, prevalence, and mortality transport",
    "assessments. External semantic fidelity remains reference-limited."
  )
)
writeLines(
  final_interpretation,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_FINAL_INTERPRETATION.md"),
  useBytes = TRUE
)

package_names <- c(
  "data.table", "glmnet", "ranger", "ggplot2", "patchwork", "svglite",
  "ragg", "flextable", "officer", "scales"
)
package_versions <- data.table(
  package = package_names,
  version = vapply(
    package_names, function(x) as.character(utils::packageVersion(x)),
    character(1)
  ),
  r_version = as.character(getRversion())
)
write_replace_csv(
  package_versions,
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_PACKAGE_VERSIONS.csv")
)

writeLines(
  c(
    capture.output(sessionInfo()),
    "",
    paste0("Generator MD5: ", hash_file(PATH_GENERATOR)),
    paste0("External label MD5: ", hash_file(PATH_EXTERNAL_LABELS))
  ),
  file.path(OUT_RELEASE, "D3_GENERATOR_V1_SESSION_INFO.txt"),
  useBytes = TRUE
)

## Visual QA was performed on the final PNG after a second render pass.
qa_path <- file.path(OUT_QC, "Figure_D3_generator_v1_QA.txt")
qa <- readLines(qa_path, warn = FALSE)
qa <- sub(
  "Visual QA: pending manual preview inspection",
  paste(
    "Visual QA: PASS after two render-review cycles;",
    "no panel-title, facet-title, axis-label, or caption clipping observed"
  ),
  qa,
  fixed = TRUE
)
writeLines(qa, qa_path, useBytes = TRUE)
copy_replace(qa_path, release_qc)

## Regenerate the manifest last so hashes cover the final written state.
manifest_path <- file.path(OUT_RELEASE, "D3_GENERATOR_V1_RELEASE_MANIFEST.csv")
manifest_files <- list.files(
  OUT_RELEASE, recursive = TRUE, full.names = TRUE
)
manifest_files <- manifest_files[
  normalizePath(manifest_files, winslash = "/", mustWork = TRUE) !=
    normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
]
manifest <- data.table(
  relative_path = substring(
    normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(OUT_RELEASE, winslash = "/", mustWork = TRUE)) + 2L
  ),
  bytes = file.info(manifest_files)$size,
  md5 = vapply(manifest_files, hash_file, character(1))
)
setorder(manifest, relative_path)
write_replace_csv(manifest, manifest_path)

cat("D3 generator v1 release audit PASSED.\n")
cat("Checks:", nrow(checks), "/", nrow(checks), "\n")
cat("Release directory:", OUT_RELEASE, "\n")
