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
dir.create(out, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  registry = file.path(pkg, "config", "G1_G7_CANONICAL_REGISTRY.csv"),
  crosswalk = file.path(pkg, "config", "G3_G4_G6_SUBDESIGN_CROSSWALK.csv"),
  sa01_qc = file.path(pkg, "output_sa01_formal", "qc", "critical_qc.csv"),
  sa01_oc = file.path(pkg, "output_sa01_formal", "tables",
                   "D1_current_rule_operating_characteristics.csv"),
  sa01_cont = file.path(pkg, "output_sa01_formal", "tables",
                     "D1_current_rule_continuous_summary.csv"),
  sa05_qc = file.path(pkg, "output_sa05_formal", "qc", "critical_qc.csv"),
  sa05_oc = file.path(pkg, "output_sa05_formal", "tables",
                   "D1_current_rule_operating_characteristics.csv"),
  sa05_cont = file.path(pkg, "output_sa05_formal", "tables",
                     "D1_current_rule_continuous_summary.csv"),
  d2_g5 = file.path(
    root, "analysis_archive", "simulations", "locked_results_20260712",
    "02_SA02_locked", "corrected_domain2_v2", "tables",
    "domain2_delta_auc_setting_summary_corrected_12rows.csv"
  ),
  d3_oc = file.path(
    root, "analysis_archive", "simulations",
    "D3_regenerated_generator_transport_validation_20260802",
    "output_formal_v1", "tables",
    "D3_sim_operating_characteristics_by_scenario.csv"
  ),
  d3_qc = file.path(
    root, "analysis_archive", "simulations",
    "D3_regenerated_generator_transport_validation_20260802",
    "output_formal_v1", "audit", "D3_formal_second_pass_QC.csv"
  ),
  d4_oc = file.path(
    root, "analysis_archive", "simulations",
    "D4_longitudinal_G4_G6_validation_20260801",
    "output_formal_v2", "tables", "D4_sim_operating_characteristics.csv"
  ),
  d4_qc = file.path(
    root, "analysis_archive", "simulations",
    "D4_longitudinal_G4_G6_validation_20260801",
    "output_formal_v2", "logs", "D4_sim_independent_QC.csv"
  )
)

missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
if (length(missing)) {
  stop("Missing governed input(s): ", paste(missing, collapse = ", "),
       call. = FALSE)
}

read_qc_pass <- function(path) {
  x <- fread(path)
  if (!all(c("pass") %in% names(x))) return(FALSE)
  all(toupper(as.character(x$pass)) == "TRUE")
}
if (!read_qc_pass(paths$sa01_qc) || !read_qc_pass(paths$sa05_qc)) {
  stop("SA-01 or SA-05 critical QC is not fully passing.", call. = FALSE)
}

registry <- fread(paths$registry)
crosswalk <- fread(paths$crosswalk)
sa01 <- fread(paths$sa01_oc)
sa05 <- fread(paths$sa05_oc)
d2 <- fread(paths$d2_g5)
d3 <- fread(paths$d3_oc)
d4 <- fread(paths$d4_oc)

registry[, display_order := as.integer(sub("G", "", canonical_id))]
setorder(registry, display_order)
fwrite(registry, file.path(out, "Table_G1_G7_canonical_reader_registry.csv"))

domain_map <- rbindlist(list(
  data.table(
    domain = "D1 Discreteness",
    reader_scenario = c("G1", "G7", "G2"),
    canonical_family = c("G1", "G7", "G2"),
    truth_or_role = c(
      "Negative: continuous severity continuum",
      "Positive: aligned two-component discrete structure",
      "Legacy diagnostic: three true components"
    ),
    current_formal_use = c(
      "Included in current Lloyd/N1 three-state verdict",
      "Included in current Lloyd/N1 three-state verdict",
      "Excluded from current fixed-K2 verdict"
    ),
    current_source_id = c(
      "sa01_oc", "sa01_oc",
      "legacy_G2_recovery_only"
    )
  ),
  data.table(
    domain = "D2 Incremental information",
    reader_scenario = c("G1", "G2", "G3", "G4", "G5", "G7"),
    canonical_family = c("G1", "G2", "G3", "G4", "G5", "G7"),
    truth_or_role = c(
      rep("Negative: no label-only prognostic information", 4),
      "Positive: fixed noisy-oracle label",
      "Negative: discrete structure is outcome-null"
    ),
    current_formal_use = c(
      rep("Retain locked module-specific D2 result", 4),
      "Metric-level positive-control accuracy scan",
      "Orthogonal D1-positive/D2-negative control"
    ),
    current_source_id = c(rep("locked_D2_main_harness", 4),
                          "d2_g5", "locked_G7_D2")
  ),
  data.table(
    domain = "D3 Transportability",
    reader_scenario = c("G3a", "G3b", "G3c"),
    canonical_family = "G3",
    truth_or_role = c(
      "Negative: same-distribution application should not alert",
      "Positive: covariate drift should alert",
      "Positive: calibration-interface stress should alert"
    ),
    current_formal_use = "Included in current regenerated-generator D3 validation",
    current_source_id = "d3_oc"
  ),
  data.table(
    domain = "D4 Actionability / support identifiability",
    reader_scenario = c("G4-L", "G6-L", "D4-LQ"),
    canonical_family = c("G4", "G6", "G4/G6"),
    truth_or_role = c(
      "Positive: longitudinal support failure",
      "Negative: longitudinal adequate support",
      "QC only: first-event ledger"
    ),
    current_formal_use = c(
      "Included in current D4 operating characteristics",
      "Included in current D4 operating characteristics",
      "Excluded from operating-characteristic denominators"
    ),
    current_source_id = "d4_oc"
  )
), fill = TRUE)
domain_map[, domain_order := match(
  domain,
  c("D1 Discreteness", "D2 Incremental information",
    "D3 Transportability", "D4 Actionability / support identifiability")
)]
setorder(domain_map, domain_order)
fwrite(domain_map, file.path(out, "Table_G1_G7_by_domain_reader_map.csv"))
fwrite(crosswalk, file.path(out, "Table_current_subdesign_crosswalk.csv"))

sa01_reader <- sa01[, .(
  sensitivity_analysis = "SA-01 sample-size scan",
  sample_size = sprintf("%s derivation / %s evaluation", n_train, n_eval),
  reader_scenario = reader_id,
  truth,
  valid_n,
  failed_n,
  discrete_n = discrete_events,
  no_discrete_n = no_discrete_events,
  inconclusive_n = inconclusive_events,
  strict_expected_state_n = strict_correct_events,
  strict_expected_state_rate = strict_correct_rate,
  strict_wilson_low = strict_correct_low,
  strict_wilson_high = strict_correct_high
)]
fwrite(sa01_reader, file.path(out, "Table_D1_SA01_current_rule_reader_summary.csv"))

sa05_reader <- sa05[, .(
  sensitivity_analysis = "SA-05 missingness mechanism x intensity",
  mechanism = missingness_mechanism,
  target_missingness,
  reader_scenario = reader_id,
  truth,
  valid_n,
  failed_n,
  discrete_n = discrete_events,
  no_discrete_n = no_discrete_events,
  inconclusive_n = inconclusive_events,
  strict_expected_state_n = strict_correct_events,
  strict_expected_state_rate = strict_correct_rate,
  strict_wilson_low = strict_correct_low,
  strict_wilson_high = strict_correct_high
)]
fwrite(sa05_reader, file.path(out, "Table_D1_SA05_current_rule_reader_summary.csv"))

sa01_control <- sa01[, .(
  domain = "D1",
  analysis = "SA-01",
  reader_scenario = reader_id,
  truth = unique(truth),
  n = sum(valid_n),
  failed = sum(failed_n),
  discrete = sum(discrete_events),
  no_discrete = sum(no_discrete_events),
  inconclusive = sum(inconclusive_events),
  strict_expected_state = sum(strict_correct_events),
  note = if (reader_id[1L] == "G1")
    "No false DISCRETE_EVIDENCE alerts; INCONCLUSIVE is not counted as a correct negative" else
    "Aligned fixed-K2 positive control"
), by = reader_id]
sa01_control[, reader_id := NULL]
sa05_control <- sa05[, .(
  domain = "D1",
  analysis = "SA-05",
  reader_scenario = reader_id,
  truth = unique(truth),
  n = sum(valid_n),
  failed = sum(failed_n),
  discrete = sum(discrete_events),
  no_discrete = sum(no_discrete_events),
  inconclusive = sum(inconclusive_events),
  strict_expected_state = sum(strict_correct_events),
  note = if (reader_id[1L] == "G1")
    "No false DISCRETE_EVIDENCE alerts; INCONCLUSIVE varies by missingness condition" else
    "Aligned fixed-K2 positive control"
), by = reader_id]
sa05_control[, reader_id := NULL]

d3_reader <- d3[, .(
  domain = "D3",
  analysis = "Current regenerated-generator controls",
  reader_scenario = fifelse(
    scenario == "G3R_same_distribution_no_alert", "G3a",
    fifelse(scenario == "G3R_covariate_drift_alert", "G3b", "G3c")
  ),
  truth = ifelse(expected_alert, "positive alert", "negative no-alert"),
  n = analyzed_n,
  failed = failed_n,
  discrete = NA_integer_, no_discrete = NA_integer_, inconclusive = NA_integer_,
  strict_expected_state = correct_n,
  note = "D3 alert/no-alert verdict; not a D1 discreteness state"
)]
d4_use <- d4[scenario != "All applicable D4 controls"]
d4_reader <- d4_use[, .(
  domain = "D4",
  analysis = "Current longitudinal support controls",
  reader_scenario = ifelse(grepl("G4", scenario), "G4-L", "G6-L"),
  truth = ifelse(grepl("G4", scenario), "positive alert", "negative no-alert"),
  n = n_valid,
  failed = n_failed,
  discrete = NA_integer_, no_discrete = NA_integer_, inconclusive = NA_integer_,
  strict_expected_state = correct,
  note = "D4 support verdict; not a D1 discreteness state"
)]

cross_domain <- rbindlist(
  list(sa01_control, sa05_control, d3_reader, d4_reader), fill = TRUE
)
cross_domain[, domain_order := match(domain, c("D1", "D3", "D4"))]
cross_domain[, scenario_order := fifelse(
  grepl("^[G]([0-9]+)", reader_scenario),
  as.integer(sub("^G([0-9]+).*$", "\\1", reader_scenario)), 99L
)]
setorder(cross_domain, domain_order, analysis, scenario_order, reader_scenario)
cross_domain[, c("domain_order", "scenario_order") := NULL]
fwrite(cross_domain,
       file.path(out, "Table_current_formal_control_results_cross_domain.csv"))

# Reader-package snapshots; source hashes remain authoritative in the manifest.
fwrite(d2, file.path(out, "Table_D2_G5_oracle_accuracy_summary_locked.csv"))
fwrite(d3, file.path(out, "Table_D3_current_operating_characteristics.csv"))
fwrite(d4, file.path(out, "Table_D4_current_operating_characteristics.csv"))

manifest <- rbindlist(lapply(names(paths), function(id) {
  p <- normalizePath(paths[[id]], winslash = "/", mustWork = TRUE)
  data.table(
    source_id = id,
    absolute_path = p,
    sha256 = digest(p, algo = "sha256", file = TRUE, serialize = FALSE),
    size_bytes = file.info(p)$size,
    modified = format(file.info(p)$mtime, "%Y-%m-%d %H:%M:%S")
  )
}))
fwrite(manifest, file.path(out, "CURRENT_RESULT_SOURCE_MANIFEST.csv"))

readme <- c(
  "# G1-G7 reader numbering and current evidence scope",
  "",
  "The canonical family identifier is always G1-G7. Domain-specific current",
  "implementations use suffixes only when the scientific generator changed:",
  "G3a-G3c for transport controls and G4-L/G6-L for longitudinal D4 controls.",
  "Supplementary figure/table labels must retain their own prefixes (Figure S",
  "or Table S) and must never be shortened to bare S1-S7.",
  "",
  "## Domain order",
  "",
  "1. D1: G1 negative control and G7 positive control. G2 is legacy only.",
  "2. D2: G1-G4 and G7 are negative roles; G5 is the positive control.",
  "3. D3: G3a same-distribution no-alert; G3b drift alert; G3c calibration alert.",
  "4. D4: G4-L support alert; G6-L no support alert; D4-LQ is QC only.",
  "",
  "No cross-domain overall accuracy is calculated here. The domains use",
  "different claim-specific states, and D1 is explicitly three-state."
)
writeLines(readme, file.path(out, "README_G1_G7_READER_NUMBERING.md"))

checklist <- c(
  "# Manuscript and supplement replacement checklist",
  "",
  "- Replace the old SA-01 D1 binary/ARI+SigClust results with",
  "  Table_D1_SA01_current_rule_reader_summary.csv.",
  "- Replace the old SA-05 D1 binary results with",
  "  Table_D1_SA05_current_rule_reader_summary.csv; retain the locked SA-05",
  "  D2 results because this rerun changed only D1.",
  "- State that G1 produced no false DISCRETE_EVIDENCE alerts but frequently",
  "  remained INCONCLUSIVE; do not relabel INCONCLUSIVE as a correct negative.",
  "- Use G7, not G2, as the aligned fixed-K2 D1 positive control.",
  "- Label G2 as a legacy three-component recovery diagnostic outside the",
  "  current fixed-K2 D1 verdict.",
  "- Use G3a/G3b/G3c for current D3 controls and G4-L/G6-L for current D4",
  "  controls. Keep legacy G3/G4/G6 output only in provenance/deprecation text.",
  "- Any old Figure 7 or operating-characteristics table that assumes binary",
  "  D1 correctness or includes G2 as the current D1 positive control is",
  "  obsolete and must not be submitted without regeneration.",
  "- Do not report a single pooled four-domain accuracy from these mixed",
  "  claim-specific states."
)
writeLines(checklist, file.path(out, "MANUSCRIPT_REPLACEMENT_CHECKLIST.md"))

cat("Integrated reader outputs written to: ", out, "\n", sep = "")
