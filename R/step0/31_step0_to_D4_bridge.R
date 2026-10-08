#!/usr/bin/env Rscript

## WP10-B: Domain 1 identity to Domain 4 RRT positivity/overlap bridge.
## Descriptive identifiability audit only. No treatment-effect model is fitted.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(pROC)
  library(ggplot2)
  library(patchwork)
})

PROJECT_ROOT <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
AUDIT_ROOT <- file.path(
  PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"
)
SCRIPT_PATH <- file.path(AUDIT_ROOT, "scripts", "37_domain1_to_domain4_bridge.R")
PRERUN_SPEC <- file.path(
  AUDIT_ROOT, "provenance", "36_37_WP10_prerun_spec_FROZEN_20260718.md"
)
WP10_ROOT <- file.path(AUDIT_ROOT, "wp10_domain2_domain4_bridge_20260718")
RUN_ROOT <- file.path(WP10_ROOT, "domain4")
TABLE_DIR <- file.path(RUN_ROOT, "tables")
FIGURE_DIR <- file.path(RUN_ROOT, "figures")
LOG_DIR <- file.path(RUN_ROOT, "logs")
PROVENANCE_DIR <- file.path(RUN_ROOT, "provenance")
for (d in c(RUN_ROOT, TABLE_DIR, FIGURE_DIR, LOG_DIR, PROVENANCE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

INPUT <- c(
  script = SCRIPT_PATH,
  frozen_spec = PRERUN_SPEC,
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  baseline = file.path(PROJECT_ROOT, "data", "baseline_covars.csv"),
  rrt_first = file.path(PROJECT_ROOT, "output", "rrt_first.csv"),
  rrt_base_times = file.path(PROJECT_ROOT, "output", "rrt_hte_base_times.csv"),
  rrt_sql = file.path(PROJECT_ROOT, "sql", "00_extract_rrt_hte_features_abs.sql"),
  historical_source_map = file.path(AUDIT_ROOT, "provenance", "D1_source_map.csv"),
  missing_criterion_audit = file.path(
    AUDIT_ROOT, "tables", "Table_D1_AKI_criterion_by_cluster.csv"
  ),
  missing_criterion_log = file.path(
    AUDIT_ROOT, "logs", "29a_domain1_aki_alignment_log.md"
  )
)
for (p in INPUT) {
  if (!file.exists(p)) stop("Missing required WP10-D4 input: ", p, call. = FALSE)
}

spec_text <- readLines(PRERUN_SPEC, warn = FALSE, encoding = "UTF-8")
if (!any(spec_text == "Status: **FROZEN - AUTHORIZED FOR FORMAL EXECUTION**")) {
  stop("Frozen WP10 specification lacks formal authorization marker.", call. = FALSE)
}

COMPLETION <- file.path(RUN_ROOT, "37_WP10_D4_run_completed_with_declared_missing.ok")
if (file.exists(COMPLETION)) {
  stop("WP10-D4 output is already complete and will not be overwritten: ", RUN_ROOT,
       call. = FALSE)
}

EXPECTED_SHA256 <- c(
  labels = "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
  final_full = "CE0B38D2D3D6263A5C4615A6AEA1E78B19840BD4B35046DBFA65EE6BFBBA633B",
  baseline = "3604A7AD95B0F3E2AF4CAA632C753C11E0C56CF0CC3BB5B6AB8DE984D85DB107",
  rrt_first = "5D50E32F9874F44F0D67E4D06331236162B056642B50C8ED16514429F7102F2F",
  rrt_base_times = "0E7AD66D2D7D9BF25209AAA9047445E1D2EB374B99D8B7B7C5B520A1BBB41C58"
)

sha256_file <- function(path) {
  toupper(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(x, tmp, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path, call. = FALSE)
}

atomic_write_lines <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic text write failed: ", path, call. = FALSE)
}

require_columns <- function(dt, cols, label) {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) {
    stop(label, " lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

input_manifest <- rbindlist(lapply(names(INPUT), function(nm) {
  p <- INPUT[[nm]]
  hash <- sha256_file(p)
  data.table(
    input_name = nm,
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    sha256 = hash,
    expected_sha256 = if (nm %in% names(EXPECTED_SHA256)) EXPECTED_SHA256[[nm]] else NA_character_,
    hash_match = if (nm %in% names(EXPECTED_SHA256)) hash == EXPECTED_SHA256[[nm]] else NA,
    bytes = file.info(p)$size
  )
}))
locked_rows <- input_manifest[input_name %in% names(EXPECTED_SHA256)]
if (any(!locked_rows$hash_match)) {
  stop(
    "At least one locked WP10-D4 input hash does not match:\n",
    paste(capture.output(print(locked_rows[hash_match == FALSE])), collapse = "\n"),
    call. = FALSE
  )
}
atomic_fwrite(
  input_manifest,
  file.path(PROVENANCE_DIR, "37_WP10_D4_input_sha256_manifest.csv")
)

# Preserve, rather than silently edit, the historical truncated hash entry.
source_map <- fread(INPUT[["historical_source_map"]])
old_rrt_hash <- source_map[source_id == "rrt_first", sha256]
if (length(old_rrt_hash) != 1L) old_rrt_hash <- NA_character_
hash_note <- c(
  "# RRT hash provenance correction note",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("Historical D1_source_map.csv value (", nchar(old_rrt_hash), " characters): `", old_rrt_hash, "`"),
  paste0("Current verified SHA-256 (64 characters): `", sha256_file(INPUT[["rrt_first"]]), "`"),
  "",
  paste(
    "The historical entry omitted the terminal `2F`. The source map is retained",
    "unchanged for provenance; this WP10 manifest records the complete hash."
  )
)
atomic_write_lines(hash_note, file.path(PROVENANCE_DIR, "37_RRT_hash_correction_note.md"))

sql_text <- paste(readLines(INPUT[["rrt_sql"]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
sql_time_origin_verified <- grepl(
  "charttime\\s*-\\s*c\\.intime|r\\.charttime\\s*-\\s*c\\.intime",
  sql_text, ignore.case = TRUE, perl = TRUE
)
if (!sql_time_origin_verified) {
  stop("Could not verify first_rrt_hours relative to ICU intime in SQL.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Assemble the locked cohort and independently audit RRT timing
# -------------------------------------------------------------------------

labels <- as.data.table(readRDS(INPUT[["labels"]]))
final_needed <- c(
  "stay_id", "gender", "age", "sofa_score", "creatinine_max", "bun_max",
  "potassium_max", "bicarbonate_min", "ph_min", "lactate_max",
  "urine_output_24h_ml", "mbp_min", "heart_rate_max", "resp_rate_max", "gcs_min"
)
final <- fread(INPUT[["final_full"]], select = final_needed)
baseline <- fread(INPUT[["baseline"]])
rrt <- fread(INPUT[["rrt_first"]])
times <- fread(INPUT[["rrt_base_times"]])

require_columns(labels, c("stay_id", "cluster_k2"), "Locked labels")
require_columns(final, final_needed, "final_full")
require_columns(
  baseline,
  c("stay_id", "charlson_comorbidity_index", "sapsii", "aki_stage_0_24h"),
  "baseline_covars"
)
require_columns(rrt, c("stay_id", "first_rrt_hours"), "rrt_first")
require_columns(
  times,
  c("stay_id", "first_rrt_hours", "death_time_hours"),
  "rrt_hte_base_times"
)

objects <- list(labels = labels, final = final, baseline = baseline, rrt = rrt, times = times)
for (nm in names(objects)) {
  x <- objects[[nm]]
  if (anyDuplicated(x$stay_id)) stop(nm, " contains duplicated stay_id.", call. = FALSE)
}
if (nrow(labels) != 20049L || nrow(final) != 20049L || nrow(times) != 20049L) {
  stop("Expected 20,049 rows in labels/final/time anchors.", call. = FALSE)
}
if (!all(labels$stay_id %in% baseline$stay_id)) {
  stop("The upstream baseline table does not cover every locked stay_id.", call. = FALSE)
}

setnames(rrt, "first_rrt_hours", "first_rrt_hours_rrt_file")
setnames(times, "first_rrt_hours", "first_rrt_hours_time_file")
dt <- merge(labels, final, by = "stay_id", all.x = TRUE, sort = FALSE)
dt <- merge(
  dt,
  baseline[, .(stay_id, charlson_comorbidity_index, sapsii, aki_stage_0_24h)],
  by = "stay_id", all.x = TRUE, sort = FALSE
)
dt <- merge(
  dt,
  times[, .(stay_id, first_rrt_hours_time_file, death_time_hours)],
  by = "stay_id", all.x = TRUE, sort = FALSE
)
dt <- merge(
  dt,
  rrt[, .(stay_id, first_rrt_hours_rrt_file)],
  by = "stay_id", all.x = TRUE, sort = FALSE
)
if (nrow(dt) != 20049L || anyNA(dt$cluster_k2)) {
  stop("Cohort merge failed the locked 20,049-row label check.", call. = FALSE)
}

for (v in c("first_rrt_hours_time_file", "first_rrt_hours_rrt_file", "death_time_hours")) {
  set(dt, j = v, value = suppressWarnings(as.numeric(dt[[v]])))
}
rrt_compare <- dt[!is.na(first_rrt_hours_rrt_file), .(
  n = .N,
  missing_in_time_file = sum(is.na(first_rrt_hours_time_file)),
  max_abs_difference = max(abs(first_rrt_hours_rrt_file - first_rrt_hours_time_file), na.rm = TRUE)
)]
if (rrt_compare$missing_in_time_file != 0L ||
    !is.finite(rrt_compare$max_abs_difference) ||
    rrt_compare$max_abs_difference > 1e-8) {
  stop("RRT timing differs between rrt_first.csv and base-time audit file.", call. = FALSE)
}
dt[, first_rrt_hours := first_rrt_hours_time_file]
dt[, phenotype := factor(
  cluster_k2,
  levels = c(1L, 2L),
  labels = c("C1 higher-risk", "C2 lower-risk")
)]
dt[, aki_stage_0_24h := suppressWarnings(as.integer(aki_stage_0_24h))]

# -------------------------------------------------------------------------
# RRT rates by combined KDIGO stage
# -------------------------------------------------------------------------

build_population <- function(data, mode, eligibility) {
  origin <- if (mode == "window_start") 24 else 72
  d0 <- data[is.na(first_rrt_hours) | first_rrt_hours >= 0]
  d1 <- d0[is.na(first_rrt_hours) | first_rrt_hours >= 24]
  d2 <- d1[is.na(death_time_hours) | death_time_hours >= origin]
  d3 <- if (eligibility == "KDIGO_2_3") {
    d2[aki_stage_0_24h %in% c(2L, 3L)]
  } else {
    copy(d2)
  }
  d3[, treat_24_72 := as.integer(
    !is.na(first_rrt_hours) & first_rrt_hours >= 24 & first_rrt_hours < 72
  )]
  flow <- data.table(
    risk_set_mode = mode,
    eligibility = eligibility,
    risk_origin_h = origin,
    n_locked_cohort = nrow(data),
    excluded_pre_icu_rrt = nrow(data) - nrow(d0),
    excluded_rrt_before_24h = nrow(d0) - nrow(d1),
    excluded_death_before_origin = nrow(d1) - nrow(d2),
    excluded_by_eligibility = nrow(d2) - nrow(d3),
    n_analyzed = nrow(d3)
  )
  list(data = d3, flow = flow)
}

stage_label <- function(x) {
  fifelse(is.na(x), "Missing", paste0("KDIGO stage ", x))
}

descriptive <- dt[is.na(first_rrt_hours) | first_rrt_hours >= 0]
descriptive[, `:=`(
  rrt_any_post_icu = as.integer(!is.na(first_rrt_hours) & first_rrt_hours >= 0),
  rrt_0_24 = as.integer(!is.na(first_rrt_hours) & first_rrt_hours >= 0 & first_rrt_hours < 24),
  rrt_24_72 = as.integer(!is.na(first_rrt_hours) & first_rrt_hours >= 24 & first_rrt_hours < 72),
  kdigo_label = stage_label(aki_stage_0_24h)
)]
rrt_stage_descriptive <- descriptive[, .(
  n = .N,
  rrt_any_post_icu_n = sum(rrt_any_post_icu),
  rrt_any_post_icu_rate = mean(rrt_any_post_icu),
  rrt_0_24_n = sum(rrt_0_24),
  rrt_0_24_rate = mean(rrt_0_24),
  rrt_24_72_n = sum(rrt_24_72),
  rrt_24_72_rate = mean(rrt_24_72)
), by = .(kdigo_label, aki_stage_0_24h)]
rrt_stage_descriptive[, `:=`(
  analysis_population = "descriptive_no_preICU_RRT",
  risk_set_mode = "not_applicable",
  risk_origin_h = NA_real_
)]

population_grid <- CJ(
  risk_set_mode = c("window_start", "window_end"),
  eligibility = c("full", "KDIGO_2_3"),
  sorted = FALSE
)
populations <- lapply(seq_len(nrow(population_grid)), function(i) {
  build_population(
    dt,
    population_grid$risk_set_mode[i],
    population_grid$eligibility[i]
  )
})
names(populations) <- paste(
  population_grid$risk_set_mode,
  population_grid$eligibility,
  sep = "__"
)
flow_table <- rbindlist(lapply(populations, `[[`, "flow"), fill = TRUE)
atomic_fwrite(flow_table, file.path(TABLE_DIR, "37_WP10_D4_cohort_flow.csv"))

rrt_stage_risk <- rbindlist(lapply(names(populations), function(nm) {
  d <- populations[[nm]]$data
  d[, kdigo_label := stage_label(aki_stage_0_24h)]
  ans <- d[, .(
    n = .N,
    rrt_any_post_icu_n = NA_integer_,
    rrt_any_post_icu_rate = NA_real_,
    rrt_0_24_n = NA_integer_,
    rrt_0_24_rate = NA_real_,
    rrt_24_72_n = sum(treat_24_72),
    rrt_24_72_rate = mean(treat_24_72)
  ), by = .(kdigo_label, aki_stage_0_24h)]
  parts <- strsplit(nm, "__", fixed = TRUE)[[1L]]
  ans[, `:=`(
    analysis_population = parts[2L],
    risk_set_mode = parts[1L],
    risk_origin_h = ifelse(parts[1L] == "window_start", 24, 72)
  )]
  ans
}), fill = TRUE)

rrt_by_kdigo <- rbindlist(list(rrt_stage_descriptive, rrt_stage_risk), fill = TRUE)
setcolorder(rrt_by_kdigo, c(
  "analysis_population", "risk_set_mode", "risk_origin_h",
  "aki_stage_0_24h", "kdigo_label", "n",
  "rrt_any_post_icu_n", "rrt_any_post_icu_rate",
  "rrt_0_24_n", "rrt_0_24_rate", "rrt_24_72_n", "rrt_24_72_rate"
))
atomic_fwrite(rrt_by_kdigo, file.path(TABLE_DIR, "Table_D1_D4_RRT_by_KDIGO.csv"))

criterion_table <- data.table(
  aki_criterion = c("creatinine_only", "urine_only", "both", "neither"),
  status = "MISSING_NOT_RECONSTRUCTED",
  n = NA_integer_,
  rrt_24_72_n = NA_integer_,
  rrt_24_72_rate = NA_real_,
  source = "No authoritative 0-24 h component-criterion export",
  reason = paste(
    "Combined KDIGO stage cannot be decomposed without criterion-level source fields;",
    "proxy reconstruction is prohibited."
  )
)
atomic_fwrite(
  criterion_table,
  file.path(TABLE_DIR, "Table_D1_D4_RRT_by_AKI_criterion.csv")
)

# -------------------------------------------------------------------------
# Propensity and overlap diagnostics
# -------------------------------------------------------------------------

PS_COVARIATES <- c(
  "age", "gender", "sofa_score", "aki_stage_0_24h", "sapsii",
  "creatinine_max", "bun_max", "potassium_max", "bicarbonate_min",
  "ph_min", "lactate_max", "urine_output_24h_ml", "mbp_min",
  "heart_rate_max", "resp_rate_max", "gcs_min"
)

median_mode_impute <- function(data, vars) {
  d <- copy(data)
  for (v in vars) {
    x <- d[[v]]
    if (is.numeric(x) || is.integer(x)) {
      med <- suppressWarnings(stats::median(x, na.rm = TRUE))
      if (!is.finite(med)) stop("Cannot impute numeric covariate: ", v, call. = FALSE)
      x[is.na(x)] <- med
    } else {
      tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
      if (!length(tab)) stop("Cannot impute categorical covariate: ", v, call. = FALSE)
      x[is.na(x) | !nzchar(as.character(x))] <- names(tab)[1L]
      x <- factor(x)
    }
    set(d, j = v, value = x)
  }
  d
}

fit_overlap <- function(pop, mode, eligibility) {
  d <- copy(pop)
  covars <- PS_COVARIATES[vapply(PS_COVARIATES, function(v) {
    v %in% names(d) && uniqueN(d[[v]], na.rm = TRUE) > 1L
  }, logical(1))]
  if (length(covars) < 3L) stop("Too few nonconstant PS covariates.", call. = FALSE)
  if (sum(d$treat_24_72 == 1L) < 5L || sum(d$treat_24_72 == 0L) < 5L) {
    stop("Insufficient treated or untreated observations for propensity model.", call. = FALSE)
  }
  dm <- median_mode_impute(d, covars)
  warning_text <- character()
  fit <- withCallingHandlers(
    glm(
      reformulate(covars, response = "treat_24_72"),
      data = dm, family = binomial()
    ),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(fit$converged)) stop("Propensity model did not converge.", call. = FALSE)
  ps_raw <- as.numeric(stats::predict(fit, type = "response"))
  eps <- 1e-6
  ps <- pmin(pmax(ps_raw, eps), 1 - eps)
  dm[, ps := ps]
  roc_obj <- pROC::roc(
    response = dm$treat_24_72,
    predictor = dm$ps,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  c_stat <- as.numeric(pROC::auc(roc_obj))
  treated_ps <- dm[treat_24_72 == 1L, ps]
  untreated_ps <- dm[treat_24_72 == 0L, ps]
  support_low <- max(min(treated_ps), min(untreated_ps))
  support_high <- min(max(treated_ps), max(untreated_ps))
  dm[, outside_overlap := ps < support_low | ps > support_high]
  dm[, ipw_ate := fifelse(treat_24_72 == 1L, 1 / ps, 1 / (1 - ps))]
  outside <- mean(dm$outside_overlap)
  ess <- sum(dm$ipw_ate)^2 / sum(dm$ipw_ate^2)
  ess_ratio <- ess / nrow(dm)
  extreme <- mean(dm$ipw_ate > 10)
  cells <- merge(
    CJ(cluster_k2 = c(1L, 2L), treat_24_72 = c(0L, 1L)),
    dm[, .N, by = .(cluster_k2, treat_24_72)],
    by = c("cluster_k2", "treat_24_72"), all.x = TRUE
  )
  cells[is.na(N), N := 0L]
  min_cell <- min(cells$N)
  counts <- dm[, .(n = .N, treated_n = sum(treat_24_72)), by = cluster_k2]
  c1 <- counts[cluster_k2 == 1L]
  c2 <- counts[cluster_k2 == 2L]
  alerts <- c(
    outside_overlap = outside > 0.10,
    ess_ratio = ess_ratio < 0.25,
    minimum_cell = min_cell < 30,
    extreme_weight = extreme > 0.01
  )
  summary <- data.table(
    risk_set_mode = mode,
    risk_origin_h = ifelse(mode == "window_start", 24, 72),
    eligibility = eligibility,
    n = nrow(dm),
    treated_n = sum(dm$treat_24_72),
    untreated_n = sum(dm$treat_24_72 == 0L),
    treatment_prevalence = mean(dm$treat_24_72),
    c1_n = c1$n,
    c1_treated_n = c1$treated_n,
    c1_treatment_prevalence = c1$treated_n / c1$n,
    c2_n = c2$n,
    c2_treated_n = c2$treated_n,
    c2_treatment_prevalence = c2$treated_n / c2$n,
    ps_c_statistic = c_stat,
    common_support_low = support_low,
    common_support_high = support_high,
    outside_overlap_proportion = outside,
    ate_ipw_ess = ess,
    ess_ratio = ess_ratio,
    minimum_phenotype_treatment_cell = min_cell,
    extreme_weight_threshold = 10,
    extreme_weight_proportion = extreme,
    outside_overlap_alert = alerts[["outside_overlap"]],
    ess_ratio_alert = alerts[["ess_ratio"]],
    minimum_cell_alert = alerts[["minimum_cell"]],
    extreme_weight_alert = alerts[["extreme_weight"]],
    identifiability = if (any(alerts)) "non_identifiable" else "adequate_by_prespecified_gates",
    ps_covariates = paste(covars, collapse = ";"),
    ps_n_covariates = length(covars),
    ps_converged = fit$converged,
    ps_warning_count = length(warning_text),
    ps_warning_text = paste(unique(warning_text), collapse = " | "),
    effect_model_fitted = FALSE
  )
  patient <- dm[, .(
    stay_id, risk_set_mode = mode, eligibility, cluster_k2,
    phenotype, treat_24_72, aki_stage_0_24h, ps, outside_overlap, ipw_ate
  )]
  list(summary = summary, patient = patient, cells = cells)
}

overlap_results <- lapply(names(populations), function(nm) {
  parts <- strsplit(nm, "__", fixed = TRUE)[[1L]]
  fit_overlap(populations[[nm]]$data, parts[1L], parts[2L])
})

overlap_summary <- rbindlist(lapply(overlap_results, `[[`, "summary"), fill = TRUE)
overlap_patients <- rbindlist(lapply(overlap_results, `[[`, "patient"), fill = TRUE)
cell_table <- rbindlist(lapply(seq_along(overlap_results), function(i) {
  x <- copy(overlap_results[[i]]$cells)
  x[, `:=`(
    risk_set_mode = overlap_results[[i]]$summary$risk_set_mode,
    eligibility = overlap_results[[i]]$summary$eligibility
  )]
  x
}), fill = TRUE)

overlap_summary[, `:=`(
  mode_order = match(risk_set_mode, c("window_start", "window_end")),
  eligibility_order = match(eligibility, c("full", "KDIGO_2_3"))
)]
setorder(overlap_summary, mode_order, eligibility_order)
overlap_summary[, c("mode_order", "eligibility_order") := NULL]
atomic_fwrite(
  overlap_summary,
  file.path(TABLE_DIR, "Table_D1_D4_restricted_overlap.csv")
)
atomic_fwrite(
  overlap_patients,
  file.path(TABLE_DIR, "37_WP10_D4_propensity_patient_level.csv")
)
atomic_fwrite(
  cell_table,
  file.path(TABLE_DIR, "37_WP10_D4_phenotype_treatment_cells.csv")
)

# -------------------------------------------------------------------------
# Figure: treatment prevalence, PS distributions, gates, and PS C-statistic
# -------------------------------------------------------------------------

phen_long <- melt(
  overlap_summary,
  id.vars = c("risk_set_mode", "eligibility"),
  measure.vars = c("c1_treatment_prevalence", "c2_treatment_prevalence"),
  variable.name = "phenotype", value.name = "estimate"
)
phen_long[, phenotype := factor(
  phenotype,
  levels = c("c1_treatment_prevalence", "c2_treatment_prevalence"),
  labels = c("C1 higher-risk", "C2 lower-risk")
)]
phen_long[, population := factor(
  paste(risk_set_mode, eligibility, sep = "__"),
  levels = c(
    "window_start__full", "window_start__KDIGO_2_3",
    "window_end__full", "window_end__KDIGO_2_3"
  ),
  labels = c(
    "24-h risk set\nFull", "24-h risk set\nKDIGO 2-3",
    "72-h survivors\nFull", "72-h survivors\nKDIGO 2-3"
  )
)]

ps_plot <- overlap_patients[risk_set_mode == "window_start"]
ps_plot[, treatment := factor(
  treat_24_72, levels = c(0, 1), labels = c("No RRT", "RRT 24-72 h")
)]
ps_plot[, eligibility := factor(
  eligibility, levels = c("full", "KDIGO_2_3"),
  labels = c("Full cohort", "KDIGO stage 2-3")
)]

gate_long <- rbindlist(list(
  overlap_summary[, .(
    risk_set_mode, eligibility, metric = "Outside overlap",
    estimate = outside_overlap_proportion, threshold = 0.10,
    alert_ratio = outside_overlap_proportion / 0.10
  )],
  overlap_summary[, .(
    risk_set_mode, eligibility, metric = "ESS/N",
    estimate = ess_ratio, threshold = 0.25,
    alert_ratio = 0.25 / ess_ratio
  )],
  overlap_summary[, .(
    risk_set_mode, eligibility, metric = "Minimum cell",
    estimate = minimum_phenotype_treatment_cell, threshold = 30,
    alert_ratio = 30 / minimum_phenotype_treatment_cell
  )],
  overlap_summary[, .(
    risk_set_mode, eligibility, metric = "Extreme weights",
    estimate = extreme_weight_proportion, threshold = 0.01,
    alert_ratio = extreme_weight_proportion / 0.01
  )]
), fill = TRUE)
gate_long[alert_ratio == 0, alert_ratio := 1e-3]
gate_long[, alert := alert_ratio > 1]
gate_long[, population := factor(
  paste(risk_set_mode, eligibility, sep = "__"),
  levels = c(
    "window_start__full", "window_start__KDIGO_2_3",
    "window_end__full", "window_end__KDIGO_2_3"
  ),
  labels = c(
    "24-h/full", "24-h/KDIGO 2-3", "72-h/full", "72-h/KDIGO 2-3"
  )
)]

cstat_plot <- overlap_summary[, .(
  risk_set_mode, eligibility, ps_c_statistic,
  population = factor(
    paste(risk_set_mode, eligibility, sep = "__"),
    levels = c(
      "window_start__full", "window_start__KDIGO_2_3",
      "window_end__full", "window_end__KDIGO_2_3"
    ),
    labels = c(
      "24-h/full", "24-h/KDIGO 2-3", "72-h/full", "72-h/KDIGO 2-3"
    )
  )
)]

figure_source <- rbindlist(list(
  phen_long[, .(
    panel = "A", risk_set_mode, eligibility,
    group = as.character(phenotype), x = as.character(population),
    estimate, threshold = NA_real_, alert_ratio = NA_real_
  )],
  ps_plot[, .(
    panel = "B", risk_set_mode, eligibility = as.character(eligibility),
    group = paste(phenotype, treatment, sep = ";"), x = as.character(stay_id),
    estimate = ps, threshold = NA_real_, alert_ratio = NA_real_
  )],
  gate_long[, .(
    panel = "C", risk_set_mode, eligibility, group = metric,
    x = as.character(population), estimate, threshold, alert_ratio
  )],
  cstat_plot[, .(
    panel = "D", risk_set_mode, eligibility, group = "PS C-statistic",
    x = as.character(population), estimate = ps_c_statistic,
    threshold = NA_real_, alert_ratio = NA_real_
  )]
), fill = TRUE)
atomic_fwrite(
  figure_source,
  file.path(TABLE_DIR, "37_Figure_D1_D4_restricted_overlap_source_data.csv")
)

theme_pub <- function() {
  theme_classic(base_family = "sans", base_size = 9.5) +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, colour = "grey35"),
      plot.tag = element_text(size = 11, face = "bold"),
      strip.text = element_text(size = 8.5, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin = margin(6, 8, 6, 6)
    )
}

phenotype_colours <- c(
  "C1 higher-risk" = "#B4553F",
  "C2 lower-risk" = "#3E6E93"
)

p_a <- ggplot(phen_long, aes(x = population, y = estimate, fill = phenotype)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = phenotype_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Subsequent RRT prevalence",
    subtitle = "First RRT from 24 h to <72 h",
    x = NULL, y = "Treatment prevalence"
  ) + theme_pub()

p_b <- ggplot(
  ps_plot,
  aes(x = ps, colour = phenotype, linetype = treatment)
) +
  geom_density(linewidth = 0.65, adjust = 1.1) +
  facet_wrap(~eligibility, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = phenotype_colours) +
  scale_linetype_manual(values = c("No RRT" = "solid", "RRT 24-72 h" = "22")) +
  labs(
    title = "Propensity-score distributions",
    subtitle = "Primary 24-h risk set; curves are descriptive",
    x = "Estimated treatment propensity", y = "Density"
  ) + theme_pub()

p_c <- ggplot(
  gate_long,
  aes(x = population, y = alert_ratio, shape = metric, colour = alert)
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  geom_point(size = 2.1, position = position_jitter(width = 0.10, height = 0)) +
  scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "#3A5A8C", `TRUE` = "#C77A30"),
                      labels = c(`FALSE` = "Within gate", `TRUE` = "Alert")) +
  labs(
    title = "Prespecified overlap gates",
    subtitle = "Values above 1 trigger the corresponding alert",
    x = NULL, y = "Alert ratio (log scale)"
  ) + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_d <- ggplot(cstat_plot, aes(x = population, y = ps_c_statistic)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2.2, colour = "#2E2E2E") +
  scale_y_continuous(limits = c(0.5, 1), breaks = c(0.5, 0.75, 1)) +
  labs(
    title = "Treatment-model discrimination",
    subtitle = "Higher values indicate stronger indication-driven assignment",
    x = NULL, y = "PS C-statistic"
  ) + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_base <- file.path(FIGURE_DIR, "Figure_D1_D4_restricted_overlap")
ggsave(paste0(fig_base, ".pdf"), fig, width = 183 / 25.4, height = 145 / 25.4,
       device = cairo_pdf, family = "sans", bg = "white")
ggsave(paste0(fig_base, ".png"), fig, width = 183 / 25.4, height = 145 / 25.4,
       units = "in", dpi = 600, bg = "white")
grDevices::tiff(
  paste0(fig_base, ".tiff"), width = 183 / 25.4, height = 145 / 25.4,
  units = "in", res = 600, compression = "lzw", type = "cairo", bg = "white"
)
print(fig)
grDevices::dev.off()

# -------------------------------------------------------------------------
# QC, log, manifests, and declared-missing completion
# -------------------------------------------------------------------------

criterion_source <- fread(INPUT[["missing_criterion_audit"]])
criterion_missing_confirmed <- nrow(criterion_source) == 4L &&
  any(grepl("MISSING_NOT_RECONSTRUCTED", unlist(criterion_source), fixed = TRUE))

qc <- rbindlist(list(
  data.table(check = "locked_input_hashes_match", pass = all(locked_rows$hash_match),
             observed = sum(locked_rows$hash_match), expected = nrow(locked_rows)),
  data.table(check = "historical_rrt_hash_discrepancy_preserved",
             pass = nchar(old_rrt_hash) == 62L &&
               nchar(sha256_file(INPUT[["rrt_first"]])) == 64L,
             observed = paste(nchar(old_rrt_hash), 64L, sep = "/"), expected = "62/64"),
  data.table(check = "rrt_sql_time_origin_verified", pass = sql_time_origin_verified,
             observed = sql_time_origin_verified, expected = TRUE),
  data.table(check = "rrt_files_agree", pass = rrt_compare$missing_in_time_file == 0L &&
               rrt_compare$max_abs_difference <= 1e-8,
             observed = paste(rrt_compare$missing_in_time_file,
                              signif(rrt_compare$max_abs_difference, 4), sep = "/"),
             expected = "0/<=1e-8"),
  data.table(check = "locked_cohort_and_labels", pass = nrow(dt) == 20049L &&
               all(table(dt$cluster_k2) == c(`1` = 3992L, `2` = 16057L)),
             observed = paste(nrow(dt), paste(table(dt$cluster_k2), collapse = "/"), sep = ";"),
             expected = "20049;3992/16057"),
  data.table(check = "four_population_rows", pass = nrow(overlap_summary) == 4L,
             observed = nrow(overlap_summary), expected = 4L),
  data.table(check = "propensity_models_converged", pass = all(overlap_summary$ps_converged),
             observed = sum(overlap_summary$ps_converged), expected = 4L),
  data.table(check = "finite_overlap_metrics", pass = all(is.finite(unlist(
    overlap_summary[, .(
      ps_c_statistic, outside_overlap_proportion, ess_ratio,
      minimum_phenotype_treatment_cell, extreme_weight_proportion
    )]
  ))), observed = "all four populations", expected = "finite"),
  data.table(check = "no_effect_model_fitted", pass = all(!overlap_summary$effect_model_fitted),
             observed = sum(overlap_summary$effect_model_fitted), expected = 0L),
  data.table(check = "criterion_missing_confirmed_and_not_reconstructed",
             pass = criterion_missing_confirmed &&
               all(criterion_table$status == "MISSING_NOT_RECONSTRUCTED") &&
               all(is.na(criterion_table$rrt_24_72_rate)),
             observed = paste(criterion_missing_confirmed,
                              all(criterion_table$status == "MISSING_NOT_RECONSTRUCTED"),
                              all(is.na(criterion_table$rrt_24_72_rate)), sep = "/"),
             expected = "TRUE/TRUE/TRUE"),
  data.table(check = "figure_exports_exist",
             pass = all(file.exists(paste0(fig_base, c(".pdf", ".png", ".tiff")))),
             observed = sum(file.exists(paste0(fig_base, c(".pdf", ".png", ".tiff")))),
             expected = 3L)
), fill = TRUE)
atomic_fwrite(qc, file.path(LOG_DIR, "37_WP10_D4_QC.csv"))

log_lines <- c(
  "# WP10 Domain 1 to Domain 4 bridge log",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Time ordering",
  "- Phenotype/KDIGO features: ICU hours 0-24.",
  "- Treatment: first RRT from 24 h inclusive to 72 h exclusive.",
  "- Primary population: alive and RRT-free at 24 h (`window_start`).",
  paste(
    "- Sensitivity population: alive at 72 h and RRT-free at the 24-h window",
    "start (`window_end` survivor-to-landmark comparison)."
  ),
  "",
  "## Interpretation boundary",
  "- This is a positivity/identifiability audit only.",
  "- No treatment-effect model was fitted.",
  "- The analysis cannot establish presence or absence of true HTE.",
  "",
  "## Declared missing component",
  paste(
    "- Creatinine-only/urine-only/both/neither status is unavailable in an",
    "authoritative export and was not reconstructed or proxied."
  ),
  "",
  "## RRT hash correction",
  paste0("- Historical source-map value: ", old_rrt_hash),
  paste0("- Current verified full SHA-256: ", sha256_file(INPUT[["rrt_first"]]))
)
atomic_write_lines(log_lines, file.path(LOG_DIR, "37_D4_bridge_log.md"))
capture.output(sessionInfo(), file = file.path(PROVENANCE_DIR, "37_WP10_D4_sessionInfo.txt"))

if (!all(qc$pass)) {
  atomic_write_lines(c(
    "WP10-D4 QC FAILURE - DO NOT USE",
    paste0("failed_checks=", paste(qc[pass == FALSE, check], collapse = ";"))
  ), file.path(RUN_ROOT, "37_WP10_D4_QC_FAILURE.txt"))
  stop("WP10-D4 completed computation but failed structural QC.", call. = FALSE)
}

output_files <- list.files(RUN_ROOT, recursive = TRUE, full.names = TRUE)
output_files <- output_files[file.info(output_files)$isdir == FALSE]
output_files <- output_files[!grepl(
  "37_WP10_D4_output_sha256_manifest.csv$|37_WP10_D4_run_completed_with_declared_missing.ok$",
  output_files
)]
output_manifest <- rbindlist(lapply(output_files, function(p) {
  data.table(
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    sha256 = sha256_file(p),
    bytes = file.info(p)$size
  )
}))
atomic_fwrite(
  output_manifest,
  file.path(PROVENANCE_DIR, "37_WP10_D4_output_sha256_manifest.csv")
)

atomic_write_lines(c(
  "RUN COMPLETED WITH DECLARED MISSING INPUT",
  "work_package=WP10_Domain4_bridge",
  "run_mode=formal_deterministic_postprocess",
  "treatment_window=24_to_less_than_72_hours",
  "primary_risk_set=alive_and_RRT_free_at_24h",
  "sensitivity_risk_set=alive_at_72h_and_RRT_free_at_24h",
  "effect_model_fitted=FALSE",
  "criterion_decomposition=MISSING_NOT_RECONSTRUCTED",
  paste0("n_overlap_rows=", nrow(overlap_summary)),
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0(
    "output_manifest_sha256=",
    sha256_file(file.path(PROVENANCE_DIR, "37_WP10_D4_output_sha256_manifest.csv"))
  )
), COMPLETION)

print(overlap_summary)
