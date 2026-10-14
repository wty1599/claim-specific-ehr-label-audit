source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  "R", "domain4", "landmark",
  "38_d4_landmark_support_common_v1.R"
))

stage <- "38a_source_audit"
if (file.exists(d4_stage_marker(stage))) {
  stop("Stage already completed: ", stage, call. = FALSE)
}

d4_msg("Starting ", stage, " in ", D4_RUN_MODE, " mode")

manifest <- rbindlist(lapply(names(D4_INPUT), function(nm) {
  p <- D4_INPUT[[nm]]
  observed <- d4_sha256(p)
  expected <- if (nm %chin% names(D4_EXPECTED_SHA256)) {
    D4_EXPECTED_SHA256[[nm]]
  } else {
    NA_character_
  }
  data.table(
    input_name = nm,
    path = normalizePath(p, winslash = "/", mustWork = TRUE),
    bytes = file.info(p)$size,
    sha256 = observed,
    expected_sha256 = expected,
    hash_match = if (is.na(expected)) NA else observed == expected
  )
}))
locked_manifest <- manifest[!is.na(expected_sha256)]
if (any(!locked_manifest$hash_match)) {
  stop("Locked D4 input hash mismatch.", call. = FALSE)
}
d4_atomic_fwrite(
  manifest,
  file.path(D4_PROVENANCE_DIR, "D4_input_sha256_manifest.csv")
)

labels <- as.data.table(readRDS(D4_INPUT[["labels"]]))
final_needed <- c(
  "stay_id", "gender", "age", "sofa_score", "mortality_30d",
  "creatinine_max", "bun_max", "potassium_max", "bicarbonate_min",
  "ph_min", "lactate_max", "urine_output_24h_ml", "mbp_min",
  "heart_rate_max", "resp_rate_max", "gcs_min"
)
final <- fread(D4_INPUT[["final_full"]], select = final_needed)
baseline <- fread(D4_INPUT[["baseline"]])
rrt_first <- fread(D4_INPUT[["rrt_first"]])
times <- fread(D4_INPUT[["rrt_base_times"]])

d4_require_columns(labels, c("stay_id", "cluster_k2"), "locked labels")
d4_require_columns(final, final_needed, "final_full")
d4_require_columns(
  baseline, c("stay_id", "sapsii", "aki_stage_0_24h"), "baseline_covars"
)
d4_require_columns(
  times,
  c(
    "stay_id", "subject_id", "hadm_id", "intime", "outtime",
    "first_rrt_hours", "death_time_hours", "mortality_30d_from_time"
  ),
  "rrt_hte_base_times"
)
d4_require_columns(rrt_first, c("stay_id", "first_rrt_hours"), "rrt_first")

objects <- list(
  labels = labels, final = final, baseline = baseline,
  rrt_first = rrt_first, times = times
)
duplicate_audit <- rbindlist(lapply(names(objects), function(nm) {
  x <- objects[[nm]]
  data.table(
    object = nm,
    rows = nrow(x),
    unique_stays = uniqueN(x$stay_id),
    duplicate_stays = sum(duplicated(x$stay_id))
  )
}))
if (any(duplicate_audit$duplicate_stays > 0L)) {
  stop("Duplicated stay_id in a D4 source object.", call. = FALSE)
}

setnames(rrt_first, "first_rrt_hours", "first_rrt_hours_rrt_file")
setnames(times, "first_rrt_hours", "first_rrt_hours_time_file")
dt <- merge(labels, final, by = "stay_id", all = FALSE, sort = FALSE)
dt <- merge(
  dt,
  baseline[, .(stay_id, sapsii, aki_stage_0_24h)],
  by = "stay_id", all.x = TRUE, sort = FALSE
)
dt <- merge(
  dt,
  times[, .(
    stay_id, subject_id, hadm_id, intime, outtime,
    first_rrt_hours_time_file, death_time_hours,
    mortality_30d_from_time
  )],
  by = "stay_id", all.x = TRUE, sort = FALSE
)
dt <- merge(
  dt,
  rrt_first[, .(stay_id, first_rrt_hours_rrt_file)],
  by = "stay_id", all.x = TRUE, sort = FALSE
)

if (nrow(dt) != 20049L || uniqueN(dt$stay_id) != 20049L) {
  stop("Locked D4 merge is not one row per 20,049 stays.", call. = FALSE)
}
label_counts <- dt[, .N, by = cluster_k2][order(cluster_k2)]
if (!identical(label_counts$cluster_k2, c(1L, 2L)) ||
    !identical(label_counts$N, c(3992L, 16057L))) {
  stop("Locked K2 counts differ from 3,992/16,057.", call. = FALSE)
}

for (v in c(
  "first_rrt_hours_time_file", "first_rrt_hours_rrt_file",
  "death_time_hours"
)) {
  set(dt, j = v, value = suppressWarnings(as.numeric(dt[[v]])))
}
dt[, intime_parsed := as.POSIXct(intime, tz = "UTC")]
dt[, outtime_parsed := as.POSIXct(outtime, tz = "UTC")]
if (anyNA(dt$intime_parsed) || anyNA(dt$outtime_parsed)) {
  stop("Failed to parse an ICU intime/outtime value.", call. = FALSE)
}
dt[, outtime_hours := as.numeric(difftime(
  outtime_parsed, intime_parsed, units = "hours"
))]
if (any(!is.finite(dt$outtime_hours)) || any(dt$outtime_hours < 0)) {
  stop("Invalid ICU outtime relative to intime.", call. = FALSE)
}
dt[, first_rrt_hours := first_rrt_hours_time_file]
dt[, cluster_k2 := as.integer(cluster_k2)]
dt[, mortality_30d := as.integer(mortality_30d)]
dt[, aki_stage_0_24h := as.integer(aki_stage_0_24h)]

rrt_compare <- dt[!is.na(first_rrt_hours_rrt_file), .(
  n_rrt_rows = .N,
  missing_in_time_file = sum(is.na(first_rrt_hours_time_file)),
  max_abs_difference = max(
    abs(first_rrt_hours_rrt_file - first_rrt_hours_time_file),
    na.rm = TRUE
  )
)]
if (rrt_compare$missing_in_time_file != 0L ||
    rrt_compare$max_abs_difference > 1e-8) {
  stop("RRT timing files disagree.", call. = FALSE)
}

old0 <- copy(dt)
old1 <- old0[is.na(first_rrt_hours) | first_rrt_hours >= 0]
old2 <- old1[is.na(first_rrt_hours) | first_rrt_hours >= D4_LANDMARK_H]
old3 <- old2[is.na(death_time_hours) | death_time_hours >= D4_LANDMARK_H]
old4 <- old3[outtime_hours >= D4_LANDMARK_H]
old4_strict <- old3[outtime_hours > D4_LANDMARK_H]

reconciliation <- data.table(
  interface = c(
    "locked_cohort",
    "exclude_pre_ICU_RRT",
    "exclude_RRT_0_to_24h",
    "exclude_death_before_24h",
    "exclude_ICU_exit_before_24h",
    "strict_exit_after_24h_sensitivity"
  ),
  n_remaining = c(
    nrow(old0), nrow(old1), nrow(old2), nrow(old3),
    nrow(old4), nrow(old4_strict)
  ),
  removed_at_step = c(
    0L,
    nrow(old0) - nrow(old1),
    nrow(old1) - nrow(old2),
    nrow(old2) - nrow(old3),
    nrow(old3) - nrow(old4),
    nrow(old4) - nrow(old4_strict)
  )
)

source_audit <- rbindlist(list(
  duplicate_audit[, .(
    check = paste0("unique_stay_", object),
    observed = as.character(unique_stays),
    expected = ifelse(object %chin% c("labels", "final", "times"), "20049", "unique"),
    pass = duplicate_stays == 0L
  )],
  data.table(
    check = c(
      "locked_cohort_rows", "locked_C1", "locked_C2",
      "rrt_files_missing", "rrt_files_max_abs_diff",
      "timestamp_parse", "outtime_nonnegative"
    ),
    observed = c(
      nrow(dt), label_counts[cluster_k2 == 1L, N],
      label_counts[cluster_k2 == 2L, N],
      rrt_compare$missing_in_time_file,
      format(rrt_compare$max_abs_difference, scientific = TRUE),
      sum(!is.na(dt$intime_parsed) & !is.na(dt$outtime_parsed)),
      sum(dt$outtime_hours >= 0)
    ),
    expected = c("20049", "3992", "16057", "0", "<=1e-8", "20049", "20049"),
    pass = c(
      nrow(dt) == 20049L,
      label_counts[cluster_k2 == 1L, N] == 3992L,
      label_counts[cluster_k2 == 2L, N] == 16057L,
      rrt_compare$missing_in_time_file == 0L,
      rrt_compare$max_abs_difference <= 1e-8,
      sum(!is.na(dt$intime_parsed) & !is.na(dt$outtime_parsed)) == 20049L,
      sum(dt$outtime_hours >= 0) == 20049L
    )
  )
), fill = TRUE)
if (any(!source_audit$pass)) stop("D4 source audit failed.", call. = FALSE)

d4_atomic_fwrite(
  source_audit,
  file.path(D4_TABLE_DIR, "D4_source_and_hash_audit.csv")
)
d4_atomic_fwrite(
  reconciliation,
  file.path(D4_TABLE_DIR, "D4_old_denominator_reconciliation.csv")
)
d4_atomic_save_rds(
  dt,
  file.path(D4_PRIVATE_DIR, "D4_locked_source_private.rds")
)
d4_mark_stage(stage)
d4_msg("Completed ", stage)

