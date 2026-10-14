options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", args), winslash = "/")))
private <- file.path(root, "private_not_for_release")
out <- file.path(root, "outputs")
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
spec_sha <- "774AE3F8ED25DB137C1FDA4BB86B7EC99B165A7CE896D905BE909DEA88620F1B"
stopifnot(sha(file.path(root, "01_PREQUERY_SPEC.md")) == spec_sha)

project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
risk_path <- file.path(project, "d4_landmark_support_interface_20260729",
  "output_formal_v1", "private_not_for_release", "D4_landmark_patient_ledger_private.rds")
stopifnot(sha(risk_path) == "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7")
source_path <- file.path(private, "loop_source_by_stay_private.psv")
stopifnot(file.exists(source_path), file.info(source_path)$size > 0)

names_source <- c(
  "stay_id", "intime_sql", "outtime_sql", "first_start_sql",
  "first_start_24_72_sql", "qualifying_rows", "pre24_rows",
  "preicu_overlap_rows", "starts_24_72_rows", "chf_code",
  "charlson_rows"
)
src <- fread(source_path, sep = "|", header = FALSE, col.names = names_source,
             na.strings = c("", "NA"), colClasses = "character")
stopifnot(nrow(src) == 20049L, uniqueN(src$stay_id) == nrow(src),
          !anyNA(src$stay_id))
risk <- as.data.table(readRDS(risk_path))
stopifnot(nrow(risk) == 17525L, uniqueN(risk$stay_id) == nrow(risk),
          risk[cluster_k2 == 1L, .N] == 3151L,
          risk[cluster_k2 == 2L, .N] == 14374L,
          risk[cluster_k2 == 1L, sum(rrt_first_24_72)] == 265L,
          risk[cluster_k2 == 2L, sum(rrt_first_24_72)] == 49L)
idx <- match(risk$stay_id, as.integer(src$stay_id))
stopifnot(!anyNA(idx), !anyDuplicated(idx))
z <- copy(src[idx])
time_parse <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
z[, `:=`(
  intime_sql = time_parse(intime_sql),
  outtime_sql = time_parse(outtime_sql),
  first_start_sql = time_parse(first_start_sql),
  first_start_24_72_sql = time_parse(first_start_24_72_sql)
)]
for (v in c("qualifying_rows", "pre24_rows", "preicu_overlap_rows",
            "starts_24_72_rows", "chf_code", "charlson_rows")) {
  z[, (v) := as.integer(get(v))]
}
stopifnot(!anyNA(z$intime_sql), !anyNA(z$outtime_sql),
          all(z$outtime_sql >= z$intime_sql),
          all(z$qualifying_rows >= 0L & z$pre24_rows >= 0L),
          all(z$pre24_rows <= z$qualifying_rows),
          all(z$preicu_overlap_rows <= z$pre24_rows),
          all(z$starts_24_72_rows <= z$qualifying_rows),
          all(z$pre24_rows == 0L | !is.na(z$first_start_sql)),
          all(z$starts_24_72_rows == 0L | !is.na(z$first_start_24_72_sql)))

intime_r <- as.POSIXct(risk$intime_parsed, tz = "UTC")
outtime_r <- as.POSIXct(risk$outtime_parsed, tz = "UTC")
stopifnot(all(abs(as.numeric(difftime(z$intime_sql, intime_r, units = "secs"))) < 1),
          all(abs(as.numeric(difftime(z$outtime_sql, outtime_r, units = "secs"))) < 1))
z[, first_hours := as.numeric(difftime(first_start_sql, intime_sql, units = "hours"))]
z[, first_24_72_hours := as.numeric(difftime(first_start_24_72_sql, intime_sql,
                                             units = "hours"))]
stopifnot(all(is.na(z$first_24_72_hours) |
                (z$first_24_72_hours >= 24 & z$first_24_72_hours < 72)),
          all((z$pre24_rows > 0L) ==
                (!is.na(z$first_hours) & z$first_hours < 24)))
z[, `:=`(cluster_k2 = risk$cluster_k2,
         rrt_first_24_72 = risk$rrt_first_24_72,
         mortality_30d = risk$mortality_30d)]
tol <- 1e-8
rrt_h <- fifelse(is.na(risk$first_rrt_hours), Inf, risk$first_rrt_hours)
death_h <- fifelse(is.na(risk$death_time_hours), Inf, risk$death_time_hours)
exit_h <- risk$outtime_hours
z[, `:=`(
  first_opportunity_valid = !is.na(first_hours) & first_hours < 72 &
    first_hours < rrt_h - tol & first_hours <= death_h + tol &
    first_hours <= exit_h + tol,
  first_record_24_72_after_or_tied_with_RRT = !is.na(first_hours) &
    first_hours >= 24 & first_hours < 72 & first_hours >= rrt_h - tol,
  after24_opportunity_valid = !is.na(first_24_72_hours) &
    first_24_72_hours < rrt_h - tol &
    first_24_72_hours <= death_h + tol &
    first_24_72_hours <= exit_h + tol
)]

timing <- rbindlist(lapply(c(NA_integer_, 1L, 2L), function(group) {
  d <- if (is.na(group)) z else z[cluster_k2 == group]
  early <- sum(d$first_opportunity_valid & d$first_hours < 24)
  late <- sum(d$first_opportunity_valid & d$first_hours >= 24 & d$first_hours < 72)
  all_before_72 <- early + late
  data.table(
    group = if (is.na(group)) "Overall" else paste0("C", group),
    anchor_n = nrow(d),
    preicu_overlapping_record_n = sum(d$preicu_overlap_rows > 0L),
    first_record_24_72_after_or_tied_with_RRT_n = sum(
      d$first_record_24_72_after_or_tied_with_RRT
    ),
    first_record_before_72_but_competing_event_precluded_n = sum(
      !is.na(d$first_hours) & d$first_hours < 72 & !d$first_opportunity_valid
    ),
    first_record_before_24_n = early,
    first_record_24_72_n = late,
    first_record_before_72_n = all_before_72,
    first_record_before_24_fraction = if (all_before_72) early / all_before_72 else NA_real_,
    pre24_IV_loop_excluded_n = sum(d$pre24_rows > 0L),
    option_A_eligible_n = sum(d$pre24_rows == 0L),
    option_A_eligible_with_window_record_n = sum(
      d$pre24_rows == 0L & d$after24_opportunity_valid
    ),
    chf_code_positive_n = sum(d$chf_code == 1L, na.rm = TRUE),
    chf_code_unknown_n = sum(is.na(d$chf_code))
  )
}))
stopifnot(timing[group == "Overall", anchor_n] == 17525L,
          timing[group == "C1", anchor_n] + timing[group == "C2", anchor_n] == 17525L,
          timing[group == "Overall", option_A_eligible_n] ==
            sum(timing[group != "Overall", option_A_eligible_n]))

fwrite(timing, file.path(out, "timing_and_denominator_feasibility.csv"))
qc <- data.table(
  check = c("original_risk_set_reproduced", "source_cohort_unique",
            "all_risk_stays_joined", "ICU_times_match", "no_patient_file_in_outputs"),
  pass = rep(TRUE, 5L),
  observed = c("17525; C1 3151; C2 14374; RRT C1 265; C2 49",
               as.character(nrow(src)), as.character(nrow(z)),
               "within 1 second", "aggregate tables only")
)
fwrite(qc, file.path(out, "timing_source_QC.csv"))
fwrite(data.table(
  input = c("prequery_spec", "original_risk_ledger", "read_only_sql",
            "timing_analysis_script", "private_medication_extract"),
  sha256 = c(spec_sha, sha(risk_path),
             sha(file.path(root, "scripts/02_extract_loop_source_readonly.sql")),
             sha(file.path(root, "scripts/03_check_timing.R")),
             sha(source_path))
), file.path(out, "timing_input_hashes.csv"))

first <- timing[group == "Overall"]
if (first$first_record_before_72_n == 0L) {
  writeLines("STOP: no first qualifying loop administration by hour 72",
             file.path(out, "timing_decision.txt"))
  cat("STOP_NO_STARTS\n")
} else if (first$first_record_before_24_fraction > 0.5) {
  writeLines("STOP: more than half of first qualifying records by hour 72 began before hour 24; no support model may be fitted",
             file.path(out, "timing_decision.txt"))
  cat("STOP_MAJORITY_BEFORE_24\n")
} else {
  writeLines("PROCEED: no majority-before-24 stop; first-event/support analysis may be run under the locked specification",
             file.path(out, "timing_decision.txt"))
  cat("PROCEED_TO_SUPPORT\n")
}
