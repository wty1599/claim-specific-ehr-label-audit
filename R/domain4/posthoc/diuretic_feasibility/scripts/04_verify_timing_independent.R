options(stringsAsFactors = FALSE)
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", args), winslash = "/")))
source_file <- file.path(root, "private_not_for_release", "loop_source_by_stay_private.psv")
project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
risk_file <- file.path(project, "d4_landmark_support_interface_20260729",
  "output_formal_v1", "private_not_for_release", "D4_landmark_patient_ledger_private.rds")
summary_file <- file.path(root, "outputs", "timing_and_denominator_feasibility.csv")

columns <- c("stay_id", "intime", "outtime", "first_start", "first_window_start",
             "rows", "pre24_rows", "preicu_rows", "window_rows", "chf", "chf_rows")
d <- read.table(source_file, sep = "|", header = FALSE, col.names = columns,
                na.strings = "", colClasses = "character", quote = "",
                comment.char = "", fill = FALSE)
risk <- readRDS(risk_file)
match_index <- match(as.character(risk$stay_id), d$stay_id)
stopifnot(nrow(d) == 20049L, nrow(risk) == 17525L,
          !anyNA(match_index), !anyDuplicated(match_index))
d <- d[match_index, , drop = FALSE]
parse_time <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
start_h <- as.numeric(difftime(parse_time(d$first_start), parse_time(d$intime),
                               units = "hours"))
rrt_h <- risk$first_rrt_hours
rrt_h[is.na(rrt_h)] <- Inf
death_h <- risk$death_time_hours
death_h[is.na(death_h)] <- Inf
tol <- 1e-8
eligible_start <- !is.na(start_h) & start_h < 72 &
  start_h < rrt_h - tol & start_h <= death_h + tol &
  start_h <= risk$outtime_hours + tol
early <- eligible_start & start_h < 24
late <- eligible_start & start_h >= 24
reported <- read.csv(summary_file, stringsAsFactors = FALSE)

for (g in c("Overall", "C1", "C2")) {
  selected <- if (g == "Overall") rep(TRUE, nrow(risk)) else
    risk$cluster_k2 == as.integer(sub("C", "", g))
  row <- reported[reported$group == g, , drop = FALSE]
  stopifnot(nrow(row) == 1L,
            sum(selected & early) == row$first_record_before_24_n,
            sum(selected & late) == row$first_record_24_72_n,
            sum(selected & early) + sum(selected & late) ==
              row$first_record_before_72_n,
            sum(selected & as.integer(d$pre24_rows) == 0L) ==
              row$option_A_eligible_n)
}
stopifnot(reported$first_record_before_24_fraction[reported$group == "Overall"] > 0.5)
write.csv(data.frame(check = "Independent base-R timing and denominator replication",
                     pass = TRUE, details = "Overall, C1, and C2 counts match"),
          file.path(root, "outputs", "independent_timing_verification.csv"),
          row.names = FALSE)
cat("INDEPENDENT_TIMING_PASS\n")
