args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 01_snapshot_consistency.R eicu_external.csv L1_external_predictions_local.csv joined_holdout_records.csv")
}

read_input <- function(path, columns) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  x
}

e <- read_input(args[1], c("patientunitstayid", "temperature_min", "gcs_min"))
p <- read_input(args[2], c("patientunitstayid", "l1_k2", "outcome", "outcome_value"))
h <- read_input(args[3], c("patientunitstayid", "side", "l1_k2", "y"))
m <- p[p$outcome == "mortality", , drop = FALSE]

if (anyDuplicated(e$patientunitstayid) || anyDuplicated(m$patientunitstayid) ||
    anyDuplicated(h$patientunitstayid)) {
  stop("A stay identifier is duplicated within a cohort or mortality file")
}

matched <- match(m$patientunitstayid, h$patientunitstayid)
cat("eICU extract stays:", nrow(e), "\n")
cat("Temperature observed:", sum(!is.na(e$temperature_min)), "\n")
cat("GCS observed:", sum(!is.na(e$gcs_min)), "\n")
cat("L1 mortality-evaluable stays:", nrow(m), "\n")
cat("Held-out split stays:", nrow(h), "\n")
cat("Unmatched mortality stays:", sum(is.na(matched)), "\n")
cat("L1 label mismatches on matched stays:",
    sum(as.character(m$l1_k2) != as.character(h$l1_k2[matched]), na.rm = TRUE), "\n")
cat("Outcome mismatches on matched stays:",
    sum(as.character(m$outcome_value) != as.character(h$y[matched]), na.rm = TRUE), "\n")

if (anyNA(matched)) stop("The mortality files do not have the same stays")
counts <- aggregate(cbind(stays = rep(1L, nrow(h)),
                          deaths = as.integer(h$y == 1)),
                    by = list(side = h$side, l1_k2 = h$l1_k2), FUN = sum)
print(counts, row.names = FALSE)
