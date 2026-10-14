options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
out <- file.path(root, "outputs")
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
paths <- c(
  risk = file.path(project, "d4_landmark_support_interface_20260729/output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds"),
  final = file.path(project, "data/final_full.csv"),
  parallel_sofa_export = file.path(project, "data/SENECA_MIMIC_2025_ICU_24h_common_base_cohort.csv"),
  ablation_labels = file.path(root, "private_not_for_release/A4_first_seed_labels_private.rds")
)
stopifnot(all(file.exists(paths)))
before <- vapply(paths, sha, character(1))
stopifnot(before[["risk"]] == "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
          before[["final"]] == "CE0B38D2D3D6263A5C4615A6AEA1E78B19840BD4B35046DBFA65EE6BFBBA633B",
          before[["parallel_sofa_export"]] == "05E92E69E60A060EA3A23E72E4184DFE9C4A25953A0579E8C178DE02223AE0F8")

risk <- as.data.table(readRDS(paths[["risk"]]))[, .(stay_id, cluster_k2,
  sofa_score, intime, rrt_first_24_72)]
final <- fread(paths[["final"]], select = c("stay_id", "sepsis_onset_time"))
parallel <- fread(paths[["parallel_sofa_export"]],
                  select = c("stay_id", "icu_intime", "sofa_time", "sofa_score",
                             "suspected_infection_time"))
labels <- as.data.table(readRDS(paths[["ablation_labels"]]))[, .(stay_id, ablation_group)]
stopifnot(uniqueN(risk$stay_id) == 17525L, uniqueN(final$stay_id) == nrow(final),
          uniqueN(parallel$stay_id) == nrow(parallel),
          uniqueN(labels$stay_id) == nrow(labels))
joined <- merge(risk, final, by = "stay_id", all.x = TRUE, sort = FALSE)
joined <- merge(joined, parallel, by = "stay_id", all.x = TRUE, sort = FALSE)
joined <- merge(joined, labels, by = "stay_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(joined) == 17525L, !anyNA(joined$sofa_time),
          !anyNA(joined$ablation_group),
          all(joined$sofa_score.x == joined$sofa_score.y),
          all(joined$sepsis_onset_time == joined$suspected_infection_time),
          all(joined$intime == joined$icu_intime))

joined[, sofa_hours_after_icu := as.numeric(difftime(
  as.POSIXct(sofa_time, tz = "UTC"),
  as.POSIXct(icu_intime, tz = "UTC"), units = "hours"))]
stopifnot(all(is.finite(joined$sofa_hours_after_icu)))
joined[, sofa_after_landmark := sofa_hours_after_icu > 24]
summary <- rbindlist(lapply(c("K2", "A4"), function(partition) {
  group_var <- if (partition == "K2") "cluster_k2" else "ablation_group"
  joined[, .(
    n = .N,
    post_landmark_sofa_n = sum(sofa_after_landmark),
    post_landmark_sofa_record_positive_n = sum(sofa_after_landmark & rrt_first_24_72 == 1L),
    max_sofa_hours_after_icu = max(sofa_hours_after_icu)
  ), by = .(group = get(group_var))][, partition := partition][]
}), use.names = TRUE)
setcolorder(summary, c("partition", "group", "n", "post_landmark_sofa_n",
                       "post_landmark_sofa_record_positive_n",
                       "max_sofa_hours_after_icu"))
stopifnot(summary[partition == "K2", sum(post_landmark_sofa_n)] == 151L,
          summary[partition == "A4", sum(post_landmark_sofa_n)] == 151L,
          summary[, sum(post_landmark_sofa_record_positive_n)] == 0L)

after <- vapply(paths, sha, character(1))
stopifnot(identical(before, after))
fwrite(summary, file.path(out, "sofa_time_boundary_aggregate.csv"))
fwrite(data.table(role = names(paths), path = unname(paths), sha256 = unname(before),
                  unchanged = TRUE),
       file.path(out, "sofa_time_audit_source_manifest.csv"))
cat("SOFA time audit: 151/17525 reconstructed parallel SOFA times exceed hour 24; all source joins exact.\n")
