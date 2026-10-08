## =====================================================================
## Internal fidelity, sample-split generalization, and invariance audit.
## The regenerated generator is compared with the locked partition.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
require_namespace("data.table")
library(data.table)

set.seed(MASTER_SEED + 1L)
if (!file.exists(PATH_GENERATOR)) {
  stop("Frozen generator is missing. Run script 01 first.")
}
generator_file_md5_before <- hash_file(PATH_GENERATOR)
object <- readRDS(PATH_GENERATOR)
if (!inherits(object, "ehr_audit_k2_deployable_surrogate_generator")) {
  stop("Frozen object has the wrong class.")
}

locked <- as.data.table(readRDS(PATH_LOCKED_LABELS))
locked[, cluster_k2 := as.character(cluster_k2)]
M <- fread(
  PATH_FINAL_FULL,
  select = c(
    "stay_id", "subject_id", "hospital_mortality", "mortality_30d", "make30",
    FEATURES33
  ),
  showProgress = FALSE
)
if (nrow(M) != EXPECTED_MIMIC_N || anyDuplicated(M$stay_id)) {
  stop("MIMIC internal-validation denominator check failed.")
}
idx <- match(locked$stay_id, M$stay_id)
if (anyNA(idx)) stop("Locked IDs are not fully represented in final_full.")
M <- M[idx]
if (!identical(as.integer(M$stay_id), as.integer(locked$stay_id))) {
  stop("MIMIC-to-label row alignment failed.")
}

pred <- predict_ehr_audit_k2(object, M, id_col = "stay_id")
reference <- locked$cluster_k2
assigned <- as.character(pred$assigned_cluster)
apparent <- as.data.table(fidelity_metrics(reference, assigned))
apparent[, analysis := "apparent_full_development_cohort"]
setcolorder(apparent, c("analysis", setdiff(names(apparent), "analysis")))

confusion <- as.data.table(
  as.data.frame.matrix(table(
    locked = factor(reference, levels = SOURCE_LABEL_LEVELS),
    regenerated = factor(assigned, levels = SOURCE_LABEL_LEVELS)
  )),
  keep.rownames = "locked_cluster"
)
setnames(confusion, SOURCE_LABEL_LEVELS, paste0("assigned_", SOURCE_LABEL_LEVELS))

if (apparent$agreement_n != EXPECTED_APPARENT_AGREEMENT_N ||
    apparent$discordant_n != EXPECTED_APPARENT_DISCORDANT_N ||
    round(apparent$ari, 3) != EXPECTED_APPARENT_ARI_ROUNDED) {
  write_replace_csv(
    apparent, file.path(OUT_QC, "generator_v1_regression_discrepancy.csv")
  )
  stop(
    "Internal apparent fidelity differs from the locked regression reference."
  )
}

## ---------- Bootstrap uncertainty ----------
subject_id <- as.character(M$subject_id)
if (anyNA(subject_id)) stop("Missing subject_id in MIMIC validation cohort.")
rows_by_subject <- split(seq_len(nrow(M)), subject_id)
subjects <- names(rows_by_subject)
boot <- vector("list", BOOT_B)
set.seed(MASTER_SEED + 101L)
for (b in seq_len(BOOT_B)) {
  sampled <- sample(subjects, length(subjects), replace = TRUE)
  ii <- unlist(rows_by_subject[sampled], use.names = FALSE)
  fm <- fidelity_metrics(reference[ii], assigned[ii])
  boot[[b]] <- data.table(
    bootstrap_id = b,
    n = fm$n,
    agreement = fm$agreement,
    ari = fm$ari,
    prevalence_difference = fm$prevalence_difference
  )
}
boot <- rbindlist(boot)
boot_summary <- rbindlist(lapply(
  c("agreement", "ari", "prevalence_difference"),
  function(metric) {
    x <- boot[[metric]]
    ci <- percentile_interval(x)
    data.table(
      metric = metric,
      estimate = apparent[[metric]],
      bootstrap_mean = mean(x, na.rm = TRUE),
      bootstrap_sd = sd(x, na.rm = TRUE),
      bootstrap_mcse = sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))),
      ci_low = ci[1],
      ci_high = ci[2],
      bootstrap_n = sum(is.finite(x)),
      interval_method = "patient-level percentile bootstrap"
    )
  }
))

## ---------- Repeated stratified five-fold surrogate fidelity ----------
fold_rows <- list()
repeat_rows <- list()
row_counter <- 0L
for (r in seq_len(CV_REPEATS)) {
  set.seed(MASTER_SEED + 1000L + r)
  fold_id <- integer(nrow(M))
  for (k in SOURCE_LABEL_LEVELS) {
    ii <- which(reference == k)
    fold_id[ii] <- sample(rep_len(seq_len(CV_FOLDS), length(ii)))
  }
  repeat_assigned <- rep(NA_character_, nrow(M))
  for (f in seq_len(CV_FOLDS)) {
    train <- fold_id != f
    test <- fold_id == f
    fold_recipe <- derive_frozen_recipe(M[train], FEATURES33)
    train_prep <- apply_frozen_recipe(M[train], fold_recipe, FEATURES33)
    fold_centroids <- derive_locked_centroids(
      train_prep$X, reference[train], SOURCE_LABEL_LEVELS
    )
    test_prep <- apply_frozen_recipe(M[test], fold_recipe, FEATURES33)
    fold_pred <- assign_nearest_centroid(test_prep$X, fold_centroids)
    fold_assigned <- as.character(fold_pred$assigned_cluster)
    repeat_assigned[test] <- fold_assigned
    fm <- as.data.table(fidelity_metrics(reference[test], fold_assigned))
    fm[, `:=`(repeat_id = r, fold_id = f)]
    row_counter <- row_counter + 1L
    fold_rows[[row_counter]] <- fm
  }
  if (anyNA(repeat_assigned)) stop("Cross-fit prediction is incomplete.")
  fm_repeat <- as.data.table(fidelity_metrics(reference, repeat_assigned))
  fm_repeat[, repeat_id := r]
  repeat_rows[[r]] <- fm_repeat
}
crossfit_fold <- rbindlist(fold_rows)
crossfit_repeat <- rbindlist(repeat_rows)
crossfit_summary <- rbindlist(lapply(
  c(
    "agreement", "ari", "recall_c1", "recall_c2", "precision_c1",
    "precision_c2", "prevalence_difference"
  ),
  function(metric) {
    x <- crossfit_repeat[[metric]]
    data.table(
      metric = metric,
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      mcse = sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))),
      empirical_p025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
      empirical_p975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
      repeat_n = sum(is.finite(x))
    )
  }
))

## ---------- Missingness, boundary, and discordance ----------
M[, locked_cluster := reference]
M[, regenerated_cluster := assigned]
M[, concordant := locked_cluster == regenerated_cluster]
M[, missing_feature_n := pred$missing_feature_n]
M[, clipped_feature_n := pred$clipped_feature_n]
M[, absolute_margin := pred$absolute_margin]
M[, signed_margin_c1 := pred$signed_margin_c1]
M[, missingness_stratum := cut(
  missing_feature_n,
  breaks = c(-Inf, 0, 2, 5, Inf),
  labels = c("0", "1-2", "3-5", "6+"),
  right = TRUE
)]

missingness_fidelity <- M[, {
  fm <- fidelity_metrics(locked_cluster, regenerated_cluster)
  as.list(fm[1, ])
}, by = missingness_stratum]

summary_variables <- c(
  "hospital_mortality", "mortality_30d", "make30", "missing_feature_n",
  "clipped_feature_n", "absolute_margin", "creatinine_max", "bun_max",
  "bicarbonate_min", "aniongap_max", "lactate_max"
)
discordance_summary <- M[, c(
  list(n = .N),
  lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
  setNames(
    lapply(.SD, function(x) stats::median(as.numeric(x), na.rm = TRUE)),
    paste0(names(.SD), "_median")
  )
), by = .(concordant, locked_cluster, regenerated_cluster),
.SDcols = summary_variables]

private_discordant <- M[
  concordant == FALSE,
  c(
    "stay_id", "subject_id", "locked_cluster", "regenerated_cluster",
    "missing_feature_n", "clipped_feature_n", "absolute_margin",
    "signed_margin_c1", "hospital_mortality", "mortality_30d", "make30"
  ),
  with = FALSE
]

outcome_profiles <- rbindlist(list(
  M[, .(
    n = .N,
    hospital_mortality = mean(hospital_mortality, na.rm = TRUE),
    mortality_30d = mean(mortality_30d, na.rm = TRUE),
    make30 = mean(make30, na.rm = TRUE)
  ), by = .(cluster = locked_cluster)][, label_object := "locked_partition"],
  M[, .(
    n = .N,
    hospital_mortality = mean(hospital_mortality, na.rm = TRUE),
    mortality_30d = mean(mortality_30d, na.rm = TRUE),
    make30 = mean(make30, na.rm = TRUE)
  ), by = .(cluster = regenerated_cluster)][
    , label_object := "regenerated_generator"
  ]
), use.names = TRUE)
setcolorder(outcome_profiles, c("label_object", "cluster", "n",
                               "hospital_mortality", "mortality_30d", "make30"))

## ---------- Deterministic and batch-invariance tests ----------
set.seed(MASTER_SEED + 200L)
sample_idx <- sample(seq_len(nrow(M)), 1000L)
sample_data <- M[sample_idx]
sample_base <- predict_ehr_audit_k2(object, sample_data)
sample_repeat <- predict_ehr_audit_k2(object, sample_data)

set.seed(MASTER_SEED + 201L)
perm <- sample(seq_len(nrow(sample_data)))
sample_shuffled <- predict_ehr_audit_k2(object, sample_data[perm])
sample_shuffled <- sample_shuffled[order(perm), ]

single_labels <- vapply(seq_len(nrow(sample_data)), function(i) {
  as.character(predict_ehr_audit_k2(
    object, sample_data[i], return_diagnostics = FALSE
  )$assigned_cluster)
}, character(1))

split_points <- split(
  seq_len(nrow(sample_data)),
  rep(seq_len(7L), length.out = nrow(sample_data))
)
batch_pred <- vector("list", length(split_points))
for (i in seq_along(split_points)) {
  ii <- split_points[[i]]
  z <- predict_ehr_audit_k2(object, sample_data[ii])
  z$original_row <- ii
  batch_pred[[i]] <- z
}
batch_pred <- rbindlist(batch_pred)
setorder(batch_pred, original_row)

reloaded <- readRDS(PATH_GENERATOR)
sample_reloaded <- predict_ehr_audit_k2(reloaded, sample_data)

eq_numeric <- function(a, b, tol = 1e-12) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
}
invariance <- data.table(
  test = c(
    "repeated_call_labels", "repeated_call_distances",
    "row_order_labels", "row_order_distances",
    "single_patient_vs_batch_labels", "arbitrary_batch_split_labels",
    "arbitrary_batch_split_distances", "serialization_reload_labels",
    "serialization_reload_distances", "generator_file_hash_unchanged"
  ),
  pass = c(
    identical(
      as.character(sample_base$assigned_cluster),
      as.character(sample_repeat$assigned_cluster)
    ),
    eq_numeric(sample_base$distance_c1, sample_repeat$distance_c1) &&
      eq_numeric(sample_base$distance_c2, sample_repeat$distance_c2),
    identical(
      as.character(sample_base$assigned_cluster),
      as.character(sample_shuffled$assigned_cluster)
    ),
    eq_numeric(sample_base$distance_c1, sample_shuffled$distance_c1) &&
      eq_numeric(sample_base$distance_c2, sample_shuffled$distance_c2),
    identical(as.character(sample_base$assigned_cluster), single_labels),
    identical(
      as.character(sample_base$assigned_cluster),
      as.character(batch_pred$assigned_cluster)
    ),
    eq_numeric(sample_base$distance_c1, batch_pred$distance_c1) &&
      eq_numeric(sample_base$distance_c2, batch_pred$distance_c2),
    identical(
      as.character(sample_base$assigned_cluster),
      as.character(sample_reloaded$assigned_cluster)
    ),
    eq_numeric(sample_base$distance_c1, sample_reloaded$distance_c1) &&
      eq_numeric(sample_base$distance_c2, sample_reloaded$distance_c2),
    identical(generator_file_md5_before, hash_file(PATH_GENERATOR))
  )
)
if (any(!invariance$pass)) {
  write_replace_csv(
    invariance, file.path(OUT_QC, "generator_v1_invariance_tests.csv")
  )
  stop("At least one generator invariance test failed.")
}

denominator_audit <- data.table(
  item = c(
    "mimic_rows", "unique_stay_id", "unique_subject_id",
    "locked_labels", "apparent_agreement_n", "apparent_discordant_n",
    "bootstrap_replicates", "crossfit_repeats", "crossfit_folds"
  ),
  value = c(
    nrow(M), uniqueN(M$stay_id), uniqueN(M$subject_id), nrow(locked),
    apparent$agreement_n, apparent$discordant_n, nrow(boot),
    CV_REPEATS, CV_FOLDS
  )
)

write_replace_csv(
  apparent, file.path(OUT_TABLES, "internal_fidelity_summary.csv")
)
write_replace_csv(
  confusion, file.path(OUT_TABLES, "internal_fidelity_confusion_matrix.csv")
)
write_replace_csv(
  boot, file.path(OUT_TABLES, "internal_fidelity_bootstrap.csv")
)
write_replace_csv(
  boot_summary, file.path(OUT_TABLES, "internal_fidelity_bootstrap_summary.csv")
)
write_replace_csv(
  crossfit_fold,
  file.path(OUT_TABLES, "internal_fidelity_crossfit_by_fold.csv")
)
write_replace_csv(
  crossfit_repeat,
  file.path(OUT_TABLES, "internal_fidelity_crossfit_by_repeat.csv")
)
write_replace_csv(
  crossfit_summary,
  file.path(OUT_TABLES, "internal_fidelity_crossfit_summary.csv")
)
write_replace_csv(
  missingness_fidelity,
  file.path(OUT_TABLES, "internal_fidelity_by_missingness.csv")
)
write_replace_csv(
  discordance_summary,
  file.path(OUT_TABLES, "internal_discordant_patient_summary.csv")
)
write_replace_csv(
  outcome_profiles,
  file.path(OUT_TABLES, "internal_locked_vs_regenerated_outcomes.csv")
)
write_replace_csv(
  private_discordant,
  file.path(OUT_PRIVATE, "internal_discordant_records_local.csv")
)
write_replace_csv(
  invariance, file.path(OUT_QC, "generator_v1_invariance_tests.csv")
)
write_replace_csv(
  denominator_audit,
  file.path(OUT_QC, "generator_v1_internal_denominator_audit.csv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_LOGS, "generator_v1_internal_validation_sessionInfo.txt"),
  useBytes = TRUE
)

cat("Internal fidelity validation completed.\n")
print(apparent)
cat("\nBootstrap summary:\n")
print(boot_summary)
cat("\nRepeated five-fold summary:\n")
print(crossfit_summary)
cat("\nAll invariance tests passed.\n")
