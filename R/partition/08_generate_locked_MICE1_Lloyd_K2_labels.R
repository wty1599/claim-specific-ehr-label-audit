## =====================================================================
## 23c_save_mice_matrix_labels.R
## Save MICE primary standardized matrix + K=2 labels after mice_primary.rds exists.
## Purpose: fix the final saving/alignment step after 23_mice_primary_pooling.R
##          successfully completed ARI/discreteness/B2 but failed at rank-correlation check.
##
## Outputs:
##   output/model/X_primary33_std_mice.rds
##   output/model/labels_primary_mice.rds
##   output/qc/mice_primary_label_summary.csv
##   output/qc/mice_primary_feature_alignment_check.csv
## =====================================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({
  need <- c("data.table", "mice")
  for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(mice)
})

set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PRIMARY_IMP <- 1L

## Prefer authoritative primary feature list from feature_sets.rds
feature_rds <- file.path(DIR_MODEL, "feature_sets.rds")
if (file.exists(feature_rds)) {
  fs <- readRDS(feature_rds)
  if (!"primary" %in% names(fs)) stop("feature_sets.rds exists but has no $primary element")
  FEATURES33 <- fs$primary
  message("Using authoritative feature_sets.rds$primary, n = ", length(FEATURES33))
} else {
  FEATURES33 <- c(
    "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
    "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
    "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml", "glucose_max",
    "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min", "sodium_min",
    "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min", "chloride_min",
    "chloride_max", "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
    "spo2_min", "heart_rate_max", "mbp_min")
  warning("feature_sets.rds not found; using manual 33-feature list")
}
stopifnot(length(FEATURES33) == 33)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
winsor <- function(v, p = c(.01, .99)) {
  q <- quantile(v, p, na.rm = TRUE)
  v[v < q[1]] <- q[1]
  v[v > q[2]] <- q[2]
  v
}
KM <- function(m, k, ns = 25) {
  kmeans(m, k, nstart = ns, iter.max = 100, algorithm = "Lloyd")
}
parse_sex <- function(g) as.integer(toupper(substr(as.character(g), 1, 1)) == "M")

## ---------- recreate the same RAW frame used in 23_mice_primary_pooling ----------
d <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id"), names(d))[1]
if (is.na(mid)) stop("No stay_id or hadm_id found in final_full")
bc <- fread(file.path(DIR_DATA, "baseline_covars.csv"))
need <- c(mid, "age", "gender", "sofa_score", "mortality_30d", "make30", FEATURES33)
miss <- setdiff(need, names(d))
if (length(miss)) stop("PATH_FINAL_FULL missing columns: ", paste(miss, collapse = ", "))
core <- d[, ..need]
core[, sex := parse_sex(gender)]
if (!"aki_stage_0_24h" %in% names(bc)) stop("baseline_covars.csv missing aki_stage_0_24h")
core <- merge(core, bc[, c(mid, "aki_stage_0_24h"), with = FALSE], by = mid, all.x = TRUE)
core <- core[complete.cases(core[, c("age", "sex", "sofa_score", "aki_stage_0_24h", "mortality_30d", "make30"), with = FALSE])]
RAW <- copy(core)
for (j in FEATURES33) RAW[[j]] <- winsor(as.numeric(RAW[[j]]))
cat(sprintf("Recreated RAW frame: n=%d, features=%d, id=%s\n", nrow(RAW), length(FEATURES33), mid))

## ---------- load existing MICE object ----------
mice_path <- file.path(DIR_MODEL, "mice_primary.rds")
if (!file.exists(mice_path)) stop("Missing ", mice_path, ". Run 23_mice_primary_pooling first.")
imp <- readRDS(mice_path)
if (imp$m < PRIMARY_IMP) stop("mice_primary.rds has fewer imputations than PRIMARY_IMP")
completed <- as.data.table(complete(imp, PRIMARY_IMP))
if (!all(FEATURES33 %in% names(completed))) {
  stop("Completed MICE data missing features: ", paste(setdiff(FEATURES33, names(completed)), collapse = ", "))
}
if (nrow(completed) != nrow(RAW)) {
  stop("Row mismatch: completed MICE n=", nrow(completed), " vs RAW n=", nrow(RAW))
}

## ---------- robust alignment check: observed values should be unchanged after MICE ----------
check <- rbindlist(lapply(FEATURES33, function(j) {
  obs <- which(!is.na(RAW[[j]]))
  mx <- if (length(obs)) max(abs(as.numeric(completed[[j]][obs]) - as.numeric(RAW[[j]][obs])), na.rm = TRUE) else NA_real_
  data.table(feature = j, n_observed = length(obs), max_abs_diff_observed = mx)
}))
fwrite(check, file.path(OUT_DIR, "mice_primary_feature_alignment_check.csv"))
bad <- check[is.finite(max_abs_diff_observed) & max_abs_diff_observed > 1e-8]
if (nrow(bad)) {
  print(bad)
  stop("Alignment check failed: observed values changed after MICE extraction")
}
cat("Observed-value alignment check passed for all features.\n")

## ---------- standardize imputed matrix and cluster ----------
Xmat <- as.matrix(completed[, ..FEATURES33])
storage.mode(Xmat) <- "double"
Xp <- scale(Xmat)
colnames(Xp) <- FEATURES33
if (anyNA(Xp)) stop("Standardized MICE matrix contains NA")

km <- KM(Xp, 2)
cl <- km$cluster
y30 <- RAW$mortality_30d
hi <- which.max(tapply(y30, cl, mean))
lab <- ifelse(cl == hi, 1L, 2L)  # 1 = higher 30-day mortality risk

## ---------- save outputs ----------
Xsave <- as.data.table(as.data.frame(Xp))
Xsave[, stay_id := as.integer(RAW[[mid]])]
setcolorder(Xsave, c("stay_id", FEATURES33))
labels <- data.table(stay_id = as.integer(RAW[[mid]]), cluster_k2 = lab)

if (anyDuplicated(Xsave$stay_id)) warning("Duplicate stay_id values detected in saved matrix")
if (!identical(Xsave$stay_id, labels$stay_id)) stop("Xsave and labels stay_id mismatch")

saveRDS(Xsave, file.path(DIR_MODEL, "X_primary33_std_mice.rds"))
saveRDS(labels, file.path(DIR_MODEL, "labels_primary_mice.rds"))

summ <- labels[, .(n = .N), by = cluster_k2]
summ <- merge(summ, data.table(cluster_k2 = lab, mortality_30d = y30, make30 = RAW$make30)[,
  .(mortality_30d = mean(mortality_30d), make30 = mean(make30)), by = cluster_k2],
  by = "cluster_k2")
summ[, `:=`(mortality_30d_pct = round(100 * mortality_30d, 2), make30_pct = round(100 * make30, 2))]
fwrite(summ, file.path(OUT_DIR, "mice_primary_label_summary.csv"))

cat("\nSaved:\n")
cat(" - ", file.path(DIR_MODEL, "X_primary33_std_mice.rds"), "\n", sep = "")
cat(" - ", file.path(DIR_MODEL, "labels_primary_mice.rds"), "\n", sep = "")
cat(" - ", file.path(OUT_DIR, "mice_primary_label_summary.csv"), "\n", sep = "")
cat("\nLabel summary (cluster_k2=1 is higher 30-day mortality):\n")
print(summ)
cat("\nNext: rerun 20_positive_controls_mice.R and 21_external_paired_mice.R.\n")
