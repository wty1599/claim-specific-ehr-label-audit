## ============================================================




## ============================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
source(file.path(PROJ_ROOT, "R", "common", "00b_cluster_utils.R"))
need_pkg("data.table"); library(data.table)

USE_MICE   <- FALSE
MICE_M     <- 5
MICE_MAXIT <- 5
set.seed(20250101)

## ---------- STEP 1: three feature sets from the clustering specification ----------
dt  <- fread(PATH_FINAL_FULL)
dec <- fread(file.path(DIR_QC, "qc_variable_decision.csv"))

high_miss      <- dec[suggest == "consider_exclusion_high_missingness", variable]
vars_inclusive <- setdiff(dec$variable, high_miss)               # S2b: 42

DROP_REDUNDANT <- c(
  "pt_max", "creatinine_min", "ast_max", "wbc_min", "abs_neutrophils_max",
  "po2_min", "baseexcess_min", "sbp_min", "mbp_mean")
vars_primary <- setdiff(vars_inclusive, DROP_REDUNDANT)          # S1: 33

vars_restricted <- c(
  "pao2fio2ratio_min", "platelets_min", "inr_max", "bilirubin_total_max",
  "mbp_min", "lactate_max", "heart_rate_max", "gcs_min",
  "creatinine_max", "bun_max", "urine_output_24h_ml",
  "wbc_max", "hemoglobin_min", "bicarbonate_min",
  "sodium_min", "potassium_max", "resp_rate_max", "temperature_max")  # S2a: 18

feature_sets <- list(primary = vars_primary,
                     sens_restricted = vars_restricted,
                     sens_inclusive  = vars_inclusive)
for (nm in names(feature_sets)) cat(sprintf("%-16s : %d variables\n", nm, length(feature_sets[[nm]])))
stopifnot(length(vars_primary) == 33, length(vars_restricted) == 18,
          length(vars_inclusive) == 42, all(vars_restricted %in% vars_inclusive))


Xf <- as.data.frame(dt[, vars_inclusive, with = FALSE])
cl <- apply_phys_ranges(Xf, vars_inclusive)
Xf <- cl$df
fwrite(cl$report, file.path(DIR_OUTPUT, "range_clean_report.csv"))
log_msg("Features with most values set to NA by physiological bounds:"); print(utils::head(cl$report, 8))


wb <- lapply(Xf, winsor_fit)
for (v in vars_inclusive) Xf[[v]] <- winsor_apply(Xf[[v]], wb[[v]])


if (USE_MICE) {
  need_pkg("mice")
  ckpt <- file.path(CKPT_DIR, "mice_imp.rds"); imp <- load_ckpt(ckpt)
  if (is.null(imp)) {
    log_msg("MICE m=", MICE_M, " maxit=", MICE_MAXIT, " (printFlag reports progress)")
    imp <- mice::mice(Xf, m = MICE_M, method = "pmm", maxit = MICE_MAXIT,
                      seed = 123, printFlag = TRUE)
    save_ckpt(imp, ckpt)
  } else log_msg("Loading MICE checkpoint")
  Xf_done <- mice::complete(imp, 1)
  for (v in vars_inclusive) Xf_done[[v]] <- winsor_apply(Xf_done[[v]], wb[[v]])
} else {
  Xf_done <- Xf
  for (v in vars_inclusive) { na <- is.na(Xf_done[[v]]); if (any(na)) Xf_done[[v]][na] <- median(Xf_done[[v]], na.rm = TRUE) }
}


Xs <- lapply(feature_sets, function(v) scale(as.matrix(Xf_done[, v, drop = FALSE])))
stopifnot(!any(vapply(Xs, anyNA, logical(1))))


saveRDS(Xs,           file.path(DIR_MODEL, "Xs_sets.rds"))
saveRDS(feature_sets, file.path(DIR_MODEL, "feature_sets.rds"))

OUTCOME_VARS <- c("make30","rrt_30d","persistent_rd_30d","aki_kdigo_7d","aki_stage_kdigo_7d",
                  "rrt_7d","icu_mortality","hospital_mortality",
                  "mortality_28d","mortality_30d","mortality_90d","icu_los_days","hospital_los_days")
id_col <- intersect(c("stay_id","hadm_id","subject_id"), names(dt))[1]
out <- as.data.frame(dt[, intersect(c(id_col, OUTCOME_VARS), names(dt)), with = FALSE])
saveRDS(list(id_col = id_col, outcomes = out), file.path(DIR_MODEL, "outcomes_local.rds"))

log_msg("Preprocessing complete. Samples=", nrow(Xs$primary), "; feature-set dimensions=",
        paste(sapply(Xs, ncol), collapse = "/"), "; imputation=", ifelse(USE_MICE, "MICE(1)", "median"))
