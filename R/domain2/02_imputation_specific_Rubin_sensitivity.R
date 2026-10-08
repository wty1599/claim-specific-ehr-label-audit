## =====================================================================
## 23b_mice_b2_rubin_only.R
## Standalone B2 analysis under MICE: Rubin-pooled decisive contrast
##   Question: Does the K=2 label add prognostic information beyond the
##             same 33 first-24h physiological features modelled directly?
##
## This script intentionally does NOT save X_primary33_std_mice.rds or labels.
## It only runs the B2 benchmark, avoiding the alignment check in
## 23_mice_primary_pooling.R.
##
## Inputs:
##   - 00_paths.R
##   - output/model/mice_primary.rds  (preferred; created by 23_mice_primary_pooling.R)
##   - final_full and baseline_covars.csv
##
## Outputs:
##   - output/qc/mice_b2_rubin_only.csv
##   - output/qc/mice_b2_per_imputation.csv
##   - output/qc/mice_b2_bootvar_byimp.csv
## =====================================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({
  need <- c("data.table", "mice", "glmnet")
  for (p in need) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(data.table)
  library(mice)
  library(glmnet)
})

set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- user-adjustable settings ----------------
MICE_RDS   <- file.path(DIR_MODEL, "mice_primary.rds")
MICE_M_USE <- 5L          # use first 5 imputations from mice_primary.rds
K_FOLDS    <- 10L         # outer CV folds for AUC estimation
INNER_FOLDS <- 5L         # inner CV folds for glmnet lambda selection
ENET_ALPHA <- 0.5
BOOT_B     <- 500L        # bootstrap variance for Rubin pooling CI
NSTART_KM  <- 25L

## ---------------- authoritative feature list ----------------
MANUAL_FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max",
  "ptt_max", "inr_max", "temperature_min", "temperature_max",
  "urine_output_24h_ml", "glucose_max", "aniongap_max", "potassium_min",
  "potassium_max", "hemoglobin_min", "sodium_min", "sodium_max", "wbc_max",
  "platelets_min", "bicarbonate_min", "chloride_min", "chloride_max", "bun_max",
  "creatinine_max", "resp_rate_max", "gcs_min", "spo2_min", "heart_rate_max",
  "mbp_min")

FEATURES33 <- MANUAL_FEATURES33
fset_path <- file.path(DIR_MODEL, "feature_sets.rds")
if (file.exists(fset_path)) {
  fs <- readRDS(fset_path)
  if (!is.null(fs$primary) && length(fs$primary) == 33L) {
    FEATURES33 <- fs$primary
    cat("Using authoritative feature list from feature_sets.rds$primary\n")
  } else {
    warning("feature_sets.rds found, but $primary is missing or not length 33; using manual FEATURES33.")
  }
} else {
  warning("feature_sets.rds not found; using manual FEATURES33.")
}

BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")

## ---------------- helper functions ----------------
winsor <- function(v, p = c(.01, .99)) {
  q <- quantile(v, p, na.rm = TRUE)
  v[v < q[1]] <- q[1]
  v[v > q[2]] <- q[2]
  v
}

parse_sex <- function(g) {
  if (is.numeric(g)) return(as.integer(g == 1))
  as.integer(toupper(substr(as.character(g), 1, 1)) == "M")
}

fauc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_stratified_folds <- function(y, k = 10L, seed = 1L) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cl in c(0, 1)) {
    ix <- which(y == cl)
    fold[ix] <- sample(rep_len(seq_len(k), length(ix)))
  }
  fold
}

KM2_label_highrisk <- function(Xc, y30) {
  km <- kmeans(Xc, centers = 2, nstart = NSTART_KM, iter.max = 100, algorithm = "Lloyd")
  cl <- km$cluster
  hi <- which.max(tapply(y30, cl, mean))
  as.integer(ifelse(cl == hi, 1L, 2L))  # 1 = higher 30-day mortality risk
}

fit_cv_glmnet_predict <- function(D, y, fold, inner_folds = INNER_FOLDS, alpha = ENET_ALPHA) {
  p <- numeric(length(y))
  for (f in sort(unique(fold))) {
    tr <- fold != f
    te <- fold == f
    set.seed(10000 + f)
    m <- cv.glmnet(
      x = D[tr, , drop = FALSE],
      y = y[tr],
      family = "binomial",
      alpha = alpha,
      nfolds = inner_folds,
      type.measure = "deviance"
    )
    p[te] <- as.numeric(predict(m, D[te, , drop = FALSE], s = "lambda.min", type = "response"))
  }
  p
}

## ---------------- reconstruct the analytic frame ----------------
if (!file.exists(MICE_RDS)) stop("MICE object not found: ", MICE_RDS,
                                 "\nRun 23_mice_primary_pooling.R first.")

imp <- readRDS(MICE_RDS)
if (imp$m < MICE_M_USE) stop("mice_primary.rds contains only ", imp$m, " imputations, but MICE_M_USE=", MICE_M_USE)

cat(sprintf("EHR audit project loaded. Using MICE object: %s\n", MICE_RDS))
cat(sprintf("Using first %d imputations; outer CV=%d; inner glmnet CV=%d; BOOT_B=%d\n",
            MICE_M_USE, K_FOLDS, INNER_FOLDS, BOOT_B))

D0 <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id"), names(D0))[1]
if (is.na(mid)) stop("No stay_id or hadm_id found in final_full.")

bc <- fread(file.path(DIR_DATA, "baseline_covars.csv"))
need <- c(mid, "age", "gender", "sofa_score", "mortality_30d", "make30", FEATURES33)
miss <- setdiff(need, names(D0))
if (length(miss)) stop("final_full missing: ", paste(miss, collapse = ", "))
if (!all(c(mid, "aki_stage_0_24h") %in% names(bc)))
  stop("baseline_covars.csv missing ", mid, " or aki_stage_0_24h")

RAW <- D0[, ..need]
RAW[, sex := parse_sex(gender)]
RAW <- merge(RAW, bc[, c(mid, "aki_stage_0_24h"), with = FALSE], by = mid, all.x = TRUE)
RAW <- RAW[complete.cases(RAW[, c("age", "sex", "sofa_score", "aki_stage_0_24h", "mortality_30d", "make30"), with = FALSE])]
for (j in FEATURES33) RAW[[j]] <- winsor(as.numeric(RAW[[j]]))

if (nrow(RAW) != nrow(complete(imp, 1))) {
  stop(sprintf("Row mismatch: RAW n=%d but completed MICE data n=%d. Recreate mice_primary.rds with the same filtering.",
               nrow(RAW), nrow(complete(imp, 1))))
}

## Optional observed-value alignment check against imp$data, not against imputed/scaled matrix.
common_obs_check <- intersect(FEATURES33, names(imp$data))
if (length(common_obs_check) > 0) {
  chk <- common_obs_check[1]
  obs <- !is.na(imp$data[[chk]]) & !is.na(RAW[[chk]])
  if (sum(obs) > 100) {
    diff_med <- median(abs(imp$data[[chk]][obs] - RAW[[chk]][obs]), na.rm = TRUE)
    cat(sprintf("Observed-value alignment check using %s: median absolute difference = %.6g\n", chk, diff_med))
  }
}

BASE <- as.matrix(RAW[, ..BASE_NUM])
storage.mode(BASE) <- "double"

## ---------------- B2 Rubin pooling ----------------
run_one_outcome <- function(outcome) {
  y <- as.integer(RAW[[outcome]])
  y30 <- as.integer(RAW$mortality_30d)
  fold <- make_stratified_folds(y, K_FOLDS, seed = ifelse(outcome == "mortality_30d", 1101L, 2202L))

  per_imp <- vector("list", MICE_M_USE)
  boot_var <- vector("list", MICE_M_USE)

  cat(sprintf("\n--- B2 outcome: %s ---\n", outcome))
  for (mm in seq_len(MICE_M_USE)) {
    cat(sprintf("  imputation %d/%d\n", mm, MICE_M_USE))
    comp <- complete(imp, mm)
    if (!all(FEATURES33 %in% names(comp))) stop("Completed imputation missing features.")
    Xc <- scale(as.matrix(comp[, FEATURES33]))
    lab <- factor(KM2_label_highrisk(Xc, y30))

    Dbase <- model.matrix(~ ., data = data.frame(BASE))[, -1, drop = FALSE]
    Dbp   <- model.matrix(~ ., data = data.frame(BASE, lab = lab))[, -1, drop = FALSE]
    Dfeat <- cbind(BASE, Xc)
    Dflab <- cbind(BASE, Xc, model.matrix(~ lab)[, -1, drop = FALSE])

    p_base <- fit_cv_glmnet_predict(Dbase, y, fold)
    p_bp   <- fit_cv_glmnet_predict(Dbp,   y, fold)
    p_feat <- fit_cv_glmnet_predict(Dfeat, y, fold)
    p_flab <- fit_cv_glmnet_predict(Dflab, y, fold)

    auc_base <- fauc(y, p_base)
    auc_bp   <- fauc(y, p_bp)
    auc_feat <- fauc(y, p_feat)
    auc_flab <- fauc(y, p_flab)
    d_decisive <- auc_flab - auc_feat

    ## bootstrap within imputation for variance of decisive dAUC
    set.seed(9000 + mm + ifelse(outcome == "mortality_30d", 0L, 100L))
    dboot <- replicate(BOOT_B, {
      ii <- sample.int(length(y), length(y), replace = TRUE)
      fauc(y[ii], p_flab[ii]) - fauc(y[ii], p_feat[ii])
    })
    vboot <- var(dboot, na.rm = TRUE)

    per_imp[[mm]] <- data.table(
      outcome = outcome,
      imputation = mm,
      auc_base = auc_base,
      auc_base_phen = auc_bp,
      auc_feat = auc_feat,
      auc_feat_lab = auc_flab,
      incr_base_to_basephen = auc_bp - auc_base,
      decisive_dAUC_feat_plus_label = d_decisive
    )
    boot_var[[mm]] <- data.table(
      outcome = outcome,
      imputation = mm,
      boot_var_decisive_dAUC = vboot,
      boot_sd_decisive_dAUC = sqrt(vboot)
    )
  }

  pi <- rbindlist(per_imp)
  bv <- rbindlist(boot_var)

  theta <- pi$decisive_dAUC_feat_plus_label
  W <- mean(bv$boot_var_decisive_dAUC, na.rm = TRUE)
  Bv <- var(theta, na.rm = TRUE)
  Tt <- W + (1 + 1 / MICE_M_USE) * Bv
  se <- sqrt(Tt)
  theta_bar <- mean(theta, na.rm = TRUE)

  pooled <- data.table(
    outcome = outcome,
    m = MICE_M_USE,
    auc_base = round(mean(pi$auc_base), 4),
    auc_base_phen = round(mean(pi$auc_base_phen), 4),
    auc_feat = round(mean(pi$auc_feat), 4),
    auc_feat_lab = round(mean(pi$auc_feat_lab), 4),
    incr_base_to_basephen = round(mean(pi$incr_base_to_basephen), 4),
    decisive_dAUC_feat_plus_label = round(theta_bar, 5),
    dAUC_lo = round(theta_bar - 1.96 * se, 5),
    dAUC_hi = round(theta_bar + 1.96 * se, 5),
    within_imp_var_W = signif(W, 4),
    between_imp_var_B = signif(Bv, 4),
    total_var_T = signif(Tt, 4)
  )

  list(per_imp = pi, boot_var = bv, pooled = pooled)
}

res_mort <- run_one_outcome("mortality_30d")
res_make <- run_one_outcome("make30")

per_imp_all <- rbind(res_mort$per_imp, res_make$per_imp)
bootvar_all <- rbind(res_mort$boot_var, res_make$boot_var)
pooled_all <- rbind(res_mort$pooled, res_make$pooled)

cat("\n===== B2 Rubin-pooled decisive contrast =====\n")
print(pooled_all)
cat("\n===== Per-imputation B2 estimates =====\n")
print(per_imp_all)

fwrite(pooled_all, file.path(OUT_DIR, "mice_b2_rubin_only.csv"))
fwrite(per_imp_all, file.path(OUT_DIR, "mice_b2_per_imputation.csv"))
fwrite(bootvar_all, file.path(OUT_DIR, "mice_b2_bootvar_byimp.csv"))

cat("\nDone. Wrote:\n")
cat(" - ", file.path(OUT_DIR, "mice_b2_rubin_only.csv"), "\n", sep = "")
cat(" - ", file.path(OUT_DIR, "mice_b2_per_imputation.csv"), "\n", sep = "")
cat(" - ", file.path(OUT_DIR, "mice_b2_bootvar_byimp.csv"), "\n", sep = "")

cat("\nInterpretation:\n")
cat("If decisive_dAUC_feat_plus_label is approximately 0 and its CI includes 0,\n")
cat("the K=2 label remains redundant beyond the same 33 physiological features under MICE.\n")
## =====================================================================
