## =====================================================================
## 24_crossfit_k2_b2_fig_v4.R
## Strict cross-fitted K=2 label + B2 raw-feature benchmark + figure output
##
## Purpose:
##   Sensitivity analysis for Figure/Methods: address possible leakage from
##   full-cohort clustering. For every outer fold, preprocessing, MICE,
##   scaling, K=2 centroids, and all prediction models are learned from the
##   TRAINING fold only. The held-out fold is imputed via mice::ignore,
##   standardized with training parameters, assigned to K=2 by training
##   centroids, and predicted by training-fitted models.
##
## Models:
##   base      = age + sex + SOFA + AKI stage
##   base_phen = base + cross-fitted K=2 label
##   feat      = base + 33 first-24h physiological features, elastic-net logistic
##   feat_lab  = feat + cross-fitted K=2 label
##
## Key contrast:
##   feat_lab - feat. If ~0 with CI including 0, the K=2 label remains redundant
##   beyond the same raw physiology even under strict cross-fitting.
##
## Design choices in this figure script:
##   - Always runs both outcomes: mortality_30d and make30.
##   - No outcome-selection switch, to keep figure data complete.
##   - Fold-level checkpoints and resume support.
##   - Safe scaling for zero/near-zero SD variables.
##   - MICE warnings/loggedEvents captured to fold logs.
##   - Outputs summary CSV, prediction CSVs, and two figure files.
## =====================================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

suppressPackageStartupMessages({
  need <- c("data.table", "glmnet", "mice", "ggplot2")
  for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(glmnet)
  library(mice)
  library(ggplot2)
})

set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FIG_DIR <- file.path(DIR_OUTPUT, "fig")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- user controls ----------------
## Formal figure run: keep QUICK_RUN <- FALSE.
## Smoke test only: set QUICK_RUN <- TRUE. It still runs both outcomes, but with
## fewer folds/imputations/bootstraps; do not use QUICK_RUN output in manuscript.
QUICK_RUN <- FALSE

K_FOLDS    <- if (QUICK_RUN) 3L else 10L
MICE_M     <- if (QUICK_RUN) 2L else 5L
MICE_MAXIT <- if (QUICK_RUN) 3L else 10L
BOOT_B     <- if (QUICK_RUN) 100L else 500L
ENET_ALPHA <- 0.5
INNER_CV   <- 5L
NSTART_KM  <- 25L
RESUME_FROM_CHECKPOINT <- TRUE
OUTCOMES <- c("mortality_30d", "make30")

FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
  "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml", "glucose_max",
  "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min", "sodium_min",
  "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min", "chloride_min",
  "chloride_max", "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
PRED_COLS <- c("base", "base_phen", "feat", "feat_lab")

## ---------------- helpers ----------------
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

boot_ci_delta <- function(y, p_aug, p_ref, B = BOOT_B) {
  n <- length(y)
  d <- replicate(B, {
    i <- sample.int(n, replace = TRUE)
    fauc(y[i], p_aug[i]) - fauc(y[i], p_ref[i])
  })
  quantile(d, c(.025, .975), na.rm = TRUE)
}

safe_winq <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 5L) {
    med <- median(x, na.rm = TRUE)
    return(c(med, med))
  }
  as.numeric(quantile(x, c(.01, .99), na.rm = TRUE, names = FALSE))
}

safe_clip <- function(x, q) {
  x <- as.numeric(x)
  pmin(pmax(x, q[1]), q[2])
}

safe_scale_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr, na.rm = TRUE)
  sdv <- apply(Xtr, 2, sd, na.rm = TRUE)
  bad <- !is.finite(sdv) | is.na(sdv) | sdv < 1e-8
  if (any(bad)) {
    warning("Zero/near-zero SD in training fold for: ", paste(colnames(Xtr)[bad], collapse = ", "))
    sdv[bad] <- 1
  }
  Ztr <- scale(Xtr, center = mu, scale = sdv)
  Zte <- scale(Xte, center = mu, scale = sdv)
  Ztr[!is.finite(Ztr)] <- 0
  Zte[!is.finite(Zte)] <- 0
  list(Ztr = Ztr, Zte = Zte, mu = mu, sd = sdv, zero_sd = names(sdv)[bad])
}

make_stratified_folds <- function(y, K) {
  folds <- integer(length(y))
  for (cl in c(0, 1)) {
    ix <- which(y == cl)
    folds[ix] <- sample(rep_len(seq_len(K), length(ix)))
  }
  folds
}

d2_euclidean <- function(A, C) {
  aa <- rowSums(A^2)
  bb <- rowSums(C^2)
  outer(aa, bb, "+") - 2 * A %*% t(C)
}

read_or_make_folds <- function(y, yname) {
  fold_file <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_folds.csv", yname))
  if (RESUME_FROM_CHECKPOINT && file.exists(fold_file)) {
    ff <- fread(fold_file)
    if (nrow(ff) == length(y) && "fold" %in% names(ff) && max(ff$fold, na.rm = TRUE) == K_FOLDS) {
      cat(sprintf("  [%s] reusing existing fold assignment: %s\n", yname, fold_file))
      return(ff$fold)
    }
    cat(sprintf("  [%s] existing fold file incompatible with current K=%d; regenerating\n", yname, K_FOLDS))
  }
  folds <- make_stratified_folds(y, K_FOLDS)
  fwrite(data.table(row_id = seq_along(y), outcome = yname, y = y, fold = folds, K_folds = K_FOLDS), fold_file)
  folds
}

read_or_init_predictions <- function(yname, n) {
  pred_file <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_predictions.csv", yname))
  if (RESUME_FROM_CHECKPOINT && file.exists(pred_file)) {
    pp <- fread(pred_file)
    settings_ok <- all(c("K_folds", "MICE_M", "MICE_MAXIT") %in% names(pp)) &&
      all(pp$K_folds == K_FOLDS) && all(pp$MICE_M == MICE_M) && all(pp$MICE_MAXIT == MICE_MAXIT)
    if (nrow(pp) == n && all(PRED_COLS %in% names(pp)) && settings_ok) {
      cat(sprintf("  [%s] reusing prediction checkpoint: %s\n", yname, pred_file))
      P <- as.matrix(pp[, ..PRED_COLS])
      storage.mode(P) <- "double"
      return(P)
    }
    cat(sprintf("  [%s] existing prediction checkpoint incompatible with current settings; starting fresh\n", yname))
  }
  matrix(NA_real_, n, length(PRED_COLS), dimnames = list(NULL, PRED_COLS))
}

save_predictions <- function(yname, y, folds, P) {
  pred_file <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_predictions.csv", yname))
  out <- data.table(
    row_id = seq_along(y), outcome = yname, y = y, fold = folds,
    K_folds = K_FOLDS, MICE_M = MICE_M, MICE_MAXIT = MICE_MAXIT,
    QUICK_RUN = QUICK_RUN
  )
  out <- cbind(out, as.data.table(P))
  fwrite(out, pred_file)
}

append_fold_log <- function(row) {
  log_file <- file.path(OUT_DIR, "crossfit_k2_b2_fold_log.csv")
  fwrite(row, log_file, append = file.exists(log_file))
}

## Fold-wise imputation: test rows are included only as rows to impute, not for fitting imputation models.
impute_train_test <- function(Xtr, Xte, seed) {
  Xall <- rbind(Xtr, Xte)
  ignore_vec <- c(rep(FALSE, nrow(Xtr)), rep(TRUE, nrow(Xte)))
  pred <- make.predictorMatrix(Xall)
  diag(pred) <- 0
  warn_vec <- character(0)
  im <- withCallingHandlers(
    mice(
      Xall,
      m = MICE_M,
      maxit = MICE_MAXIT,
      method = "pmm",
      predictorMatrix = pred,
      ignore = ignore_vec,
      printFlag = FALSE,
      seed = seed
    ),
    warning = function(w) {
      warn_vec <<- c(warn_vec, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  completed <- lapply(seq_len(MICE_M), function(k) {
    full <- complete(im, k)
    list(
      tr = full[seq_len(nrow(Xtr)), , drop = FALSE],
      te = full[(nrow(Xtr) + 1L):nrow(Xall), , drop = FALSE]
    )
  })
  list(
    imps = completed,
    warnings = unique(warn_vec),
    logged_events = if (!is.null(im$loggedEvents)) as.data.table(im$loggedEvents) else data.table()
  )
}

## ---------------- data ----------------
d <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id"), names(d))[1]
if (is.na(mid)) stop("No stay_id/hadm_id found in PATH_FINAL_FULL")

bc <- fread(file.path(DIR_DATA, "baseline_covars.csv"))
need <- c(mid, "age", "gender", "sofa_score", "mortality_30d", "make30", FEATURES33)
miss <- setdiff(need, names(d))
if (length(miss)) stop("final_full missing columns: ", paste(miss, collapse = ", "))
if (!all(c(mid, "aki_stage_0_24h") %in% names(bc))) {
  stop("baseline_covars.csv must contain ", mid, " and aki_stage_0_24h")
}

M <- d[, ..need]
M[, sex := as.integer(toupper(substr(as.character(gender), 1, 1)) == "M")]
M <- merge(M, bc[, c(mid, "aki_stage_0_24h"), with = FALSE], by = mid, all.x = TRUE)
for (j in c(BASE_NUM, "mortality_30d", "make30", FEATURES33)) M[[j]] <- as.numeric(M[[j]])
M <- M[complete.cases(M[, c(BASE_NUM, "mortality_30d", "make30"), with = FALSE])]

cat(sprintf(
  "n=%d ; strict fold-wise MICE for figure (K=%d, m=%d, maxit=%d; outcomes=mortality_30d+make30)\n",
  nrow(M), K_FOLDS, MICE_M, MICE_MAXIT
))

## ---------------- core analysis ----------------
run_outcome <- function(yname) {
  y <- as.integer(M[[yname]])
  if (!all(y %in% c(0L, 1L))) stop("Outcome must be binary 0/1: ", yname)

  folds <- read_or_make_folds(y, yname)
  P <- read_or_init_predictions(yname, length(y))

  for (f in seq_len(K_FOLDS)) {
    te <- which(folds == f)
    already_done <- all(complete.cases(P[te, , drop = FALSE]))
    if (RESUME_FROM_CHECKPOINT && already_done) {
      cat(sprintf("  [%s] fold %d/%d already complete; skipping\n", yname, f, K_FOLDS))
      next
    }

    cat(sprintf("  [%s] fold %d/%d\n", yname, f, K_FOLDS))
    tr <- which(folds != f)
    ytr <- y[tr]

    Btr <- as.matrix(M[tr, ..BASE_NUM]); storage.mode(Btr) <- "double"
    Bte <- as.matrix(M[te, ..BASE_NUM]); storage.mode(Bte) <- "double"

    ## Train-only winsor bounds.
    Xtr0 <- as.data.frame(M[tr, ..FEATURES33])
    Xte0 <- as.data.frame(M[te, ..FEATURES33])
    for (j in FEATURES33) {
      q <- safe_winq(Xtr0[[j]])
      Xtr0[[j]] <- safe_clip(Xtr0[[j]], q)
      Xte0[[j]] <- safe_clip(Xte0[[j]], q)
    }

    imp_seed <- 20240601 + 1000L * match(yname, OUTCOMES) + f
    imp_res <- impute_train_test(Xtr0, Xte0, seed = imp_seed)
    imps <- imp_res$imps

    Pm <- array(
      NA_real_, dim = c(length(te), length(PRED_COLS), MICE_M),
      dimnames = list(NULL, PRED_COLS, paste0("imp", seq_len(MICE_M)))
    )
    fold_zero_sd <- character(0)

    for (k in seq_len(MICE_M)) {
      Xtr <- as.matrix(imps[[k]]$tr); storage.mode(Xtr) <- "double"
      Xte <- as.matrix(imps[[k]]$te); storage.mode(Xte) <- "double"

      sc <- safe_scale_train_test(Xtr, Xte)
      Ztr <- sc$Ztr; Zte <- sc$Zte
      fold_zero_sd <- union(fold_zero_sd, sc$zero_sd)

      set.seed(20240601 + 100L * f + k)
      km <- kmeans(Ztr, 2, nstart = NSTART_KM, iter.max = 100, algorithm = "Lloyd")
      hi_cl <- which.max(tapply(ytr, km$cluster, mean))
      cen <- km$centers
      ltr <- factor(ifelse(km$cluster == hi_cl, 1L, 2L), levels = c(1, 2))
      lte <- factor(ifelse(max.col(-d2_euclidean(Zte, cen), ties.method = "first") == hi_cl, 1L, 2L), levels = c(1, 2))

      ## Fit models in training fold only.
      df0 <- data.frame(ytr = ytr, Btr, check.names = FALSE)
      df1 <- data.frame(ytr = ytr, Btr, lab = ltr, check.names = FALSE)
      g0 <- suppressWarnings(glm(ytr ~ ., data = df0, family = binomial))
      g1 <- suppressWarnings(glm(ytr ~ ., data = df1, family = binomial))

      Df  <- cbind(Btr, Ztr)
      Dfe <- cbind(Bte, Zte)
      lab_tr_mat <- model.matrix(~ ltr)[, -1, drop = FALSE]
      lab_te_mat <- model.matrix(~ lte)[, -1, drop = FALSE]

      cvf <- cv.glmnet(Df, ytr, family = "binomial", alpha = ENET_ALPHA, nfolds = INNER_CV)
      cvfl <- cv.glmnet(cbind(Df, lab_tr_mat), ytr, family = "binomial", alpha = ENET_ALPHA, nfolds = INNER_CV)

      Pm[, "base", k] <- as.numeric(predict(g0, newdata = data.frame(Bte, check.names = FALSE), type = "response"))
      Pm[, "base_phen", k] <- as.numeric(predict(g1, newdata = data.frame(Bte, lab = lte, check.names = FALSE), type = "response"))
      Pm[, "feat", k] <- as.numeric(predict(cvf, Dfe, s = "lambda.min", type = "response"))
      Pm[, "feat_lab", k] <- as.numeric(predict(cvfl, cbind(Dfe, lab_te_mat), s = "lambda.min", type = "response"))
    }

    ## Probability-level averaging over m completed training/test imputations.
    P[te, ] <- apply(Pm, c(1, 2), mean, na.rm = TRUE)
    save_predictions(yname, y, folds, P)

    fold_log <- data.table(
      outcome = yname,
      fold = f,
      n_train = length(tr),
      n_test = length(te),
      events_train = sum(y[tr] == 1),
      events_test = sum(y[te] == 1),
      mice_warnings_n = length(imp_res$warnings),
      mice_logged_events_n = nrow(imp_res$logged_events),
      zero_sd_features = if (length(fold_zero_sd)) paste(fold_zero_sd, collapse = ";") else "",
      K_folds = K_FOLDS,
      MICE_M = MICE_M,
      MICE_MAXIT = MICE_MAXIT,
      QUICK_RUN = QUICK_RUN,
      completed_at = as.character(Sys.time())
    )
    append_fold_log(fold_log)

    if (length(imp_res$warnings)) {
      warning_file <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_fold%d_mice_warnings.txt", yname, f))
      writeLines(imp_res$warnings, warning_file)
    }
    if (nrow(imp_res$logged_events)) {
      events_file <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_fold%d_mice_loggedEvents.csv", yname, f))
      fwrite(imp_res$logged_events, events_file)
    }
  }

  if (any(!complete.cases(P))) {
    stop("Not all predictions complete for ", yname, ". Check checkpoint files and rerun with RESUME_FROM_CHECKPOINT=TRUE.")
  }

  dec <- fauc(y, P[, "feat_lab"]) - fauc(y, P[, "feat"])
  ci <- boot_ci_delta(y, P[, "feat_lab"], P[, "feat"], B = BOOT_B)

  res <- data.table(
    outcome = yname,
    outcome_label = fifelse(yname == "mortality_30d", "30-day mortality", "MAKE-30"),
    imputation = sprintf("strict fold-wise MICE probability-averaged (m=%d)", MICE_M),
    K_folds = K_FOLDS,
    auc_base = round(fauc(y, P[, "base"]), 4),
    auc_base_phen = round(fauc(y, P[, "base_phen"]), 4),
    incr_base_to_basephen = round(fauc(y, P[, "base_phen"]) - fauc(y, P[, "base"]), 4),
    auc_feat = round(fauc(y, P[, "feat"]), 4),
    auc_feat_lab = round(fauc(y, P[, "feat_lab"]), 4),
    decisive_dAUC = round(dec, 4),
    dAUC_lo = round(as.numeric(ci[1]), 4),
    dAUC_hi = round(as.numeric(ci[2]), 4)
  )
  fwrite(res, file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_summary.csv", yname)))
  res
}

res_list <- lapply(OUTCOMES, run_outcome)
res <- rbindlist(res_list, fill = TRUE)

cat("\n===== CROSS-FITTED STRICT FOLD-WISE MICE K2 LABEL B2 =====\n")
print(res)
cat("\nReference full-cohort-label MICE Rubin: mortality decisive ≈ 0.00001; make30 ≈ -0.00008\n")
cat("Interpretation: if decisive_dAUC remains ~0 and CI includes 0, redundancy is not explained by full-cohort label leakage or single-matrix imputation.\n")

fwrite(res, file.path(OUT_DIR, "crossfit_k2_b2_summary.csv"))

pred_all <- rbindlist(lapply(OUTCOMES, function(yname) {
  f <- file.path(OUT_DIR, sprintf("crossfit_k2_b2_%s_predictions.csv", yname))
  if (file.exists(f)) fread(f) else NULL
}), fill = TRUE)
if (nrow(pred_all)) fwrite(pred_all, file.path(OUT_DIR, "crossfit_k2_b2_all_predictions.csv"))

## ---------------- figure outputs ----------------
auc_long <- melt(
  res,
  id.vars = c("outcome", "outcome_label"),
  measure.vars = c("auc_base", "auc_base_phen", "auc_feat", "auc_feat_lab"),
  variable.name = "model", value.name = "AUC"
)
auc_long[, model := factor(model,
  levels = c("auc_base", "auc_base_phen", "auc_feat", "auc_feat_lab"),
  labels = c("Baseline", "Baseline + K2", "Raw33 elastic-net", "Raw33 elastic-net + K2")
)]
fwrite(auc_long, file.path(OUT_DIR, "crossfit_k2_b2_figdata_auc_long.csv"))

fig_auc <- ggplot(auc_long, aes(x = model, y = AUC, group = outcome_label)) +
  geom_point(size = 2.4) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ outcome_label, nrow = 1) +
  coord_cartesian(ylim = c(max(0.50, min(auc_long$AUC, na.rm = TRUE) - 0.03), min(0.95, max(auc_long$AUC, na.rm = TRUE) + 0.03))) +
  labs(x = NULL, y = "Cross-fitted AUC", title = "Strict fold-wise MICE cross-fitted benchmark") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))

dauc <- res[, .(outcome, outcome_label, decisive_dAUC, dAUC_lo, dAUC_hi)]
fwrite(dauc, file.path(OUT_DIR, "crossfit_k2_b2_figdata_decisive_dauc.csv"))
fig_dauc <- ggplot(dauc, aes(y = outcome_label, x = decisive_dAUC, xmin = dAUC_lo, xmax = dAUC_hi)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  labs(x = "ΔAUC: raw33 elastic-net + K2 minus raw33 elastic-net", y = NULL,
       title = "Decisive contrast under strict cross-fitting") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

for (ext in c("png", "pdf")) {
  ggsave(file.path(FIG_DIR, paste0("crossfit_k2_b2_auc_ladder.", ext)), fig_auc, width = 8.5, height = 4.2, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0("crossfit_k2_b2_decisive_dauc_forest.", ext)), fig_dauc, width = 7.2, height = 3.2, dpi = 300)
}

cat("Done. Wrote:\n")
cat(" - output/qc/crossfit_k2_b2_summary.csv\n")
cat(" - output/qc/crossfit_k2_b2_<outcome>_summary.csv\n")
cat(" - output/qc/crossfit_k2_b2_<outcome>_predictions.csv\n")
cat(" - output/qc/crossfit_k2_b2_fold_log.csv\n")
cat(" - output/qc/crossfit_k2_b2_figdata_auc_long.csv\n")
cat(" - output/qc/crossfit_k2_b2_figdata_decisive_dauc.csv\n")
cat(" - output/fig/crossfit_k2_b2_auc_ladder.(png/pdf)\n")
cat(" - output/fig/crossfit_k2_b2_decisive_dauc_forest.(png/pdf)\n")
