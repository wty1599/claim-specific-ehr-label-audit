## =====================================================================
## 26_landmark_24h_sensitivity.R
## Landmark sensitivity for first-24h physiology models
##
## Rationale:
##   If a patient dies early, leaves ICU early, or receives RRT during the first
##   24h, the 0-24h feature window may already contain part of the outcome
##   process. This script performs a 24h landmark sensitivity analysis.
##
## Landmark cohort:
##   Core landmark: alive at 24h AND ICU length of stay >= 24h.
##   RRT/MAKE landmark: core landmark AND no RRT in first 24h by default.
##
## Post-landmark outcomes:
##   1) death_after24_30d : 30-day mortality among patients alive at 24h
##   2) rrt_after24_30d   : RRT after 24h through 30d among patients without RRT_24h
##   3) make_after24_30d  : death_after24_30d OR rrt_after24_30d OR persistent RD30,
##                          among patients without RRT_24h by default
##
## Models evaluated by 10-fold CV:
##   base      = age + sex + SOFA + AKI stage
##   base_phen = base + MICE K=2 label
##   feat_pen  = base + 33 first-24h physiological features, elastic-net logistic
##   feat_phen = feat_pen + MICE K=2 label
##
## Key contrasts:
##   base_phen - base      : apparent gain over coarse clinical baseline
##   feat_phen - feat_pen  : decisive raw-feature benchmark contrast
##
## Inputs expected:
##   output/model/X_primary33_std_mice.rds
##   output/model/labels_primary_mice.rds
##   data/final_full via PATH_FINAL_FULL
##   data/baseline_covars.csv
##
## Optional SQL fallback:
##   If final_full does not contain required landmark columns, set USE_SQL=TRUE
##   and configure DB connection variables below or through environment variables.
## =====================================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({
  need <- c("data.table", "glmnet")
  for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(glmnet)
})

set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
FIG_DIR <- if (exists("DIR_FIG")) DIR_FIG else file.path(DIR_OUTPUT, "fig")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- user controls ----------------
K_FOLDS <- 10L
BOOT_B  <- 1000L
ENET_ALPHA <- 0.5
INNER_CV <- 5L
USE_NO_RRT24_FOR_RRT_AND_MAKE <- TRUE
USE_SQL_IF_NEEDED <- FALSE
SQL_SCHEMA <- Sys.getenv("EHR_AUDIT_SQL_SCHEMA", "project_sa_aki")
SQL_TABLE  <- Sys.getenv("EHR_AUDIT_SQL_TABLE", "final_full")

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
PRED_COLS <- c("base", "base_phen", "feat_pen", "feat_phen")

## ---------------- helpers ----------------
fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok]); p <- as.numeric(p[ok])
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  mean((as.numeric(p[ok]) - as.integer(y[ok]))^2)
}

clamp <- function(p, eps = 1e-6) pmin(pmax(as.numeric(p), eps), 1 - eps)

cal_slope_int <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok]); p <- clamp(p[ok])
  if (length(unique(y)) < 2) return(c(slope = NA_real_, intercept = NA_real_))
  lp <- qlogis(p)
  slope <- tryCatch(unname(coef(glm(y ~ lp, family = binomial))[2]), error = function(e) NA_real_)
  intercept <- tryCatch(unname(coef(glm(y ~ 1, offset = lp, family = binomial))[1]), error = function(e) NA_real_)
  c(slope = slope, intercept = intercept)
}

boot_ci_delta <- function(y, p_aug, p_ref, B = BOOT_B) {
  ok <- !is.na(y) & !is.na(p_aug) & !is.na(p_ref)
  y <- as.integer(y[ok]); p_aug <- as.numeric(p_aug[ok]); p_ref <- as.numeric(p_ref[ok])
  n <- length(y)
  if (n == 0 || length(unique(y)) < 2) return(c(lo = NA_real_, hi = NA_real_))
  d <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    fast_auc(y[i], p_aug[i]) - fast_auc(y[i], p_ref[i])
  })
  qs <- quantile(d, c(.025, .975), na.rm = TRUE)
  c(lo = as.numeric(qs[1]), hi = as.numeric(qs[2]))
}

stratified_folds <- function(y, K = K_FOLDS) {
  y <- as.integer(y)
  folds <- integer(length(y))
  for (cl in sort(unique(y))) {
    ix <- which(y == cl)
    folds[ix] <- sample(rep_len(seq_len(K), length(ix)))
  }
  folds
}

parse_sex <- function(g) {
  if (is.numeric(g)) return(as.integer(g == 1))
  as.integer(toupper(substr(as.character(g), 1, 1)) == "M")
}

coalesce_cols <- function(dt, choices) {
  hit <- intersect(choices, names(dt))
  if (!length(hit)) return(rep(NA_real_, nrow(dt)))
  out <- dt[[hit[1]]]
  if (length(hit) > 1) {
    for (h in hit[-1]) {
      idx <- is.na(out)
      out[idx] <- dt[[h]][idx]
    }
  }
  out
}

read_mice_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing MICE primary matrix: ", path,
         "\nRun 23c_save_mice_matrix_labels.R or 23_mice_primary_pooling first.")
  }
  x <- readRDS(path)
  if (is.matrix(x)) {
    X <- as.data.table(x)
    id_col <- attr(x, "id_col")
    if (!is.null(attr(x, "stay_id"))) X[, stay_id := attr(x, "stay_id")]
  } else {
    X <- as.data.table(x)
  }
  id_candidates <- intersect(c("stay_id", "hadm_id", "subject_id"), names(X))
  if (!length(id_candidates)) stop("MICE matrix must contain stay_id/hadm_id/subject_id column.")
  id_col <- id_candidates[1]
  miss <- setdiff(FEATURES33, names(X))
  if (length(miss)) stop("MICE matrix missing features: ", paste(miss, collapse = ", "))
  setnames(X, id_col, "id_key")
  X[, id_key := as.character(id_key)]
  X[, c("id_key", FEATURES33), with = FALSE]
}

read_mice_labels <- function(path) {
  if (!file.exists(path)) {
    stop("Missing MICE primary labels: ", path,
         "\nRun 23c_save_mice_matrix_labels.R or 23_mice_primary_pooling first.")
  }
  lab <- as.data.table(readRDS(path))
  id_candidates <- intersect(c("stay_id", "hadm_id", "subject_id"), names(lab))
  if (!length(id_candidates)) stop("labels_primary_mice.rds must contain stay_id/hadm_id/subject_id.")
  id_col <- id_candidates[1]
  if (!"cluster_k2" %in% names(lab)) stop("labels_primary_mice.rds must contain cluster_k2.")
  setnames(lab, id_col, "id_key")
  lab[, id_key := as.character(id_key)]
  lab <- unique(lab[, .(id_key, cluster_k2)], by = "id_key")
  lab[, cluster_k2 := factor(as.character(cluster_k2), levels = sort(unique(as.character(cluster_k2))))]
  lab
}

## Optional SQL extraction. This is intentionally conservative; final_full CSV is preferred.
read_final_full_or_sql <- function() {
  ff <- fread(PATH_FINAL_FULL)
  needed <- c("icu_los_days", "days_to_death", "rrt_24h", "rrt_30d", "persistent_rd_30d",
              "mortality_30d", "make30")
  missing_needed <- setdiff(needed, names(ff))
  if (!length(missing_needed) || !USE_SQL_IF_NEEDED) return(ff)

  message("final_full missing landmark columns: ", paste(missing_needed, collapse = ", "))
  message("Trying SQL extraction because USE_SQL_IF_NEEDED=TRUE ...")
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RPostgres", quietly = TRUE)) {
    stop("DBI/RPostgres not installed. Install or set USE_SQL_IF_NEEDED=FALSE and rebuild final_full.")
  }
  con <- DBI::dbConnect(
    RPostgres::Postgres(),
    dbname = Sys.getenv("PGDATABASE", "mimiciv"),
    host = Sys.getenv("PGHOST", "localhost"),
    port = as.integer(Sys.getenv("PGPORT", "5432")),
    user = Sys.getenv("PGUSER", Sys.info()[["user"]]),
    password = Sys.getenv("PGPASSWORD", "")
  )
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  sql <- sprintf("SELECT * FROM %s.%s", SQL_SCHEMA, SQL_TABLE)
  as.data.table(DBI::dbGetQuery(con, sql))
}

## ---------------- data ----------------
ff <- read_final_full_or_sql()
id_col <- intersect(c("stay_id", "hadm_id", "subject_id"), names(ff))[1]
if (is.na(id_col)) stop("No id column found in final_full.")

bc_path <- file.path(DIR_DATA, "baseline_covars.csv")
if (!file.exists(bc_path)) stop("Missing baseline_covars.csv: ", bc_path)
bc <- fread(bc_path)
if (!all(c(id_col, "aki_stage_0_24h") %in% names(bc))) {
  stop("baseline_covars.csv must contain ", id_col, " and aki_stage_0_24h.")
}

needed_ff <- c(id_col, "age", "gender", "sofa_score", "icu_los_days", "days_to_death",
               "mortality_30d", "make30", "rrt_24h", "rrt_30d", "rrt_7d", "persistent_rd_30d")
needed_ff <- intersect(needed_ff, names(ff))
D <- ff[, ..needed_ff]
D[, sex := parse_sex(gender)]
D <- merge(D, bc[, c(id_col, "aki_stage_0_24h"), with = FALSE], by = id_col, all.x = TRUE)
setnames(D, id_col, "id_key")
D[, id_key := as.character(id_key)]

## RRT fallback if rrt_30d unavailable.
if (!"rrt_30d" %in% names(D)) {
  if ("rrt_7d" %in% names(D)) {
    warning("rrt_30d not found; using rrt_7d as fallback for post-24h RRT sensitivity.")
    D[, rrt_30d := rrt_7d]
  } else {
    warning("Neither rrt_30d nor rrt_7d found; RRT and MAKE-after24 RRT component cannot be evaluated.")
    D[, rrt_30d := NA_integer_]
  }
}
for (v in c(BASE_NUM, "icu_los_days", "days_to_death", "mortality_30d", "make30", "rrt_24h", "rrt_30d", "persistent_rd_30d")) {
  if (v %in% names(D)) D[[v]] <- suppressWarnings(as.numeric(D[[v]]))
}

X <- read_mice_matrix(file.path(DIR_MODEL, "X_primary33_std_mice.rds"))
lab <- read_mice_labels(file.path(DIR_MODEL, "labels_primary_mice.rds"))
D <- merge(D, X, by = "id_key", all.x = FALSE)
D <- merge(D, lab, by = "id_key", all.x = FALSE)

## ---------------- landmark definitions ----------------
D[, alive_at_24h := fifelse(is.na(days_to_death), TRUE, days_to_death > 1)]
D[, icu_ge_24h := !is.na(icu_los_days) & icu_los_days >= 1]
D[, no_rrt_24h := is.na(rrt_24h) | rrt_24h == 0]
D[, landmark_core := alive_at_24h & icu_ge_24h]
D[, landmark_no_rrt24 := landmark_core & no_rrt_24h]

## Post-landmark outcomes. Since landmark_core requires alive at 24h, mortality_30d
## among landmark patients is by construction death after 24h through day 30.
D[, death_after24_30d := as.integer(mortality_30d == 1)]
D[, rrt_after24_30d := as.integer(rrt_30d == 1 & no_rrt_24h)]
D[, make_after24_30d := as.integer((death_after24_30d == 1) |
                                     (rrt_after24_30d == 1) |
                                     (persistent_rd_30d == 1))]

cohort_audit <- rbindlist(list(
  data.table(step = "Merged with MICE matrix and labels", n = nrow(D)),
  data.table(step = "Alive at 24h", n = sum(D$alive_at_24h, na.rm = TRUE)),
  data.table(step = "ICU LOS >=24h", n = sum(D$icu_ge_24h, na.rm = TRUE)),
  data.table(step = "Landmark core: alive at 24h and ICU LOS >=24h", n = sum(D$landmark_core, na.rm = TRUE)),
  data.table(step = "Landmark no early RRT: core and no RRT in first 24h", n = sum(D$landmark_no_rrt24, na.rm = TRUE))
))
fwrite(cohort_audit, file.path(OUT_DIR, "landmark24_cohort_audit.csv"))
cat("\n===== 24h landmark cohort audit =====\n")
print(cohort_audit)

## ---------------- model evaluation ----------------
fit_predict_cv <- function(dat, outcome, tag) {
  y <- as.integer(dat[[outcome]])
  keep <- complete.cases(dat[, c(BASE_NUM, "cluster_k2", FEATURES33, outcome), with = FALSE]) & !is.na(y)
  dat <- dat[keep]
  y <- as.integer(dat[[outcome]])
  if (length(unique(y)) < 2) stop("Outcome has <2 classes after filtering: ", outcome)

  folds <- stratified_folds(y, K_FOLDS)
  P <- matrix(NA_real_, nrow(dat), length(PRED_COLS), dimnames = list(NULL, PRED_COLS))

  for (f in seq_len(K_FOLDS)) {
    tr <- folds != f
    te <- folds == f
    ytr <- y[tr]
    Btr <- as.matrix(dat[tr, ..BASE_NUM]); storage.mode(Btr) <- "double"
    Bte <- as.matrix(dat[te, ..BASE_NUM]); storage.mode(Bte) <- "double"
    Xtr <- as.matrix(dat[tr, ..FEATURES33]); storage.mode(Xtr) <- "double"
    Xte <- as.matrix(dat[te, ..FEATURES33]); storage.mode(Xte) <- "double"
    lab_tr <- factor(dat$cluster_k2[tr])
    lab_te <- factor(dat$cluster_k2[te], levels = levels(lab_tr))

    g0 <- suppressWarnings(glm(ytr ~ ., data = data.frame(ytr = ytr, Btr, check.names = FALSE), family = binomial))
    g1 <- suppressWarnings(glm(ytr ~ ., data = data.frame(ytr = ytr, Btr, lab = lab_tr, check.names = FALSE), family = binomial))

    Df <- cbind(Btr, Xtr)
    Dfe <- cbind(Bte, Xte)
    lab_tr_mat <- model.matrix(~ lab_tr)[, -1, drop = FALSE]
    lab_te_mat <- model.matrix(~ lab_te)[, -1, drop = FALSE]

    cvf <- cv.glmnet(Df, ytr, family = "binomial", alpha = ENET_ALPHA, nfolds = INNER_CV, type.measure = "auc")
    cvfl <- cv.glmnet(cbind(Df, lab_tr_mat), ytr, family = "binomial", alpha = ENET_ALPHA, nfolds = INNER_CV, type.measure = "auc")

    P[te, "base"] <- as.numeric(predict(g0, newdata = data.frame(Bte, check.names = FALSE), type = "response"))
    P[te, "base_phen"] <- as.numeric(predict(g1, newdata = data.frame(Bte, lab = lab_te, check.names = FALSE), type = "response"))
    P[te, "feat_pen"] <- as.numeric(predict(cvf, Dfe, s = "lambda.min", type = "response"))
    P[te, "feat_phen"] <- as.numeric(predict(cvfl, cbind(Dfe, lab_te_mat), s = "lambda.min", type = "response"))
  }

  perf <- rbindlist(lapply(PRED_COLS, function(m) {
    cs <- cal_slope_int(y, P[, m])
    data.table(
      analysis = tag,
      outcome = outcome,
      model = m,
      n = length(y),
      events = sum(y == 1),
      event_rate = round(mean(y == 1), 4),
      auc = round(fast_auc(y, P[, m]), 4),
      brier = round(brier(y, P[, m]), 4),
      calib_slope = round(cs["slope"], 3),
      calib_intercept = round(cs["intercept"], 3)
    )
  }))

  comps <- rbindlist(list(
    data.table(analysis = tag, outcome = outcome, comparison = "base_phen vs base",
               auc_ref = fast_auc(y, P[, "base"]), auc_aug = fast_auc(y, P[, "base_phen"]),
               dAUC = fast_auc(y, P[, "base_phen"]) - fast_auc(y, P[, "base"]),
               t(boot_ci_delta(y, P[, "base_phen"], P[, "base"]))),
    data.table(analysis = tag, outcome = outcome, comparison = "feat_phen vs feat_pen",
               auc_ref = fast_auc(y, P[, "feat_pen"]), auc_aug = fast_auc(y, P[, "feat_phen"]),
               dAUC = fast_auc(y, P[, "feat_phen"]) - fast_auc(y, P[, "feat_pen"]),
               t(boot_ci_delta(y, P[, "feat_phen"], P[, "feat_pen"])))
  ), fill = TRUE)
  setnames(comps, c("lo", "hi"), c("dAUC_lo", "dAUC_hi"), skip_absent = TRUE)
  for (j in c("auc_ref", "auc_aug", "dAUC", "dAUC_lo", "dAUC_hi")) comps[, (j) := round(get(j), 4)]

  pred <- data.table(id_key = dat$id_key, analysis = tag, outcome = outcome, y = y, fold = folds)
  pred <- cbind(pred, as.data.table(P))

  list(performance = perf, comparisons = comps, predictions = pred)
}

analyses <- list(
  list(outcome = "death_after24_30d", tag = "Landmark core: alive and ICU_LOS>=24h", data = D[landmark_core == TRUE]),
  list(outcome = "rrt_after24_30d", tag = "Landmark no early RRT: post-24h RRT", data = D[landmark_no_rrt24 == TRUE]),
  list(outcome = "make_after24_30d", tag = "Landmark no early RRT: post-24h MAKE", data = D[if (USE_NO_RRT24_FOR_RRT_AND_MAKE) landmark_no_rrt24 == TRUE else landmark_core == TRUE])
)

results <- lapply(analyses, function(a) {
  cat("\n===== Running landmark sensitivity: ", a$outcome, " | ", a$tag, " =====\n", sep = "")
  fit_predict_cv(a$data, a$outcome, a$tag)
})

perf_all <- rbindlist(lapply(results, `[[`, "performance"), fill = TRUE)
comp_all <- rbindlist(lapply(results, `[[`, "comparisons"), fill = TRUE)
pred_all <- rbindlist(lapply(results, `[[`, "predictions"), fill = TRUE)

fwrite(perf_all, file.path(OUT_DIR, "landmark24_model_performance.csv"))
fwrite(comp_all, file.path(OUT_DIR, "landmark24_incremental_contrasts.csv"))
fwrite(pred_all, file.path(OUT_DIR, "landmark24_predictions.csv"))

cat("\n===== Landmark 24h model performance =====\n")
print(perf_all)
cat("\n===== Landmark 24h incremental contrasts =====\n")
print(comp_all)

cat("\nInterpretation guide:\n")
cat("1) base_phen vs base estimates whether the K=2 label still adds apparent value over a coarse baseline after excluding early deaths/early ICU discharge.\n")
cat("2) feat_phen vs feat_pen is the decisive benchmark. If dAUC remains near 0 and CI includes 0, the K=2 label remains redundant beyond the same 33 raw physiological features even under the 24h landmark design.\n")
cat("3) RRT/MAKE analyses are run in the no-early-RRT landmark set by default, reducing the risk that 0-24h RRT is already part of the outcome process.\n")

## Optional simple forest plot for decisive dAUC.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  forest <- comp_all[comparison == "feat_phen vs feat_pen"]
  forest[, label := fifelse(outcome == "death_after24_30d", "30-day death after 24h",
                            fifelse(outcome == "rrt_after24_30d", "RRT after 24h", "MAKE after 24h"))]
  p <- ggplot(forest, aes(x = dAUC, y = label)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = dAUC_lo, xmax = dAUC_hi), height = 0.16) +
    geom_point(size = 2.4) +
    labs(x = "Delta AUC: feat + K=2 label minus feat-only model",
         y = NULL,
         title = "24h landmark sensitivity: incremental value beyond raw physiology") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG_DIR, "landmark24_decisive_dauc_forest.png"), p, width = 7.2, height = 3.8, dpi = 300)
  ggsave(file.path(FIG_DIR, "landmark24_decisive_dauc_forest.pdf"), p, width = 7.2, height = 3.8)
}

cat("\nDone. Wrote:\n")
cat(" - output/qc/landmark24_cohort_audit.csv\n")
cat(" - output/qc/landmark24_model_performance.csv\n")
cat(" - output/qc/landmark24_incremental_contrasts.csv\n")
cat(" - output/qc/landmark24_predictions.csv\n")
cat(" - output/fig/landmark24_decisive_dauc_forest.[png/pdf] if ggplot2 is available\n")
