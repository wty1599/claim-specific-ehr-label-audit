## =====================================================================
## 22_sensitivity_aggregate_mice_safe.R
## Robustness aggregation for the MICE primary framework:
##   (a) discreteness across available feature matrices;
##   (b) B2 redundancy across available clustering-label variants.
##
## Changes from older 22_sensitivity_aggregate.R:
##   - Requires X_primary33_std_mice.rds for the primary33 matrix.
##   - Uses labels_primary_mice.rds for kmeans_k2.
##   - Does not assume cluster_k3 exists in labels_primary_mice.rds; it is added
##     only when present, otherwise skipped cleanly.
##   - Writes a provenance/audit table of all matrices and label files.
##   - Avoids rbindlist failures when all optional sensitivity inputs are absent.
## =====================================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({
  need <- c("data.table", "cluster", "diptest", "MASS", "glmnet")
  for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(cluster)
  library(diptest)
  library(MASS)
  library(glmnet)
})
set.seed(20240601)
OUT_DIR <- file.path(DIR_OUTPUT, "qc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- configuration ----------------
PRIMARY_MICE <- file.path(DIR_MODEL, "X_primary33_std_mice.rds")
if (!file.exists(PRIMARY_MICE)) {
  stop("Missing MICE primary matrix: ", PRIMARY_MICE,
       "\nRun 23_mice_primary_pooling_v2_fixed.R or 23c_save_mice_matrix_labels.R first.")
}

MANUAL_FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
  "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml", "glucose_max",
  "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min", "sodium_min",
  "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min", "chloride_min",
  "chloride_max", "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)

load_primary_features <- function() {
  fset_path <- file.path(DIR_MODEL, "feature_sets.rds")
  if (file.exists(fset_path)) {
    fs <- readRDS(fset_path)
    if (!is.null(fs$primary) && length(fs$primary) == 33L) {
      message("Using authoritative feature_sets.rds$primary")
      return(as.character(fs$primary))
    }
    warning("feature_sets.rds exists but $primary is missing or not length 33; using manual FEATURES33.")
  } else {
    warning("feature_sets.rds not found; using manual FEATURES33.")
  }
  MANUAL_FEATURES33
}

FEATURES33 <- load_primary_features()
stopifnot(length(FEATURES33) == 33L)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")

FEATURE_MATRICES <- list(
  primary33_mice = PRIMARY_MICE,
  restricted18   = file.path(DIR_MODEL, "X_restricted18_std.rds"),
  inclusive42    = file.path(DIR_MODEL, "X_inclusive42_std.rds")
)

KM <- function(m, k, ns = 15) kmeans(m, k, nstart = ns, iter.max = 100, algorithm = "Lloyd")
fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]
  p <- p[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
parse_sex <- function(g) {
  if (is.numeric(g)) return(as.integer(g == 1))
  x <- toupper(substr(trimws(as.character(g)), 1, 1))
  fifelse(x == "M", 1L, fifelse(x == "F", 0L, NA_integer_))
}

rbind_nonempty <- function(x, fill = TRUE) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(data.table())
  rbindlist(x, fill = fill, use.names = TRUE)
}

check_unique_id <- function(dt, id_col, name) {
  if (!(id_col %in% names(dt))) stop(name, " missing ID column: ", id_col)
  nd <- anyDuplicated(dt[[id_col]])
  if (nd) stop(name, " contains duplicated ", id_col, "; first duplicated index = ", nd)
}

normalize_id_col <- function(dt, id_col, name) {
  id_candidates <- intersect(c(id_col, "stay_id", "hadm_id", "subject_id"), names(dt))
  if (!(id_col %in% names(dt))) {
    if (!length(id_candidates)) stop(name, " lacks an ID column compatible with ", id_col)
    setnames(dt, id_candidates[1], id_col)
  }
  dt
}

## ---------- label-file registry: k3 only if actually present ----------
make_label_registry <- function(mid) {
  registry <- list()

  mice_label_file <- file.path(DIR_MODEL, "labels_primary_mice.rds")
  if (!file.exists(mice_label_file)) {
    stop("Missing MICE label file: ", mice_label_file,
         "\nRun 23_mice_primary_pooling_v2_fixed.R or 23c_save_mice_matrix_labels.R first.")
  }
  lab_mice <- as.data.table(readRDS(mice_label_file))
  lab_mice <- normalize_id_col(lab_mice, mid, "labels_primary_mice.rds")
  check_unique_id(lab_mice, mid, "labels_primary_mice.rds")
  if (!"cluster_k2" %in% names(lab_mice)) stop("labels_primary_mice.rds lacks cluster_k2.")

  registry$kmeans_k2_mice <- list(file = mice_label_file, col = "cluster_k2", required = TRUE)
  if ("cluster_k3" %in% names(lab_mice)) {
    registry$kmeans_k3_mice <- list(file = mice_label_file, col = "cluster_k3", required = FALSE)
  } else {
    message("labels_primary_mice.rds has no cluster_k3; kmeans_k3_mice sensitivity will be skipped.")
  }

  optional <- list(
    gmm_k2  = list(file = file.path(DIR_MODEL, "labels_gmm.rds"),  col = "cluster", required = FALSE),
    ward_k2 = list(file = file.path(DIR_MODEL, "labels_ward.rds"), col = "cluster", required = FALSE)
  )
  c(registry, optional)
}

## ---------- provenance audit ----------
write_initial_audit <- function(label_files) {
  matrix_audit <- rbindlist(lapply(names(FEATURE_MATRICES), function(nm) {
    data.table(type = "matrix", name = nm, path = FEATURE_MATRICES[[nm]], exists = file.exists(FEATURE_MATRICES[[nm]]))
  }))
  label_audit <- rbindlist(lapply(names(label_files), function(nm) {
    data.table(type = "label", name = nm, path = label_files[[nm]]$file,
               column = label_files[[nm]]$col, exists = file.exists(label_files[[nm]]$file))
  }), fill = TRUE)
  out <- rbind(matrix_audit, label_audit, fill = TRUE)
  fwrite(out, file.path(OUT_DIR, "sensitivity_mice_provenance_audit.csv"))
  invisible(out)
}

## ---------- (a) discreteness per available feature matrix ----------
disc <- rbind_nonempty(lapply(names(FEATURE_MATRICES), function(nm) {
  f <- FEATURE_MATRICES[[nm]]
  if (!file.exists(f)) {
    cat("[skip disc]", nm, "missing file\n")
    return(NULL)
  }
  Z <- as.data.table(readRDS(f))
  id_cols <- intersect(c("stay_id", "hadm_id", "subject_id"), names(Z))
  feats <- setdiff(names(Z), id_cols)
  if (nm == "primary33_mice") {
    miss <- setdiff(FEATURES33, names(Z))
    if (length(miss)) stop("Primary MICE matrix missing features: ", paste(miss, collapse = ", "))
    feats <- FEATURES33
  }
  m <- as.matrix(Z[, ..feats])
  storage.mode(m) <- "double"
  if (anyNA(m)) {
    cat("[skip disc]", nm, "contains NA after loading\n")
    return(NULL)
  }
  m <- m[complete.cases(m), , drop = FALSE]
  if (nrow(m) < 100 || ncol(m) < 2) {
    cat("[skip disc]", nm, "too small after cleaning\n")
    return(NULL)
  }

  km <- KM(m, 2)
  axis <- as.numeric(m %*% (km$centers[1, ] - km$centers[2, ]))
  m_gap <- m[sample.int(nrow(m), min(3000, nrow(m))), , drop = FALSE]
  cg <- clusGap(m_gap, function(x, k) KM(x, k, 8), K.max = 12, B = 50,
                d.power = 2, spaceH0 = "scaledPCA")
  ci <- km$tot.withinss / sum(sweep(m, 2, colMeans(m))^2)
  Sig <- cov(m)
  cn <- replicate(80, {
    s <- MASS::mvrnorm(nrow(m), colMeans(m), Sig)
    k2 <- KM(s, 2)
    k2$tot.withinss / sum(sweep(s, 2, colMeans(s))^2)
  })
  data.table(
    feature_set = nm,
    n = nrow(m),
    d = ncol(m),
    dip_axis_p = signif(diptest::dip.test(axis)$p.value, 3),
    gap_selected_K = cluster::maxSE(cg$Tab[, "gap"], cg$Tab[, "SE.sim"], "firstSEmax"),
    gap_rises_to_max = as.logical(which.max(cg$Tab[, "gap"]) == 12),
    sig_effect_pct = round(100 * (mean(cn) - ci) / mean(cn), 3),
    matrix_source = basename(f)
  )
}))
cat("\n===== (a) Discreteness robustness across feature sets =====\n")
if (nrow(disc)) print(disc) else cat("No discreteness sensitivity result generated.\n")

## ---------- base table for B2 ----------
mimic <- fread(PATH_FINAL_FULL)
mid <- intersect(c("stay_id", "hadm_id", "subject_id"), names(mimic))[1]
if (is.na(mid)) stop("No stay_id/hadm_id/subject_id found in PATH_FINAL_FULL.")
check_unique_id(mimic, mid, "PATH_FINAL_FULL")

bc <- fread(file.path(DIR_DATA, "baseline_covars.csv"))
check_unique_id(bc, mid, "baseline_covars.csv")
if (!"aki_stage_0_24h" %in% names(bc)) stop("baseline_covars.csv missing aki_stage_0_24h")

base <- mimic[, .SD, .SDcols = c(mid, "age", "gender", "sofa_score", "mortality_30d", "make30")]
base[, sex := parse_sex(gender)]
base <- merge(base, bc[, .SD, .SDcols = c(mid, "aki_stage_0_24h")], by = mid, all.x = TRUE)

Xp <- as.data.table(readRDS(PRIMARY_MICE))
Xp <- normalize_id_col(Xp, mid, "X_primary33_std_mice.rds")
check_unique_id(Xp, mid, "X_primary33_std_mice.rds")
stopifnot(all(FEATURES33 %in% names(Xp)))
base <- merge(base, Xp[, .SD, .SDcols = c(mid, FEATURES33)], by = mid, all.x = TRUE)

LABEL_FILES <- make_label_registry(mid)
write_initial_audit(LABEL_FILES)

b2 <- rbind_nonempty(lapply(names(LABEL_FILES), function(nm) {
  L <- LABEL_FILES[[nm]]
  if (!file.exists(L$file)) {
    cat("[skip B2]", nm, "missing file\n")
    return(NULL)
  }
  lab <- as.data.table(readRDS(L$file))
  lab <- normalize_id_col(lab, mid, basename(L$file))
  if (!(L$col %in% names(lab))) {
    cat("[skip B2]", nm, "missing column", L$col, "\n")
    return(NULL)
  }
  check_unique_id(lab, mid, basename(L$file))
  lab <- lab[, .SD, .SDcols = c(mid, L$col)]
  setnames(lab, L$col, "lab")

  D <- merge(base, lab, by = mid, all.x = TRUE)
  D <- D[complete.cases(D[, c("age", "sex", "sofa_score", "aki_stage_0_24h",
                              "mortality_30d", "make30", "lab", FEATURES33), with = FALSE])]
  if (nrow(D) < 100 || length(unique(D$lab)) < 2) {
    cat("[skip B2]", nm, "insufficient rows or label levels after filtering\n")
    return(NULL)
  }
  D[, lab := factor(lab)]

  one <- function(yn) {
    y <- as.integer(D[[yn]])
    if (anyNA(y) || !all(y %in% c(0, 1))) stop("Outcome must be binary without NA: ", yn)
    XF <- as.matrix(D[, .SD, .SDcols = c("age", "sex", "sofa_score", "aki_stage_0_24h", FEATURES33)])
    storage.mode(XF) <- "double"
    XL <- model.matrix(~ lab, D)[, -1, drop = FALSE]
    folds <- integer(length(y))
    for (cl in c(0, 1)) {
      ix <- which(y == cl)
      folds[ix] <- sample(rep_len(seq_len(10), length(ix)))
    }
    cvp <- function(dz) {
      p <- numeric(length(y))
      for (f in seq_len(10)) {
        tr <- folds != f
        te <- folds == f
        m <- cv.glmnet(dz[tr, , drop = FALSE], y[tr], family = "binomial", alpha = .5, nfolds = 5)
        p[te] <- as.numeric(predict(m, dz[te, , drop = FALSE], s = "lambda.min", type = "response"))
      }
      p
    }
    p0 <- cvp(XF)
    p1 <- cvp(cbind(XF, XL))
    c(auc_feat = fast_auc(y, p0), auc_feat_lab = fast_auc(y, p1), dAUC = fast_auc(y, p1) - fast_auc(y, p0))
  }

  m30 <- one("mortality_30d")
  mk <- one("make30")
  if (m30["auc_feat"] < 0.75) {
    stop(sprintf("[validation failed] %s: mort_auc_feat=%.4f is far below expected (~0.806). Check PRIMARY_MICE alignment.",
                 nm, m30["auc_feat"]))
  }
  data.table(
    clustering = nm,
    label_source = basename(L$file),
    label_column = L$col,
    n = nrow(D),
    mort_auc_feat = round(m30["auc_feat"], 4),
    mort_auc_feat_lab = round(m30["auc_feat_lab"], 4),
    mort_dAUC = round(m30["dAUC"], 4),
    make_auc_feat = round(mk["auc_feat"], 4),
    make_auc_feat_lab = round(mk["auc_feat_lab"], 4),
    make_dAUC = round(mk["dAUC"], 4)
  )
}))

cat("\n===== (b) B2 redundancy across clustering variants (does ANY label add over 33 feats?) =====\n")
if (nrow(b2)) print(b2) else cat("No B2 sensitivity result generated.\n")

if (nrow(disc)) fwrite(disc, file.path(OUT_DIR, "sens_discreteness.csv"))
if (nrow(b2)) fwrite(b2, file.path(OUT_DIR, "sens_b2_variants.csv"))
cat("\nValidation target: mort_auc_feat should be around 0.806 rather than baseline-only values.\n")
cat("cluster_k3 is now optional: it is analyzed only if labels_primary_mice.rds actually contains cluster_k3.\n")
cat("Done.\n")
