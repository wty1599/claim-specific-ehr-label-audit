options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(MASS)
  library(mice)
})

if (!exists("D3R_SCRIPT_DIR", inherits = FALSE)) {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(this_file) || !nzchar(this_file)) {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(arg)) stop("Cannot resolve D3 simulation script directory.")
    this_file <- sub("^--file=", "", arg[1])
  }
  D3R_SCRIPT_DIR <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))
}
D3R_ROOT <- dirname(D3R_SCRIPT_DIR)
D3R_PROJECT_ROOT <- dirname(dirname(D3R_ROOT))
D3R_SPEC <- file.path(D3R_ROOT, "spec", "D3_SIM_PRERUN_SPEC_FROZEN_20260802.md")
D3R_EXPECTED_HASHES <- file.path(D3R_ROOT, "spec", "D3_SIM_PUBLIC_SOURCE_HASHES_v1.csv")

D3R_RUN_MODE <- Sys.getenv("D3_SIM_RUN_MODE", unset = "smoke")
if (!D3R_RUN_MODE %in% c("smoke", "formal")) stop("D3_SIM_RUN_MODE must be smoke or formal.")
D3R_OUTPUT_ROOT <- Sys.getenv(
  "D3_SIM_OUTPUT_ROOT",
  unset = if (D3R_RUN_MODE == "formal")
    file.path(D3R_PROJECT_ROOT, "outputs", "simulation", "D3_G4abc_formal") else
    file.path(D3R_PROJECT_ROOT, "outputs", "simulation", paste0("D3_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S")))
)
D3R_TABLE_DIR <- file.path(D3R_OUTPUT_ROOT, "tables")
D3R_FIGURE_DIR <- file.path(D3R_OUTPUT_ROOT, "figures")
D3R_LOG_DIR <- file.path(D3R_OUTPUT_ROOT, "logs")
D3R_PROVENANCE_DIR <- file.path(D3R_OUTPUT_ROOT, "provenance")
D3R_CHECKPOINT_DIR <- file.path(D3R_OUTPUT_ROOT, "checkpoints")

D3R_BASE_SEED <- 20260802L
D3R_N_TRAIN <- if (D3R_RUN_MODE == "formal") 5000L else 1000L
D3R_N_EXTERNAL <- if (D3R_RUN_MODE == "formal") 3000L else 600L
D3R_N_REP <- if (D3R_RUN_MODE == "formal") 100L else 2L
D3R_MICE_M <- 5L
D3R_MICE_MAXIT <- 5L
D3R_PRIMARY_IMPUTATION <- 1L
D3R_CV_FOLDS <- 5L
D3R_WORKERS <- as.integer(Sys.getenv(
  "D3_SIM_WORKERS", unset = if (D3R_RUN_MODE == "formal") "4" else "2"
))

D3R_SCENARIOS <- c(
  same = "G3R_same_distribution_no_alert",
  drift = "G3R_covariate_drift_alert",
  cal = "D3CQ_calibration_stress_control"
)
D3R_EXPECTED_ALERT <- c(
  G3R_same_distribution_no_alert = FALSE,
  G3R_covariate_drift_alert = TRUE,
  D3CQ_calibration_stress_control = TRUE
)

D3R_PREVALENCE_GATE <- 0.10
D3R_ASSIGNMENT_ARI_GATE <- 0.80
D3R_CALIBRATION_SLOPE_RANGE <- c(0.80, 1.20)
D3R_CALIBRATION_INTERCEPT_GATE <- 0.20
D3R_WINSOR_PROBS <- c(0.01, 0.99)
D3R_K <- 2L
D3R_KMEANS_NSTART <- 50L
D3R_KMEANS_ITER_MAX <- 100L

D3R_FEATURES <- c(
  "aniongap_max", "creatinine_max", "bun_max", "lactate_max",
  "bicarbonate_min", "ph_min", "urine_output_24h", "mbp_min",
  "pao2fio2_min", "spo2_min", "heart_rate_max", "resp_rate_max",
  "wbc_max", "platelets_min", "bilirubin_max", "inr_max", "ptt_max",
  "potassium_max", "sodium_max", "temperature_max", "gcs_min",
  "hemoglobin_min", "glucose_max", "chloride_max", "alt_max", "ast_max",
  "albumin_min", "calcium_min", "magnesium_max", "vasopressor_dose_max",
  "shock_index_max", "baseexcess_min", "fio2_max"
)
D3R_BASELINE <- c("age", "sex", "sofa", "aki_stage")
stopifnot(length(D3R_FEATURES) == 33L, !anyDuplicated(D3R_FEATURES))

D3R_SEVERITY_LOADING <- c(
  aniongap_max = 0.80, creatinine_max = 0.95, bun_max = 0.90,
  lactate_max = 0.85, bicarbonate_min = -0.85, ph_min = -0.75,
  urine_output_24h = -0.95, mbp_min = -0.70, pao2fio2_min = -0.65,
  spo2_min = -0.35, heart_rate_max = 0.55, resp_rate_max = 0.50,
  wbc_max = 0.35, platelets_min = -0.40, bilirubin_max = 0.45,
  inr_max = 0.45, ptt_max = 0.40, potassium_max = 0.35,
  sodium_max = 0.10, temperature_max = 0.25, gcs_min = -0.60,
  hemoglobin_min = -0.25, glucose_max = 0.25, chloride_max = 0.15,
  alt_max = 0.30, ast_max = 0.35, albumin_min = -0.45,
  calcium_min = -0.20, magnesium_max = 0.20,
  vasopressor_dose_max = 0.75, shock_index_max = 0.70,
  baseexcess_min = -0.80, fio2_max = 0.55
)

D3R_FEATURE_SCALE <- data.table(
  feature = D3R_FEATURES,
  mean = c(16, 2.0, 38, 2.8, 20, 7.31, 1200, 65, 230, 95, 115, 26,
           16, 155, 2.0, 1.4, 42, 4.6, 140, 37.8, 12, 10, 165, 104,
           90, 110, 2.8, 8.0, 2.2, 0.08, 1.0, -4.0, 0.45),
  sd = c(5, 1.2, 25, 2.2, 5, 0.10, 850, 14, 90, 4, 22, 7, 9, 80,
         2.5, 0.5, 20, 0.8, 5, 1.0, 3, 2.0, 80, 6, 300, 500, 0.7,
         0.7, 0.5, 0.15, 0.5, 5.0, 0.25)
)
D3R_POSITIVE_FEATURES <- setdiff(
  D3R_FEATURES, c("ph_min", "temperature_max", "baseexcess_min")
)
D3R_BOUNDED_LOW <- c(ph_min = 6.70, spo2_min = 50, gcs_min = 3, calcium_min = 4)
D3R_BOUNDED_HIGH <- c(
  ph_min = 7.60, spo2_min = 100, gcs_min = 15,
  temperature_max = 43, fio2_max = 1
)

d3r_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
d3r_clamp_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
d3r_logit <- function(p) qlogis(d3r_clamp_prob(p))
d3r_safe_plogis <- function(x) plogis(pmin(pmax(x, -30), 30))
d3r_sha256 <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
d3r_hash_object <- function(x) digest(x, algo = "sha256", serialize = TRUE)

d3r_init_dirs <- function() {
  if (D3R_RUN_MODE == "formal" && file.exists(file.path(D3R_OUTPUT_ROOT, "D3_SIMULATION_COMPLETED.ok"))) {
    stop("Formal output is already complete; refusing to overwrite: ", D3R_OUTPUT_ROOT)
  }
  for (d in c(D3R_OUTPUT_ROOT, D3R_TABLE_DIR, D3R_FIGURE_DIR, D3R_LOG_DIR,
              D3R_PROVENANCE_DIR, D3R_CHECKPOINT_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

d3r_atomic_fwrite <- function(x, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(as.data.table(x), tmp, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path)
  invisible(path)
}

d3r_atomic_write_lines <- function(x, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic write failed: ", path)
  invisible(path)
}

d3r_seed <- function(scenario, repeat_id) {
  D3R_BASE_SEED + match(scenario, unname(D3R_SCENARIOS)) * 100000L + as.integer(repeat_id) * 1000L
}

d3r_make_covariance <- function(p = length(D3R_FEATURES)) {
  groups <- rep(1:6, length.out = p)
  Sigma <- matrix(0.15, p, p)
  diag(Sigma) <- 1
  for (g in unique(groups)) {
    idx <- which(groups == g)
    Sigma[idx, idx] <- 0.45
    diag(Sigma)[idx] <- 1
  }
  Sigma
}

d3r_to_clinical_scale <- function(Z) {
  X <- as.data.table(Z)
  setnames(X, D3R_FEATURES)
  for (v in D3R_FEATURES) {
    m <- D3R_FEATURE_SCALE[feature == v, mean]
    s <- D3R_FEATURE_SCALE[feature == v, sd]
    X[[v]] <- m + s * X[[v]]
    if (v %in% D3R_POSITIVE_FEATURES) X[[v]] <- pmax(X[[v]], 0.001)
    if (v %in% names(D3R_BOUNDED_LOW)) X[[v]] <- pmax(X[[v]], D3R_BOUNDED_LOW[[v]])
    if (v %in% names(D3R_BOUNDED_HIGH)) X[[v]] <- pmin(X[[v]], D3R_BOUNDED_HIGH[[v]])
  }
  X[]
}

d3r_simulate_latent <- function(n, severity_mean = 0, severity_sd = 1,
                                external_shift = NULL) {
  eps <- MASS::mvrnorm(n, rep(0, length(D3R_FEATURES)), d3r_make_covariance())
  severity <- rnorm(n, severity_mean, severity_sd)
  renal <- 0.65 * severity + rnorm(n, 0, 0.75)
  acid <- 0.60 * severity + rnorm(n, 0, 0.75)
  shock <- 0.70 * severity + rnorm(n, 0, 0.75)
  resp <- 0.55 * severity + rnorm(n, 0, 0.80)
  inflam <- 0.40 * severity + rnorm(n, 0, 0.90)
  load <- matrix(0, nrow = length(D3R_FEATURES), ncol = 5,
                 dimnames = list(D3R_FEATURES, c("renal", "acid", "shock", "resp", "inflam")))
  load[c("creatinine_max", "bun_max", "urine_output_24h", "potassium_max"), "renal"] <- c(0.8, 0.7, -0.8, 0.4)
  load[c("lactate_max", "bicarbonate_min", "ph_min", "aniongap_max", "baseexcess_min"), "acid"] <- c(0.7, -0.7, -0.6, 0.6, -0.7)
  load[c("mbp_min", "vasopressor_dose_max", "shock_index_max", "heart_rate_max"), "shock"] <- c(-0.7, 0.8, 0.8, 0.5)
  load[c("pao2fio2_min", "spo2_min", "fio2_max", "resp_rate_max"), "resp"] <- c(-0.8, -0.4, 0.7, 0.4)
  load[c("wbc_max", "platelets_min", "bilirubin_max", "inr_max", "ptt_max", "temperature_max"), "inflam"] <- c(0.5, -0.4, 0.4, 0.4, 0.3, 0.3)
  latent <- cbind(renal, acid, shock, resp, inflam)
  Z <- matrix(0, n, length(D3R_FEATURES), dimnames = list(NULL, D3R_FEATURES))
  for (j in seq_along(D3R_FEATURES)) {
    v <- D3R_FEATURES[j]
    Z[, j] <- D3R_SEVERITY_LOADING[[v]] * severity +
      as.numeric(latent %*% load[v, ]) + 0.65 * eps[, j]
  }
  if (!is.null(external_shift)) {
    for (v in names(external_shift)) if (v %in% colnames(Z)) Z[, v] <- Z[, v] + external_shift[[v]]
    for (v in c("creatinine_max", "lactate_max", "bicarbonate_min", "mbp_min")) {
      Z[, v] <- Z[, v] + rnorm(n, 0, 0.35)
    }
  }
  list(Z = Z, severity = severity, renal = renal, shock = shock)
}

d3r_baseline <- function(severity) {
  n <- length(severity)
  age <- pmin(pmax(round(rnorm(n, 64 + 4 * severity, 13)), 18), 95)
  sex <- rbinom(n, 1, 0.55)
  sofa <- pmin(pmax(round(7 + 3 * severity + rnorm(n, 0, 2)), 0), 24)
  p3 <- d3r_safe_plogis(-1.2 + 1.25 * severity)
  p2 <- d3r_safe_plogis(-0.5 + 0.75 * severity) * (1 - p3)
  u <- runif(n)
  aki_stage <- ifelse(u < p3, 3L, ifelse(u < p3 + p2, 2L, 1L))
  data.table(age, sex, sofa, aki_stage)
}

d3r_outcome <- function(severity, base, calibration_stress = FALSE) {
  intercept <- if (calibration_stress) -1.85 else -2.35
  severity_beta <- if (calibration_stress) 1.35 else 0.90
  lp <- intercept + severity_beta * severity + 0.018 * (base$age - 65) +
    0.06 * (base$sofa - 7) + 0.18 * (base$aki_stage - 1)
  mortality_30d <- rbinom(length(severity), 1, d3r_safe_plogis(lp))
  lp_make <- -1.10 + 1.05 * severity + 0.04 * (base$sofa - 7) +
    0.45 * (base$aki_stage - 1)
  make30 <- rbinom(length(severity), 1, d3r_safe_plogis(lp_make))
  data.table(mortality_30d, make30, true_mortality_lp = lp)
}

d3r_inject_missing <- function(dt, severity, mechanism = c("MAR", "MNAR"), intensity = 0.08) {
  mechanism <- match.arg(mechanism)
  out <- copy(dt)
  for (v in D3R_FEATURES) {
    if (mechanism == "MAR") {
      p <- d3r_safe_plogis(d3r_logit(intensity) - 0.65 * severity)
    } else {
      z <- as.numeric(scale(out[[v]])); z[!is.finite(z)] <- 0
      p <- d3r_safe_plogis(d3r_logit(intensity) + 0.35 * severity + 0.25 * z)
    }
    out[rbinom(nrow(out), 1, p) == 1, (v) := NA_real_]
  }
  out[]
}

d3r_generate_dataset <- function(scenario, cohort = c("derivation", "external"), n, repeat_id, seed) {
  cohort <- match.arg(cohort)
  set.seed(seed)
  drift <- scenario == D3R_SCENARIOS[["drift"]] && cohort == "external"
  cal_stress <- scenario == D3R_SCENARIOS[["cal"]] && cohort == "external"
  shift <- if (drift) c(creatinine_max = 0.45, lactate_max = 0.40,
                        bicarbonate_min = -0.35, mbp_min = -0.35,
                        urine_output_24h = -0.25, pao2fio2_min = -0.25) else NULL
  latent <- d3r_simulate_latent(
    n, severity_mean = if (drift) 0.55 else 0,
    severity_sd = if (drift) 1.10 else 1, external_shift = shift
  )
  features_complete <- d3r_to_clinical_scale(latent$Z)
  base <- d3r_baseline(latent$severity)
  outcomes <- d3r_outcome(latent$severity, base, calibration_stress = cal_stress)
  core <- data.table(
    patient_id = sprintf("%s_%03d_%05d", ifelse(cohort == "derivation", "D", "V"), repeat_id, seq_len(n)),
    cohort, scenario, repeat_id, severity_z = latent$severity,
    latent_renal = latent$renal, latent_shock = latent$shock
  )
  complete <- cbind(core, copy(features_complete), base, outcomes)
  incomplete_features <- d3r_inject_missing(
    features_complete, latent$severity,
    mechanism = if (drift) "MNAR" else "MAR",
    intensity = if (drift) 0.18 else 0.08
  )
  raw <- cbind(core, incomplete_features, base, outcomes)
  list(raw = raw, complete = complete)
}

d3r_mice <- function(dt, seed) {
  vars <- c(D3R_FEATURES, D3R_BASELINE)
  imp_data <- as.data.frame(dt[, ..vars])
  ini <- mice::mice(imp_data, maxit = 0, printFlag = FALSE)
  imp <- mice::mice(
    imp_data, m = D3R_MICE_M, maxit = D3R_MICE_MAXIT,
    method = ini$method, predictorMatrix = ini$predictorMatrix,
    seed = seed, printFlag = FALSE
  )
  lapply(seq_len(D3R_MICE_M), function(i) {
    comp <- as.data.table(mice::complete(imp, i))
    out <- copy(dt)
    for (v in vars) set(out, j = v, value = comp[[v]])
    out
  })
}

d3r_winsor_fit <- function(dt) rbindlist(lapply(D3R_FEATURES, function(v) {
  q <- quantile(dt[[v]], D3R_WINSOR_PROBS, na.rm = TRUE, names = FALSE, type = 7)
  data.table(feature = v, low = q[1], high = q[2])
}))

d3r_winsor_apply <- function(dt, bounds) {
  out <- copy(dt)
  for (i in seq_len(nrow(bounds))) {
    v <- bounds$feature[i]
    out[[v]] <- pmin(pmax(out[[v]], bounds$low[i]), bounds$high[i])
  }
  out
}

d3r_fit_discovery <- function(dt) {
  if (anyNA(dt[, ..D3R_FEATURES])) stop("Discovery input contains missing features.")
  bounds <- d3r_winsor_fit(dt)
  w <- d3r_winsor_apply(dt, bounds)
  X <- as.matrix(w[, ..D3R_FEATURES])
  center <- colMeans(X)
  scalev <- apply(X, 2, sd)
  if (any(!is.finite(scalev) | scalev <= 0)) stop("Discovery scale failure.")
  Xz <- sweep(sweep(X, 2, center, "-"), 2, scalev, "/")
  km <- kmeans(Xz, centers = 2, nstart = D3R_KMEANS_NSTART,
               iter.max = D3R_KMEANS_ITER_MAX, algorithm = "Lloyd")
  means <- tapply(dt$severity_z, km$cluster, mean)
  high <- as.integer(names(which.max(means)))
  labels <- ifelse(km$cluster == high, "C1", "C2")
  other <- setdiff(1:2, high)
  centroids <- rbind(C1 = km$centers[high, ], C2 = km$centers[other, ])
  list(labels = labels, bounds = bounds, center = center, scale = scalev,
       centroids = centroids, Xz = Xz)
}

d3r_apply_discovery_complete <- function(object, dt) {
  w <- d3r_winsor_apply(dt, object$bounds)
  X <- as.matrix(w[, ..D3R_FEATURES])
  Xz <- sweep(sweep(X, 2, object$center, "-"), 2, object$scale, "/")
  d <- sapply(1:2, function(i) rowSums((sweep(Xz, 2, object$centroids[i, ], "-"))^2))
  rownames(object$centroids)[max.col(-d, ties.method = "first")]
}

d3r_derive_recipe <- function(dt) {
  rbindlist(lapply(D3R_FEATURES, function(v) {
    x <- as.numeric(dt[[v]])
    q <- quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE, type = 7)
    xw <- pmin(pmax(x, q[1]), q[2])
    data.table(
      feature = v, winsor_p01 = q[1], winsor_p99 = q[2],
      imputation_median = median(xw, na.rm = TRUE),
      reference_mean = mean(xw, na.rm = TRUE),
      reference_sd = sd(xw, na.rm = TRUE)
    )
  }))
}

d3r_apply_recipe <- function(dt, recipe) {
  X <- matrix(NA_real_, nrow(dt), length(D3R_FEATURES), dimnames = list(NULL, D3R_FEATURES))
  missing_n <- clipped_n <- integer(nrow(dt))
  for (k in seq_along(D3R_FEATURES)) {
    v <- D3R_FEATURES[k]; r <- recipe[k]
    x <- as.numeric(dt[[v]]); miss <- !is.finite(x)
    missing_n <- missing_n + miss
    clipped <- !miss & (x < r$winsor_p01 | x > r$winsor_p99)
    clipped_n <- clipped_n + clipped
    x <- pmin(pmax(x, r$winsor_p01), r$winsor_p99)
    x[miss] <- r$imputation_median
    X[, k] <- (x - r$reference_mean) / r$reference_sd
  }
  if (any(!is.finite(X))) stop("Regenerated preprocessing produced non-finite values.")
  list(X = X, missing_n = missing_n, clipped_n = clipped_n)
}

d3r_build_generator <- function(raw, locked_labels) {
  recipe <- d3r_derive_recipe(raw)
  recipe_numeric <- as.matrix(
    recipe[, setdiff(names(recipe), "feature"), with = FALSE]
  )
  if (any(!is.finite(recipe_numeric)) || any(recipe$reference_sd <= 0)) {
    stop("Invalid regenerated recipe.")
  }
  prep <- d3r_apply_recipe(raw, recipe)
  centroids <- t(vapply(c("C1", "C2"), function(k) {
    colMeans(prep$X[locked_labels == k, , drop = FALSE])
  }, numeric(length(D3R_FEATURES))))
  rownames(centroids) <- c("C1", "C2"); colnames(centroids) <- D3R_FEATURES
  list(recipe = recipe, centroids = centroids,
       feature_names = D3R_FEATURES, distance = "squared_euclidean")
}

d3r_predict_generator <- function(object, dt) {
  prep <- d3r_apply_recipe(dt, object$recipe)
  d <- sapply(1:2, function(i) rowSums((sweep(prep$X, 2, object$centroids[i, ], "-"))^2))
  idx <- max.col(-d, ties.method = "first")
  data.table(
    assigned = rownames(object$centroids)[idx], distance_c1 = d[, 1],
    distance_c2 = d[, 2], absolute_margin = abs(d[, 2] - d[, 1]),
    missing_n = prep$missing_n, clipped_n = prep$clipped_n
  )
}

d3r_ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  choose2 <- function(z) z * (z - 1) / 2
  index <- sum(choose2(tab)); a <- rowSums(tab); b <- colSums(tab)
  expected <- sum(choose2(a)) * sum(choose2(b)) / choose2(n)
  denom <- 0.5 * (sum(choose2(a)) + sum(choose2(b))) - expected
  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(NA_real_)
  (index - expected) / denom
}

d3r_fidelity <- function(reference, assigned) {
  data.table(
    n = length(reference), agreement = mean(reference == assigned),
    ari = d3r_ari(reference, assigned),
    reference_c1 = mean(reference == "C1"), assigned_c1 = mean(assigned == "C1"),
    prevalence_difference = mean(assigned == "C1") - mean(reference == "C1")
  )
}

d3r_crossfit_fidelity <- function(raw, locked, seed) {
  set.seed(seed)
  fold <- integer(nrow(raw))
  for (k in c("C1", "C2")) {
    ii <- which(locked == k)
    fold[ii] <- sample(rep_len(seq_len(D3R_CV_FOLDS), length(ii)))
  }
  pred <- rep(NA_character_, nrow(raw))
  for (f in seq_len(D3R_CV_FOLDS)) {
    train <- fold != f; test <- fold == f
    gen <- d3r_build_generator(raw[train], locked[train])
    pred[test] <- d3r_predict_generator(gen, raw[test])$assigned
  }
  if (anyNA(pred)) stop("Cross-fit generator predictions incomplete.")
  d3r_fidelity(locked, pred)
}

d3r_fit_base_k2 <- function(dt, labels) {
  d <- copy(dt)
  d[, generated_k2 := factor(labels, levels = c("C2", "C1"))]
  glm(mortality_30d ~ age + sex + sofa + aki_stage + generated_k2,
      data = d, family = binomial())
}

d3r_auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); y <- as.integer(y[ok]); p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

d3r_calibration <- function(y, p) {
  p <- d3r_clamp_prob(p); lp <- qlogis(p)
  slope_fit <- glm(y ~ lp, family = binomial())
  int_fit <- glm(y ~ 1, offset = lp, family = binomial())
  c(slope = unname(coef(slope_fit)[2]), intercept = unname(coef(int_fit)[1]))
}

d3r_run_repeat <- function(scenario, repeat_id) {
  seed <- d3r_seed(scenario, repeat_id)
  stage <- "generation"
  tryCatch({
    deriv <- d3r_generate_dataset(scenario, "derivation", D3R_N_TRAIN, repeat_id, seed + 1L)
    external <- d3r_generate_dataset(scenario, "external", D3R_N_EXTERNAL, repeat_id, seed + 2L)

    stage <- "mice"
    imps <- d3r_mice(deriv$raw, seed + 101L)
    stage <- "discovery"
    discoveries <- lapply(imps, d3r_fit_discovery)
    primary <- discoveries[[D3R_PRIMARY_IMPUTATION]]
    locked <- primary$labels
    mice_ari <- vapply(discoveries, function(z) d3r_ari(locked, z$labels), numeric(1))

    stage <- "generator_construction"
    generator <- d3r_build_generator(deriv$raw, locked)
    dev_pred <- d3r_predict_generator(generator, deriv$raw)
    internal <- d3r_fidelity(locked, dev_pred$assigned)
    crossfit <- d3r_crossfit_fidelity(deriv$raw, locked, seed + 202L)

    stage <- "external_assignment"
    ext_pred <- d3r_predict_generator(generator, external$raw)
    complete_reference <- d3r_apply_discovery_complete(primary, external$complete)
    ext_fidelity <- d3r_fidelity(complete_reference, ext_pred$assigned)

    # Determinism and batch invariance use the same frozen object.
    ncheck <- min(100L, nrow(external$raw))
    check <- external$raw[seq_len(ncheck)]
    p1 <- d3r_predict_generator(generator, check)
    p2 <- d3r_predict_generator(generator, check)
    perm <- rev(seq_len(ncheck))
    p3 <- d3r_predict_generator(generator, check[perm])[rev(seq_len(ncheck))]
    invariant <- identical(p1$assigned, p2$assigned) && identical(p1$assigned, p3$assigned) &&
      isTRUE(all.equal(p1$distance_c1, p2$distance_c1, tolerance = 0))

    stage <- "outcome_model"
    model <- d3r_fit_base_k2(deriv$raw, dev_pred$assigned)
    ext_model <- copy(external$raw)
    ext_model[, generated_k2 := factor(ext_pred$assigned, levels = c("C2", "C1"))]
    prob <- as.numeric(predict(model, newdata = ext_model, type = "response"))
    cal <- d3r_calibration(ext_model$mortality_30d, prob)
    auc <- d3r_auc(ext_model$mortality_30d, prob)
    brier <- mean((ext_model$mortality_30d - prob)^2)

    prevalence_drift <- mean(ext_pred$assigned == "C1") - mean(dev_pred$assigned == "C1")
    prevalence_alert <- abs(prevalence_drift) > D3R_PREVALENCE_GATE
    assignment_alert <- is.finite(ext_fidelity$ari) && ext_fidelity$ari < D3R_ASSIGNMENT_ARI_GATE
    calibration_alert <- abs(cal[["intercept"]]) > D3R_CALIBRATION_INTERCEPT_GATE ||
      cal[["slope"]] < D3R_CALIBRATION_SLOPE_RANGE[1] ||
      cal[["slope"]] > D3R_CALIBRATION_SLOPE_RANGE[2]
    combined_alert <- prevalence_alert || assignment_alert || calibration_alert

    data.table(
      scenario, repeat_id, expected_alert = unname(D3R_EXPECTED_ALERT[[scenario]]),
      failure = FALSE, failure_stage = NA_character_, failure_reason = NA_character_,
      n_train = nrow(deriv$raw), n_external = nrow(external$raw),
      internal_agreement = internal$agreement, internal_ari = internal$ari,
      internal_prevalence_difference = internal$prevalence_difference,
      crossfit_agreement = crossfit$agreement, crossfit_ari = crossfit$ari,
      mice_ari_mean = mean(mice_ari), mice_ari_min = min(mice_ari),
      external_reference_agreement = ext_fidelity$agreement,
      external_reference_ari = ext_fidelity$ari,
      derivation_c1_prevalence = mean(dev_pred$assigned == "C1"),
      external_c1_prevalence = mean(ext_pred$assigned == "C1"),
      prevalence_drift = prevalence_drift,
      calibration_slope = cal[["slope"]], calibration_intercept = cal[["intercept"]],
      external_auc = auc, external_brier = brier,
      external_assignment_complete = all(!is.na(ext_pred$assigned)),
      batch_invariance_pass = invariant,
      external_missing_feature_mean = mean(ext_pred$missing_n),
      external_clipped_feature_mean = mean(ext_pred$clipped_n),
      external_margin_median = median(ext_pred$absolute_margin),
      prevalence_alert, assignment_alert, calibration_alert, combined_alert,
      correct_classification = combined_alert == unname(D3R_EXPECTED_ALERT[[scenario]])
    )
  }, error = function(e) {
    data.table(
      scenario, repeat_id, expected_alert = unname(D3R_EXPECTED_ALERT[[scenario]]),
      failure = TRUE, failure_stage = stage, failure_reason = conditionMessage(e)
    )
  })
}

d3r_wilson <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, ctr - half), high = min(1, ctr + half))
}
