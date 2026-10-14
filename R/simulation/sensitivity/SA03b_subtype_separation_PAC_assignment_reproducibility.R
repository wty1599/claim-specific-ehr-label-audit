
################################################################################
## COMMON STANDALONE SIMULATION UTILITIES
## Generated for the EHR audit / EHR-derived subphenotype four-domain framework.
##
## Design principles:
##   - no source() dependency on older simulation files;
##   - MICE m = 5, completed dataset 1 defines the primary clustering result;
##   - derivation preprocessing and centroids are frozen for held-out/external data;
##   - heavy discreteness diagnostics are run on completed dataset 1 only;
##   - every diagnostic returns a named list/data.table contract;
##   - task-level RDS checkpoints support restart;
##   - all-task failure stops before empty-table post-processing.
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(warn = 1)

## =============================================================================
## Package checks
## =============================================================================

required_pkgs <- c(
  "data.table", "mice", "mclust", "diptest", "cluster",
  "glmnet", "pROC", "future", "future.apply", "ggplot2"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(mclust)
  library(diptest)
  library(cluster)
  library(glmnet)
  library(pROC)
  library(future)
  library(future.apply)
  library(ggplot2)
})

## =============================================================================
## Generic utilities
## =============================================================================

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  message(msg)
  if (exists("LOG_FILE", inherits = TRUE)) {
    cat(msg, "\n", file = get("LOG_FILE", inherits = TRUE), append = TRUE)
  }
  invisible(msg)
}

safe_num <- function(x, default = NA_real_) {
  if (length(x) == 0L || all(is.na(x))) return(default)
  as.numeric(x[1L])
}

safe_int <- function(x, default = NA_integer_) {
  if (length(x) == 0L || all(is.na(x))) return(default)
  as.integer(x[1L])
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  if (sum(is.finite(x)) < 2L) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, p) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

safe_rbindlist <- function(x) {
  x <- Filter(
    function(z) {
      !is.null(z) &&
        (is.data.frame(z) || data.table::is.data.table(z)) &&
        nrow(z) > 0L
    },
    x
  )
  if (!length(x)) return(data.table())
  x <- lapply(x, function(z) {
    z <- as.data.table(z)
    if (anyDuplicated(names(z))) setnames(z, make.unique(names(z), sep = "_dup"))
    z
  })
  rbindlist(x, fill = TRUE, use.names = TRUE)
}


write_csv_safe <- function(x, path) {
  x <- tryCatch(as.data.table(x), error = function(e) data.table())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (ncol(x) == 0L) {
    writeLines("", con = path, useBytes = TRUE)
  } else {
    data.table::fwrite(x, path)
  }
  invisible(path)
}

mcse_mean <- function(x) {
  x <- as.numeric(x)
  n <- sum(is.finite(x))
  if (n < 2L) return(NA_real_)
  stats::sd(x, na.rm = TRUE) / sqrt(n)
}

mcse_prop <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 1L) return(NA_real_)
  p <- mean(x)
  sqrt(p * (1 - p) / n)
}

summarise_numeric <- function(dt, group_cols, value_col, metric_name = value_col) {
  if (!nrow(dt) || !value_col %in% names(dt)) return(data.table())
  dt[, {
    z <- as.numeric(get(value_col))
    z <- z[is.finite(z)]
    n <- length(z)
    list(
      metric = metric_name,
      n = n,
      mean = if (n) mean(z) else NA_real_,
      sd = if (n > 1L) sd(z) else NA_real_,
      mcse = if (n > 1L) sd(z) / sqrt(n) else NA_real_,
      median = if (n) median(z) else NA_real_,
      empirical_p025 = if (n) safe_quantile(z, 0.025) else NA_real_,
      empirical_p975 = if (n) safe_quantile(z, 0.975) else NA_real_
    )
  }, by = group_cols]
}

summarise_binary <- function(dt, group_cols, value_col, metric_name = value_col) {
  if (!nrow(dt) || !value_col %in% names(dt)) return(data.table())
  dt[, {
    z <- as.numeric(get(value_col))
    z <- z[is.finite(z)]
    n <- length(z)
    p <- if (n) mean(z) else NA_real_
    list(
      metric = metric_name,
      n = n,
      proportion = p,
      mcse = if (n) sqrt(p * (1 - p) / n) else NA_real_,
      exact_low = if (n) stats::binom.test(sum(z), n)$conf.int[1L] else NA_real_,
      exact_high = if (n) stats::binom.test(sum(z), n)$conf.int[2L] else NA_real_
    )
  }, by = group_cols]
}

write_session_info <- function(path) {
  zz <- file(path, open = "wt", encoding = "UTF-8")
  ## Revert the sink BEFORE closing the connection. Closing an active sink
  ## target first leaves sink() pointing at an invalid connection, which throws
  ## "invalid connection" on exit.
  on.exit({
    sink()
    close(zz)
  }, add = TRUE)
  sink(zz)
  print(sessionInfo())
}

## =============================================================================
## Strict diagnostic contracts
## =============================================================================

new_diag <- function(
  ok = FALSE,
  stage = NA_character_,
  error = NA_character_,
  table = data.table(),
  summary = list()
) {
  structure(
    list(
      ok = isTRUE(ok),
      stage = as.character(stage)[1L],
      error = as.character(error)[1L],
      table = as.data.table(table),
      summary = summary
    ),
    class = c("ehr_audit_diag", "list")
  )
}

validate_diag <- function(x, stage) {
  required <- c("ok", "stage", "error", "table", "summary")
  if (!is.list(x) || !all(required %in% names(x))) {
    return(new_diag(
      ok = FALSE,
      stage = stage,
      error = paste0(
        "Invalid diagnostic return: class=",
        paste(class(x), collapse = "/"),
        "; names=",
        paste(names(x), collapse = ",")
      )
    ))
  }
  x$table <- tryCatch(as.data.table(x$table), error = function(e) data.table())
  x
}

safe_stage <- function(stage, expr) {
  tryCatch(
    {
      out <- force(expr)
      validate_diag(out, stage)
    },
    error = function(e) {
      new_diag(
        ok = FALSE,
        stage = stage,
        error = conditionMessage(e)
      )
    }
  )
}

## =============================================================================
## Feature metadata and data-generating mechanisms
## =============================================================================

feature_names <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min",
  "abs_lymphocytes_min", "lactate_max", "ph_min", "pco2_max",
  "calcium_min", "calcium_max", "ptt_max", "inr_max",
  "temperature_min", "temperature_max", "urine_output_24h_ml",
  "glucose_max", "aniongap_max", "potassium_min", "potassium_max",
  "hemoglobin_min", "sodium_min", "sodium_max", "wbc_max",
  "platelets_min", "bicarbonate_min", "chloride_min", "chloride_max",
  "bun_max", "creatinine_max", "resp_rate_max", "gcs_min",
  "spo2_min", "heart_rate_max", "mbp_min"
)

baseline_vars <- c("age", "sex", "sofa", "aki_stage")
outcome_vars <- c("make30", "mortality_30d")
imputation_vars <- c(feature_names, baseline_vars)

severity_loading <- c(
  bilirubin_total_max = 0.45,
  alt_max = 0.25,
  pao2fio2ratio_min = -0.65,
  abs_lymphocytes_min = -0.35,
  lactate_max = 0.85,
  ph_min = -0.75,
  pco2_max = 0.30,
  calcium_min = -0.20,
  calcium_max = 0.10,
  ptt_max = 0.40,
  inr_max = 0.45,
  temperature_min = -0.15,
  temperature_max = 0.20,
  urine_output_24h_ml = -0.95,
  glucose_max = 0.25,
  aniongap_max = 0.80,
  potassium_min = -0.10,
  potassium_max = 0.35,
  hemoglobin_min = -0.25,
  sodium_min = -0.10,
  sodium_max = 0.10,
  wbc_max = 0.35,
  platelets_min = -0.40,
  bicarbonate_min = -0.85,
  chloride_min = -0.10,
  chloride_max = 0.15,
  bun_max = 0.90,
  creatinine_max = 0.95,
  resp_rate_max = 0.50,
  gcs_min = -0.60,
  spo2_min = -0.35,
  heart_rate_max = 0.55,
  mbp_min = -0.70
)

stopifnot(
  length(feature_names) == 33L,
  identical(names(severity_loading), feature_names)
)

make_subtype_patterns <- function() {
  p <- length(feature_names)
  out <- matrix(0, nrow = 3L, ncol = p)
  colnames(out) <- feature_names
  rownames(out) <- paste0("Subtype", 1:3)

  out[1L, c(
    "lactate_max", "aniongap_max", "heart_rate_max", "resp_rate_max",
    "wbc_max", "temperature_max"
  )] <- 1
  out[1L, c(
    "mbp_min", "bicarbonate_min", "ph_min", "platelets_min"
  )] <- -1

  out[2L, c(
    "creatinine_max", "bun_max", "potassium_max", "inr_max", "ptt_max"
  )] <- 1
  out[2L, c(
    "urine_output_24h_ml", "calcium_min", "hemoglobin_min"
  )] <- -1

  out[3L, c(
    "bilirubin_total_max", "alt_max", "pco2_max", "chloride_max"
  )] <- 1
  out[3L, c(
    "pao2fio2ratio_min", "spo2_min", "gcs_min",
    "abs_lymphocytes_min"
  )] <- -1

  out <- out - matrix(colMeans(out), nrow = 3L, ncol = p, byrow = TRUE)
  row_norm <- sqrt(rowMeans(out^2))
  out <- out / row_norm
  out
}

subtype_patterns <- make_subtype_patterns()

generate_factor_loadings <- function(p, q = 5L, seed = 91001L) {
  set.seed(seed)
  L <- matrix(rnorm(p * q, sd = 0.25), nrow = p, ncol = q)
  L
}

factor_loadings <- generate_factor_loadings(length(feature_names))

generate_cohort <- function(
  n,
  delta = 0,
  subtype_prob = c(1/3, 1/3, 1/3),
  severity_mean = 0,
  external_drift = FALSE,
  seed
) {
  set.seed(seed)

  subtype <- sample.int(3L, n, replace = TRUE, prob = subtype_prob)
  severity <- rnorm(n, mean = severity_mean, sd = 1)

  F <- matrix(rnorm(n * ncol(factor_loadings)), nrow = n)
  correlated_noise <- F %*% t(factor_loadings)
  independent_noise <- matrix(
    rnorm(n * length(feature_names), sd = 0.75),
    nrow = n
  )

  X <- severity %o% severity_loading +
    correlated_noise +
    independent_noise +
    delta * subtype_patterns[subtype, , drop = FALSE]

  if (isTRUE(external_drift)) {
    shift <- setNames(rep(0, length(feature_names)), feature_names)
    shift[c(
      "lactate_max", "creatinine_max", "bun_max",
      "resp_rate_max", "heart_rate_max"
    )] <- c(0.35, 0.30, 0.30, 0.20, 0.20)
    shift[c(
      "mbp_min", "pao2fio2ratio_min", "urine_output_24h_ml"
    )] <- c(-0.30, -0.35, -0.25)

    X <- sweep(X, 2L, shift, FUN = "+")
    scale_mult <- setNames(rep(1, length(feature_names)), feature_names)
    scale_mult[c(
      "lactate_max", "creatinine_max", "pao2fio2ratio_min",
      "urine_output_24h_ml"
    )] <- c(1.20, 1.15, 1.20, 1.15)
    X <- sweep(X, 2L, scale_mult, FUN = "*")
  }

  age <- pmin(95, pmax(18, round(65 + 11 * rnorm(n) + 2.0 * severity)))
  sex <- rbinom(n, 1, plogis(-0.05 + 0.08 * severity))
  sofa <- pmin(24, pmax(0, round(6 + 2.8 * severity + rnorm(n, sd = 2))))
  aki_stage <- pmin(
    3L,
    pmax(
      0L,
      as.integer(cut(
        severity + rnorm(n, sd = 0.8),
        breaks = c(-Inf, -0.6, 0.2, 1.0, Inf),
        labels = FALSE
      )) - 1L
    )
  )

  lp_mort <- -2.35 +
    0.70 * severity +
    0.018 * (age - 65) +
    0.08 * sex +
    0.065 * sofa +
    0.16 * aki_stage

  lp_make <- -1.30 +
    0.85 * severity +
    0.012 * (age - 65) +
    0.055 * sofa +
    0.28 * aki_stage

  mortality_30d <- rbinom(n, 1, plogis(lp_mort))
  make30 <- rbinom(n, 1, plogis(lp_make))

  dt <- as.data.table(X)
  setnames(dt, feature_names)
  dt[, patient_id := seq_len(.N)]
  dt[, latent_subtype := subtype]
  dt[, severity_true := severity]
  dt[, age := age]
  dt[, sex := sex]
  dt[, sofa := sofa]
  dt[, aki_stage := aki_stage]
  dt[, mortality_30d := mortality_30d]
  dt[, make30 := make30]

  setcolorder(
    dt,
    c(
      "patient_id", "latent_subtype", "severity_true",
      baseline_vars, feature_names, outcome_vars
    )
  )
  dt
}

inject_mar_missingness <- function(
  dt,
  seed,
  base_low = 0.05,
  base_high = 0.28,
  external_extra = 0
) {
  set.seed(seed)
  out <- copy(dt)
  p <- length(feature_names)
  base_rate <- seq(base_low, base_high, length.out = p)
  base_rate <- base_rate[sample.int(p)]

  sev <- as.numeric(scale(out$severity_true))
  age_z <- as.numeric(scale(out$age))

  for (j in seq_along(feature_names)) {
    v <- feature_names[j]
    xz <- as.numeric(scale(out[[v]]))
    xz[!is.finite(xz)] <- 0
    lin <- qlogis(pmin(0.80, base_rate[j] + external_extra)) +
      0.30 * sev +
      0.12 * age_z +
      0.10 * xz
    prob <- pmin(0.85, pmax(0.005, plogis(lin)))
    miss <- rbinom(nrow(out), 1, prob) == 1L
    set(out, which(miss), v, NA_real_)
  }

  for (v in c("sofa", "aki_stage")) {
    prob <- pmin(
      0.20,
      pmax(
        0.01,
        0.03 + external_extra / 2 + 0.03 * plogis(sev)
      )
    )
    miss <- rbinom(nrow(out), 1, prob) == 1L
    set(out, which(miss), v, NA_real_)
  }

  out
}

split_derivation <- function(dt, train_prop = 0.70, seed) {
  set.seed(seed)
  idx <- sample.int(nrow(dt))
  n_train <- floor(train_prop * nrow(dt))
  list(
    analysis_train = copy(dt[idx[seq_len(n_train)]]),
    internal_test = copy(dt[idx[(n_train + 1L):length(idx)]])
  )
}

mice_complete_sets <- function(dt, m, maxit, seed) {
  imp_cols <- intersect(imputation_vars, names(dt))
  imp_data <- as.data.frame(dt[, ..imp_cols])

  method <- rep("pmm", length(imp_cols))
  names(method) <- imp_cols
  no_missing <- vapply(imp_data, function(x) !anyNA(x), logical(1))
  method[no_missing] <- ""

  predictor_matrix <- mice::make.predictorMatrix(imp_data)
  diag(predictor_matrix) <- 0

  fit <- mice::mice(
    imp_data,
    m = m,
    maxit = maxit,
    method = method,
    predictorMatrix = predictor_matrix,
    seed = seed,
    printFlag = FALSE
  )

  lapply(seq_len(m), function(i) {
    completed <- mice::complete(fit, action = i)
    out <- copy(dt)
    for (v in imp_cols) {
      set(out, j = v, value = completed[[v]])
    }
    out[, imputation_id := i]
    out
  })
}

## =============================================================================
## Preprocessing, clustering, matching, and assignment reproducibility
## =============================================================================

fit_preprocess <- function(dt, features, probs = c(0.01, 0.99)) {
  bounds <- lapply(features, function(v) {
    as.numeric(quantile(dt[[v]], probs = probs, na.rm = TRUE, names = FALSE))
  })
  names(bounds) <- features

  Xw <- sapply(features, function(v) {
    b <- bounds[[v]]
    pmin(pmax(as.numeric(dt[[v]]), b[1L]), b[2L])
  })
  Xw <- as.matrix(Xw)
  colnames(Xw) <- features

  center <- colMeans(Xw)
  scale_value <- apply(Xw, 2L, sd)
  scale_value[!is.finite(scale_value) | scale_value < 1e-8] <- 1

  list(
    bounds = bounds,
    center = center,
    scale = scale_value,
    features = features
  )
}

apply_preprocess <- function(dt, prep) {
  Xw <- sapply(prep$features, function(v) {
    b <- prep$bounds[[v]]
    pmin(pmax(as.numeric(dt[[v]]), b[1L]), b[2L])
  })
  Xw <- as.matrix(Xw)
  colnames(Xw) <- prep$features
  Xz <- sweep(Xw, 2L, prep$center, FUN = "-")
  Xz <- sweep(Xz, 2L, prep$scale, FUN = "/")
  Xz[!is.finite(Xz)] <- 0
  Xz
}

safe_kmeans <- function(X, k, seed, nstart = 50L, iter.max = 500L) {
  X <- as.matrix(X)
  if (nrow(X) <= k || ncol(X) < 1L || any(!is.finite(X))) {
    return(list(
      ok = FALSE,
      labels = rep(NA_integer_, nrow(X)),
      centers = matrix(NA_real_, nrow = k, ncol = ncol(X)),
      tot_withinss = NA_real_,
      algorithm_used = NA_character_,
      warning = NA_character_,
      error = "Invalid matrix for k-means"
    ))
  }

  warning_text <- character()

  set.seed(seed)
  fit <- tryCatch(
    withCallingHandlers(
      stats::kmeans(
        X,
        centers = k,
        nstart = nstart,
        iter.max = iter.max,
        algorithm = "Hartigan-Wong"
      ),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  algorithm_used <- "Hartigan-Wong"
  quick_transfer_warning <- any(
    grepl(
      "Quick-TRANSfer|Quick-Transfer|maximum",
      warning_text,
      ignore.case = TRUE
    )
  )

  if (is.null(fit) || quick_transfer_warning) {
    warning_text_lloyd <- character()
    set.seed(seed + 17L)
    fit_lloyd <- tryCatch(
      withCallingHandlers(
        stats::kmeans(
          X,
          centers = k,
          nstart = nstart,
          iter.max = iter.max,
          algorithm = "Lloyd"
        ),
        warning = function(w) {
          warning_text_lloyd <<- c(
            warning_text_lloyd,
            conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_lloyd)) {
      fit <- fit_lloyd
      algorithm_used <- "Lloyd_fallback"
      warning_text <- c(warning_text, warning_text_lloyd)
    }
  }

  if (is.null(fit) || !is.list(fit) ||
      !all(c("cluster", "centers", "tot.withinss") %in% names(fit))) {
    return(list(
      ok = FALSE,
      labels = rep(NA_integer_, nrow(X)),
      centers = matrix(NA_real_, nrow = k, ncol = ncol(X)),
      tot_withinss = NA_real_,
      algorithm_used = algorithm_used,
      warning = paste(unique(warning_text), collapse = " | "),
      error = "k-means failed or returned an invalid object"
    ))
  }

  list(
    ok = TRUE,
    labels = as.integer(fit$cluster),
    centers = as.matrix(fit$centers),
    tot_withinss = as.numeric(fit$tot.withinss),
    algorithm_used = algorithm_used,
    warning = if (length(warning_text)) {
      paste(unique(warning_text), collapse = " | ")
    } else {
      NA_character_
    },
    error = NA_character_
  )
}

relabel_k2_by_burden <- function(fit) {
  if (!is.list(fit) || !isTRUE(fit$ok) || nrow(fit$centers) != 2L) return(fit)
  loading <- severity_loading[colnames(fit$centers)]
  score <- as.numeric(fit$centers %*% loading)
  high_old <- which.max(score)
  map <- integer(2L)
  map[high_old] <- 1L
  map[setdiff(1:2, high_old)] <- 2L
  fit$labels <- map[fit$labels]
  fit$centers <- fit$centers[c(high_old, setdiff(1:2, high_old)), , drop = FALSE]
  fit
}

assign_centroids <- function(X, centers) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)
  d <- sapply(seq_len(nrow(centers)), function(k) {
    rowSums((X - matrix(
      centers[k, ],
      nrow = nrow(X),
      ncol = ncol(X),
      byrow = TRUE
    ))^2)
  })
  if (is.null(dim(d))) d <- matrix(d, ncol = 1L)
  max.col(-d, ties.method = "first")
}

adjusted_rand <- function(a, b) {
  if (length(a) != length(b) || !length(a) ||
      anyNA(a) || anyNA(b)) return(NA_real_)
  as.numeric(mclust::adjustedRandIndex(a, b))
}


truth_structure_metrics <- function(labels, truth) {
  if (length(labels) != length(truth) || anyNA(labels) || anyNA(truth)) {
    return(list(
      purity = NA_real_,
      weighted_cluster_entropy = NA_real_,
      normalized_mutual_information = NA_real_
    ))
  }

  tab <- table(labels, truth)
  n <- sum(tab)
  purity <- sum(apply(tab, 1L, max)) / n

  row_prob <- rowSums(tab) / n
  entropy_by_cluster <- vapply(seq_len(nrow(tab)), function(i) {
    p <- tab[i, ] / sum(tab[i, ])
    p <- p[p > 0]
    -sum(p * log(p))
  }, numeric(1))
  weighted_entropy <- sum(row_prob * entropy_by_cluster)

  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  mi <- 0
  for (i in seq_len(nrow(pxy))) {
    for (j in seq_len(ncol(pxy))) {
      if (pxy[i, j] > 0) {
        mi <- mi + pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  nmi <- if (hx > 0 && hy > 0) mi / sqrt(hx * hy) else NA_real_

  list(
    purity = as.numeric(purity),
    weighted_cluster_entropy = as.numeric(weighted_entropy),
    normalized_mutual_information = as.numeric(nmi)
  )
}

truth_composition_table <- function(
  labels,
  truth,
  task_id,
  repeat_id,
  setting_name,
  imputation_id,
  algorithm,
  k,
  sample_role = "derivation"
) {
  tab <- as.data.table(
    as.data.frame(
      table(
        cluster = factor(labels, levels = seq_len(k)),
        latent_subtype = factor(truth, levels = 1:3)
      )
    )
  )
  setnames(tab, "Freq", "n")
  tab[, cluster_total := sum(n), by = cluster]
  tab[, subtype_total := sum(n), by = latent_subtype]
  tab[, proportion_within_cluster := n / cluster_total]
  tab[, proportion_within_subtype := n / subtype_total]
  tab[, `:=`(
    task_id = task_id,
    repeat_id = repeat_id,
    setting = setting_name,
    imputation_id = imputation_id,
    algorithm = algorithm,
    k = k,
    sample_role = sample_role
  )]
  setcolorder(
    tab,
    c(
      "task_id", "repeat_id", "setting", "imputation_id",
      "algorithm", "k", "sample_role", "cluster", "latent_subtype", "n",
      "proportion_within_cluster", "proportion_within_subtype"
    )
  )
  tab
}

all_permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
  do.call(
    rbind,
    lapply(v, function(x) {
      rest <- all_permutations(v[v != x])
      cbind(x, rest)
    })
  )
}

best_label_match <- function(reference, candidate, k) {
  if (length(reference) != length(candidate) ||
      anyNA(reference) || anyNA(candidate)) {
    return(list(
      agreement = NA_real_,
      mapped = rep(NA_integer_, length(candidate)),
      map = rep(NA_integer_, k)
    ))
  }

  perms <- all_permutations(seq_len(k))
  agreement <- apply(perms, 1L, function(p) {
    mean(reference == p[candidate])
  })
  best <- which.max(agreement)
  list(
    agreement = agreement[best],
    mapped = perms[best, candidate],
    map = as.integer(perms[best, ])
  )
}

cluster_prevalence <- function(labels, k) {
  tab <- table(factor(labels, levels = seq_len(k)))
  as.numeric(tab) / sum(tab)
}

centroids_from_labels <- function(X, labels, k) {
  out <- matrix(NA_real_, nrow = k, ncol = ncol(X))
  colnames(out) <- colnames(X)
  for (j in seq_len(k)) {
    idx <- which(labels == j)
    if (length(idx)) out[j, ] <- colMeans(X[idx, , drop = FALSE])
  }
  out
}

profile_correlation_after_match <- function(
  derivation_centers,
  external_centers,
  map_external_to_derivation
) {
  k <- nrow(derivation_centers)
  mapped_external <- matrix(
    NA_real_,
    nrow = k,
    ncol = ncol(external_centers)
  )
  for (old in seq_len(k)) {
    new <- map_external_to_derivation[old]
    mapped_external[new, ] <- external_centers[old, ]
  }
  cors <- vapply(seq_len(k), function(j) {
    suppressWarnings(cor(
      derivation_centers[j, ],
      mapped_external[j, ],
      use = "pairwise.complete.obs"
    ))
  }, numeric(1))
  safe_mean(cors)
}

assignment_reproducibility_metrics <- function(
  frozen_labels,
  recluster_labels,
  derivation_labels,
  X_external,
  derivation_centers,
  k
) {
  match <- best_label_match(frozen_labels, recluster_labels, k)
  mapped_recluster <- match$mapped

  prev_der <- cluster_prevalence(derivation_labels, k)
  prev_ext <- cluster_prevalence(mapped_recluster, k)
  ext_centers <- centroids_from_labels(X_external, recluster_labels, k)

  data.table(
    k = k,
    frozen_vs_reclustered_ari = adjusted_rand(
      frozen_labels,
      recluster_labels
    ),
    matched_agreement = match$agreement,
    mean_abs_prevalence_difference = mean(abs(prev_ext - prev_der)),
    max_abs_prevalence_difference = max(abs(prev_ext - prev_der)),
    centroid_profile_correlation = profile_correlation_after_match(
      derivation_centers,
      ext_centers,
      match$map
    )
  )
}

outcome_separation <- function(dt, labels, outcome = "mortality_30d") {
  if (length(labels) != nrow(dt) || anyNA(labels)) return(NA_real_)
  y <- as.numeric(dt[[outcome]])
  safe_mean(y[labels == 1L]) - safe_mean(y[labels == 2L])
}

## =============================================================================
## Consensus PAC and consensus k-means
## =============================================================================

sample_pair_index <- function(n, n_pairs, seed) {
  set.seed(seed)
  i <- sample.int(n, n_pairs, replace = TRUE)
  j <- sample.int(n, n_pairs, replace = TRUE)
  same <- i == j
  while (any(same)) {
    j[same] <- sample.int(n, sum(same), replace = TRUE)
    same <- i == j
  }
  lo <- pmin(i, j)
  hi <- pmax(i, j)
  data.table(i = lo, j = hi)
}

consensus_pac_sampled_pairs <- function(
  X,
  k_values = 2:4,
  reps = 50L,
  subsample_fraction = 0.80,
  subsample_n = 800L,
  n_pairs = 50000L,
  lower = 0.10,
  upper = 0.90,
  seed = 1001L
) {
  tryCatch({
    X <- as.matrix(X)
    n_sub <- min(nrow(X), as.integer(subsample_n))
    if (n_sub < max(k_values) + 10L) {
      stop("Insufficient rows for PAC")
    }

    set.seed(seed)
    idx_sub <- sort(sample.int(nrow(X), n_sub))
    Xs <- X[idx_sub, , drop = FALSE]
    pairs <- sample_pair_index(n_sub, n_pairs, seed + 11L)

    rows <- lapply(k_values, function(k) {
      co_observed <- integer(nrow(pairs))
      co_clustered <- integer(nrow(pairs))

      for (b in seq_len(reps)) {
        set.seed(seed + 10000L * k + b)
        idx <- sort(sample.int(
          n_sub,
          max(k + 2L, floor(subsample_fraction * n_sub)),
          replace = FALSE
        ))
        fit <- safe_kmeans(
          Xs[idx, , drop = FALSE],
          k = k,
          seed = seed + 20000L * k + b,
          nstart = 20L
        )
        if (!isTRUE(fit$ok)) next

        included <- logical(n_sub)
        included[idx] <- TRUE
        lab <- rep(NA_integer_, n_sub)
        lab[idx] <- fit$labels

        observed_pair <- included[pairs$i] & included[pairs$j]
        obs_idx <- which(observed_pair)
        if (!length(obs_idx)) next

        co_observed[obs_idx] <- co_observed[obs_idx] + 1L
        same_cluster <- lab[pairs$i[obs_idx]] == lab[pairs$j[obs_idx]]
        co_clustered[obs_idx] <- co_clustered[obs_idx] +
          as.integer(same_cluster)
      }

      usable <- co_observed > 0L
      consensus <- co_clustered[usable] / co_observed[usable]
      if (!length(consensus)) {
        return(data.table(
          k = k,
          pac = NA_real_,
          n_pairs_observed = 0L,
          pair_coverage = 0,
          mean_consensus = NA_real_,
          cdf_area = NA_real_
        ))
      }

      data.table(
        k = k,
        pac = mean(consensus > lower & consensus < upper),
        n_pairs_observed = length(consensus),
        pair_coverage = mean(usable),
        mean_consensus = mean(consensus),
        cdf_area = 1 - mean(consensus)
      )
    })

    tab <- rbindlist(rows)
    setorder(tab, k)
    tab[, delta_area := c(NA_real_, diff(cdf_area))]

    if (!any(is.finite(tab$pac))) {
      stop("PAC was unavailable for every candidate K")
    }

    finite_pac <- which(is.finite(tab$pac))
    selected_min_pac <- if (length(finite_pac)) {
      tab$k[finite_pac[which.min(tab$pac[finite_pac])]]
    } else {
      NA_integer_
    }

    delta_candidates <- which(is.finite(tab$delta_area))
    selected_max_delta <- if (length(delta_candidates)) {
      tab$k[delta_candidates[which.max(tab$delta_area[delta_candidates])]]
    } else {
      NA_integer_
    }

    new_diag(
      ok = TRUE,
      stage = "consensus_pac",
      table = tab,
      summary = list(
        selected_k_min_pac = as.integer(selected_min_pac),
        selected_k_max_delta_area = as.integer(selected_max_delta)
      )
    )
  }, error = function(e) {
    new_diag(
      ok = FALSE,
      stage = "consensus_pac",
      error = conditionMessage(e)
    )
  })
}

consensus_kmeans_fit <- function(
  X,
  k,
  reps = 25L,
  subsample_fraction = 0.80,
  consensus_n = 600L,
  seed = 2001L
) {
  tryCatch({
    X <- as.matrix(X)
    n_sub <- min(nrow(X), as.integer(consensus_n))
    if (n_sub <= k + 5L) stop("Insufficient rows for consensus fit")

    set.seed(seed)
    idx_sub <- sort(sample.int(nrow(X), n_sub))
    Xs <- X[idx_sub, , drop = FALSE]

    co <- matrix(0, nrow = n_sub, ncol = n_sub)
    obs <- matrix(0, nrow = n_sub, ncol = n_sub)
    n_success <- 0L

    for (b in seq_len(reps)) {
      set.seed(seed + b)
      idx <- sort(sample.int(
        n_sub,
        max(k + 2L, floor(subsample_fraction * n_sub)),
        replace = FALSE
      ))
      fit <- safe_kmeans(
        Xs[idx, , drop = FALSE],
        k = k,
        seed = seed + 10000L + b,
        nstart = 20L
      )
      if (!isTRUE(fit$ok)) next

      same <- outer(fit$labels, fit$labels, FUN = "==")
      co[idx, idx] <- co[idx, idx] + same
      obs[idx, idx] <- obs[idx, idx] + 1
      n_success <- n_success + 1L
    }

    if (n_success < max(5L, floor(0.50 * reps))) {
      stop("Too few successful consensus resamples")
    }

    consensus <- co / pmax(obs, 1)
    consensus[obs == 0] <- 0.5
    diag(consensus) <- 1

    hc <- hclust(as.dist(1 - consensus), method = "average")
    sub_labels <- as.integer(cutree(hc, k = k))
    centers <- centroids_from_labels(Xs, sub_labels, k)

    if (any(!is.finite(centers))) {
      stop("Consensus-derived centroids contain non-finite values")
    }

    labels_all <- assign_centroids(X, centers)
    out <- list(
      ok = TRUE,
      labels = as.integer(labels_all),
      centers = centers,
      consensus_subsample_index = idx_sub,
      consensus_subsample_labels = sub_labels,
      n_success = n_success,
      fallback = FALSE,
      error = NA_character_
    )

    if (k == 2L) out <- relabel_k2_by_burden(out)
    out
  }, error = function(e) {
    fallback <- safe_kmeans(
      X,
      k = k,
      seed = seed + 999999L,
      nstart = 100L
    )
    if (k == 2L) fallback <- relabel_k2_by_burden(fallback)
    list(
      ok = isTRUE(fallback$ok),
      labels = fallback$labels,
      centers = fallback$centers,
      consensus_subsample_index = integer(),
      consensus_subsample_labels = integer(),
      n_success = 0L,
      fallback = TRUE,
      error = paste0(
        "Consensus fit failed; plain k-means fallback used. ",
        conditionMessage(e)
      )
    )
  })
}

## =============================================================================
## Additional Domain 1 diagnostics
## =============================================================================

dip_diagnostic <- function(X, labels_k2, seed, subsample_n = 1200L) {
  tryCatch({
    X <- as.matrix(X)
    n_sub <- min(nrow(X), as.integer(subsample_n))
    set.seed(seed)
    idx <- sort(sample.int(nrow(X), n_sub))
    Xs <- X[idx, , drop = FALSE]
    labs <- labels_k2[idx]

    centers <- centroids_from_labels(Xs, labs, 2L)
    direction <- centers[1L, ] - centers[2L, ]
    if (sum(direction^2) < 1e-12) direction[1L] <- 1
    discriminant_axis <- as.numeric(Xs %*% direction)

    pc <- prcomp(Xs, center = FALSE, scale. = FALSE)
    max_pc <- min(5L, ncol(pc$x))

    p_values <- c(
      discriminant_axis = diptest::dip.test(discriminant_axis)$p.value,
      vapply(seq_len(max_pc), function(j) {
        diptest::dip.test(pc$x[, j])$p.value
      }, numeric(1))
    )

    names(p_values)[-1L] <- paste0("PC", seq_len(max_pc))
    q_values <- p.adjust(p_values, method = "BH")

    tab <- data.table(
      axis = names(p_values),
      p_value = as.numeric(p_values),
      p_fdr = as.numeric(q_values)
    )

    new_diag(
      ok = TRUE,
      stage = "dip",
      table = tab,
      summary = list(
        discriminant_p = safe_num(
          tab[axis == "discriminant_axis", p_value]
        ),
        discriminant_p_fdr = safe_num(
          tab[axis == "discriminant_axis", p_fdr]
        ),
        any_axis_fdr_lt_0_05 = any(tab$p_fdr < 0.05, na.rm = TRUE)
      )
    )
  }, error = function(e) {
    new_diag(
      ok = FALSE,
      stage = "dip",
      error = conditionMessage(e)
    )
  })
}

cluster_index <- function(X, labels) {
  overall <- colMeans(X)
  total_ss <- sum((X - matrix(
    overall,
    nrow = nrow(X),
    ncol = ncol(X),
    byrow = TRUE
  ))^2)
  within_ss <- 0
  for (g in sort(unique(labels))) {
    idx <- which(labels == g)
    center <- colMeans(X[idx, , drop = FALSE])
    within_ss <- within_ss + sum((X[idx, , drop = FALSE] - matrix(
      center,
      nrow = length(idx),
      ncol = ncol(X),
      byrow = TRUE
    ))^2)
  }
  within_ss / total_ss
}

sigclust_like <- function(
  X,
  seed,
  B = 199L,
  subsample_n = 1000L
) {
  tryCatch({
    X <- as.matrix(X)
    n_sub <- min(nrow(X), as.integer(subsample_n))
    set.seed(seed)
    idx <- sort(sample.int(nrow(X), n_sub))
    Xs <- X[idx, , drop = FALSE]

    obs_fit <- safe_kmeans(Xs, 2L, seed = seed + 1L, nstart = 50L)
    if (!isTRUE(obs_fit$ok)) stop(obs_fit$error)
    obs_ci <- cluster_index(Xs, obs_fit$labels)

    sv <- svd(Xs, nu = 0L, nv = min(ncol(Xs), nrow(Xs) - 1L))
    eig <- (sv$d^2) / max(1, nrow(Xs) - 1L)
    eig <- pmax(eig, 1e-8)
    q <- length(eig)

    null_ci <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      set.seed(seed + 1000L + b)
      Z <- matrix(rnorm(n_sub * q), nrow = n_sub)
      Z <- sweep(Z, 2L, sqrt(eig), FUN = "*")
      fit_b <- safe_kmeans(
        Z,
        2L,
        seed = seed + 50000L + b,
        nstart = 20L
      )
      if (isTRUE(fit_b$ok)) {
        null_ci[b] <- cluster_index(Z, fit_b$labels)
      }
    }

    null_ci <- null_ci[is.finite(null_ci)]
    if (length(null_ci) < max(5L, floor(0.60 * B))) {
      stop("Too few successful SigClust-like null simulations")
    }

    p_value <- (1 + sum(null_ci <= obs_ci)) / (length(null_ci) + 1)
    null_mean <- mean(null_ci)
    null_sd <- sd(null_ci)
    effect_sd <- if (is.finite(null_sd) && null_sd > 0) {
      (null_mean - obs_ci) / null_sd
    } else {
      NA_real_
    }

    tab <- data.table(
      observed_cluster_index = obs_ci,
      null_mean_cluster_index = null_mean,
      null_sd_cluster_index = null_sd,
      p_value = p_value,
      effect_sd = effect_sd,
      relative_effect_pct = 100 * (null_mean - obs_ci) / null_mean,
      B_requested = B,
      B_successful = length(null_ci),
      mcse_p = sqrt(p_value * (1 - p_value) / (length(null_ci) + 1))
    )

    new_diag(
      ok = TRUE,
      stage = "sigclust_like",
      table = tab,
      summary = as.list(tab[1L])
    )
  }, error = function(e) {
    new_diag(
      ok = FALSE,
      stage = "sigclust_like",
      error = conditionMessage(e)
    )
  })
}

silhouette_diagnostic <- function(
  X,
  k_values = 2:4,
  seed,
  subsample_n = 1200L
) {
  tryCatch({
    X <- as.matrix(X)
    n_sub <- min(nrow(X), as.integer(subsample_n))
    set.seed(seed)
    idx <- sort(sample.int(nrow(X), n_sub))
    Xs <- X[idx, , drop = FALSE]
    d <- dist(Xs)

    tab <- rbindlist(lapply(k_values, function(k) {
      fit <- safe_kmeans(
        Xs,
        k = k,
        seed = seed + 1000L * k,
        nstart = 50L
      )
      sil <- if (isTRUE(fit$ok)) {
        safe_mean(cluster::silhouette(fit$labels, d)[, "sil_width"])
      } else {
        NA_real_
      }
      data.table(k = k, mean_silhouette = sil)
    }))

    finite_silhouette <- which(is.finite(tab$mean_silhouette))
    if (!length(finite_silhouette)) {
      stop("Silhouette was unavailable for every candidate K")
    }
    selected <- if (length(finite_silhouette)) {
      tab$k[
        finite_silhouette[
          which.max(tab$mean_silhouette[finite_silhouette])
        ]
      ]
    } else {
      NA_integer_
    }

    new_diag(
      ok = TRUE,
      stage = "silhouette",
      table = tab,
      summary = list(selected_k_max_silhouette = as.integer(selected))
    )
  }, error = function(e) {
    new_diag(
      ok = FALSE,
      stage = "silhouette",
      error = conditionMessage(e)
    )
  })
}

## =============================================================================
## Prediction and Rubin-style pooling
## =============================================================================

safe_auc <- function(y, p) {
  if (length(unique(y[is.finite(y)])) < 2L ||
      any(!is.finite(p))) return(NA_real_)
  as.numeric(
    suppressMessages(
      pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))
    )
  )
}

delong_delta_variance <- function(y, p0, p1) {
  tryCatch({
    r0 <- pROC::roc(y, p0, quiet = TRUE, direction = "<")
    r1 <- pROC::roc(y, p1, quiet = TRUE, direction = "<")
    v <- as.numeric(
      pROC::var(r1, method = "delong") +
        pROC::var(r0, method = "delong") -
        2 * pROC::cov(r1, r0, method = "delong")
    )
    max(v, 0)
  }, error = function(e) NA_real_)
}

fit_glm_predictions <- function(train_df, eval_df, outcome, predictors) {
  form <- reformulate(predictors, response = outcome)
  fit <- suppressWarnings(glm(
    form,
    data = train_df,
    family = binomial()
  ))
  as.numeric(predict(fit, newdata = eval_df, type = "response"))
}

fit_en_predictions <- function(
  x_train,
  y_train,
  x_eval,
  seed,
  alpha = 0.5,
  nfolds = 5L
) {
  set.seed(seed)
  folds <- min(nfolds, max(3L, floor(length(y_train) / 40L)))
  fit <- glmnet::cv.glmnet(
    x = as.matrix(x_train),
    y = y_train,
    family = "binomial",
    alpha = alpha,
    nfolds = folds,
    type.measure = "deviance",
    standardize = TRUE
  )
  as.numeric(predict(
    fit,
    newx = as.matrix(x_eval),
    s = "lambda.min",
    type = "response"
  ))
}

prediction_metrics_one_outcome <- function(
  train_dt,
  test_dt,
  external_dt,
  X_train,
  X_test,
  X_external,
  label_train,
  label_test,
  label_external,
  outcome,
  seed,
  algorithm = "kmeans"
) {
  tryCatch({
    tr <- as.data.frame(train_dt[, c(baseline_vars, outcome), with = FALSE])
    te <- as.data.frame(test_dt[, c(baseline_vars, outcome), with = FALSE])
    ex <- as.data.frame(external_dt[, c(baseline_vars, outcome), with = FALSE])

    tr$K2_label <- factor(label_train)
    te$K2_label <- factor(label_test, levels = levels(tr$K2_label))
    ex$K2_label <- factor(label_external, levels = levels(tr$K2_label))

    base_predictors <- baseline_vars
    base_k2_predictors <- c(baseline_vars, "K2_label")

    p_base_test <- fit_glm_predictions(
      tr, te, outcome, base_predictors
    )
    p_base_k2_test <- fit_glm_predictions(
      tr, te, outcome, base_k2_predictors
    )
    p_base_ext <- fit_glm_predictions(
      tr, ex, outcome, base_predictors
    )
    p_base_k2_ext <- fit_glm_predictions(
      tr, ex, outcome, base_k2_predictors
    )

    base_mm_tr <- model.matrix(
      ~ age + sex + sofa + aki_stage - 1,
      data = tr
    )
    base_mm_te <- model.matrix(
      ~ age + sex + sofa + aki_stage - 1,
      data = te
    )
    base_mm_ex <- model.matrix(
      ~ age + sex + sofa + aki_stage - 1,
      data = ex
    )

    raw_tr <- cbind(base_mm_tr, X_train)
    raw_te <- cbind(base_mm_te, X_test)
    raw_ex <- cbind(base_mm_ex, X_external)

    raw_k2_tr <- cbind(raw_tr, K2_label = as.numeric(label_train == 1L))
    raw_k2_te <- cbind(raw_te, K2_label = as.numeric(label_test == 1L))
    raw_k2_ex <- cbind(raw_ex, K2_label = as.numeric(label_external == 1L))

    y_train <- as.numeric(train_dt[[outcome]])
    y_test <- as.numeric(test_dt[[outcome]])
    y_ext <- as.numeric(external_dt[[outcome]])

    p_raw_test <- fit_en_predictions(
      raw_tr, y_train, raw_te, seed = seed + 10L
    )
    p_raw_k2_test <- fit_en_predictions(
      raw_k2_tr, y_train, raw_k2_te, seed = seed + 20L
    )
    p_raw_ext <- fit_en_predictions(
      raw_tr, y_train, raw_ex, seed = seed + 10L
    )
    p_raw_k2_ext <- fit_en_predictions(
      raw_k2_tr, y_train, raw_k2_ex, seed = seed + 20L
    )

    out <- rbindlist(list(
      data.table(
        validation = "internal_test",
        contrast = "Base+K2_minus_Base",
        auc_reference = safe_auc(y_test, p_base_test),
        auc_augmented = safe_auc(y_test, p_base_k2_test),
        delta_auc = safe_auc(y_test, p_base_k2_test) -
          safe_auc(y_test, p_base_test),
        within_imputation_variance = delong_delta_variance(
          y_test,
          p_base_test,
          p_base_k2_test
        )
      ),
      data.table(
        validation = "internal_test",
        contrast = "Raw-EN+K2_minus_Raw-EN",
        auc_reference = safe_auc(y_test, p_raw_test),
        auc_augmented = safe_auc(y_test, p_raw_k2_test),
        delta_auc = safe_auc(y_test, p_raw_k2_test) -
          safe_auc(y_test, p_raw_test),
        within_imputation_variance = delong_delta_variance(
          y_test,
          p_raw_test,
          p_raw_k2_test
        )
      ),
      data.table(
        validation = "external",
        contrast = "Base+K2_minus_Base",
        auc_reference = safe_auc(y_ext, p_base_ext),
        auc_augmented = safe_auc(y_ext, p_base_k2_ext),
        delta_auc = safe_auc(y_ext, p_base_k2_ext) -
          safe_auc(y_ext, p_base_ext),
        within_imputation_variance = delong_delta_variance(
          y_ext,
          p_base_ext,
          p_base_k2_ext
        )
      ),
      data.table(
        validation = "external",
        contrast = "Raw-EN+K2_minus_Raw-EN",
        auc_reference = safe_auc(y_ext, p_raw_ext),
        auc_augmented = safe_auc(y_ext, p_raw_k2_ext),
        delta_auc = safe_auc(y_ext, p_raw_k2_ext) -
          safe_auc(y_ext, p_raw_ext),
        within_imputation_variance = delong_delta_variance(
          y_ext,
          p_raw_ext,
          p_raw_k2_ext
        )
      )
    ))

    out[, outcome := outcome]
    out[, algorithm := algorithm]
    out
  }, error = function(e) {
    data.table(
      validation = c("internal_test", "internal_test", "external", "external"),
      contrast = rep(
        c("Base+K2_minus_Base", "Raw-EN+K2_minus_Raw-EN"),
        2L
      ),
      auc_reference = NA_real_,
      auc_augmented = NA_real_,
      delta_auc = NA_real_,
      within_imputation_variance = NA_real_,
      outcome = outcome,
      algorithm = algorithm,
      prediction_error = conditionMessage(e)
    )
  })
}

rubin_pool_delta_auc <- function(dt) {
  if (!nrow(dt)) return(data.table())
  dt[, {
    q <- as.numeric(delta_auc)
    u <- as.numeric(within_imputation_variance)
    keep <- is.finite(q)
    q <- q[keep]
    u <- u[keep]
    m <- length(q)

    if (!m) {
      list(
        m = 0L,
        pooled_delta_auc = NA_real_,
        total_variance = NA_real_,
        standard_error = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        ci_excludes_zero = NA
      )
    } else {
      qbar <- mean(q)
      ubar <- if (any(is.finite(u))) mean(u[is.finite(u)]) else 0
      b <- if (m > 1L) var(q) else 0
      total <- ubar + (1 + 1 / m) * b
      se <- sqrt(max(total, 0))
      list(
        m = m,
        pooled_delta_auc = qbar,
        total_variance = total,
        standard_error = se,
        ci_low = qbar - 1.96 * se,
        ci_high = qbar + 1.96 * se,
        ci_excludes_zero = (qbar - 1.96 * se) > 0
      )
    }
  }, by = c(
    "task_id", "repeat_id", "validation",
    "contrast", "outcome", "algorithm"
  )]
}

## =============================================================================
## Plotting helpers
## =============================================================================

save_plot_both <- function(p, filename_stem, width = 8, height = 5) {
  ggplot2::ggsave(
    filename = file.path(FIG_DIR, paste0(filename_stem, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
  ggplot2::ggsave(
    filename = file.path(FIG_DIR, paste0(filename_stem, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

theme_publication <- function() {
  theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title.position = "plot"
    )
}


################################################################################
## 03b_subtype_separation_pac_assignment_reproducibility_ARCHIVE.R
##
## Sensitivity analysis 03b:
##   subtype-separation dose-response under a true K = 3 data-generating model.
##
## Primary questions:
##   1. At what separation strength does the true K = 3 structure become
##      recoverable?
##   2. Can stable/prognostic K = 2 partitions appear when the latent classes
##      are not recoverable?
##   3. How do PAC, dip, SigClust-like effect, K2/K3 ARI, MICE stability,
##      and selected-K behavior change across separation levels?
##   4. Does frozen-centroid assignment reproducibility improve with
##      separation, and is K = 3 more reproducible than the forced K = 2
##      compression? NOTE: no external drift is applied here. The external
##      cohort is an independent draw from the SAME data-generating mechanism,
##      so this quantifies same-distribution assignment reproducibility, not
##      generalizability under drift (the latter is addressed in SA-06 / S3).
##   5. Does K2 add information beyond Raw-EN as separation changes?
##
## PAC is intentionally evaluated here rather than in the main simulation.
################################################################################

## =============================================================================
## User settings
## =============================================================================

PROJECT_DIR <- Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = "")
ANALYSIS_ID <- "sensitivity_03b_subtype_separation_pac_assignment_reproducibility_FINAL"
SIM_ROOT <- file.path(PROJECT_DIR, ANALYSIS_ID)

USE_TEST_MODE <- FALSE
N_REP_FORMAL <- 100L          # aligned with the pre-registered plan (B = 100 primary)
N_REP_TEST <- 1L
N_ACTIVE_REP <- if (USE_TEST_MODE) N_REP_TEST else N_REP_FORMAL

DELTA_GRID <- c(0.00, 0.75, 1.50, 2.25, 3.20, 4.00)

N_DERIVATION <- if (USE_TEST_MODE) 600L else 5000L
N_EXTERNAL <- if (USE_TEST_MODE) 400L else 3000L
TRAIN_PROP <- 0.70

MICE_M <- if (USE_TEST_MODE) 2L else 5L
MICE_MAXIT <- if (USE_TEST_MODE) 2L else 5L

PAC_K_VALUES <- 2:4
PAC_REPS <- if (USE_TEST_MODE) 5L else 50L
PAC_SUBSAMPLE_N <- if (USE_TEST_MODE) 250L else 800L
PAC_PAIR_N <- if (USE_TEST_MODE) 3000L else 50000L
PAC_SUBSAMPLE_FRACTION <- 0.80
PAC_LOWER <- 0.10
PAC_UPPER <- 0.90

SIGCLUST_B <- if (USE_TEST_MODE) 19L else 199L

## Parallelism tuned for a 16-core / 32 GB workstation.
## These tasks are CPU-bound and memory-light (each worker holds only small
## cohorts plus transient consensus/sigclust objects, well under ~1 GB), so the
## binding constraint is cores, not RAM. We leave 2 logical cores for the master
## R session + OS/RStudio and cap at 12 to avoid oversubscription. On Windows,
## future::multisession is the correct plan (multicore/fork is unavailable).
N_WORKERS <- if (USE_TEST_MODE) {
  2L
} else {
  max(1L, min(12L, future::availableCores() - 2L))
}
GLOBAL_SEED <- 2026071203L
RESUME_COMPLETED_TASKS <- TRUE
OVERWRITE_OUTPUT_TABLES <- TRUE

RUN_DOMAIN2 <- TRUE

RUN_MODE <- if (USE_TEST_MODE) "test" else "formal"

CONFIG_SIGNATURE <- paste(
  ANALYSIS_ID,
  RUN_MODE,
  paste(DELTA_GRID, collapse = ","),
  N_ACTIVE_REP,
  N_DERIVATION,
  N_EXTERNAL,
  TRAIN_PROP,
  MICE_M,
  MICE_MAXIT,
  PAC_REPS,
  PAC_SUBSAMPLE_N,
  PAC_PAIR_N,
  PAC_SUBSAMPLE_FRACTION,
  PAC_LOWER,
  PAC_UPPER,
  SIGCLUST_B,
  sep = "|"
)

OUT_DIR <- file.path(SIM_ROOT, paste0("output_", RUN_MODE))
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
RAW_DIR <- file.path(OUT_DIR, "raw_task_rds")

for (d in c(SIM_ROOT, OUT_DIR, TAB_DIR, FIG_DIR, LOG_DIR, RAW_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(
  LOG_DIR,
  paste0("simulation03b_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

set.seed(GLOBAL_SEED)

log_msg("Sensitivity 03b started.")
log_msg("Run mode: ", RUN_MODE)
log_msg("N_ACTIVE_REP: ", N_ACTIVE_REP)
log_msg("Delta settings: ", paste(DELTA_GRID, collapse = ", "))
log_msg("MICE m/maxit: ", MICE_M, "/", MICE_MAXIT)
log_msg("PAC reps/subsample/pairs: ", PAC_REPS, "/", PAC_SUBSAMPLE_N, "/", PAC_PAIR_N)
log_msg("SigClust-like B: ", SIGCLUST_B)
log_msg("Workers: ", N_WORKERS)

## =============================================================================
## Task-level simulation
## =============================================================================

run_03b_task <- function(task_row) {
  task_id <- as.character(task_row$task_id)
  repeat_id <- as.integer(task_row$repeat_id)
  delta <- as.numeric(task_row$delta)
  task_seed <- as.integer(task_row$task_seed)

  derivation_full <- generate_cohort(
    n = N_DERIVATION,
    delta = delta,
    subtype_prob = c(1/3, 1/3, 1/3),
    severity_mean = 0,
    external_drift = FALSE,
    seed = task_seed + 1L
  )

  external_full <- generate_cohort(
    n = N_EXTERNAL,
    delta = delta,
    subtype_prob = c(1/3, 1/3, 1/3),
    severity_mean = 0,
    external_drift = FALSE,
    seed = task_seed + 2L
  )

  split <- split_derivation(
    derivation_full,
    train_prop = TRAIN_PROP,
    seed = task_seed + 3L
  )

  train_missing <- inject_mar_missingness(
    split$analysis_train,
    seed = task_seed + 10L
  )
  test_missing <- inject_mar_missingness(
    split$internal_test,
    seed = task_seed + 11L
  )
  external_missing <- inject_mar_missingness(
    external_full,
    seed = task_seed + 12L
  )

  train_sets <- mice_complete_sets(
    train_missing,
    m = MICE_M,
    maxit = MICE_MAXIT,
    seed = task_seed + 100L
  )
  test_sets <- mice_complete_sets(
    test_missing,
    m = MICE_M,
    maxit = MICE_MAXIT,
    seed = task_seed + 200L
  )
  external_sets <- mice_complete_sets(
    external_missing,
    m = MICE_M,
    maxit = MICE_MAXIT,
    seed = task_seed + 300L
  )

  d1_list <- list()
  d2_list <- list()
  d3_list <- list()
  truth_composition_list <- list()
  k2_labels_by_imp <- list()
  k3_labels_by_imp <- list()

  pac_diag <- new_diag(stage = "consensus_pac")
  dip_diag <- new_diag(stage = "dip")
  sig_diag <- new_diag(stage = "sigclust_like")
  silhouette_diag <- new_diag(stage = "silhouette")

  for (imp in seq_len(MICE_M)) {
    train_dt <- train_sets[[imp]]
    test_dt <- test_sets[[imp]]
    ext_dt <- external_sets[[imp]]

    prep <- fit_preprocess(train_dt, feature_names)
    X_train <- apply_preprocess(train_dt, prep)
    X_test <- apply_preprocess(test_dt, prep)
    X_ext <- apply_preprocess(ext_dt, prep)

    ## Independent external reclustering uses external-local preprocessing.
    ## Frozen assignment and external prediction retain derivation preprocessing.
    prep_ext_local <- fit_preprocess(ext_dt, feature_names)
    X_ext_local <- apply_preprocess(ext_dt, prep_ext_local)

    fit_k2 <- safe_kmeans(
      X_train,
      k = 2L,
      seed = task_seed + 1000L + imp,
      nstart = 100L
    )
    fit_k2 <- relabel_k2_by_burden(fit_k2)

    fit_k3 <- safe_kmeans(
      X_train,
      k = 3L,
      seed = task_seed + 2000L + imp,
      nstart = 100L
    )

    if (!isTRUE(fit_k2$ok) || !isTRUE(fit_k3$ok)) {
      stop(
        "Primary k-means failed in task ", task_id,
        ", imputation ", imp,
        call. = FALSE
      )
    }

    label_k2_test <- assign_centroids(X_test, fit_k2$centers)
    label_k2_ext_frozen <- assign_centroids(X_ext, fit_k2$centers)
    label_k3_test <- assign_centroids(X_test, fit_k3$centers)
    label_k3_ext_frozen <- assign_centroids(X_ext, fit_k3$centers)

    ext_fit_k2 <- safe_kmeans(
      X_ext_local,
      k = 2L,
      seed = task_seed + 3000L + imp,
      nstart = 100L
    )
    ext_fit_k2 <- relabel_k2_by_burden(ext_fit_k2)

    ext_fit_k3 <- safe_kmeans(
      X_ext_local,
      k = 3L,
      seed = task_seed + 4000L + imp,
      nstart = 100L
    )

    if (!isTRUE(ext_fit_k2$ok) || !isTRUE(ext_fit_k3$ok)) {
      stop(
        "External reclustering failed in task ", task_id,
        ", imputation ", imp,
        call. = FALSE
      )
    }

    k2_labels_by_imp[[imp]] <- fit_k2$labels
    k3_labels_by_imp[[imp]] <- fit_k3$labels

    truth_k2 <- truth_structure_metrics(
      fit_k2$labels,
      train_dt$latent_subtype
    )
    truth_k3 <- truth_structure_metrics(
      fit_k3$labels,
      train_dt$latent_subtype
    )

    if (imp == 1L) {
      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = fit_k2$labels,
          truth = train_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 2L
        )
      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = fit_k3$labels,
          truth = train_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 3L,
          sample_role = "derivation"
        )

      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = label_k2_ext_frozen,
          truth = ext_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 2L,
          sample_role = "external_frozen"
        )
      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = ext_fit_k2$labels,
          truth = ext_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 2L,
          sample_role = "external_reclustered"
        )
      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = label_k3_ext_frozen,
          truth = ext_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 3L,
          sample_role = "external_frozen"
        )
      truth_composition_list[[length(truth_composition_list) + 1L]] <-
        truth_composition_table(
          labels = ext_fit_k3$labels,
          truth = ext_dt$latent_subtype,
          task_id = task_id,
          repeat_id = repeat_id,
          setting_name = as.character(delta),
          imputation_id = imp,
          algorithm = "plain_kmeans",
          k = 3L,
          sample_role = "external_reclustered"
        )
    }

    d1_list[[length(d1_list) + 1L]] <- rbindlist(list(
      data.table(
        task_id = task_id,
        repeat_id = repeat_id,
        delta = delta,
        imputation_id = imp,
        k = 2L,
        ari_vs_latent_subtype = adjusted_rand(
          fit_k2$labels,
          train_dt$latent_subtype
        ),
        truth_purity = truth_k2$purity,
        truth_weighted_cluster_entropy =
          truth_k2$weighted_cluster_entropy,
        truth_normalized_mutual_information =
          truth_k2$normalized_mutual_information,
        outcome_separation_C1_minus_C2 = outcome_separation(
          train_dt,
          fit_k2$labels,
          outcome = "mortality_30d"
        ),
        cluster_1_prevalence = mean(fit_k2$labels == 1L),
        cluster_2_prevalence = mean(fit_k2$labels == 2L),
        cluster_3_prevalence = NA_real_
      ),
      data.table(
        task_id = task_id,
        repeat_id = repeat_id,
        delta = delta,
        imputation_id = imp,
        k = 3L,
        ari_vs_latent_subtype = adjusted_rand(
          fit_k3$labels,
          train_dt$latent_subtype
        ),
        truth_purity = truth_k3$purity,
        truth_weighted_cluster_entropy =
          truth_k3$weighted_cluster_entropy,
        truth_normalized_mutual_information =
          truth_k3$normalized_mutual_information,
        outcome_separation_C1_minus_C2 = NA_real_,
        cluster_1_prevalence = mean(fit_k3$labels == 1L),
        cluster_2_prevalence = mean(fit_k3$labels == 2L),
        cluster_3_prevalence = mean(fit_k3$labels == 3L)
      )
    ))

    tm_k2 <- assignment_reproducibility_metrics(
      frozen_labels = label_k2_ext_frozen,
      recluster_labels = ext_fit_k2$labels,
      derivation_labels = fit_k2$labels,
      X_external = X_ext,
      derivation_centers = fit_k2$centers,
      k = 2L
    )
    tm_k2[, `:=`(
      task_id = task_id,
      repeat_id = repeat_id,
      delta = delta,
      imputation_id = imp,
      ari_frozen_vs_external_truth = adjusted_rand(
        label_k2_ext_frozen,
        ext_dt$latent_subtype
      ),
      ari_reclustered_vs_external_truth = adjusted_rand(
        ext_fit_k2$labels,
        ext_dt$latent_subtype
      )
    )]

    tm_k3 <- assignment_reproducibility_metrics(
      frozen_labels = label_k3_ext_frozen,
      recluster_labels = ext_fit_k3$labels,
      derivation_labels = fit_k3$labels,
      X_external = X_ext,
      derivation_centers = fit_k3$centers,
      k = 3L
    )
    tm_k3[, `:=`(
      task_id = task_id,
      repeat_id = repeat_id,
      delta = delta,
      imputation_id = imp,
      ari_frozen_vs_external_truth = adjusted_rand(
        label_k3_ext_frozen,
        ext_dt$latent_subtype
      ),
      ari_reclustered_vs_external_truth = adjusted_rand(
        ext_fit_k3$labels,
        ext_dt$latent_subtype
      )
    )]

    d3_list[[imp]] <- rbindlist(list(tm_k2, tm_k3), fill = TRUE)

    if (isTRUE(RUN_DOMAIN2)) {
      d2_imp <- rbindlist(lapply(
        outcome_vars,
        function(outcome) {
          prediction_metrics_one_outcome(
            train_dt = train_dt,
            test_dt = test_dt,
            external_dt = ext_dt,
            X_train = X_train,
            X_test = X_test,
            X_external = X_ext,
            label_train = fit_k2$labels,
            label_test = label_k2_test,
            label_external = label_k2_ext_frozen,
            outcome = outcome,
            seed = task_seed + 5000L + 100L * imp +
              match(outcome, outcome_vars),
            algorithm = "plain_kmeans"
          )
        }
      ))
      d2_imp[, `:=`(
        task_id = task_id,
        repeat_id = repeat_id,
        delta = delta,
        imputation_id = imp
      )]
      d2_list[[imp]] <- d2_imp
    }

    if (imp == 1L) {
      pac_diag <- consensus_pac_sampled_pairs(
        X_train,
        k_values = PAC_K_VALUES,
        reps = PAC_REPS,
        subsample_fraction = PAC_SUBSAMPLE_FRACTION,
        subsample_n = PAC_SUBSAMPLE_N,
        n_pairs = PAC_PAIR_N,
        lower = PAC_LOWER,
        upper = PAC_UPPER,
        seed = task_seed + 6000L
      )

      dip_diag <- dip_diagnostic(
        X_train,
        labels_k2 = fit_k2$labels,
        seed = task_seed + 7000L
      )

      sig_diag <- sigclust_like(
        X_train,
        seed = task_seed + 8000L,
        B = SIGCLUST_B
      )

      silhouette_diag <- silhouette_diagnostic(
        X_train,
        k_values = PAC_K_VALUES,
        seed = task_seed + 9000L
      )

    }
  }

  stability_rows <- rbindlist(lapply(c(2L, 3L), function(k) {
    labs <- if (k == 2L) k2_labels_by_imp else k3_labels_by_imp
    ## Exclude the trivial imputation-1-vs-itself comparison (ARI = 1), which
    ## otherwise inflates mean stability and the >= 0.90 proportion.
    compare_ids <- setdiff(seq_len(MICE_M), 1L)
    ari <- if (length(compare_ids)) {
      vapply(compare_ids, function(imp) {
        adjusted_rand(labs[[1L]], labs[[imp]])
      }, numeric(1))
    } else {
      NA_real_
    }
    finite_ari <- ari[is.finite(ari)]
    data.table(
      task_id = task_id,
      repeat_id = repeat_id,
      delta = delta,
      k = k,
      m = MICE_M,
      n_stability_comparisons = length(compare_ids),
      ari_mean_vs_primary = if (length(finite_ari)) mean(finite_ari) else NA_real_,
      ari_median_vs_primary = if (length(finite_ari)) median(finite_ari) else NA_real_,
      ari_min_vs_primary = if (length(finite_ari)) min(finite_ari) else NA_real_,
      ari_prop_ge_0_90 = if (length(finite_ari)) mean(finite_ari >= 0.90) else NA_real_
    )
  }))

  diagnostic_summary <- data.table(
    task_id = task_id,
    repeat_id = repeat_id,
    delta = delta,
    selected_k_pac = safe_int(
      pac_diag$summary$selected_k_min_pac
    ),
    selected_k_delta_area = safe_int(
      pac_diag$summary$selected_k_max_delta_area
    ),
    selected_k_silhouette = safe_int(
      silhouette_diag$summary$selected_k_max_silhouette
    ),
    dip_discriminant_p = safe_num(
      dip_diag$summary$discriminant_p
    ),
    dip_discriminant_p_fdr = safe_num(
      dip_diag$summary$discriminant_p_fdr
    ),
    dip_any_axis_fdr_lt_0_05 = isTRUE(
      dip_diag$summary$any_axis_fdr_lt_0_05
    ),
    sigclust_p = safe_num(sig_diag$summary$p_value),
    sigclust_effect_sd = safe_num(sig_diag$summary$effect_sd),
    sigclust_relative_effect_pct = safe_num(
      sig_diag$summary$relative_effect_pct
    ),
    pac_ok = isTRUE(pac_diag$ok),
    dip_ok = isTRUE(dip_diag$ok),
    sigclust_ok = isTRUE(sig_diag$ok),
    silhouette_ok = isTRUE(silhouette_diag$ok)
  )

  add_task_columns <- function(tab) {
    if (!nrow(tab)) return(tab)
    tab[, `:=`(
      task_id = task_id,
      repeat_id = repeat_id,
      delta = delta,
      imputation_id = 1L
    )]
    tab
  }

  list(
    ok = TRUE,
    config_signature = CONFIG_SIGNATURE,
    task_id = task_id,
    repeat_id = repeat_id,
    delta = delta,
    domain1_by_imputation = safe_rbindlist(d1_list),
    domain2_by_imputation = safe_rbindlist(d2_list),
    domain3_by_imputation = safe_rbindlist(d3_list),
    truth_composition_primary = safe_rbindlist(truth_composition_list),
    mice_stability = stability_rows,
    diagnostic_summary = diagnostic_summary,
    pac_by_k = add_task_columns(copy(pac_diag$table)),
    dip_by_axis = add_task_columns(copy(dip_diag$table)),
    sigclust = add_task_columns(copy(sig_diag$table)),
    silhouette_by_k = add_task_columns(copy(silhouette_diag$table)),
    diagnostic_errors = data.table(
      task_id = task_id,
      repeat_id = repeat_id,
      delta = delta,
      pac_error = pac_diag$error,
      dip_error = dip_diag$error,
      sigclust_error = sig_diag$error,
      silhouette_error = silhouette_diag$error
    )
  )
}

## =============================================================================
## Parallel execution with task-level checkpoints
## =============================================================================

tasks <- CJ(
  delta = DELTA_GRID,
  repeat_id = seq_len(N_ACTIVE_REP)
)
tasks[, delta_label := gsub("\\.", "p", sprintf("%.2f", delta))]
tasks[, task_id := sprintf("delta_%s_rep_%03d", delta_label, repeat_id)]
tasks[, task_seed := GLOBAL_SEED +
  as.integer(factor(delta, levels = DELTA_GRID)) * 100000L +
  repeat_id]

fwrite(
  tasks,
  file.path(TAB_DIR, "sensitivity03b_task_manifest.csv")
)

config <- data.table(
  parameter = c(
    "ANALYSIS_ID", "RUN_MODE", "N_REP", "N_DERIVATION", "N_EXTERNAL",
    "TRAIN_PROP", "MICE_M", "MICE_MAXIT", "PAC_REPS",
    "PAC_SUBSAMPLE_N", "PAC_PAIR_N", "PAC_SUBSAMPLE_FRACTION",
    "PAC_LOWER", "PAC_UPPER", "SIGCLUST_B",
    "N_WORKERS", "GLOBAL_SEED"
  ),
  value = c(
    ANALYSIS_ID, RUN_MODE, N_ACTIVE_REP, N_DERIVATION, N_EXTERNAL,
    TRAIN_PROP, MICE_M, MICE_MAXIT, PAC_REPS,
    PAC_SUBSAMPLE_N, PAC_PAIR_N, PAC_SUBSAMPLE_FRACTION,
    PAC_LOWER, PAC_UPPER, SIGCLUST_B,
    N_WORKERS, GLOBAL_SEED
  )
)
fwrite(config, file.path(TAB_DIR, "sensitivity03b_configuration.csv"))

if (N_WORKERS > 1L) {
  ## Guard against the default 500 MiB globals limit when exporting task
  ## metadata and helper closures to multisession workers.
  options(future.globals.maxSize = 2 * 1024^3)
  future::plan(future::multisession, workers = N_WORKERS)
} else {
  future::plan(future::sequential)
}

log_msg("Total tasks: ", nrow(tasks))
log_msg("Parallel workers: ", N_WORKERS)

worker_fun <- function(i) {
  data.table::setDTthreads(1L)
  task_row <- tasks[i]
  message(
    "[", timestamp(), "] Running ",
    task_row$task_id,
    " | delta=", task_row$delta,
    " | repeat ", task_row$repeat_id, "/", N_ACTIVE_REP
  )
  checkpoint <- file.path(
    RAW_DIR,
    paste0(task_row$task_id, ".rds")
  )

  if (isTRUE(RESUME_COMPLETED_TASKS) && file.exists(checkpoint)) {
    old <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
    if (
      is.list(old) &&
      isTRUE(old$ok) &&
      identical(old$config_signature, CONFIG_SIGNATURE)
    ) {
      return(old)
    }
  }

  result <- tryCatch(
    run_03b_task(task_row),
    error = function(e) {
      list(
        ok = FALSE,
        config_signature = CONFIG_SIGNATURE,
        task_id = as.character(task_row$task_id),
        repeat_id = as.integer(task_row$repeat_id),
        delta = as.numeric(task_row$delta),
        error = conditionMessage(e)
      )
    }
  )

  saveRDS(result, checkpoint)
  message(
    "[", timestamp(), "] Finished ",
    task_row$task_id,
    " | status=", if (isTRUE(result$ok)) "ok" else "error"
  )
  result
}

results <- future.apply::future_lapply(
  seq_len(nrow(tasks)),
  worker_fun,
  future.seed = TRUE,
  future.scheduling = 1
)

future::plan(future::sequential)
invisible(gc(full = TRUE))

error_log <- rbindlist(lapply(results, function(z) {
  if (is.list(z) && !isTRUE(z$ok)) {
    data.table(
      task_id = z$task_id,
      repeat_id = z$repeat_id,
      delta = z$delta,
      error = z$error
    )
  } else {
    NULL
  }
}), fill = TRUE)

write_csv_safe(
  error_log,
  file.path(LOG_DIR, "sensitivity03b_parallel_error_log.csv")
)

success <- Filter(function(z) is.list(z) && isTRUE(z$ok), results)

if (!length(success)) {
  stop(
    "All sensitivity 03b tasks failed. See ",
    file.path(LOG_DIR, "sensitivity03b_parallel_error_log.csv"),
    call. = FALSE
  )
}

log_msg(
  "Successful tasks: ", length(success),
  " / ", nrow(tasks),
  "; failed: ", nrow(error_log)
)

## =============================================================================
## Aggregation
## =============================================================================

d1_imp <- safe_rbindlist(lapply(success, `[[`, "domain1_by_imputation"))
d2_imp <- safe_rbindlist(lapply(success, `[[`, "domain2_by_imputation"))
d3_imp <- safe_rbindlist(lapply(success, `[[`, "domain3_by_imputation"))
truth_composition_primary <- safe_rbindlist(
  lapply(success, `[[`, "truth_composition_primary")
)
mice_stability <- safe_rbindlist(lapply(success, `[[`, "mice_stability"))
diag_summary <- safe_rbindlist(lapply(success, `[[`, "diagnostic_summary"))
pac_by_k <- safe_rbindlist(lapply(success, `[[`, "pac_by_k"))
dip_by_axis <- safe_rbindlist(lapply(success, `[[`, "dip_by_axis"))
sigclust_dt <- safe_rbindlist(lapply(success, `[[`, "sigclust"))
silhouette_by_k <- safe_rbindlist(lapply(success, `[[`, "silhouette_by_k"))
diagnostic_errors <- safe_rbindlist(
  lapply(success, `[[`, "diagnostic_errors")
)

d1_repeat <- d1_imp[, .(
  ari_vs_latent_subtype = safe_mean(ari_vs_latent_subtype),
  truth_purity = safe_mean(truth_purity),
  truth_weighted_cluster_entropy = safe_mean(
    truth_weighted_cluster_entropy
  ),
  truth_normalized_mutual_information = safe_mean(
    truth_normalized_mutual_information
  ),
  outcome_separation_C1_minus_C2 = safe_mean(
    outcome_separation_C1_minus_C2
  ),
  cluster_1_prevalence = safe_mean(cluster_1_prevalence)
), by = .(task_id, repeat_id, delta, k)]

d3_repeat <- d3_imp[, .(
  frozen_vs_reclustered_ari = safe_mean(
    frozen_vs_reclustered_ari
  ),
  matched_agreement = safe_mean(matched_agreement),
  mean_abs_prevalence_difference = safe_mean(
    mean_abs_prevalence_difference
  ),
  max_abs_prevalence_difference = safe_mean(
    max_abs_prevalence_difference
  ),
  centroid_profile_correlation = safe_mean(
    centroid_profile_correlation
  ),
  ari_frozen_vs_external_truth = safe_mean(
    ari_frozen_vs_external_truth
  ),
  ari_reclustered_vs_external_truth = safe_mean(
    ari_reclustered_vs_external_truth
  )
), by = .(task_id, repeat_id, delta, k)]


d1_k_wide <- dcast(
  d1_repeat,
  task_id + repeat_id + delta ~ k,
  value.var = "ari_vs_latent_subtype"
)
setnames(
  d1_k_wide,
  old = intersect(c("2", "3"), names(d1_k_wide)),
  new = c("ari_k2", "ari_k3")[seq_along(
    intersect(c("2", "3"), names(d1_k_wide))
  )]
)
d1_k_comparison <- if (all(c("ari_k2", "ari_k3") %in% names(d1_k_wide))) {
  d1_k_wide[, .(
    task_id,
    repeat_id,
    delta,
    ari_k2,
    ari_k3,
    ari_k3_minus_k2 = ari_k3 - ari_k2
  )]
} else {
  data.table()
}

d3_k_wide <- dcast(
  d3_repeat,
  task_id + repeat_id + delta ~ k,
  value.var = c(
    "frozen_vs_reclustered_ari",
    "matched_agreement",
    "mean_abs_prevalence_difference",
    "centroid_profile_correlation"
  )
)

d3_k_comparison <- if (all(c(
  "frozen_vs_reclustered_ari_2",
  "frozen_vs_reclustered_ari_3"
) %in% names(d3_k_wide))) {
  d3_k_wide[, .(
    task_id,
    repeat_id,
    delta,
    assignment_reproducibility_ari_k2 = frozen_vs_reclustered_ari_2,
    assignment_reproducibility_ari_k3 = frozen_vs_reclustered_ari_3,
    assignment_reproducibility_ari_k3_minus_k2 =
      frozen_vs_reclustered_ari_3 -
      frozen_vs_reclustered_ari_2,
    matched_agreement_k3_minus_k2 =
      matched_agreement_3 - matched_agreement_2,
    prevalence_error_k3_minus_k2 =
      mean_abs_prevalence_difference_3 -
      mean_abs_prevalence_difference_2,
    profile_correlation_k3_minus_k2 =
      centroid_profile_correlation_3 -
      centroid_profile_correlation_2
  )]
} else {
  data.table()
}

diagnostic_failure_summary <- diag_summary[, .(
  n_repeats = .N,
  pac_failure_rate = mean(!pac_ok),
  dip_failure_rate = mean(!dip_ok),
  sigclust_failure_rate = mean(!sigclust_ok),
  silhouette_failure_rate = mean(!silhouette_ok)
), by = delta]

d2_rubin <- if (nrow(d2_imp)) {
  rubin_pool_delta_auc(d2_imp)
} else {
  data.table()
}
if (nrow(d2_rubin)) {
  d2_rubin <- merge(
    d2_rubin,
    unique(d2_imp[, .(task_id, delta)]),
    by = "task_id",
    all.x = TRUE
  )
}

make_selected_k_distribution_03 <- function(values, method_name) {
  out <- data.table(
    delta = diag_summary$delta,
    selected_k = as.integer(values)
  )
  out <- out[is.finite(selected_k), .(n = .N), by = .(delta, selected_k)]
  out[, proportion := n / sum(n), by = delta]
  out[, method := method_name]
  setcolorder(out, c("delta", "method", "selected_k", "n", "proportion"))
  out
}

selected_k_distribution <- rbindlist(list(
  make_selected_k_distribution_03(
    diag_summary$selected_k_pac,
    "PAC_minimum"
  ),
  make_selected_k_distribution_03(
    diag_summary$selected_k_delta_area,
    "PAC_delta_area"
  ),
  make_selected_k_distribution_03(
    diag_summary$selected_k_silhouette,
    "Silhouette"
  )
), fill = TRUE)

## =============================================================================
## Operating characteristics and Monte Carlo error
## =============================================================================

op_list <- list()

add_op <- function(source, value, group_cols, metric) {
  tmp <- copy(source)
  expr <- parse(text = value)[[1L]]
  tmp[, .op_value := as.numeric(
    eval(expr, envir = as.list(tmp), enclos = parent.frame())
  )]
  summarise_binary(tmp, group_cols, ".op_value", metric)
}

op_list[[length(op_list) + 1L]] <- add_op(
  d1_repeat[k == 3L],
  "ari_vs_latent_subtype >= 0.90",
  "delta",
  "K3_true_label_recovery_ARI_ge_0.90"
)

op_list[[length(op_list) + 1L]] <- add_op(
  d1_repeat[k == 2L],
  "ari_vs_latent_subtype >= 0.40 & ari_vs_latent_subtype <= 0.70",
  "delta",
  "K2_compression_ARI_between_0.40_and_0.70"
)

op_list[[length(op_list) + 1L]] <- add_op(
  diag_summary,
  "selected_k_pac == 3L",
  "delta",
  "PAC_selects_K3"
)

op_list[[length(op_list) + 1L]] <- add_op(
  diag_summary,
  "dip_discriminant_p_fdr < 0.05",
  "delta",
  "Dip_discriminant_FDR_rejection"
)

op_list[[length(op_list) + 1L]] <- add_op(
  diag_summary,
  "sigclust_p < 0.05",
  "delta",
  "SigClust_like_rejection"
)

op_list[[length(op_list) + 1L]] <- add_op(
  d3_repeat[k == 3L],
  "frozen_vs_reclustered_ari >= 0.80",
  "delta",
  ## No external drift is applied in SA-03b; this is frozen-assignment
  ## reproducibility on an independent same-distribution replicate.
  "K3_frozen_assignment_reproducibility_ARI_ge_0.80"
)

op_list[[length(op_list) + 1L]] <- add_op(
  mice_stability[k == 3L],
  "ari_mean_vs_primary >= 0.90",
  "delta",
  "K3_MICE_stability_ARI_ge_0.90"
)

if (nrow(d2_rubin)) {
  ## Replaces the previous abs(dAUC) < 0.001 flag, which was too tight to fire
  ## given the Monte Carlo SD of pooled dAUC and so carried no information.
  ## We instead report (i) an equivalence band and (ii) a decisive-increment
  ## rate based on whether the Rubin 95% CI excludes zero.
  op_list[[length(op_list) + 1L]] <- add_op(
    d2_rubin[contrast == "Raw-EN+K2_minus_Raw-EN"],
    "abs(pooled_delta_auc) < 0.005",
    c("delta", "validation", "outcome"),
    "Raw_EN_increment_within_equivalence_band_abs_dAUC_lt_0.005"
  )
  op_list[[length(op_list) + 1L]] <- add_op(
    d2_rubin[contrast == "Raw-EN+K2_minus_Raw-EN"],
    "ci_excludes_zero",
    c("delta", "validation", "outcome"),
    "Raw_EN_increment_decisive_Rubin_CI_excludes_zero"
  )
}

operating_characteristics <- safe_rbindlist(op_list)

mcse_summary <- safe_rbindlist(list(
  summarise_numeric(
    d1_repeat,
    c("delta", "k"),
    "ari_vs_latent_subtype",
    "ARI_vs_latent_subtype"
  ),
  summarise_numeric(
    d1_repeat,
    c("delta", "k"),
    "truth_purity",
    "Truth_purity"
  ),
  summarise_numeric(
    d1_repeat,
    c("delta", "k"),
    "truth_normalized_mutual_information",
    "Truth_normalized_mutual_information"
  ),
  summarise_numeric(
    d1_repeat[k == 2L],
    "delta",
    "outcome_separation_C1_minus_C2",
    "K2_outcome_separation"
  ),
  summarise_numeric(
    d3_repeat,
    c("delta", "k"),
    "frozen_vs_reclustered_ari",
    "Frozen_vs_reclustered_ARI"
  ),
  summarise_numeric(
    mice_stability,
    c("delta", "k"),
    "ari_mean_vs_primary",
    "MICE_label_stability_ARI"
  ),
  if (nrow(d2_rubin)) {
    summarise_numeric(
      d2_rubin,
      c("delta", "validation", "outcome", "contrast"),
      "pooled_delta_auc",
      "Rubin_pooled_delta_AUC"
    )
  } else {
    data.table()
  }
))

## =============================================================================
## ADEMP table
## =============================================================================

ademp <- data.table(
  component = c(
    "Aim",
    "Data-generating mechanism",
    "Estimands",
    "Methods",
    "Performance measures",
    "Repetitions and Monte Carlo error",
    "Primary interpretation rule"
  ),
  specification = c(
    paste(
      "Evaluate dose-response of discreteness, K recovery, PAC,",
      "same-distribution frozen-assignment reproducibility, and predictive",
      "redundancy as true K=3 subtype separation increases."
    ),
    paste0(
      "Three latent classes; delta in {",
      paste(DELTA_GRID, collapse = ", "),
      "}; continuous severity remains outcome-generating; derivation/external",
      " cohorts share the same DGM."
    ),
    paste(
      "K2/K3 ARI to latent truth; PAC for K=2-4; selected K;",
      "dip and SigClust-like diagnostics; MICE ARI;",
      "frozen-vs-reclustered ARI; prevalence/profile correspondence;",
      "internal/external delta AUC."
    ),
    paste(
      "MICE m=5; 1%/99% winsorization; derivation standardization;",
      "plain k-means with frozen centroids; sampled-pair PAC;",
      "held-out prediction with Base, Base+K2, Raw-EN, Raw-EN+K2."
    ),
    paste(
      "Mean, SD, empirical 2.5%-97.5% interval, MCSE;",
      "selected-K proportion; recovery/assignment-reproducibility/stability proportions;",
      "task failure and diagnostic failure rates."
    ),
    paste0(
      N_ACTIVE_REP,
      " repeats per delta in this run; MCSE reported for all principal",
      " means and proportions."
    ),
    paste(
      "No single index proves discreteness. Evidence requires concordance",
      "among K3 recovery, PAC/selected K, multimodality diagnostics,",
      "MICE stability, and same-distribution assignment reproducibility."
    )
  )
)

## =============================================================================
## Save tables
## =============================================================================

tables_to_write <- list(
  sensitivity03b_domain1_by_imputation = d1_imp,
  sensitivity03b_domain1_by_repeat = d1_repeat,
  sensitivity03b_truth_composition_primary = truth_composition_primary,
  sensitivity03b_k2_k3_recovery_comparison = d1_k_comparison,
  sensitivity03b_pac_by_k = pac_by_k,
  sensitivity03b_dip_by_axis = dip_by_axis,
  sensitivity03b_sigclust_like = sigclust_dt,
  sensitivity03b_silhouette_by_k = silhouette_by_k,
  sensitivity03b_diagnostic_summary = diag_summary,
  sensitivity03b_selected_k_distribution = selected_k_distribution,
  sensitivity03b_mice_label_stability = mice_stability,
  sensitivity03b_domain2_by_imputation = d2_imp,
  sensitivity03b_domain2_rubin_by_repeat = d2_rubin,
  sensitivity03b_same_distribution_assignment_reproducibility_by_imputation = d3_imp,
  sensitivity03b_same_distribution_assignment_reproducibility_by_repeat = d3_repeat,
  sensitivity03b_k2_k3_assignment_reproducibility_comparison = d3_k_comparison,
  sensitivity03b_operating_characteristics = operating_characteristics,
  sensitivity03b_diagnostic_failure_summary = diagnostic_failure_summary,
  sensitivity03b_mcse_summary = mcse_summary,
  sensitivity03b_ADEMP_table = ademp,
  sensitivity03b_diagnostic_errors = diagnostic_errors,
  sensitivity03b_task_failure_log = error_log
)

for (nm in names(tables_to_write)) {
  write_csv_safe(
    tables_to_write[[nm]],
    file.path(TAB_DIR, paste0(nm, ".csv"))
  )
}

## =============================================================================
## Figures
## =============================================================================

if (nrow(pac_by_k)) {
  pac_plot_dt <- pac_by_k[, .(
    mean = safe_mean(pac),
    low = safe_quantile(pac, 0.025),
    high = safe_quantile(pac, 0.975)
  ), by = .(delta, k)]

  p_pac <- ggplot(
    pac_plot_dt,
    aes(
      x = delta,
      y = mean,
      group = factor(k),
      linetype = factor(k),
      shape = factor(k)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.08
    ) +
    labs(
      x = "Subtype separation (delta)",
      y = "PAC",
      linetype = "K",
      shape = "K",
      title = "Consensus PAC across subtype-separation levels"
    ) +
    theme_publication()

  save_plot_both(
    p_pac,
    "Figure_S03b_PAC_by_separation",
    width = 8,
    height = 5
  )
}

if (nrow(d1_repeat)) {
  ari_plot_dt <- d1_repeat[, .(
    mean = safe_mean(ari_vs_latent_subtype),
    low = safe_quantile(ari_vs_latent_subtype, 0.025),
    high = safe_quantile(ari_vs_latent_subtype, 0.975)
  ), by = .(delta, k)]

  p_ari <- ggplot(
    ari_plot_dt,
    aes(
      x = delta,
      y = mean,
      group = factor(k),
      linetype = factor(k),
      shape = factor(k)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.08
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Subtype separation (delta)",
      y = "Adjusted Rand index versus latent truth",
      linetype = "K",
      shape = "K",
      title = "Recovery of the latent three-class structure"
    ) +
    theme_publication()

  save_plot_both(
    p_ari,
    "Figure_S03b_true_label_recovery",
    width = 8,
    height = 5
  )
}

if (nrow(d3_repeat)) {
  tr_plot_dt <- d3_repeat[, .(
    mean = safe_mean(frozen_vs_reclustered_ari),
    low = safe_quantile(frozen_vs_reclustered_ari, 0.025),
    high = safe_quantile(frozen_vs_reclustered_ari, 0.975)
  ), by = .(delta, k)]

  p_tr <- ggplot(
    tr_plot_dt,
    aes(
      x = delta,
      y = mean,
      group = factor(k),
      linetype = factor(k),
      shape = factor(k)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.08
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Subtype separation (delta)",
      y = "Frozen-versus-reclustered ARI",
      linetype = "K",
      shape = "K",
      title = "Frozen-assignment reproducibility (independent same-distribution replicate, no drift)"
    ) +
    theme_publication()

  save_plot_both(
    p_tr,
    "Figure_S03b_same_distribution_assignment_reproducibility",
    width = 8,
    height = 5
  )
}

if (nrow(selected_k_distribution)) {
  p_k <- ggplot(
    selected_k_distribution[
      method %in% c("PAC_minimum", "Silhouette")
    ],
    aes(
      x = factor(selected_k),
      y = proportion,
      fill = factor(selected_k)
    )
  ) +
    geom_col() +
    facet_grid(method ~ delta, labeller = label_both) +
    labs(
      x = "Selected K",
      y = "Proportion of repeats",
      fill = "Selected K",
      title = "Selected-K distributions"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  save_plot_both(
    p_k,
    "Figure_S03b_selected_K_distribution",
    width = 11,
    height = 7
  )
}

if (nrow(d2_rubin)) {
  dauc_plot_dt <- d2_rubin[
    contrast == "Raw-EN+K2_minus_Raw-EN",
    .(
      mean = safe_mean(pooled_delta_auc),
      low = safe_quantile(pooled_delta_auc, 0.025),
      high = safe_quantile(pooled_delta_auc, 0.975)
    ),
    by = .(delta, validation, outcome)
  ]

  p_dauc <- ggplot(
    dauc_plot_dt,
    aes(
      x = delta,
      y = mean,
      shape = outcome,
      linetype = outcome
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.08
    ) +
    facet_wrap(~ validation, scales = "free_y") +
    labs(
      x = "Subtype separation (delta)",
      y = "Raw-EN+K2 minus Raw-EN delta AUC",
      shape = "Outcome",
      linetype = "Outcome",
      title = "Incremental information beyond raw physiology"
    ) +
    theme_publication()

  save_plot_both(
    p_dauc,
    "Figure_S03b_decisive_delta_AUC",
    width = 9,
    height = 5
  )
}

## =============================================================================
## Completion log
## =============================================================================

completion <- data.table(
  analysis_id = ANALYSIS_ID,
  run_mode = RUN_MODE,
  expected_tasks = nrow(tasks),
  successful_tasks = length(success),
  failed_tasks = nrow(error_log),
  completed_at = timestamp()
)

fwrite(
  completion,
  file.path(LOG_DIR, "sensitivity03b_run_completion.csv")
)
if (nrow(error_log) == 0L && length(success) == nrow(tasks)) {
  writeLines(
    paste0(
      "completed_at=", timestamp(), "\n",
      "run_mode=", RUN_MODE, "\n",
      "expected_tasks=", nrow(tasks), "\n",
      "successful_tasks=", length(success), "\n",
      "failed_tasks=0\n"
    ),
    file.path(LOG_DIR, "run_completed.ok")
  )
} else {
  writeLines(
    paste0(
      "partial_at=", timestamp(), "\n",
      "run_mode=", RUN_MODE, "\n",
      "expected_tasks=", nrow(tasks), "\n",
      "successful_tasks=", length(success), "\n",
      "failed_tasks=", nrow(error_log), "\n"
    ),
    file.path(LOG_DIR, "run_partial_with_failures.txt")
  )
}

write_session_info(
  file.path(LOG_DIR, "sessionInfo.txt")
)

log_msg("Sensitivity 03b completed.")
log_msg("Tables: ", TAB_DIR)
log_msg("Figures: ", FIG_DIR)
