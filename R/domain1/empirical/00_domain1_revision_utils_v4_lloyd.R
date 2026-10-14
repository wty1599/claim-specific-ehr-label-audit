## Shared utilities for the governed Domain 1 Lloyd-aligned package.
##
## This file is a new implementation snapshot. It does not modify the v3 source.
## The formal K2 contract is Lloyd, nstart=25, iter.max=100, with no
## clustering-algorithm fallback.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

assert_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop("Missing required input [", label, "]: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

first_existing_column <- function(dt, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(dt)]
  if (!length(hit)) {
    if (required) {
      stop("None of the required columns exist: ", paste(candidates, collapse = ", "), call. = FALSE)
    }
    return(NULL)
  }
  hit[1L]
}

wilson_interval <- function(x, n, conf.level = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0 || x < 0 || x > n) {
    return(c(low = NA_real_, high = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(low = max(0, centre - half), high = min(1, centre + half))
}

bh_discriminant_p <- function(dt, suffix = "") {
  stems <- c(
    "dip_discriminant_p", "dip_pc1_p", "dip_pc2_p", "dip_pc3_p",
    "dip_pc4_p", "dip_pc5_p"
  )
  cols <- paste0(stems, suffix)
  if (!all(cols %in% names(dt))) {
    return(rep(NA_real_, nrow(dt)))
  }
  vapply(seq_len(nrow(dt)), function(i) {
    p <- safe_numeric(unlist(dt[i, ..cols], use.names = FALSE))
    if (length(p) != 6L || any(!is.finite(p))) return(NA_real_)
    stats::p.adjust(p, method = "BH")[1L]
  }, numeric(1))
}

classify_domain1 <- function(dip_p_fdr, mahalanobis_delta,
                             mahalanobis_threshold,
                             dip_alpha = D1_DIP_FDR_ALPHA) {
  missing <- !is.finite(dip_p_fdr) | !is.finite(mahalanobis_delta) |
    !is.finite(mahalanobis_threshold)
  positive <- !missing & dip_p_fdr < dip_alpha &
    mahalanobis_delta >= mahalanobis_threshold
  negative <- !missing & dip_p_fdr >= dip_alpha &
    mahalanobis_delta < mahalanobis_threshold
  out <- rep("inconclusive", length(dip_p_fdr))
  out[negative] <- "no_discrete_evidence"
  out[positive] <- "discrete_evidence"
  out[missing] <- "failed"
  out
}

summarise_numeric_vector <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  data.table::data.table(
    n = n,
    mean = if (n) mean(x) else NA_real_,
    sd = if (n > 1L) stats::sd(x) else NA_real_,
    mcse = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_,
    median = if (n) stats::median(x) else NA_real_,
    p025 = if (n) as.numeric(stats::quantile(x, 0.025, names = FALSE)) else NA_real_,
    p975 = if (n) as.numeric(stats::quantile(x, 0.975, names = FALSE)) else NA_real_
  )
}

safe_kmeans <- function(X, centers = 2L, seed = 1L,
                        nstart = 25L, iter.max = 100L) {
  failure <- function(message, warnings = character()) {
    structure(
      list(
        error = message,
        algorithm_used = "Lloyd",
        fallback_used = FALSE,
        fallback_reason = NA_character_,
        warnings = warnings,
        nstart = 25L,
        iter.max = 100L,
        requested_nstart = nstart,
        requested_iter.max = iter.max,
        requested_centers = centers
      ),
      class = "d1_kmeans_error"
    )
  }

  if (!is.numeric(nstart) || length(nstart) != 1L || !is.finite(nstart) ||
      nstart != 25) {
    return(failure("Governed Lloyd contract requires nstart=25."))
  }
  if (!is.numeric(iter.max) || length(iter.max) != 1L || !is.finite(iter.max) ||
      iter.max != 100) {
    return(failure("Governed Lloyd contract requires iter.max=100."))
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    return(failure("seed must be one finite integer-like numeric value."))
  }
  if (!is.numeric(centers) || length(centers) != 1L || !is.finite(centers) ||
      as.integer(centers) != centers || centers != 2L) {
    return(failure("Governed Domain 1 contract requires centers=2."))
  }

  X <- tryCatch(as.matrix(X), error = function(e) NULL)
  if (is.null(X) || !is.numeric(X) || nrow(X) < as.integer(centers) ||
      ncol(X) < 1L || any(!is.finite(X))) {
    return(failure(
      "X must be a finite numeric matrix with at least centers rows and one column."
    ))
  }

  warning_text <- character()
  error_text <- NA_character_
  set.seed(as.integer(seed))
  fit <- tryCatch(
    withCallingHandlers(
      stats::kmeans(
        X,
        centers = as.integer(centers),
        nstart = 25L,
        iter.max = 100L,
        algorithm = "Lloyd"
      ),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    return(failure(paste0("Lloyd k-means failed: ", error_text), warning_text))
  }
  if (length(warning_text)) {
    return(failure(
      paste0(
        "Lloyd k-means emitted a warning and is invalid under the governed contract: ",
        paste(unique(warning_text), collapse = " | ")
      ),
      warning_text
    ))
  }
  if (!is.null(fit$ifault) && length(fit$ifault) == 1L &&
      is.finite(fit$ifault) && fit$ifault != 0L) {
    return(failure(paste0("Lloyd k-means returned nonzero ifault=", fit$ifault, ".")))
  }
  if (length(unique(fit$cluster)) != as.integer(centers) ||
      length(fit$size) != as.integer(centers) || any(fit$size <= 0L)) {
    return(failure("Lloyd k-means returned an empty or missing cluster."))
  }
  attr(fit, "algorithm_used") <- "Lloyd"
  attr(fit, "fallback_used") <- FALSE
  attr(fit, "fallback_reason") <- NA_character_
  attr(fit, "warnings") <- warning_text
  attr(fit, "nstart") <- 25L
  attr(fit, "iter.max") <- 100L
  fit
}

kmeans_audit_record <- function(fit, prefix) {
  if (is.null(fit)) {
    values <- list(
      algorithm_used = "not_run", fallback_used = FALSE,
      fallback_reason = NA_character_, warning_count = 0L,
      warning_text = "", nstart = 25L, itermax = 100L, failed = FALSE
    )
  } else if (inherits(fit, "d1_kmeans_error")) {
    warning_text <- fit$warnings
    if (is.null(warning_text)) warning_text <- character()
    values <- list(
      algorithm_used = if (is.null(fit$algorithm_used)) "failed" else fit$algorithm_used,
      fallback_used = if (is.null(fit$fallback_used)) TRUE else fit$fallback_used,
      fallback_reason = if (is.null(fit$fallback_reason)) NA_character_ else fit$fallback_reason,
      warning_count = length(warning_text),
      warning_text = paste(unique(warning_text), collapse = " | "),
      nstart = if (is.null(fit$nstart)) 25L else as.integer(fit$nstart),
      itermax = if (is.null(fit$iter.max)) 100L else as.integer(fit$iter.max),
      failed = TRUE
    )
  } else {
    warning_text <- attr(fit, "warnings")
    if (is.null(warning_text)) warning_text <- character()
    fallback_value <- attr(fit, "fallback_used")
    if (!is.logical(fallback_value) || length(fallback_value) != 1L ||
        is.na(fallback_value)) fallback_value <- NA
    values <- list(
      algorithm_used = attr(fit, "algorithm_used"),
      fallback_used = fallback_value,
      fallback_reason = attr(fit, "fallback_reason"),
      warning_count = length(warning_text),
      warning_text = paste(unique(warning_text), collapse = " | "),
      nstart = as.integer(attr(fit, "nstart")),
      itermax = as.integer(attr(fit, "iter.max")),
      failed = FALSE
    )
  }
  names(values) <- paste0(prefix, "_", names(values))
  data.table::as.data.table(values)
}

cluster_index <- function(X, labels) {
  X <- as.matrix(X)
  overall <- colMeans(X)
  total_ss <- sum(sweep(X, 2L, overall, FUN = "-")^2)
  if (!is.finite(total_ss) || total_ss <= 0) return(NA_real_)
  within_ss <- 0
  for (g in sort(unique(labels))) {
    Xi <- X[labels == g, , drop = FALSE]
    if (!nrow(Xi)) next
    within_ss <- within_ss + sum(sweep(Xi, 2L, colMeans(Xi), FUN = "-")^2)
  }
  within_ss / total_ss
}

nearest_centroid_labels <- function(X, centers) {
  dmat <- vapply(seq_len(nrow(centers)), function(g) {
    rowSums(sweep(X, 2L, centers[g, ], FUN = "-")^2)
  }, numeric(nrow(X)))
  max.col(-as.matrix(dmat), ties.method = "first")
}

fixed_center_cluster_index <- function(X, labels, centers) {
  X <- as.matrix(X)
  overall <- colMeans(X)
  total_ss <- sum(sweep(X, 2L, overall, FUN = "-")^2)
  if (!is.finite(total_ss) || total_ss <= 0) return(NA_real_)
  within_ss <- sum(vapply(seq_len(nrow(X)), function(i) {
    sum((X[i, ] - centers[labels[i], ])^2)
  }, numeric(1)))
  within_ss / total_ss
}

pooled_within_mahalanobis <- function(X, labels, centers,
                                      ridge_multiplier = D1_MAHALANOBIS_RIDGE_MULTIPLIER) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)
  if (!is.numeric(X) || !is.numeric(centers) || any(!is.finite(X)) ||
      any(!is.finite(centers))) {
    stop("Mahalanobis inputs must be finite numeric matrices.", call. = FALSE)
  }
  if (length(labels) != nrow(X) || any(!labels %in% seq_len(nrow(centers)))) {
    stop("Mahalanobis labels do not match X/centers.", call. = FALSE)
  }
  if (nrow(X) <= nrow(centers)) {
    stop("Mahalanobis denominator requires n > K.", call. = FALSE)
  }
  residuals <- X - centers[labels, , drop = FALSE]
  denom <- nrow(X) - nrow(centers)
  pooled_cov <- crossprod(residuals) / denom
  ridge <- ridge_multiplier * mean(diag(pooled_cov))
  if (!is.finite(ridge) || ridge <= 0) {
    stop("Mahalanobis ridge is non-positive or non-finite.", call. = FALSE)
  }
  pooled_cov_ridge <- pooled_cov + diag(ridge, ncol(X))
  difference <- centers[1L, ] - centers[2L, ]
  solved <- tryCatch(
    solve(pooled_cov_ridge, difference),
    error = function(e) qr.solve(pooled_cov_ridge, difference)
  )
  list(
    delta = sqrt(max(0, sum(difference * solved))),
    fisher_direction = as.numeric(solved),
    ridge = ridge
  )
}

domain1_external_diagnostics <- function(X_train, X_evaluation,
                                         evaluation_labels = NULL, seed) {
  X_train <- as.matrix(X_train)
  X_evaluation <- as.matrix(X_evaluation)
  if (ncol(X_train) != ncol(X_evaluation)) {
    stop("Training/evaluation column mismatch.", call. = FALSE)
  }
  if (any(!is.finite(X_train)) || any(!is.finite(X_evaluation))) {
    stop("D1 matrices contain non-finite values.", call. = FALSE)
  }
  if (!is.null(evaluation_labels) &&
      nrow(X_evaluation) != length(evaluation_labels)) {
    stop("Evaluation matrix/label row mismatch.", call. = FALSE)
  }

  observed_fit <- safe_kmeans(X_train, centers = 2L, seed = seed + 2L)
  primary_kmeans_audit <- kmeans_audit_record(observed_fit, "kmeans_primary")
  if (inherits(observed_fit, "d1_kmeans_error")) {
    stop("Training K2 fit failed: ", observed_fit$error, call. = FALSE)
  }
  centres <- observed_fit$centers
  holdout_labels <- nearest_centroid_labels(X_evaluation, centres)

  separation <- pooled_within_mahalanobis(
    X_train, observed_fit$cluster, centres
  )
  direction <- separation$fisher_direction
  if (sum(direction^2) <= 1e-12) direction[1L] <- 1
  discriminant_axis <- as.numeric(X_evaluation %*% direction)

  pc <- stats::prcomp(X_train, center = TRUE, scale. = FALSE)
  holdout_pc <- stats::predict(pc, newdata = X_evaluation)
  if (ncol(holdout_pc) < 5L) {
    stop("BH6 shape contract requires at least five held-out PC projections.",
         call. = FALSE)
  }
  p_raw <- c(
    discriminant_axis = diptest::dip.test(discriminant_axis)$p.value,
    vapply(seq_len(5L), function(j) {
      diptest::dip.test(holdout_pc[, j])$p.value
    }, numeric(1))
  )
  names(p_raw)[-1L] <- paste0("PC", seq_len(5L))
  if (length(p_raw) != 6L || any(!is.finite(p_raw))) {
    stop("BH6 raw shape P values are incomplete or non-finite.", call. = FALSE)
  }
  p_fdr <- stats::p.adjust(p_raw, method = "BH")
  if (length(p_fdr) != 6L || any(!is.finite(p_fdr))) {
    stop("BH6 adjusted shape P values are incomplete or non-finite.", call. = FALSE)
  }

  observed_ci <- fixed_center_cluster_index(
    X_evaluation, holdout_labels, centres
  )
  input_label_agreement <- if (!is.null(evaluation_labels)) {
    tryCatch(
      mclust::adjustedRandIndex(
        as.integer(as.factor(evaluation_labels)), holdout_labels
      ),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  cbind(data.table::data.table(
    n_training = nrow(X_train),
    n_evaluated = nrow(X_evaluation),
    partition_source = "Lloyd K2 fitted on training data; nearest-centroid assignment on held-out evaluation",
    shape_source = "heldout_fisher_discriminant_bh6",
    reference_scope = "processed_space_fixed_k2_n_specific",
    kmeans_algorithm = "Lloyd",
    kmeans_nstart = 25L,
    kmeans_itermax = 100L,
    dip_evaluation_source = "independent/held-out observations only",
    dip_discriminant_p = unname(p_raw["discriminant_axis"]),
    dip_discriminant_p_fdr = unname(p_fdr["discriminant_axis"]),
    dip_pc1_p = unname(p_raw["PC1"]),
    dip_pc2_p = unname(p_raw["PC2"]),
    dip_pc3_p = unname(p_raw["PC3"]),
    dip_pc4_p = unname(p_raw["PC4"]),
    dip_pc5_p = unname(p_raw["PC5"]),
    dip_pc1_p_fdr = unname(p_fdr["PC1"]),
    dip_pc2_p_fdr = unname(p_fdr["PC2"]),
    dip_pc3_p_fdr = unname(p_fdr["PC3"]),
    dip_pc4_p_fdr = unname(p_fdr["PC4"]),
    dip_pc5_p_fdr = unname(p_fdr["PC5"]),
    dip_any_axis_fdr_lt_0_05 = any(p_fdr < D1_DIP_FDR_ALPHA),
    mahalanobis_centroid_delta = separation$delta,
    mahalanobis_ridge = separation$ridge,
    observed_cluster_index = observed_ci,
    refit_vs_input_label_ari_descriptive = input_label_agreement,
    oracle_ari_used_in_verdict = FALSE
  ), primary_kmeans_audit)
}

domain1_diagnostics <- function(X, labels = NULL, seed,
                                evaluation_n = D1_EVALUATION_N,
                                train_fraction = D1_TRAIN_FRACTION, ...) {
  X <- as.matrix(X)
  n_eval_total <- min(nrow(X), as.integer(evaluation_n))
  set.seed(seed)
  idx <- sort(sample.int(nrow(X), n_eval_total))
  n_train <- max(2L, min(n_eval_total - 2L, floor(n_eval_total * train_fraction)))
  set.seed(seed + 1L)
  train_local <- sort(sample.int(n_eval_total, n_train))
  evaluation_local <- setdiff(seq_len(n_eval_total), train_local)
  domain1_external_diagnostics(
    X_train = X[idx[train_local], , drop = FALSE],
    X_evaluation = X[idx[evaluation_local], , drop = FALSE],
    evaluation_labels = if (!is.null(labels)) labels[idx[evaluation_local]] else NULL,
    seed = seed
  )
}

saveRDS_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, tmp)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Could not finalize RDS: ", path, call. = FALSE)
  invisible(path)
}

checkpoint_signature <- function(...) {
  paste(..., sep = "|", collapse = "|")
}

code_bundle_signature <- function(paths) {
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
  hashes <- unname(tools::md5sum(paths))
  paste(paste(paths, hashes, sep = "="), collapse = ";")
}

read_checkpoint_if_valid <- function(path, signature) {
  if (!file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || is.null(attr(obj, "checkpoint_signature")) ||
      !identical(attr(obj, "checkpoint_signature"), signature)) {
    return(NULL)
  }
  obj
}

write_checkpoint <- function(object, path, signature) {
  attr(object, "checkpoint_signature") <- signature
  saveRDS_atomic(object, path)
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (anyDuplicated(names(x))) {
    duplicated_names <- unique(names(x)[duplicated(names(x))])
    stop(
      "Refusing to write duplicate CSV columns: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(data.table::as.data.table(x), tmp)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Could not finalize output: ", path, call. = FALSE)
  invisible(path)
}

hash_manifest <- function(paths, roles, output_path) {
  stopifnot(length(paths) == length(roles))
  exists <- file.exists(paths)
  md5 <- rep(NA_character_, length(paths))
  md5[exists] <- unname(tools::md5sum(paths[exists]))
  out <- data.table::data.table(
    role = roles,
    path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    exists = exists,
    size_bytes = ifelse(exists, file.info(paths)$size, NA_real_),
    modified_time = ifelse(exists, as.character(file.info(paths)$mtime), NA_character_),
    md5 = md5
  )
  write_csv_atomic(out, output_path)
  out
}
