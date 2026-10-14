## =====================================================================
## Pure helper functions for D3 deployable surrogate generator v1.
## No file is read and no output is written when this file is sourced.
## =====================================================================

require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

hash_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  tolower(unname(tools::md5sum(path)))
}

hash_object <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf, version = 3)
  hash_file(tf)
}

write_atomic_csv <- function(x, path) {
  require_namespace("data.table")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tf <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  data.table::fwrite(x, tf)
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  if (!file.rename(tf, path)) {
    unlink(tf)
    stop("Atomic write failed: ", path, call. = FALSE)
  }
  invisible(path)
}

write_replace_csv <- function(x, path) {
  require_namespace("data.table")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path)
  invisible(path)
}

safe_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    stop("Refusing to overwrite existing model object: ", path, call. = FALSE)
  }
  saveRDS(x, path, version = 3)
  invisible(path)
}

validate_feature_schema <- function(df, features, object_name = "data") {
  missing <- setdiff(features, names(df))
  if (length(missing)) {
    stop(
      object_name, " is missing feature(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  if (anyDuplicated(features)) {
    stop("Feature schema contains duplicates.", call. = FALSE)
  }
  invisible(TRUE)
}

numeric_feature_frame <- function(df, features) {
  out <- df[, features, drop = FALSE]
  for (j in features) {
    out[[j]] <- suppressWarnings(as.numeric(out[[j]]))
  }
  out
}

derive_frozen_recipe <- function(df, features) {
  require_namespace("data.table")
  validate_feature_schema(df, features)
  rows <- lapply(features, function(j) {
    x <- suppressWarnings(as.numeric(df[[j]]))
    if (!any(is.finite(x))) {
      stop("Feature has no finite development values: ", j, call. = FALSE)
    }
    q <- stats::quantile(
      x, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE, type = 7
    )
    xw <- pmin(pmax(x, q[1]), q[2])
    med <- stats::median(xw, na.rm = TRUE)
    mu <- mean(xw, na.rm = TRUE)
    s <- stats::sd(xw, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) {
      stop("Feature has non-positive post-winsor SD: ", j, call. = FALSE)
    }
    data.table::data.table(
      feature = j,
      winsor_p01 = q[1],
      winsor_p99 = q[2],
      imputation_median = med,
      reference_mean = mu,
      reference_sd = s
    )
  })
  data.table::rbindlist(rows)
}

validate_recipe <- function(recipe, features) {
  required <- c(
    "feature", "winsor_p01", "winsor_p99", "imputation_median",
    "reference_mean", "reference_sd"
  )
  missing <- setdiff(required, names(recipe))
  if (length(missing)) {
    stop("Recipe is missing column(s): ", paste(missing, collapse = ", "))
  }
  if (!identical(as.character(recipe$feature), as.character(features))) {
    stop("Recipe feature order does not match locked feature order.")
  }
  numeric_cols <- setdiff(required, "feature")
  recipe_df <- as.data.frame(recipe)
  if (any(!is.finite(as.matrix(recipe_df[, numeric_cols, drop = FALSE])))) {
    stop("Recipe contains non-finite parameters.")
  }
  if (any(recipe$reference_sd <= 0)) {
    stop("Recipe contains non-positive standard deviations.")
  }
  invisible(TRUE)
}

apply_frozen_recipe <- function(df, recipe, features) {
  validate_feature_schema(df, features)
  validate_recipe(recipe, features)
  n <- nrow(df)
  p <- length(features)
  X <- matrix(NA_real_, nrow = n, ncol = p)
  colnames(X) <- features
  missing_n <- integer(n)
  clipped_n <- integer(n)

  for (k in seq_along(features)) {
    j <- features[k]
    r <- recipe[k, , drop = FALSE]
    x <- suppressWarnings(as.numeric(df[[j]]))
    miss <- !is.finite(x)
    missing_n <- missing_n + as.integer(miss)
    low <- !miss & x < r$winsor_p01
    high <- !miss & x > r$winsor_p99
    clipped_n <- clipped_n + as.integer(low | high)
    x <- pmin(pmax(x, r$winsor_p01), r$winsor_p99)
    x[miss] <- r$imputation_median
    X[, k] <- (x - r$reference_mean) / r$reference_sd
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("Frozen preprocessing produced non-finite standardized values.")
  }
  list(X = X, missing_n = missing_n, clipped_n = clipped_n)
}

derive_locked_centroids <- function(X, labels, levels = c("1", "2")) {
  labels <- as.character(labels)
  if (!all(labels %in% levels) || length(unique(labels)) != 2L) {
    stop("Locked labels do not contain exactly the expected two levels.")
  }
  cent <- t(vapply(
    levels,
    function(k) colMeans(X[labels == k, , drop = FALSE]),
    numeric(ncol(X))
  ))
  rownames(cent) <- levels
  colnames(cent) <- colnames(X)
  if (any(!is.finite(cent))) stop("Centroids contain non-finite values.")
  cent
}

squared_distance_matrix <- function(X, centroids) {
  aa <- rowSums(X^2)
  bb <- rowSums(centroids^2)
  d <- outer(aa, bb, "+") - 2 * X %*% t(centroids)
  pmax(d, 0)
}

assign_nearest_centroid <- function(X, centroids) {
  d <- squared_distance_matrix(X, centroids)
  idx <- max.col(-d, ties.method = "first")
  labels <- rownames(centroids)[idx]
  d1 <- d[, 1]
  d2 <- d[, 2]
  data.frame(
    assigned_cluster = labels,
    distance_c1 = d1,
    distance_c2 = d2,
    signed_margin_c1 = d2 - d1,
    absolute_margin = abs(d2 - d1),
    stringsAsFactors = FALSE
  )
}

predict_ehr_audit_k2 <- function(
    object, newdata, id_col = NULL, strict_schema = TRUE,
    return_diagnostics = TRUE) {
  required <- c("feature_names", "recipe", "centroids", "distance_metric")
  if (!is.list(object) || !all(required %in% names(object)) ||
      !identical(object$distance_metric, "squared_euclidean") ||
      length(object$feature_names) != 33L ||
      !identical(colnames(object$centroids), as.character(object$feature_names)) ||
      !identical(rownames(object$centroids), c("1", "2"))) {
    stop("Object is not a recognized EHR audit deployable surrogate generator.")
  }
  if (isTRUE(strict_schema)) {
    validate_feature_schema(newdata, object$feature_names, "newdata")
  }
  prep <- apply_frozen_recipe(
    newdata, object$recipe, object$feature_names
  )
  ans <- assign_nearest_centroid(prep$X, object$centroids)
  ans$missing_feature_n <- prep$missing_n
  ans$clipped_feature_n <- prep$clipped_n
  ans$generator_version <- object$generator_version
  ans$generator_content_hash <- object$generator_content_hash
  if (!is.null(id_col)) {
    if (!id_col %in% names(newdata)) {
      stop("Requested ID column is absent from newdata: ", id_col)
    }
    id_frame <- as.data.frame(newdata)[, id_col, drop = FALSE]
    ans <- cbind(id_frame, ans)
  }
  if (!isTRUE(return_diagnostics)) {
    keep <- c(if (!is.null(id_col)) id_col else character(), "assigned_cluster")
    ans <- ans[, keep, drop = FALSE]
  }
  ans
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  a <- rowSums(tab)
  b <- colSums(tab)
  index <- sum(choose2(tab))
  expected <- sum(choose2(a)) * sum(choose2(b)) / choose2(n)
  max_index <- 0.5 * (sum(choose2(a)) + sum(choose2(b)))
  denom <- max_index - expected
  if (abs(denom) < .Machine$double.eps) return(NA_real_)
  (index - expected) / denom
}

fidelity_metrics <- function(reference, assigned) {
  reference <- as.character(reference)
  assigned <- as.character(assigned)
  ok <- !is.na(reference) & !is.na(assigned)
  reference <- reference[ok]
  assigned <- assigned[ok]
  lev <- sort(unique(c(reference, assigned)))
  if (!all(c("1", "2") %in% lev)) {
    stop("Fidelity metrics require labels 1 and 2.")
  }
  tab <- table(
    reference = factor(reference, levels = c("1", "2")),
    assigned = factor(assigned, levels = c("1", "2"))
  )
  n <- sum(tab)
  agreement_n <- sum(diag(tab))
  recall_c1 <- tab["1", "1"] / sum(tab["1", ])
  recall_c2 <- tab["2", "2"] / sum(tab["2", ])
  precision_c1 <- tab["1", "1"] / sum(tab[, "1"])
  precision_c2 <- tab["2", "2"] / sum(tab[, "2"])
  data.frame(
    n = n,
    agreement_n = agreement_n,
    discordant_n = n - agreement_n,
    agreement = agreement_n / n,
    ari = adjusted_rand_index(reference, assigned),
    recall_c1 = recall_c1,
    recall_c2 = recall_c2,
    precision_c1 = precision_c1,
    precision_c2 = precision_c2,
    reference_c1_prevalence = mean(reference == "1"),
    assigned_c1_prevalence = mean(assigned == "1"),
    prevalence_difference = mean(assigned == "1") - mean(reference == "1")
  )
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

clamp_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

calibration_slope_intercept <- function(y, p) {
  lp <- stats::qlogis(clamp_probability(p))
  slope <- tryCatch(
    unname(stats::coef(stats::glm(y ~ lp, family = stats::binomial()))[2]),
    error = function(e) NA_real_
  )
  intercept <- tryCatch(
    unname(stats::coef(stats::glm(
      y ~ 1, offset = lp, family = stats::binomial()
    ))[1]),
    error = function(e) NA_real_
  )
  c(slope = slope, intercept = intercept)
}

cluster_bootstrap_indices <- function(cluster_id, B, seed) {
  set.seed(seed)
  cluster_id <- as.character(cluster_id)
  clusters <- unique(cluster_id)
  rows <- split(seq_along(cluster_id), cluster_id)
  lapply(seq_len(B), function(i) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    unlist(rows[sampled], use.names = FALSE)
  })
}

percentile_interval <- function(x, probs = c(0.025, 0.975)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
}

parse_sex <- function(x) {
  if (is.numeric(x)) return(as.integer(x == 1))
  y <- toupper(substr(trimws(as.character(x)), 1, 1))
  ifelse(y == "M", 1L, ifelse(y == "F", 0L, NA_integer_))
}

parse_age <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- gsub(">\\s*89", "90", as.character(x))
  suppressWarnings(as.numeric(y))
}
