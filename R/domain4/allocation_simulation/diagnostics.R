arm_diagnostics <- function(w) {
  n <- length(w)
  missing_result <- c(ess = NA_real_, retention = NA_real_,
                      top1_share = NA_real_)
  if (n == 0L || any(!is.finite(w)) || any(w <= 0)) {
    return(missing_result)
  }
  total <- sum(w)
  squared <- sum(w^2)
  if (!is.finite(total) || total <= 0 || !is.finite(squared) || squared <= 0) {
    return(missing_result)
  }
  k <- as.integer(ceiling(0.01 * n))
  first_top <- n - k + 1L
  partially_sorted <- sort.int(w, partial = first_top)
  ess <- total^2 / squared
  c(ess = ess, retention = ess / n,
    top1_share = sum(partially_sorted[first_top:n]) / total)
}

pooled_ess <- function(w) {
  if (length(w) == 0L || any(!is.finite(w)) || any(w <= 0)) {
    return(NA_real_)
  }
  total <- sum(w)
  squared <- sum(w^2)
  if (!is.finite(total) || total <= 0 || !is.finite(squared) || squared <= 0) {
    return(NA_real_)
  }
  total^2 / squared
}

target_weights <- function(A, e, target) {
  stopifnot(length(A) == length(e), all(A %in% c(0L, 1L)),
            all(is.finite(e)), all(e > 0 & e < 1))
  switch(target,
         ATE = ifelse(A == 1L, 1 / e, 1 / (1 - e)),
         ATT = ifelse(A == 1L, 1, e / (1 - e)),
         ATO = ifelse(A == 1L, 1 - e, e),
         stop("Unknown weighting target: ", target))
}

absolute_smd <- function(x, A, w = rep(1, length(x))) {
  stopifnot(length(x) == length(A), length(w) == length(A))
  treated <- A == 1L
  control <- A == 0L
  if (sum(treated) < 2L || sum(control) < 2L ||
      any(!is.finite(w)) || any(w <= 0)) {
    return(NA_real_)
  }
  denominator <- sqrt((stats::var(x[treated]) + stats::var(x[control])) / 2)
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  abs(stats::weighted.mean(x[treated], w[treated]) -
        stats::weighted.mean(x[control], w[control])) / denominator
}

empirical_overlap <- function(A, e) {
  result <- c(lower = NA_real_, upper = NA_real_, empty = NA_real_,
              below_treated = NA_real_, above_treated = NA_real_,
              below_control = NA_real_, above_control = NA_real_)
  if (!any(A == 1L) || !any(A == 0L)) {
    return(result)
  }
  e1 <- e[A == 1L]
  e0 <- e[A == 0L]
  lower <- max(min(e1), min(e0))
  upper <- min(max(e1), max(e0))
  result[c("lower", "upper", "empty")] <- c(lower, upper, as.numeric(lower > upper))
  if (lower <= upper) {
    result[c("below_treated", "above_treated", "below_control", "above_control")] <-
      c(mean(e1 < lower), mean(e1 > upper), mean(e0 < lower), mean(e0 > upper))
  }
  result
}
