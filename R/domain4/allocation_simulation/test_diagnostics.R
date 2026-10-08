args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(file_arg))
source(file.path(script_dir, "diagnostics.R"))

near <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12))
tests <- character()
check <- function(name, condition) {
  if (!isTRUE(condition)) stop("Unit test failed: ", name)
  tests <<- c(tests, name)
}

A <- c(rep(1L, 3L), rep(0L, 7L))
p <- mean(A)
e <- rep(p, length(A))
w_ate <- target_weights(A, e, "ATE")
check("constant-score ATE pooled ESS identity",
      near(pooled_ess(w_ate) / length(A), 4 * p * (1 - p)))
stabilized <- ifelse(A == 1L, p / e, (1 - p) / (1 - e))
check("constant-score stabilized weights are one; treatment count unchanged",
      all(stabilized == 1) && sum(A) == 3L)

w <- c(1, 2, 4, 3, 5, 6, 9)
arms <- c(rep(1L, 3L), rep(0L, 4L))
rescaled <- w
rescaled[arms == 1L] <- 5 * rescaled[arms == 1L]
check("within-arm rescaling preserves treated ESS",
      near(arm_diagnostics(w[arms == 1L])["ess"],
           arm_diagnostics(rescaled[arms == 1L])["ess"]))
check("within-arm rescaling preserves treated top-1-percent share",
      near(arm_diagnostics(w[arms == 1L])["top1_share"],
           arm_diagnostics(rescaled[arms == 1L])["top1_share"]))
check("within-arm rescaling changes pooled ESS",
      !near(pooled_ess(w), pooled_ess(rescaled)))
check("equal within-arm weights yield arm size",
      near(arm_diagnostics(rep(7, 12))["ess"], 12))
check("positive weights cannot exceed arm size in ESS",
      arm_diagnostics(c(1, 2, 3, 5, 8))["ess"] <= 5)
check("empty treated arm yields uncomputable ESS",
      is.na(arm_diagnostics(numeric(0))["ess"]))
check("zero and invalid weights yield uncomputable ESS",
      all(is.na(c(arm_diagnostics(c(0, 0))["ess"],
                  arm_diagnostics(c(1, 0))["ess"],
                  arm_diagnostics(c(1, Inf))["ess"],
                  arm_diagnostics(c(1, NA_real_))["ess"]))))
check("target-specific weights match definitions",
      near(target_weights(A, e, "ATT"), ifelse(A == 1L, 1, e / (1 - e))) &&
        near(target_weights(A, e, "ATO"), ifelse(A == 1L, 1 - e, e)))

out_dir <- file.path(Sys.getenv("EHR_AUDIT_OUTPUT_ROOT", unset = tempdir()),
                     "domain4_allocation_tests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(data.frame(test = tests, passed = TRUE),
                 file.path(out_dir, "unit_test_results.csv"), row.names = FALSE)
cat(length(tests), "unit tests passed\n")
