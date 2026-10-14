args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[[1]] else Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/domain3/generator/00_d3_generator_utils_v1.R"))
required <- c("apply_frozen_recipe", "assign_nearest_centroid", "parse_age",
              "parse_sex", "fast_auc", "calibration_slope_intercept",
              "cluster_bootstrap_indices", "percentile_interval")
stopifnot(all(vapply(required, function(name) is.function(get(name)), logical(1))))
checks <- character()
check <- function(name, condition) {
  stopifnot(isTRUE(condition))
  checks <<- c(checks, name)
}
# Small artificial fixtures only; no clinical input or fitted object is loaded.
features <- c("a", "b")
recipe <- data.frame(feature = features, winsor_p01 = c(0, 0),
                     winsor_p99 = c(10, 10), imputation_median = c(5, 5),
                     reference_mean = c(5, 5), reference_sd = c(5, 5))
fixture <- data.frame(a = c(-1, NA, 10), b = c(0, 5, Inf))
prepared <- apply_frozen_recipe(fixture, recipe, features)
expected <- matrix(c(-1, 0, 1, -1, 0, 0), ncol = 2,
                   dimnames = list(NULL, features))
check("fixed recipe clipping, missing-value imputation and scaling",
      identical(prepared$X, expected) &&
      identical(prepared$missing_n, c(0L, 1L, 1L)) &&
      identical(prepared$clipped_n, c(1L, 0L, 0L)))
centroids <- rbind(c(-1, -1), c(1, 1))
colnames(centroids) <- features
rownames(centroids) <- c("1", "2")
assigned <- assign_nearest_centroid(prepared$X, centroids)
check("nearest-centroid assignment and deterministic tie",
      identical(as.integer(assigned$assigned_cluster), c(1L, 1L, 2L)))
check("AUC ordering and ties", identical(fast_auc(c(0, 1), c(0.1, 0.9)), 1) &&
      identical(fast_auc(c(0, 1), c(0.5, 0.5)), 0.5))
indices <- cluster_bootstrap_indices(c("a", "a", "b", "c"), 3L, 123L)
check("cluster-bootstrap fixed seed",
      identical(indices, cluster_bootstrap_indices(c("a", "a", "b", "c"), 3L, 123L)))
check("percentile interval", isTRUE(all.equal(unname(percentile_interval(1:100)),
      unname(quantile(1:100, c(0.025, 0.975))), tolerance = 0)))
inventory <- read.csv(file.path(repo, "R_PACKAGE_REQUIREMENTS.csv"), stringsAsFactors = FALSE)
inventory_paths <- strsplit(inventory$files, "; ", fixed = TRUE)
check("package inventory references existing released files",
      all(file.exists(file.path(repo, unlist(inventory_paths)))) &&
      identical(as.integer(lengths(inventory_paths)), as.integer(inventory$n_scripts)))
raw_entry <- "R/domain3/raw_temperature_completion/scripts/01_restore_and_evaluate_raw.R"
raw_packages <- c("data.table", "digest", "jsonlite", "glmnet", "ranger")
check("Raw completion dependencies are mapped to their entry point",
      all(vapply(raw_packages, function(package) {
        row <- which(inventory$package == package)
        length(row) == 1L && raw_entry %in% inventory_paths[[row]]
      }, logical(1))))
cat(length(checks), "public dependency utility checks passed; artificial fixtures only.\n")
