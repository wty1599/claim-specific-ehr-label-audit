args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) args[[1]] else Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/domain3/generator/00_d3_generator_utils_v1.R"))
raw <- parse(file.path(repo,
  "R/domain3/raw_temperature_completion/scripts/01_restore_and_evaluate_raw.R"),
  keep.source = FALSE)
definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("prepare_old_inputs")), as.list(raw))
stopifnot(length(definition) == 1L)
eval(definition[[1]])
calls <- codetools::findGlobals(prepare_old_inputs, merge = FALSE)$functions
stopifnot("predict_ehr_audit_k2" %in% calls,
          !"predict_saaki_k2" %in% calls,
          is.function(predict_ehr_audit_k2))
# Artificial 33-feature fixture; no clinical records or fitted objects are read.
features <- paste0("feature_", seq_len(33L))
recipe <- data.frame(feature = features, winsor_p01 = -100, winsor_p99 = 100,
  imputation_median = 0, reference_mean = 0, reference_sd = 1)
centroids <- rbind(rep(-1, 33L), rep(1, 33L))
dimnames(centroids) <- list(c("1", "2"), features)
object <- list(feature_names = features, recipe = recipe, centroids = centroids,
  distance_metric = "squared_euclidean", generator_version = "artificial",
  generator_content_hash = "artificial")
fixture <- as.data.frame(rbind(rep(-1, 33L), rep(1, 33L), rep(NA_real_, 33L)))
names(fixture) <- features
fixture$stay_id <- c("artificial_a", "artificial_b", "artificial_c")
assigned <- predict_ehr_audit_k2(object, fixture, id_col = "stay_id")
stopifnot(identical(as.integer(assigned$assigned_cluster), c(1L, 2L, 1L)),
          identical(assigned$stay_id, fixture$stay_id),
          identical(assigned$missing_feature_n, c(0L, 0L, 33L)))
cat("Raw public function interface and 33-feature artificial assignment checks passed.\n")
