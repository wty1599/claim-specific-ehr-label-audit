options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(digest)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
output_root <- Sys.getenv("EHR_AUDIT_OUTPUT_ROOT")
base <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT")
stopifnot(nzchar(base), nzchar(repo))
source(file.path(repo, "R/common/01_assert_external_output.R"))
ehr_audit_assert_external_output(output_root, repo)
root <- file.path(output_root, "domain3_parameter_recovery")
out <- file.path(root, "outputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

model <- file.path(output_root, "model")
inputs <- c(
  matrix = file.path(model, "X_primary33_std_mice.rds"),
  labels = file.path(model, "labels_primary_mice.rds"),
  feature_sets = file.path(model, "feature_sets.rds"),
  mice = file.path(model, "mice_primary.rds"),
  final_full = file.path(base, "final_full.csv"),
  baseline = file.path(base, "baseline_covars.csv")
)
stopifnot(all(file.exists(inputs)))
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
expected_sha <- c(
  matrix = "9896AD4906207CA910E54EB187507037106BFF212299586141908F387CFE6279",
  labels = "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
  feature_sets = "8BC1AABF3ACB88BD6C37F8A857A308ABE92955C65264249D7F31F4827F73C8DF",
  mice = "B9EF638A9BBD3AA6672D15F788187A70765197B11279A2C1A4A721B0BE87EA80"
)
observed_sha <- vapply(inputs, sha, character(1))
stopifnot(identical(unname(observed_sha[names(expected_sha)]),
                    unname(expected_sha)))
expected_md5 <- c(
  matrix = "15a6f073e7e6d6c2d611b0b2d04b8f99",
  labels = "e741e9cf6bc39bf8075127a76dab5b2a",
  mice = "9da5ae19f60efaa10d4b54845d897f4e",
  final_full = "a3c6582cc68922787b432632f3aff8d7",
  baseline = "972e5b37ab179908d926c2cd34618129"
)
observed_md5 <- unname(tools::md5sum(inputs))
names(observed_md5) <- names(inputs)
stopifnot(identical(tolower(observed_md5[names(expected_md5)]),
                    expected_md5))

provenance <- data.table(
  role = names(inputs), path = unname(inputs),
  bytes = as.numeric(file.info(inputs)$size),
  modified = format(file.info(inputs)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  sha256 = unname(observed_sha), md5 = unname(observed_md5)
)
fwrite(provenance, file.path(out, "input_provenance.csv"))

X <- as.data.table(readRDS(inputs[["matrix"]]))
locked <- as.data.table(readRDS(inputs[["labels"]]))
fs <- readRDS(inputs[["feature_sets"]])
features <- as.character(fs[["primary"]])
imp <- readRDS(inputs[["mice"]])
stopifnot(length(features) == 33L, !anyDuplicated(features),
          identical(names(X), c("stay_id", features)),
          identical(names(locked), c("stay_id", "cluster_k2")),
          nrow(X) == 20049L, nrow(locked) == 20049L,
          !anyDuplicated(X$stay_id), !anyDuplicated(locked$stay_id),
          identical(X$stay_id, locked$stay_id),
          identical(sort(unique(locked$cluster_k2)), 1:2),
          inherits(imp, "mids"), imp$m >= 1L,
          nrow(imp$data) == 20049L,
          identical(names(imp$data), features))
Z <- as.matrix(X[, ..features])
storage.mode(Z) <- "double"
stopifnot(all(is.finite(Z)))
labels <- as.integer(locked$cluster_k2)
centroids <- rbind(
  `1` = colMeans(Z[labels == 1L, , drop = FALSE]),
  `2` = colMeans(Z[labels == 2L, , drop = FALSE])
)
sq1 <- rowSums(sweep(Z, 2L, centroids["1", ], "-")^2)
sq2 <- rowSums(sweep(Z, 2L, centroids["2", ], "-")^2)
reassigned <- ifelse(sq1 <= sq2, 1L, 2L)
margin <- abs(sq1 - sq2)
fixed_label_check <- data.table(
  n = nrow(Z), matched = sum(reassigned == labels),
  mismatched = sum(reassigned != labels),
  exact_distance_ties = sum(sq1 == sq2),
  near_ties_1e_10 = sum(margin < 1e-10),
  minimum_absolute_distance_margin = min(margin),
  c1_n = sum(labels == 1L), c2_n = sum(labels == 2L)
)
fwrite(fixed_label_check, file.path(out, "fixed_label_recovery_check.csv"))
if (fixed_label_check$mismatched != 0L) {
  stop("Fixed-label centroids did not reproduce all locked labels; stop before coordinate recovery.")
}

# Rebuild only the historical row mapping and pre-imputation winsorized frame.
# The saved mids object is extracted, never re-imputed.
d <- fread(inputs[["final_full"]], showProgress = FALSE)
b <- fread(inputs[["baseline"]], showProgress = FALSE)
needed <- c("stay_id", "age", "gender", "sofa_score", "mortality_30d",
            "make30", features)
stopifnot(all(needed %in% names(d)),
          all(c("stay_id", "aki_stage_0_24h") %in% names(b)),
          !anyDuplicated(d$stay_id), !anyDuplicated(b$stay_id))
raw <- d[, ..needed]
raw[, sex := as.integer(toupper(substr(as.character(gender), 1L, 1L)) == "M")]
raw <- merge(raw, b[, .(stay_id, aki_stage_0_24h)], by = "stay_id", all.x = TRUE)
raw <- raw[complete.cases(raw[, c("age", "sex", "sofa_score",
                                   "aki_stage_0_24h", "mortality_30d",
                                   "make30"), with = FALSE])]
stopifnot(nrow(raw) == 20049L, identical(raw$stay_id, X$stay_id))

limits <- matrix(NA_real_, nrow = length(features), ncol = 2L,
                 dimnames = list(features, c("winsor_p01", "winsor_p99")))
for (i in seq_along(features)) {
  j <- features[i]
  v <- as.numeric(raw[[j]])
  q <- as.numeric(stats::quantile(v, c(0.01, 0.99), na.rm = TRUE))
  v[v < q[1]] <- q[1]
  v[v > q[2]] <- q[2]
  raw[[j]] <- v
  limits[i, ] <- q
}

missing_pattern_match <- identical(is.na(as.matrix(raw[, ..features])),
                                   is.na(as.matrix(imp$data[, features])))
max_observed_difference <- max(vapply(features, function(j) {
  obs <- !is.na(raw[[j]])
  if (!any(obs)) return(0)
  max(abs(as.numeric(raw[[j]][obs]) - as.numeric(imp$data[[j]][obs])))
}, numeric(1)))
stopifnot(missing_pattern_match, max_observed_difference <= 1e-10)
completed <- mice::complete(imp, action = 1L)
stopifnot(nrow(completed) == nrow(raw), identical(names(completed), features))
completed_matrix <- as.matrix(completed[, features, drop = FALSE])
storage.mode(completed_matrix) <- "double"
stopifnot(all(is.finite(completed_matrix)))
z_from_completion <- scale(completed_matrix)
max_matrix_difference <- max(abs(z_from_completion - Z))
center <- attr(z_from_completion, "scaled:center")
spread <- attr(z_from_completion, "scaled:scale")
coordinate_check <- data.table(
  route = "B_saved_MICE_completion_1",
  row_mapping_identical = TRUE,
  missing_pattern_identical = missing_pattern_match,
  max_observed_raw_vs_mids = max_observed_difference,
  max_saved_vs_rescaled_matrix = max_matrix_difference,
  scale_attribute_present_in_saved_table =
    all(c("scaled:center", "scaled:scale") %in% names(attributes(X))),
  feature_n = length(features), n = nrow(raw)
)
fwrite(coordinate_check, file.path(out, "coordinate_recovery_check.csv"))
if (max_matrix_difference > 1e-10 || any(!is.finite(center)) ||
    any(!is.finite(spread)) || any(spread <= 0)) {
  stop("Saved matrix was not reproduced from MICE completion 1; original coordinate system unqualified.")
}

recipe <- data.table(
  feature_order = seq_along(features), feature = features,
  winsor_p01 = limits[, 1L], winsor_p99 = limits[, 2L],
  imputation_median = vapply(features, function(j) median(raw[[j]], na.rm = TRUE), numeric(1)),
  original_center = unname(center[features]),
  original_scale = unname(spread[features]),
  centroid_c1 = unname(centroids["1", features]),
  centroid_c2 = unname(centroids["2", features])
)
numeric_recipe_cols <- setdiff(names(recipe), "feature")
stopifnot(nrow(recipe) == 33L,
          all(is.finite(as.matrix(recipe[, ..numeric_recipe_cols]))))
fwrite(recipe, file.path(out, "recovered_aggregate_parameters.csv"))
aggregate_object <- list(
  object_type = "recovered_original_centroids_with_inductive_median_adapter",
  recovery_date = "2026-09-25",
  route = "B: extract saved MICE completion 1; no imputation refit",
  feature_names = features,
  winsor_limits = limits,
  original_center = center[features],
  original_scale = spread[features],
  original_centroids = centroids,
  inductive_median_from_pre_imputation_winsorized_observed =
    recipe$imputation_median,
  source_sha256 = observed_sha,
  verification = list(
    fixed_label_matches = fixed_label_check$matched,
    max_saved_vs_rescaled_matrix = max_matrix_difference,
    max_observed_raw_vs_mids = max_observed_difference
  )
)
saveRDS(aggregate_object,
        file.path(out, "original_parameters_aggregate_v1.rds"), version = 3)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("RECOVERY_PASS", fixed_label_check$matched, "of", nrow(Z),
    "matrix max abs difference", format(max_matrix_difference, digits = 8), "\n")
