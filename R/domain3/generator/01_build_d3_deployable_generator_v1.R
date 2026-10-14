## =====================================================================
## Build and freeze the MIMIC-only regenerated deployable surrogate.
## This script must not read eICU data or any external outcome.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))
require_namespace("data.table")
library(data.table)

set.seed(MASTER_SEED)

build_inputs <- c(
  PATH_FINAL_FULL, PATH_LOCKED_LABELS,
  PATH_REFERENCE_RECIPE, PATH_REFERENCE_CENTROIDS
)
missing_inputs <- build_inputs[!file.exists(build_inputs)]
if (length(missing_inputs)) {
  stop("Missing authoritative input(s): ", paste(missing_inputs, collapse = "; "))
}

fingerprints <- data.table(
  role = c(
    "final_full", "locked_labels", "reference_recipe",
    "reference_centroids"
  ),
  path = c(
    PATH_FINAL_FULL, PATH_LOCKED_LABELS,
    PATH_REFERENCE_RECIPE, PATH_REFERENCE_CENTROIDS
  ),
  expected_md5 = unname(LOCKED_MD5[c(
    "final_full", "locked_labels", "reference_recipe", "reference_centroids"
  )])
)
fingerprints[, observed_md5 := vapply(path, hash_file, character(1))]
fingerprints[, match := observed_md5 == expected_md5]
write_replace_csv(
  fingerprints, file.path(OUT_QC, "generator_v1_input_fingerprints.csv")
)
if (any(!fingerprints$match)) {
  stop(
    "Authoritative input fingerprint mismatch: ",
    paste(fingerprints[!match, role], collapse = ", ")
  )
}

locked <- as.data.table(readRDS(PATH_LOCKED_LABELS))
if (!all(c("stay_id", "cluster_k2") %in% names(locked))) {
  stop("Locked label artifact lacks stay_id or cluster_k2.")
}
if (nrow(locked) != EXPECTED_MIMIC_N || anyDuplicated(locked$stay_id)) {
  stop("Locked label denominator or uniqueness check failed.")
}
locked[, cluster_k2 := as.character(cluster_k2)]
if (!identical(sort(unique(locked$cluster_k2)), SOURCE_LABEL_LEVELS)) {
  stop("Locked label levels differ from the protected 1/2 coding.")
}

M <- fread(
  PATH_FINAL_FULL,
  select = c("stay_id", "subject_id", FEATURES33),
  showProgress = FALSE
)
if (nrow(M) != EXPECTED_MIMIC_N || anyDuplicated(M$stay_id)) {
  stop("MIMIC development denominator or stay_id uniqueness check failed.")
}
idx <- match(locked$stay_id, M$stay_id)
if (anyNA(idx)) stop("Some locked-label stay IDs are absent from final_full.")
M <- M[idx]
if (!identical(as.integer(M$stay_id), as.integer(locked$stay_id))) {
  stop("MIMIC rows could not be aligned exactly to locked labels.")
}
validate_feature_schema(M, FEATURES33, "MIMIC development data")

recipe <- derive_frozen_recipe(M, FEATURES33)
prep <- apply_frozen_recipe(M, recipe, FEATURES33)
centroids <- derive_locked_centroids(
  prep$X, locked$cluster_k2, SOURCE_LABEL_LEVELS
)

feature_schema <- data.table(
  feature_order = seq_along(FEATURES33),
  feature = FEATURES33,
  feature_order_hash = hash_object(FEATURES33)
)

build_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
generator_core <- list(
  object_class = "ehr_audit_k2_deployable_surrogate_generator",
  generator_name = GENERATOR_NAME,
  generator_version = GENERATOR_VERSION,
  specification_date = "2026-07-28",
  post_hoc_status = GENERATOR_STATUS,
  source_label_path = normalizePath(
    PATH_LOCKED_LABELS, winslash = "/", mustWork = TRUE
  ),
  source_label_md5 = hash_file(PATH_LOCKED_LABELS),
  source_data_path = normalizePath(
    PATH_FINAL_FULL, winslash = "/", mustWork = TRUE
  ),
  source_data_md5 = hash_file(PATH_FINAL_FULL),
  training_id_hash = hash_object(as.integer(M$stay_id)),
  feature_names = FEATURES33,
  feature_order_hash = hash_object(FEATURES33),
  recipe = recipe,
  recipe_hash = hash_object(recipe),
  centroids = centroids,
  centroid_hash = hash_object(centroids),
  label_levels = SOURCE_LABEL_LEVELS,
  label_orientation = list(
    higher_risk_label = SOURCE_HIGHER_RISK_LABEL,
    rule = "preserved from locked cluster_k2 coding; not re-estimated"
  ),
  missingness_policy = MISSINGNESS_POLICY,
  distance_metric = DISTANCE_METRIC,
  assignment_rule = paste(
    "frozen MIMIC 1st/99th winsorization;",
    "frozen post-winsor medians; frozen means/SDs;",
    "nearest locked-label centroid"
  ),
  source_object_note = paste(
    "Post hoc regenerated deployable surrogate.",
    "Not the exact original MICE-imputation-1 generator."
  ),
  master_seed = MASTER_SEED,
  r_version = R.version.string,
  package_versions = c(
    data.table = as.character(utils::packageVersion("data.table"))
  )
)
generator_core$generator_content_hash <- hash_object(generator_core)
class(generator_core) <- c(
  "ehr_audit_k2_deployable_surrogate_generator", "list"
)

reference_recipe <- fread(PATH_REFERENCE_RECIPE)
reference_recipe <- reference_recipe[match(FEATURES33, feature)]
reference_centroids <- fread(PATH_REFERENCE_CENTROIDS)
if (!"source_label" %in% names(reference_centroids)) {
  stop("Reference centroid file lacks source_label.")
}
reference_centroids <- reference_centroids[
  match(SOURCE_LABEL_LEVELS, as.character(source_label))
]
reference_centroid_matrix <- as.matrix(reference_centroids[, ..FEATURES33])
rownames(reference_centroid_matrix) <- SOURCE_LABEL_LEVELS
storage.mode(reference_centroid_matrix) <- "double"

recipe_numeric <- c(
  "winsor_p01", "winsor_p99", "imputation_median", "reference_mean",
  "reference_sd"
)
recipe_diff <- max(abs(
  as.matrix(recipe[, ..recipe_numeric]) -
    as.matrix(reference_recipe[, ..recipe_numeric])
))
centroid_diff <- max(abs(centroids - reference_centroid_matrix))

regression_comparison <- data.table(
  component = c("recipe", "centroids"),
  max_abs_difference = c(recipe_diff, centroid_diff),
  tolerance = c(1e-10, 1e-10),
  pass = c(recipe_diff <= 1e-10, centroid_diff <= 1e-10),
  reference_path = c(PATH_REFERENCE_RECIPE, PATH_REFERENCE_CENTROIDS)
)
write_replace_csv(
  regression_comparison,
  file.path(OUT_QC, "generator_v1_regression_comparison.csv")
)
if (any(!regression_comparison$pass)) {
  stop(
    "New generator does not reproduce the existing surrogate recipe/centroids."
  )
}

if (file.exists(PATH_GENERATOR)) {
  stop("Generator output already exists; refusing to overwrite: ", PATH_GENERATOR)
}
saveRDS(generator_core, PATH_GENERATOR, version = 3)

generator_file_md5 <- hash_file(PATH_GENERATOR)
manifest <- data.table(
  generator_name = GENERATOR_NAME,
  generator_version = GENERATOR_VERSION,
  generator_status = GENERATOR_STATUS,
  generator_path = PATH_GENERATOR,
  generator_file_md5 = generator_file_md5,
  generator_content_hash = generator_core$generator_content_hash,
  source_label_md5 = generator_core$source_label_md5,
  source_data_md5 = generator_core$source_data_md5,
  training_id_hash = generator_core$training_id_hash,
  feature_order_hash = generator_core$feature_order_hash,
  recipe_hash = generator_core$recipe_hash,
  centroid_hash = generator_core$centroid_hash,
  n_training = nrow(M),
  feature_n = length(FEATURES33),
  external_data_read = FALSE,
  external_outcome_read = FALSE,
  build_utc = build_utc
)

write_replace_csv(recipe, file.path(OUT_QC, "generator_v1_recipe.csv"))
write_replace_csv(
  as.data.table(centroids, keep.rownames = "source_label"),
  file.path(OUT_QC, "generator_v1_centroids.csv")
)
write_replace_csv(
  feature_schema, file.path(OUT_QC, "generator_v1_feature_schema.csv")
)
write_replace_csv(
  manifest, file.path(OUT_QC, "generator_v1_build_manifest.csv")
)

yml <- c(
  paste0("generator_name: ", GENERATOR_NAME),
  paste0("generator_version: ", GENERATOR_VERSION),
  paste0("status: ", GENERATOR_STATUS),
  paste0("master_seed: ", MASTER_SEED),
  paste0("mimic_n: ", EXPECTED_MIMIC_N),
  paste0("eicu_stays_expected: ", EXPECTED_EICU_STAYS),
  paste0("feature_n: ", length(FEATURES33)),
  paste0("source_label_md5: ", generator_core$source_label_md5),
  paste0("source_data_md5: ", generator_core$source_data_md5),
  paste0("generator_file_md5: ", generator_file_md5),
  paste0("generator_content_hash: ", generator_core$generator_content_hash),
  "external_outcomes_used_for_construction: false",
  "external_reclustering: false",
  "replaces_locked_labels: false"
)
writeLines(
  yml, file.path(OUT_ROOT, "D3_deployable_generator_v1_locked.yml"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_LOGS, "generator_v1_build_sessionInfo.txt"),
  useBytes = TRUE
)

cat("Generator built and frozen:\n", PATH_GENERATOR, "\n")
cat("Generator file MD5:", generator_file_md5, "\n")
cat("Recipe max difference:", format(recipe_diff, scientific = TRUE), "\n")
cat("Centroid max difference:", format(centroid_diff, scientific = TRUE), "\n")
