options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(mclust)
  library(digest)
})

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(file_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/"))
root <- dirname(script_dir)
private_dir <- file.path(root, "private_not_for_release")
output_dir <- file.path(root, "outputs")
dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

spec <- file.path(root, "01_LOCKED_POSTHOC_SPEC.md")
locked_spec_sha <- "95A3A1FF8ACA21B7992AB7438E90E65C5ABFCDB36F4AD4393A9B429953C7115C"
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
stopifnot(sha(spec) == locked_spec_sha)

project <- Sys.getenv("EHR_AUDIT_RESTRICTED_DATA_ROOT", unset = "")
stopifnot(nzchar(project), dir.exists(project))
paths <- c(
  matrix = file.path(project, "output/model/X_primary33_std_mice.rds"),
  fixed_labels = file.path(project, "output/model/labels_primary_mice.rds"),
  final_full = file.path(project, "data/final_full.csv"),
  risk_ledger = file.path(project, "d4_landmark_support_interface_20260729",
                          "output_formal_v1/private_not_for_release/D4_landmark_patient_ledger_private.rds"),
  registry = file.path(project, "domain1_renal_identity_and_method_closure_20260717",
                       "provenance/D1_feature_domain_registry.csv"),
  seed_metrics = file.path(project, "domain1_renal_identity_and_method_closure_20260717",
                           "tables/Table_D1_ablation_seed_stability.csv")
)
expected_sha <- c(
  matrix = "9896AD4906207CA910E54EB187507037106BFF212299586141908F387CFE6279",
  fixed_labels = "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
  final_full = "CE0B38D2D3D6263A5C4615A6AEA1E78B19840BD4B35046DBFA65EE6BFBBA633B",
  risk_ledger = "233EBC8593B6789B04A7F84C66BDBC19BBCEC503D08573D4892C8FAB7C0486A7",
  registry = "A20D2EAC3C698654E5724E3AE5C1C49BAABC58C1C08E72C401963D8894631673",
  seed_metrics = "888DE81D3C67897D0A501907938B249F060B2D731830D09548B12CCCFB02B6FC"
)
stopifnot(all(file.exists(paths)))
before_sha <- vapply(paths, sha, character(1))
stopifnot(identical(before_sha, expected_sha))

matrix_dt <- as.data.table(readRDS(paths[["matrix"]]))
fixed_dt <- as.data.table(readRDS(paths[["fixed_labels"]]))
full_dt <- fread(paths[["final_full"]], select = c("stay_id", "mortality_30d"))
risk_dt <- as.data.table(readRDS(paths[["risk_ledger"]]))
registry <- fread(paths[["registry"]])
seed_metrics <- fread(paths[["seed_metrics"]])[variant_id == "A4"]

stopifnot(nrow(matrix_dt) == 20049L, ncol(matrix_dt) == 34L,
          nrow(fixed_dt) == 20049L, nrow(risk_dt) == 17525L,
          uniqueN(matrix_dt$stay_id) == nrow(matrix_dt),
          uniqueN(fixed_dt$stay_id) == nrow(fixed_dt),
          uniqueN(full_dt$stay_id) == nrow(full_dt),
          uniqueN(risk_dt$stay_id) == nrow(risk_dt),
          identical(matrix_dt$stay_id, fixed_dt$stay_id),
          all(fixed_dt$cluster_k2 %in% 1:2),
          sum(fixed_dt$cluster_k2 == 1L) == 3992L)

features33 <- setdiff(names(matrix_dt), "stay_id")
stopifnot(length(features33) == 33L, setequal(features33, registry$feature),
          all(vapply(matrix_dt[, ..features33], is.numeric, logical(1))))
strict_renal <- registry[strict_renal_core == 1L, feature]
strict_acid <- registry[strict_acid_base == 1L, feature]
excluded <- union(strict_renal, strict_acid)
features_a4 <- setdiff(features33, excluded)
stopifnot(length(excluded) == 7L, length(features_a4) == 26L,
          setequal(excluded, c("urine_output_24h_ml", "bun_max", "creatinine_max",
                               "ph_min", "pco2_max", "aniongap_max", "bicarbonate_min")))

x_all <- as.matrix(matrix_dt[, ..features33])
x_a4 <- as.matrix(matrix_dt[, ..features_a4])
stopifnot(all(is.finite(x_all)), all(is.finite(x_a4)))
reference <- as.integer(fixed_dt$cluster_k2)
seed <- 20240601L
seed_set <- c(seed, 2026071730L + 0:18)
stopifnot(nrow(seed_metrics) == 20L, setequal(seed_metrics$seed, seed_set),
          seed_metrics[seed == 20240601L, .N] == 1L)

set.seed(seed)
a0 <- kmeans(x_all, centers = 2L, nstart = 25L, iter.max = 100L, algorithm = "Lloyd")
a0_aligned <- if (mean(a0$cluster == reference) >= mean((3L - a0$cluster) == reference))
  as.integer(a0$cluster) else as.integer(3L - a0$cluster)
stopifnot(identical(a0_aligned, reference),
          abs(mclust::adjustedRandIndex(a0_aligned, reference) - 1) < 1e-15)

set.seed(seed)
a4 <- kmeans(x_a4, centers = 2L, nstart = 25L, iter.max = 100L, algorithm = "Lloyd")
raw_a4 <- as.integer(a4$cluster)
ari <- mclust::adjustedRandIndex(reference, raw_a4)
archive_ari <- seed_metrics[seed == 20240601L, ari]
stopifnot(length(archive_ari) == 1L, abs(ari - archive_ari) < 1e-12)

full_idx <- match(matrix_dt$stay_id, full_dt$stay_id)
stopifnot(!anyNA(full_idx))
death <- full_dt$mortality_30d[full_idx]
stopifnot(!anyNA(death), all(death %in% 0:1))
death_raw1 <- mean(death[raw_a4 == 1L])
death_raw2 <- mean(death[raw_a4 == 2L])
stopifnot(is.finite(death_raw1), is.finite(death_raw2), death_raw1 != death_raw2)
a4_oriented <- if (death_raw1 > death_raw2) raw_a4 else 3L - raw_a4

private_labels <- data.table(stay_id = matrix_dt$stay_id,
                             cluster_k2 = reference,
                             ablation_group = a4_oriented)
setkey(private_labels, stay_id)
risk_idx <- match(risk_dt$stay_id, private_labels$stay_id)
stopifnot(!anyNA(risk_idx),
          all(as.integer(risk_dt$cluster_k2) == private_labels$cluster_k2[risk_idx]))
risk_groups <- private_labels$ablation_group[risk_idx]
stopifnot(all(tabulate(risk_groups, nbins = 2L) >= 500L))

ari_rank <- rank(seed_metrics$ari, ties.method = "min")
seed_rank <- ari_rank[which(seed_metrics$seed == seed)]
median_ari <- median(seed_metrics$ari)
selection <- data.table(
  ablation_variant = "A4 strict renal-plus-acid-base block removed",
  selection_rule = "Original first seed; patient-level per-seed labels not in located Step 0 archive",
  seed = seed, n_seeds_archived = nrow(seed_metrics),
  ari_against_locked_k2 = ari, archived_ari = archive_ari,
  ari_ascending_rank_of_20 = seed_rank, ari_median_of_20 = median_ari,
  n_features = length(features_a4),
  orientation = if (death_raw1 > death_raw2) "raw1=A1" else "raw2=A1",
  death_rate_A1_full = mean(death[a4_oriented == 1L]),
  death_rate_A2_full = mean(death[a4_oriented == 2L])
)
feature_audit <- registry[, .(feature, strict_renal_core, strict_acid_base)]
feature_audit[, included_in_a4 := feature %chin% features_a4]
feature_audit[, original_matrix_position := match(feature, features33)]
setorder(feature_audit, original_matrix_position)

full_cross <- as.data.table(table(K2 = reference, A4 = a4_oriented))
full_cross[, population := "Full 20,049 cohort"]
risk_cross <- as.data.table(table(K2 = private_labels$cluster_k2[risk_idx], A4 = risk_groups))
risk_cross[, population := "Hour-24 risk set"]
cross <- rbindlist(list(full_cross, risk_cross), use.names = TRUE)
setcolorder(cross, c("population", "K2", "A4", "N"))
group_counts <- rbindlist(list(
  data.table(population = "Full 20,049 cohort", group = c("A1", "A2"),
             n = as.integer(tabulate(a4_oriented, 2L))),
  data.table(population = "Hour-24 risk set", group = c("A1", "A2"),
             n = as.integer(tabulate(risk_groups, 2L)))
))

after_sha <- vapply(paths, sha, character(1))
stopifnot(identical(before_sha, after_sha))
source_manifest <- data.table(role = names(paths), path = unname(paths),
                              sha256_before = unname(before_sha),
                              sha256_after = unname(after_sha), unchanged = TRUE)

saveRDS(private_labels, file.path(private_dir, "A4_first_seed_labels_private.rds"))
fwrite(selection, file.path(output_dir, "A4_seed_selection.csv"))
fwrite(feature_audit, file.path(output_dir, "A4_feature_audit.csv"))
fwrite(group_counts, file.path(output_dir, "A4_group_counts.csv"))
fwrite(cross, file.path(output_dir, "K2_A4_cross_tab.csv"))
fwrite(source_manifest, file.path(output_dir, "source_manifest_preanalysis.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_session_label_construction.txt"))
cat("A4 first-seed label construction passed; ARI=", format(ari, digits = 16),
    "; risk group sizes=", paste(tabulate(risk_groups, 2L), collapse = "/"), "\n", sep = "")
