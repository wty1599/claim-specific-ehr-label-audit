options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
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
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
aggregate_path <- file.path(out, "original_parameters_aggregate_v1.rds")
stopifnot(file.exists(aggregate_path))
old_path <- file.path(output_root, "domain3_deployable_generator_v1_20260728/model",
                      "ehr_audit_k2_deployable_surrogate_generator_v1.rds")
stopifnot(tolower(unname(tools::md5sum(old_path))) ==
            "4722d24103151ccc5bcbcc737772f249")
old <- readRDS(old_path)
original <- readRDS(aggregate_path)
stopifnot(identical(original$object_type,
                    "recovered_original_centroids_with_inductive_median_adapter"),
          identical(original$verification$fixed_label_matches, 20049L),
          identical(unname(original$source_sha256[c("matrix", "labels", "feature_sets", "mice")]),
                    c("9896AD4906207CA910E54EB187507037106BFF212299586141908F387CFE6279",
                      "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
                      "8BC1AABF3ACB88BD6C37F8A857A308ABE92955C65264249D7F31F4827F73C8DF",
                      "B9EF638A9BBD3AA6672D15F788187A70765197B11279A2C1A4A721B0BE87EA80")))
features <- original$feature_names
stopifnot(length(features) == 33L,
          identical(as.character(old$feature_names), features),
          identical(as.character(old$recipe$feature), features))

X <- as.data.table(readRDS(file.path(output_root, "model/X_primary33_std_mice.rds")))
L <- as.data.table(readRDS(file.path(output_root, "model/labels_primary_mice.rds")))
M <- fread(file.path(base, "final_full.csv"),
           select = c("stay_id", features), showProgress = FALSE)
stopifnot(nrow(X) == 20049L, nrow(L) == 20049L, nrow(M) == 20049L,
          !anyDuplicated(M$stay_id), identical(X$stay_id, L$stay_id))
ii <- match(L$stay_id, M$stay_id)
stopifnot(!anyNA(ii))
M <- M[ii]
stopifnot(identical(M$stay_id, L$stay_id))
locked <- as.integer(L$cluster_k2)

stopifnot(max(abs(as.numeric(old$recipe$winsor_p01) -
                  original$winsor_limits[, 1L])) == 0,
          max(abs(as.numeric(old$recipe$winsor_p99) -
                  original$winsor_limits[, 2L])) == 0,
          max(abs(as.numeric(old$recipe$imputation_median) -
                  original$inductive_median_from_pre_imputation_winsorized_observed)) == 0)

assign_fixed <- function(Z, centroids) {
  stopifnot(identical(colnames(Z), colnames(centroids)),
            identical(rownames(centroids), c("1", "2")),
            all(is.finite(Z)), all(is.finite(centroids)))
  a <- rowSums(sweep(Z, 2L, centroids["1", ], "-")^2)
  b <- rowSums(sweep(Z, 2L, centroids["2", ], "-")^2)
  ifelse(a <= b, 1L, 2L)
}
apply_rule <- function(df, features, lower, upper, median_value,
                       center, spread, centroids) {
  stopifnot(identical(names(df)[-1L], features),
            length(lower) == length(features),
            length(upper) == length(features),
            length(median_value) == length(features),
            length(center) == length(features),
            length(spread) == length(features), all(spread > 0))
  Z <- matrix(NA_real_, nrow = nrow(df), ncol = length(features),
              dimnames = list(NULL, features))
  for (j in seq_along(features)) {
    v <- suppressWarnings(as.numeric(df[[features[j]]]))
    absent <- !is.finite(v)
    v <- pmin(pmax(v, lower[j]), upper[j])
    v[absent] <- median_value[j]
    Z[, j] <- (v - center[j]) / spread[j]
  }
  assign_fixed(Z, centroids)
}

Z0 <- as.matrix(X[, ..features])
L0 <- assign_fixed(Z0, original$original_centroids)
stopifnot(identical(as.integer(L0), locked))
L1 <- apply_rule(
  M, features, original$winsor_limits[, 1L], original$winsor_limits[, 2L],
  original$inductive_median_from_pre_imputation_winsorized_observed,
  unname(original$original_center), unname(original$original_scale),
  original$original_centroids
)
L2 <- apply_rule(
  M, features, old$recipe$winsor_p01, old$recipe$winsor_p99,
  old$recipe$imputation_median, old$recipe$reference_mean,
  old$recipe$reference_sd, old$centroids
)
stopifnot(sum(L2 == locked) == 19767L, sum(L2 != locked) == 282L)

archived_discordance_path <- file.path(
  output_root, "domain3_deployable_generator_v1_20260728/private_local",
  "internal_discordant_records_local.csv"
)
stopifnot(file.exists(archived_discordance_path))
archived_discordance <- fread(archived_discordance_path, showProgress = FALSE)
stopifnot(all(c("stay_id", "locked_cluster", "regenerated_cluster") %in%
                names(archived_discordance)),
          !anyDuplicated(archived_discordance$stay_id))
archived_idx <- match(archived_discordance$stay_id, L$stay_id)
stopifnot(!anyNA(archived_idx))
archive_check <- data.table(
  archived_record_n = nrow(archived_discordance),
  reconstructed_discordant_n = sum(L2 != locked),
  discordant_id_sets_identical = setequal(
    archived_discordance$stay_id, L$stay_id[L2 != locked]
  ),
  archived_locked_labels_match = sum(
    as.integer(archived_discordance$locked_cluster) == locked[archived_idx]
  ),
  archived_generated_labels_match = sum(
    as.integer(archived_discordance$regenerated_cluster) == L2[archived_idx]
  ),
  archive_sha256 = sha(archived_discordance_path)
)
stopifnot(archive_check$archived_record_n == 282L,
          archive_check$discordant_id_sets_identical,
          archive_check$archived_locked_labels_match == 282L,
          archive_check$archived_generated_labels_match == 282L)
fwrite(archive_check, file.path(out, "L2_archived_record_reconciliation.csv"))

ari <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(n) n * (n - 1) / 2
  total_pairs <- choose2(sum(tab))
  observed <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- (row_pairs + col_pairs) / 2
  (observed - expected) / (maximum - expected)
}
fidelity <- rbindlist(lapply(list(L0 = L0, L1 = L1, L2 = L2),
                            function(pred) {
  data.table(n = length(pred), concordant_n = sum(pred == locked),
             discordant_n = sum(pred != locked),
             concordance = mean(pred == locked), ari = ari(locked, pred),
             c1_n = sum(pred == 1L), c1_prevalence = mean(pred == 1L))
}), idcol = "layer")
fwrite(fidelity, file.path(out, "L0_L1_L2_fidelity.csv"))

flip01 <- L0 != L1
flip12 <- L1 != L2
transition <- data.table(
  classification = c("neither_step", "L0_to_L1_only", "L1_to_L2_only",
                     "both_steps_cancel", "all_L0_to_L1_flips",
                     "all_L1_to_L2_flips", "final_L0_to_L2_discordance"),
  n = c(sum(!flip01 & !flip12), sum(flip01 & !flip12),
        sum(!flip01 & flip12), sum(flip01 & flip12),
        sum(flip01), sum(flip12), sum(L0 != L2))
)
stopifnot(transition[classification %in% c("L0_to_L1_only", "L1_to_L2_only"),
                     sum(n)] == 282L,
          transition[classification %in% c("neither_step", "L0_to_L1_only",
                                          "L1_to_L2_only", "both_steps_cancel"),
                     sum(n)] == 20049L)
fwrite(transition, file.path(out, "L0_L1_L2_transition_accounting.csv"))

repeat_rows <- vector("list", 10L)
for (r in seq_len(10L)) {
  set.seed(20260728L + 1000L + r)
  fold <- integer(nrow(M))
  for (k in 1:2) {
    members <- which(locked == k)
    fold[members] <- sample(rep_len(1:5, length(members)))
  }
  split_prediction <- integer(nrow(M))
  for (f in 1:5) {
    train <- fold != f
    test <- fold == f
    low <- high <- med <- numeric(length(features))
    for (j in seq_along(features)) {
      x <- suppressWarnings(as.numeric(M[[features[j]]][train]))
      q <- as.numeric(stats::quantile(x, c(0.01, 0.99), na.rm = TRUE))
      x <- pmin(pmax(x, q[1]), q[2])
      low[j] <- q[1]
      high[j] <- q[2]
      med[j] <- stats::median(x, na.rm = TRUE)
    }
    split_prediction[test] <- apply_rule(
      M[test], features, low, high, med,
      unname(original$original_center), unname(original$original_scale),
      original$original_centroids
    )
  }
  stopifnot(all(split_prediction %in% 1:2))
  repeat_rows[[r]] <- data.table(
    repeat_id = r,
    agreement_with_locked = mean(split_prediction == locked),
    ari_with_locked = ari(locked, split_prediction),
    agreement_with_full_L1 = mean(split_prediction == L1),
    ari_with_full_L1 = ari(L1, split_prediction),
    c1_n = sum(split_prediction == 1L)
  )
}
split_check <- rbindlist(repeat_rows)
fwrite(split_check, file.path(out, "split_recipe_perturbation_repeats.csv"))
split_summary <- rbindlist(lapply(setdiff(names(split_check), "repeat_id"),
                                 function(metric) {
  v <- split_check[[metric]]
  data.table(metric = metric, mean = mean(v), sd = stats::sd(v),
             min = min(v), max = max(v))
}))
fwrite(split_summary, file.path(out, "split_recipe_perturbation_summary.csv"))

writeLines(capture.output(sessionInfo()), file.path(out, "decomposition_sessionInfo.txt"))
cat("DECOMPOSITION_PASS L1 discordance", sum(L1 != locked),
    "L2 discordance", sum(L2 != locked), "\n")
