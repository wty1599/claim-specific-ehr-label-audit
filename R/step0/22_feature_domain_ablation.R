## =====================================================================
## 29b_domain1_feature_domain_ablation.R
## Structural dependence of the locked K2 partition on prespecified
## feature domains. Outcomes are descriptive only and never used to align
## or select an ablation result.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "mclust", "ggplot2", "digest")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
  library(ggplot2)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), winslash = "/", mustWork = TRUE)
OUT_ROOT <- normalizePath(file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"), winslash = "/", mustWork = TRUE)
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"), figures = file.path(OUT_ROOT, "figures"),
  logs = file.path(OUT_ROOT, "logs"), checkpoints = file.path(OUT_ROOT, "checkpoints"),
  provenance = file.path(OUT_ROOT, "provenance")
)
INPUT <- list(
  matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  final_full = file.path(PROJECT_ROOT, "data", "final_full.csv"),
  registry = file.path(DIRS$provenance, "D1_feature_domain_registry.csv"),
  label_builder = file.path(PROJECT_ROOT, "r", "23c_save_mice_matrix_labels.R")
)
if (any(!file.exists(unlist(INPUT)))) stop("Missing input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))

LOCKED_SEED <- 20240601L
SEEDS <- c(LOCKED_SEED, 2026071730L + 0:18)
NSTART <- 25L
ITER_MAX <- 100L
ALGORITHM <- "Lloyd"

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, bom = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}
write_lines_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not write ", path)
}
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)

entropy <- function(counts) {
  p <- counts[counts > 0] / sum(counts)
  -sum(p * log(p))
}

adjusted_mutual_information <- function(reference, candidate) {
  tab <- table(reference, candidate)
  n <- sum(tab)
  a <- rowSums(tab)
  b <- colSums(tab)
  mi <- 0
  for (i in seq_along(a)) for (j in seq_along(b)) {
    nij <- tab[i, j]
    if (nij > 0) mi <- mi + (nij / n) * log((n * nij) / (a[i] * b[j]))
  }
  emi <- 0
  for (ai in a) for (bj in b) {
    lower <- max(1L, ai + bj - n)
    upper <- min(ai, bj)
    if (lower <= upper) {
      nij <- lower:upper
      probability <- dhyper(nij, bj, n - bj, ai)
      term <- (nij / n) * log((n * nij) / (ai * bj))
      emi <- emi + sum(probability * term)
    }
  }
  h_ref <- entropy(a)
  h_candidate <- entropy(b)
  denominator <- mean(c(h_ref, h_candidate)) - emi
  if (abs(denominator) < 1e-15) return(ifelse(abs(mi - emi) < 1e-15, 1, NA_real_))
  (mi - emi) / denominator
}

kappa_binary <- function(reference, candidate) {
  observed <- mean(reference == candidate)
  p_ref <- mean(reference == 1L)
  p_new <- mean(candidate == 1L)
  expected <- p_ref * p_new + (1 - p_ref) * (1 - p_new)
  if (abs(1 - expected) < 1e-12) NA_real_ else (observed - expected) / (1 - expected)
}

jaccard <- function(a, b) safe_div(sum(a & b), sum(a | b))

align_to_reference <- function(raw_label, reference) {
  direct <- mean(raw_label == reference)
  swapped <- mean((3L - raw_label) == reference)
  if (direct >= swapped) as.integer(raw_label) else as.integer(3L - raw_label)
}

run_ablation <- function(X_variant, X_all, reference, mortality, make30, seed, variant_id, variant_label, feature_names) {
  set.seed(seed)
  km <- kmeans(
    X_variant, centers = 2L, nstart = NSTART, iter.max = ITER_MAX,
    algorithm = ALGORITHM
  )
  aligned <- align_to_reference(km$cluster, reference)
  ref_centroids <- do.call(rbind, lapply(1:2, function(g) colMeans(X_all[reference == g, , drop = FALSE])))
  new_centroids <- do.call(rbind, lapply(1:2, function(g) colMeans(X_all[aligned == g, , drop = FALSE])))
  centroid_difference <- as.numeric(new_centroids - ref_centroids)
  data.table(
    variant_id = variant_id,
    variant_label = variant_label,
    seed = seed,
    locked_seed = seed == LOCKED_SEED,
    n_features = ncol(X_variant),
    included_features = paste(feature_names, collapse = ";"),
    ari = mclust::adjustedRandIndex(reference, aligned),
    adjusted_mutual_information = adjusted_mutual_information(reference, aligned),
    kappa = kappa_binary(reference, aligned),
    label_switch_proportion = mean(reference != aligned),
    c1_prevalence = mean(aligned == 1L),
    jaccard_c1 = jaccard(reference == 1L, aligned == 1L),
    jaccard_c2 = jaccard(reference == 2L, aligned == 2L),
    centroid_profile_rmse_all33 = sqrt(mean(centroid_difference^2)),
    centroid_profile_correlation_all33 = cor(as.numeric(ref_centroids), as.numeric(new_centroids)),
    mortality_c1 = mean(mortality[aligned == 1L]),
    mortality_c2 = mean(mortality[aligned == 2L]),
    mortality_risk_difference_c1_minus_c2 = mean(mortality[aligned == 1L]) - mean(mortality[aligned == 2L]),
    make30_c1 = mean(make30[aligned == 1L]),
    make30_c2 = mean(make30[aligned == 2L]),
    make30_risk_difference_c1_minus_c2 = mean(make30[aligned == 1L]) - mean(make30[aligned == 2L]),
    outcome_use = "descriptive only; not used for label alignment or structural decisions"
  )
}

export_plot <- function(plot, stem, width = 8.0, height = 5.2) {
  cairo_pdf(file.path(DIRS$figures, paste0(stem, ".pdf")), width = width, height = height, family = "sans")
  print(plot)
  dev.off()
  png(file.path(DIRS$figures, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300, type = "cairo")
  print(plot)
  dev.off()
  tiff(file.path(DIRS$figures, paste0(stem, ".tiff")), width = width, height = height, units = "in", res = 600, compression = "lzw")
  print(plot)
  dev.off()
}

X_dt <- as.data.table(readRDS(INPUT$matrix))
labels <- as.data.table(readRDS(INPUT$labels))
final_full <- fread(INPUT$final_full, showProgress = FALSE)
registry <- fread(INPUT$registry)
if (!identical(X_dt$stay_id, labels$stay_id)) stop("Locked matrix/label order mismatch")
if (nrow(X_dt) != 20049L || sum(labels$cluster_k2 == 1L) != 3992L) stop("Locked denominator changed")
features33 <- setdiff(names(X_dt), "stay_id")
if (!setequal(features33, registry$feature)) stop("Registry mismatch")
X <- as.matrix(X_dt[, ..features33])
reference <- as.integer(labels$cluster_k2)
idx <- match(labels$stay_id, final_full$stay_id)
if (anyNA(idx)) stop("Missing final_full stay_id")
mortality <- as.numeric(final_full$mortality_30d[idx])
make30 <- as.numeric(final_full$make30[idx])

renal_core <- registry[strict_renal_core == 1L, feature]
acid_base <- registry[strict_acid_base == 1L, feature]
renal_acid <- union(renal_core, acid_base)
nonrenal_sofa_proxy <- registry[nonrenal_sofa_proxy == 1L, feature]

variants <- list(
  A0 = list(label = "A0: all 33 features", features = features33),
  A1 = list(label = "A1: remove creatinine + BUN", features = setdiff(features33, c("creatinine_max", "bun_max"))),
  A2 = list(label = "A2: remove renal core", features = setdiff(features33, renal_core)),
  A3 = list(label = "A3: remove acid-base block", features = setdiff(features33, acid_base)),
  A4 = list(label = "A4: remove renal + acid-base", features = setdiff(features33, renal_acid)),
  A5 = list(label = "A5: remove non-renal SOFA proxy block", features = setdiff(features33, nonrenal_sofa_proxy)),
  A6 = list(label = "A6: renal core only", features = renal_core),
  A7 = list(label = "A7: renal + acid-base only", features = renal_acid)
)

cat("A0 locked-seed reproduction check...\n")
a0 <- run_ablation(X, X, reference, mortality, make30, LOCKED_SEED, "A0", variants$A0$label, features33)
if (!isTRUE(all.equal(a0$ari, 1)) || a0$label_switch_proportion != 0) {
  stop("A0 did not reproduce locked labels; ablation analysis stopped")
}

results <- vector("list", length(variants) * length(SEEDS))
position <- 0L
for (variant_id in names(variants)) {
  spec <- variants[[variant_id]]
  X_variant <- X[, spec$features, drop = FALSE]
  cat(sprintf("Running %s with %d features across %d seeds...\n", variant_id, ncol(X_variant), length(SEEDS)))
  for (seed in SEEDS) {
    position <- position + 1L
    results[[position]] <- run_ablation(
      X_variant, X, reference, mortality, make30, seed,
      variant_id, spec$label, spec$features
    )
  }
  saveRDS(results[seq_len(position)], file.path(DIRS$checkpoints, "29b_ablation_checkpoint.rds"))
}
seed_results <- rbindlist(results, fill = TRUE)
primary <- seed_results[locked_seed == TRUE]
summary <- seed_results[, .(
  n_seeds = .N,
  n_features = first(n_features),
  included_features = first(included_features),
  ari_primary_seed = ari[seed == LOCKED_SEED],
  ari_mean = mean(ari), ari_sd = sd(ari), ari_p025 = quantile(ari, 0.025, type = 6), ari_p975 = quantile(ari, 0.975, type = 6),
  ami_mean = mean(adjusted_mutual_information), ami_sd = sd(adjusted_mutual_information),
  kappa_mean = mean(kappa), kappa_sd = sd(kappa),
  label_switch_mean = mean(label_switch_proportion),
  label_switch_p025 = quantile(label_switch_proportion, 0.025, type = 6),
  label_switch_p975 = quantile(label_switch_proportion, 0.975, type = 6),
  c1_prevalence_mean = mean(c1_prevalence),
  c1_prevalence_p025 = quantile(c1_prevalence, 0.025, type = 6),
  c1_prevalence_p975 = quantile(c1_prevalence, 0.975, type = 6),
  jaccard_c1_mean = mean(jaccard_c1), jaccard_c2_mean = mean(jaccard_c2),
  centroid_profile_rmse_mean = mean(centroid_profile_rmse_all33),
  centroid_profile_correlation_mean = mean(centroid_profile_correlation_all33),
  mortality_risk_difference_mean = mean(mortality_risk_difference_c1_minus_c2),
  make30_risk_difference_mean = mean(make30_risk_difference_c1_minus_c2),
  interpretation = "continuous descriptive metrics; no post-hoc pass/fail threshold"
), by = .(variant_id, variant_label)]

write_csv_atomic(primary, file.path(DIRS$tables, "Table_D1_ablation_vs_reference.csv"))
write_csv_atomic(seed_results, file.path(DIRS$tables, "Table_D1_ablation_seed_stability.csv"))
write_csv_atomic(summary, file.path(DIRS$tables, "29b_Table_D1_ablation_seed_summary.csv"))

plot_data <- copy(seed_results)
plot_data[, variant_label := factor(variant_label, levels = rev(unique(variant_label)))]
p_ari <- ggplot(plot_data, aes(x = ari, y = variant_label)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey92", colour = "grey45", linewidth = 0.45) +
  geom_point(data = plot_data[locked_seed == TRUE], colour = "#2E2E2E", size = 2.0) +
  labs(
    x = "Adjusted Rand index versus locked K2 labels", y = NULL,
    title = "Feature-domain ablation of the locked empirical partition",
    subtitle = "Boxes summarize 20 prespecified k-means seeds; dark points show the locked seed"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"), axis.text = element_text(size = 9))
export_plot(p_ari, "Figure_D1_ablation_ARI")

p_prev <- ggplot(plot_data, aes(x = c1_prevalence, y = variant_label)) +
  geom_vline(xintercept = mean(reference == 1L), linetype = 2, colour = "grey70", linewidth = 0.5) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey92", colour = "grey45", linewidth = 0.45) +
  geom_point(data = plot_data[locked_seed == TRUE], colour = "#2E2E2E", size = 2.0) +
  scale_x_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(
    x = "Aligned C1 prevalence", y = NULL,
    title = "Cluster prevalence after feature-domain ablation",
    subtitle = "Dashed line marks the locked C1 prevalence (19.91%)"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"), axis.text = element_text(size = 9))
export_plot(p_prev, "Figure_D1_ablation_cluster_prevalence")

log_lines <- c(
  "# 29b feature-domain ablation log",
  "",
  paste0("- Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- Locked seed: ", LOCKED_SEED, "; additional prespecified seeds: ", length(SEEDS) - 1L, "."),
  paste0("- k-means: K=2, Lloyd, nstart=", NSTART, ", iter.max=", ITER_MAX, "."),
  "- A0 exactly reproduced the locked labels (ARI=1; switch proportion=0) before any ablation was run.",
  "- Labels were aligned to the locked reference by maximum patient agreement; mortality and MAKE30 were not used for alignment.",
  "- The non-renal SOFA proxy block contains the direct 33-feature counterparts of the five non-renal SOFA domains: PaO2/FiO2, platelets, bilirubin, MAP, and GCS.",
  "- Mortality and MAKE30 separations are descriptive only.",
  "- No post-hoc ARI or other pass/fail threshold was used.",
  "",
  "## Seed summary",
  paste(capture.output(print(summary)), collapse = "\n")
)
write_lines_atomic(log_lines, file.path(DIRS$logs, "29b_domain1_ablation_log.md"))
capture.output(sessionInfo(), file = file.path(DIRS$logs, "29b_sessionInfo.txt"))
write_lines_atomic("RUN COMPLETED", file.path(OUT_ROOT, "29b_run_completed.ok"))

cat("WP3 ablation completed.\n")
print(summary)

