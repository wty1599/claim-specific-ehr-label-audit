## =====================================================================
## 29c_domain1_loading_formalization.R
## Formalize the locked 33-feature Fisher direction and descriptive
## sum(w^2) contributions. Contributions are not causal or independent.
## =====================================================================

suppressPackageStartupMessages({
  required <- c("data.table", "ggplot2", "digest")
  unavailable <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(unavailable)) stop("Missing R packages: ", paste(unavailable, collapse = ", "))
  library(data.table)
  library(ggplot2)
})

PROJECT_ROOT <- normalizePath(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), winslash = "/", mustWork = TRUE)
OUT_ROOT <- normalizePath(file.path(PROJECT_ROOT, "domain1_renal_identity_and_method_closure_20260717"), winslash = "/", mustWork = TRUE)
DIRS <- list(
  tables = file.path(OUT_ROOT, "tables"), figures = file.path(OUT_ROOT, "figures"),
  logs = file.path(OUT_ROOT, "logs"), provenance = file.path(OUT_ROOT, "provenance")
)
INPUT <- list(
  loadings = file.path(PROJECT_ROOT, "analysis_archive", "simulations", "Domain1_copula_ladder_step0_step1_20260717", "output", "formal", "step0_fisher_loadings.csv"),
  matrix = file.path(PROJECT_ROOT, "output", "model", "X_primary33_std_mice.rds"),
  labels = file.path(PROJECT_ROOT, "output", "model", "labels_primary_mice.rds"),
  registry = file.path(DIRS$provenance, "D1_feature_domain_registry.csv"),
  labels_dictionary = file.path(PROJECT_ROOT, "r", "feature_labels.R"),
  wp1 = file.path(DIRS$tables, "Table_D1_identity_models_oof.csv"),
  wp3 = file.path(DIRS$tables, "29b_Table_D1_ablation_seed_summary.csv")
)
if (any(!file.exists(unlist(INPUT)))) stop("Missing input(s): ", paste(names(INPUT)[!file.exists(unlist(INPUT))], collapse = ", "))

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
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

export_plot <- function(plot, stem, width = 8.0, height = 6.8) {
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

loading <- fread(INPUT$loadings)
registry <- fread(INPUT$registry)
X_dt <- as.data.table(readRDS(INPUT$matrix))
labels <- as.data.table(readRDS(INPUT$labels))
wp1 <- fread(INPUT$wp1)
wp3 <- fread(INPUT$wp3)
if (nrow(loading) != 33L || anyDuplicated(loading$feature)) stop("Loading file is not a unique 33-feature table")
if (!setequal(loading$feature, registry$feature)) stop("Loading/registry feature mismatch")
if (!identical(X_dt$stay_id, labels$stay_id)) stop("Locked matrix/label mismatch")
if (!all(loading$feature %in% names(X_dt))) stop("Loading feature absent from locked matrix")

label_env <- new.env(parent = baseenv())
source(INPUT$labels_dictionary, local = label_env, encoding = "UTF-8")
display_labels <- label_env$FEAT_LABEL

loading <- merge(loading, registry, by = "feature", all.x = TRUE, sort = FALSE)
setorder(loading, loading_rank)
loading[, `:=`(
  standardized_variable_name = feature,
  display_label = unname(display_labels[feature]),
  signed_w = fisher_loading_unit_norm,
  abs_w = abs(fisher_loading_unit_norm),
  w_squared = fisher_loading_unit_norm^2
)]
if (anyNA(loading$display_label)) loading[is.na(display_label), display_label := feature]
if (abs(sum(loading$w_squared) - 1) > 1e-8) stop("Unit-normalized loading squares do not sum to 1")

selected_features <- loading$feature
X <- as.matrix(X_dt[, ..selected_features])
storage.mode(X) <- "double"
w <- as.numeric(loading$signed_w)
names(w) <- loading$feature
score <- as.numeric(X %*% w)
mean_score_c1 <- mean(score[labels$cluster_k2 == 1L])
mean_score_c2 <- mean(score[labels$cluster_k2 == 2L])
direction <- if (mean_score_c1 > mean_score_c2) "positive signed w points toward C1" else "positive signed w points toward C2"

loading[, partition_group := fifelse(
  strict_renal_core == 1L, "Renal core",
  fifelse(strict_acid_base == 1L, "Acid-base",
    fifelse(nonrenal_sofa_proxy == 1L, "Non-renal SOFA direct proxy", "Remaining")
  )
)]

partition <- loading[, .(
  n_features = .N,
  sum_w_squared = sum(w_squared),
  share_percent = 100 * sum(w_squared),
  features = paste(feature, collapse = ";"),
  row_type = "mutually exclusive primary partition"
), by = .(system = partition_group)]

make_overlap_row <- function(name, features, status) {
  subset <- loading[feature %in% features]
  data.table(
    system = name,
    n_features = nrow(subset),
    sum_w_squared = sum(subset$w_squared),
    share_percent = 100 * sum(subset$w_squared),
    features = paste(subset$feature, collapse = ";"),
    row_type = status
  )
}

renal <- loading[strict_renal_core == 1L, feature]
acid <- loading[strict_acid_base == 1L, feature]
expanded <- loading[expanded_renal_acid_metabolic == 1L, feature]
nonrenal <- loading[nonrenal_sofa_proxy == 1L, feature]
overlap <- rbindlist(list(
  make_overlap_row("BUN + creatinine", c("bun_max", "creatinine_max"), "overlapping prespecified summary"),
  make_overlap_row("Strict renal core", renal, "overlapping prespecified summary"),
  make_overlap_row("Strict acid-base", acid, "overlapping prespecified summary"),
  make_overlap_row("Strict renal + acid-base", union(renal, acid), "overlapping prespecified summary"),
  make_overlap_row("Expanded renal/acid-base/metabolic", expanded, "sensitivity definition; includes lactate and potassium maximum"),
  make_overlap_row("Non-renal SOFA direct proxy", nonrenal, "overlapping prespecified summary")
))
contribution <- rbindlist(list(partition, overlap), fill = TRUE)
contribution[, interpretation := "sum(w^2) is a descriptive share of the standardized discriminant direction, not an independent or causal contribution"]

loading_output <- loading[, .(
  standardized_variable_name, display_label, clinical_domain = primary_domain,
  signed_w, abs_w, loading_rank, w_squared, partition_group,
  strict_renal_core, strict_acid_base, expanded_renal_acid_metabolic,
  nonrenal_sofa_proxy, definition_status,
  direction_check = direction,
  source_loading_file = INPUT$loadings
)]

wp1_values <- setNames(wp1$auc, wp1$model_id)
wp3_values <- setNames(wp3$ari_mean, wp3$variant_id)
cross_check <- data.table(
  evidence_id = c("WP4_loading", "WP1_prediction", "WP3_ablation", "WP3_nonrenal_control"),
  observed_quantity = c(
    contribution[system == "Strict renal + acid-base", share_percent],
    wp1_values[["renal_acid_base"]] - wp1_values[["sofa_nonrenal"]],
    wp3_values[["A4"]],
    wp3_values[["A5"]]
  ),
  quantity_definition = c(
    "strict renal+acid-base sum(w^2) share, percent",
    "OOF AUC renal+acid-base minus non-renal SOFA",
    "mean ARI after removing strict renal+acid-base",
    "mean ARI after removing non-renal SOFA proxy block"
  ),
  direction = c(
    "large loading share supports concentration in renal/acid-base variables",
    "positive difference supports stronger label reconstruction by renal/acid-base variables",
    "lower ARI indicates structural dependence on renal/acid-base variables",
    "high ARI indicates limited dependence on the non-renal SOFA proxy block"
  ),
  pass_fail_use = "none; cross-work-package consistency is descriptive"
)

write_csv_atomic(loading_output, file.path(DIRS$tables, "Table_D1_discriminant_loadings_all33.csv"))
write_csv_atomic(contribution, file.path(DIRS$tables, "Table_D1_loading_contribution_by_system.csv"))
write_csv_atomic(cross_check, file.path(DIRS$tables, "29c_Table_D1_cross_workpackage_consistency.csv"))

plot_loading <- copy(loading_output)
plot_loading[, display_label := factor(display_label, levels = rev(display_label[order(abs_w)]))]
p_load <- ggplot(plot_loading, aes(x = signed_w, y = display_label)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = signed_w, yend = display_label), colour = "grey55", linewidth = 0.55) +
  geom_point(size = 2.0, colour = "#2E2E2E") +
  labs(
    x = "Unit-normalized Fisher loading (signed toward C1)", y = NULL,
    title = "Locked discriminant direction across 33 standardized features",
    subtitle = "Loadings are descriptive; their squared values are not causal contributions"
  ) +
  theme_classic(base_family = "sans", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"), axis.text = element_text(size = 8.2))
export_plot(p_load, "Figure_D1_loading_lollipop", height = 8.2)

plot_system <- partition[order(share_percent)]
plot_system[, system := factor(system, levels = system)]
p_system <- ggplot(plot_system, aes(x = share_percent, y = system)) +
  geom_col(width = 0.62, fill = "#3A5A8C") +
  geom_text(aes(label = sprintf("%.1f%%", share_percent)), hjust = -0.12, size = 3.1, family = "sans") +
  scale_x_continuous(limits = c(0, max(plot_system$share_percent) * 1.18), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Share of sum(w^2), %", y = NULL,
    title = "Descriptive contribution of prespecified feature systems",
    subtitle = "Primary groups are mutually exclusive and sum to 100%"
  ) +
  theme_classic(base_family = "sans", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9, colour = "grey35"), axis.text = element_text(size = 9))
export_plot(p_system, "Figure_D1_loading_system_contribution", height = 4.6)

provenance <- data.table(
  source = names(INPUT), path = unlist(INPUT),
  sha256 = vapply(unlist(INPUT), sha256, character(1)), read_only = TRUE
)
write_csv_atomic(provenance, file.path(DIRS$provenance, "29c_loading_input_checksums.csv"))

strict_share <- contribution[system == "Strict renal + acid-base", share_percent]
expanded_share <- contribution[system == "Expanded renal/acid-base/metabolic", share_percent]
nonrenal_share <- contribution[
  system == "Non-renal SOFA direct proxy" & row_type == "overlapping prespecified summary",
  share_percent
]
log_lines <- c(
  "# 29c loading provenance and formalization",
  "",
  paste0("- Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- Formal loading source: `", INPUT$loadings, "`."),
  paste0("- Loading SHA-256: `", sha256(INPUT$loadings), "`."),
  paste0("- Locked matrix SHA-256: `", sha256(INPUT$matrix), "`."),
  paste0("- Locked labels SHA-256: `", sha256(INPUT$labels), "`."),
  paste0("- Score means: C1=", sprintf("%.6f", mean_score_c1), "; C2=", sprintf("%.6f", mean_score_c2), "; ", direction, "."),
  paste0("- Strict renal+acid-base share: ", sprintf("%.3f%%", strict_share), "."),
  paste0("- Expanded renal/acid-base/metabolic share: ", sprintf("%.3f%%", expanded_share), "."),
  paste0("- Non-renal SOFA direct-proxy share: ", sprintf("%.3f%%", nonrenal_share), "."),
  "- The expanded definition adds lactate_max and potassium_max and is reported only as a sensitivity definition.",
  "- sum(w^2) is descriptive and is not an independent or causal contribution.",
  "",
  "## Cross-work-package consistency",
  paste(capture.output(print(cross_check)), collapse = "\n")
)
write_lines_atomic(log_lines, file.path(DIRS$logs, "29c_loading_provenance.md"))
capture.output(sessionInfo(), file = file.path(DIRS$logs, "29c_sessionInfo.txt"))
write_lines_atomic("RUN COMPLETED", file.path(OUT_ROOT, "29c_run_completed.ok"))

cat("WP4 loading formalization completed.\n")
print(contribution)
print(cross_check)
