## ============================================================
## 10b_supp_tableS1_figure2_mice_local.R
## Generate Table 1-related supplementary tables and figures
## using MICE-primary K=2 labels.
##
## Outputs are aggregate only. No patient-level frame is exported.
## ============================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

QC_DIR  <- file.path(DIR_OUTPUT, "qc")
FIG_DIR <- file.path(DIR_OUTPUT, "fig")
dir.create(QC_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PROJECT_ROOT <- dirname(DIR_DATA)
R_DIR <- if (exists("DIR_R")) DIR_R else file.path(PROJECT_ROOT, "r")

## -------------------------------
## 1. Required paths
## -------------------------------
path_final_full <- PATH_FINAL_FULL
path_eicu       <- file.path(DIR_DATA, "eicu_external.csv")
path_lab_mice   <- file.path(DIR_MODEL, "labels_primary_mice.rds")
path_X_mice     <- file.path(DIR_MODEL, "X_primary33_std_mice.rds")
path_sum_mice   <- file.path(QC_DIR, "mice_primary_label_summary.csv")
path_eicu_lab   <- file.path(DIR_OUTPUT, "eicu_cluster_labels_mice.csv")
path_feature_sets <- file.path(DIR_MODEL, "feature_sets.rds")
path_feature_labels <- file.path(R_DIR, "feature_labels.R")

required <- data.table(
  item = c(
    "PATH_FINAL_FULL",
    "eicu_external.csv",
    "labels_primary_mice.rds",
    "X_primary33_std_mice.rds",
    "mice_primary_label_summary.csv",
    "eicu_cluster_labels_mice.csv",
    "feature_sets.rds",
    "feature_labels.R"
  ),
  path = c(
    path_final_full,
    path_eicu,
    path_lab_mice,
    path_X_mice,
    path_sum_mice,
    path_eicu_lab,
    path_feature_sets,
    path_feature_labels
  )
)
required[, exists := file.exists(path)]
required[, size_MB := ifelse(exists, round(file.info(path)$size / 1024^2, 3), NA_real_)]
required[, modified_time := ifelse(exists, as.character(file.info(path)$mtime), NA_character_)]
print(required)
fwrite(required, file.path(QC_DIR, "supp_table1_related_required_files_audit.csv"))

must_exist <- required[item %in% c(
  "PATH_FINAL_FULL", "eicu_external.csv", "labels_primary_mice.rds",
  "X_primary33_std_mice.rds", "mice_primary_label_summary.csv",
  "eicu_cluster_labels_mice.csv"
)]
if (any(!must_exist$exists)) {
  stop("Missing required files:\n", paste(must_exist[exists == FALSE, path], collapse = "\n"))
}

## -------------------------------
## 2. Helper functions
## -------------------------------
num <- function(x) suppressWarnings(as.numeric(x))

med_iqr <- function(x) {
  x <- num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return("—")
  q <- quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
}

safe_quant <- function(x, p) {
  x <- num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

smd_fun <- function(x, g) {
  x <- num(x)
  g <- as.integer(g)
  x1 <- x[g == 1]
  x2 <- x[g == 2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s1 <- sd(x1)
  s2 <- sd(x2)
  sp <- sqrt(((length(x1) - 1) * s1^2 + (length(x2) - 1) * s2^2) / (length(x1) + length(x2) - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(x1) - mean(x2)) / sp
}

first_existing <- function(dt, candidates) {
  x <- intersect(candidates, names(dt))
  if (length(x) == 0) return(NA_character_)
  x[1]
}

clean_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\bmax\\b", "max", x, ignore.case = TRUE)
  x <- gsub("\\bmin\\b", "min", x, ignore.case = TRUE)
  x <- gsub("\\b24h\\b", "24h", x, ignore.case = TRUE)
  x
}

normalise_names <- function(nm) {
  tolower(gsub("[^a-zA-Z0-9]+", "_", nm))
}

## -------------------------------
## 3. Feature list and label dictionary
## -------------------------------
X_mice <- as.data.table(readRDS(path_X_mice))
lab_mice <- as.data.table(readRDS(path_lab_mice))

id_mimic <- intersect(c("stay_id", "hadm_id", "subject_id"), intersect(names(X_mice), names(lab_mice)))[1]
if (is.na(id_mimic)) stop("Cannot identify shared MIMIC ID column between X_primary33_std_mice.rds and labels_primary_mice.rds.")
if (!"cluster_k2" %in% names(lab_mice)) stop("labels_primary_mice.rds does not contain cluster_k2.")

if (file.exists(path_feature_sets)) {
  fs <- readRDS(path_feature_sets)
  if (is.list(fs) && "primary" %in% names(fs)) {
    FEATURES <- as.character(fs$primary)
    message("Using authoritative feature_sets.rds$primary")
  } else {
    FEATURES <- setdiff(names(X_mice), id_mimic)
    message("feature_sets.rds found but no $primary; using columns in X_primary33_std_mice.rds")
  }
} else {
  FEATURES <- setdiff(names(X_mice), id_mimic)
  message("feature_sets.rds not found; using columns in X_primary33_std_mice.rds")
}
FEATURES <- FEATURES[FEATURES %in% names(X_mice)]
if (length(FEATURES) != 33) warning("Expected 33 primary features, found ", length(FEATURES), ".")

load_feature_dict <- function(features) {
  base <- data.table(
    feature = features,
    feature_label = clean_label(features),
    unit = "",
    system = "Unclassified",
    display_order = seq_along(features)
  )
  if (!file.exists(path_feature_labels)) return(base)
  env <- new.env(parent = emptyenv())
  tryCatch(source(path_feature_labels, local = env), error = function(e) NULL)
  objs <- mget(ls(env), envir = env, ifnotfound = list(NULL))
  dfs <- Filter(function(z) is.data.frame(z) || is.data.table(z), objs)
  if (length(dfs) == 0) return(base)
  for (z in dfs) {
    z <- as.data.table(z)
    setnames(z, names(z), normalise_names(names(z)))
    feature_col <- first_existing(z, c("feature", "variable", "var", "name", "feature_name"))
    if (is.na(feature_col)) next
    if (!any(as.character(z[[feature_col]]) %in% features)) next
    label_col <- first_existing(z, c("feature_label", "label", "display_label", "display", "label_en", "pretty_name"))
    unit_col  <- first_existing(z, c("unit", "units"))
    sys_col   <- first_existing(z, c("system", "domain", "organ_system", "category", "feature_group", "group"))
    out <- data.table(feature = as.character(z[[feature_col]]))
    out[, feature_label := if (!is.na(label_col)) as.character(z[[label_col]]) else clean_label(feature)]
    out[, unit := if (!is.na(unit_col)) as.character(z[[unit_col]]) else ""]
    out[, system := if (!is.na(sys_col)) as.character(z[[sys_col]]) else "Unclassified"]
    out <- out[feature %in% features]
    out <- unique(out, by = "feature")
    ans <- merge(base[, .(feature, display_order)], out, by = "feature", all.x = TRUE, sort = FALSE)
    ans[is.na(feature_label), feature_label := clean_label(feature)]
    ans[is.na(unit), unit := ""]
    ans[is.na(system) | system == "", system := "Unclassified"]
    setorder(ans, display_order)
    return(ans)
  }
  base
}

feat_dict <- load_feature_dict(FEATURES)
fwrite(feat_dict, file.path(QC_DIR, "supp_feature_dictionary_used.csv"))

## -------------------------------
## 4. Load and label MIMIC/eICU frames
## -------------------------------
Mraw <- fread(path_final_full)
Eraw <- fread(path_eicu)
Elab <- fread(path_eicu_lab)

if ("cluster_k2" %in% names(Mraw)) Mraw[, cluster_k2 := NULL]
if ("cluster_k2" %in% names(Eraw)) Eraw[, cluster_k2 := NULL]

if (!id_mimic %in% names(Mraw)) stop("PATH_FINAL_FULL does not contain MIMIC ID column: ", id_mimic)
if (anyDuplicated(lab_mice[[id_mimic]]) > 0) stop("Duplicate IDs in labels_primary_mice.rds.")

M <- merge(
  Mraw,
  lab_mice[, .SD, .SDcols = c(id_mimic, "cluster_k2")],
  by = id_mimic,
  all.y = TRUE
)
M[, cohort := "MIMIC-IV"]
M[, cohort_key := "MIMIC"]
M[, phenotype_short := fifelse(cluster_k2 == 1, "C1", "C2")]
M[, phenotype := fifelse(cluster_k2 == 1, "C1 higher-risk", "C2 lower-risk")]
M[, group_key := paste0(cohort_key, "_", phenotype_short)]

id_eicu <- intersect(c("patientunitstayid", "stay_id", "hadm_id", "subject_id", "uniquepid"), intersect(names(Eraw), names(Elab)))[1]
if (is.na(id_eicu)) stop("Cannot identify shared eICU ID column between eicu_external.csv and eicu_cluster_labels_mice.csv.")
if (!"cluster_k2" %in% names(Elab)) stop("eicu_cluster_labels_mice.csv does not contain cluster_k2.")
if (anyDuplicated(Elab[[id_eicu]]) > 0) stop("Duplicate IDs in eicu_cluster_labels_mice.csv.")

E <- merge(
  Eraw,
  Elab[, .SD, .SDcols = c(id_eicu, "cluster_k2")],
  by = id_eicu,
  all.y = TRUE
)
E[, cohort := "eICU"]
E[, cohort_key := "eICU"]
E[, phenotype_short := fifelse(cluster_k2 == 1, "C1", "C2")]
E[, phenotype := fifelse(cluster_k2 == 1, "C1 higher-risk", "C2 lower-risk")]
E[, group_key := paste0(cohort_key, "_", phenotype_short)]

label_dist <- rbindlist(list(
  M[, .N, by = .(cohort, cohort_key, cluster_k2, phenotype_short, phenotype)],
  E[, .N, by = .(cohort, cohort_key, cluster_k2, phenotype_short, phenotype)]
), fill = TRUE)
label_dist[, pct := round(100 * N / sum(N), 2), by = cohort]
setorder(label_dist, cohort_key, cluster_k2)
fwrite(label_dist, file.path(QC_DIR, "supp_label_distribution_mice.csv"))

prev_drift <- dcast(label_dist[cluster_k2 == 1], . ~ cohort_key, value.var = "pct")
prev_drift[, `:=`(
  comparison = "High-risk phenotype prevalence: eICU minus MIMIC",
  MIMIC_pct_C1 = MIMIC,
  eICU_pct_C1 = eICU,
  absolute_difference_percentage_points = round(eICU - MIMIC, 2),
  relative_ratio = round(eICU / MIMIC, 3)
)]
prev_drift <- prev_drift[, .(comparison, MIMIC_pct_C1, eICU_pct_C1, absolute_difference_percentage_points, relative_ratio)]
fwrite(prev_drift, file.path(QC_DIR, "supp_prevalence_drift_mice.csv"))

## -------------------------------
## 5. Supplementary Table S1: raw 33-feature summaries
## -------------------------------
summarise_raw_feature <- function(D, feature, cohort_key) {
  if (!feature %in% names(D)) {
    return(data.table(
      feature = feature,
      cohort_key = cohort_key,
      group_key = paste0(cohort_key, "_", c("C1", "C2")),
      summary = "—",
      n_nonmissing = 0L,
      missing_pct = NA_real_
    ))
  }
  tmp <- data.table(
    group_key = D$group_key,
    value = num(D[[feature]])
  )
  tmp[, .(
    summary = med_iqr(value),
    n_nonmissing = sum(is.finite(value)),
    missing_pct = round(100 * mean(!is.finite(value)), 2)
  ), by = group_key][, `:=`(feature = feature, cohort_key = cohort_key)][]
}

raw_long <- rbindlist(c(
  lapply(FEATURES, function(f) summarise_raw_feature(M, f, "MIMIC")),
  lapply(FEATURES, function(f) summarise_raw_feature(E, f, "eICU"))
), fill = TRUE)

raw_wide <- dcast(raw_long, feature ~ group_key, value.var = "summary")
for (cn in c("MIMIC_C1", "MIMIC_C2", "eICU_C1", "eICU_C2")) {
  if (!cn %in% names(raw_wide)) raw_wide[, (cn) := "—"]
}

missing_wide <- dcast(raw_long, feature ~ group_key, value.var = "missing_pct")
setnames(missing_wide, setdiff(names(missing_wide), "feature"), paste0(setdiff(names(missing_wide), "feature"), "_missing_pct"))

## -------------------------------
## 6. Standardized profile and SMDs
## -------------------------------
## MIMIC z-profile uses the MICE-primary standardized matrix directly.
ZM <- merge(
  X_mice[, .SD, .SDcols = c(id_mimic, FEATURES)],
  lab_mice[, .SD, .SDcols = c(id_mimic, "cluster_k2")],
  by = id_mimic,
  all.x = TRUE
)
ZM[, cohort := "MIMIC-IV"]
ZM[, cohort_key := "MIMIC"]
ZM[, phenotype_short := fifelse(cluster_k2 == 1, "C1", "C2")]
ZM[, phenotype := fifelse(cluster_k2 == 1, "C1 higher-risk", "C2 lower-risk")]

## eICU z-profile: standardize eICU raw features using MIMIC raw 1%-99% winsorization, median imputation, and MIMIC mean/sd.
standardize_eicu_to_mimic <- function(feature) {
  if (!feature %in% names(M) || !feature %in% names(E)) return(NULL)
  m <- num(M[[feature]])
  e <- num(E[[feature]])
  coverage_e <- mean(is.finite(e))
  if (!is.finite(coverage_e) || coverage_e < 0.30) return(NULL)
  qlo <- safe_quant(m, 0.01)
  qhi <- safe_quant(m, 0.99)
  med <- safe_quant(m, 0.50)
  if (!is.finite(qlo) || !is.finite(qhi) || !is.finite(med) || qlo >= qhi) return(NULL)
  mw <- pmin(pmax(m, qlo), qhi)
  mw[!is.finite(mw)] <- med
  mu <- mean(mw, na.rm = TRUE)
  sig <- sd(mw, na.rm = TRUE)
  if (!is.finite(sig) || sig == 0) return(NULL)
  ew <- pmin(pmax(e, qlo), qhi)
  ew[!is.finite(ew)] <- med
  z <- (ew - mu) / sig
  list(feature = feature, z = z, eicu_coverage = coverage_e, qlo = qlo, qhi = qhi, median = med, mean_ref = mu, sd_ref = sig)
}

std_list <- lapply(FEATURES, standardize_eicu_to_mimic)
std_list <- std_list[!vapply(std_list, is.null, logical(1))]
E_COMMON_FEATURES <- vapply(std_list, function(z) z$feature, character(1))

ZE <- data.table(id_eicu_value = E[[id_eicu]])
for (obj in std_list) ZE[, (obj$feature) := obj$z]
ZE[, cluster_k2 := E$cluster_k2]
ZE[, cohort := "eICU"]
ZE[, cohort_key := "eICU"]
ZE[, phenotype_short := fifelse(cluster_k2 == 1, "C1", "C2")]
ZE[, phenotype := fifelse(cluster_k2 == 1, "C1 higher-risk", "C2 lower-risk")]

std_audit <- rbindlist(lapply(std_list, function(obj) {
  data.table(
    feature = obj$feature,
    eicu_coverage = round(obj$eicu_coverage, 3),
    mimic_winsor_p01 = obj$qlo,
    mimic_winsor_p99 = obj$qhi,
    mimic_imputation_median = obj$median,
    mimic_reference_mean = obj$mean_ref,
    mimic_reference_sd = obj$sd_ref
  )
}), fill = TRUE)

std_audit <- merge(feat_dict, std_audit, by = "feature", all.x = TRUE, sort = FALSE)
std_audit[, used_for_eicu_profile := !is.na(eicu_coverage)]
fwrite(std_audit, file.path(QC_DIR, "supp_eicu_mimic_coordinate_standardization_audit.csv"))

profile_one <- function(D, feature) {
  if (!feature %in% names(D)) return(NULL)
  D[, .(
    mean_z = mean(num(get(feature)), na.rm = TRUE),
    median_z = median(num(get(feature)), na.rm = TRUE),
    q1_z = safe_quant(get(feature), 0.25),
    q3_z = safe_quant(get(feature), 0.75),
    n_nonmissing = sum(is.finite(num(get(feature))))
  ), by = .(cohort, cohort_key, cluster_k2, phenotype_short, phenotype)][, feature := feature][]
}

profile_long <- rbindlist(c(
  lapply(FEATURES, function(f) profile_one(ZM, f)),
  lapply(E_COMMON_FEATURES, function(f) profile_one(ZE, f))
), fill = TRUE)
profile_long <- merge(feat_dict, profile_long, by = "feature", all.y = TRUE, sort = FALSE)
setorder(profile_long, cohort_key, display_order, cluster_k2)
fwrite(profile_long, file.path(QC_DIR, "figure2_z_profile_mice_labels_data.csv"))

profile_diff <- dcast(
  profile_long,
  cohort + cohort_key + feature + feature_label + unit + system + display_order ~ phenotype_short,
  value.var = "mean_z"
)
for (cn in c("C1", "C2")) if (!cn %in% names(profile_diff)) profile_diff[, (cn) := NA_real_]
profile_diff[, difference_C1_minus_C2 := C1 - C2]
setorder(profile_diff, cohort_key, display_order)
fwrite(profile_diff, file.path(QC_DIR, "figure2_z_profile_mice_labels_difference.csv"))

smd_one <- function(D, feature, cohort_key) {
  if (!feature %in% names(D)) return(data.table(feature = feature, cohort_key = cohort_key, smd_C1_minus_C2 = NA_real_))
  data.table(
    feature = feature,
    cohort_key = cohort_key,
    smd_C1_minus_C2 = smd_fun(D[[feature]], D$cluster_k2)
  )
}

smd_long <- rbindlist(c(
  lapply(FEATURES, function(f) smd_one(ZM, f, "MIMIC")),
  lapply(FEATURES, function(f) smd_one(ZE, f, "eICU"))
), fill = TRUE)
smd_wide <- dcast(smd_long, feature ~ cohort_key, value.var = "smd_C1_minus_C2")
setnames(smd_wide, intersect(c("MIMIC", "eICU"), names(smd_wide)), paste0(intersect(c("MIMIC", "eICU"), names(smd_wide)), "_SMD_C1_minus_C2"))

## Supplementary Table S1 final
TableS1 <- Reduce(function(x, y) merge(x, y, by = "feature", all.x = TRUE, sort = FALSE), list(
  feat_dict,
  raw_wide,
  smd_wide,
  missing_wide
))
setorder(TableS1, display_order)

## More readable column order
preferred_cols <- c(
  "system", "feature", "feature_label", "unit",
  "MIMIC_C1", "MIMIC_C2", "MIMIC_SMD_C1_minus_C2",
  "eICU_C1", "eICU_C2", "eICU_SMD_C1_minus_C2",
  "MIMIC_C1_missing_pct", "MIMIC_C2_missing_pct",
  "eICU_C1_missing_pct", "eICU_C2_missing_pct",
  "display_order"
)
setcolorder(TableS1, intersect(preferred_cols, names(TableS1)))
fwrite(TableS1, file.path(QC_DIR, "tableS1_33_features_by_mice_subphenotype.csv"))

## Missingness audit in long and wide formats
missing_long <- raw_long[, .(feature, cohort_key, group_key, n_nonmissing, missing_pct)]
missing_long <- merge(feat_dict, missing_long, by = "feature", all.y = TRUE, sort = FALSE)
setorder(missing_long, cohort_key, display_order, group_key)
fwrite(missing_long, file.path(QC_DIR, "tableS1_33_features_missingness_audit.csv"))

## -------------------------------
## 7. Figures
## -------------------------------
## Figure 2: standardized C1-C2 profile.
fig2_dat <- profile_diff[!is.na(difference_C1_minus_C2)]
fig2_dat[, feature_plot := factor(feature_label, levels = rev(unique(feat_dict$feature_label[order(feat_dict$display_order)])))]

p_fig2 <- ggplot(fig2_dat, aes(x = feature_plot, y = difference_C1_minus_C2)) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_segment(aes(xend = feature_plot, y = 0, yend = difference_C1_minus_C2), linewidth = 0.35) +
  geom_point(size = 1.7) +
  coord_flip() +
  facet_wrap(~ cohort, nrow = 1) +
  labs(
    x = NULL,
    y = "Mean standardized difference (C1 − C2)",
    title = "Physiological profile difference by MICE-primary K=2 subphenotype",
    subtitle = "C1 denotes the higher-risk phenotype; eICU features are standardized to the MIMIC raw-feature coordinate system."
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    axis.text.y = element_text(size = 6.7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "figure2_z_profile_mice_labels.png"), p_fig2, width = 9.5, height = 8.2, dpi = 320)
ggsave(file.path(FIG_DIR, "figure2_z_profile_mice_labels.pdf"), p_fig2, width = 9.5, height = 8.2, device = cairo_pdf)

## Figure S: phenotype prevalence distribution.
dist_plot <- copy(label_dist)
dist_plot[, phenotype := factor(phenotype, levels = c("C1 higher-risk", "C2 lower-risk"))]

p_dist <- ggplot(dist_plot, aes(x = cohort, y = pct, fill = phenotype)) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_text(aes(label = paste0(N, "\n", pct, "%")), position = position_stack(vjust = 0.5), size = 3.1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(
    x = NULL,
    y = "Patients (%)",
    fill = NULL,
    title = "Distribution of MICE-primary K=2 subphenotypes",
    subtitle = "High-risk phenotype prevalence is higher in eICU than MIMIC-IV."
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    legend.position = "top"
  )

ggsave(file.path(FIG_DIR, "figureS_label_distribution_mice.png"), p_dist, width = 5.8, height = 4.2, dpi = 320)
ggsave(file.path(FIG_DIR, "figureS_label_distribution_mice.pdf"), p_dist, width = 5.8, height = 4.2, device = cairo_pdf)

## -------------------------------
## 8. Optional DOCX export for Table S1
## -------------------------------
if (requireNamespace("officer", quietly = TRUE) && requireNamespace("flextable", quietly = TRUE)) {
  tryCatch({
    TableS1_doc <- copy(TableS1)
    ## Trim to main columns for Word readability; full CSV contains missingness columns.
    doc_cols <- intersect(c(
      "system", "feature_label", "unit",
      "MIMIC_C1", "MIMIC_C2", "MIMIC_SMD_C1_minus_C2",
      "eICU_C1", "eICU_C2", "eICU_SMD_C1_minus_C2"
    ), names(TableS1_doc))
    TableS1_doc <- TableS1_doc[, ..doc_cols]
    ft <- flextable::flextable(TableS1_doc)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft, size = 8, part = "all")
    ft <- flextable::autofit(ft)
    doc <- officer::read_docx()
    doc <- officer::body_add_par(
      doc,
      "Supplementary Table S1. First-24h physiological features by MICE-primary K=2 subphenotype.",
      style = "heading 1"
    )
    doc <- flextable::body_add_flextable(doc, value = ft)
    doc <- officer::body_add_par(
      doc,
      "Values are median [IQR]. SMD is calculated as C1 minus C2 using standardized feature values. C1 denotes the higher-risk phenotype.",
      style = "Normal"
    )
    print(doc, target = file.path(QC_DIR, "tableS1_33_features_by_mice_subphenotype.docx"))
  }, error = function(e) {
    message("DOCX export failed; CSV outputs remain available. Error: ", conditionMessage(e))
  })
} else {
  message("Packages officer/flextable not installed; DOCX export skipped. CSV and figure outputs were generated.")
}

## -------------------------------
## 9. Output audit
## -------------------------------
outputs <- data.table(
  output = c(
    "tableS1_33_features_by_mice_subphenotype.csv",
    "tableS1_33_features_by_mice_subphenotype.docx",
    "tableS1_33_features_missingness_audit.csv",
    "figure2_z_profile_mice_labels_data.csv",
    "figure2_z_profile_mice_labels_difference.csv",
    "figure2_z_profile_mice_labels.png",
    "figure2_z_profile_mice_labels.pdf",
    "supp_label_distribution_mice.csv",
    "supp_prevalence_drift_mice.csv",
    "figureS_label_distribution_mice.png",
    "figureS_label_distribution_mice.pdf",
    "supp_eicu_mimic_coordinate_standardization_audit.csv",
    "supp_feature_dictionary_used.csv"
  ),
  path = c(
    file.path(QC_DIR, "tableS1_33_features_by_mice_subphenotype.csv"),
    file.path(QC_DIR, "tableS1_33_features_by_mice_subphenotype.docx"),
    file.path(QC_DIR, "tableS1_33_features_missingness_audit.csv"),
    file.path(QC_DIR, "figure2_z_profile_mice_labels_data.csv"),
    file.path(QC_DIR, "figure2_z_profile_mice_labels_difference.csv"),
    file.path(FIG_DIR, "figure2_z_profile_mice_labels.png"),
    file.path(FIG_DIR, "figure2_z_profile_mice_labels.pdf"),
    file.path(QC_DIR, "supp_label_distribution_mice.csv"),
    file.path(QC_DIR, "supp_prevalence_drift_mice.csv"),
    file.path(FIG_DIR, "figureS_label_distribution_mice.png"),
    file.path(FIG_DIR, "figureS_label_distribution_mice.pdf"),
    file.path(QC_DIR, "supp_eicu_mimic_coordinate_standardization_audit.csv"),
    file.path(QC_DIR, "supp_feature_dictionary_used.csv")
  )
)
outputs[, exists := file.exists(path)]
outputs[, size_MB := ifelse(exists, round(file.info(path)$size / 1024^2, 3), NA_real_)]
fwrite(outputs, file.path(QC_DIR, "supp_table1_related_outputs_audit.csv"))
print(outputs)

cat("\nFinished supplementary Table 1-related outputs.\n")
cat("Key files:\n")
cat(" - ", file.path(QC_DIR, "tableS1_33_features_by_mice_subphenotype.csv"), "\n", sep = "")
cat(" - ", file.path(FIG_DIR, "figure2_z_profile_mice_labels.png"), "\n", sep = "")
cat(" - ", file.path(FIG_DIR, "figureS_label_distribution_mice.png"), "\n", sep = "")
cat(" - ", file.path(QC_DIR, "supp_prevalence_drift_mice.csv"), "\n", sep = "")
