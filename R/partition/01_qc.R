# ============================================================




# ============================================================

source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

library(data.table)

# ---------------------------------------------------------

# ---------------------------------------------------------

if (!file.exists(PATH_FINAL_FULL)) {
  cat("final_full.csv not found; CSV files in the data directory:\n")
  print(list.files(DIR_DATA, pattern = "\\.csv$", ignore.case = TRUE))
  stop("Check PATH_FINAL_FULL.")
}

dt <- fread(PATH_FINAL_FULL)
N  <- nrow(dt)
cat(sprintf("Loaded: %d rows x %d columns (expected 20049 x 86)\n", N, ncol(dt)))


clust_vars <- c(
  "wbc_max","wbc_min","platelets_min","hemoglobin_min",
  "abs_neutrophils_max","abs_lymphocytes_min","bands_max",
  "bilirubin_total_max","alt_max","ast_max","albumin_min",
  "inr_max","pt_max","ptt_max","fibrinogen_min",
  "creatinine_max","creatinine_min","bun_max","aniongap_max",
  "bicarbonate_min","glucose_max","sodium_min","sodium_max",
  "potassium_min","potassium_max","calcium_min","calcium_max",
  "chloride_min","chloride_max","lactate_max","ph_min",
  "po2_min","pco2_max","so2_min","pao2fio2ratio_min","baseexcess_min",
  "heart_rate_max","sbp_min","mbp_min","mbp_mean","resp_rate_max",
  "temperature_min","temperature_max","spo2_min","gcs_min",
  "urine_output_24h_ml"
)

miss_in_data <- setdiff(clust_vars, names(dt))
if (length(miss_in_data)) {
  cat("Warning: clustering features not found:\n"); print(miss_in_data)
}
clust_vars <- intersect(clust_vars, names(dt))
cat(sprintf("Clustering features included in QC: %d\n", length(clust_vars)))


# ---------------------------------------------------------

# ---------------------------------------------------------

miss_tbl <- data.table(
  variable  = names(dt),
  n_missing = vapply(dt, function(x) sum(is.na(x)), numeric(1))
)
miss_tbl[, miss_rate_pct := round(n_missing / N * 100, 2)]
miss_tbl[, is_cluster_var := variable %in% clust_vars]
miss_tbl[, tier := fifelse(miss_rate_pct < 10,  "1_usable_<10%",
                    fifelse(miss_rate_pct <= 40, "2_MICE_10-40%",
                                                 "3_drop_>40%"))]
setorder(miss_tbl, -miss_rate_pct)
miss_clust <- miss_tbl[is_cluster_var == TRUE]

cat("\n---- Clustering-feature missingness strata ----\n"); print(miss_clust[, .N, by = tier])
cat("\n---- >40% (consider exclusion) ----\n"); print(miss_clust[tier == "3_drop_>40%", .(variable, miss_rate_pct)])
cat("\n---- 10-40% (MICE candidates) ----\n"); print(miss_clust[tier == "2_MICE_10-40%", .(variable, miss_rate_pct)])

fwrite(miss_tbl, file.path(DIR_QC, "qc_missing_report.csv"))


# ---------------------------------------------------------

# ---------------------------------------------------------

num_clust <- clust_vars[vapply(dt[, ..clust_vars], is.numeric, logical(1))]

dist_tbl <- rbindlist(lapply(num_clust, function(v) {
  xo <- dt[[v]][!is.na(dt[[v]])]
  q  <- quantile(xo, c(0.01, 0.25, 0.5, 0.75, 0.99), names = FALSE)
  data.table(variable = v, n_valid = length(xo),
             mean = round(mean(xo),3), sd = round(sd(xo),3),
             min = round(min(xo),3), p1 = round(q[1],3), p25 = round(q[2],3),
             median = round(q[3],3), p75 = round(q[4],3), p99 = round(q[5],3),
             max = round(max(xo),3))
}))
cat("\n---- Distribution summary ----\n"); print(dist_tbl)
fwrite(dist_tbl, file.path(DIR_QC, "qc_distribution_report.csv"))


# ---------------------------------------------------------

# ---------------------------------------------------------

phys_ranges <- list(
  wbc_max = c(0.1,100),  wbc_min = c(0.1,100),
  platelets_min = c(5,1500),  hemoglobin_min = c(2,25),
  abs_neutrophils_max = c(0,100),  abs_lymphocytes_min = c(0,50),
  bands_max = c(0,100),  bilirubin_total_max = c(0,60),
  alt_max = c(0,15000),  ast_max = c(0,15000),  albumin_min = c(0.5,7),
  inr_max = c(0.5,25),  pt_max = c(5,150),  ptt_max = c(10,250),
  fibrinogen_min = c(20,1200),
  creatinine_max = c(0.1,30),  creatinine_min = c(0.1,30),  bun_max = c(1,300),
  aniongap_max = c(0,50),  bicarbonate_min = c(3,60),  glucose_max = c(10,2000),
  sodium_min = c(100,190),  sodium_max = c(100,190),
  potassium_min = c(1,12),  potassium_max = c(1,12),
  calcium_min = c(3,20),  calcium_max = c(3,20),
  chloride_min = c(60,160), chloride_max = c(60,160),
  lactate_max = c(0.1,40),  ph_min = c(6.5,8.0),
  po2_min = c(10,700),  pco2_max = c(5,200),  so2_min = c(20,100),
  pao2fio2ratio_min = c(20,700),  baseexcess_min = c(-40,40),
  heart_rate_max = c(10,350),  sbp_min = c(20,300),
  mbp_min = c(10,250),  mbp_mean = c(10,250),  resp_rate_max = c(2,80),
  temperature_min = c(25,45),  temperature_max = c(25,45),
  spo2_min = c(30,100),  gcs_min = c(3,15),  urine_output_24h_ml = c(0,15000)
)

range_tbl <- rbindlist(lapply(intersect(names(phys_ranges), num_clust), function(v) {
  x <- dt[[v]]; lo <- phys_ranges[[v]][1]; hi <- phys_ranges[[v]][2]
  n_lo <- sum(x < lo, na.rm = TRUE); n_hi <- sum(x > hi, na.rm = TRUE)
  data.table(variable = v, lower = lo, upper = hi,
             n_below = n_lo, n_above = n_hi, n_out = n_lo + n_hi,
             pct_out = round((n_lo + n_hi)/N*100, 3))
}))
setorder(range_tbl, -pct_out)
cat("\n---- Values outside physiological bounds ----\n"); print(range_tbl)
fwrite(range_tbl, file.path(DIR_QC, "qc_range_check.csv"))


# ---------------------------------------------------------

# ---------------------------------------------------------

cor_mat <- cor(dt[, ..num_clust], method = "spearman", use = "pairwise.complete.obs")
fwrite(as.data.table(round(cor_mat,3), keep.rownames = "variable"),
       file.path(DIR_QC, "qc_spearman_matrix.csv"))

up <- which(upper.tri(cor_mat), arr.ind = TRUE)
high_corr <- data.table(
  var1 = rownames(cor_mat)[up[,1]], var2 = colnames(cor_mat)[up[,2]],
  rho  = round(cor_mat[up], 3))[abs(rho) > 0.8]
high_corr <- high_corr[order(-abs(rho))]
cat("\n---- |rho| > 0.8 variable pairs ----\n")
if (nrow(high_corr)) print(high_corr) else cat("No highly correlated pairs.\n")
fwrite(high_corr, file.path(DIR_QC, "qc_high_corr_pairs.csv"))


# ---------------------------------------------------------

# ---------------------------------------------------------

drop_high_miss <- miss_clust[tier == "3_drop_>40%", variable]
mice_vars      <- miss_clust[tier == "2_MICE_10-40%", variable]
keep_direct    <- miss_clust[tier == "1_usable_<10%", variable]

decision <- data.table(
  variable = clust_vars,
  miss_rate_pct = miss_tbl[match(clust_vars, variable), miss_rate_pct],
  suggest = fifelse(clust_vars %in% drop_high_miss, "consider_exclusion_high_missingness",
             fifelse(clust_vars %in% mice_vars,      "retain_with_MICE", "ready_without_imputation")))
if (nrow(high_corr))
  decision[variable %in% unique(c(high_corr$var1, high_corr$var2)),
           note := "Highly correlated pair; consider retaining one variable"]
setorder(decision, -miss_rate_pct)
fwrite(decision, file.path(DIR_QC, "qc_variable_decision.csv"))

cat("\n============ QC complete ============\n")
cat("Output directory: ", DIR_QC, "\n")
cat(sprintf("ready_without_imputation %d / MICE %d / consider exclusion %d\n",
            length(keep_direct), length(mice_vars), length(drop_high_miss)))
