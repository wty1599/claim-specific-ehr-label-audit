## ============================================================







## ============================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({ library(data.table); library(pheatmap) })
set.seed(20250101)

Xs  <- readRDS(file.path(DIR_MODEL, "Xs_sets.rds"))
OL  <- readRDS(file.path(DIR_MODEL, "outcomes_local.rds"))
oc  <- as.data.table(OL$outcomes); id_col <- OL$id_col

Xp   <- Xs$primary
feat <- colnames(Xp); if (is.null(feat)) feat <- names(attr(Xp, "scaled:center"))
colnames(Xp) <- feat
stopifnot(nrow(Xp) == nrow(oc))

run_K <- function(K){
  km  <- kmeans(Xp, centers = K, nstart = 25, iter.max = 100)
  raw <- km$cluster

  ord <- names(sort(tapply(oc$mortality_30d, raw, mean), decreasing = TRUE))
  remap <- setNames(seq_along(ord), ord)
  lab <- remap[as.character(raw)]


  prop <- as.data.table(table(cluster = lab))[, .(cluster, n = N,
             pct = round(100*N/sum(N), 2))]
  cat(sprintf("\n===== K = %d =====\n", K)); cat("[Proportions]\n"); print(prop)


  prof <- sapply(sort(unique(lab)), function(k) colMeans(Xp[lab==k,,drop=FALSE]))
  colnames(prof) <- paste0("C", sort(unique(lab)))
  fwrite(data.table(feature = rownames(prof), prof),
         file.path(DIR_OUTPUT, sprintf("profile_zscore_primary_K%d.csv", K)))
  rng <- max(abs(prof))
  pheatmap(prof, cluster_cols = FALSE, cluster_rows = TRUE,
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(51),
           breaks = seq(-rng, rng, length.out = 52),
           display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6,
           main = sprintf("Standardized profile (primary, K=%d)", K),
           filename = file.path(DIR_FIG, sprintf("profile_heatmap_primary_K%d.png", K)),
           width = 4 + 0.5*K, height = 9)


  d <- cbind(cluster = lab, oc)
  outc <- d[, .(n = .N,
                make30_pct   = round(100*mean(make30, na.rm=TRUE), 1),
                mort30_pct   = round(100*mean(mortality_30d, na.rm=TRUE), 1),
                mort90_pct   = round(100*mean(mortality_90d, na.rm=TRUE), 1),
                aki_any_pct  = round(100*mean(aki_kdigo_7d, na.rm=TRUE), 1),
                aki_stage3_pct = round(100*mean(aki_stage_kdigo_7d>=3, na.rm=TRUE), 1),
                rrt7d_pct    = round(100*mean(rrt_7d, na.rm=TRUE), 1),
                persist_rd_pct = round(100*mean(persistent_rd_30d, na.rm=TRUE), 1),
                icu_los_med  = round(median(icu_los_days, na.rm=TRUE), 2)),
              by = cluster][order(cluster)]
  cat("[Outcomes]\n"); print(outc)
  fwrite(outc, file.path(DIR_OUTPUT, sprintf("outcomes_by_cluster_primary_K%d.csv", K)))


  data.table(setNames(list(oc[[id_col]], lab), c(id_col, sprintf("cluster_k%d", K))))
}

lab2 <- run_K(2)
lab3 <- run_K(3)


labs <- merge(lab2, lab3, by = id_col)
fwrite(labs, file.path(DIR_OUTPUT, "cluster_labels_primary_k2_k3.csv"))
saveRDS(labs, file.path(DIR_MODEL, "labels_primary.rds"))

cat("\nPartition diagnostics complete:\n")
cat("  Cluster proportions (smallest cluster >~5%?), profile heatmaps, and outcome separation\n")
cat("  Figures: output/fig/profile_heatmap_primary_K2.png / _K3.png\n")
