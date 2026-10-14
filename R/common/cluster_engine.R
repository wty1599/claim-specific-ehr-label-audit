## ============================================================




## ============================================================
source("00b_cluster_utils.R")
need_pkg(c("cluster", "fpc", "mclust", "data.table"))
library(data.table)
library(mclust)

if (!exists("REPS"))       REPS       <- 300
if (!exists("MAXK"))       MAXK       <- 6
if (!exists("P_ITEM"))     P_ITEM     <- 0.80
if (!exists("P_FEAT"))     P_FEAT     <- 1.00    # Feature resampling fraction (CCP pFeature=1).
if (!exists("CONS_N"))     CONS_N     <- 5000
if (!exists("CKPT_EVERY")) CKPT_EVERY <- 25
if (!exists("N_SUB"))      N_SUB      <- 5000
if (!exists("SEED"))       SEED       <- 20250101


consensus_one_K <- function(Xc, K, reps, p_item, p_feat, seed, ckpt_file, ckpt_every) {
  ns <- nrow(Xc); p <- ncol(Xc)
  reps_store <- load_ckpt(ckpt_file)
  
  if (is.null(reps_store)) {
    reps_store <- vector("list", reps)
  } else {
    if (!is.list(reps_store)) {
      stop("Invalid checkpoint file format: ", ckpt_file)
    }
    
    old_len <- length(reps_store)
    
    if (old_len < reps) {
      length(reps_store) <- reps
      log_msg(sprintf("  Checkpoint expanded: %d -> %d reps", old_len, reps))
      save_ckpt(reps_store, ckpt_file)
    }
    
    if (old_len > reps) {
      reps_store <- reps_store[seq_len(reps)]
      log_msg(sprintf("  Checkpoint truncated: %d -> %d reps", old_len, reps))
      save_ckpt(reps_store, ckpt_file)
    }
  }
  done <- sum(!vapply(reps_store, is.null, logical(1)))
  if (done > 0) log_msg(sprintf("  K=%d resuming; completed %d/%d", K, done, reps))
  if (done < reps) {
    pb <- txtProgressBar(min = 0, max = reps, initial = done, style = 3)
    for (r in seq_len(reps)) {
      if (!is.null(reps_store[[r]])) { setTxtProgressBar(pb, r); next }
      set.seed((seed + K * 100003L + r) %% 2147483647L)
      items <- sort(sample.int(ns, max(2L, floor(p_item * ns))))
      feats <- sort(sample.int(p,  max(2L, floor(p_feat * p))))
      km <- tryCatch(kmeans(Xc[items, feats, drop = FALSE], centers = K,
                            nstart = 1, iter.max = 50, algorithm = "Lloyd"),
                     error = function(e) NULL)
      reps_store[[r]] <- if (is.null(km)) list(idx = integer(0), lab = integer(0))
                         else list(idx = items, lab = as.integer(km$cluster))
      if (r %% ckpt_every == 0 || r == reps) save_ckpt(reps_store, ckpt_file)
      setTxtProgressBar(pb, r)
    }
    close(pb)
  }
  M <- matrix(0L, ns, ns); I <- matrix(0L, ns, ns)
  for (r in seq_len(reps)) {
    rr <- reps_store[[r]]; if (length(rr$idx) < 2) next
    I[rr$idx, rr$idx] <- I[rr$idx, rr$idx] + 1L
    M[rr$idx, rr$idx] <- M[rr$idx, rr$idx] + outer(rr$lab, rr$lab, "==")
  }
  C <- matrix(0, ns, ns); nz <- I > 0; C[nz] <- M[nz] / I[nz]; diag(C) <- 1
  ut <- upper.tri(C); vals <- C[ut & nz]
  pac <- mean(vals > 0.1 & vals < 0.9)
  d   <- as.dist(1 - C)
  cls <- cutree(hclust(d, method = "average"), k = K)
  sil <- mean(cluster::silhouette(cls, d)[, "sil_width"])
  auc <- mean(ecdf(vals)(seq(0, 1, 0.01)))
  rm(M, I, C); gc()
  list(K = K, pac = pac, sil = sil, auc = auc)
}


run_consensus_selection <- function(Xset, tag) {
  set.seed(SEED)
  n <- nrow(Xset)
  cons_idx <- if (n > CONS_N) sort(sample.int(n, CONS_N)) else seq_len(n)
  Xc <- Xset[cons_idx, , drop = FALSE]; ns <- nrow(Xc)
  res <- list()
  for (K in 2:MAXK) {
    log_msg(sprintf("=== [%s] Consensus K=%d (reps=%d, ns=%d, pFeat=%.2f) ===", tag, K, REPS, ns, P_FEAT))
    t0 <- Sys.time()
    ck <- file.path(CKPT_DIR, sprintf("cons_%s_K%d_ns%d_seed%d.rds", tag, K, ns, SEED))
    res[[as.character(K)]] <- consensus_one_K(Xc, K, REPS, P_ITEM, P_FEAT, SEED, ck, CKPT_EVERY)
    log_msg(sprintf("  Elapsed %.1f minutes", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  ksel <- data.frame(K = 2:MAXK, PAC = sapply(res, `[[`, "pac"),
                     silhouette = sapply(res, `[[`, "sil"), cdf_auc = sapply(res, `[[`, "auc"))
  ksel$delta_area <- c(NA, diff(ksel$cdf_auc) / head(ksel$cdf_auc, -1))
  fwrite(ksel, file.path(DIR_OUTPUT, sprintf("k_selection_%s.csv", tag)))
  png(file.path(DIR_FIG, sprintf("k_selection_%s.png", tag)), 1400, 600, res = 130)
  op <- par(mfrow = c(1, 2))
  plot(ksel$K, ksel$PAC, type = "b", pch = 19, xlab = "K", ylab = "PAC",
       main = sprintf("[%s] PAC (lower is better)", tag)); grid()
  plot(ksel$K, ksel$silhouette, type = "b", pch = 19, xlab = "K", ylab = "silhouette",
       main = sprintf("[%s] Silhouette", tag)); grid()
  par(op); dev.off()
  log_msg(sprintf("[%s] k selection table:", tag)); print(ksel)
  ksel
}


finalize_set <- function(Xset, tag, k) {
  n <- nrow(Xset)
  log_msg(sprintf("=== [%s] Part B: full-cohort k-means K=%d ===", tag, k))
  set.seed(SEED)
  lab <- kmeans(Xset, centers = k, nstart = 50, iter.max = 100)$cluster


  sz <- as.data.frame(table(cluster = lab)); sz$pct <- round(100 * sz$Freq / n, 2)
  sz$flag_small <- sz$pct < 5
  fwrite(sz, file.path(DIR_OUTPUT, sprintf("cluster_sizes_%s.csv", tag))); print(sz)


  prof <- t(sapply(sort(unique(lab)), function(kk) colMeans(Xset[lab == kk, , drop = FALSE])))
  rownames(prof) <- paste0("C", sort(unique(lab)))
  fwrite(data.frame(cluster = rownames(prof), prof),
         file.path(DIR_OUTPUT, sprintf("zprofile_%s.csv", tag)))
  png(file.path(DIR_FIG, sprintf("heatmap_%s.png", tag)), 1500, 900, res = 130)
  heatmap(prof, Colv = NA, Rowv = NA, scale = "none",
          col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
          margins = c(10, 5), main = sprintf("[%s] Cluster x feature z", tag))
  dev.off()


  idx <- if (is.na(N_SUB)) seq_len(n) else { set.seed(2); sort(sample.int(n, min(N_SUB, n))) }
  sub <- Xset[idx, ]
  lab_gmm <- mclust::Mclust(sub, G = k, verbose = FALSE)$classification
  lab_hc  <- cutree(hclust(dist(sub), method = "ward.D2"), k = k)
  ari <- data.frame(pair = c("KM_vs_GMM", "KM_vs_HC", "GMM_vs_HC"),
                    ARI = round(c(mclust::adjustedRandIndex(lab[idx], lab_gmm),
                                  mclust::adjustedRandIndex(lab[idx], lab_hc),
                                  mclust::adjustedRandIndex(lab_gmm, lab_hc)), 3))
  fwrite(ari, file.path(DIR_OUTPUT, sprintf("method_agreement_%s.csv", tag))); print(ari)


  log_msg(sprintf("[%s] bootstrap stability (B=100)…", tag))
  cb <- fpc::clusterboot(Xset, B = 100, bootmethod = "boot",
                         clustermethod = fpc::kmeansCBI, krange = k,
                         seed = SEED, count = TRUE)
  stab <- data.frame(cluster = seq_along(cb$bootmean),
                     jaccard = round(cb$bootmean, 3), pass_0.75 = cb$bootmean > 0.75)
  fwrite(stab, file.path(DIR_OUTPUT, sprintf("stability_%s.csv", tag))); print(stab)


  oc <- readRDS(file.path(DIR_MODEL, "outcomes_local.rds"))$outcomes
  oc$cluster <- factor(lab)
  cat_v <- intersect(c("make30","rrt_7d","rrt_30d","persistent_rd_30d","icu_mortality",
                       "hospital_mortality","mortality_28d","mortality_30d","mortality_90d",
                       "aki_stage_kdigo_7d"), names(oc))
  ct <- rbindlist(lapply(cat_v, function(v) {
    tb <- table(oc$cluster, oc[[v]]); pr <- round(prop.table(tb, 1), 3)
    p  <- tryCatch(chisq.test(tb)$p.value, error = function(e) NA)
    data.table(var = v, cluster = rownames(pr), as.data.frame.matrix(pr), p_chisq = round(p, 4))
  }), fill = TRUE)
  fwrite(ct, file.path(DIR_OUTPUT, sprintf("outcome_categorical_%s.csv", tag)))
  num_v <- intersect(c("icu_los_days", "hospital_los_days"), names(oc))
  if (length(num_v)) {
    nt <- rbindlist(lapply(num_v, function(v) {
      ag <- aggregate(oc[[v]], list(oc$cluster), median, na.rm = TRUE)
      p  <- tryCatch(kruskal.test(oc[[v]] ~ oc$cluster)$p.value, error = function(e) NA)
      data.table(var = v, cluster = ag[[1]], median = round(ag[[2]], 2), p_kruskal = round(p, 4))
    }))
    fwrite(nt, file.path(DIR_OUTPUT, sprintf("outcome_numeric_%s.csv", tag)))
  }

  saveRDS(lab, file.path(DIR_MODEL, sprintf("labels_%s.rds", tag)))
  log_msg(sprintf("[%s] Part B complete", tag))
  invisible(lab)
}
