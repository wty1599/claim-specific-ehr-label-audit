## ============================================================


## ============================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
source(file.path(PROJ_ROOT, "R", "common", "00b_cluster_utils.R"))
need_pkg(c("data.table", "mclust")); library(data.table)

lab_p <- readRDS(file.path(DIR_MODEL, "labels_primary.rds"))
lab_r <- readRDS(file.path(DIR_MODEL, "labels_sens_restricted.rds"))
lab_i <- readRDS(file.path(DIR_MODEL, "labels_sens_inclusive.rds"))
fsets <- readRDS(file.path(DIR_MODEL, "feature_sets.rds"))

cmp <- rbindlist(list(
  data.table(sensitivity = "sens_restricted", n_vars = length(fsets$sens_restricted),
             ARI_vs_primary = round(mclust::adjustedRandIndex(lab_p, lab_r), 3)),
  data.table(sensitivity = "sens_inclusive",  n_vars = length(fsets$sens_inclusive),
             ARI_vs_primary = round(mclust::adjustedRandIndex(lab_p, lab_i), 3))
))
fwrite(cmp, file.path(DIR_OUTPUT, "cluster_sensitivity_vs_primary.csv"))
cat("\nPrimary analysis (33 features) robustness comparison:\n"); print(cmp)
cat("\nCross-tabulation primary × restricted:\n"); print(table(primary = lab_p, restricted = lab_r))
cat("\nCross-tabulation primary × inclusive:\n");  print(table(primary = lab_p, inclusive  = lab_i))
cat("\n(ARI>0.7 is the historical feature-selection robustness criterion)\n")


oc <- readRDS(file.path(DIR_MODEL, "outcomes_local.rds"))
assign_tbl <- data.table(id = oc$outcomes[[oc$id_col]],
                         cluster_primary = lab_p,
                         cluster_restricted = lab_r,
                         cluster_inclusive = lab_i)
setnames(assign_tbl, "id", oc$id_col)
fwrite(assign_tbl, file.path(DIR_OUTPUT, "cluster_assignments.csv"))
cat("\nAssignments saved locally: ", file.path(DIR_OUTPUT, "cluster_assignments.csv"), "\n")
