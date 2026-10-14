## ============================================================


## ============================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))

REPS    <- 1000
CONS_N  <- 5000
FINAL_K <- 2

source(file.path(PROJ_ROOT, "R", "common", "cluster_engine.R"))

TAG  <- "sens_inclusive"
Xset <- readRDS(file.path(DIR_MODEL, "Xs_sets.rds"))[[TAG]]

run_consensus_selection(Xset, TAG)
if (!is.na(FINAL_K)) finalize_set(Xset, TAG, FINAL_K)
  
