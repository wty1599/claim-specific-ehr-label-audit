


source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages(library(data.table))
mimic <- fread(PATH_FINAL_FULL)
OL <- readRDS(file.path(DIR_MODEL,"outcomes_local.rds")); oc <- as.data.table(OL$outcomes)
cat("=== MIMIC final cohort ===\n")
cat("Final EHR audit cohort n =", nrow(mimic), "\n")
cat("MAKE-30 events:", sum(oc$make30), sprintf("(%.1f%%)\n", 100*mean(oc$make30)))
cat("30-day mortalityevents:", sum(oc$mortality_30d), sprintf("(%.1f%%)\n", 100*mean(oc$mortality_30d)))
cat("90-day mortalityevents:", sum(oc$mortality_90d), sprintf("(%.1f%%)\n", 100*mean(oc$mortality_90d)))
cat("\n[Upstream cohort counts are defined by the project_sa_aki source SQL]\n")
cat("  Source SQL: SELECT COUNT(*) FROM project_sa_aki.<intermediate_table>;\n")
