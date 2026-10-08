if (!exists("D3R_ROOT")) source(file.path(D3R_SCRIPT_DIR, "41_d3r_sim_common_v1.R"))

tasks <- CJ(scenario = unname(D3R_SCENARIOS), repeat_id = seq_len(D3R_N_REP), sorted = FALSE)
tasks[, task_id := .I]
d3r_atomic_fwrite(tasks, file.path(D3R_PROVENANCE_DIR, "D3_sim_task_manifest.csv"), overwrite = TRUE)

run_task <- function(i) {
  scenario <- tasks$scenario[i]; repeat_id <- tasks$repeat_id[i]
  ck <- file.path(D3R_CHECKPOINT_DIR, sprintf("%s_rep%03d.rds", scenario, repeat_id))
  if (file.exists(ck)) return(readRDS(ck))
  ans <- d3r_run_repeat(scenario, repeat_id)
  tmp <- paste0(ck, ".tmp_", Sys.getpid())
  saveRDS(ans, tmp, version = 3)
  if (!file.rename(tmp, ck)) stop("Checkpoint write failed: ", ck)
  ans
}

if (D3R_WORKERS > 1L) {
  suppressPackageStartupMessages(library(future.apply))
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::multisession, workers = D3R_WORKERS)
  options(future.globals.maxSize = 2 * 1024^3)
  data.table::setDTthreads(1L)
  res <- future.apply::future_lapply(seq_len(nrow(tasks)), run_task, future.seed = TRUE)
  future::plan(future::sequential)
} else {
  res <- lapply(seq_len(nrow(tasks)), run_task)
}
repeat_level <- rbindlist(res, fill = TRUE)
setorder(repeat_level, scenario, repeat_id)
d3r_atomic_fwrite(repeat_level, file.path(D3R_TABLE_DIR, "D3_sim_repeat_level.csv"), overwrite = TRUE)
d3r_msg("Repeat execution complete: ", nrow(repeat_level), " rows.")
