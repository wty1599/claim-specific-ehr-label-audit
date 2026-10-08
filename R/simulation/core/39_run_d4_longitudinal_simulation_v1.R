options(stringsAsFactors = FALSE, warn = 1)

root <- file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation"
)
scripts <- file.path(root, "core")

source(file.path(scripts, "39a_d4_sim_freeze_spec_v1.R"))
source(file.path(scripts, "39d_d4_sim_apply_support_engine_v1.R"))
source(file.path(scripts, "39e_d4_sim_operating_characteristics_v1.R"))
source(file.path(scripts, "39f_d4_sim_figures_tables_v1.R"))
source(file.path(scripts, "39g_d4_sim_independent_qc_v1.R"))

sim_atomic_write_lines(
  capture.output(sessionInfo()),
  file.path(SIM_LOG_DIR, "sessionInfo.txt")
)
sim_msg("D4 longitudinal simulation completed: ", SIM_OUTPUT_ROOT)
