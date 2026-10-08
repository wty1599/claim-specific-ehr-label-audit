#!/usr/bin/env Rscript

Sys.setenv(D1_SENS_ANALYSIS = "SA01")
project_root <- Sys.getenv(
  "MIMIC_IV_ROOT", unset = Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd())
)
source(
  file.path(
    project_root, "d1_current_rule_sensitivity_20260803", "scripts",
    "00_d1_current_sensitivity_engine.R"
  ),
  encoding = "UTF-8"
)

