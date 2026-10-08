options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/source_correction/00_correction_paths.R"))
agg <- aggregate
src <- figure_source
stopifnot(dir.exists(src))

read_src <- function(name) fread(file.path(src, name), showProgress = FALSE)
write_src <- function(x, name) fwrite(x, file.path(src, name))

# Update only Base + K2. Raw-EN/Raw-RF remain temperature-containing historical
# references; this script does not establish that their predictions are unchanged.
performance <- read_src("Figure4_panelC_L1_transport.csv")
external <- fread(file.path(agg, "L1_base_k2_external_performance.csv"))
stopifnot(nrow(performance) == 8L, nrow(external) == 2L)
for (outcome_name in c("mortality", "make")) {
  i <- which(performance$outcome == outcome_name & performance$model == "base_gen")
  z <- external[outcome == outcome_name]
  stopifnot(length(i) == 1L, nrow(z) == 1L,
            performance$n_external_stays[i] == z$n_external_stays)
  for (field in c("external_auc", "external_auc_ci_low", "external_auc_ci_high",
                  "external_auc_bootstrap_mcse", "calibration_slope",
                  "calibration_intercept", "brier_score")) {
    performance[[field]][i] <- z[[field]]
  }
  performance$transport_drop_internal_minus_external[i] <-
    performance$internal_cv_auc[i] - z$external_auc
  performance$model_label[i] <- "Base + K2"
  performance$emphasis[i] <- "Post hoc external temperature correction"
}
write_src(performance, "Figure4_panelC_L1_transport.csv")

calibration <- read_src("Figure4_panelC_L1_calibration.csv")
new_cal <- fread(file.path(agg, "L1_base_k2_calibration_bootstrap.csv"))
stopifnot(nrow(calibration) == 4L, nrow(new_cal) == 4L,
          all(calibration$model == "base_gen"))
for (i in seq_len(nrow(calibration))) {
  z <- new_cal[outcome == calibration$outcome[i] &
                 parameter == calibration$parameter[i]]
  stopifnot(nrow(z) == 1L, z$bootstrap_valid == 1000L,
            z$bootstrap_failed == 0L)
  for (field in c("estimate", "bootstrap_mean", "bootstrap_sd",
                  "bootstrap_mcse", "ci_low", "ci_high",
                  "bootstrap_requested", "bootstrap_valid", "bootstrap_failed",
                  "bootstrap_seed")) {
    calibration[[field]][i] <- z[[field]]
  }
}
calibration[, `:=`(prediction_input_md5 =
  unname(tools::md5sum(file.path(root, "private_local/corrected_external_predictions_local.csv"))),
  performance_input_md5 = unname(tools::md5sum(
    file.path(agg, "L1_base_k2_external_performance.csv"))),
  state_rule = "post hoc source correction; original reference range retained")]
write_src(calibration, "Figure4_panelC_L1_calibration.csv")

prev <- read_src("Figure4_pooled_L1_prevalence.csv")
i <- which(prev$cohort == "eICU")
stopifnot(length(i) == 1L, prev$n[i] == 17465L)
prev$c1_n[i] <- 3319L
prev$c1_prevalence[i] <- 3319 / 17465
write_src(prev, "Figure4_pooled_L1_prevalence.csv")

for (name in c("heldout_performance_point.csv",
               "heldout_performance_intervals.csv",
               "heldout_calibration_curves.csv",
               "hospital_holdout_split_summary.csv")) {
  ok <- file.copy(file.path(agg, name), file.path(src, name), overwrite = TRUE)
  stopifnot(isTRUE(ok))
}

risk <- fread(file.path(agg, "heldout_risk_distributions.csv"))
risk <- risk[model %in% c("base", "base_l1_k2") &
               method %in% c("no_update", "intercept_slope")]
screened <- risk[bin_low < 0.65]
tail <- risk[bin_low >= 0.65,
             .(bin_low = 0.65, bin_high = 1, stays = sum(stays)),
             by = .(model, method)]
screened <- rbindlist(list(screened, tail), use.names = TRUE)
setorder(screened, model, method, bin_low)
stopifnot(nrow(screened) == 56L, all(screened$stays >= 10L),
          all(screened[, sum(stays), by = .(model, method)]$V1 == 6478L))
write_src(screened, "heldout_risk_distribution_screened.csv")

cat("FIGURE4_SOURCE_READY: corrected L1 metrics, 96 screened hospitals, ",
    "four 14-bin screened risk distributions\n", sep = "")
