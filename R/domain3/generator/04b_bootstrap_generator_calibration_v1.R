## =====================================================================
## D3 regenerated deployable surrogate generator v1
## Patient-cluster bootstrap intervals for external base_gen calibration.
##
## This is a repeat-level post-processing step. It reads frozen external
## predictions, does not refit any outcome model, and does not change labels
## or the prespecified point-estimate transport gates.
## =====================================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot resolve script path from Rscript --file.")
script_path <- sub("^--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_d3_generator_config_v1.R"))
source(file.path(script_dir, "00_d3_generator_utils_v1.R"))

for (pkg in c("data.table")) require_namespace(pkg)
library(data.table)

prediction_path <- file.path(
  OUT_PRIVATE, "d3_generator_v1_external_predictions_local.csv"
)
performance_path <- file.path(
  OUT_TABLES, "d3_generator_v1_model_performance.csv"
)
output_path <- file.path(
  OUT_TABLES, "d3_generator_v1_base_gen_calibration_bootstrap.csv"
)

if (!file.exists(prediction_path)) {
  stop("Frozen external prediction file is missing: ", prediction_path)
}
if (!file.exists(performance_path)) {
  stop("External model performance file is missing: ", performance_path)
}

predictions <- fread(prediction_path)
performance <- fread(performance_path)

required_columns <- c(
  "uniquepid", "outcome", "outcome_value", "model", "prediction"
)
missing_columns <- setdiff(required_columns, names(predictions))
if (length(missing_columns)) {
  stop(
    "Frozen prediction file lacks required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

base_gen <- predictions[
  model == "base_gen" &
    !is.na(uniquepid) &
    !is.na(outcome_value) &
    !is.na(prediction)
]
if (!nrow(base_gen)) stop("No complete base_gen predictions were found.")

outcomes <- sort(unique(base_gen$outcome))
expected_outcomes <- c("make", "mortality")
if (!setequal(outcomes, expected_outcomes)) {
  stop(
    "Unexpected base_gen outcomes: ",
    paste(outcomes, collapse = ", "),
    "; expected: ",
    paste(expected_outcomes, collapse = ", ")
  )
}

bootstrap_rows <- list()
row_counter <- 0L

for (outcome_index in seq_along(outcomes)) {
  outcome_name <- outcomes[outcome_index]
  D <- base_gen[outcome == outcome_name]
  if (uniqueN(D$outcome_value) != 2L) {
    stop("Outcome is not binary in base_gen predictions: ", outcome_name)
  }

  point <- calibration_slope_intercept(D$outcome_value, D$prediction)
  bootstrap_seed <- MASTER_SEED + 9100L + outcome_index
  indices <- cluster_bootstrap_indices(
    D$uniquepid,
    B = MODEL_BOOT_B,
    seed = bootstrap_seed
  )

  bootstrap_values <- vapply(indices, function(ii) {
    calibration_slope_intercept(
      D$outcome_value[ii],
      D$prediction[ii]
    )
  }, numeric(2))

  for (parameter_name in c("slope", "intercept")) {
    values <- as.numeric(bootstrap_values[parameter_name, ])
    valid <- is.finite(values)
    interval <- if (any(valid)) {
      percentile_interval(values[valid])
    } else {
      c(NA_real_, NA_real_)
    }

    reference_estimate <- performance[
      outcome == outcome_name & model == "base_gen",
      if (parameter_name == "slope") calibration_slope else calibration_intercept
    ]
    if (length(reference_estimate) != 1L) {
      stop(
        "Expected one base_gen performance row for outcome: ",
        outcome_name
      )
    }
    if (!isTRUE(all.equal(
      unname(point[parameter_name]),
      unname(reference_estimate),
      tolerance = 1e-12
    ))) {
      stop(
        "Point estimate does not match the locked performance table for ",
        outcome_name, " ", parameter_name
      )
    }

    row_counter <- row_counter + 1L
    bootstrap_rows[[row_counter]] <- data.table(
      outcome = outcome_name,
      model = "base_gen",
      parameter = parameter_name,
      estimate = unname(point[parameter_name]),
      bootstrap_mean = mean(values[valid]),
      bootstrap_sd = stats::sd(values[valid]),
      bootstrap_mcse = stats::sd(values[valid]) / sqrt(sum(valid)),
      ci_low = interval[1],
      ci_high = interval[2],
      bootstrap_requested = MODEL_BOOT_B,
      bootstrap_valid = sum(valid),
      bootstrap_failed = sum(!valid),
      n_stays = nrow(D),
      n_patients = uniqueN(D$uniquepid),
      n_events = sum(D$outcome_value == 1L),
      resampling_unit = "uniquepid",
      bootstrap_seed = bootstrap_seed,
      prediction_input_md5 = hash_file(prediction_path),
      performance_input_md5 = hash_file(performance_path),
      prognostic_model_refit = FALSE,
      calibration_model_refit_each_bootstrap = TRUE,
      labels_regenerated = FALSE,
      state_rule = "prespecified_point_estimate_gate_unchanged"
    )
  }
}

result <- rbindlist(bootstrap_rows, use.names = TRUE)
result[, parameter_order := match(parameter, c("slope", "intercept"))]
setorder(result, outcome, parameter_order)
result[, parameter_order := NULL]
write_replace_csv(result, output_path)

writeLines(
  c(
    capture.output(sessionInfo()),
    "",
    paste0("Input prediction MD5: ", hash_file(prediction_path)),
    paste0("Input performance MD5: ", hash_file(performance_path)),
    paste0("Bootstrap replicates requested: ", MODEL_BOOT_B),
    "Resampling unit: uniquepid",
    "Outcome models refitted: FALSE",
    "Labels regenerated: FALSE",
    "Prespecified point-estimate state rule changed: FALSE"
  ),
  file.path(OUT_LOGS, "04b_bootstrap_generator_calibration_v1_sessionInfo.txt"),
  useBytes = TRUE
)

cat("Calibration bootstrap completed.\n")
cat("Output:", output_path, "\n")
print(result)
