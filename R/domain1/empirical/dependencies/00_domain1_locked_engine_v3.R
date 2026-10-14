## Read-only loader and wrappers for the locked S1/S2/S7 simulation engine.
## The archived script is parsed in two definition-only blocks; its batch run,
## future plan, directories, and S7 generation gate are never executed.

load_locked_main_s7_engine <- function(path = D1_INPUTS$main_s1_s2_s7_engine) {
  path <- assert_file(path, "locked S1/S2/S7 engine")
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  settings_start <- grep("^GLOBAL_SEED <-", lines)[1L]
  settings_stop <- grep("^## 1\\. Packages and logging", lines)[1L] - 1L
  definitions_start <- grep("^## ---- Robust table binding", lines)[1L]
  definitions_stop <- grep("^## 10\\. Batch run", lines)[1L] - 1L
  marks <- c(settings_start, settings_stop, definitions_start, definitions_stop)
  if (any(!is.finite(marks)) || settings_start >= settings_stop ||
      definitions_start >= definitions_stop) {
    stop("Could not isolate definition-only blocks in locked engine.", call. = FALSE)
  }

  env <- new.env(parent = globalenv())
  env$log_msg <- function(...) invisible(NULL)
  ## The settings block constructs archived candidate paths from PROJECT_DIR,
  ## but the definition-only loader intentionally skips the earlier directory
  ## creation lines. Supply a harmless read-only base solely so those constants
  ## can be parsed; none of the candidate files are read by this package.
  env$PROJECT_DIR <- dirname(path)
  eval(parse(text = paste(lines[settings_start:settings_stop], collapse = "\n")), env)
  ## The revised rerun never calls the archived S7 generation gate.
  env$S7_REQUIRE_EXISTING_D1_COMBINED_RULE <- FALSE
  eval(parse(text = paste(lines[definitions_start:definitions_stop], collapse = "\n")), env)

  required <- c(
    "feature_names", "generate_dataset", "make_mice_imputations",
    "simulate_latent_features", "simulate_latent_features_s7",
    "to_clinical_scale", "winsor_fit", "winsor_apply",
    "standardize_train", "standardize_apply"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(missing)) {
    stop("Locked engine extraction missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(env$feature_names) != 33L) {
    stop("Locked main/S7 engine does not define 33 features.", call. = FALSE)
  }
  env
}

preprocess_train_evaluation <- function(engine, train_dt, evaluation_dt) {
  vars <- engine$feature_names
  wb <- engine$winsor_fit(train_dt, vars)
  train_w <- engine$winsor_apply(train_dt, wb)
  evaluation_w <- engine$winsor_apply(evaluation_dt, wb)
  std <- engine$standardize_train(train_w, vars)
  X_evaluation <- engine$standardize_apply(
    evaluation_w, vars, std$center, std$scale
  )
  list(
    X_train = std$Xz,
    X_evaluation = X_evaluation,
    winsor_bounds = wb,
    center = std$center,
    scale = std$scale
  )
}

run_locked_mice_with_warning_audit <- function(engine, ...) {
  warning_text <- character()
  imputations <- withCallingHandlers(
    engine$make_mice_imputations(...),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    imputations = imputations,
    warning_count = length(warning_text),
    warning_text = paste(unique(warning_text), collapse = " | ")
  )
}

simulate_complete_d1_pair <- function(engine, dgm_family, n_train, n_evaluation,
                                      delta, seed) {
  dgm_family <- match.arg(dgm_family, c("S2_like_K3", "S7_like_K2", "S1_continuum"))

  make_one <- function(n, local_seed) {
    set.seed(local_seed)
    if (dgm_family == "S1_continuum") {
      true_subtype <- rep("continuum", n)
      latent <- engine$simulate_latent_features(
        n = n, severity_mean = 0, severity_sd = 1,
        true_subtype = true_subtype, subtype_delta = 0,
        external_shift = NULL, severity_signal_multiplier = 1,
        organ_signal_multiplier = 1, noise_sd = 0.65,
        incremental_group = NULL, incremental_delta = 0
      )
    } else if (dgm_family == "S2_like_K3") {
      true_subtype <- sample(
        c("metabolic", "renal_dysfunction", "shock"),
        size = n, replace = TRUE, prob = c(0.35, 0.35, 0.30)
      )
      latent <- engine$simulate_latent_features(
        n = n, severity_mean = 0, severity_sd = 1,
        true_subtype = true_subtype, subtype_delta = delta,
        external_shift = NULL,
        severity_signal_multiplier = engine$S2_SEVERITY_SIGNAL_MULTIPLIER,
        organ_signal_multiplier = engine$S2_ORGAN_SIGNAL_MULTIPLIER,
        noise_sd = engine$S2_NOISE_SD,
        incremental_group = NULL, incremental_delta = 0
      )
    } else {
      old_delta <- engine$S7_CLUSTER_DELTA
      on.exit(assign("S7_CLUSTER_DELTA", old_delta, envir = engine), add = TRUE)
      assign("S7_CLUSTER_DELTA", delta, envir = engine)
      latent <- engine$simulate_latent_features_s7(n)
      true_subtype <- latent$true_subtype
    }
    features <- data.table::copy(data.table::as.data.table(
      engine$to_clinical_scale(latent$Z)
    ))
    features[, `:=`(
      true_subtype = true_subtype,
      severity_z = latent$severity_z
    )]
    features
  }

  train <- make_one(n_train, seed + 1L)
  evaluation <- make_one(n_evaluation, seed + 2L)
  prep <- preprocess_train_evaluation(engine, train, evaluation)
  c(prep, list(train_dt = train, evaluation_dt = evaluation))
}

simulate_mice_d1_pair <- function(engine, scenario, repeat_id, seed,
                                  n_train = 5000L,
                                  n_evaluation = D1_EVALUATION_N) {
  scenario <- match.arg(scenario, c(
    "S1_pure_severity_continuum",
    "S2_true_discrete_subtypes",
    "S7_discrete_structure_outcome_null"
  ))
  train_raw <- engine$generate_dataset(
    scenario = scenario, cohort = "derivation", n = n_train,
    repeat_id = repeat_id, seed = seed + 1L
  )
  ## Use derivation mode for the independent evaluation cohort so the D1 test
  ## assesses same-distribution structure rather than Domain 3 drift.
  evaluation_raw <- engine$generate_dataset(
    scenario = scenario, cohort = "derivation", n = n_evaluation,
    repeat_id = repeat_id, seed = seed + 2L
  )
  train_mice <- run_locked_mice_with_warning_audit(
    engine,
    train_raw, scenario, repeat_id, cohort = "derivation",
    m = D1_MICE_M, maxit = D1_MICE_MAXIT, seed = seed + 101L
  )
  evaluation_mice <- run_locked_mice_with_warning_audit(
    engine,
    evaluation_raw, scenario, repeat_id, cohort = "same_distribution_evaluation",
    m = D1_MICE_M, maxit = D1_MICE_MAXIT, seed = seed + 202L
  )
  train_imp <- train_mice$imputations
  evaluation_imp <- evaluation_mice$imputations
  train <- train_imp[[D1_PRIMARY_IMPUTATION_ID]]
  evaluation <- evaluation_imp[[D1_PRIMARY_IMPUTATION_ID]]
  prep <- preprocess_train_evaluation(engine, train, evaluation)
  c(prep, list(
    train_dt = train,
    evaluation_dt = evaluation,
    imputation_m = D1_MICE_M,
    primary_imputation_id = D1_PRIMARY_IMPUTATION_ID,
    mice_warning_count = train_mice$warning_count + evaluation_mice$warning_count,
    mice_warning_text = paste(
      c(train_mice$warning_text, evaluation_mice$warning_text)[
        nzchar(c(train_mice$warning_text, evaluation_mice$warning_text))
      ],
      collapse = " | "
    )
  ))
}
