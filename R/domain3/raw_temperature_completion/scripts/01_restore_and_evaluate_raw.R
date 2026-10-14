options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})
for (pkg in c("glmnet", "ranger")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
if (as.character(getRversion()) != "4.6.0") stop("R 4.6.0 is required.")
for (cat in c("LC_COLLATE", "LC_CTYPE", "LC_MONETARY", "LC_TIME")) {
  Sys.setlocale(cat, "Chinese (Simplified)_China.utf8")
}
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(file_arg) == 1L)
entry <- normalizePath(sub("^--file=", "", file_arg), winslash = "/")
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(entry), "../../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/raw_temperature_completion/00_raw_paths.R"))
args <- commandArgs(TRUE)
stage <- sub("^--stage=", "", grep("^--stage=", args, value = TRUE))
stopifnot(length(stage) == 1L, stage %in% c("preflight", "recover", "evaluate"))
SPEC_SHA <- "2CA190DCAF14FEE54A76606C73CFA80D9A872E13259B408A30BEE5799789E905"
stopifnot(all(file.exists(PATHS)))
sha <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
stopifnot(sha(PATHS[["specification"]]) == SPEC_SHA)
stopifnot(sha(file.path(RAW_CODE_ROOT, "spec/CALIBRATION_DEFINITION_ADDENDUM_20261008.md")) ==
  "7CF08AD3B936CE087DC3BAAFA800FDE9871BA251B8B33670E5553442914C4C8D")
safe_csv <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  fwrite(x, path)
}
safe_rds <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  saveRDS(x, path, version = 3)
}
safe_json <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  write_json(x, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
}
now <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
record_path <- file.path(ROOT, "provenance/pre_run_record.json")
manifest_path <- file.path(ROOT, "provenance/input_manifest.csv")
input_snapshot <- function() {
  f <- file.info(PATHS)
  data.table(role = names(PATHS), path = unname(PATHS), bytes = f$size,
             modified_utc = format(f$mtime, tz = "UTC", usetz = TRUE),
             sha256 = vapply(PATHS, sha, character(1)))
}
check_inputs <- function() {
  old <- fread(manifest_path)
  actual <- input_snapshot()
  stopifnot(identical(old$role, actual$role), identical(old$sha256, actual$sha256))
  invisible(TRUE)
}

# Original constants and pure helpers; the original configuration is never sourced.
MASTER_SEED <- 20260728L
MODEL_BOOT_B <- 1000L
MODEL_CV_FOLDS <- 10L
RF_TREES <- 300L
ENET_ALPHA <- 0.5
SOURCE_LABEL_LEVELS <- c("1", "2")
FEATURES33 <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min",
  "abs_lymphocytes_min", "lactate_max", "ph_min", "pco2_max",
  "calcium_min", "calcium_max", "ptt_max", "inr_max",
  "temperature_min", "temperature_max", "urine_output_24h_ml",
  "glucose_max", "aniongap_max", "potassium_min", "potassium_max",
  "hemoglobin_min", "sodium_min", "sodium_max", "wbc_max",
  "platelets_min", "bicarbonate_min", "chloride_min", "chloride_max",
  "bun_max", "creatinine_max", "resp_rate_max", "gcs_min", "spo2_min",
  "heart_rate_max", "mbp_min"
)
BASE_NUM <- c("age", "sex", "sofa_score", "aki_stage_0_24h")
MODEL_KINDS <- c("base", "base_gen", "feat_pen", "feat_rf")
RAW_KINDS <- c("feat_pen", "feat_rf")
OUTCOMES <- list(
  mortality = list(mimic = "hospital_mortality", eicu = "hosp_mortality"),
  make = list(mimic = "make30", eicu = "make_hosp")
)
source(PATHS[["original_utils"]])

fit_model <- function(kind, y, base, X, cluster) {
  switch(
    kind,
    base = stats::glm(
      y ~ ., data = data.frame(y = y, base), family = stats::binomial()
    ),
    base_gen = stats::glm(
      y ~ .,
      data = data.frame(y = y, base, cluster_k2 = factor(cluster)),
      family = stats::binomial()
    ),
    feat_pen = glmnet::cv.glmnet(
      cbind(base, X), y, family = "binomial", alpha = ENET_ALPHA,
      nfolds = 5
    ),
    feat_rf = ranger::ranger(
      x = as.data.frame(cbind(base, X)),
      y = factor(y),
      probability = TRUE,
      num.trees = RF_TREES,
      num.threads = 0,
      seed = 1L
    ),
    stop("Unknown model kind: ", kind)
  )
}

predict_model <- function(kind, fit, base, X, cluster) {
  switch(
    kind,
    base = as.numeric(stats::predict(
      fit, newdata = data.frame(base), type = "response"
    )),
    base_gen = as.numeric(stats::predict(
      fit,
      newdata = data.frame(base, cluster_k2 = factor(
        cluster, levels = levels(fit$model$cluster_k2)
      )),
      type = "response"
    )),
    feat_pen = as.numeric(stats::predict(
      fit, newx = cbind(base, X), s = "lambda.min", type = "response"
    )),
    feat_rf = as.numeric(stats::predict(
      fit, data = as.data.frame(cbind(base, X))
    )$predictions[, "1"]),
    stop("Unknown model kind: ", kind)
  )
}

make_stratified_folds <- function(y, k) {
  folds <- integer(length(y))
  for (class_value in c(0L, 1L)) {
    ii <- which(y == class_value)
    folds[ii] <- sample(rep_len(seq_len(k), length(ii)))
  }
  folds
}

internal_cv_predictions <- function(y, base, X, cluster, model_kinds) {
  folds <- make_stratified_folds(y, MODEL_CV_FOLDS)
  out <- matrix(NA_real_, nrow = length(y), ncol = length(model_kinds))
  colnames(out) <- model_kinds
  for (kind in model_kinds) {
    for (f in seq_len(MODEL_CV_FOLDS)) {
      train <- folds != f
      test <- folds == f
      fit <- fit_model(kind, y[train], base[train, , drop = FALSE],
                       X[train, , drop = FALSE], cluster[train])
      out[test, kind] <- predict_model(
        kind, fit, base[test, , drop = FALSE], X[test, , drop = FALSE],
        cluster[test]
      )
      cat("CV_PROGRESS", kind, f, now(), "\n")
      flush.console()
    }
  }
  if (anyNA(out)) stop("Internal cross-validated predictions are incomplete.")
  attr(out, "folds") <- folds
  out
}

bootstrap_auc_ci <- function(y, p, cluster_id, B) {
  cluster_id <- as.character(cluster_id)
  clusters <- unique(cluster_id)
  rows_by_cluster <- split(seq_along(cluster_id), cluster_id)
  values <- vapply(seq_len(B), function(b) {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    ii <- unlist(rows_by_cluster[sampled], use.names = FALSE)
    fast_auc(y[ii], p[ii])
  }, numeric(1))
  ci <- percentile_interval(values)
  c(low = ci[1], high = ci[2], mcse = sd(values, na.rm = TRUE) /
      sqrt(sum(is.finite(values))))
}

# Parse, but never evaluate, the original script to compare copied function bodies.
original_definition <- function(name) {
  e <- parse(PATHS[["original_script"]])
  hits <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                   identical(x[[2]], as.name(name)), as.list(e))
  stopifnot(length(hits) == 1L)
  eval(hits[[1]][[3]], envir = new.env(parent = globalenv()))
}
for (name in c("fit_model", "predict_model", "make_stratified_folds", "bootstrap_auc_ci")) {
  stopifnot(identical(body(get(name)), body(original_definition(name))))
}

if (stage == "preflight") {
  stopifnot(as.character(packageVersion("glmnet")) == "5.0",
            as.character(packageVersion("ranger")) == "0.18.0",
            as.character(packageVersion("data.table")) == "1.18.4")
  safe_csv(input_snapshot(), manifest_path)
  safe_json(list(
    recorded_utc = now(), analysis_date_local = "2026-10-08",
    specification_sha256 = SPEC_SHA, script_sha256 = sha(entry),
    post_hoc = TRUE, raw_recipe = "archived object$recipe (L2 coordinates)",
    initial_source_seed = 20240601L, original_master_seed = MASTER_SEED,
    rng_kind = RNGkind(), original_model_order = MODEL_KINDS,
    original_cv_folds = MODEL_CV_FOLDS, elastic_net_alpha = ENET_ALPHA,
    elastic_net_inner_folds = 5L, elastic_net_prediction_lambda = "lambda.min",
    random_forest_trees = RF_TREES, random_forest_seed = 1L,
    random_forest_threads = 0L, prediction_tolerance = 1e-12,
    bootstrap_replicates = 1000L, bootstrap_unit = "uniquepid",
    paired_seeds = c(mortality = MASTER_SEED + 5002L, make = MASTER_SEED + 5006L),
    uncertainty = "patient percentile bootstrap conditional on fixed recovered source fits",
    calibration_slope = "joint glm(y ~ logit(p)); probability clamp 1e-6",
    calibration_intercept = "offset-only glm(y ~ 1, offset=logit(p)); original display (CITL)",
    joint_calibration_intercept = "separately named intercept from glm(y ~ logit(p))",
    external_model_fitting = FALSE, clustering = FALSE, mice = FALSE,
    R_version = R.version.string,
    packages = lapply(c("data.table", "glmnet", "ranger", "digest", "jsonlite"), function(p) {
      list(package = p, version = as.character(packageVersion(p)))
    })
  ), record_path)
  writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs/preflight_sessionInfo.txt"))
  cat("PREFLIGHT_RECORD_COMPLETE", now(), "\n")
  quit(status = 0L)
}

stopifnot(file.exists(record_path), file.exists(manifest_path))
record <- read_json(record_path, simplifyVector = TRUE)
stopifnot(record$script_sha256 == sha(entry), record$specification_sha256 == SPEC_SHA)
check_inputs()

prepare_old_inputs <- function() {
  object <- readRDS(PATHS[["generator"]])
  stopifnot(identical(as.character(object$feature_names), FEATURES33))
  freeze <- fread(PATHS[["label_freeze"]])
  stopifnot(hash_file(PATHS[["generator"]]) == freeze$generator_file_md5[1],
            hash_file(PATHS[["external_labels"]]) == freeze$label_file_md5[1],
            hash_file(PATHS[["external"]]) == freeze$eicu_input_md5[1],
            !isTRUE(freeze$external_outcomes_read[1]), !isTRUE(freeze$eicu_reclustered[1]),
            hash_file(PATHS[["baseline"]]) == "972e5b37ab179908d926c2cd34618129")
  M <- fread(PATHS[["development"]], select = c(
    "stay_id", "subject_id", "age", "gender", "sofa_score",
    "hospital_mortality", "make30", FEATURES33), showProgress = FALSE)
  bc <- fread(PATHS[["baseline"]], select = c("stay_id", "aki_stage_0_24h"), showProgress = FALSE)
  stopifnot(!anyDuplicated(M$stay_id), !anyDuplicated(bc$stay_id))
  M <- merge(M, bc, by = "stay_id", all.x = TRUE)
  M[, sex := parse_sex(gender)]
  M <- M[complete.cases(M[, c(BASE_NUM, "hospital_mortality", "make30"), with = FALSE])]
  stopifnot(nrow(M) == 20049L)
  source_assignment <- as.data.table(predict_ehr_audit_k2(object, M, id_col = "stay_id"))
  stopifnot(identical(source_assignment$stay_id, M$stay_id))
  M[, cluster_k2 := factor(source_assignment$assigned_cluster, levels = SOURCE_LABEL_LEVELS)]
  Xm <- apply_frozen_recipe(M, object$recipe, FEATURES33)$X
  BASEm <- as.matrix(M[, ..BASE_NUM])
  storage.mode(BASEm) <- "double"
  E <- fread(PATHS[["external"]], showProgress = FALSE)
  labels <- fread(PATHS[["external_labels"]], showProgress = FALSE)
  E <- merge(E, labels[, .(
    patientunitstayid, assigned_cluster, distance_c1, distance_c2,
    signed_margin_c1, absolute_margin, missing_feature_n,
    clipped_feature_n, generator_version, generator_content_hash)],
    by = "patientunitstayid", all.x = FALSE, all.y = TRUE)
  stopifnot(nrow(E) == 17465L, uniqueN(E$uniquepid) == 16212L,
            uniqueN(E$hospitalid) == 199L, !anyDuplicated(E$patientunitstayid))
  E[, age := parse_age(age)]
  E[, sex := parse_sex(sex)]
  setnames(E, "sofa", "sofa_score", skip_absent = TRUE)
  E[, aki_stage_0_24h := as.numeric(aki_stage_0_24h)]
  E[, hosp_mortality := as.integer(hosp_mortality)]
  E[, make_hosp := as.integer(hosp_mortality == 1L |
      (!is.na(scr_discharge) & !is.na(scr_baseline) & scr_discharge >= 1.5 * scr_baseline))]
  E[, cluster_k2 := factor(as.character(assigned_cluster), levels = SOURCE_LABEL_LEVELS)]
  Xe <- apply_frozen_recipe(E, object$recipe, FEATURES33)$X
  BASEe <- as.matrix(E[, ..BASE_NUM])
  storage.mode(BASEe) <- "double"
  list(object = object, M = M, E = E, Xm = Xm, Xe = Xe, BASEm = BASEm, BASEe = BASEe,
       complete = complete.cases(BASEe) & !is.na(E$cluster_k2))
}

if (stage == "recover") {
  stopifnot(!file.exists(file.path(ROOT, "provenance/recovery_verified.json")))
  D <- prepare_old_inputs()
  set.seed(20240601L)
  recovered <- list()
  cv_store <- list()
  prediction_rows <- list()
  summary_rows <- list()
  row_counter <- 0L
  rng_checkpoints <- list()
  for (outcome_name in names(OUTCOMES)) {
    cat("RECOVERY_OUTCOME_BEGIN", outcome_name, now(), "\n")
    rng_checkpoints[[paste0(outcome_name, "_before_cv")]] <- .Random.seed
    spec <- OUTCOMES[[outcome_name]]
    y_m <- as.integer(D$M[[spec$mimic]])
    cv_pred <- internal_cv_predictions(y_m, D$BASEm, D$Xm, D$M$cluster_k2, MODEL_KINDS)
    fits <- lapply(MODEL_KINDS, function(kind) {
      fit_model(kind, y_m, D$BASEm, D$Xm, D$M$cluster_k2)
    })
    names(fits) <- MODEL_KINDS
    recovered[[outcome_name]] <- fits
    cv_store[[outcome_name]] <- cv_pred
    rng_checkpoints[[paste0(outcome_name, "_after_full_fits")]] <- .Random.seed
    sel <- which(D$complete & !is.na(D$E[[spec$eicu]]))
    y_e <- as.integer(D$E[[spec$eicu]][sel])
    stopifnot(length(sel) == if (outcome_name == "mortality") 17284L else 17333L,
              sum(y_e == 1L) == if (outcome_name == "mortality") 2940L else 5714L)
    for (kind in MODEL_KINDS) {
      p_e <- predict_model(kind, fits[[kind]], D$BASEe[sel, , drop = FALSE],
                           D$Xe[sel, , drop = FALSE], D$E$cluster_k2[sel])
      auc_ci <- bootstrap_auc_ci(y_e, p_e, D$E$uniquepid[sel], MODEL_BOOT_B)
      cal <- calibration_slope_intercept(y_e, p_e)
      row_counter <- row_counter + 1L
      prediction_rows[[row_counter]] <- data.table(
        patientunitstayid = D$E$patientunitstayid[sel], uniquepid = D$E$uniquepid[sel],
        hospitalid = D$E$hospitalid[sel], outcome = outcome_name,
        outcome_value = y_e, model = kind, prediction = p_e)
      summary_rows[[row_counter]] <- data.table(
        outcome = outcome_name, model = kind, n_external_stays = length(sel),
        n_external_patients = uniqueN(D$E$uniquepid[sel]), n_external_events = sum(y_e == 1L),
        internal_cv_auc = fast_auc(y_m, cv_pred[, kind]), external_auc = fast_auc(y_e, p_e),
        external_auc_ci_low = auc_ci["low"], external_auc_ci_high = auc_ci["high"],
        external_auc_bootstrap_mcse = auc_ci["mcse"], calibration_slope = cal["slope"],
        calibration_intercept = cal["intercept"], brier_score = mean((p_e - y_e)^2))
      cat("RECOVERY_EXTERNAL_COMPLETE", outcome_name, kind, now(), "\n")
    }
    rng_checkpoints[[paste0(outcome_name, "_after_external_bootstrap")]] <- .Random.seed
    safe_rds(list(fits = fits, cv = cv_pred, rng = .Random.seed),
      file.path(ROOT, "private_local", paste0("recovery_checkpoint_", outcome_name, ".rds")))
  }
  new <- rbindlist(prediction_rows)
  old <- fread(PATHS[["archived_predictions"]])
  stopifnot(!anyDuplicated(new[, .(outcome, model, patientunitstayid)]),
            !anyDuplicated(old[, .(outcome, model, patientunitstayid)]))
  setorder(new, outcome, model, patientunitstayid)
  setorder(old, outcome, model, patientunitstayid)
  stopifnot(nrow(new) == nrow(old),
            identical(new$outcome, old$outcome), identical(new$model, old$model),
            identical(new$patientunitstayid, old$patientunitstayid),
            identical(as.character(new$uniquepid), as.character(old$uniquepid)),
            identical(new$outcome_value, as.integer(old$outcome_value)))
  new[, archived_prediction := old$prediction]
  new[, absolute_error := abs(prediction - archived_prediction)]
  qc <- new[, .(n_stays = .N, max_absolute_error = max(absolute_error),
                mean_absolute_error = mean(absolute_error),
                failures_above_tolerance = sum(!is.finite(absolute_error) | absolute_error > 1e-12),
                tolerance = 1e-12), by = .(outcome, model)]
  safe_csv(qc, file.path(ROOT, "aggregate_outputs/old_prediction_reproduction.csv"))
  safe_csv(new, file.path(ROOT, "private_local/reconstructed_old_predictions_local.csv"))
  safe_rds(recovered, file.path(ROOT, "private_local/recovered_source_fits_local.rds"))
  safe_rds(list(cv = cv_store, rng_checkpoints = rng_checkpoints),
           file.path(ROOT, "private_local/recovered_cv_and_rng_local.rds"))
  safe_csv(rbindlist(summary_rows), file.path(ROOT, "aggregate_outputs/reconstructed_old_performance.csv"))
  print(qc)
  check_inputs()
  if (any(qc$failures_above_tolerance > 0L)) {
    safe_json(list(status = "FAILED_OLD_PREDICTION_REPRODUCTION", time_utc = now(),
      corrected_inputs_evaluated = FALSE), file.path(ROOT, "provenance/recovery_failed.json"))
    stop("Archived predictions not reproduced; corrected evaluation is prohibited.")
  }
  safe_json(list(status = "VERIFIED", time_utc = now(),
    n_prediction_rows = nrow(new), max_absolute_error = max(new$absolute_error),
    tolerance = 1e-12, specification_sha256 = SPEC_SHA,
    fit_file_sha256 = sha(file.path(ROOT, "private_local/recovered_source_fits_local.rds"))),
    file.path(ROOT, "provenance/recovery_verified.json"))
  cat("ALL_OLD_PREDICTIONS_VERIFIED", now(), "\n")
  quit(status = 0L)
}

if (stage == "evaluate") {
  verified_path <- file.path(ROOT, "provenance/recovery_verified.json")
  stopifnot(file.exists(verified_path))
  verified <- read_json(verified_path, simplifyVector = TRUE)
  fit_path <- file.path(ROOT, "private_local/recovered_source_fits_local.rds")
  stopifnot(verified$status == "VERIFIED", verified$max_absolute_error <= 1e-12,
            verified$fit_file_sha256 == sha(fit_path))
  D <- prepare_old_inputs()
  fits <- readRDS(fit_path)
  T <- fread(PATHS[["corrected_temperature"]], showProgress = FALSE)
  stopifnot(nrow(T) == nrow(D$E), !anyDuplicated(T$patientunitstayid))
  ti <- match(D$E$patientunitstayid, T$patientunitstayid)
  stopifnot(!anyNA(ti))
  T <- T[ti]
  stopifnot(identical(D$E$patientunitstayid, T$patientunitstayid))
  temperature_names <- c("temperature_min", "temperature_max")
  corrected_names <- paste0(temperature_names, "_corrected")
  stopifnot(all(corrected_names %in% names(T)))
  corrected <- copy(D$E)
  for (field in temperature_names) corrected[[field]] <- T[[paste0(field, "_corrected")]]
  other_names <- setdiff(names(corrected), temperature_names)
  stopifnot(identical(corrected[, ..other_names], D$E[, ..other_names]))
  Xcorrected <- apply_frozen_recipe(corrected, D$object$recipe, FEATURES33)$X
  stopifnot(identical(D$Xe[, setdiff(FEATURES33, temperature_names), drop = FALSE],
                     Xcorrected[, setdiff(FEATURES33, temperature_names), drop = FALSE]))
  archived <- fread(PATHS[["archived_predictions"]])
  current <- fread(PATHS[["corrected_base_k2_predictions"]])
  prediction_rows <- list()
  performance_rows <- list()
  contrast_rows <- list()
  bootstrap_rows <- list()
  failure_rows <- list()
  for (outcome_name in names(OUTCOMES)) {
    spec <- OUTCOMES[[outcome_name]]
    sel <- which(D$complete & !is.na(D$E[[spec$eicu]]))
    y <- as.integer(D$E[[spec$eicu]][sel])
    ids <- D$E$patientunitstayid[sel]
    patient <- as.character(D$E$uniquepid[sel])
    cur <- current[outcome == outcome_name]
    stopifnot(!anyDuplicated(cur$patientunitstayid), nrow(cur) == length(ids))
    ci <- match(ids, cur$patientunitstayid)
    stopifnot(!anyNA(ci), identical(y, as.integer(cur$outcome_value[ci])),
              identical(patient, as.character(cur$uniquepid[ci])))
    p_current <- cur$prediction[ci]
    pred <- list()
    for (kind in RAW_KINDS) {
      old <- archived[outcome == outcome_name & model == kind]
      oi <- match(ids, old$patientunitstayid)
      stopifnot(!anyNA(oi))
      p_new <- predict_model(kind, fits[[outcome_name]][[kind]], D$BASEe[sel, , drop = FALSE],
                             Xcorrected[sel, , drop = FALSE], D$E$cluster_k2[sel])
      stopifnot(all(is.finite(p_new)), all(p_new >= 0 & p_new <= 1))
      pred[[kind]] <- list(old = old$prediction[oi], corrected = p_new)
      prediction_rows[[paste(outcome_name, kind)]] <- data.table(
        patientunitstayid = ids, uniquepid = patient, outcome = outcome_name,
        outcome_value = y, model = kind, archived_prediction = old$prediction[oi],
        corrected_prediction = p_new, current_base_k2_prediction = p_current)
    }
    seed <- MASTER_SEED + if (outcome_name == "mortality") 5002L else 5006L
    ii_list <- cluster_bootstrap_indices(patient, 1000L, seed)
    for (kind in RAW_KINDS) {
      for (version in c("old", "corrected")) {
        p <- pred[[kind]][[version]]
        auc_vals <- rep(NA_real_, 1000L)
        for (b in seq_len(1000L)) {
          ans <- tryCatch(fast_auc(y[ii_list[[b]]], p[ii_list[[b]]]), error = identity)
          if (inherits(ans, "error") || !is.finite(ans)) {
            failure_rows[[length(failure_rows) + 1L]] <- data.table(
              outcome = outcome_name, model = kind, version = version, replicate = b,
              metric = "auc", reason = if (inherits(ans, "error")) conditionMessage(ans) else "nonfinite_auc")
          } else auc_vals[b] <- ans
        }
        valid <- is.finite(auc_vals)
        stopifnot(any(valid))
        interval <- percentile_interval(auc_vals[valid])
        cal <- calibration_slope_intercept(y, p)
        joint <- glm(y ~ qlogis(clamp_probability(p)), family = binomial())
        performance_rows[[length(performance_rows) + 1L]] <- data.table(
          outcome = outcome_name, model = kind, input_version = version,
          n_external_stays = length(ids), n_external_patients = uniqueN(patient),
          n_external_events = sum(y == 1L), external_auc = fast_auc(y, p),
          external_auc_ci_low = interval[1], external_auc_ci_high = interval[2],
          external_auc_bootstrap_mcse = sd(auc_vals[valid]) / sqrt(sum(valid)),
          calibration_slope = unname(cal["slope"]),
          calibration_intercept = unname(cal["intercept"]),
          calibration_joint_intercept = unname(coef(joint)[1]),
          brier_score = mean((p - y)^2), bootstrap_requested = 1000L,
          bootstrap_valid = sum(valid), bootstrap_failed = sum(!valid), bootstrap_seed = seed,
          intercept_definition = "offset-only CITL; archived Figure 4 convention")
        bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.table(
          outcome = outcome_name, model = kind, input_version = version,
          metric = "auc", replicate = seq_len(1000L), value = auc_vals, bootstrap_seed = seed)
      }
      paired_vals <- vapply(ii_list, function(ii) {
        fast_auc(y[ii], pred[[kind]]$corrected[ii]) - fast_auc(y[ii], p_current[ii])
      }, numeric(1))
      valid <- is.finite(paired_vals)
      if (any(!valid)) {
        for (b in which(!valid)) failure_rows[[length(failure_rows) + 1L]] <- data.table(
          outcome = outcome_name, model = kind, version = "corrected", replicate = b,
          metric = "paired_delta_auc", reason = "nonfinite_paired_auc")
      }
      stopifnot(any(valid))
      interval <- percentile_interval(paired_vals[valid])
      contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
        outcome = outcome_name, model_a = kind, model_b = "base_l1_k2_corrected",
        contrast = paste0(kind, "_minus_corrected_base_l1_k2"),
        delta_auc = fast_auc(y, pred[[kind]]$corrected) - fast_auc(y, p_current),
        ci_low = interval[1], ci_high = interval[2], bootstrap_mean = mean(paired_vals[valid]),
        bootstrap_mcse = sd(paired_vals[valid]) / sqrt(sum(valid)), bootstrap_requested = 1000L,
        bootstrap_valid = sum(valid), bootstrap_failed = sum(!valid), bootstrap_seed = seed,
        resampling_unit = "uniquepid", uncertainty = "fixed source fits and predictions")
      bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.table(
        outcome = outcome_name, model = kind, input_version = "corrected",
        metric = "paired_delta_auc", replicate = seq_len(1000L), value = paired_vals, bootstrap_seed = seed)
    }
    cat("CORRECTED_OUTCOME_COMPLETE", outcome_name, now(), "\n")
  }
  perf <- rbindlist(performance_rows)
  pair <- rbindlist(contrast_rows)
  failures <- if (length(failure_rows)) rbindlist(failure_rows) else data.table(
    outcome = character(), model = character(), version = character(), replicate = integer(),
    metric = character(), reason = character())
  safe_csv(rbindlist(prediction_rows), file.path(ROOT, "private_local/corrected_raw_predictions_local.csv"))
  safe_csv(perf, file.path(ROOT, "aggregate_outputs/raw_old_corrected_performance.csv"))
  safe_csv(pair, file.path(ROOT, "aggregate_outputs/raw_vs_corrected_base_k2_paired_auc.csv"))
  safe_csv(rbindlist(bootstrap_rows), file.path(ROOT, "aggregate_outputs/patient_bootstrap_replicate_metrics.csv"))
  safe_csv(failures, file.path(ROOT, "aggregate_outputs/bootstrap_failures.csv"))
  oldp <- perf[input_version == "old"]
  newp <- perf[input_version == "corrected"]
  comparison <- merge(oldp, newp, by = c("outcome", "model"), suffixes = c("_old", "_corrected"))
  for (metric in c("external_auc", "brier_score", "calibration_slope", "calibration_intercept")) {
    comparison[[paste0(metric, "_change")]] <- comparison[[paste0(metric, "_corrected")]] - comparison[[paste0(metric, "_old")]]
  }
  safe_csv(comparison, file.path(ROOT, "aggregate_outputs/raw_old_corrected_comparison.csv"))
  figure <- fread(PATHS[["figure_transport"]])
  stopifnot(nrow(figure) == 8L)
  for (i in seq_len(nrow(newp))) {
    z <- newp[i]
    j <- which(figure$outcome == z$outcome & figure$model == z$model)
    stopifnot(length(j) == 1L, figure$n_external_stays[j] == z$n_external_stays,
              figure$n_external_patients[j] == z$n_external_patients,
              figure$n_external_events[j] == z$n_external_events)
    for (field in c("external_auc", "external_auc_ci_low", "external_auc_ci_high",
      "external_auc_bootstrap_mcse", "calibration_slope", "calibration_intercept", "brier_score")) {
      figure[[field]][j] <- z[[field]]
    }
    figure$transport_drop_internal_minus_external[j] <- figure$internal_cv_auc[j] - z$external_auc
    figure$emphasis[j] <- "Corrected-temperature external evaluation; recovered frozen Raw model"
  }
  stopifnot(identical(figure[!model %in% RAW_KINDS], fread(PATHS[["figure_transport"]])[!model %in% RAW_KINDS]))
  safe_csv(figure, file.path(ROOT, "figure_source/Figure4_panelC_L1_transport.csv"))
  calibration <- fread(PATHS[["figure_calibration"]])
  stopifnot(nrow(calibration) == 4L, all(calibration$model == "base_gen"))
  safe_csv(calibration, file.path(ROOT, "figure_source/Figure4_panelC_L1_calibration.csv"))
  check_inputs()
  safe_json(list(status = "COMPLETE", time_utc = now(),
    old_predictions_verified = TRUE, bootstrap_failure_rows = nrow(failures),
    patient_bootstrap_requested_per_statistic = 1000L,
    corrected_fields_only = temperature_names,
    raw_recipe = "unchanged archived L2 object$recipe",
    current_base_k2_refitted = FALSE, current_hospital_recalibration_refitted = FALSE,
    original_inputs_unchanged = TRUE), file.path(ROOT, "provenance/completion_record.json"))
  writeLines(capture.output(sessionInfo()), file.path(ROOT, "logs/evaluation_sessionInfo.txt"))
  print(perf)
  print(pair)
  cat("RAW_TEMPERATURE_COMPLETION_FINISHED", now(), "\n")
}
