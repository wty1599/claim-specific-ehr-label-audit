options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

D4_PROJECT_ROOT <- normalizePath(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)
D4_RESTRICTED_DATA_ROOT <- Sys.getenv(
  "EHR_AUDIT_RESTRICTED_DATA_ROOT",
  unset = file.path(D4_PROJECT_ROOT, "data_restricted")
)
D4_PACKAGE_ROOT <- file.path(D4_PROJECT_ROOT, "R", "domain4")
D4_SCRIPT_DIR <- file.path(D4_PACKAGE_ROOT, "landmark")
D4_SPEC <- file.path(
  D4_PACKAGE_ROOT, "spec", "D4_LANDMARK_SUPPORT_PRERUN_SPEC_FROZEN_20260730.md"
)

D4_RUN_MODE <- Sys.getenv("D4_RUN_MODE", unset = "smoke")
if (!D4_RUN_MODE %chin% c("smoke", "formal")) {
  stop("D4_RUN_MODE must be smoke or formal.", call. = FALSE)
}

default_output <- if (D4_RUN_MODE == "formal") {
  file.path(D4_PACKAGE_ROOT, "output_formal_v1")
} else {
  file.path(
    D4_PACKAGE_ROOT,
    paste0("output_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
}
D4_OUTPUT_ROOT <- Sys.getenv("D4_OUTPUT_ROOT", unset = default_output)
D4_TABLE_DIR <- file.path(D4_OUTPUT_ROOT, "tables")
D4_FIGURE_DIR <- file.path(D4_OUTPUT_ROOT, "figures")
D4_LOG_DIR <- file.path(D4_OUTPUT_ROOT, "logs")
D4_PROVENANCE_DIR <- file.path(D4_OUTPUT_ROOT, "provenance")
D4_PRIVATE_DIR <- file.path(D4_OUTPUT_ROOT, "private_not_for_release")

D4_INPUT <- c(
  labels = file.path(D4_RESTRICTED_DATA_ROOT, "model", "labels_primary_mice.rds"),
  final_full = file.path(D4_RESTRICTED_DATA_ROOT, "final_full.csv"),
  baseline = file.path(D4_RESTRICTED_DATA_ROOT, "baseline_covars.csv"),
  rrt_first = file.path(D4_RESTRICTED_DATA_ROOT, "rrt_first.csv"),
  rrt_base_times = file.path(D4_RESTRICTED_DATA_ROOT, "rrt_base_times.csv"),
  rrt_sql = file.path(D4_PROJECT_ROOT, "sql", "domain4", "01_extract_rrt_timing_inputs.sql"),
  frozen_spec = D4_SPEC
)

D4_EXPECTED_SHA256 <- c(
  labels = "159576D0B2796E9786DB0EB4C8747AF987BE5E53AF200419E3E37567C1D875F9",
  final_full = "CE0B38D2D3D6263A5C4615A6AEA1E78B19840BD4B35046DBFA65EE6BFBBA633B",
  baseline = "3604A7AD95B0F3E2AF4CAA632C753C11E0C56CF0CC3BB5B6AB8DE984D85DB107",
  rrt_first = "5D50E32F9874F44F0D67E4D06331236162B056642B50C8ED16514429F7102F2F",
  rrt_base_times = "0E7AD66D2D7D9BF25209AAA9047445E1D2EB374B99D8B7B7C5B520A1BBB41C58",
  rrt_sql = "8926BDF1C3362D786AEAE93FECBBE579B9B1A3724FFADE0CA7DE8EBCB7776B3A"
)

D4_LANDMARK_H <- 24
D4_WINDOW_END_H <- 72
D4_TIE_TOL <- 1e-8
D4_PRIMARY_GATE <- "B"
D4_GATES <- list(
  A = list(frac = 0.005, events = 5L),
  B = list(frac = 0.010, events = 10L),
  C = list(frac = 0.020, events = 10L),
  D = list(frac = 0.010, events = 20L)
)
D4_SUPPORT_THRESHOLDS <- list(
  outside_overlap_max = 0.10,
  ess_ratio_min = 0.25,
  minimum_cell_min = 30L,
  extreme_weight_cut = 10,
  extreme_weight_prop_max = 0.01
)
D4_PS_COVARIATES <- c(
  "cluster_k2", "age", "gender", "sofa_score", "aki_stage_0_24h",
  "sapsii", "creatinine_max", "bun_max", "potassium_max",
  "bicarbonate_min", "ph_min", "lactate_max", "urine_output_24h_ml",
  "mbp_min", "heart_rate_max", "resp_rate_max", "gcs_min"
)
D4_FORBIDDEN_PS_TERMS <- c(
  "mortality_30d", "death_time_hours", "first_rrt_hours",
  "rrt_first_24_72", "first_event_state", "make30", "rrt_30d",
  "persistent_rd_30d", "terminal_scr"
)

d4_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n",
      sep = "")
}

d4_sha256 <- function(path) {
  toupper(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

d4_require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing D4 files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
}

d4_require_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(label, " lacks columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

d4_atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  fwrite(as.data.table(x), tmp, na = "")
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Atomic CSV write failed: ", path, call. = FALSE)
  }
  invisible(path)
}

d4_atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, tmp, compress = "xz")
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Atomic RDS write failed: ", path, call. = FALSE)
  }
  invisible(path)
}

d4_atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    stop("Atomic text write failed: ", path, call. = FALSE)
  }
  invisible(path)
}

d4_stage_marker <- function(stage) {
  file.path(D4_LOG_DIR, paste0(stage, ".ok"))
}

d4_require_stage <- function(stage) {
  p <- d4_stage_marker(stage)
  if (!file.exists(p)) stop("Required prior stage missing: ", p, call. = FALSE)
}

d4_mark_stage <- function(stage) {
  d4_atomic_write_lines(
    c(
      paste0("stage=", stage),
      paste0("run_mode=", D4_RUN_MODE),
      paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    ),
    d4_stage_marker(stage)
  )
}

d4_binary_auc <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- score[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

d4_impute_data <- function(data, covariates, mode = c("median_mode", "complete_case")) {
  mode <- match.arg(mode)
  d <- copy(data)
  audit <- rbindlist(lapply(covariates, function(v) {
    x <- d[[v]]
    data.table(
      variable = v,
      class = paste(class(x), collapse = "/"),
      missing_n = sum(is.na(x) | (is.character(x) & !nzchar(x))),
      missing_prop = mean(is.na(x) | (is.character(x) & !nzchar(x)))
    )
  }))
  if (mode == "complete_case") {
    keep <- complete.cases(d[, ..covariates])
    audit[, `:=`(imputation_method = "none_complete_case", fill_value = NA_character_)]
    return(list(data = d[keep], audit = audit, n_removed = sum(!keep)))
  }
  fill <- character(length(covariates))
  for (i in seq_along(covariates)) {
    v <- covariates[i]
    x <- d[[v]]
    if (is.numeric(x) || is.integer(x)) {
      value <- suppressWarnings(stats::median(x, na.rm = TRUE))
      if (!is.finite(value)) stop("Cannot impute numeric covariate: ", v)
      x[is.na(x)] <- value
      fill[i] <- format(value, digits = 16, scientific = FALSE, trim = TRUE)
    } else {
      x <- as.character(x)
      tab <- sort(table(x[!is.na(x) & nzchar(x)]), decreasing = TRUE)
      if (!length(tab)) stop("Cannot impute categorical covariate: ", v)
      value <- names(tab)[1L]
      x[is.na(x) | !nzchar(x)] <- value
      d[[v]] <- factor(x)
      fill[i] <- value
      next
    }
    d[[v]] <- x
  }
  audit[, `:=`(imputation_method = "median_mode", fill_value = fill)]
  list(data = d, audit = audit, n_removed = 0L)
}

d4_fit_support <- function(data, imputation_mode = "median_mode") {
  d4_require_columns(
    data,
    c(D4_PS_COVARIATES, "rrt_first_24_72", "mortality_30d",
      "first_event_state", "stay_id"),
    "D4 risk-set data"
  )
  if (any(D4_FORBIDDEN_PS_TERMS %chin% D4_PS_COVARIATES)) {
    stop("Forbidden term present in the propensity covariates.", call. = FALSE)
  }
  imp <- d4_impute_data(data, D4_PS_COVARIATES, imputation_mode)
  dm <- as.data.table(imp$data)
  set(dm, j = "cluster_k2", value = factor(dm$cluster_k2, levels = c(1L, 2L)))
  set(dm, j = "gender", value = factor(dm$gender))
  if (nrow(dm) < 1L || uniqueN(dm$rrt_first_24_72) < 2L) {
    stop("Allocation model lacks both treatment levels.", call. = FALSE)
  }
  form <- reformulate(D4_PS_COVARIATES, response = "rrt_first_24_72")
  warning_text <- character()
  fit <- withCallingHandlers(
    glm(form, data = dm, family = binomial()),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (!isTRUE(fit$converged)) {
    stop("Primary unpenalized propensity model did not converge.", call. = FALSE)
  }
  ps_raw <- as.numeric(predict(fit, type = "response"))
  eps <- 1e-6
  set(dm, j = "ps", value = pmin(pmax(ps_raw, eps), 1 - eps))
  set(
    dm,
    j = "ipw_ate",
    value = fifelse(dm$rrt_first_24_72 == 1L, 1 / dm$ps, 1 / (1 - dm$ps))
  )

  summarize_group <- function(z, label) {
    treated <- z[rrt_first_24_72 == 1L]
    control <- z[rrt_first_24_72 == 0L]
    if (!nrow(treated) || !nrow(control)) {
      return(data.table(
        phenotype = label, n = nrow(z), n_treated = nrow(treated),
        n_control = nrow(control), treatment_prevalence = mean(z$rrt_first_24_72),
        common_support_low = NA_real_, common_support_high = NA_real_,
        outside_overlap_proportion = NA_real_, ate_ipw_ess = NA_real_,
        ess_ratio = NA_real_, minimum_cell = min(nrow(treated), nrow(control)),
        extreme_weight_proportion = mean(z$ipw_ate > D4_SUPPORT_THRESHOLDS$extreme_weight_cut),
        ps_auc = d4_binary_auc(z$rrt_first_24_72, z$ps)
      ))
    }
    support_low <- max(min(treated$ps), min(control$ps))
    support_high <- min(max(treated$ps), max(control$ps))
    z[, outside_overlap := ps < support_low | ps > support_high]
    ess <- sum(z$ipw_ate)^2 / sum(z$ipw_ate^2)
    data.table(
      phenotype = label,
      n = nrow(z),
      n_treated = nrow(treated),
      n_control = nrow(control),
      treatment_prevalence = mean(z$rrt_first_24_72),
      common_support_low = support_low,
      common_support_high = support_high,
      outside_overlap_proportion = mean(z$outside_overlap),
      ate_ipw_ess = ess,
      ess_ratio = ess / nrow(z),
      minimum_cell = min(nrow(treated), nrow(control)),
      extreme_weight_proportion =
        mean(z$ipw_ate > D4_SUPPORT_THRESHOLDS$extreme_weight_cut),
      ps_auc = d4_binary_auc(z$rrt_first_24_72, z$ps)
    )
  }

  overall <- summarize_group(copy(dm), "Overall")
  by_k2 <- rbindlist(lapply(c("1", "2"), function(k) {
    summarize_group(
      copy(dm[as.character(cluster_k2) == k]),
      if (k == "1") "C1 higher-risk" else "C2 lower-risk"
    )
  }))
  support <- rbindlist(list(overall, by_k2), fill = TRUE)

  for (phen in c("C1 higher-risk", "C2 lower-risk")) {
    row <- support[phenotype == phen]
    idx <- if (phen == "C1 higher-risk") {
      as.character(dm$cluster_k2) == "1"
    } else {
      as.character(dm$cluster_k2) == "2"
    }
    dm[idx, outside_overlap_within_k2 :=
         ps < row$common_support_low | ps > row$common_support_high]
  }

  cells <- dm[, .(
    n = .N,
    mortality_30d_events = sum(mortality_30d == 1L, na.rm = TRUE),
    mortality_30d_missing = sum(is.na(mortality_30d))
  ), by = .(
    cluster_k2 = as.character(cluster_k2),
    rrt_first_24_72
  )]
  cells[, phenotype := fifelse(
    cluster_k2 == "1", "C1 higher-risk", "C2 lower-risk"
  )]
  setcolorder(
    cells,
    c("cluster_k2", "phenotype", "rrt_first_24_72", "n",
      "mortality_30d_events", "mortality_30d_missing")
  )

  model_qc <- data.table(
    imputation_mode = imputation_mode,
    n_input = nrow(data),
    n_model = nrow(dm),
    n_removed = imp$n_removed,
    n_treated = sum(dm$rrt_first_24_72),
    formula = paste(deparse(form), collapse = ""),
    cluster_k2_in_formula = "cluster_k2" %chin% all.vars(form),
    forbidden_term_in_formula =
      any(D4_FORBIDDEN_PS_TERMS %chin% attr(terms(form), "term.labels")),
    converged = isTRUE(fit$converged),
    warning_count = length(warning_text),
    warning_text = paste(unique(warning_text), collapse = " | "),
    coefficient_count = length(coef(fit)),
    nonfinite_coefficient_count = sum(!is.finite(coef(fit))),
    fitted_probability_min = min(ps_raw),
    fitted_probability_max = max(ps_raw),
    clipped_probability_n = sum(ps_raw < eps | ps_raw > 1 - eps),
    pooled_ps_auc = d4_binary_auc(dm$rrt_first_24_72, dm$ps),
    effect_model_fitted = FALSE
  )

  list(
    data = dm,
    support = support,
    cells = cells,
    imputation = imp$audit,
    model_qc = model_qc
  )
}

d4_apply_state <- function(support, cells) {
  groups <- c("C1 higher-risk", "C2 lower-risk")
  s <- support[phenotype %chin% groups]
  if (nrow(s) != 2L || anyNA(s[, .(
    outside_overlap_proportion, ess_ratio, minimum_cell,
    extreme_weight_proportion
  )])) {
    return(list(
      state = "INTERFACE_INCOMPLETE",
      triggers = "missing phenotype-specific support metric",
      gate_table = data.table(),
      support_alerts = data.table()
    ))
  }
  gate_rows <- rbindlist(lapply(names(D4_GATES), function(g) {
    rule <- D4_GATES[[g]]
    rbindlist(lapply(c("1", "2"), function(k) {
      a <- cells[cluster_k2 == k & rrt_first_24_72 == 1L]
      b <- cells[cluster_k2 == k & rrt_first_24_72 == 0L]
      n_t <- if (nrow(a)) a$n else 0L
      n_c <- if (nrow(b)) b$n else 0L
      ev_t <- if (nrow(a)) a$mortality_30d_events else 0L
      ev_c <- if (nrow(b)) b$mortality_30d_events else 0L
      frac <- n_t / (n_t + n_c)
      data.table(
        gate = g,
        phenotype = if (k == "1") "C1 higher-risk" else "C2 lower-risk",
        treated_fraction_threshold = rule$frac,
        event_threshold = rule$events,
        n_treated = n_t,
        n_control = n_c,
        treated_fraction = frac,
        events_treated = ev_t,
        events_control = ev_c,
        pass = frac >= rule$frac &&
          ev_t >= rule$events &&
          ev_c >= rule$events
      )
    }))
  }))
  support_alerts <- s[, .(
    phenotype,
    outside_overlap_proportion,
    outside_overlap_alert =
      outside_overlap_proportion >
        D4_SUPPORT_THRESHOLDS$outside_overlap_max,
    ess_ratio,
    ess_ratio_alert = ess_ratio < D4_SUPPORT_THRESHOLDS$ess_ratio_min,
    minimum_cell,
    minimum_cell_alert =
      minimum_cell < D4_SUPPORT_THRESHOLDS$minimum_cell_min,
    extreme_weight_proportion,
    extreme_weight_alert =
      extreme_weight_proportion >
        D4_SUPPORT_THRESHOLDS$extreme_weight_prop_max
  )]
  gate_b_fail <- gate_rows[gate == D4_PRIMARY_GATE & !pass]
  alert_cols <- c(
    "outside_overlap_alert", "ess_ratio_alert",
    "minimum_cell_alert", "extreme_weight_alert"
  )
  support_fail <- support_alerts[
    rowSums(as.matrix(support_alerts[, ..alert_cols])) > 0L
  ]
  if (nrow(gate_b_fail) || nrow(support_fail)) {
    pieces <- c(
      if (nrow(gate_b_fail)) paste0(
        "Gate B failed: ",
        paste(gate_b_fail$phenotype, collapse = "; ")
      ),
      if (nrow(support_fail)) paste0(
        "support alert: ",
        paste(support_fail$phenotype, collapse = "; ")
      )
    )
    state <- "SUPPORT_INADEQUATE"
    triggers <- paste(pieces, collapse = " | ")
  } else {
    state <- "SUPPORT_ADEQUATE_FOR_MODEL_DEVELOPMENT"
    triggers <- "none"
  }
  list(
    state = state,
    triggers = triggers,
    gate_table = gate_rows,
    support_alerts = support_alerts
  )
}

d4_init_dirs <- function() {
  for (d in c(
    D4_OUTPUT_ROOT, D4_TABLE_DIR, D4_FIGURE_DIR, D4_LOG_DIR,
    D4_PROVENANCE_DIR, D4_PRIVATE_DIR
  )) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

d4_require_files(D4_INPUT)
spec_text <- readLines(D4_SPEC, warn = FALSE, encoding = "UTF-8")
if (!any(spec_text == "Status: **FROZEN - AUTHORIZED FOR FORMAL EXECUTION**")) {
  stop("D4 frozen specification lacks the authorization marker.", call. = FALSE)
}
d4_init_dirs()
