################################################################################
# 05_missingness_mechanism_intensity_FINAL.R
#
# Complete standalone rerun of Sensitivity Analysis 05.
#
# Locked design:
#   S1 pure severity continuum; S2 true K=3 subtypes; S5 Domain 2 positive control
#   MCAR / MAR / MNAR x target missingness 0 / 0.15 / 0.30
#   formal 50 repeats per setting = 1350 tasks
#
# Main additions:
#   realized-missingness fidelity QC; variable-level rates; mechanism-signal QC;
#   paired complete-data contrasts; MICE completion/stability/model variability;
#   Domain 1 and Domain 2 operating characteristics; MC-SE; ADEMP.
#
# Deliberate exclusions: PAC, Gap, Domain 3, Domain 4, GMM, Ward, external S5.
################################################################################

options(stringsAsFactors = FALSE, warn = 1)

## ---- user settings -----------------------------------------------------------
PROJECT_DIR <- file.path(Sys.getenv("EHR_AUDIT_WORK_ROOT", unset = ""), "analysis_archive", "simulations", "reproduction_runs")
ANALYSIS_ID <- "sensitivity_05_missingness_mechanism_intensity_FINAL"
SIM_ROOT <- file.path(PROJECT_DIR, ANALYSIS_ID)

USE_TEST_MODE <- FALSE
N_REP_FORMAL <- 50L
N_REP_TEST <- 1L
N_ACTIVE_REP <- if (USE_TEST_MODE) N_REP_TEST else N_REP_FORMAL

SCENARIOS <- c(
  "S1_pure_severity_continuum",
  "S2_true_discrete_subtypes",
  "S5_domain2_positive_control"
)
MECHANISMS <- c("MCAR", "MAR", "MNAR")
TARGET_GRID <- c(0, 0.15, 0.30)

N_TOTAL <- if (USE_TEST_MODE) 900L else 5000L
TRAIN_PROP <- 0.70
MICE_M <- if (USE_TEST_MODE) 2L else 5L
MICE_MAXIT <- if (USE_TEST_MODE) 2L else 5L
PRIMARY_IMPUTATION_ID <- 1L
N_WORKERS <- if (USE_TEST_MODE) 2L else 14L
GLOBAL_SEED <- 2026071205L

WINSOR_PROBS <- c(0.01, 0.99)
GLMNET_ALPHA <- 0.5
GLMNET_NFOLD <- 5L
SIGCLUST_B <- if (USE_TEST_MODE) 9L else 49L
SIGCLUST_N <- if (USE_TEST_MODE) 300L else 1000L
DIP_N <- if (USE_TEST_MODE) 400L else 1500L

DIP_ALPHA <- 0.05
SIGCLUST_EFFECT_CUTOFF <- 3.0
TRUE_K3_ARI_CUTOFF <- 0.80
MICE_ARI_CUTOFF <- 0.90
D2_PRACTICAL_NULL <- 0.001

OVERALL_ABS_ERROR_TOL <- 0.015
OVERALL_Q975_ABS_ERROR_TOL <- 0.030
MECHANISM_RANGE_TOL <- 0.020
MECHANISM_Q975_RANGE_TOL <- 0.040
MAR_SIGNAL_AUC_MIN <- 0.58
MNAR_SIGNAL_AUC_MIN <- 0.58
MCAR_ABS_COR_MAX <- 0.05
ZERO_EQUIVALENCE_TOL <- 1e-10

S5_ORACLE_ACCURACY <- 0.90
S5_FEATURE_SIGNAL <- 0.20
S5_MORTALITY_LOGOR <- 1.05
S5_MAKE_LOGOR <- 0.85

RESUME_COMPLETED_TASKS <- TRUE
ALLOW_PARTIAL_OUTPUTS <- TRUE
STOP_ON_CRITICAL_QC <- TRUE

RUN_MODE <- if (USE_TEST_MODE) "test" else "formal"
OUT_DIR <- file.path(SIM_ROOT, paste0("output_", RUN_MODE))
TAB_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
QC_DIR <- file.path(OUT_DIR, "qc")
LOG_DIR <- file.path(OUT_DIR, "logs")
CKPT_DIR <- file.path(OUT_DIR, "task_checkpoints")
for (d in c(SIM_ROOT, OUT_DIR, TAB_DIR, FIG_DIR, QC_DIR, LOG_DIR, CKPT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
LOG_FILE <- file.path(LOG_DIR, paste0("sensitivity05_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
CONFIG_SIGNATURE <- paste(
  ANALYSIS_ID, RUN_MODE, paste(SCENARIOS, collapse = "|"),
  paste(MECHANISMS, collapse = "|"), paste(TARGET_GRID, collapse = "|"),
  N_ACTIVE_REP, N_TOTAL, TRAIN_PROP, MICE_M, MICE_MAXIT,
  SIGCLUST_B, S5_ORACLE_ACCURACY, S5_FEATURE_SIGNAL,
  S5_MORTALITY_LOGOR, S5_MAKE_LOGOR, GLOBAL_SEED, sep = "__"
)

## ---- packages ----------------------------------------------------------------
required_pkgs <- c(
  "data.table", "mice", "mclust", "diptest", "cluster", "glmnet",
  "pROC", "future", "future.apply", "ggplot2"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(data.table); library(mice); library(mclust); library(diptest)
  library(cluster); library(glmnet); library(pROC)
  library(future); library(future.apply); library(ggplot2)
})
options(future.globals.maxSize = 12 * 1024^3)
data.table::setDTthreads(1L)

## ---- generic utilities -------------------------------------------------------
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  x <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  message(x); cat(x, "\n", file = LOG_FILE, append = TRUE); invisible(x)
}
safe_mean <- function(x) { x <- as.numeric(x); if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) }
safe_sd <- function(x) { x <- as.numeric(x); x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else sd(x) }
safe_q <- function(x, p) { x <- as.numeric(x); x <- x[is.finite(x)]; if (!length(x)) NA_real_ else as.numeric(quantile(x, p, names = FALSE, type = 7)) }
safe_num <- function(x, default = NA_real_) if (!length(x) || all(is.na(x))) default else as.numeric(x[1L])
safe_rbind <- function(x) {
  x <- Filter(function(z) !is.null(z) && (is.data.frame(z) || is.data.table(z)) && nrow(z), x)
  if (!length(x)) return(data.table())
  rbindlist(lapply(x, as.data.table), fill = TRUE, use.names = TRUE)
}
summary_numeric <- function(dt, by, value, metric) {
  if (!nrow(dt) || !value %in% names(dt)) return(data.table())
  dt[, { z <- as.numeric(get(value)); z <- z[is.finite(z)]; n <- length(z)
    list(metric = metric, nonmissing_repeats = n,
         mean = if (n) mean(z) else NA_real_, empirical_sd = if (n > 1L) sd(z) else NA_real_,
         mcse_mean = if (n > 1L) sd(z)/sqrt(n) else NA_real_,
         median = if (n) median(z) else NA_real_, q025 = if (n) safe_q(z, .025) else NA_real_,
         q25 = if (n) safe_q(z, .25) else NA_real_, q75 = if (n) safe_q(z, .75) else NA_real_,
         q975 = if (n) safe_q(z, .975) else NA_real_)
  }, by = by]
}
summary_binary <- function(dt, by, value, metric, denominator) {
  if (!nrow(dt) || !value %in% names(dt)) return(data.table())
  dt[, { z <- as.numeric(get(value)); num <- sum(z == 1, na.rm = TRUE); obs <- sum(is.finite(z)); den <- as.integer(denominator)
    p <- num/den; ci <- binom.test(num, den)$conf.int
    list(metric = metric, numerator = num, denominator = den, observed_nonmissing = obs,
         missing_metric_repeats = den - obs, proportion = p,
         mcse = sqrt(p*(1-p)/den), exact_low = ci[1], exact_high = ci[2])
  }, by = by]
}
stdz <- function(x) { x <- as.numeric(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s < 1e-8) rep(0, length(x)) else (x - mean(x, na.rm = TRUE))/s }
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 5L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}
rank_auc <- function(y, score) {
  ok <- is.finite(y) & is.finite(score); y <- as.integer(y[ok]); score <- score[ok]
  if (!length(y) || length(unique(y)) < 2L) return(NA_real_)
  n1 <- sum(y == 1L); n0 <- sum(y == 0L); r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1*(n1+1)/2)/(n1*n0)
}
adjusted_rand <- function(a, b) if (length(a) != length(b) || anyNA(a) || anyNA(b)) NA_real_ else as.numeric(mclust::adjustedRandIndex(a, b))
write_session_info <- function(path) {
  info <- capture.output(utils::sessionInfo())
  writeLines(
    text = info,
    con = path,
    useBytes = TRUE
  )
  invisible(path)
}
theme_publication <- function() theme_classic(base_size = 11) + theme(
  legend.position = "bottom", strip.background = element_blank(),
  strip.text = element_text(face = "bold"), plot.title.position = "plot"
)
save_plot_both <- function(p, stem, width = 8, height = 5) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 300)
}

## ---- feature metadata and DGM ------------------------------------------------
feature_names <- c(
  "bilirubin_total_max", "alt_max", "pao2fio2ratio_min", "abs_lymphocytes_min",
  "lactate_max", "ph_min", "pco2_max", "calcium_min", "calcium_max", "ptt_max",
  "inr_max", "temperature_min", "temperature_max", "urine_output_24h_ml",
  "glucose_max", "aniongap_max", "potassium_min", "potassium_max", "hemoglobin_min",
  "sodium_min", "sodium_max", "wbc_max", "platelets_min", "bicarbonate_min",
  "chloride_min", "chloride_max", "bun_max", "creatinine_max", "resp_rate_max",
  "gcs_min", "spo2_min", "heart_rate_max", "mbp_min"
)
baseline_vars <- c("age", "sex", "sofa", "aki_stage")
outcome_vars <- c("make30", "mortality_30d")
imputation_vars <- c(feature_names, baseline_vars)
severity_loading <- c(
  .45,.25,-.65,-.35,.85,-.75,.30,-.20,.10,.40,.45,-.15,.20,-.95,.25,.80,
  -.10,.35,-.25,-.10,.10,.35,-.40,-.85,-.10,.15,.90,.95,.50,-.60,-.35,.55,-.70
)
names(severity_loading) <- feature_names

make_patterns <- function() {
  out <- matrix(0, 3, length(feature_names), dimnames = list(paste0("Subtype",1:3), feature_names))
  out[1, c("lactate_max","aniongap_max","heart_rate_max","resp_rate_max","wbc_max","temperature_max")] <- 1
  out[1, c("mbp_min","bicarbonate_min","ph_min","platelets_min")] <- -1
  out[2, c("creatinine_max","bun_max","potassium_max","inr_max","ptt_max")] <- 1
  out[2, c("urine_output_24h_ml","calcium_min","hemoglobin_min")] <- -1
  out[3, c("bilirubin_total_max","alt_max","pco2_max","chloride_max")] <- 1
  out[3, c("pao2fio2ratio_min","spo2_min","gcs_min","abs_lymphocytes_min")] <- -1
  out <- out - matrix(colMeans(out), 3, ncol(out), byrow = TRUE)
  out/sqrt(rowMeans(out^2))
}
subtype_patterns <- make_patterns()
set.seed(91001L); factor_loadings <- matrix(rnorm(length(feature_names)*5, sd = .25), nrow = length(feature_names))
set.seed(91002L); missing_weights <- sample(seq(.60,1.40,length.out = length(feature_names))); missing_weights <- missing_weights/mean(missing_weights); names(missing_weights) <- feature_names

scenario_pars <- function(scenario) switch(scenario,
  S1_pure_severity_continuum = list(delta = 0, severity_mult = 1, noise = .75, hidden_signal = 0),
  S2_true_discrete_subtypes = list(delta = 3.2, severity_mult = .35, noise = .45, hidden_signal = 0),
  S5_domain2_positive_control = list(delta = 0, severity_mult = 1, noise = .75, hidden_signal = S5_FEATURE_SIGNAL),
  stop("Unknown scenario: ", scenario, call. = FALSE)
)

generate_complete <- function(scenario, repeat_id, seed) {
  p <- scenario_pars(scenario); set.seed(seed)
  subtype <- sample.int(3L, N_TOTAL, replace = TRUE); severity <- rnorm(N_TOTAL); hidden <- rbinom(N_TOTAL,1,.5)
  F <- matrix(rnorm(N_TOTAL*ncol(factor_loadings)), nrow = N_TOTAL)
  noise <- matrix(rnorm(N_TOTAL*length(feature_names), sd = p$noise), nrow = N_TOTAL)
  hidden_pattern <- setNames(rep(0,length(feature_names)), feature_names)
  hidden_pattern[c("lactate_max","bicarbonate_min","bun_max","gcs_min","mbp_min")] <- c(1,-1,1,-1,-1)
  X <- p$severity_mult*(severity %o% severity_loading) + F %*% t(factor_loadings) + noise +
       p$delta*subtype_patterns[subtype,,drop=FALSE] + p$hidden_signal*((hidden-.5) %o% hidden_pattern)
  age <- pmin(95,pmax(18,round(65+11*rnorm(N_TOTAL)+2*severity)))
  sex <- rbinom(N_TOTAL,1,plogis(-.05+.08*severity))
  sofa <- pmin(24,pmax(0,round(6+2.8*severity+rnorm(N_TOTAL,sd=2))))
  aki <- pmin(3L,pmax(0L,as.integer(cut(severity+rnorm(N_TOTAL,sd=.8), c(-Inf,-.6,.2,1,Inf), labels=FALSE))-1L))
  h <- if (scenario == "S5_domain2_positive_control") hidden else 0
  mort <- rbinom(N_TOTAL,1,plogis(-2.35+.70*severity+.018*(age-65)+.08*sex+.065*sofa+.16*aki+S5_MORTALITY_LOGOR*h))
  make <- rbinom(N_TOTAL,1,plogis(-1.30+.85*severity+.012*(age-65)+.055*sofa+.28*aki+S5_MAKE_LOGOR*h))
  oracle <- ifelse(rbinom(N_TOTAL,1,1-S5_ORACLE_ACCURACY)==1, 1-hidden, hidden)
  set.seed(seed+100L); ord <- sample.int(N_TOTAL); ntr <- floor(TRAIN_PROP*N_TOTAL); role <- rep("test",N_TOTAL); role[ord[seq_len(ntr)]] <- "train"
  dt <- as.data.table(X); setnames(dt, feature_names)
  dt[, `:=`(
    patient_id = sprintf("%s_R%03d_%06d", substr(scenario,1,2), repeat_id, seq_len(.N)),
    split_role = role, latent_subtype = subtype, severity_true = severity,
    hidden_state = hidden, oracle_label = oracle, age = age, sex = sex,
    sofa = sofa, aki_stage = aki, make30 = make, mortality_30d = mort
  )]
  setcolorder(dt, c("patient_id","split_role","latent_subtype","severity_true","hidden_state","oracle_label",baseline_vars,feature_names,outcome_vars))
  dt
}

log_msg("Script initialized; configuration signature: ", CONFIG_SIGNATURE)

## =============================================================================
## 4. Calibrated missingness injection and realized-fidelity QC
## =============================================================================

calibrate_intercept <- function(score, target) {
  if (target <= 0) return(-Inf)
  uniroot(function(a) mean(plogis(a + score)) - target, c(-40,40), tol = 1e-12)$root
}

target_by_variable <- function(target) {
  if (target == 0) return(setNames(rep(0,length(feature_names)), feature_names))
  p <- pmin(target*missing_weights, .60); p <- p*target/mean(p); p <- pmin(p,.65); p <- p*target/mean(p)
  names(p) <- feature_names; p
}

missing_scores <- function(dt, variable, mechanism, scenario) {
  mar <- .55*stdz(dt$sofa) + .25*stdz(dt$age) + .20*stdz(dt$aki_stage) + .15*(dt$sex-mean(dt$sex))
  own <- stdz(dt[[variable]])
  latent <- if (scenario == "S5_domain2_positive_control") .25*stdz(dt$hidden_state) else
            if (scenario == "S2_true_discrete_subtypes") .20*stdz(dt$latent_subtype) else 0
  mnar <- .80*own + .25*stdz(dt$severity_true) + latent
  selected <- switch(mechanism, MCAR = rep(0,nrow(dt)), MAR = mar, MNAR = mnar,
                     stop("Unknown mechanism: ", mechanism, call. = FALSE))
  list(selected = selected, mar = mar, own = own, mnar = mnar)
}

inject_one_cohort <- function(dt, scenario, mechanism, target, repeat_id, cohort, seed_base) {
  out <- copy(dt); targets <- target_by_variable(target)
  masks <- matrix(FALSE,nrow(out),length(feature_names),dimnames=list(NULL,feature_names))
  rows <- vector("list",length(feature_names))
  for (j in seq_along(feature_names)) {
    v <- feature_names[j]; tv <- targets[[v]]; sc <- missing_scores(dt,v,mechanism,scenario)
    if (tv == 0) { alpha <- -Inf; prob <- rep(0,nrow(dt)); mask <- rep(FALSE,nrow(dt)) } else {
      alpha <- calibrate_intercept(sc$selected,tv); prob <- plogis(alpha+sc$selected)
      # intensity is excluded from the seed: 15% and 30% masks are nested.
      set.seed(seed_base + 1000L*j); mask <- runif(nrow(dt)) < prob
    }
    masks[,j] <- mask; if (any(mask)) set(out,which(mask),v,NA_real_)
    rows[[j]] <- data.table(
      scenario=scenario, missingness_mechanism=mechanism, target_missingness=target,
      repeat_id=repeat_id, cohort=cohort, variable=v,
      target_variable_missingness=tv, expected_probability_mean=mean(prob),
      realized_variable_missingness=mean(mask), target_minus_realized=tv-mean(mask),
      absolute_target_error=abs(tv-mean(mask)), calibration_intercept=alpha,
      selected_signal_auc=if(tv>0) rank_auc(as.integer(mask),sc$selected) else NA_real_,
      mar_signal_auc=if(tv>0) rank_auc(as.integer(mask),sc$mar) else NA_real_,
      own_value_signal_auc=if(tv>0) rank_auc(as.integer(mask),sc$own) else NA_real_,
      selected_signal_correlation=safe_cor(as.integer(mask),sc$selected),
      mar_signal_correlation=safe_cor(as.integer(mask),sc$mar),
      own_value_signal_correlation=safe_cor(as.integer(mask),sc$own),
      missing_cells=sum(mask), eligible_cells=length(mask)
    )
  }
  overall <- data.table(
    scenario=scenario, missingness_mechanism=mechanism, target_missingness=target,
    repeat_id=repeat_id, cohort=cohort, realized_overall_missingness=mean(masks),
    target_minus_realized=target-mean(masks), absolute_target_error=abs(target-mean(masks)),
    missing_cells=sum(masks), eligible_cells=length(masks)
  )
  list(data=out, mask=masks, overall=overall, variable=rbindlist(rows))
}

combine_missing_qc <- function(a,b) {
  o <- rbindlist(list(a$overall,b$overall))
  overall <- o[, .(
    realized_overall_missingness=sum(missing_cells)/sum(eligible_cells),
    target_minus_realized=unique(target_missingness)-sum(missing_cells)/sum(eligible_cells),
    absolute_target_error=abs(unique(target_missingness)-sum(missing_cells)/sum(eligible_cells)),
    missing_cells=sum(missing_cells), eligible_cells=sum(eligible_cells),
    train_realized_missingness=realized_overall_missingness[cohort=="train"],
    test_realized_missingness=realized_overall_missingness[cohort=="test"]
  ), by=.(scenario,missingness_mechanism,target_missingness,repeat_id)]
  v <- rbindlist(list(a$variable,b$variable))
  variable <- v[, .(
    target_variable_missingness=unique(target_variable_missingness),
    expected_probability_mean=weighted.mean(expected_probability_mean,eligible_cells),
    realized_variable_missingness=sum(missing_cells)/sum(eligible_cells),
    target_minus_realized=unique(target_variable_missingness)-sum(missing_cells)/sum(eligible_cells),
    absolute_target_error=abs(unique(target_variable_missingness)-sum(missing_cells)/sum(eligible_cells)),
    selected_signal_auc=weighted.mean(selected_signal_auc,eligible_cells,na.rm=TRUE),
    mar_signal_auc=weighted.mean(mar_signal_auc,eligible_cells,na.rm=TRUE),
    own_value_signal_auc=weighted.mean(own_value_signal_auc,eligible_cells,na.rm=TRUE),
    selected_signal_correlation=weighted.mean(selected_signal_correlation,eligible_cells,na.rm=TRUE),
    mar_signal_correlation=weighted.mean(mar_signal_correlation,eligible_cells,na.rm=TRUE),
    own_value_signal_correlation=weighted.mean(own_value_signal_correlation,eligible_cells,na.rm=TRUE),
    missing_cells=sum(missing_cells), eligible_cells=sum(eligible_cells)
  ), by=.(scenario,missingness_mechanism,target_missingness,repeat_id,variable)]
  for (nm in c("selected_signal_auc","mar_signal_auc","own_value_signal_auc",
               "selected_signal_correlation","mar_signal_correlation","own_value_signal_correlation")) {
    variable[!is.finite(get(nm)), (nm):=NA_real_]
  }
  list(overall=overall,variable=variable)
}

## =============================================================================
## 5. MICE completion QC
## =============================================================================

mice_sets <- function(dt,m,maxit,seed,cohort) {
  cols <- intersect(imputation_vars,names(dt)); x <- as.data.frame(dt[,..cols]); n_before <- sum(is.na(x))
  if (n_before == 0L) {
    sets <- lapply(seq_len(m),function(i){z<-copy(dt);z[,imputation_id:=i];z})
    return(list(ok=TRUE,sets=sets,expected=m,completed=m,n_before=0L,n_after=0L,logged=0L,warnings=""))
  }
  method <- setNames(rep("pmm",length(cols)),cols); method[vapply(x,function(z)!anyNA(z),logical(1))] <- ""
  pred <- mice::make.predictorMatrix(x); diag(pred) <- 0; warns <- character()
  fit <- tryCatch(withCallingHandlers(
    mice::mice(x,m=m,maxit=maxit,method=method,predictorMatrix=pred,seed=seed,printFlag=FALSE),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}
  ),error=function(e) stop("MICE failed in ",cohort,": ",conditionMessage(e),call.=FALSE))
  sets <- lapply(seq_len(m),function(i){comp<-mice::complete(fit,i);z<-copy(dt);for(v in cols)set(z,j=v,value=comp[[v]]);z[,imputation_id:=i];z})
  n_after <- sum(vapply(sets,function(z)sum(is.na(z[,..cols])),numeric(1)))
  if(length(sets)!=m || n_after>0) stop("Incomplete MICE output in ",cohort,call.=FALSE)
  list(ok=TRUE,sets=sets,expected=m,completed=length(sets),n_before=n_before,n_after=n_after,
       logged=if(is.null(fit$loggedEvents))0L else nrow(fit$loggedEvents),
       warnings=paste(unique(warns),collapse=" | "))
}

## =============================================================================
## 6. Preprocessing, k-means, Domain 1 diagnostics
## =============================================================================

fit_prep <- function(dt) {
  bounds <- setNames(lapply(feature_names,function(v)as.numeric(quantile(dt[[v]],WINSOR_PROBS,names=FALSE,type=7))),feature_names)
  X <- sapply(feature_names,function(v){b<-bounds[[v]];pmin(pmax(as.numeric(dt[[v]]),b[1]),b[2])})
  X <- as.matrix(X); colnames(X)<-feature_names; ctr<-colMeans(X); scl<-apply(X,2,sd);scl[!is.finite(scl)|scl<1e-8]<-1
  list(bounds=bounds,center=ctr,scale=scl)
}
apply_prep <- function(dt,prep) {
  X<-sapply(feature_names,function(v){b<-prep$bounds[[v]];pmin(pmax(as.numeric(dt[[v]]),b[1]),b[2])})
  X<-as.matrix(X);colnames(X)<-feature_names;X<-sweep(X,2,prep$center,"-");X<-sweep(X,2,prep$scale,"/");X[!is.finite(X)]<-0;X
}

safe_kmeans <- function(X,k,seed,nstart=100L) {
  warns<-character();set.seed(seed)
  fit<-tryCatch(withCallingHandlers(kmeans(X,k,nstart=nstart,iter.max=500,algorithm="Hartigan-Wong"),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)NULL)
  fallback<-is.null(fit)||any(grepl("Quick-TRANSfer|maximum number of steps",warns,ignore.case=TRUE))
  if(fallback){set.seed(seed+100000L);fit<-tryCatch(suppressWarnings(kmeans(X,k,nstart=nstart,iter.max=500,algorithm="Lloyd")),error=function(e)NULL)}
  if(is.null(fit)||length(fit$cluster)!=nrow(X)||anyNA(fit$cluster)||length(unique(fit$cluster))!=k)
    stop("Invalid k-means K=",k,call.=FALSE)
  list(labels=as.integer(fit$cluster),centers=as.matrix(fit$centers),fallback=fallback,warnings=paste(unique(warns),collapse=" | "))
}
k2_map <- function(centers){score<-as.numeric(centers%*%severity_loading[colnames(centers)]);hi<-which.max(score);map<-integer(2);map[hi]<-1L;map[setdiff(1:2,hi)]<-2L;map}
map_labels <- function(x,map) as.integer(map[x])
assign_centroids <- function(X,C){D<-sapply(seq_len(nrow(C)),function(g)rowSums((X-matrix(C[g,],nrow(X),ncol(X),byrow=TRUE))^2));if(is.null(dim(D)))D<-matrix(D,ncol=1);as.integer(max.col(-D,ties.method="first"))}
outcome_sep <- function(dt,lab,outcome) safe_mean(dt[[outcome]][lab==1])-safe_mean(dt[[outcome]][lab==2])

run_dip <- function(X,labels,seed) {
  n<-min(nrow(X),DIP_N);set.seed(seed);idx<-sort(sample.int(nrow(X),n));Xs<-X[idx,,drop=FALSE];l<-labels[idx]
  C<-rbind(colMeans(Xs[l==1,,drop=FALSE]),colMeans(Xs[l==2,,drop=FALSE]));d<-C[1,]-C[2,];if(sum(d^2)<1e-12)d[1]<-1
  axis<-as.numeric(Xs%*%d);pc<-prcomp(Xs,center=FALSE,scale.=FALSE);mp<-min(5,ncol(pc$x))
  p<-c(discriminant_axis=diptest::dip.test(axis)$p.value,vapply(seq_len(mp),function(j)diptest::dip.test(pc$x[,j])$p.value,numeric(1)))
  names(p)[-1]<-paste0("PC",seq_len(mp));q<-p.adjust(p,"BH");tab<-data.table(axis=names(p),p_value=as.numeric(p),p_fdr=as.numeric(q))
  list(table=tab,disc_p=tab[axis=="discriminant_axis",p_value][1],disc_q=tab[axis=="discriminant_axis",p_fdr][1],any_fdr=any(tab$p_fdr<.05))
}
cluster_index <- function(X,l){ctr<-colMeans(X);tot<-sum((X-matrix(ctr,nrow(X),ncol(X),byrow=TRUE))^2);within<-0;for(g in unique(l)){idx<-which(l==g);c<-colMeans(X[idx,,drop=FALSE]);within<-within+sum((X[idx,,drop=FALSE]-matrix(c,length(idx),ncol(X),byrow=TRUE))^2)};within/tot}
run_sig <- function(X,seed) {
  n<-min(nrow(X),SIGCLUST_N);set.seed(seed);idx<-sort(sample.int(nrow(X),n));Xs<-X[idx,,drop=FALSE]
  obsfit<-safe_kmeans(Xs,2,seed+1,50);obs<-cluster_index(Xs,obsfit$labels)
  sv<-svd(Xs,nu=0,nv=min(ncol(Xs),nrow(Xs)-1));eig<-pmax(sv$d^2/max(1,nrow(Xs)-1),1e-8);q<-length(eig)
  nul<-rep(NA_real_,SIGCLUST_B);for(b in seq_len(SIGCLUST_B)){set.seed(seed+1000+b);Z<-sweep(matrix(rnorm(n*q),nrow=n),2,sqrt(eig),"*");f<-safe_kmeans(Z,2,seed+50000+b,20);nul[b]<-cluster_index(Z,f$labels)}
  nul<-nul[is.finite(nul)];if(length(nul)<max(5,floor(.6*SIGCLUST_B)))stop("Too few SigClust null fits",call.=FALSE)
  mu<-mean(nul);s<-sd(nul);effect<-if(is.finite(s)&&s>0)(mu-obs)/s else NA_real_;p<-(1+sum(nul<=obs))/(length(nul)+1)
  data.table(observed_index=obs,null_mean=mu,null_sd=s,standardized_effect_sd=effect,
             relative_effect_pct=100*(mu-obs)/mu,auxiliary_p=p,B_successful=length(nul))
}

## =============================================================================
## 7. Domain 2 prediction and Rubin pooling
## =============================================================================

safe_auc <- function(y,p){ok<-is.finite(y)&is.finite(p);if(sum(ok)<10||length(unique(y[ok]))<2)return(NA_real_);as.numeric(suppressMessages(pROC::auc(pROC::roc(y[ok],p[ok],quiet=TRUE,direction="<"))))}
delong_var <- function(y,p0,p1) tryCatch({ok<-is.finite(y)&is.finite(p0)&is.finite(p1);r0<-pROC::roc(y[ok],p0[ok],quiet=TRUE,direction="<");r1<-pROC::roc(y[ok],p1[ok],quiet=TRUE,direction="<");max(as.numeric(pROC::var.roc(r1,method="delong")+pROC::var.roc(r0,method="delong")-2*pROC::cov.roc(r1,r0,method="delong")),0)},error=function(e)NA_real_)
continuous_nri <- function(y,p0,p1){ok<-is.finite(y)&is.finite(p0)&is.finite(p1);y<-as.integer(y[ok]);p0<-p0[ok];p1<-p1[ok];if(length(unique(y))<2)return(NA_real_);up<-p1>p0;down<-p1<p0;(mean(up[y==1])-mean(down[y==1]))+(mean(down[y==0])-mean(up[y==0]))}
idi <- function(y,p0,p1){ok<-is.finite(y)&is.finite(p0)&is.finite(p1);y<-as.integer(y[ok]);p0<-p0[ok];p1<-p1[ok];if(length(unique(y))<2)return(NA_real_);(mean(p1[y==1])-mean(p1[y==0]))-(mean(p0[y==1])-mean(p0[y==0]))}
make_foldid <- function(y,seed){nf<-min(GLMNET_NFOLD,max(3L,floor(length(y)/50)));set.seed(seed);f<-integer(length(y));for(v in unique(y)){idx<-which(y==v);f[idx]<-sample(rep(seq_len(nf),length.out=length(idx)))};f}
fit_en <- function(X,y,Xnew,foldid){fit<-glmnet::cv.glmnet(as.matrix(X),as.numeric(y),family="binomial",alpha=GLMNET_ALPHA,foldid=foldid,type.measure="deviance",standardize=TRUE);as.numeric(predict(fit,as.matrix(Xnew),s="lambda.min",type="response"))}
fit_glm <- function(train,test,outcome,predictors){f<-reformulate(predictors,response=outcome);fit<-suppressWarnings(glm(f,data=train,family=binomial()));as.numeric(predict(fit,newdata=test,type="response"))}

prediction_one_outcome <- function(train,test,Xtr,Xte,label_tr,label_te,label_source,outcome,seed) {
  tr<-as.data.frame(train[,c(baseline_vars,outcome),with=FALSE]);te<-as.data.frame(test[,c(baseline_vars,outcome),with=FALSE])
  tr$label<-factor(label_tr);te$label<-factor(label_te,levels=levels(tr$label))
  p_base<-fit_glm(tr,te,outcome,baseline_vars);p_base_l<-fit_glm(tr,te,outcome,c(baseline_vars,"label"))
  Btr<-model.matrix(~age+sex+sofa+aki_stage-1,tr);Bte<-model.matrix(~age+sex+sofa+aki_stage-1,te)
  Rtr<-cbind(Btr,Xtr);Rte<-cbind(Bte,Xte);RLtr<-cbind(Rtr,label=as.numeric(label_tr==1));RLte<-cbind(Rte,label=as.numeric(label_te==1))
  ytr<-train[[outcome]];yte<-test[[outcome]];fold<-make_foldid(ytr,seed)
  p_raw<-fit_en(Rtr,ytr,Rte,fold);p_raw_l<-fit_en(RLtr,ytr,RLte,fold)
  rowfun<-function(ref,aug,p0,p1){a0<-safe_auc(yte,p0);a1<-safe_auc(yte,p1);data.table(
    outcome=outcome,label_source=label_source,reference_model=ref,augmented_model=aug,
    contrast=paste0(aug,"_minus_",ref),auc_reference=a0,auc_augmented=a1,delta_auc=a1-a0,
    within_imputation_variance=delong_var(yte,p0,p1),continuous_nri=continuous_nri(yte,p0,p1),idi=idi(yte,p0,p1),
    brier_reference=mean((yte-p0)^2),brier_augmented=mean((yte-p1)^2),delta_brier=mean((yte-p1)^2)-mean((yte-p0)^2))}
  rbindlist(list(rowfun("Base","Base+label",p_base,p_base_l),rowfun("Raw-EN","Raw-EN+label",p_raw,p_raw_l)))
}

rubin_pool <- function(dt) dt[, {
  q<-delta_auc[is.finite(delta_auc)];u<-within_imputation_variance[is.finite(delta_auc)];m<-length(q)
  if(!m) list(m=0L,pooled_auc_reference=NA_real_,pooled_auc_augmented=NA_real_,pooled_delta_auc=NA_real_,within_variance_mean=NA_real_,between_variance=NA_real_,total_variance=NA_real_,standard_error=NA_real_,ci_low=NA_real_,ci_high=NA_real_,ci_positive=NA,practical_null=NA,pooled_continuous_nri=NA_real_,pooled_idi=NA_real_,pooled_delta_brier=NA_real_) else {
    qbar<-mean(q);ubar<-if(any(is.finite(u)))mean(u[is.finite(u)]) else 0;b<-if(m>1)var(q) else 0;t<-ubar+(1+1/m)*b;se<-sqrt(max(t,0))
    list(m=m,pooled_auc_reference=safe_mean(auc_reference),pooled_auc_augmented=safe_mean(auc_augmented),pooled_delta_auc=qbar,
         within_variance_mean=ubar,between_variance=b,total_variance=t,standard_error=se,ci_low=qbar-1.96*se,ci_high=qbar+1.96*se,
         ci_positive=(qbar-1.96*se)>0,practical_null=abs(qbar)<D2_PRACTICAL_NULL,
         pooled_continuous_nri=safe_mean(continuous_nri),pooled_idi=safe_mean(idi),pooled_delta_brier=safe_mean(delta_brier))
  }
}, by=.(task_id,scenario,missingness_mechanism,target_missingness,repeat_id,outcome,label_source,reference_model,augmented_model,contrast)]

## =============================================================================
## 8. One task = one scenario/mechanism/intensity/repeat
## =============================================================================

run_task <- function(task) {
  task_id<-as.character(task$task_id);scenario<-as.character(task$scenario);mechanism<-as.character(task$missingness_mechanism)
  target<-as.numeric(task$target_missingness);rep_id<-as.integer(task$repeat_id);base_seed<-as.integer(task$scenario_repeat_seed);mask_seed<-as.integer(task$mask_seed_base)
  complete<-generate_complete(scenario,rep_id,base_seed);tr0<-complete[split_role=="train"];te0<-complete[split_role=="test"]
  trm<-inject_one_cohort(tr0,scenario,mechanism,target,rep_id,"train",mask_seed+100000L)
  tem<-inject_one_cohort(te0,scenario,mechanism,target,rep_id,"test",mask_seed+200000L)
  realized<-combine_missing_qc(trm,tem)
  mice_seed_offset<-as.integer(round(target*1000))+10000L*match(mechanism,MECHANISMS)
  tri<-mice_sets(trm$data,MICE_M,MICE_MAXIT,base_seed+300000L+mice_seed_offset,"train")
  tei<-mice_sets(tem$data,MICE_M,MICE_MAXIT,base_seed+400000L+mice_seed_offset,"test")
  mice_status<-rbindlist(list(
    data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,cohort="train",expected_imputations=tri$expected,completed_imputations=tri$completed,n_missing_before=tri$n_before,n_missing_after=tri$n_after,n_logged_events=tri$logged,warning_count=if(nzchar(tri$warnings))length(strsplit(tri$warnings," \\| ")[[1]]) else 0L,warnings=tri$warnings,mice_success=tri$ok),
    data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,cohort="test",expected_imputations=tei$expected,completed_imputations=tei$completed,n_missing_before=tei$n_before,n_missing_after=tei$n_after,n_logged_events=tei$logged,warning_count=if(nzchar(tei$warnings))length(strsplit(tei$warnings," \\| ")[[1]]) else 0L,warnings=tei$warnings,mice_success=tei$ok)
  ))
  k2labs<-vector("list",MICE_M);k3labs<-vector("list",MICE_M);d1<-list();d2<-list();dipaxes<-NULL;dip<-NULL;sig<-NULL
  for(imp in seq_len(MICE_M)){
    tr<-tri$sets[[imp]];te<-tei$sets[[imp]];prep<-fit_prep(tr);Xtr<-apply_prep(tr,prep);Xte<-apply_prep(te,prep)
    f2<-safe_kmeans(Xtr,2,base_seed+500000L+1000L*imp);f3<-safe_kmeans(Xtr,3,base_seed+600000L+1000L*imp)
    map<-k2_map(f2$centers);lab2<-map_labels(f2$labels,map);C2<-f2$centers[order(map),,drop=FALSE];lab2te<-assign_centroids(Xte,C2)
    k2labs[[imp]]<-lab2;k3labs[[imp]]<-f3$labels
    d1[[imp]]<-rbindlist(list(
      data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,imputation_id=imp,k=2L,ari_vs_latent_subtype=adjusted_rand(lab2,tr$latent_subtype),mortality_30d_C1_minus_C2=outcome_sep(tr,lab2,"mortality_30d"),make30_C1_minus_C2=outcome_sep(tr,lab2,"make30"),cluster_1_prevalence=mean(lab2==1),cluster_min_prevalence=min(table(lab2))/length(lab2),lloyd_fallback=f2$fallback,kmeans_warnings=f2$warnings),
      data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,imputation_id=imp,k=3L,ari_vs_latent_subtype=adjusted_rand(f3$labels,tr$latent_subtype),mortality_30d_C1_minus_C2=NA_real_,make30_C1_minus_C2=NA_real_,cluster_1_prevalence=mean(f3$labels==1),cluster_min_prevalence=min(table(f3$labels))/length(f3$labels),lloyd_fallback=f3$fallback,kmeans_warnings=f3$warnings)
    ))
    if(scenario=="S5_domain2_positive_control"){ltr<-tr$oracle_label;lte<-te$oracle_label;lsource<-"fixed_noisy_oracle"}else{ltr<-lab2;lte<-lab2te;lsource<-"K2_cluster"}
    one<-rbindlist(lapply(seq_along(outcome_vars),function(oi)prediction_one_outcome(tr,te,Xtr,Xte,ltr,lte,lsource,outcome_vars[oi],base_seed+700000L+10000L*oi+imp)))
    one[,`:=`(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,imputation_id=imp)];d2[[imp]]<-one
    if(imp==PRIMARY_IMPUTATION_ID){dip<-run_dip(Xtr,lab2,base_seed+800000L);sig<-run_sig(Xtr,base_seed+900000L);dipaxes<-copy(dip$table);dipaxes[,`:=`(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id)]}
  }
  d1imp<-safe_rbind(d1);d2imp<-safe_rbind(d2)
  stability<-rbindlist(lapply(c(2L,3L),function(k){L<-if(k==2L)k2labs else k3labs;a<-vapply(seq_len(MICE_M),function(i)adjusted_rand(L[[1]],L[[i]]),numeric(1));data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,k=k,m=MICE_M,ari_mean_vs_primary=mean(a,na.rm=TRUE),ari_median_vs_primary=median(a,na.rm=TRUE),ari_min_vs_primary=min(a,na.rm=TRUE),ari_prop_ge_0_90=mean(a>=MICE_ARI_CUTOFF,na.rm=TRUE))}))
  d1rep<-d1imp[,.(ari_vs_latent_subtype=safe_mean(ari_vs_latent_subtype),mortality_30d_C1_minus_C2=safe_mean(mortality_30d_C1_minus_C2),make30_C1_minus_C2=safe_mean(make30_C1_minus_C2),cluster_1_prevalence=safe_mean(cluster_1_prevalence),cluster_min_prevalence=safe_mean(cluster_min_prevalence),lloyd_fallback_rate=mean(as.numeric(lloyd_fallback),na.rm=TRUE)),by=.(task_id,scenario,missingness_mechanism,target_missingness,repeat_id,k)]
  diag<-data.table(task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,dip_discriminant_p=dip$disc_p,dip_discriminant_p_fdr=dip$disc_q,dip_any_axis_fdr_lt_0_05=dip$any_fdr,sigclust_observed_index=sig$observed_index,sigclust_null_mean=sig$null_mean,sigclust_null_sd=sig$null_sd,sigclust_standardized_effect_sd=sig$standardized_effect_sd,sigclust_relative_effect_pct=sig$relative_effect_pct,sigclust_auxiliary_p=sig$auxiliary_p,sigclust_B_successful=sig$B_successful)
  d2rep<-rubin_pool(d2imp)
  variability<-d2imp[,.(auc_reference_across_imputation_sd=safe_sd(auc_reference),auc_augmented_across_imputation_sd=safe_sd(auc_augmented),delta_auc_across_imputation_sd=safe_sd(delta_auc),continuous_nri_across_imputation_sd=safe_sd(continuous_nri),idi_across_imputation_sd=safe_sd(idi)),by=.(task_id,scenario,missingness_mechanism,target_missingness,repeat_id,outcome,label_source,contrast)]
  list(ok=TRUE,config_signature=CONFIG_SIGNATURE,task_id=task_id,scenario=scenario,missingness_mechanism=mechanism,target_missingness=target,repeat_id=rep_id,
       realized_overall=realized$overall,realized_variable=realized$variable,mice_status=mice_status,mice_stability=stability,model_variability=variability,
       domain1_imp=d1imp,domain1_repeat=d1rep,domain1_diag=diag,dip_axes=dipaxes,domain2_imp=d2imp,domain2_repeat=d2rep)
}

## =============================================================================
## 9. Manifest, multisession execution, checkpointing
## =============================================================================

tasks <- CJ(scenario=SCENARIOS,missingness_mechanism=MECHANISMS,target_missingness=TARGET_GRID,repeat_id=seq_len(N_ACTIVE_REP))
tasks[,`:=`(scenario_index=match(scenario,SCENARIOS),mechanism_index=match(missingness_mechanism,MECHANISMS))]
tasks[,scenario_repeat_seed:=GLOBAL_SEED+scenario_index*10000000L+repeat_id*10000L]
tasks[,mask_seed_base:=GLOBAL_SEED+scenario_index*10000000L+repeat_id*10000L+mechanism_index*1000000L]
tasks[,task_id:=sprintf("%s_%s_miss_%02d_rep_%03d",c(S1_pure_severity_continuum="S1",S2_true_discrete_subtypes="S2",S5_domain2_positive_control="S5")[scenario],missingness_mechanism,as.integer(round(target_missingness*100)),repeat_id)]
fwrite(tasks,file.path(TAB_DIR,"sensitivity05_task_manifest.csv"))
config <- data.table(parameter=c("ANALYSIS_ID","RUN_MODE","SCENARIOS","MECHANISMS","TARGET_GRID","N_REP_PER_SETTING","N_TOTAL","TRAIN_PROP","MICE_M","MICE_MAXIT","SIGCLUST_B","S5_ORACLE_ACCURACY","N_WORKERS","GLOBAL_SEED","CONFIG_SIGNATURE"),
                     value=c(ANALYSIS_ID,RUN_MODE,paste(SCENARIOS,collapse="|"),paste(MECHANISMS,collapse="|"),paste(TARGET_GRID,collapse="|"),N_ACTIVE_REP,N_TOTAL,TRAIN_PROP,MICE_M,MICE_MAXIT,SIGCLUST_B,S5_ORACLE_ACCURACY,N_WORKERS,GLOBAL_SEED,CONFIG_SIGNATURE))
fwrite(config,file.path(TAB_DIR,"sensitivity05_configuration.csv"))

if (N_WORKERS > 1L) future::plan(future::multisession, workers = N_WORKERS) else future::plan(future::sequential)
log_msg("Total tasks: ",nrow(tasks))
worker_fun <- function(i){
  suppressPackageStartupMessages({library(data.table);library(mice);library(mclust);library(diptest);library(cluster);library(glmnet);library(pROC)})
  data.table::setDTthreads(1L);task<-tasks[i];path<-file.path(CKPT_DIR,paste0(task$task_id,".rds"))
  if(RESUME_COMPLETED_TASKS&&file.exists(path)){old<-tryCatch(readRDS(path),error=function(e)NULL);if(is.list(old)&&isTRUE(old$ok)&&identical(old$config_signature,CONFIG_SIGNATURE))return(old)}
  start<-Sys.time();res<-tryCatch({z<-run_task(task);z$started_at<-format(start,"%Y-%m-%d %H:%M:%S");z$finished_at<-timestamp();z$elapsed_minutes<-as.numeric(difftime(Sys.time(),start,units="mins"));z},
    error=function(e)list(ok=FALSE,config_signature=CONFIG_SIGNATURE,task_id=as.character(task$task_id),scenario=as.character(task$scenario),missingness_mechanism=as.character(task$missingness_mechanism),target_missingness=as.numeric(task$target_missingness),repeat_id=as.integer(task$repeat_id),started_at=format(start,"%Y-%m-%d %H:%M:%S"),finished_at=timestamp(),elapsed_minutes=as.numeric(difftime(Sys.time(),start,units="mins")),error=conditionMessage(e)))
  saveRDS(res,path,compress="xz");res
}
results <- future.apply::future_lapply(seq_len(nrow(tasks)),worker_fun,future.seed=TRUE,future.scheduling=1,
  future.packages=c("data.table","mice","mclust","diptest","cluster","glmnet","pROC"))
future::plan(future::sequential);invisible(gc(full=TRUE))

## ---- completion QC -----------------------------------------------------------
errors <- safe_rbind(lapply(results,function(z)if(is.list(z)&&!isTRUE(z$ok))data.table(task_id=z$task_id,scenario=z$scenario,missingness_mechanism=z$missingness_mechanism,target_missingness=z$target_missingness,repeat_id=z$repeat_id,started_at=z$started_at,finished_at=z$finished_at,elapsed_minutes=z$elapsed_minutes,error=z$error)else NULL))
if (!nrow(errors)) {
  errors <- data.table(
    task_id=character(), scenario=character(), missingness_mechanism=character(),
    target_missingness=numeric(), repeat_id=integer(), started_at=character(),
    finished_at=character(), elapsed_minutes=numeric(), error=character()
  )
}
fwrite(errors,file.path(LOG_DIR,"sensitivity05_task_failure_log.csv"))
success <- Filter(function(z)is.list(z)&&isTRUE(z$ok)&&identical(z$config_signature,CONFIG_SIGNATURE),results)
if(!length(success))stop("All SA-05 tasks failed; see failure log.",call.=FALSE)
obs <- safe_rbind(lapply(success,function(z)data.table(scenario=z$scenario,missingness_mechanism=z$missingness_mechanism,target_missingness=z$target_missingness,repeat_id=z$repeat_id)))
dup <- obs[,.N,by=.(scenario,missingness_mechanism,target_missingness,repeat_id)][N>1]
expected <- tasks[,.(scenario,missingness_mechanism,target_missingness,repeat_id)];missing <- fsetdiff(expected,unique(obs))
completion <- tasks[,.(expected_repeats=uniqueN(repeat_id)),by=.(scenario,missingness_mechanism,target_missingness)]
completion <- merge(completion,obs[,.(observed_repeats=uniqueN(repeat_id)),by=.(scenario,missingness_mechanism,target_missingness)],all.x=TRUE)
completion <- merge(completion,errors[,.(failed_repeats=uniqueN(repeat_id)),by=.(scenario,missingness_mechanism,target_missingness)],all.x=TRUE)
completion <- merge(completion,dup[,.(duplicate_repeats=uniqueN(repeat_id)),by=.(scenario,missingness_mechanism,target_missingness)],all.x=TRUE)
for(v in c("observed_repeats","failed_repeats","duplicate_repeats"))completion[is.na(get(v)),(v):=0L]
completion[,`:=`(missing_repeats=expected_repeats-observed_repeats,completion_proportion=observed_repeats/expected_repeats)]
completion[,structural_qc_pass:=observed_repeats==expected_repeats&failed_repeats==0&duplicate_repeats==0&missing_repeats==0]
fwrite(completion,file.path(QC_DIR,"05_run_completion_qc.csv"));fwrite(dup,file.path(QC_DIR,"05_duplicate_repeat_records.csv"));fwrite(missing,file.path(QC_DIR,"05_missing_repeat_records.csv"))
log_msg("Successful tasks: ",length(success)," / ",nrow(tasks),"; failed: ",nrow(errors))

## ---- collect repeat-level objects -------------------------------------------
realized_overall <- safe_rbind(lapply(success,`[[`,"realized_overall"))
realized_variable <- safe_rbind(lapply(success,`[[`,"realized_variable"))
mice_status <- safe_rbind(lapply(success,`[[`,"mice_status"))
mice_stability <- safe_rbind(lapply(success,`[[`,"mice_stability"))
model_variability <- safe_rbind(lapply(success,`[[`,"model_variability"))
d1imp <- safe_rbind(lapply(success,`[[`,"domain1_imp"));d1rep <- safe_rbind(lapply(success,`[[`,"domain1_repeat"));d1diag <- safe_rbind(lapply(success,`[[`,"domain1_diag"));dip_axes <- safe_rbind(lapply(success,`[[`,"dip_axes"))
d2imp <- safe_rbind(lapply(success,`[[`,"domain2_imp"));d2rep <- safe_rbind(lapply(success,`[[`,"domain2_repeat"))

## =============================================================================
## 10. Realized-missingness fidelity QC
## =============================================================================

realized_summary <- realized_overall[,.(repeats=.N,realized_mean=safe_mean(realized_overall_missingness),realized_sd=safe_sd(realized_overall_missingness),realized_median=median(realized_overall_missingness),realized_q025=safe_q(realized_overall_missingness,.025),realized_q975=safe_q(realized_overall_missingness,.975),mcse_mean=safe_sd(realized_overall_missingness)/sqrt(.N),mean_target_minus_realized=safe_mean(target_minus_realized),mean_absolute_target_error=safe_mean(absolute_target_error),q975_absolute_target_error=safe_q(absolute_target_error,.975)),by=.(scenario,missingness_mechanism,target_missingness)]
variable_summary <- realized_variable[,.(repeats=.N,target_variable_missingness=unique(target_variable_missingness),realized_mean=safe_mean(realized_variable_missingness),realized_sd=safe_sd(realized_variable_missingness),realized_q025=safe_q(realized_variable_missingness,.025),realized_q975=safe_q(realized_variable_missingness,.975),mcse_mean=safe_sd(realized_variable_missingness)/sqrt(.N),mean_absolute_target_error=safe_mean(absolute_target_error),mean_selected_signal_auc=safe_mean(selected_signal_auc),mean_mar_signal_auc=safe_mean(mar_signal_auc),mean_own_value_signal_auc=safe_mean(own_value_signal_auc),mean_abs_selected_signal_correlation=safe_mean(abs(selected_signal_correlation))),by=.(scenario,missingness_mechanism,target_missingness,variable)]
mech_fidelity_rep <- realized_variable[target_missingness>0,.(
  mean_selected_signal_auc=safe_mean(selected_signal_auc),
  mean_abs_selected_signal_correlation=safe_mean(abs(selected_signal_correlation)),
  mean_mar_signal_auc=safe_mean(mar_signal_auc),
  mean_own_value_signal_auc=safe_mean(own_value_signal_auc),
  mcar_mean_abs_nonconstant_signal_correlation=safe_mean(c(abs(mar_signal_correlation),abs(own_value_signal_correlation)))
),by=.(scenario,missingness_mechanism,target_missingness,repeat_id)]
mech_fidelity <- mech_fidelity_rep[,.(
  repeats=.N,mean_selected_signal_auc=safe_mean(mean_selected_signal_auc),
  q025_selected_signal_auc=safe_q(mean_selected_signal_auc,.025),q975_selected_signal_auc=safe_q(mean_selected_signal_auc,.975),
  mean_abs_selected_signal_correlation=safe_mean(mean_abs_selected_signal_correlation),
  mean_mar_signal_auc=safe_mean(mean_mar_signal_auc),mean_own_value_signal_auc=safe_mean(mean_own_value_signal_auc),
  mcar_mean_abs_nonconstant_signal_correlation=safe_mean(mcar_mean_abs_nonconstant_signal_correlation)
),by=.(scenario,missingness_mechanism,target_missingness)]
mech_fidelity[,fidelity_pass:=fifelse(
  missingness_mechanism=="MCAR",mcar_mean_abs_nonconstant_signal_correlation<=MCAR_ABS_COR_MAX,
  fifelse(missingness_mechanism=="MAR",mean_selected_signal_auc>=MAR_SIGNAL_AUC_MIN,mean_selected_signal_auc>=MNAR_SIGNAL_AUC_MIN)
)]
balance_rep <- realized_overall[,.(mechanism_min=min(realized_overall_missingness),mechanism_max=max(realized_overall_missingness),mechanism_range=max(realized_overall_missingness)-min(realized_overall_missingness)),by=.(scenario,target_missingness,repeat_id)]
balance_summary <- balance_rep[,.(repeats=.N,mean_mechanism_range=safe_mean(mechanism_range),q975_mechanism_range=safe_q(mechanism_range,.975),mechanism_balance_pass=safe_mean(mechanism_range)<=MECHANISM_RANGE_TOL&safe_q(mechanism_range,.975)<=MECHANISM_Q975_RANGE_TOL),by=.(scenario,target_missingness)]

mon <- dcast(realized_overall,scenario+missingness_mechanism+repeat_id~target_missingness,value.var="realized_overall_missingness")
if(all(c("0","0.15","0.3")%in%names(mon))){setnames(mon,c("0","0.15","0.3"),c("realized_0","realized_15","realized_30"));mon[,monotonic_pass:=realized_0<=realized_15&realized_15<=realized_30]}else mon[,monotonic_pass:=FALSE]
mon_summary <- mon[,.(repeats=.N,monotonic_pass_rate=mean(monotonic_pass),all_repeats_monotonic=all(monotonic_pass)),by=.(scenario,missingness_mechanism)]

## zero-intensity equivalence across mechanisms
z <- dcast(d1rep[target_missingness==0],scenario+repeat_id+k~missingness_mechanism,value.var=c("ari_vs_latent_subtype","mortality_30d_C1_minus_C2","cluster_1_prevalence"))
zero_eq <- z[,.(scenario,repeat_id,k,max_abs_difference=max(abs(c(
  ari_vs_latent_subtype_MCAR-ari_vs_latent_subtype_MAR,ari_vs_latent_subtype_MCAR-ari_vs_latent_subtype_MNAR,
  mortality_30d_C1_minus_C2_MCAR-mortality_30d_C1_minus_C2_MAR,mortality_30d_C1_minus_C2_MCAR-mortality_30d_C1_minus_C2_MNAR,
  cluster_1_prevalence_MCAR-cluster_1_prevalence_MAR,cluster_1_prevalence_MCAR-cluster_1_prevalence_MNAR)),na.rm=TRUE))]
zero_eq[!is.finite(max_abs_difference),max_abs_difference:=0];zero_eq[,equivalence_pass:=max_abs_difference<=ZERO_EQUIVALENCE_TOL]

z2 <- dcast(
  d2rep[target_missingness == 0],
  scenario + repeat_id + outcome + label_source + contrast ~ missingness_mechanism,
  value.var = "pooled_delta_auc"
)
if (all(c("MCAR","MAR","MNAR") %in% names(z2))) {
  zero_eq_d2 <- z2[, .(
    scenario, repeat_id, outcome, label_source, contrast,
    max_abs_difference = pmax(
      abs(MCAR - MAR), abs(MCAR - MNAR), abs(MAR - MNAR),
      na.rm = TRUE
    )
  )]
  zero_eq_d2[!is.finite(max_abs_difference), max_abs_difference := 0]
  zero_eq_d2[, equivalence_pass := max_abs_difference <= ZERO_EQUIVALENCE_TOL]
} else {
  zero_eq_d2 <- data.table(equivalence_pass = FALSE)
}

setting_fidelity <- merge(realized_summary,mech_fidelity,by=c("scenario","missingness_mechanism","target_missingness"),all.x=TRUE)
setting_fidelity[,overall_target_fidelity_pass:=mean_absolute_target_error<=OVERALL_ABS_ERROR_TOL&q975_absolute_target_error<=OVERALL_Q975_ABS_ERROR_TOL]
setting_fidelity[target_missingness==0,fidelity_pass:=TRUE]
setting_fidelity[,setting_fidelity_pass:=overall_target_fidelity_pass&fidelity_pass]

## =============================================================================
## 11. MICE, Domain 1, Domain 2 summaries and contrasts
## =============================================================================

mice_completion <- mice_status[,.(expected_imputations=sum(expected_imputations),completed_imputations=sum(completed_imputations),mice_failure_rate=mean(!mice_success),convergence_warning_rate=mean(warning_count>0|n_logged_events>0),incomplete_chain_rate=mean(n_missing_after>0|completed_imputations<expected_imputations),mean_logged_events=safe_mean(n_logged_events),residual_missing_cells=sum(n_missing_after)),by=.(scenario,missingness_mechanism,target_missingness,cohort)]
mice_stability_summary <- mice_stability[,.(nonmissing_repeats=sum(is.finite(ari_mean_vs_primary)),mean_ari=safe_mean(ari_mean_vs_primary),empirical_sd=safe_sd(ari_mean_vs_primary),median=median(ari_mean_vs_primary,na.rm=TRUE),q025=safe_q(ari_mean_vs_primary,.025),q975=safe_q(ari_mean_vs_primary,.975),proportion_ari_ge_0_90=mean(ari_mean_vs_primary>=MICE_ARI_CUTOFF,na.rm=TRUE),mcse_mean=safe_sd(ari_mean_vs_primary)/sqrt(sum(is.finite(ari_mean_vs_primary)))),by=.(scenario,missingness_mechanism,target_missingness,k)]
model_variability_summary <- model_variability[,.(repeats=.N,mean_auc_reference_across_imputation_sd=safe_mean(auc_reference_across_imputation_sd),mean_auc_augmented_across_imputation_sd=safe_mean(auc_augmented_across_imputation_sd),mean_delta_auc_across_imputation_sd=safe_mean(delta_auc_across_imputation_sd),mean_continuous_nri_across_imputation_sd=safe_mean(continuous_nri_across_imputation_sd),mean_idi_across_imputation_sd=safe_mean(idi_across_imputation_sd)),by=.(scenario,missingness_mechanism,target_missingness,outcome,label_source,contrast)]

d1 <- merge(d1rep,d1diag,by=c("task_id","scenario","missingness_mechanism","target_missingness","repeat_id"),all.x=TRUE)
d1 <- merge(d1,mice_stability[,.(task_id,scenario,missingness_mechanism,target_missingness,repeat_id,k,mice_ari_mean_vs_primary=ari_mean_vs_primary,mice_ari_min_vs_primary=ari_min_vs_primary)],by=c("task_id","scenario","missingness_mechanism","target_missingness","repeat_id","k"),all.x=TRUE)

d1_summary <- safe_rbind(list(
  summary_numeric(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],c("scenario","missingness_mechanism","target_missingness","k"),"ari_vs_latent_subtype","ARI_vs_latent_subtype"),
  summary_numeric(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")&k==2L],c("scenario","missingness_mechanism","target_missingness"),"mortality_30d_C1_minus_C2","K2_mortality_separation"),
  summary_numeric(unique(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes"),.(scenario,missingness_mechanism,target_missingness,repeat_id,sigclust_standardized_effect_sd)]),c("scenario","missingness_mechanism","target_missingness"),"sigclust_standardized_effect_sd","SigClust_like_standardized_effect_SD"),
  summary_numeric(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],c("scenario","missingness_mechanism","target_missingness","k"),"mice_ari_mean_vs_primary","MICE_label_stability_ARI")
))

s1diag <- unique(d1[scenario=="S1_pure_severity_continuum",.(scenario,missingness_mechanism,target_missingness,repeat_id,dip_discriminant_p_fdr,sigclust_standardized_effect_sd)])
s2diag <- unique(d1[scenario=="S2_true_discrete_subtypes",.(scenario,missingness_mechanism,target_missingness,repeat_id,dip_discriminant_p_fdr,sigclust_standardized_effect_sd)])
d1_oc <- safe_rbind(list(
  summary_binary(s1diag[,flag:=dip_discriminant_p_fdr<DIP_ALPHA],c("scenario","missingness_mechanism","target_missingness"),"flag","S1_Dip_false_positive_rate",N_ACTIVE_REP),
  summary_binary(s1diag[,flag2:=sigclust_standardized_effect_sd>=SIGCLUST_EFFECT_CUTOFF],c("scenario","missingness_mechanism","target_missingness"),"flag2","S1_standardized_effect_false_positive_rate",N_ACTIVE_REP),
  summary_binary(s2diag[,flag:=dip_discriminant_p_fdr<DIP_ALPHA],c("scenario","missingness_mechanism","target_missingness"),"flag","S2_Dip_detection_rate",N_ACTIVE_REP),
  summary_binary(s2diag[,flag2:=sigclust_standardized_effect_sd>=SIGCLUST_EFFECT_CUTOFF],c("scenario","missingness_mechanism","target_missingness"),"flag2","S2_standardized_effect_detection_rate",N_ACTIVE_REP),
  summary_binary(d1[scenario=="S2_true_discrete_subtypes"&k==3L][,flag:=ari_vs_latent_subtype>=TRUE_K3_ARI_CUTOFF],c("scenario","missingness_mechanism","target_missingness","k"),"flag","S2_K3_true_structure_recovery_rate",N_ACTIVE_REP)
))

d2_summary <- summary_numeric(d2rep,c("scenario","missingness_mechanism","target_missingness","outcome","label_source","contrast"),"pooled_delta_auc","Rubin_pooled_delta_AUC")
d2_secondary <- safe_rbind(list(
  summary_numeric(d2rep,c("scenario","missingness_mechanism","target_missingness","outcome","label_source","contrast"),"pooled_continuous_nri","Continuous_NRI_secondary"),
  summary_numeric(d2rep,c("scenario","missingness_mechanism","target_missingness","outcome","label_source","contrast"),"pooled_idi","IDI_secondary")
))
decisive <- d2rep[contrast=="Raw-EN+label_minus_Raw-EN"]
d2_oc <- safe_rbind(list(
  summary_binary(decisive[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")][,flag:=ci_positive],c("scenario","missingness_mechanism","target_missingness","outcome"),"flag","Null_control_positive_CI_false_positive_rate",N_ACTIVE_REP),
  summary_binary(decisive[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")][,flag2:=practical_null],c("scenario","missingness_mechanism","target_missingness","outcome"),"flag2","Null_control_practical_null_rate",N_ACTIVE_REP),
  summary_binary(decisive[scenario=="S5_domain2_positive_control"][,flag:=ci_positive],c("scenario","missingness_mechanism","target_missingness","outcome"),"flag","S5_positive_control_detection_rate",N_ACTIVE_REP)
))
d2_failure <- d2rep[,.(expected_repeats=N_ACTIVE_REP,observed_repeat_metrics=sum(is.finite(pooled_delta_auc)),missing_repeat_metrics=N_ACTIVE_REP-sum(is.finite(pooled_delta_auc)),metric_failure_rate=(N_ACTIVE_REP-sum(is.finite(pooled_delta_auc)))/N_ACTIVE_REP),by=.(scenario,missingness_mechanism,target_missingness,outcome,label_source,contrast)]

## ---- complete-data contrasts -------------------------------------------------
ref_d1 <- d1[missingness_mechanism=="MCAR"&target_missingness==0,.(scenario,repeat_id,k,complete_ari=ari_vs_latent_subtype,complete_mortality_sep=mortality_30d_C1_minus_C2,complete_mice_ari=mice_ari_mean_vs_primary,complete_sig_effect=sigclust_standardized_effect_sd)]
d1_vs_complete <- merge(d1,ref_d1,by=c("scenario","repeat_id","k"),all.x=TRUE)
d1_vs_complete[,`:=`(ari_change_vs_complete=ari_vs_latent_subtype-complete_ari,mortality_separation_change_vs_complete=mortality_30d_C1_minus_C2-complete_mortality_sep,mice_ari_change_vs_complete=mice_ari_mean_vs_primary-complete_mice_ari,sigclust_effect_change_vs_complete=sigclust_standardized_effect_sd-complete_sig_effect)]
ref_d2 <- d2rep[missingness_mechanism=="MCAR"&target_missingness==0,.(scenario,repeat_id,outcome,label_source,contrast,complete_delta_auc=pooled_delta_auc,complete_nri=pooled_continuous_nri,complete_idi=pooled_idi)]
d2_vs_complete <- merge(d2rep,ref_d2,by=c("scenario","repeat_id","outcome","label_source","contrast"),all.x=TRUE)
d2_vs_complete[,`:=`(delta_auc_change_vs_complete=pooled_delta_auc-complete_delta_auc,continuous_nri_change_vs_complete=pooled_continuous_nri-complete_nri,idi_change_vs_complete=pooled_idi-complete_idi)]
ref_mv <- model_variability[missingness_mechanism=="MCAR"&target_missingness==0,.(scenario,repeat_id,outcome,label_source,contrast,complete_delta_auc_sd=delta_auc_across_imputation_sd)]
mv_vs_complete <- merge(model_variability,ref_mv,by=c("scenario","repeat_id","outcome","label_source","contrast"),all.x=TRUE)
mv_vs_complete[,delta_auc_variability_change_vs_complete:=delta_auc_across_imputation_sd-complete_delta_auc_sd]

## ---- paired mechanism and intensity contrasts -------------------------------
mechanism_contrast <- function(dt,value,id,metric){w<-dcast(dt,as.formula(paste(paste(id,collapse=" + "),"~ missingness_mechanism")),value.var=value);if(!all(MECHANISMS%in%names(w)))return(data.table());w[,`:=`(metric=metric,MAR_minus_MCAR=MAR-MCAR,MNAR_minus_MCAR=MNAR-MCAR,MNAR_minus_MAR=MNAR-MAR)];w}
intensity_contrast <- function(dt,value,id,metric){w<-dcast(dt,as.formula(paste(paste(id,collapse=" + "),"~ target_missingness")),value.var=value);if(!all(c("0","0.15","0.3")%in%names(w)))return(data.table());setnames(w,c("0","0.15","0.3"),c("value_0","value_15","value_30"));w[,`:=`(metric=metric,change_15_minus_0=value_15-value_0,change_30_minus_0=value_30-value_0,change_30_minus_15=value_30-value_15)];w}
d1_mech <- safe_rbind(list(
  mechanism_contrast(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],"ari_vs_latent_subtype",c("scenario","target_missingness","repeat_id","k"),"ARI_vs_latent_subtype"),
  mechanism_contrast(unique(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes"),.(scenario,missingness_mechanism,target_missingness,repeat_id,sigclust_standardized_effect_sd)]),"sigclust_standardized_effect_sd",c("scenario","target_missingness","repeat_id"),"SigClust_like_standardized_effect")
))
d2_mech <- mechanism_contrast(decisive,"pooled_delta_auc",c("scenario","target_missingness","repeat_id","outcome","label_source","contrast"),"Rubin_pooled_delta_AUC")
d1_trend <- safe_rbind(list(
  intensity_contrast(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],"ari_vs_latent_subtype",c("scenario","missingness_mechanism","repeat_id","k"),"ARI_vs_latent_subtype"),
  intensity_contrast(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],"mice_ari_mean_vs_primary",c("scenario","missingness_mechanism","repeat_id","k"),"MICE_label_stability_ARI")
))
d2_trend <- intensity_contrast(decisive,"pooled_delta_auc",c("scenario","missingness_mechanism","repeat_id","outcome","label_source","contrast"),"Rubin_pooled_delta_AUC")


contrast_summary <- function(dt, group_cols, contrast_cols, summary_name) {
  if (!nrow(dt)) return(data.table())
  cols <- intersect(contrast_cols, names(dt))
  if (!length(cols)) return(data.table())
  long <- melt(dt, id.vars = group_cols, measure.vars = cols,
               variable.name = "contrast_name", value.name = "contrast_value")
  out <- summary_numeric(long, c(group_cols, "contrast_name"),
                         "contrast_value", summary_name)
  out
}

d1_mech_summary <- contrast_summary(
  d1_mech,
  intersect(c("scenario","target_missingness","k","metric"), names(d1_mech)),
  c("MAR_minus_MCAR","MNAR_minus_MCAR","MNAR_minus_MAR"),
  "Paired_missingness_mechanism_contrast"
)
d2_mech_summary <- contrast_summary(
  d2_mech,
  intersect(c("scenario","target_missingness","outcome","label_source","contrast","metric"), names(d2_mech)),
  c("MAR_minus_MCAR","MNAR_minus_MCAR","MNAR_minus_MAR"),
  "Paired_missingness_mechanism_contrast"
)
d1_trend_summary <- contrast_summary(
  d1_trend,
  intersect(c("scenario","missingness_mechanism","k","metric"), names(d1_trend)),
  c("change_15_minus_0","change_30_minus_0","change_30_minus_15"),
  "Paired_missingness_intensity_contrast"
)
d2_trend_summary <- contrast_summary(
  d2_trend,
  intersect(c("scenario","missingness_mechanism","outcome","label_source","contrast","metric"), names(d2_trend)),
  c("change_15_minus_0","change_30_minus_0","change_30_minus_15"),
  "Paired_missingness_intensity_contrast"
)

## ---- MC-SE and ADEMP ----------------------------------------------------------
mcse <- safe_rbind(list(
  summary_numeric(realized_overall,c("scenario","missingness_mechanism","target_missingness"),"realized_overall_missingness","Realized_overall_missingness"),
  summary_numeric(mice_stability,c("scenario","missingness_mechanism","target_missingness","k"),"ari_mean_vs_primary","MICE_label_stability_ARI"),
  summary_numeric(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes")],c("scenario","missingness_mechanism","target_missingness","k"),"ari_vs_latent_subtype","Domain1_ARI_vs_truth"),
  summary_numeric(unique(d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes"),.(scenario,missingness_mechanism,target_missingness,repeat_id,sigclust_standardized_effect_sd)]),c("scenario","missingness_mechanism","target_missingness"),"sigclust_standardized_effect_sd","SigClust_like_standardized_effect"),
  summary_numeric(decisive,c("scenario","missingness_mechanism","target_missingness","outcome","label_source"),"pooled_delta_auc","Decisive_Rubin_pooled_delta_AUC"),
  summary_numeric(model_variability,c("scenario","missingness_mechanism","target_missingness","outcome","label_source","contrast"),"delta_auc_across_imputation_sd","Across_imputation_delta_AUC_SD")
))
ademp <- data.table(component=c("Aim","Data-generating mechanisms","Missingness mechanisms","Intensities","Estimands","Methods","Performance measures","Repetitions and MC error","Realized-missingness fidelity","Interpretation boundary"),specification=c(
  "Assess robustness of Domain 1 and Domain 2 conclusions to missingness mechanism and intensity.",
  "S1 continuum; S2 true K=3; S5 hidden state plus fixed 90%-accurate oracle label.",
  "MCAR constant probabilities; MAR depends on fully observed baseline covariates; MNAR depends on each feature's own true value plus latent severity/state; variable/cohort intercepts are calibrated.",
  "Target feature-cell missingness 0%, 15%, and 30%.",
  "Realized missingness; MICE stability; S1 false positives; S2 K3 recovery; Raw-EN+label minus Raw-EN Delta AUC; paired complete-data and mechanism contrasts.",
  "Train/test imputed separately with MICE m=5 maxit=5; post-imputation 1st/99th winsorization; training-fitted scaling and centroids; held-out elastic-net prediction; Rubin Delta AUC pooling.",
  "Means, empirical intervals, MC-SE, exact-binomial operating characteristics, task/model failure, mechanism fidelity, and intensity trends.",
  paste0(N_ACTIVE_REP," repeats per setting in this run; all formal summaries use repeat-level results."),
  "Overall/variable realized rates, target error, mechanism-signal AUC/correlation, mechanism balance, zero equivalence, and monotonicity are explicit QC outputs.",
  "MICE stability is reproducibility evidence, not discreteness; PAC, Gap, transportability, Domain 4, and alternative algorithms are excluded."
))

## =============================================================================
## 12. Critical QC, provenance, outputs, and figures
## =============================================================================

critical_qc <- data.table(
  qc_item=c(
    "All expected tasks completed","No duplicate repeat records","Target-zero missingness is exactly zero",
    "Overall realized missingness matches target","Realized rates are balanced across mechanisms",
    "Missingness increases monotonically","Mechanism implementation fidelity",
    "MICE completed with no residual missingness","Zero-setting analyses are identical across mechanisms",
    "Critical result tables are nonempty"
  ),
  pass=c(
    isTRUE(all(completion$structural_qc_pass)),nrow(dup)==0L,
    isTRUE(all(realized_overall[target_missingness==0,realized_overall_missingness]==0)),
    isTRUE(all(setting_fidelity$overall_target_fidelity_pass)),isTRUE(all(balance_summary$mechanism_balance_pass)),
    isTRUE(all(mon_summary$all_repeats_monotonic)),isTRUE(all(setting_fidelity$fidelity_pass)),
    isTRUE(all(mice_completion$residual_missing_cells==0 & mice_completion$completed_imputations==mice_completion$expected_imputations)),
    isTRUE(all(zero_eq$equivalence_pass)) && isTRUE(all(zero_eq_d2$equivalence_pass)),
    all(c(nrow(realized_overall)>0,nrow(realized_variable)>0,nrow(d1)>0,nrow(d2rep)>0,nrow(mice_stability)>0))
  ),
  critical=TRUE,
  note=c(
    paste0(length(success),"/",nrow(tasks)," successful tasks"),paste0(nrow(dup)," duplicate keys"),
    "All target-zero masks must be empty",
    paste0("mean abs error <= ",OVERALL_ABS_ERROR_TOL,"; q97.5 <= ",OVERALL_Q975_ABS_ERROR_TOL),
    paste0("mean range <= ",MECHANISM_RANGE_TOL,"; q97.5 <= ",MECHANISM_Q975_RANGE_TOL),
    "0% <= 15% <= 30% in every paired repeat",
    paste0("MAR/MNAR AUC >= ",MAR_SIGNAL_AUC_MIN,"; MCAR mean abs correlation <= ",MCAR_ABS_COR_MAX),
    "Expected completed imputations and zero residual cells",paste0("tolerance ",ZERO_EQUIVALENCE_TOL),
    "Realized missingness, MICE, Domain 1, and Domain 2 outputs"
  )
)
critical_qc[is.na(pass), pass := FALSE]
if(USE_TEST_MODE) critical_qc[qc_item%in%c("Overall realized missingness matches target","Realized rates are balanced across mechanisms","Mechanism implementation fidelity"),critical:=FALSE]
critical_failures <- critical_qc[critical==TRUE & pass==FALSE]

if (nrow(errors)) {
  errors[, error_stage := fcase(
    grepl("MICE", error, ignore.case = TRUE), "MICE",
    grepl("k-means|cluster", error, ignore.case = TRUE), "clustering",
    grepl("glmnet|AUC|prediction|model", error, ignore.case = TRUE), "prediction",
    default = "other"
  )]
  failure_summary <- errors[, .(
    count = .N,
    proportion_of_expected_repeats = .N / N_ACTIVE_REP
  ), by = .(
    scenario, missingness_mechanism, target_missingness,
    error_stage, error
  )]
} else {
  errors[, error_stage := character()]
  failure_summary <- data.table()
}

mice_task_failure <- if (nrow(errors)) {
  errors[, .(
    expected_repeats = N_ACTIVE_REP,
    mice_failed_repeats = sum(error_stage == "MICE"),
    mice_task_failure_rate = sum(error_stage == "MICE") / N_ACTIVE_REP
  ), by = .(scenario, missingness_mechanism, target_missingness)]
} else {
  unique(tasks[, .(
    scenario, missingness_mechanism, target_missingness,
    expected_repeats = N_ACTIVE_REP,
    mice_failed_repeats = 0L,
    mice_task_failure_rate = 0
  )])
}
provenance <- data.table(
  output_group=c("realized_missingness","MICE","Domain1","Domain2","operating_characteristics","figures"),
  source_level="task/repeat-level rerun",
  derivation=c(
    "Calibrated feature-cell masks before imputation",
    "Separate train/test MICE; outcomes and latent fields excluded",
    "Training-set K2/K3 plus primary-imputation dip and SigClust-like effect",
    "Held-out Base/Raw-EN comparisons pooled across five imputations",
    "Expected-repeat denominators; no conditioning on successful metric estimation",
    "Generated directly from final CSV data"
  ),
  deprecated_input_used=FALSE,
  config_signature=CONFIG_SIGNATURE
)
metric_map <- data.table(
  metric=c("realized_overall_missingness","selected_signal_auc","ari_vs_latent_subtype","sigclust_standardized_effect_sd","mice_ari_mean_vs_primary","pooled_delta_auc","ci_positive","practical_null"),
  definition=c(
    "Missing feature cells divided by all eligible feature cells",
    "AUC of the prespecified mechanism score for the realized missingness indicator",
    "ARI between fitted labels and known latent class; S1 latent class is noninformative",
    "Null-mean minus observed clustering index divided by null SD",
    "Mean label ARI versus completed dataset 1",
    "Rubin-pooled paired AUC difference across completed datasets",
    "Lower bound of Rubin 95% CI is greater than zero",
    paste0("Absolute pooled Delta AUC < ",D2_PRACTICAL_NULL)
  ),
  primary_or_secondary=c("QC","QC","primary D1","primary auxiliary D1","reproducibility","primary D2","operating characteristic","operating characteristic")
)

## Save all tables.
tables <- list(
  "05_realized_missingness_by_repeat"=realized_overall,
  "05_realized_missingness_summary"=realized_summary,
  "05_variable_specific_missingness_by_repeat"=realized_variable,
  "05_variable_specific_missingness_summary"=variable_summary,
  "05_missingness_mechanism_fidelity_by_repeat"=mech_fidelity_rep,
  "05_missingness_mechanism_fidelity_summary"=mech_fidelity,
  "05_missingness_balance_across_mechanisms_by_repeat"=balance_rep,
  "05_missingness_balance_across_mechanisms"=balance_summary,
  "05_missingness_intensity_monotonicity_by_repeat"=mon,
  "05_missingness_intensity_monotonicity_summary"=mon_summary,
  "05_zero_missingness_domain1_equivalence_qc"=zero_eq,
  "05_zero_missingness_domain2_equivalence_qc"=zero_eq_d2,
  "05_mice_status_by_task"=mice_status,
  "05_mice_completion_summary"=mice_completion,
  "05_mice_task_failure_summary"=mice_task_failure,
  "05_mice_label_stability_by_repeat"=mice_stability,
  "05_mice_label_stability_summary"=mice_stability_summary,
  "05_mice_stability_contrast_vs_complete"=d1_vs_complete[,.(task_id,scenario,missingness_mechanism,target_missingness,repeat_id,k,mice_ari_mean_vs_primary,complete_mice_ari,mice_ari_change_vs_complete)],
  "05_mice_model_variability_by_repeat"=model_variability,
  "05_mice_model_variability_summary"=model_variability_summary,
  "05_mice_model_variability_contrast_vs_complete"=mv_vs_complete,
  "05_domain1_by_imputation"=d1imp,
  "05_domain1_repeat_level_standardized"=d1,
  "05_domain1_dip_axes"=dip_axes,
  "05_domain1_summary"=d1_summary,
  "05_domain1_operating_characteristics"=d1_oc,
  "05_domain1_contrast_vs_complete"=d1_vs_complete,
  "05_domain1_mechanism_contrasts"=d1_mech,
  "05_domain1_mechanism_contrast_summary"=d1_mech_summary,
  "05_domain1_intensity_trends"=d1_trend,
  "05_domain1_intensity_trend_summary"=d1_trend_summary,
  "05_domain2_by_imputation"=d2imp,
  "05_domain2_by_repeat_Rubin"=d2rep,
  "05_domain2_delta_auc_summary"=d2_summary,
  "05_domain2_operating_characteristics"=d2_oc,
  "05_domain2_contrast_vs_complete"=d2_vs_complete,
  "05_domain2_mechanism_contrasts"=d2_mech,
  "05_domain2_mechanism_contrast_summary"=d2_mech_summary,
  "05_domain2_intensity_trends"=d2_trend,
  "05_domain2_intensity_trend_summary"=d2_trend_summary,
  "05_domain2_nri_idi_secondary"=d2_secondary,
  "05_domain2_failure_summary"=d2_failure,
  "05_mcse_summary"=mcse,
  "05_ADEMP_table"=ademp,
  "05_task_failure_summary"=failure_summary,
  "05_result_provenance"=provenance,
  "05_metric_definition_map"=metric_map
)
for(nm in names(tables)) fwrite(tables[[nm]],file.path(TAB_DIR,paste0(nm,".csv")))
fwrite(setting_fidelity,file.path(QC_DIR,"05_realized_missingness_setting_fidelity_qc.csv"))
fwrite(critical_qc,file.path(QC_DIR,"05_critical_qc_summary.csv"))

## ---- publication-oriented figures -------------------------------------------
fd_miss <- realized_summary[,.(scenario,missingness_mechanism,target_missingness,realized_mean,realized_q025,realized_q975)]
fwrite(fd_miss,file.path(TAB_DIR,"Figure_S05A_realized_missingness_data.csv"))
p1 <- ggplot(fd_miss,aes(target_missingness,realized_mean,shape=missingness_mechanism))+geom_abline(intercept=0,slope=1,linetype=2,linewidth=.4)+geom_point(position=position_dodge(.015),size=2)+geom_errorbar(aes(ymin=realized_q025,ymax=realized_q975),width=.008,position=position_dodge(.015))+facet_wrap(~scenario)+labs(x="Target feature-cell missingness",y="Realized feature-cell missingness",shape="Mechanism",title="Realized-missingness fidelity")+theme_publication()
save_plot_both(p1,"Figure_S05A_realized_missingness_fidelity",10,5)

fd_d1 <- d1[scenario%in%c("S1_pure_severity_continuum","S2_true_discrete_subtypes"),.(mean=safe_mean(ari_vs_latent_subtype),low=safe_q(ari_vs_latent_subtype,.025),high=safe_q(ari_vs_latent_subtype,.975)),by=.(scenario,missingness_mechanism,target_missingness,k)]
fwrite(fd_d1,file.path(TAB_DIR,"Figure_S05B_domain1_recovery_data.csv"))
p2 <- ggplot(fd_d1,aes(factor(target_missingness),mean,shape=missingness_mechanism))+geom_point(position=position_dodge(.35),size=2)+geom_errorbar(aes(ymin=low,ymax=high),width=.14,position=position_dodge(.35))+facet_grid(k~scenario,labeller=label_both)+coord_cartesian(ylim=c(-.05,1))+labs(x="Target missingness",y="ARI versus latent truth",shape="Mechanism",title="Domain 1 robustness to missingness")+theme_publication()
save_plot_both(p2,"Figure_S05B_domain1_missingness_robustness",10,6)

fd_mice <- mice_stability[,.(mean=safe_mean(ari_mean_vs_primary),low=safe_q(ari_mean_vs_primary,.025),high=safe_q(ari_mean_vs_primary,.975)),by=.(scenario,missingness_mechanism,target_missingness,k)]
fwrite(fd_mice,file.path(TAB_DIR,"Figure_S05C_mice_stability_data.csv"))
p3 <- ggplot(fd_mice,aes(factor(target_missingness),mean,shape=missingness_mechanism))+geom_point(position=position_dodge(.35),size=2)+geom_errorbar(aes(ymin=low,ymax=high),width=.14,position=position_dodge(.35))+facet_grid(k~scenario,labeller=label_both)+coord_cartesian(ylim=c(0,1))+labs(x="Target missingness",y="ARI versus completed dataset 1",shape="Mechanism",title="Cross-imputation label stability")+theme_publication()
save_plot_both(p3,"Figure_S05C_mice_label_stability",11,6)

fd_d2 <- decisive[,.(mean=safe_mean(pooled_delta_auc),low=safe_q(pooled_delta_auc,.025),high=safe_q(pooled_delta_auc,.975)),by=.(scenario,missingness_mechanism,target_missingness,outcome,label_source)]
fwrite(fd_d2,file.path(TAB_DIR,"Figure_S05D_domain2_delta_auc_data.csv"))
p4 <- ggplot(fd_d2,aes(factor(target_missingness),mean,shape=missingness_mechanism))+geom_hline(yintercept=0,linewidth=.4)+geom_point(position=position_dodge(.35),size=2)+geom_errorbar(aes(ymin=low,ymax=high),width=.14,position=position_dodge(.35))+facet_grid(outcome~scenario,scales="free_y")+labs(x="Target missingness",y="Raw-EN+label minus Raw-EN Delta AUC",shape="Mechanism",title="Domain 2 robustness to missingness")+theme_publication()
save_plot_both(p4,"Figure_S05D_domain2_missingness_robustness",12,6)

## ---- completion markers ------------------------------------------------------
all_tasks_successful <- length(success)==nrow(tasks)&&nrow(errors)==0L
all_critical_qc_pass <- nrow(critical_failures)==0L
completion_summary <- data.table(analysis_id=ANALYSIS_ID,run_mode=RUN_MODE,expected_tasks=nrow(tasks),successful_tasks=length(success),failed_tasks=nrow(errors),all_tasks_successful=all_tasks_successful,all_critical_qc_pass=all_critical_qc_pass,completed_at=timestamp(),config_signature=CONFIG_SIGNATURE)
fwrite(completion_summary,file.path(LOG_DIR,"sensitivity05_run_completion_summary.csv"))
completed_marker<-file.path(LOG_DIR,"run_completed.ok");partial_marker<-file.path(LOG_DIR,"run_partial_or_qc_failure.txt")
if(file.exists(completed_marker))file.remove(completed_marker);if(file.exists(partial_marker))file.remove(partial_marker)
if(all_tasks_successful&&all_critical_qc_pass){writeLines(paste0("completed_at=",timestamp(),"\nrun_mode=",RUN_MODE,"\nexpected_tasks=",nrow(tasks),"\nsuccessful_tasks=",length(success),"\nfailed_tasks=0\nall_critical_qc_pass=TRUE\nconfig_signature=",CONFIG_SIGNATURE),completed_marker)}else if(ALLOW_PARTIAL_OUTPUTS){writeLines(paste0("partial_at=",timestamp(),"\nrun_mode=",RUN_MODE,"\nexpected_tasks=",nrow(tasks),"\nsuccessful_tasks=",length(success),"\nfailed_tasks=",nrow(errors),"\nall_critical_qc_pass=",all_critical_qc_pass,"\nReview QC and rerun failed checkpoints.\nconfig_signature=",CONFIG_SIGNATURE),partial_marker)}
write_session_info(
  file.path(LOG_DIR, "sessionInfo.txt")
)

warning_text <- capture.output(warnings())
writeLines(
  text = warning_text,
  con = file.path(LOG_DIR, "warnings_after_run.txt"),
  useBytes = TRUE
)
log_msg("SA-05 finished. all tasks successful: ",all_tasks_successful,"; all critical QC passed: ",all_critical_qc_pass)
cat("\n============================================================\nSensitivity 05 finished.\nRun mode:",RUN_MODE,"\nExpected tasks:",nrow(tasks),"\nSuccessful tasks:",length(success),"\nFailed tasks:",nrow(errors),"\nAll critical QC passed:",all_critical_qc_pass,"\nOutput:\n",OUT_DIR,"\n============================================================\n")
if(STOP_ON_CRITICAL_QC&&nrow(critical_failures)>0L)stop("Sensitivity 05 completed but critical QC failed; review 05_critical_qc_summary.csv.",call.=FALSE)
