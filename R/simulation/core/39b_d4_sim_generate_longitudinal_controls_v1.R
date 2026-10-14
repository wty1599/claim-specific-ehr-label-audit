if (!exists("SIM_ROOT")) source(file.path(
  Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()), "R", "simulation", "core",
  "39_d4_sim_common_v1.R"
))

simulate_d4_longitudinal_cohort <- function(scenario, repeat_id, n = SIM_N) {
  stopifnot(scenario %in% unname(SIM_SCENARIOS), n >= 500L)
  set.seed(sim_seed(scenario, repeat_id))

  severity <- rnorm(n)
  renal <- 0.55 * severity + sqrt(1 - 0.55^2) * rnorm(n)
  shock <- 0.50 * severity + sqrt(1 - 0.50^2) * rnorm(n)
  partition_score <- 0.85 * renal + 0.55 * shock + 0.25 * rnorm(n)
  cluster_k2 <- ifelse(
    partition_score > quantile(partition_score, 0.80), 1L, 2L
  )

  # The label is fixed before the landmark and is never refitted in D4.
  age <- sim_clamp(round(64 + 5 * severity + rnorm(n, 0, 12)), 18, 95)
  gender <- rbinom(n, 1, 0.55)
  sofa_score <- sim_clamp(round(7 + 2.5 * severity + rnorm(n, 0, 2)), 0, 24)
  sapsii <- sim_clamp(round(38 + 11 * severity + rnorm(n, 0, 8)), 0, 120)
  creatinine_max <- pmax(0.2, exp(log(1.2) + 0.40 * renal + 0.12 * severity))
  bun_max <- pmax(2, exp(log(24) + 0.35 * renal + 0.12 * severity))
  potassium_max <- sim_clamp(4.2 + 0.28 * renal + rnorm(n, 0, 0.35), 2, 8)
  bicarbonate_min <- sim_clamp(22 - 2.4 * severity - 1.2 * shock + rnorm(n, 0, 2), 5, 40)
  ph_min <- sim_clamp(7.36 - 0.045 * severity - 0.025 * shock + rnorm(n, 0, 0.035), 6.7, 7.7)
  lactate_max <- pmax(0.2, exp(log(1.8) + 0.35 * shock + 0.15 * severity))
  urine_output_24h_ml <- pmax(0, 1800 - 330 * renal - 180 * severity + rnorm(n, 0, 450))
  mbp_min <- sim_clamp(67 - 7 * shock - 3 * severity + rnorm(n, 0, 6), 20, 130)
  heart_rate_max <- sim_clamp(105 + 12 * severity + 8 * shock + rnorm(n, 0, 10), 30, 240)
  resp_rate_max <- sim_clamp(27 + 5 * severity + rnorm(n, 0, 4), 4, 70)
  gcs_min <- sim_clamp(round(14 - 1.5 * severity + rnorm(n, 0, 1.2)), 3, 15)
  aki_lp <- -0.5 + 0.9 * renal + 0.4 * severity
  aki_u <- runif(n)
  aki_stage_0_24h <- fifelse(
    aki_u < plogis(aki_lp - 1.2), 3L,
    fifelse(aki_u < plogis(aki_lp - 0.2), 2L,
            fifelse(aki_u < plogis(aki_lp + 0.8), 1L, 0L))
  )

  dt <- data.table(
    stay_id = sprintf("%s_%03d_%06d", scenario, repeat_id, seq_len(n)),
    cluster_k2, age, gender, sofa_score, aki_stage_0_24h, sapsii,
    creatinine_max, bun_max, potassium_max, bicarbonate_min, ph_min,
    lactate_max, urine_output_24h_ml, mbp_min, heart_rate_max,
    resp_rate_max, gcs_min,
    latent_severity = severity, latent_renal = renal, latent_shock = shock
  )

  # Prespecified pre-landmark exclusions. D4-LQ increases each rate.
  stress <- scenario == SIM_SCENARIOS[["D4LQ"]]
  p_preicu_rrt <- if (stress) 0.010 else 0.0002
  p_pre24_rrt <- if (stress) 0.080 else plogis(-4.10 + 0.35 * renal)
  p_pre24_death <- if (stress) 0.080 else plogis(-4.15 + 0.55 * severity)
  p_pre24_exit <- if (stress) 0.180 else plogis(-2.45 - 0.30 * severity)

  preicu_rrt <- rbinom(n, 1, p_preicu_rrt) == 1L
  pre24_rrt <- !preicu_rrt & rbinom(n, 1, p_pre24_rrt) == 1L
  pre24_death <- rbinom(n, 1, p_pre24_death) == 1L
  pre24_exit <- rbinom(n, 1, p_pre24_exit) == 1L

  dt[, first_rrt_hours := NA_real_]
  dt[preicu_rrt, first_rrt_hours := -runif(.N, 0.1, 48)]
  dt[pre24_rrt, first_rrt_hours := runif(.N, 0, 23.999)]
  dt[, death_time_hours := NA_real_]
  dt[pre24_death, death_time_hours := runif(.N, 0, 23.999)]
  dt[, outtime_hours := runif(.N, 72, 168)]
  dt[pre24_exit, outtime_hours := runif(.N, 2, 23.999)]

  eligible <- is.na(dt$first_rrt_hours) &
    (is.na(dt$death_time_hours) | dt$death_time_hours >= 24) &
    dt$outtime_hours >= 24

  # Center allocation drivers within the fixed phenotype so that the
  # phenotype-specific G4-L baselines retain their prespecified meaning.
  renal_centered <- renal - ave(renal, cluster_k2, FUN = mean)
  shock_centered <- shock - ave(shock, cluster_k2, FUN = mean)

  if (scenario == SIM_SCENARIOS[["G4L"]]) {
    base <- ifelse(dt$cluster_k2 == 1L, qlogis(0.100), qlogis(0.0045))
    p_rrt <- plogis(base + 0.65 * renal_centered + 0.25 * shock_centered)
  } else {
    p_rrt <- plogis(qlogis(0.22) + 0.18 * severity +
                      0.12 * renal + 0.10 * (dt$cluster_k2 == 1L))
    p_rrt <- if (stress) {
      sim_clamp(p_rrt, 0.12, 0.40)
    } else {
      sim_clamp(p_rrt, 0.18, 0.35)
    }
  }
  dt[, true_rrt_probability := p_rrt]

  # Generate independent candidate times for RRT, death, and ICU exit. The
  # first-event ledger below, rather than the generator, determines event order.
  idx <- which(eligible)
  candidate_wait <- function(p) {
    p <- sim_clamp(p, 1e-8, 1 - 1e-8)
    w <- rexp(length(p), rate = -log1p(-p) / 48)
    fifelse(w < 48, w, NA_real_)
  }
  p_death <- plogis(-3.75 + 0.55 * severity +
                      0.35 * (dt$cluster_k2 == 1L))
  p_exit <- plogis(-0.55 - 0.35 * severity)
  rrt_wait <- candidate_wait(p_rrt[idx])
  death_wait <- candidate_wait(p_death[idx])
  exit_wait <- candidate_wait(p_exit[idx])
  dt[idx[!is.na(rrt_wait)], first_rrt_hours := 24 + rrt_wait[!is.na(rrt_wait)]]
  dt[idx[!is.na(death_wait)], death_time_hours := 24 + death_wait[!is.na(death_wait)]]
  dt[idx[!is.na(exit_wait)], outtime_hours := 24 + exit_wait[!is.na(exit_wait)]]

  # D4-LQ deliberately creates exact RRT/death ties; the ledger must select RRT.
  if (stress) {
    tie_pool <- idx[!is.na(dt$first_rrt_hours[idx]) &
                       dt$first_rrt_hours[idx] >= 24 & dt$first_rrt_hours[idx] < 72]
    tie_n <- min(length(tie_pool), max(1L, floor(0.01 * length(idx))))
    if (tie_n > 0L) {
      tie_idx <- sample(tie_pool, tie_n)
      dt[tie_idx, death_time_hours := first_rrt_hours]
    }
  }

  mortality_lp <- -2.05 + 0.65 * severity + 0.55 * (cluster_k2 == 1L) +
    0.015 * (age - 65) + 0.035 * (sofa_score - 7)
  dt[, mortality_30d := rbinom(.N, 1, plogis(mortality_lp))]
  dt[!is.na(death_time_hours) & death_time_hours < 720, mortality_30d := 1L]

  # EHR-like MAR missingness in a subset of PS covariates. The audit engine
  # applies the same median/mode handling as the empirical D4 analysis.
  miss_vars <- c(
    "sapsii", "creatinine_max", "bun_max", "bicarbonate_min", "ph_min",
    "lactate_max", "urine_output_24h_ml", "mbp_min", "gcs_min"
  )
  for (j in seq_along(miss_vars)) {
    v <- miss_vars[j]
    p_miss <- plogis(qlogis(0.06) + 0.18 * severity + 0.05 * j / length(miss_vars))
    z <- runif(n) < p_miss
    set(dt, which(z), v, NA)
  }
  dt[]
}
