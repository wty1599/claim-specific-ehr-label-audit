## ============================================================
## Shared clustering utilities used by the preprocessing and clustering scripts.
## ============================================================

need_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install required packages: install.packages(c(",
                         paste0("'", miss, "'", collapse = ", "), "))", call. = FALSE)
  invisible(TRUE)
}
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")


CKPT_DIR <- if (exists("DIR_OUTPUT")) file.path(DIR_OUTPUT, "checkpoint") else "checkpoint"
dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
save_ckpt <- function(obj, file) saveRDS(obj, file)
load_ckpt <- function(file) if (file.exists(file)) readRDS(file) else NULL


winsor_fit   <- function(x, p = c(0.01, 0.99)) as.numeric(quantile(x, p, na.rm = TRUE, type = 7))
winsor_apply <- function(x, b) pmin(pmax(x, b[1]), b[2])



.pr <- function(v, lo, hi) data.frame(var = v, lo = lo, hi = hi, stringsAsFactors = FALSE)
PHYS_RANGE <- do.call(rbind, list(
  .pr("wbc_max", 0, 300),            .pr("wbc_min", 0, 300),
  .pr("platelets_min", 0, 2000),     .pr("hemoglobin_min", 2, 25),
  .pr("abs_neutrophils_max", 0, 250),.pr("abs_lymphocytes_min", 0, 100),
  .pr("bands_max", 0, 100),          .pr("bilirubin_total_max", 0, 80),
  .pr("alt_max", 0, 20000),          .pr("ast_max", 0, 20000),
  .pr("albumin_min", 0.5, 7),        .pr("inr_max", 0.5, 30),
  .pr("pt_max", 5, 200),             .pr("ptt_max", 10, 300),
  .pr("fibrinogen_min", 0, 1500),    .pr("creatinine_max", 0, 30),
  .pr("creatinine_min", 0, 30),      .pr("bun_max", 0, 300),
  .pr("aniongap_max", 0, 60),        .pr("bicarbonate_min", 0, 60),
  .pr("glucose_max", 10, 2000),      .pr("sodium_min", 90, 200),
  .pr("sodium_max", 90, 200),        .pr("potassium_min", 1, 12),
  .pr("potassium_max", 1, 12),       .pr("calcium_min", 2, 20),
  .pr("calcium_max", 2, 20),         .pr("chloride_min", 50, 160),
  .pr("chloride_max", 50, 160),      .pr("lactate_max", 0, 40),
  .pr("ph_min", 6.5, 7.9),           .pr("po2_min", 10, 700),
  .pr("pco2_max", 5, 250),           .pr("so2_min", 0, 100),
  .pr("pao2fio2ratio_min", 0, 700),  .pr("baseexcess_min", -40, 40),
  .pr("heart_rate_max", 0, 350),     .pr("sbp_min", 0, 300),
  .pr("mbp_min", 0, 250),            .pr("mbp_mean", 0, 250),
  .pr("resp_rate_max", 0, 90),       .pr("temperature_min", 25, 45),
  .pr("temperature_max", 25, 45),    .pr("spo2_min", 0, 100),
  .pr("gcs_min", 3, 15),             .pr("urine_output_24h_ml", 0, 20000)
))


apply_phys_ranges <- function(df, vars, range_tbl = PHYS_RANGE) {
  rep <- data.frame(var = character(), n_out = integer(), pct_out = numeric())
  for (v in intersect(vars, range_tbl$var)) {
    lo <- range_tbl$lo[range_tbl$var == v]; hi <- range_tbl$hi[range_tbl$var == v]
    x  <- df[[v]]; bad <- !is.na(x) & (x < lo | x > hi)
    if (any(bad)) df[[v]][bad] <- NA
    rep <- rbind(rep, data.frame(var = v, n_out = sum(bad), pct_out = round(100 * mean(bad), 3)))
  }
  miss <- setdiff(vars, range_tbl$var)
  if (length(miss)) log_msg("Features without physiological bounds (cleaning skipped): ", paste(miss, collapse = ", "))
  list(df = df, report = rep[order(-rep$pct_out), ])
}
