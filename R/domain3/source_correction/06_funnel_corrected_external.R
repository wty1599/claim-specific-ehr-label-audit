args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/source_correction/00_correction_paths.R"))
if (!length(args)) {
  args <- c(file.path(private, "corrected_hospital_counts_local.csv"),
            file.path(aggregate, "funnel"), file.path(spec_root, "FUNNEL_POSTHOC_SPEC.md"))
}
if (length(args) != 3L) {
  stop("Usage: Rscript 06_funnel_corrected_external.R [restricted_csv output_dir public_spec]")
}
input <- normalizePath(args[[1]], mustWork = TRUE)
ehr_audit_assert_external_output(dirname(args[[1]]), repo)
ehr_audit_assert_external_output(dirname(args[[2]]), repo)
dir.create(args[[2]], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[[2]], mustWork = TRUE)
spec <- normalizePath(args[[3]], mustWork = TRUE)
if (!requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("lme4", quietly = TRUE)) {
  stop("digest and lme4 are required")
}

expected_input_sha <- "6b8ea20e9c5fc5136ca6d26cefba8f35af7b2d7e82cd5ae836990ce8c6940c28"
expected_spec_sha <- "87f5f92f4f9c750f7fa117e6d68f2987f3777013a5feab657006485182ca045f"
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
stopifnot(identical(sha(input), expected_input_sha),
          identical(sha(spec), expected_spec_sha))

raw <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("hospitalid", "stays", "corrected_c1") %in% names(raw)),
          nrow(raw) == 199L, !anyNA(raw[c("hospitalid", "stays", "corrected_c1")]),
          !anyDuplicated(raw$hospitalid))
d <- raw[c("stays", "corrected_c1")]
names(d) <- c("n_stays", "l1_c1_n")
rm(raw)
stopifnot(all(d$n_stays > 0), all(d$n_stays == floor(d$n_stays)),
          all(d$l1_c1_n == floor(d$l1_c1_n)),
          all(d$l1_c1_n >= 0), all(d$l1_c1_n <= d$n_stays),
          sum(d$n_stays) == 17465L, sum(d$l1_c1_n) == 3319L)

p_mimic <- 3902 / 20049
p_eicu <- 3319 / 17465
d$p <- d$l1_c1_n / d$n_stays
d$point_eligible <- d$l1_c1_n >= 10L &
  (d$n_stays - d$l1_c1_n) >= 10L
d$drift_10pp <- abs(d$p - p_mimic) > 0.10
stopifnot(sum(d$point_eligible) + sum(!d$point_eligible) == 199L)

size_breaks <- c(0, 4, 9, 19, 39, 79, Inf)
size_labels <- c("1-4", "5-9", "10-19", "20-39", "40-79", "80+")
d$size_band <- cut(d$n_stays, breaks = size_breaks, labels = size_labels)
size <- do.call(rbind, lapply(size_labels, function(label) {
  x <- d[d$size_band == label, ]
  data.frame(size_band = label, hospitals = nrow(x), stays = sum(x$n_stays),
             c1_stays = sum(x$l1_c1_n), individually_displayed = sum(x$point_eligible),
             suppressed = sum(!x$point_eligible),
             outside_10pp = sum(x$drift_10pp))
}))
public_size_bands <- c(size_labels[1:4], "40+")
size_public <- do.call(rbind, lapply(public_size_bands, function(label) {
  x <- if (label == "40+") d[d$n_stays >= 40L, ] else
    d[d$size_band == label, ]
  data.frame(size_band = label, hospitals = nrow(x), stays = sum(x$n_stays),
             c1_stays = sum(x$l1_c1_n),
             individually_displayed = sum(x$point_eligible),
             suppressed = sum(!x$point_eligible),
             outside_10pp = sum(x$drift_10pp))
}))
stopifnot(sum(size_public$hospitals) == 199L,
          all(size_public$suppressed == 0L |
                size_public$suppressed >= 2L))

q <- as.numeric(quantile(d$n_stays, probs = c(0, 0.25, 0.5, 0.75, 1),
                         names = FALSE))
overall <- data.frame(hospitals = nrow(d), total_stays = sum(d$n_stays),
                      total_c1 = sum(d$l1_c1_n),
                      n_min = q[1], n_q25 = q[2], n_median = q[3],
                      n_q75 = q[4], n_max = q[5],
                      individually_displayed = sum(d$point_eligible),
                      suppressed = sum(!d$point_eligible),
                      outside_10pp_total = sum(d$drift_10pp),
                      outside_10pp_displayed = sum(d$drift_10pp & d$point_eligible),
                      outside_10pp_suppressed = sum(d$drift_10pp & !d$point_eligible))

limits <- expand.grid(center = c("MIMIC-IV", "eICU"),
                      nominal_coverage = c(0.95, 0.998),
                      stringsAsFactors = FALSE)
funnel <- do.call(rbind, lapply(seq_len(nrow(limits)), function(i) {
  center <- limits$center[i]
  coverage <- limits$nominal_coverage[i]
  p <- if (center == "MIMIC-IV") p_mimic else p_eicu
  alpha <- 1 - coverage
  low <- qbinom(alpha / 2, size = d$n_stays, prob = p)
  high <- qbinom(1 - alpha / 2, size = d$n_stays, prob = p)
  below <- d$l1_c1_n < low
  above <- d$l1_c1_n > high
  do.call(rbind, lapply(c("all", "displayed", "suppressed"), function(group) {
    take <- switch(group, all = rep(TRUE, nrow(d)),
                   displayed = d$point_eligible,
                   suppressed = !d$point_eligible)
    data.frame(center = center, p = p, coverage = coverage,
               display_group = group, hospitals = sum(take),
               below = sum(below & take), above = sum(above & take),
               outside = sum((below | above) & take))
  }))
}))

model_status <- "fit_failed"
model_message <- NA_character_
sd_point <- NA_real_
sd_low <- NA_real_
sd_high <- NA_real_
fit_warnings <- character()
d$hospital <- factor(seq_len(nrow(d)))
fit <- tryCatch(withCallingHandlers(
  lme4::glmer(cbind(l1_c1_n, n_stays - l1_c1_n) ~ 1 + (1 | hospital),
              data = d, family = binomial(), nAGQ = 9L,
              control = lme4::glmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 200000L))),
  warning = function(w) fit_warnings <<- c(fit_warnings, conditionMessage(w))),
  error = function(e) { model_message <<- conditionMessage(e); NULL })

if (!is.null(fit)) {
  conv <- fit@optinfo$conv
  bad_conv <- (!is.null(conv$opt) && conv$opt != 0L) ||
    length(conv$lme4$messages) > 0L ||
    any(grepl("converg|singular|unidentif", fit_warnings, ignore.case = TRUE))
  if (bad_conv || lme4::isSingular(fit)) {
    model_status <- "nonconverged_or_singular"
    model_message <- paste(c(fit_warnings, conv$lme4$messages), collapse = " | ")
  } else {
    sd_point <- as.numeric(attr(lme4::VarCorr(fit)$hospital, "stddev"))
    ci <- tryCatch(suppressMessages(
      confint(fit, parm = "theta_", method = "profile")),
      error = function(e) { model_message <<- conditionMessage(e); NULL })
    if (!is.null(ci) && nrow(ci) == 1L &&
        all(is.finite(ci[1, ])) && ci[1, 1] <= sd_point && sd_point <= ci[1, 2]) {
      sd_low <- ci[1, 1]
      sd_high <- ci[1, 2]
      model_status <- "converged"
    } else {
      model_status <- "profile_ci_failed"
    }
  }
}

glmm <- data.frame(model = "C1 ~ 1 + (1 | hospital)",
                   n_hospitals = nrow(d), sd_logit = sd_point,
                   ci95_low = sd_low, ci95_high = sd_high,
                   status = model_status, message = model_message,
                   lme4_version = as.character(packageVersion("lme4")))
primary <- subset(funnel, center == "MIMIC-IV" & coverage == 0.998 &
                    display_group == "all")
editorial_stop <- primary$outside <= 2L || model_status != "converged"

write.csv(overall, file.path(output, "hospital_size_overall.csv"), row.names = FALSE)
dir.create(file.path(output, "private_local"), showWarnings = FALSE)
write.csv(size, file.path(output, "private_local", "hospital_size_bands_restricted.csv"),
          row.names = FALSE)
write.csv(size_public, file.path(output, "hospital_size_bands_public.csv"),
          row.names = FALSE)
write.csv(funnel, file.path(output, "exact_binomial_funnel_counts.csv"), row.names = FALSE)
write.csv(glmm, file.path(output, "random_intercept_model.csv"), row.names = FALSE)
writeLines(c(paste0("spec_sha256=", sha(spec)),
             paste0("private_input_sha256=", sha(input)),
             paste0("editorial_stop=", editorial_stop),
             paste0("primary_99_8_outside=", primary$outside),
             paste0("model_status=", model_status)),
           file.path(output, "analysis_status.txt"))
cat("Primary 99.8% outside:", primary$outside, "\n")
cat("Random-intercept model:", model_status, "\n")
cat("Editorial stop:", editorial_stop, "\n")
