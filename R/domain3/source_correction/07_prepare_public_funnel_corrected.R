args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- Sys.getenv("EHR_AUDIT_REPO_ROOT", unset =
  normalizePath(file.path(dirname(script_path), "../../.."), winslash = "/"))
source(file.path(repo, "R/domain3/source_correction/00_correction_paths.R"))
if (!length(args)) {
  args <- c(file.path(private, "corrected_hospital_counts_local.csv"),
            file.path(aggregate, "funnel"), figure_source,
            file.path(spec_root, "FUNNEL_POSTHOC_SPEC.md"))
}
if (length(args) != 4L) {
  stop("Usage: Rscript 07_prepare_public_funnel_corrected.R [restricted_csv analysis_dir figure_source_dir public_spec]")
}
input <- normalizePath(args[[1]], mustWork = TRUE)
analysis <- normalizePath(args[[2]], mustWork = TRUE)
public_dir <- normalizePath(args[[3]], mustWork = TRUE)
ehr_audit_assert_external_output(public_dir, repo)
spec <- normalizePath(args[[4]], mustWork = TRUE)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
stopifnot(identical(sha(input),
  "6b8ea20e9c5fc5136ca6d26cefba8f35af7b2d7e82cd5ae836990ce8c6940c28"),
  identical(sha(spec),
  "87f5f92f4f9c750f7fa117e6d68f2987f3777013a5feab657006485182ca045f"))
status <- readLines(file.path(analysis, "analysis_status.txt"))
stopifnot("editorial_stop=FALSE" %in% status)

raw <- read.csv(input, stringsAsFactors = FALSE)
d <- raw[c("stays", "corrected_c1")]
names(d) <- c("n_stays", "l1_c1_n")
rm(raw)
stopifnot(nrow(d) == 199L, sum(d$n_stays) == 17465L,
          sum(d$l1_c1_n) == 3319L)
d$eligible <- d$l1_c1_n >= 10L & d$n_stays - d$l1_c1_n >= 10L
stopifnot(sum(d$eligible) == 96L)

ind <- d[d$eligible, ]
individual <- data.frame(point_type = "hospital", n_stays = ind$n_stays,
  c1_prevalence = ind$l1_c1_n / ind$n_stays,
  plot_x = ind$n_stays, hospitals_aggregated = 1L,
  n_band = "individual")

hidden <- d[!d$eligible, ]
breaks <- c(0, 4, 9, 19, 39, 79, Inf)
labels <- c("1-4", "5-9", "10-19", "20-39", "40-79", "80+")
hidden$band <- cut(hidden$n_stays, breaks = breaks, labels = labels)
safe_group <- function(x) nrow(x) >= 2L && sum(x$n_stays) >= 10L &&
  sum(x$l1_c1_n) >= 10L && sum(x$n_stays - x$l1_c1_n) >= 10L
groups <- list()
pending <- hidden[FALSE, ]
pending_labels <- character()
for (label in labels) {
  x <- hidden[hidden$band == label, ]
  if (!nrow(x)) next
  pending <- rbind(pending, x)
  pending_labels <- c(pending_labels, label)
  if (safe_group(pending)) {
    groups[[length(groups) + 1L]] <- list(data = pending, labels = pending_labels)
    pending <- hidden[FALSE, ]
    pending_labels <- character()
  }
}
if (nrow(pending)) {
  if (!length(groups)) stop("Cannot construct a safe aggregate")
  last <- length(groups)
  groups[[last]]$data <- rbind(groups[[last]]$data, pending)
  groups[[last]]$labels <- c(groups[[last]]$labels, pending_labels)
}
stopifnot(all(vapply(groups, function(g) safe_group(g$data), logical(1))),
          sum(vapply(groups, function(g) nrow(g$data), integer(1))) == 103L)
pooled <- do.call(rbind, lapply(groups, function(g) {
  x <- g$data
  data.frame(point_type = "pooled", n_stays = sum(x$n_stays),
    c1_prevalence = sum(x$l1_c1_n) / sum(x$n_stays),
    plot_x = median(x$n_stays), hospitals_aggregated = nrow(x),
    n_band = paste(g$labels, collapse = " + "))
}))
public <- rbind(individual, pooled)
stopifnot(nrow(public) == 96L + length(groups),
          sum(public$hospitals_aggregated) == 199L,
          sum(public$n_stays) == 17465L,
          abs(sum(public$n_stays * public$c1_prevalence) - 3319) < 1e-7,
          all(public$n_stays >= 10L),
          all(public$hospitals_aggregated[public$point_type == "pooled"] >= 2L),
          all(public$n_stays[public$point_type == "hospital"] >= 20L))

summary <- read.csv(file.path(analysis, "hospital_size_overall.csv"))
funnel <- read.csv(file.path(analysis, "exact_binomial_funnel_counts.csv"))
primary <- subset(funnel, center == "MIMIC-IV" & coverage == 0.998 &
                    display_group == "all")
stopifnot(nrow(primary) == 1L, primary$outside == 9L)
summary$pooled_markers <- nrow(pooled)
summary$outside_99_8_mimic <- primary$outside
write.csv(public, file.path(public_dir, "Figure4_panelB_public_funnel_points.csv"),
          row.names = FALSE)
write.csv(summary, file.path(public_dir, "Figure4_panelB_funnel_summary.csv"),
          row.names = FALSE)
cat("Public individual hospitals:", nrow(individual), "\n")
cat("Suppressed hospitals pooled in", nrow(pooled), "markers\n")
