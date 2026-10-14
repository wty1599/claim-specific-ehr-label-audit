args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(args)) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(dirname(dirname(sub("^--file=", "", script_arg[[1]]))), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

r_files <- list.files(
  repo_root,
  pattern = "\\.[Rr]$",
  recursive = TRUE,
  full.names = TRUE
)

results <- lapply(r_files, function(path) {
  err <- NA_character_
  ok <- tryCatch(
    {
      parse(file = path, encoding = "UTF-8")
      TRUE
    },
    error = function(e) {
      err <<- conditionMessage(e)
      FALSE
    }
  )

  path_norm <- gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = TRUE))
  relative_path <- substring(path_norm, nchar(repo_root) + 2L)

  data.frame(
    file = relative_path,
    parse_ok = ok,
    error = err,
    stringsAsFactors = FALSE
  )
})

results <- do.call(rbind, results)

if (any(!results$parse_ok)) {
  print(results[!results$parse_ok, , drop = FALSE])
  stop("One or more R scripts failed static parsing.", call. = FALSE)
}

cat("All", nrow(results), "R scripts parsed successfully.\n")
