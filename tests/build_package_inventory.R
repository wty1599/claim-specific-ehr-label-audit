#!/usr/bin/env Rscript
# Inventory literal package calls in the files actually included in this copy.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(file.path(dirname(sub("^--file=", "", args)), ".."),
                      winslash = "/")
files <- list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
packages <- list()
walk <- function(node, found = character()) {
  if (is.expression(node) || is.pairlist(node)) {
    parts <- as.list(node)
    for (i in seq_along(parts)) {
      found <- tryCatch(walk(parts[[i]], found), error = function(e) found)
    }
  } else if (is.call(node)) {
    fun <- node[[1L]]
    if (is.call(fun)) found <- walk(fun, found)
    if (is.symbol(fun)) {
      name <- as.character(fun)
      if (name %in% c("::", ":::")) {
        pkg <- node[[2L]]
        if (is.symbol(pkg) || is.character(pkg)) found <- c(found, as.character(pkg))
      }
      if (name %in% c("library", "require", "requireNamespace") && length(node) >= 2L) {
        pkg <- node[[2L]]
        character_only <- isTRUE(as.list(node)$character.only)
        if (is.character(pkg) ||
            (name != "requireNamespace" && is.symbol(pkg) && !character_only)) {
          found <- c(found, as.character(pkg))
        }
      }
    }
    parts <- as.list(node)[-1L]
    for (i in seq_along(parts)) {
      found <- tryCatch(walk(parts[[i]], found), error = function(e) found)
    }
  }
  found
}
for (file in files) {
  relative <- substring(normalizePath(file, winslash = "/"), nchar(root) + 2L)
  for (package in unique(walk(parse(file)))) {
    packages[[package]] <- c(packages[[package]], relative)
  }
}
names <- sort(names(packages))
inventory <- data.frame(
  package = names,
  n_scripts = vapply(packages[names], function(x) length(unique(x)), integer(1)),
  files = vapply(packages[names], function(x) paste(sort(unique(x)), collapse = "; "), character(1)),
  stringsAsFactors = FALSE
)
csv <- character()
connection <- textConnection("csv", "w", local = TRUE)
write.csv(inventory, connection, row.names = FALSE, na = "")
close(connection)
writeBin(charToRaw(paste0(paste(csv, collapse = "\n"), "\n")),
         file.path(root, "R_PACKAGE_REQUIREMENTS.csv"))
message("Inventoried ", nrow(inventory), " literal package dependencies across ",
        length(files), " R files. Dynamic package names require manual review.")
