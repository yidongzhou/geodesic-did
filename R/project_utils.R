# Utilities shared by the command-line reproduction scripts.

script_file <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(argument)) return(NULL)
  path <- sub("^--file=", "", argument[1])
  normalizePath(gsub("~+~", " ", path, fixed = TRUE), mustWork = TRUE)
}

project_root <- function(script = script_file()) {
  if (is.null(script)) return(normalizePath(".", mustWork = TRUE))
  normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
}

load_gdid_code <- function(root = project_root()) {
  source(file.path(root, "R", "lcm.R"), local = .GlobalEnv)
  source(file.path(root, "R", "gdd.R"), local = .GlobalEnv)
  invisible(root)
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf("Install the following R packages before running this script: %s",
              paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

load_single_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!expected_name %in% loaded) {
    stop(sprintf("%s does not contain object '%s'.", path, expected_name), call. = FALSE)
  }
  environment[[expected_name]]
}
