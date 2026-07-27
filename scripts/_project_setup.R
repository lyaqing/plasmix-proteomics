# Shared project setup for scripts in scripts/.

release_root <- Sys.getenv("PLASMIX_RELEASE_ROOT", unset = ".")
PROJECT_ROOT <- normalizePath(release_root, mustWork = TRUE)
required_project_dirs <- c("data", "results", "figures", "tables", "cache", "utils", "scripts")
missing_project_dirs <- required_project_dirs[!dir.exists(file.path(PROJECT_ROOT, required_project_dirs))]
if (length(missing_project_dirs)) {
    stop("Run from the 01_integration project root. Missing directories: ", paste(missing_project_dirs, collapse = ", "), call. = FALSE)
}
setwd(PROJECT_ROOT)

upstream_root_env <- Sys.getenv("PLASMIX_UPSTREAM_DATA_ROOT", unset = "")
UPSTREAM_DATA_ROOT <- if (nzchar(upstream_root_env)) upstream_root_env else file.path(dirname(PROJECT_ROOT), "00_data")

upstream_path <- function(...) {
    if (!dir.exists(UPSTREAM_DATA_ROOT)) {
        stop("The upstream 00_data directory was not found. Set PLASMIX_UPSTREAM_DATA_ROOT.", call. = FALSE)
    }
    file.path(UPSTREAM_DATA_ROOT, ...)
}

check_inputs <- function(paths, label = "input files") {
    missing_paths <- unname(paths[!file.exists(paths)])
    if (length(missing_paths)) stop("Missing ", label, ":\n", paste(missing_paths, collapse = "\n"), call. = FALSE)
    invisible(paths)
}

use_packages <- function(attach, namespace = character()) {
    packages <- unique(c(attach, namespace))
    missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_packages)) {
        stop("Install required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
    }
    suppressPackageStartupMessages(invisible(lapply(attach, library, character.only = TRUE)))
}
