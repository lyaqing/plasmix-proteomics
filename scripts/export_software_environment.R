# Export the software environment used by the Plasmix proteomics analysis.
# Outputs package, R, Bioconductor, session and optional Python environment records.

source("scripts/_project_setup.R")

# The list includes attached packages and packages intentionally called with package::function.
analysis_packages <- c(
    "bio3d", "brms", "circlize", "clue", "clusterProfiler", "ComplexHeatmap", "cowplot",
    "data.table", "dplyr", "ggpattern", "ggh4x", "ggplot2", "ggplotify", "ggpp", "ggpubr",
    "ggrepel", "ggridges", "ggtext", "grid", "irr", "lightgbm", "limma", "matrixStats",
    "metap", "openxlsx", "org.Hs.eg.db", "patchwork", "pbapply", "pdp", "Peptides",
    "purrr", "randomForest", "readr", "readxl", "reticulate", "RColorBrewer", "scales",
    "showtext", "smplot2", "SomaDataIO", "stringr", "tibble", "tidyr", "tidyverse", "xgboost"
)

package_table <- data.frame(
    package = analysis_packages,
    installed = vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE),
    version = NA_character_,
    library_path = NA_character_,
    stringsAsFactors = FALSE
)
installed_index <- which(package_table$installed)
package_table$version[installed_index] <- vapply(analysis_packages[installed_index], function(pkg) as.character(utils::packageVersion(pkg)), character(1))
package_table$library_path[installed_index] <- vapply(analysis_packages[installed_index], function(pkg) dirname(find.package(pkg)), character(1))

runtime_rows <- data.frame(
    package = c("R", "Bioconductor"),
    installed = TRUE,
    version = c(
        paste(R.version$major, R.version$minor, sep = "."),
        if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else NA_character_
    ),
    library_path = c(R.home(), NA_character_),
    stringsAsFactors = FALSE
)

environment_table <- rbind(runtime_rows, package_table)
utils::write.table(environment_table, "results/software_environment.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
capture.output(sessionInfo(), file = "results/session_info.txt")

if (requireNamespace("reticulate", quietly = TRUE)) {
    capture.output(reticulate::py_config(), file = "results/python_environment.txt")
    freesasa_available <- reticulate::py_module_available("freesasa")
    freesasa_version <- if (freesasa_available) {
        tryCatch(reticulate::py_eval("__import__('importlib.metadata', fromlist=['version']).version('freesasa')"), error = function(e) NA_character_)
    } else NA_character_
    python_packages <- data.frame(module = "freesasa", installed = freesasa_available, version = freesasa_version)
    utils::write.table(python_packages, "results/python_packages.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

message("Software environment records were written to results/.")
