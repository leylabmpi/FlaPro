# r_433_nb_install_deps.R
# Dependencies to install in the R 4.3.3+ environment for Jupyter Notebook

options(warn = 1)

# 1. Installantion managers (BiocManager, remotes, devtools)
installer_pkgs <- c("BiocManager", "remotes", "devtools")
missing_installers <- installer_pkgs[!(installer_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_installers) > 0) {
  install.packages(missing_installers)
}

# 2. from Bioconductor
bioc_packages <- c(
  "phyloseq"
)

# 3. from CRAN
cran_packages <- c(
  "readxl", "zCompositions", "reshape2", "ggpubr", "data.table", 
  "ggplot2", "tibble", "vegan", "foreach", "doParallel", 
  "lme4", "textshape", "mgcv", "MASS", "ggsci", 
  "cluster", "tidyr", "readr", "broom", "lmerTest", 
  "broom.mixed", "furrr", "testthat", "ggrepel", "pheatmap", 
  "dplyr", "purrr", "tidytable", "openxlsx", "caret", 
  "vctrs", "rcurl", "partitions", "igraph", "fontawesome", 
  "systemfonts", "bestNormalize", "coda", "gridExtra", 
  "jsonlite", "VennDiagram", "IRkernel"
)

# Installation of the absent CRAN packages
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  message("Installing packages from CRAN: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

# Installation of the absent Bioconductor packages
missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(missing_bioc) > 0) {
  message("Installing packages from Bioconductor: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc)
}

# 4. Specific packages from GitHub and Bitbucket (devtools or remotes)
message("Installing packages from GitHub / Bitbucket...")
devtools::install_github("tpq/balance", quiet = TRUE)
devtools::install_github("malucalle/selbal", quiet = TRUE)
devtools::install_bitbucket("knomics/nearestbalance", quiet = TRUE)
devtools::install_github("leylabmpi/LeyLabRMisc", quiet = TRUE)

# 5. R-kernell for Jupyter
if (requireNamespace("IRkernel", quietly = TRUE)) {
  tryCatch({
    IRkernel::installspec(user = TRUE)
  }, error = function(e) {
    warning("Cannot install IRkernel")
  })
}

message("All dependencies are set up!")