#!/usr/bin/env Rscript

# Dependency Pre-installation Script
# Automates installation from CRAN/Bioconductor/GitHub/Local sources with error handling

# Configure mirror for faster installation
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# Package Source
packages <- list(
  # Standard CRAN packages
  CRAN = c(
    "BiocManager", "devtools", "data.table", "pbmcapply", "tidyverse",
    "vroom", "reshape2", "seqinr", "stringr", "DescTools", 
    "argparse", "Matrix", "Seurat", "XML", "R.utils"
  ),
  
  # GitHub-hosted packages
  GitHub = list(
    bedtoolsr = "PhanstielLab/bedtoolsr"
  ),
  
  # Bioconductor packages
  Bioconductor = c(
    "GenomicRanges", "Biostrings", "Rsamtools", 
    "GenomicAlignments", "BiocIO", "restfulr"
  ),
  
  # Local source packages
  Local = list(
    rtracklayer = "src/rtracklayer_1.68.0.tar.gz"
  )
)

#  Check Function
safe_require <- function(pkg) {
  suppressWarnings(suppressPackageStartupMessages(
    require(pkg, character.only = TRUE, quietly = TRUE)
  ))
}

# Installation
install_with_retry <- function(install_func, pkg, max_attempts = 2, ...) {
  for (attempt in 1:max_attempts) {
    message(sprintf("Install attempt %d/%d: %s", attempt, max_attempts, pkg))
    result <- tryCatch({
      install_func(pkg, quiet = TRUE)
      if (safe_require(pkg)) return(TRUE)
      FALSE
    }, error = function(e) FALSE)
    if (result) return(TRUE)
    Sys.sleep(2^attempt) # Exponential backoff (2, 4, 8 seconds)
  }
  FALSE
}

# Main
tryCatch({
  # CRAN Packages
  for (pkg in packages$CRAN) {
    if (!safe_require(pkg)) {
      message(sprintf("[CRAN] Installing %s...", pkg))
      install.packages(pkg, quiet = TRUE)
      if (!safe_require(pkg)) stop("CRAN installation failed: ", pkg)
    }
  }

  # Bioconductor Packages
  for (pkg in packages$Bioconductor) {
    if (!safe_require(pkg)) {
      message(sprintf("[Bioconductor] Installing %s...", pkg))
      if (!install_with_retry(BiocManager::install, pkg = pkg, ask = FALSE)) {
        stop("Bioconductor installation failed: ", pkg)
      }
    }
  }

  # GitHub Packages
  for (pkg in names(packages$GitHub)) {
    if (!safe_require(pkg)) {
      repo <- packages$GitHub[[pkg]]
      message(sprintf("[GitHub] Installing %s from %s...", pkg, repo))
      devtools::install_github(repo)
      if (!safe_require(pkg)) stop("GitHub installation failed: ", repo)
    }
  }

  # Local Source
  for (pkg in names(packages$Local)) {
    if (!safe_require(pkg)) {
      tar_path <- packages$Local[[pkg]]
      if (!file.exists(tar_path)) stop("Local package not found: ", tar_path)
      message(sprintf("[Local] Installing %s from %s...", pkg, tar_path))
      install.packages(tar_path, repos = NULL, type = "source")
      if (!safe_require(pkg)) stop("Local installation failed: ", pkg)
    }
  }

  # Post-installation Verification
  message("\nInstallation Verification:")
  print(sapply(unlist(packages), function(pkg) safe_require(pkg)))

  message("\nAll dependencies successfully installed!")

}, error = function(e) {
  message("\nInstallation failed: ", conditionMessage(e))
  quit(status = 1)
})
