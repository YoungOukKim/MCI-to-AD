# ==============================================================================
# packages.R
#
# Purpose : Install all R packages required to reproduce analyses in
#           Kim Y†, Heo W†, Park S†, et al. (2025). bioRxiv.
#           https://github.com/YoungOukKim/MCI-to-AD
#
# Usage   : Rscript environment/packages.R
#           OR in R console: source("environment/packages.R")
#
# R version : 4.3.2 (tested)
# ==============================================================================

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  # Data wrangling
  "data.table",   # fast data manipulation (Scripts 01, 02, 02b, 04)
  "dplyr",        # tidyverse grammar (Scripts 02, 03)
  "tidyverse",    # ggplot2, readr, tidyr, purrr, stringr (Script 03)
  "lubridate",    # date handling (Script 03)
  # Statistical analysis
  "segmented",    # segmented regression + Davies' test (Script 04)
  "pROC",         # ROC / AUC with DeLong CI (Script 03)
  "survival",     # Cox proportional hazards (Script 03)
  # Visualization
  "ggplot2",      # all plots (Script 04)
  "survminer",    # Kaplan-Meier plots (Script 03)
  "patchwork",    # combined figure panels (Script 04, optional)
  # HDF5 for SEA-AD h5ad
  "zoo"           # rolling means (alternative to manual loop)
)

# ── Bioconductor packages ──────────────────────────────────────────────────────
bioc_pkgs <- c(
  "rhdf5"         # read SEA-AD .h5ad file (Script 01)
)

# ── Install CRAN ──────────────────────────────────────────────────────────────
message("Installing CRAN packages...")
install.packages(
  setdiff(cran_pkgs, rownames(installed.packages())),
  repos = "https://cloud.r-project.org",
  quiet = TRUE
)

# ── Install Bioconductor ──────────────────────────────────────────────────────
message("Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(
  setdiff(bioc_pkgs, rownames(installed.packages())),
  ask = FALSE, update = FALSE
)

# ── Verify ────────────────────────────────────────────────────────────────────
all_pkgs <- c(cran_pkgs, bioc_pkgs)
installed_ok <- sapply(all_pkgs, requireNamespace, quietly = TRUE)

if (all(installed_ok)) {
  message("\n✓  All ", length(all_pkgs), " packages installed successfully.")
} else {
  missing <- names(installed_ok)[!installed_ok]
  warning("The following packages could not be installed:\n  ",
          paste(missing, collapse = ", "))
}

message("\nR session info saved to environment/session_info.txt")
writeLines(capture.output(sessionInfo()),
           con = "environment/session_info.txt")
