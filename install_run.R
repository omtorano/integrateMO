# Check if BiocManager is installed, if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor dependency
BiocManager::install("sva")
BiocManager::install("WGCNA")

# Install CRAN dependencies
# GLS error in IMIFA installation indicates GLS needs to be installed on root Ubuntu system (sudo apt install libgsl-dev)
install.packages("IMIFA")
