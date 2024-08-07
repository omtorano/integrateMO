# Check if BiocManager is installed, if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor dependency
BiocManager::install("sva") #sva does not install even when listed under remotes with bioc channel in DESCRIPTION
BiocManager::install("impute") #impute is a dependency of WGCNA but does not automatically install with the installation of WGCNA
BiocManager::install("preprocessCore") #preprocessCore is a depencency of WGCNA but does not automatically isntall with the installation of WGCNA
BiocManager::install("WGCNA") #WGCNA does not install even when listed under remotes with bioc channell in DESCRIPTION

# Install CRAN dependencies
# GLS error in IMIFA installation indicates GLS needs to be installed on root Ubuntu system (sudo apt install libgsl-dev)
install.packages("IMIFA") #IMIFA does not install properly when listed in imports of DESCRIPTION
