.onLoad <- function(libname, pkgname) {
  # Check if BiocManager is installed, if not, install it
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Install Bioconductor dependencies if not already installed
  BiocManager::install(c("sva"), update = FALSE, ask = FALSE)
  # Install WGCNA, IMIFA if not already installed
  if(!requireNamespace("WGCNA", quietly = TRUE)){
	install.packages("WGCNA")
  }
  if(!requireNamespace("IMIFA", quietly = TRUE)){
	install.packages("IMIFA")
  }
  
}