if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dorothea")  # installs the regulon data package

library(dorothea)

# human regulons
data(dorothea_hs, package = "dorothea")


