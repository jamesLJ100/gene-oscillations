# Install dyngen (if needed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dyngen")

library(dyngen)


# Load the realnets data
data("realnets")

# List available GRN names
names(realnets)
