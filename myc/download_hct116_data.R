if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)

proj_root <- here::here()
setwd(proj_root)

raw_file <- file.path(proj_root, "myc/data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")

if (!file.exists(raw_file)) {
  # download data -- hct116 (makeDirectory=TRUE creates the GSM4286760/ subdirectory under baseDir)
  getGEOSuppFiles("GSM4286760", baseDir = file.path(proj_root, "myc/data"))
} else {
  cat("HCT116 raw data already present at", raw_file, "- skipping download.\n")
}
