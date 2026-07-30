if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)

proj_root <- here::here()
setwd(proj_root)

geo_dir <- file.path(proj_root, "segmentation_clock/data/GSE114186")
required_files <- file.path(geo_dir, c(
  "GSE114186_mmE95_CellData.csv.gz",  "GSE114186_mmE95_GeneData.csv.gz",  "GSE114186_mmE95_X.csv.gz",
  "GSE114186_mmESC_CellData.csv.gz",  "GSE114186_mmESC_GeneData.csv.gz",  "GSE114186_mmESC_X.csv.gz",
  "GSE114186_hsIPSC_CellData.csv.gz", "GSE114186_hsIPSC_GeneData.csv.gz", "GSE114186_hsIPSC_X.csv.gz"
))

if (!all(file.exists(required_files))) {
  # download data -- GSE114186 (mouse E9.5 PM, mESC, and hIPSC segmentation clock data)
  getGEOSuppFiles("GSE114186", baseDir = file.path(proj_root, "segmentation_clock/data"))
} else {
  cat("GSE114186 data already present at", geo_dir, "- skipping download.\n")
}
