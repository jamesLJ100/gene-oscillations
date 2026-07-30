if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery")
}
library(GEOquery)

proj_root <- here::here()
setwd(proj_root)

tar_file <- file.path(proj_root, "nfkb/data/GSE162992/GSE162992_RAW.tar")
raw_dir  <- file.path(proj_root, "nfkb/data/GSE162992/GSE162992_RAW")

if (!file.exists(tar_file)) {
  # download data -- nfkb
  getGEOSuppFiles("GSE162992", baseDir = file.path(proj_root, "nfkb/data"))
} else {
  cat("NF-kB raw archive already present at", tar_file, "- skipping download.\n")
}

if (!dir.exists(raw_dir) || length(list.files(raw_dir)) == 0) {
  untar(tar_file, exdir = raw_dir)
} else {
  cat("NF-kB raw archive already extracted at", raw_dir, "- skipping untar.\n")
}
