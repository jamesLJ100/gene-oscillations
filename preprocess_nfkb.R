#Download NFKB data
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
library(data.table)

#suppl_files <- getGEOSuppFiles("GSE162992", baseDir = "data/")

# Extract all files
#untar("data/GSE162992/GSE162992_RAW.tar", exdir = "data/GSE162992/GSE162992_RAW")

# See what's inside
#files <- list.files("data/GSE162992/GSE162992_RAW", full.names = TRUE)
#print(files)

# Install Seurat if needed
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
library(Seurat)

# Create directories
dir.create("data/GSE162992/GSE162992_RAW/WT", showWarnings = FALSE)
dir.create("data/GSE162992/GSE162992_RAW/IkBaMM", showWarnings = FALSE)

# Copy WT files
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_barcodes.tsv.gz",
          "data/GSE162992/GSE162992_RAW/WT/barcodes.tsv.gz")
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_features.tsv.gz",
          "data/GSE162992/GSE162992_RAW/WT/features.tsv.gz")
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_matrix.mtx.gz",
          "data/GSE162992/GSE162992_RAW/WT/matrix.mtx.gz")

# Copy IkBaMM files
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_barcodes.tsv.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/barcodes.tsv.gz")
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_features.tsv.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/features.tsv.gz")
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_matrix.mtx.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/matrix.mtx.gz")

# Then load as above
wt <- Read10X(data.dir = "data/GSE162992/GSE162992_RAW/WT")
ikbamm <- Read10X(data.dir = "data/GSE162992/GSE162992_RAW/IkBaMM")

