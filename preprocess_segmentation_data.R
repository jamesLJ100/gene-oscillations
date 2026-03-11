library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
library(biomaRt)
library(rhdf5)

source(file.path(proj_root, "hdfrw.R"))


# --- Load SPRING data ---
groupings <- fread("data/mme95/cell_groupings.csv")
indices   <- as.integer(readLines("data/mme95/original_cell_indices.txt"))

label_names <- groupings[[1]]
labels_mat  <- as.data.frame(t(groupings[, -1, with = FALSE]))
colnames(labels_mat) <- label_names
rownames(labels_mat) <- NULL

# --- Load GEO data ---
celldata <- fread("data/GSE114186/GSE114186_mmE95_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("data/GSE114186/GSE114186_mmE95_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("data/GSE114186/GSE114186_mmE95_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

# --- Subset to Stage 1 clusters ---
stage1_clusters <- c('aNTB1','aNTB2','pNTB','NMP','pPSM','aPSM','SOM')
stage1_flag  <- celldata$louvain %in% stage1_clusters
X_sub        <- X_raw[stage1_flag, ]
celldata_sub <- celldata[stage1_flag, ]

stopifnot("SPRING indices exceed stage1 subset" = max(indices) < nrow(X_sub))

# --- Select PM cells using 0-based SPRING indices ---
X_pm    <- X_sub[indices + 1, ]
meta_pm <- celldata_sub[indices + 1, ]

meta_pm$cluster    <- labels_mat$louvain
meta_pm$library_id <- labels_mat$library_id

cat("Dimensions of PM count matrix:", dim(X_pm), "\n")
cat("\nCluster counts:\n")
print(table(meta_pm$cluster))

# --- Get mouse gene lengths from biomaRt ---
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_info <- getBM(
  attributes = c("mgi_symbol", "start_position", "end_position"),
  filters    = "mgi_symbol",
  values     = colnames(X_pm),
  mart       = mart
)

cat("Genes queried:", ncol(X_pm), "\n")
cat("Genes returned by biomaRt:", nrow(gene_info), "\n")
cat("Duplicated symbols:", sum(duplicated(gene_info$mgi_symbol)), "\n")

gene_info <- gene_info[!duplicated(gene_info$mgi_symbol), ]
gene_info$length_kb <- abs(gene_info$end_position - gene_info$start_position) / 1000

# Subset to genes with length info
genes_keep <- colnames(X_pm)[colnames(X_pm) %in% gene_info$mgi_symbol]
X_sub      <- X_pm[, genes_keep]
gene_info  <- gene_info[match(genes_keep, gene_info$mgi_symbol), ]

cat("Genes after filtering:", ncol(X_sub), "\n")

# --- UMI -> TPM (cells x genes -> genes x cells for calculation) ---
X_t  <- t(as.matrix(X_sub))             # genes x cells
rpk  <- X_t / gene_info$length_kb       # normalise by gene length
tpm  <- t(t(rpk) / colSums(rpk)) * 1e6  # normalise each cell to sum to 1e6

rownames(tpm) <- genes_keep
colnames(tpm) <- meta_pm$unique_cell_id

cat("\nTPM dimensions (genes x cells):", dim(tpm), "\n")
cat("NAs in tpm:", sum(is.na(tpm)), "\n")
cat("Column sums (should all be 1e6):\n")
print(summary(colSums(tpm)))

# --- Save each cluster as a separate .h5 file ---
dir.create("data/mme95/tpm", showWarnings = FALSE)

clusters <- unique(meta_pm$cluster)
cat("\nSaving clusters:", paste(clusters, collapse = ", "), "\n")

for (cl in clusters) {
  cl_cells <- meta_pm$unique_cell_id[meta_pm$cluster == cl]
  cl_tpm   <- tpm[, cl_cells]
  
  outfile <- file.path("data/mme95/tpm", paste0("mmE95_PM_", cl, "_tpm.h5"))
  mat2hdf(cl_tpm, outfile)
  
  cat("Saved:", outfile, "—", ncol(cl_tpm), "cells,", nrow(cl_tpm), "genes\n")
}
