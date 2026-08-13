rm(list=ls())
library(data.table)
library(Matrix)
library(rhdf5)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "common/cell_purity_filter.R"))

# ============================================================================
# mESC preprocessing, mirroring preprocess_mme95.R's approach for the mmE9.5 PM data:
# load raw GEO counts, attach Louvain cluster labels, QC-filter cells, apply the same
# neighbor-purity cell filter as preprocess_mme95.R (see common/cell_purity_filter.R),
# save raw UMI counts per cluster.
#
# Louvain clustering source: the "full" Diaz2019 mESC SPRING session —
# https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/mESC/full
# Unlike the E9.5 data, no separate fetch from SPRING is needed here: the raw
# GEO CellData's own `louvain` column already matches that SPRING session's
# labels exactly — verified directly by comparing all 21,478 cells against
# categorical_coloring_data.json fetched from the SPRING session (100% match,
# no offset/reordering). No `leiden` column exists for this dataset (unlike
# E9.5, which had a leiden/louvain mismatch requiring this exact check).
# ============================================================================

geo_files <- file.path(proj_root, "segmentation_clock/data/GSE114186", c(
  "GSE114186_mmESC_CellData.csv.gz", "GSE114186_mmESC_GeneData.csv.gz", "GSE114186_mmESC_X.csv.gz"
))
if (!all(file.exists(geo_files))) source(file.path(proj_root, "segmentation_clock/download_data.R"))

# --- Load GEO data ---
celldata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmESC_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmESC_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("segmentation_clock/data/GSE114186/GSE114186_mmESC_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

# --- No stage/cluster subsetting: the full mESC dataset is used, matching the
# SPRING "full" session directly (no curated subplot selection like mme95) ---
X_pm    <- X_raw
meta_pm <- celldata
meta_pm$cluster <- celldata$louvain

cat("Dimensions of mESC count matrix:", dim(X_pm), "\n")
cat("\nCluster counts:\n")
print(table(meta_pm$cluster))

# --- Cell QC: drop low-count cells and cells with high mitochondrial fraction ---
# Same thresholds/rationale as preprocess_mme95.R (Cyclum paper's recommendation for droplet data).
MIN_COUNTS    <- 1000
MAX_MITO_FRAC <- 0.10

mito_genes   <- grep("^mt-", colnames(X_pm), value = TRUE)
total_counts <- rowSums(X_pm)
mito_frac    <- rowSums(X_pm[, mito_genes, drop = FALSE]) / total_counts

cat("\nMitochondrial genes found:", length(mito_genes), "\n")
cat("Total counts per cell:\n"); print(summary(total_counts))
cat("Mitochondrial fraction per cell:\n"); print(summary(mito_frac))

qc_flag <- total_counts >= MIN_COUNTS & mito_frac <= MAX_MITO_FRAC
cat(sprintf(
  "\nQC filter (min_counts=%d, max_mito_frac=%.2f): keeping %d / %d cells\n",
  MIN_COUNTS, MAX_MITO_FRAC, sum(qc_flag), length(qc_flag)
))

X_pm    <- X_pm[qc_flag, ]
meta_pm <- meta_pm[qc_flag, ]

cat("\nCluster counts after QC:\n")
print(table(meta_pm$cluster))

# --- Save raw UMI counts (cells x genes -> genes x cells) ---
# Cluster names contain spaces ("d0 ESC", "d4-5 Krt") - sanitised for filenames below.
save_cluster_counts <- function(X, meta, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  X_t <- t(as.matrix(X))   # genes x cells
  rownames(X_t) <- colnames(X)
  colnames(X_t) <- meta$unique_cell_id

  clusters <- unique(meta$cluster)
  cat("\nSaving clusters to", out_dir, ":", paste(clusters, collapse = ", "), "\n")

  for (cl in clusters) {
    cl_cells  <- meta$unique_cell_id[meta$cluster == cl]
    cl_counts <- X_t[, cl_cells]

    cl_safe     <- gsub("[^A-Za-z0-9]+", "_", cl)
    counts_file <- file.path(out_dir, paste0("mmESC_", cl_safe, "_counts.h5"))
    mat2hdf(cl_counts, counts_file)

    cat("Saved:", counts_file, "—", ncol(cl_counts), "cells,", nrow(cl_counts), "genes\n")
  }
}

save_cluster_counts(X_pm, meta_pm, "segmentation_clock/data/mmESC/counts/raw")

# --- Neighbor-label-purity cell QC: drop cells that sit ambiguously between clusters ---
# See common/cell_purity_filter.R for the method (k-NN label agreement, computed directly on
# scaled HVG expression - not a PCA/tSNE/UMAP/force-atlas projection - so the result doesn't
# depend on a choice of dimensionality-reduction method or number of components).
purity_result <- neighbor_purity_filter(X_pm, meta_pm, k_neighbors = 20L, purity_min = 0.7)
meta_pm     <- purity_result$meta
purity_flag <- meta_pm$kept

save_cluster_counts(X_pm[purity_flag, ], meta_pm[purity_flag, ], "segmentation_clock/data/mmESC/counts/purity_filtered")
