rm(list=ls())
library(data.table)
library(Matrix)
library(rhdf5)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "common/cell_purity_filter.R"))

# ============================================================================
# hIPSC preprocessing, mirroring preprocess_mesc.R / preprocess_mme95.R: load raw GEO
# counts, attach Louvain cluster labels, QC-filter cells, apply the same neighbor-purity
# cell filter as preprocess_mme95.R / preprocess_mesc.R (see common/cell_purity_filter.R),
# save raw UMI counts per cluster.
#
# Louvain clustering source: the "full" Diaz2019 hIPSC SPRING session —
# https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/hIPSC/full
# As with mESC, no separate fetch from SPRING is needed: the raw GEO CellData's
# own `louvain` column already matches that SPRING session's labels exactly —
# verified directly against categorical_coloring_data.json (100% match, all
# 14,750 cells, no offset/reordering).
# ============================================================================

geo_files <- file.path(proj_root, "segmentation_clock/data/GSE114186", c(
  "GSE114186_hsIPSC_CellData.csv.gz", "GSE114186_hsIPSC_GeneData.csv.gz", "GSE114186_hsIPSC_X.csv.gz"
))
if (!all(file.exists(geo_files))) source(file.path(proj_root, "segmentation_clock/download_data.R"))

# --- Load GEO data ---
celldata <- fread("segmentation_clock/data/GSE114186/GSE114186_hsIPSC_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("segmentation_clock/data/GSE114186/GSE114186_hsIPSC_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("segmentation_clock/data/GSE114186/GSE114186_hsIPSC_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

# --- No stage/cluster subsetting: the full hIPSC dataset is used, matching the
# SPRING "full" session directly ---
X_pm    <- X_raw
meta_pm <- celldata
meta_pm$cluster <- celldata$louvain

cat("Dimensions of hIPSC count matrix:", dim(X_pm), "\n")
cat("\nCluster counts:\n")
print(table(meta_pm$cluster))

# --- Drop the "stressed" cluster and the 5 near-empty "unknown" clusters ---
# "str" is a technical stress-response signature (dissociation/handling stress), not a
# real differentiation stage, so it isn't a meaningful comparison point alongside the
# d0/d1/d2/d3-4 time-course. unk1-5 are ~20-30 cells each pre-QC - too small to be a
# reliable cluster either way.
excluded_clusters <- c("str", "unk1", "unk2", "unk3", "unk4", "unk5")
cluster_flag <- !(meta_pm$cluster %in% excluded_clusters)
cat(sprintf("\nExcluding clusters (%s): keeping %d / %d cells\n",
            paste(excluded_clusters, collapse = ", "), sum(cluster_flag), length(cluster_flag)))

X_pm    <- X_pm[cluster_flag, ]
meta_pm <- meta_pm[cluster_flag, ]

# --- Cell QC: drop low-count cells and cells with high mitochondrial fraction ---
# Same thresholds/rationale as preprocess_mme95.R / preprocess_mesc.R. This gene panel uses
# classic bare mitochondrial gene symbols (ND1, COX1, ATP6, TRNL1, ...) with no
# distinguishing prefix (unlike mouse's "mt-" prefix), and a bare "^MT" regex would
# incorrectly catch unrelated genes (e.g. the MT1A/MT2A metallothionein family) - so
# the 37 human mitochondrial genes (13 protein-coding, 2 rRNA, 22 tRNA) are listed
# explicitly instead.
MIN_COUNTS    <- 1000
MAX_MITO_FRAC <- 0.10

mito_gene_symbols <- c(
  "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
  "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
  "RNR1", "RNR2",
  "TRNA", "TRNC", "TRND", "TRNE", "TRNF", "TRNG", "TRNH", "TRNI", "TRNK",
  "TRNL1", "TRNL2", "TRNM", "TRNN", "TRNP", "TRNQ", "TRNR", "TRNS1", "TRNS2",
  "TRNT", "TRNV", "TRNW", "TRNY"
)
mito_genes <- intersect(mito_gene_symbols, colnames(X_pm))

total_counts <- rowSums(X_pm)
mito_frac    <- rowSums(X_pm[, mito_genes, drop = FALSE]) / total_counts

cat("\nMitochondrial genes found:", length(mito_genes), "/", length(mito_gene_symbols), "\n")
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
# Cluster names contain spaces/hyphens ("d3-4 aPSM") - sanitised for filenames below.
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
    counts_file <- file.path(out_dir, paste0("hIPSC_", cl_safe, "_counts.h5"))
    mat2hdf(cl_counts, counts_file)

    cat("Saved:", counts_file, "—", ncol(cl_counts), "cells,", nrow(cl_counts), "genes\n")
  }
}

save_cluster_counts(X_pm, meta_pm, "segmentation_clock/data/hIPSC/counts/raw")

# --- Neighbor-label-purity cell QC: drop cells that sit ambiguously between clusters ---
# See common/cell_purity_filter.R for the method (k-NN label agreement, computed directly on
# scaled HVG expression - not a PCA/tSNE/UMAP/force-atlas projection - so the result doesn't
# depend on a choice of dimensionality-reduction method or number of components).
purity_result <- neighbor_purity_filter(X_pm, meta_pm, k_neighbors = 20L, purity_min = 0.7)
meta_pm     <- purity_result$meta
purity_flag <- meta_pm$kept

save_cluster_counts(X_pm[purity_flag, ], meta_pm[purity_flag, ], "segmentation_clock/data/hIPSC/counts/purity_filtered")
