rm(list=ls())

library(data.table)
library(Matrix)
library(rhdf5)
library(ggplot2)
library(patchwork)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "common/cell_purity_filter.R"))

geo_files <- file.path(proj_root, "segmentation_clock/data/GSE114186", c(
  "GSE114186_mmE95_CellData.csv.gz", "GSE114186_mmE95_GeneData.csv.gz", "GSE114186_mmE95_X.csv.gz"
))
if (!all(file.exists(geo_files))) source(file.path(proj_root, "segmentation_clock/download_data.R"))

# --- Load SPRING data ---
# cell_groupings.csv / original_cell_indices.txt exported from the "mme95" subplot of the
# Diaz2019 E95PM SPRING session:
# https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/E95PM/full
# original_cell_indices.txt gives 0-based row positions into that SPRING dataset's own cell
# ordering, which corresponds to celldata[celldata$leiden %in% stage1_clusters, ] below (NOT
# celldata$louvain — verified against the SPRING data by matching raw counts for marker genes).
# These are exported by hand from the interactive SPRING viewer (a lasso selection saved to
# file) rather than served as a downloadable resource, so unlike the GEO data above they can't
# be fetched automatically -- fail early with a clear message if they're missing.
spring_files <- file.path(proj_root, "segmentation_clock/data/mme95", c(
  "cell_groupings.csv", "original_cell_indices.txt"
))
if (!all(file.exists(spring_files))) {
  stop(
    "Missing SPRING-exported file(s): ",
    paste(spring_files[!file.exists(spring_files)], collapse = ", "),
    ". These come from a manual lasso-selection export in the interactive SPRING viewer ",
    "(https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/E95PM/full) ",
    "and aren't obtainable via GEOquery or any other automated download."
  )
}

groupings <- fread("segmentation_clock/data/mme95/cell_groupings.csv", header = FALSE)
indices   <- as.integer(readLines("segmentation_clock/data/mme95/original_cell_indices.txt"))

label_names <- groupings[[1]]
labels_mat  <- as.data.frame(t(groupings[, -1, with = FALSE]))
colnames(labels_mat) <- label_names
rownames(labels_mat) <- NULL

# --- Load GEO data ---
celldata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

# --- Subset to Stage 1 clusters ---
# Needed only to resolve original_cell_indices.txt's 0-based row positions to actual cell
# identities (rownames) - those positions are into this Stage-1-subsetted ordering, not the
# raw celldata (see the comment on spring_files above).
stage1_clusters <- c('aNTB1','aNTB2','pNTB','NMP','pPSM','aPSM','SOM')
stage1_flag  <- celldata$leiden %in% stage1_clusters
X_sub        <- X_raw[stage1_flag, ]

stopifnot("SPRING indices exceed stage1 subset" = max(indices) < nrow(X_sub))

pm_subplot_ids <- rownames(X_sub)[indices + 1]

# --- Build the "mme95 PM subplot" population directly from the authors' own SPRING export ---
# cell_groupings.csv's `louvain` row already IS the authors' final cluster assignment for
# these cells (aPSM/dmSOM/MPC/NMP/pPSM/scSOM), exported from their own Louvain run on the
# broader Neural+PSM subset - https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?datasets/Diaz2019/E95PM/full.
# There's no need to re-derive a clustering ourselves (an earlier version of this script did,
# and it neither matched their cluster count nor could be mapped back to their cluster names -
# Louvain cluster IDs/counts aren't stable across reruns or population scope). Using their
# labels directly is both simpler and exactly reproduces their grouping.
X_pm    <- X_sub[pm_subplot_ids, ]
meta_pm <- celldata[pm_subplot_ids, ]

meta_pm$library_id <- labels_mat$library_id
meta_pm$cluster    <- labels_mat$louvain

cat("Dimensions of PM count matrix:", dim(X_pm), "\n")
cat("Cluster counts (from the authors' own SPRING export):\n")
print(table(meta_pm$cluster))

# --- Cell QC: drop low-count cells and cells with high mitochondrial fraction ---
# Per the Cyclum paper's recommendation for droplet data.
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

# --- Save raw UMI counts (cells x genes -> genes x cells) ---
# This is inDrops (droplet/UMI) data (GSM3137206), and per the Cyclum paper's methods, droplet
# data should use "read counts normalized and log2 transformed" rather than TPM — gene-length
# normalization corrects a bias that read-count protocols have and UMI counting doesn't. The
# paper also says not to filter out genes. Both downstream algorithms take these raw UMI counts
# directly and do their own normalization: scPrisma's runner (algorithms/run_scPrisma.py) applies
# sc.pp.normalize_total/log1p, and Cyclum's runner (algorithms/run_cyclum.py) applies CPM + log2 +
# per-gene scaling — so there's a single shared input file type instead of a separate TPM/CPM
# preprocessing step here.

# Save each cluster in `meta` as a separate genes x cells .h5 file under out_dir.
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

    counts_file <- file.path(out_dir, paste0("mmE95_PM_", cl, "_counts.h5"))
    mat2hdf(cl_counts, counts_file)

    cat("Saved:", counts_file, "—", ncol(cl_counts), "cells,", nrow(cl_counts), "genes\n")
  }
}

save_cluster_counts(X_pm, meta_pm, "segmentation_clock/data/mme95/counts/raw")

# --- Neighbor-label-purity cell QC: drop cells that sit ambiguously between clusters ---
# See common/cell_purity_filter.R for the method (k-NN label agreement, computed directly on
# scaled HVG expression - not a PCA/tSNE/UMAP/force-atlas projection - so the result doesn't
# depend on a choice of dimensionality-reduction method or number of components).
K_NEIGHBORS <- 20L
PURITY_MIN  <- 0.7   # keep cells whose neighbours agree with their own label >= 70% of the time

purity_result <- neighbor_purity_filter(X_pm, meta_pm, k_neighbors = K_NEIGHBORS, purity_min = PURITY_MIN)
meta_pm     <- purity_result$meta
purity_flag <- meta_pm$kept

save_cluster_counts(X_pm[purity_flag, ], meta_pm[purity_flag, ], "segmentation_clock/data/mme95/counts/purity_filtered")

# Save per-cell cluster/purity/kept results (keyed by rowname = cell identity) so
# plot_mme95_knn_graph.R can overlay them on the paper's own graph layout without
# recomputing QC/purity itself.
saveRDS(meta_pm, "segmentation_clock/data/mme95/counts/purity_filter_results.rds")

# --- PCA plot: clusters before/after neighbor-purity filtering ---
# For visualization only - a 2D projection of the same scaled HVG expression the purity filter
# was computed from. The filtering decision above doesn't depend on this projection.
pca_embedding <- prcomp(purity_result$scaled_expr)$x[, 1:2]

figures_dir <- file.path(proj_root, "segmentation_clock/figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

plot_df <- data.frame(
  PC1     = pca_embedding[, 1],
  PC2     = pca_embedding[, 2],
  cluster = factor(meta_pm$cluster),
  kept    = purity_flag
)

p_before <- ggplot(plot_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(title = "Before neighbor-purity filtering", color = "Cluster") +
  theme_bw(base_size = 12)

plot_df$status <- factor(
  ifelse(plot_df$kept, as.character(plot_df$cluster), "Removed"),
  levels = c(levels(plot_df$cluster), "Removed")
)
status_colors <- setNames(
  c(scales::hue_pal()(nlevels(plot_df$cluster)), "grey80"),
  levels(plot_df$status)
)

p_after <- ggplot(plot_df, aes(PC1, PC2, color = status)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  labs(title = "After neighbor-purity filtering (grey = removed)", color = "") +
  theme_bw(base_size = 12)

pca_plot <- p_before + p_after
ggsave(file.path(figures_dir, "pca_purity_filtering.pdf"), pca_plot, width = 12, height = 5, dpi = 300)
cat("\nSaved PCA before/after filtering plot to:", file.path(figures_dir, "pca_purity_filtering.pdf"), "\n")
