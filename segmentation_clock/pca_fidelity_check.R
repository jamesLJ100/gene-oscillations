rm(list=ls())
library(data.table)
library(jsonlite)
library(ggplot2)
library(cluster)

proj_root <- here::here()
setwd(proj_root)

# ============================================================================
# Diagnostic script: is our silhouette-filtering embedding a reasonable proxy
# for the original clustering, or is it meaningfully different?
#
# Per the source paper's methods:
#   - Mouse E9.5 PSM neighbour graph: 20 PCA dimensions (matches what we used).
#   - HVG selection: top 2,000 genes by a "bin-normalized overdispersion
#     metric" — NOT raw variance, which is what preprocess_mme95.R currently uses.
#   - Normalization: total-count normalize, log-normalize, scale to unit
#     variance / zero mean (i.e. z-score) before PCA.
#   - Clustering itself was on the FULL 2,340-cell pool, not our narrower
#     ~1,480-cell post-QC PM subset.
#
# Deliberately NOT changed here: cell-level QC. This script applies no
# filtering of its own (no min-genes, no mito) — the only thing being tested
# is embedding fidelity (PC count / HVG method / population scope). The "old"
# comparison arm below reproduces preprocess_mme95.R's current min_counts/mito QC
# exactly, unchanged, purely so the two methods are compared on the same cells.
#
# This script rebuilds the embedding with the corrected HVG method, computes
# silhouette scores in it, and compares against what preprocess_mme95.R's current
# (30 PC, naive-variance HVG, 1,480-cell) approach produces for the same
# cells. Exploratory only — does not feed the pipeline.
# ============================================================================

N_HVG   <- 2000
N_PCS   <- 20
N_BINS  <- 20
SIL_MIN <- 0

# --- Load raw GEO data (same as preprocess_mme95.R) ---
celldata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

# --- Full 2,340-cell Stage-1 pool (verified elsewhere to match SPRING's own
# row ordering exactly via marker-gene correlation) ---
stage1_clusters <- c('aNTB1','aNTB2','pNTB','NMP','pPSM','aPSM','SOM')
stage1_flag   <- celldata$leiden %in% stage1_clusters
X_full        <- X_raw[stage1_flag, ]
celldata_full <- celldata[stage1_flag, ]
cat("Full Stage-1 pool:", nrow(X_full), "cells,", ncol(X_full), "genes\n")
stopifnot("Expected 2340 cells matching the paper's full E95PM pool" = nrow(X_full) == 2340)

# --- Fetch Louvain labels for the full pool directly from SPRING ---
spring_url   <- "https://kleintools.hms.harvard.edu/tools/datasets/Diaz2019/E95PM/full/categorical_coloring_data.json"
full_cat     <- fromJSON(spring_url)
louvain_full <- full_cat[["louvain"]][["label_list"]]
stopifnot("Louvain label count must match the full cell pool" = length(louvain_full) == nrow(X_full))
celldata_full$cluster <- louvain_full

cat("\nCluster counts (full pool, no additional QC applied):\n")
print(table(celldata_full$cluster))

# --- Bin-normalized overdispersion HVG selection (Scanpy/Seurat "seurat" flavor) ---
# Raw variance is biased toward highly-expressed genes (Poisson-like mean-variance
# relationship). This selects genes with high variance RELATIVE to other genes at a
# similar expression level: bin by mean expression, z-score dispersion within bin.
select_hvg_dispersion <- function(logcpm, n_top, n_bins) {
  gene_mean <- colMeans(logcpm)
  gene_var  <- colMeans(logcpm^2) - gene_mean^2
  dispersion <- ifelse(gene_mean > 0, gene_var / gene_mean, 0)
  log_disp   <- log1p(dispersion)

  bins <- cut(gene_mean, breaks = n_bins, labels = FALSE)
  norm_disp <- ave(log_disp, bins, FUN = function(x) {
    s <- sd(x)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x)) / s
  })

  order(norm_disp, decreasing = TRUE)[seq_len(n_top)]
}

cpm    <- X_full / rowSums(X_full) * 1e6
logcpm <- log1p(cpm)
hvg    <- select_hvg_dispersion(logcpm, N_HVG, N_BINS)

# PCA with scale.=TRUE == "scaled to unit variance and zero mean" per the paper
pca <- prcomp(logcpm[, hvg], center = TRUE, scale. = TRUE, rank. = N_PCS)

plot_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cluster = factor(celldata_full$cluster))
p <- ggplot(plot_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(
    title    = "Full pool, 20 PCs, bin-normalized-dispersion HVGs (paper's stated method)",
    subtitle = "PCA scatter, not ForceAtlas2 - for cluster-separation comparison only",
    color    = "Cluster"
  ) +
  theme_bw(base_size = 12)

figures_dir <- file.path(proj_root, "segmentation_clock/figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
ggsave(file.path(figures_dir, "pca_full_population_20pc.pdf"), p, width = 7, height = 5.5, dpi = 300)

# --- Silhouette on the full population (all 10 categories provide context,
# matching the original clustering's actual scope) ---
cell_dist_full <- dist(pca$x)
sil_full <- silhouette(as.integer(factor(celldata_full$cluster)), cell_dist_full)
celldata_full$silhouette_new <- sil_full[, "sil_width"]

# --- Map down to the same "mme95" PM-subplot cells preprocess_mme95.R evaluates ---
groupings <- fread("segmentation_clock/data/mme95/cell_groupings.csv", header = FALSE)
indices   <- as.integer(readLines("segmentation_clock/data/mme95/original_cell_indices.txt"))
label_names <- groupings[[1]]
labels_mat  <- as.data.frame(t(groupings[, -1, with = FALSE]))
colnames(labels_mat) <- label_names

orig_ids <- rownames(celldata)[stage1_flag]
pm_ids   <- orig_ids[indices + 1]

meta_new <- celldata_full[pm_ids, ]
meta_new$cluster_pm <- labels_mat$louvain[match(pm_ids, orig_ids[indices + 1])]

# --- Reproduce preprocess_mme95.R's CURRENT approach for the same cells, for comparison ---
X_pm_old <- X_raw[stage1_flag, ][indices + 1, ]
rownames(X_pm_old) <- orig_ids[indices + 1]
meta_old <- data.frame(row.names = orig_ids[indices + 1], cluster = labels_mat$louvain)

mito_genes   <- grep("^mt-", colnames(X_pm_old), value = TRUE)
total_counts <- rowSums(X_pm_old)
mito_frac    <- rowSums(X_pm_old[, mito_genes, drop = FALSE]) / total_counts
qc_flag_old  <- total_counts >= 1000 & mito_frac <= 0.10
X_pm_old     <- X_pm_old[qc_flag_old, ]
meta_old     <- meta_old[qc_flag_old, , drop = FALSE]

cpm_old    <- X_pm_old / rowSums(X_pm_old) * 1e6
logcpm_old <- log1p(cpm_old)
gene_means_old <- colMeans(logcpm_old)
gene_vars_old  <- colMeans(logcpm_old^2) - gene_means_old^2
hvg_old <- order(gene_vars_old, decreasing = TRUE)[seq_len(N_HVG)]
pca_old <- prcomp(logcpm_old[, hvg_old], center = TRUE, scale. = TRUE, rank. = 30)
sil_old <- silhouette(as.integer(factor(meta_old$cluster)), dist(pca_old$x))
meta_old$silhouette_old <- sil_old[, "sil_width"]

# --- Compare old vs new silhouette scores for the cells present in both ---
common_ids <- intersect(rownames(meta_old), rownames(meta_new))
cat(sprintf("\nCells comparable between old and new methods: %d\n", length(common_ids)))

compare_df <- data.frame(
  cluster = meta_old[common_ids, "cluster"],
  old     = meta_old[common_ids, "silhouette_old"],
  new     = meta_new[common_ids, "silhouette_new"]
)

cat("\nCorrelation between old and new silhouette scores:",
    round(cor(compare_df$old, compare_df$new), 3), "\n")

old_keep <- compare_df$old >= SIL_MIN
new_keep <- compare_df$new >= SIL_MIN
cat(sprintf("\nOld method keeps: %d / %d\n", sum(old_keep), nrow(compare_df)))
cat(sprintf("New method keeps: %d / %d\n", sum(new_keep), nrow(compare_df)))
cat(sprintf("Cells where the keep/remove decision FLIPS: %d / %d (%.1f%%)\n",
            sum(old_keep != new_keep), nrow(compare_df), 100 * mean(old_keep != new_keep)))

cat("\nPer-cluster keep counts, old vs new:\n")
print(data.frame(
  cluster  = tapply(compare_df$cluster, compare_df$cluster, unique),
  n        = tapply(compare_df$cluster, compare_df$cluster, length),
  old_keep = tapply(old_keep, compare_df$cluster, sum),
  new_keep = tapply(new_keep, compare_df$cluster, sum)
))
