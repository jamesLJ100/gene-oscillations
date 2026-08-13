rm(list=ls())

# Reproduces the source paper's own visualization - "k-NN graph of mouse neural tube, PSM and
# somite clusters at E9.5 (2,340 cells, 20 principal component dimensions)" - by replicating
# their exact Stage-2 embedding recipe (filter_cells -> normalize_total -> log1p -> HVG500 ->
# scale -> PCA20 -> neighbors(20nn,20pcs) -> ForceAtlas2 graph layout), then overlays our
# neighbor-purity filter results from preprocess_mme95.R on that same layout, for a direct
# visual comparison against the paper's figure.
# https://github.com/wagnerde/Diaz2019/blob/master/mmE95/Diaz2019_mmE95.ipynb (cells 31-40)
#
# Exploratory/diagnostic only - does not feed preprocess_mme95.R's pipeline, which uses the
# authors' own SPRING-exported PM-subplot cluster labels directly rather than this recomputed
# embedding. Run preprocess_mme95.R first so purity_filter_results.rds exists.

# R's own linked libraries and Python's scientific stack (loaded below via reticulate) each
# bundle their own OpenMP runtime on macOS, which crashes on double-init unless this is set
# first - a known R+reticulate interop issue, not specific to this script's logic.
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")

library(data.table)
library(reticulate)
library(ggplot2)

proj_root <- here::here()
setwd(proj_root)

purity_results_file <- "segmentation_clock/data/mme95/counts/purity_filter_results.rds"
if (!file.exists(purity_results_file)) {
  stop("Missing ", purity_results_file, " - run preprocess_mme95.R first.")
}
purity_results <- readRDS(purity_results_file)

# --- Load GEO + SPRING data (same source files as preprocess_mme95.R) ---
geo_files <- file.path(proj_root, "segmentation_clock/data/GSE114186", c(
  "GSE114186_mmE95_CellData.csv.gz", "GSE114186_mmE95_GeneData.csv.gz", "GSE114186_mmE95_X.csv.gz"
))
if (!all(file.exists(geo_files))) source(file.path(proj_root, "segmentation_clock/download_data.R"))

celldata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_CellData.csv.gz") |> as.data.frame()
rownames(celldata) <- make.unique(as.character(celldata[[1]]))
celldata <- celldata[, -1, drop = FALSE]

genedata <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_GeneData.csv.gz") |> as.data.frame()
rownames(genedata) <- make.unique(as.character(genedata[[1]]))
genedata <- genedata[, -1, drop = FALSE]

X_raw <- fread("segmentation_clock/data/GSE114186/GSE114186_mmE95_X.csv.gz", header = FALSE) |> as.matrix()
rownames(X_raw) <- rownames(celldata)
colnames(X_raw) <- rownames(genedata)

stage1_clusters <- c('aNTB1','aNTB2','pNTB','NMP','pPSM','aPSM','SOM')
stage1_flag     <- celldata$leiden %in% stage1_clusters
X_sub           <- X_raw[stage1_flag, ]

cat("Full Neural+PSM (Stage-1) population:", nrow(X_sub), "cells\n")

# --- Replicate the notebook's exact Stage-2 preprocessing + graph embedding (cells 33-36) ---
use_condaenv("scPrisma_env", required = TRUE)
sc <- import("scanpy")
ad <- import("anndata")

adata <- ad$AnnData(X = X_sub)
adata$obs_names <- rownames(X_sub)
adata$var_names <- colnames(X_sub)

sc$pp$filter_cells(adata, min_genes = 250L)
sc$pp$normalize_total(adata)
sc$pp$log1p(adata)
sc$pp$highly_variable_genes(adata, n_top_genes = 500L)
sc$pp$scale(adata)
sc$tl$pca(adata, n_comps = 20L, svd_solver = "arpack")
sc$pp$neighbors(adata, n_neighbors = 20L, n_pcs = 20L, use_rep = "X_pca")
sc$tl$draw_graph(adata, layout = "fa", iterations = 400L, random_state = 0L)

cat("Cells retained after filter_cells(min_genes=250):", nrow(adata$obs), "\n")

fa_coords <- as.matrix(adata$obsm[["X_draw_graph_fa"]])
# Matches the notebook's own cosmetic rotation (cell 36), for visual consistency with the
# paper's figure / SPRING layout orientation.
fa_coords <- fa_coords[, 2:1]
fa_coords[, 1] <- -fa_coords[, 1]

plot_df <- data.frame(
  FA1     = fa_coords[, 1],
  FA2     = fa_coords[, 2],
  cell_id = as.character(reticulate::py_to_r(adata$obs_names$values)),
  stringsAsFactors = FALSE
)

# --- Overlay our results: the authors' 6 PM cluster names for PM-lineage cells (everything
# else in the Neural+PSM population is labelled "Other (neural)"), and our neighbor-purity
# filter's keep/remove call for PM cells. ---
pm_match <- match(plot_df$cell_id, rownames(purity_results))
plot_df$cluster <- ifelse(is.na(pm_match), "Other (neural)", purity_results$cluster[pm_match])
plot_df$kept    <- ifelse(is.na(pm_match), NA, purity_results$kept[pm_match])

cat("\nCluster counts on the graph:\n")
print(table(plot_df$cluster))

plot_df$status <- factor(
  ifelse(is.na(plot_df$kept), "Other (neural)",
         ifelse(plot_df$kept, plot_df$cluster, "Removed (low purity)")),
  levels = c(setdiff(unique(plot_df$cluster), "Other (neural)"), "Removed (low purity)", "Other (neural)")
)

p_clusters <- ggplot(plot_df, aes(FA1, FA2, color = cluster)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c(setNames(
    scales::hue_pal()(length(unique(plot_df$cluster)) - 1),
    setdiff(unique(plot_df$cluster), "Other (neural)")
  ), "Other (neural)" = "grey80")) +
  labs(title = "Neural tube, PSM and somite clusters at E9.5 (k-NN graph, 20 PCs)",
       color = "Cluster") +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

status_colors <- c(setNames(
  scales::hue_pal()(length(levels(plot_df$status)) - 2),
  setdiff(levels(plot_df$status), c("Removed (low purity)", "Other (neural)"))
), "Removed (low purity)" = "black", "Other (neural)" = "grey80")

p_purity <- ggplot(plot_df, aes(FA1, FA2, color = status)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  labs(title = "Neighbor-purity filter result (black = removed)", color = "") +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

figures_dir <- file.path(proj_root, "segmentation_clock/figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

ggsave(file.path(figures_dir, "mme95_knn_graph_clusters.pdf"), p_clusters, width = 8, height = 7, dpi = 300)
ggsave(file.path(figures_dir, "mme95_knn_graph_purity_filter.pdf"), p_purity, width = 8, height = 7, dpi = 300)
cat("\nSaved:\n",
    file.path(figures_dir, "mme95_knn_graph_clusters.pdf"), "\n",
    file.path(figures_dir, "mme95_knn_graph_purity_filter.pdf"), "\n")
