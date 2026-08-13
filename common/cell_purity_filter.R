# Flags cells whose local expression neighbourhood disagrees with their own cluster label -
# used to drop ambiguous/boundary cells sitting between two clusters. For each cell, finds its
# k nearest neighbours by Euclidean distance in normalized, log-transformed, HVG-selected,
# scaled expression (HVGs chosen via scanpy's mean-corrected dispersion method, not plain
# variance ranking - plain variance is confounded with mean expression in count data and mostly
# just picks highly-expressed genes rather than genes that discriminate cell identity), then
# computes what fraction of those neighbours share the cell's own cluster label. This is
# deliberately NOT computed in a PCA/tSNE/UMAP/force-atlas-reduced space, so the result doesn't
# depend on a choice of dimensionality-reduction method or number of components.
#
# X:    raw UMI counts, cells x genes, one row per cell in `meta`'s order
# meta: data.frame, one row per cell (same order as X), with a `cluster` column
# Returns a list with:
#   meta        - `meta` with `purity` (fraction of k nearest neighbours sharing the cell's own
#                 cluster label) and `kept` (purity >= purity_min) columns added
#   scaled_expr - the HVG-selected, scaled expression matrix the filter was computed from
#                 (cells x n_hvg), for callers that want to reuse it (e.g. for a visualization)
neighbor_purity_filter <- function(X, meta, k_neighbors = 20L, purity_min = 0.5,
                                    n_hvg = 500L, conda_env = "scPrisma_env") {
  # R's own linked libraries (via Matrix/rhdf5) and Python's scientific stack (loaded below via
  # reticulate) each bundle their own OpenMP runtime on macOS, which crashes on double-init
  # unless this is set first - a known R+reticulate interop issue, not specific to this filter.
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")

  total_counts <- rowSums(X)
  cpm    <- sweep(X, 1, total_counts, "/") * 1e4
  logexp <- log1p(cpm)

  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <- reticulate::import("scanpy")
  ad <- reticulate::import("anndata")

  adata <- ad$AnnData(X = logexp)
  adata$var_names <- colnames(logexp)
  sc$pp$highly_variable_genes(adata, n_top_genes = as.integer(n_hvg))
  hvg <- colnames(logexp)[adata$var[["highly_variable"]]]

  # Scaling (zero mean, unit variance per gene) is the standard step between HVG selection and
  # distance/neighbour computation in scRNA-seq pipelines - without it, genes with higher
  # absolute expression would still dominate the Euclidean distance even after HVG selection.
  scaled_expr <- scale(logexp[, hvg])

  knn_idx <- FNN::get.knn(scaled_expr, k = k_neighbors)$nn.index
  neighbor_labels <- matrix(meta$cluster[knn_idx], nrow = nrow(knn_idx))

  meta$purity <- rowMeans(neighbor_labels == meta$cluster)
  meta$kept   <- meta$purity >= purity_min

  cat("\nNeighbor-purity summary (all clusters):\n")
  print(summary(meta$purity))
  cat("\nNeighbor-purity by cluster:\n")
  print(aggregate(purity ~ cluster, data = meta, FUN = function(x) round(summary(x), 3)))
  cat(sprintf(
    "\nNeighbor-purity filter (purity >= %.2f): keeping %d / %d cells\n",
    purity_min, sum(meta$kept), nrow(meta)
  ))
  cat("Cluster counts after neighbor-purity filter:\n")
  print(table(meta$cluster[meta$kept]))

  list(meta = meta, scaled_expr = scaled_expr)
}
