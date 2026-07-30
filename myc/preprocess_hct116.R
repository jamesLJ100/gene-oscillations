library(data.table)
library(rhdf5)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))

raw_file <- file.path(proj_root, "myc/data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")
if (!file.exists(raw_file)) source(file.path(proj_root, "myc/download_hct116_data.R"))

# read the file
umi_matrix <- fread(raw_file)
# set gene names as rownames
umi_matrix <- as.data.frame(umi_matrix)
rownames(umi_matrix) <- umi_matrix$V1
umi_matrix$V1 <- NULL

cat("=== UMI matrix ===\n")
cat("Dimensions:", dim(umi_matrix), "\n")
cat("NAs in umi_matrix:", sum(is.na(umi_matrix)), "\n")

# === Dataset summary statistics (full matrix, pre-gene-filter) ===
cat("\n=== Dataset Summary Statistics (pre-filter) ===\n")

num_features        <- nrow(umi_matrix)
num_cells           <- ncol(umi_matrix)
lib_sizes           <- colSums(umi_matrix)
median_library_size <- median(lib_sizes)
dropout_rate        <- sum(umi_matrix == 0) / (num_features * num_cells)

cpm        <- t(t(umi_matrix) / lib_sizes) * 1e6
log2cpm    <- log2(cpm + 1)

gene_mean_log2cpm  <- rowMeans(log2cpm)
median_avg_log2cpm <- median(gene_mean_log2cpm)

# use matrixStats::rowVars if available, otherwise fall back to apply
if (requireNamespace("matrixStats", quietly = TRUE)) {
  gene_var_log2cpm <- matrixStats::rowVars(as.matrix(log2cpm))
} else {
  gene_var_log2cpm <- apply(log2cpm, 1, var)
}
median_var_log2cpm <- median(gene_var_log2cpm)

cat(sprintf("num_cells:                   %d\n",   num_cells))
cat(sprintf("num_features:                %d\n",   num_features))
cat(sprintf("median_library_size:         %.1f\n", median_library_size))
cat(sprintf("dropout_rate:                %.4f\n", dropout_rate))
cat(sprintf("median_average_log2_cpm:     %.4f\n", median_avg_log2cpm))
cat(sprintf("median_variance_log2_cpm:    %.4f\n", median_var_log2cpm))

# save raw UMI counts (Cyclum/scPrisma both normalize internally; gene-length/TPM correction
# is only appropriate for read-count protocols, not droplet/UMI data -- see README)
mat2hdf(as.matrix(umi_matrix), "myc/data/GSM4286760/HCT116_umi.h5")