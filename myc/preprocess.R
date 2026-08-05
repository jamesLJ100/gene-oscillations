library(data.table)
library(rhdf5)
library(matrixStats)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))
source(file.path(proj_root, "common/utils.R"))

raw_file <- file.path(proj_root, "myc/data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")
if (!file.exists(raw_file)) source(file.path(proj_root, "myc/download.R"))

umi_matrix <- fread(raw_file)
umi_matrix <- as.data.frame(umi_matrix)
rownames(umi_matrix) <- umi_matrix$V1
umi_matrix$V1 <- NULL
umi_matrix <- as.matrix(umi_matrix)

cat("Dimensions:", dim(umi_matrix), "\n")
cat("NAs in umi_matrix:", sum(is.na(umi_matrix)), "\n")

# Cell QC: drop cells with high mitochondrial fraction or a low UMI count
qc_flag    <- qc_cell_mask(umi_matrix, mito_prefix = "^MT-")
umi_matrix <- umi_matrix[, qc_flag]

num_features        <- nrow(umi_matrix)
num_cells           <- ncol(umi_matrix)
lib_sizes           <- colSums(umi_matrix)
median_library_size <- median(lib_sizes)
dropout_rate        <- sum(umi_matrix == 0) / (num_features * num_cells)

cpm     <- sweep(umi_matrix, 2, lib_sizes, "/") * 1e6
log2cpm <- log2(cpm + 1)

gene_mean_log2cpm  <- rowMeans(log2cpm)
median_avg_log2cpm <- median(gene_mean_log2cpm)

gene_var_log2cpm   <- rowVars(log2cpm)
median_var_log2cpm <- median(gene_var_log2cpm)

cat(sprintf("num_cells:                   %d\n",   num_cells))
cat(sprintf("num_features:                %d\n",   num_features))
cat(sprintf("median_library_size:         %.1f\n", median_library_size))
cat(sprintf("dropout_rate:                %.4f\n", dropout_rate))
cat(sprintf("median_average_log2_cpm:     %.4f\n", median_avg_log2cpm))
cat(sprintf("median_variance_log2_cpm:    %.4f\n", median_var_log2cpm))

# save raw UMI counts
mat2hdf(umi_matrix, "myc/data/GSM4286760/HCT116_umi.h5")