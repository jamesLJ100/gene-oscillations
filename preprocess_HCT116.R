library(data.table)
library(biomaRt)
library(rhdf5)

# read the file
umi_matrix <- fread("data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")
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

# get gene length info from ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("hgnc_symbol", "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = rownames(umi_matrix),
  mart       = mart
)

cat("\n=== biomaRt results ===\n")
cat("Genes queried:", nrow(umi_matrix), "\n")
cat("Genes returned by biomaRt:", nrow(gene_info), "\n")
cat("Duplicated symbols:", sum(duplicated(gene_info$hgnc_symbol)), "\n")

# remove duplicated gene symbols
gene_info <- gene_info[!duplicated(gene_info$hgnc_symbol), ]
cat("After dedup:", nrow(gene_info), "genes\n")

# calculate gene length in kb
gene_info$length_kb <- abs(gene_info$end_position - gene_info$start_position) / 1000

# === missing gene analysis ===
cat("\n=== Missing gene analysis ===\n")
missing_genes <- gene_info$hgnc_symbol[!gene_info$hgnc_symbol %in% rownames(umi_matrix)]
cat("Genes in gene_info but not in umi_matrix:", length(missing_genes), "\n")
if (length(missing_genes) > 0) {
  for (g in missing_genes) {
    cat("\nMissing gene:", g, "\n")
    cat("  nchar:", nchar(g), "\n")
    cat("  chartr check (trimmed):", trimws(g), "\n")
    cat("  identical to trimmed:", identical(g, trimws(g)), "\n")
    
    case_match <- rownames(umi_matrix)[tolower(rownames(umi_matrix)) == tolower(g)]
    cat("  Case-insensitive match in umi_matrix:",
        if (length(case_match) > 0) case_match else "none", "\n")
    
    partial_match <- rownames(umi_matrix)[grepl(paste0("^", g, "$"), rownames(umi_matrix), ignore.case = TRUE)]
    cat("  Partial match in umi_matrix:",
        if (length(partial_match) > 0) partial_match else "none", "\n")
  }
}

not_in_biomart <- rownames(umi_matrix)[!rownames(umi_matrix) %in% gene_info$hgnc_symbol]
cat("\nGenes in umi_matrix not returned by biomaRt:", length(not_in_biomart), "\n")
cat("First few:", head(not_in_biomart), "\n")

# remove missing genes from gene_info before subsetting
gene_info <- gene_info[gene_info$hgnc_symbol %in% rownames(umi_matrix), ]
cat("\nFinal gene_info after removing missing genes:", nrow(gene_info), "\n")

# subset matrix to genes we have length info for
umi_sub <- umi_matrix[gene_info$hgnc_symbol, ]
cat("\n=== Subsetting matrix ===\n")
cat("Dimensions of umi_sub:", dim(umi_sub), "\n")
cat("NAs in umi_sub:", sum(is.na(umi_sub)), "\n")

# === Dataset summary statistics (post-gene-filter) ===
cat("\n=== Dataset Summary Statistics (post-filter) ===\n")

num_features_sub        <- nrow(umi_sub)
num_cells_sub           <- ncol(umi_sub)
lib_sizes_sub           <- colSums(umi_sub)
median_library_size_sub <- median(lib_sizes_sub)
dropout_rate_sub        <- sum(umi_sub == 0) / (num_features_sub * num_cells_sub)

cpm_sub     <- t(t(umi_sub) / lib_sizes_sub) * 1e6
log2cpm_sub <- log2(cpm_sub + 1)

gene_mean_log2cpm_sub  <- rowMeans(log2cpm_sub)
median_avg_log2cpm_sub <- median(gene_mean_log2cpm_sub)

if (requireNamespace("matrixStats", quietly = TRUE)) {
  gene_var_log2cpm_sub <- matrixStats::rowVars(as.matrix(log2cpm_sub))
} else {
  gene_var_log2cpm_sub <- apply(log2cpm_sub, 1, var)
}
median_var_log2cpm_sub <- median(gene_var_log2cpm_sub)

cat(sprintf("num_cells:                   %d\n",   num_cells_sub))
cat(sprintf("num_features:                %d\n",   num_features_sub))
cat(sprintf("median_library_size:         %.1f\n", median_library_size_sub))
cat(sprintf("dropout_rate:                %.4f\n", dropout_rate_sub))
cat(sprintf("median_average_log2_cpm:     %.4f\n", median_avg_log2cpm_sub))
cat(sprintf("median_variance_log2_cpm:    %.4f\n", median_var_log2cpm_sub))

# convert to TPM
rpk <- umi_sub / gene_info$length_kb
tpm <- t(t(rpk) / colSums(rpk)) * 1e6

cat("\n=== TPM ===\n")
cat("Dimensions:", dim(tpm), "\n")
cat("NAs in tpm:", sum(is.na(tpm)), "\n")
cat("Sample values:\n")
print(tpm[1:5, 1:5])

# save
mat2hdf(tpm, "data/GSM4286760/HCT116_tpm.h5")