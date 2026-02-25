
# read the file
umi_matrix <- fread("data/GSM4286760/GSM4286760_scRNA_HCT116_WT_10x_umiTable.txt.gz")

# set gene names as rownames
umi_matrix <- as.data.frame(umi_matrix)
rownames(umi_matrix) <- umi_matrix$V1
umi_matrix$V1 <- NULL

cat("=== UMI matrix ===\n")
cat("Dimensions:", dim(umi_matrix), "\n")
cat("NAs in umi_matrix:", sum(is.na(umi_matrix)), "\n")

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
    
    # check if it exists in umi_matrix under different case
    case_match <- rownames(umi_matrix)[tolower(rownames(umi_matrix)) == tolower(g)]
    cat("  Case-insensitive match in umi_matrix:", 
        if (length(case_match) > 0) case_match else "none", "\n")
    
    # check for partial matches
    partial_match <- rownames(umi_matrix)[grepl(paste0("^", g, "$"), rownames(umi_matrix), ignore.case = TRUE)]
    cat("  Partial match in umi_matrix:", 
        if (length(partial_match) > 0) partial_match else "none", "\n")
  }
}

# also check the reverse - genes in umi_matrix not returned by biomaRt
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