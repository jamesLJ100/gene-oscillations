source(file.path(proj_root, "hdfrw.R"))


# Download NFKB data
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
library(data.table)

#suppl_files <- getGEOSuppFiles("GSE162992", baseDir = "data/")
#untar("data/GSE162992/GSE162992_RAW.tar", exdir = "data/GSE162992/GSE162992_RAW")
#files <- list.files("data/GSE162992/GSE162992_RAW", full.names = TRUE)
#print(files)

if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
library(Seurat)


# Create directories
dir.create("data/GSE162992/GSE162992_RAW/WT", showWarnings = FALSE, recursive = TRUE)
dir.create("data/GSE162992/GSE162992_RAW/IkBaMM", showWarnings = FALSE, recursive = TRUE)

# Copy WT files
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_barcodes.tsv.gz",
          "data/GSE162992/GSE162992_RAW/WT/barcodes.tsv.gz", overwrite = TRUE)
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_features.tsv.gz",
          "data/GSE162992/GSE162992_RAW/WT/features.tsv.gz", overwrite = TRUE)
file.copy("data/GSE162992/GSE162992_RAW/GSM4969860_WT_matrix.mtx.gz",
          "data/GSE162992/GSE162992_RAW/WT/matrix.mtx.gz", overwrite = TRUE)

# Copy IkBaMM (SS) files
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_barcodes.tsv.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/barcodes.tsv.gz", overwrite = TRUE)
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_features.tsv.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/features.tsv.gz", overwrite = TRUE)
file.copy("data/GSE162992/GSE162992_RAW/GSM4969861_IkBaMM_matrix.mtx.gz",
          "data/GSE162992/GSE162992_RAW/IkBaMM/matrix.mtx.gz", overwrite = TRUE)

# Load data
wt     <- Read10X(data.dir = "data/GSE162992/GSE162992_RAW/WT")
ss     <- Read10X(data.dir = "data/GSE162992/GSE162992_RAW/IkBaMM")

# ── Function to process a single mouse type ───────────────────────────────────

process_mouse_type <- function(data, mouse_type, prefix_pattern) {
  
  cat("\n========================================\n")
  cat("Processing:", mouse_type, "\n")
  cat("========================================\n")
  
  # ── Step 1: Create Seurat object and filter by minimum features ─────────────
  
  seurat_obj <- CreateSeuratObject(
    counts       = data$`Gene Expression`,
    project      = mouse_type,
    min.features = 1500
  )
  
  cat("Cells after >=1500 feature filter:", ncol(seurat_obj), "\n")
  
  # ── Step 2: Calculate QC metrics and filter ─────────────────────────────────
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  seurat_obj <- subset(
    seurat_obj, 
    subset = nFeature_RNA > 200 & percent.mt < 20
  )
  
  cat("Cells after QC filtering:", ncol(seurat_obj), "\n")
  
  # ── Step 3: Assign HTO labels from antibody capture ─────────────────────────
  
  ab       <- data$`Antibody Capture`
  ab_dense <- as.matrix(ab)
  
  shared_barcodes <- intersect(colnames(seurat_obj), colnames(ab_dense))
  ab_filtered     <- ab_dense[, shared_barcodes]
  
  cat("Shared barcodes between GEX and HTO:", length(shared_barcodes), "\n")
  
  total_per_cell <- colSums(ab_filtered)
  safe_totals    <- ifelse(total_per_cell == 0, 1, total_per_cell)
  
  ab_proportions <- sweep(ab_filtered, 2, safe_totals, FUN = "/")
  
  max_ab <- vapply(seq_len(ncol(ab_proportions)), function(i) {
    col <- ab_proportions[, i]
    if (all(col == 0)) return(NA_integer_)
    which.max(col)
  }, integer(1))
  
  max_prop <- vapply(seq_len(ncol(ab_proportions)), function(i) {
    col <- ab_proportions[, i]
    if (all(col == 0)) return(0)
    max(col)
  }, numeric(1))
  
  # Remove prefix pattern from antibody names
  clean_names <- gsub(prefix_pattern, "", rownames(ab_filtered))
  
  cell_labels <- ifelse(
    !is.na(max_ab) & max_prop >= 0.75,
    clean_names[max_ab],
    "Unassigned"
  )
  names(cell_labels) <- shared_barcodes
  
  cat("\nCell label distribution:\n")
  print(table(cell_labels))
  
  # ── Step 4: Add labels to Seurat object ─────────────────────────────────────
  
  seurat_obj$HTO_classification <- cell_labels[colnames(seurat_obj)]
  
  cat("\nSeurat metadata label distribution:\n")
  print(table(seurat_obj$HTO_classification))
  
  assigned_obj <- subset(seurat_obj, HTO_classification != "Unassigned")
  
  cat("\nFinal assigned cells per condition:\n")
  print(table(assigned_obj$HTO_classification))
  
  return(assigned_obj)
}

# ── Process both mouse types ──────────────────────────────────────────────────

wt_assigned <- process_mouse_type(wt, "WT", "^WT_|_TotalSeqB$")
ss_assigned <- process_mouse_type(ss, "IkBaMM", "^IkBaMM_|_TotalSeqB$")

# ── Step 5: Save to HDF5 files by condition ──────────────────────────────────

save_by_condition <- function(seurat_obj, mouse_type_label) {
  
  # Create output directory
  dir.create("data/GSE162992/processed", showWarnings = FALSE, recursive = TRUE)
  
  # Get unique conditions
  conditions <- unique(seurat_obj$HTO_classification)
  conditions <- conditions[conditions != "Unassigned"]
  
  cat("\n========================================\n")
  cat("Saving files for:", mouse_type_label, "\n")
  cat("========================================\n")
  
  for (condition in conditions) {
    # Subset cells for this condition
    subset_obj <- subset(seurat_obj, HTO_classification == condition)
    
    # Extract the count matrix
    counts_mat <- GetAssayData(subset_obj, assay = "RNA", layer = "counts")
    
    # Convert to dense matrix if sparse
    if (inherits(counts_mat, "sparseMatrix")) {
      counts_mat <- as.matrix(counts_mat)
    }
    
    # Clean condition name (remove MM_ prefix if present)
    clean_condition <- gsub("^MM_", "", condition)
    
    # Create filename
    output_file <- sprintf("data/GSE162992/processed/%s_%s.h5", 
                           tolower(mouse_type_label), 
                           tolower(clean_condition))
    
    # Save to HDF5
    mat2hdf(counts_mat, output_file)
    
    cat(sprintf("Saved %s_%s: %d genes x %d cells\n", 
                mouse_type_label, clean_condition, 
                nrow(counts_mat), ncol(counts_mat)))
  }
}

# Save files for both mouse types
save_by_condition(wt_assigned, "wt")
save_by_condition(ss_assigned, "ss")

cat("\n========================================\n")
cat("Processing complete!\n")
cat("========================================\n")