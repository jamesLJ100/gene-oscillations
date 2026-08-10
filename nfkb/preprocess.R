library(data.table)
library(Seurat)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "common/hdfrw.R"))

raw_dir <- file.path(proj_root, "nfkb/data/GSE162992/GSE162992_RAW")

required_files <- file.path(raw_dir, c(
  "GSM4969860_WT_barcodes.tsv.gz", "GSM4969860_WT_features.tsv.gz", "GSM4969860_WT_matrix.mtx.gz",
  "GSM4969861_IkBaMM_barcodes.tsv.gz", "GSM4969861_IkBaMM_features.tsv.gz", "GSM4969861_IkBaMM_matrix.mtx.gz"
))
if (!all(file.exists(required_files))) source(file.path(proj_root, "nfkb/download.R"))

# Copy each GSM's raw files into a per-mouse-type folder named the way Read10X()
# expects (barcodes/features/matrix.*, no GSM/condition prefix).
gsm_by_type <- c(WT = "GSM4969860", IkBaMM = "GSM4969861")
file_ext    <- c(barcodes = "tsv.gz", features = "tsv.gz", matrix = "mtx.gz")

for (mouse_type in names(gsm_by_type)) {
  dir.create(file.path(raw_dir, mouse_type), showWarnings = FALSE, recursive = TRUE)
  for (file_type in names(file_ext)) {
    file.copy(
      file.path(raw_dir, sprintf("%s_%s_%s.%s", gsm_by_type[[mouse_type]], mouse_type,
                                  file_type, file_ext[[file_type]])),
      file.path(raw_dir, mouse_type, sprintf("%s.%s", file_type, file_ext[[file_type]])),
      overwrite = TRUE
    )
  }
}

# Load data
wt     <- Read10X(data.dir = file.path(raw_dir, "WT"))
ss     <- Read10X(data.dir = file.path(raw_dir, "IkBaMM"))


process_mouse_type <- function(data, mouse_type) {

  # e.g. "WT" -> "^WT_|_TotalSeqB$" - derived rather than passed separately so the
  # two can't drift out of sync.
  prefix_pattern <- sprintf("^%s_|_TotalSeqB$", mouse_type)
  
  cat("Processing:", mouse_type, "\n")

  # Create Seurat object and filter by minimum features
  
  seurat_obj <- CreateSeuratObject(
    counts       = data$`Gene Expression`,
    project      = mouse_type,
    min.features = 1500
  )
  
  cat("Cells after >=1500 feature filter:", ncol(seurat_obj), "\n")
  
  # Calculate QC metrics and filter
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  seurat_obj <- subset(
    seurat_obj, 
    subset = nFeature_RNA > 200 & percent.mt < 20
  )
  
  cat("Cells after QC filtering:", ncol(seurat_obj), "\n")
  
  # Assign HTO labels from antibody capture
  
  ab       <- data$`Antibody Capture`
  ab_dense <- as.matrix(ab)
  
  shared_barcodes <- intersect(colnames(seurat_obj), colnames(ab_dense))
  ab_filtered     <- ab_dense[, shared_barcodes]
  
  cat("Shared barcodes between GEX and HTO:", length(shared_barcodes), "\n")
  
  total_per_cell <- colSums(ab_filtered)
  safe_totals    <- ifelse(total_per_cell == 0, 1, total_per_cell)
  
  ab_proportions <- sweep(ab_filtered, 2, safe_totals, FUN = "/")
  
  # max.col() operates row-wise, so transpose to get one row per cell. A cell with
  # no antibody signal at all (all-zero column) still gets an arbitrary index here,
  # but its max_prop is 0, so it's excluded below by the >= 0.75 threshold regardless
  # - no special-casing needed for that case.
  max_ab   <- max.col(t(ab_proportions), ties.method = "first")
  max_prop <- apply(ab_proportions, 2, max)

  # Remove prefix pattern from antibody names
  clean_names <- gsub(prefix_pattern, "", rownames(ab_filtered))

  cell_labels <- ifelse(
    max_prop >= 0.75,
    clean_names[max_ab],
    "Unassigned"
  )
  names(cell_labels) <- shared_barcodes
  
  cat("\nCell label distribution:\n")
  print(table(cell_labels))
  
  # Add labels to Seurat object
  
  seurat_obj$HTO_classification <- cell_labels[colnames(seurat_obj)]
  
  cat("\nSeurat metadata label distribution:\n")
  print(table(seurat_obj$HTO_classification))
  
  assigned_obj <- subset(seurat_obj, HTO_classification != "Unassigned")
  
  cat("\nFinal assigned cells per condition:\n")
  print(table(assigned_obj$HTO_classification))
  
  return(assigned_obj)
}


wt_assigned <- process_mouse_type(wt, "WT")
ss_assigned <- process_mouse_type(ss, "IkBaMM")

# Save to HDF5 files by condition
save_by_condition <- function(seurat_obj, mouse_type_label) {

  processed_dir <- file.path(proj_root, "nfkb/data/GSE162992/processed")
  dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
  
  conditions <- unique(seurat_obj$HTO_classification)
  conditions <- conditions[conditions != "Unassigned"]
  
  cat("Saving files for:", mouse_type_label, "\n")

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
    output_file <- file.path(processed_dir, sprintf("%s_%s.h5",
                                                      tolower(mouse_type_label),
                                                      tolower(clean_condition)))
    
    # Save to HDF5
    mat2hdf(counts_mat, output_file)
    
    cat(sprintf("Saved %s_%s: %d genes x %d cells\n", 
                mouse_type_label, clean_condition, 
                nrow(counts_mat), ncol(counts_mat)))
  }
}

save_by_condition(wt_assigned, "wt")
save_by_condition(ss_assigned, "ss")
