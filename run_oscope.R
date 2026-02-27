#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(Oscope)
})

set.seed(42)

# ===== OSCOPE FUNCTION =====
apply_oscope <- function(gene_expr_matrix, maxK = 10, NCThre = 100, top_pairs_quantile = 0.05) {
  # gene_expr_matrix should be genes x samples
  # top_pairs_quantile: threshold for selecting top gene pairs (default 5%)
  
  gene_expr_matrix <- as.matrix(gene_expr_matrix)
  
  # 1. Normalization: Calculate size factors
  if (sum(is.na(MedianNorm(gene_expr_matrix))) != 0) {
    Sizes <- MedianNorm(gene_expr_matrix, alternative = TRUE)
  } else {
    Sizes <- MedianNorm(gene_expr_matrix)
  }
  
  # Get normalized expression matrix
  DataNorm <- GetNormalizedMat(gene_expr_matrix, Sizes)
  
  # 2. Pre-processing: Filter for high mean/variance genes
  MV <- CalcMV(Data = gene_expr_matrix, Sizes = Sizes, NormData = FALSE, MeanCutLow = 0)
  
  if (length(MV$GeneToUse) == 0) {
    warning("No genes passed mean-variance filtering")
    return(NULL)
  }
  
  # Subset to high mean/high variance genes
  DataSubset <- DataNorm[MV$GeneToUse, , drop = FALSE]
  
  # 3. Rescaling for sine model
  DataInput <- NormForSine(DataSubset, qt1 = 0.05, qt2 = 0.95)
  
  # Remove incomplete cases
  DataInput_subset <- as.matrix(DataInput[complete.cases(DataInput), ])
  
  if (nrow(DataInput_subset) == 0) {
    warning("No complete cases after rescaling")
    return(NULL)
  }
  
  # 4. Run paired-sine model
  cat("  Running OscopeSine...\n")
  SineRes <- OscopeSine(DataInput_subset)
  
  # 5. IDENTIFY CANDIDATE OSCILLATORY GENES
  # From the paper: "candidate oscillatory genes are those genes in the top gene pairs"
  # Gene pairs are ranked by -log10(ε²), which is the sine score (SimiMat)
  
  # Get upper triangle of similarity matrix (to avoid double counting pairs)
  n_genes <- nrow(SineRes$SimiMat)
  upper_tri_indices <- which(upper.tri(SineRes$SimiMat), arr.ind = TRUE)
  
  # Extract sine scores for all gene pairs
  pair_scores <- data.frame(
    gene1 = rownames(SineRes$SimiMat)[upper_tri_indices[, 1]],
    gene2 = rownames(SineRes$SimiMat)[upper_tri_indices[, 2]],
    sine_score = SineRes$SimiMat[upper_tri_indices],
    stringsAsFactors = FALSE
  )
  
  # Determine threshold for top pairs
  threshold <- quantile(pair_scores$sine_score, probs = 1 - top_pairs_quantile, na.rm = TRUE)
  top_pairs <- pair_scores[pair_scores$sine_score >= threshold, ]
  
  # Genes appearing in top pairs are candidate oscillatory genes
  candidate_genes <- unique(c(top_pairs$gene1, top_pairs$gene2))
  
  cat(sprintf("  Identified %d candidate oscillatory genes from top %.1f%% of gene pairs\n", 
              length(candidate_genes), top_pairs_quantile * 100))
  cat(sprintf("  Sine score threshold: %.3f\n", threshold))
  
  # 6. K-medoids clustering (only on candidate genes)
  cat("  Running OscopeKM...\n")
  KMRes <- OscopeKM(SineRes, maxK = maxK)
  
  # 7. Flag clusters with low sine scores or lack of phase shifts
  if (length(KMRes) > 0) {
    cat("  Running FlagCluster...\n")
    ToRM <- FlagCluster(SineRes, KMRes, DataInput_subset)
    
    # Remove flagged clusters
    if (length(ToRM$FlagID) > 0) {
      cat(sprintf("  Removing %d flagged cluster(s) (linear relationships without phase shifts)\n", 
                  length(ToRM$FlagID)))
      KMResUse <- KMRes[-ToRM$FlagID]
    } else {
      KMResUse <- KMRes
    }
  } else {
    warning("No gene clusters identified by K-medoids")
    KMResUse <- list()
    ToRM <- list(FlagID = integer(0))
  }
  
  # 8. Extended nearest insertion to recover cell order (if we have valid clusters)
  ENIRes <- NULL
  clustered_oscillating_genes <- character(0)
  
  if (length(KMResUse) > 0) {
    cat("  Running OscopeENI...\n")
    ENIRes <- OscopeENI(KMRes = KMResUse, Data = DataInput_subset, NCThre = NCThre)
    
    # Genes in non-flagged clusters
    clustered_oscillating_genes <- unique(unlist(KMResUse))
    
    cat(sprintf("  %d genes in %d non-flagged cluster(s)\n", 
                length(clustered_oscillating_genes), length(KMResUse)))
  }
  
  # 9. Create gene classification table
  gene_classification <- data.frame(
    gene = rownames(gene_expr_matrix),
    is_candidate_oscillator = rownames(gene_expr_matrix) %in% candidate_genes,
    in_valid_cluster = rownames(gene_expr_matrix) %in% clustered_oscillating_genes,
    cluster = NA,
    stringsAsFactors = FALSE
  )
  
  # Assign cluster membership
  if (length(KMResUse) > 0) {
    for (i in seq_along(KMResUse)) {
      cluster_genes <- KMResUse[[i]]
      gene_classification$cluster[gene_classification$gene %in% cluster_genes] <- i
    }
  }
  
  # Calculate mean sine score for each gene (for ranking/reference)
  mean_sine_scores <- rowMeans(SineRes$SimiMat, na.rm = TRUE)
  gene_classification$mean_sine_score <- NA
  gene_classification$mean_sine_score[match(names(mean_sine_scores), 
                                            gene_classification$gene)] <- mean_sine_scores
  
  return(list(
    candidate_oscillating_genes = candidate_genes,  # From top gene pairs
    clustered_oscillating_genes = clustered_oscillating_genes,  # In non-flagged clusters
    gene_classification = gene_classification,
    top_pairs = top_pairs,
    top_pairs_threshold = threshold,
    clusters = KMResUse,
    cell_orders = ENIRes,
    flagged_clusters = ToRM$FlagID,
    sine_results = SineRes,
    filtered_genes = MV$GeneToUse,
    normalized_data = DataNorm
  ))
}

# ===== READ DATA FILE =====
read_data_file <- function(file_path) {
  # Read H5 file
  # Assuming data is stored in a dataset, adjust path as needed
  h5_contents <- hdf2mat(file_path)
  cat("H5 file contents:\n")
  print(h5_contents)
  
  
  # If there are row/col names, read them too
  # row_names <- h5read(file_path, "row_names")  # adjust as needed
  # col_names <- h5read(file_path, "col_names")  # adjust as needed
  
  # Convert to matrix if needed
  if (!is.matrix(df)) {
    df <- as.matrix(df)
  }
  
  return(df)
}

# ===== PROCESS SINGLE FILE =====
process_data_file <- function(file_path) {
  cat(sprintf("Processing file: %s\n", file_path))
  
  tryCatch({
    # Read data file
    df <- read_data_file(file_path)
    
    cat(sprintf("Matrix shape: %d x %d\n", nrow(df), ncol(df)))
    cat("Note: Matrix should be genes x samples for Oscope\n")
    
    # Start timing
    start_time <- Sys.time()
    
    # Apply Oscope
    result <- apply_oscope(df)
    
    # End timing
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Save results
    filepath <- dirname(file_path)
    fullflname <- basename(file_path)
    fname <- tools::file_path_sans_ext(fullflname)
    
    output_dir <- file.path(filepath, 'oscope')
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Save full results as RDS file
    output_path <- file.path(output_dir, sprintf('%s_oscope.rds', fname))
    saveRDS(result, output_path)
    cat(sprintf("Results saved to: %s\n", output_path))
    
    # Save gene classification as CSV
    if (!is.null(result) && !is.null(result$gene_classification)) {
      csv_path <- file.path(output_dir, sprintf('%s_gene_classification.csv', fname))
      
      # Sort by candidate status, cluster membership, then sine score
      gene_df <- result$gene_classification
      gene_df <- gene_df[order(-gene_df$is_candidate_oscillator,
                               -gene_df$in_valid_cluster,
                               -gene_df$mean_sine_score, 
                               na.last = TRUE), ]
      
      write.csv(gene_df, csv_path, row.names = FALSE)
      cat(sprintf("Gene classification saved to: %s\n", csv_path))
      
      # Print summary
      n_candidates <- sum(gene_df$is_candidate_oscillator, na.rm = TRUE)
      n_clustered <- sum(gene_df$in_valid_cluster, na.rm = TRUE)
      n_total <- nrow(gene_df)
      cat(sprintf("  --> %d / %d genes are candidate oscillators (%.1f%%) [from top gene pairs]\n", 
                  n_candidates, n_total, 100 * n_candidates / n_total))
      cat(sprintf("  --> %d / %d genes in valid clusters (%.1f%%) [clustered candidates with phase shifts]\n", 
                  n_clustered, n_total, 100 * n_clustered / n_total))
    }
    
    # Save candidate oscillating genes list
    if (!is.null(result) && length(result$candidate_oscillating_genes) > 0) {
      genes_path <- file.path(output_dir, sprintf('%s_candidate_oscillating_genes.txt', fname))
      writeLines(result$candidate_oscillating_genes, genes_path)
      cat(sprintf("Candidate oscillating genes saved to: %s\n", genes_path))
    }
    
    # Save clustered oscillating genes list
    if (!is.null(result) && length(result$clustered_oscillating_genes) > 0) {
      genes_path <- file.path(output_dir, sprintf('%s_clustered_oscillating_genes.txt', fname))
      writeLines(result$clustered_oscillating_genes, genes_path)
      cat(sprintf("Clustered oscillating genes saved to: %s\n", genes_path))
    }
    
    cat(sprintf("Runtime: %.2f seconds\n", elapsed))
    
    return(data.frame(
      file = fname,
      runtime_seconds = elapsed,
      n_genes = nrow(df),
      n_samples = ncol(df),
      n_candidate_oscillating = if (!is.null(result)) length(result$candidate_oscillating_genes) else 0,
      n_clustered_oscillating = if (!is.null(result)) length(result$clustered_oscillating_genes) else 0,
      n_clusters = if (!is.null(result) && !is.null(result$clusters)) length(result$clusters) else 0,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat(sprintf("Error processing file %s: %s\n", file_path, e$message))
    return(data.frame(
      file = basename(file_path),
      runtime_seconds = NA,
      n_genes = NA,
      n_samples = NA,
      n_oscillating = NA,
      n_clusters = NA,
      error = e$message,
      stringsAsFactors = FALSE
    ))
  })
}

# ===== MAIN FUNCTION =====
main <- function(input_dir) {
  args <- commandArgs(trailingOnly = TRUE)
  
  input_dir_abs <- normalizePath(input_dir, mustWork = FALSE)
  
  cat(sprintf("Current working directory: %s\n", getwd()))
  
  if (!dir.exists(input_dir_abs)) {
    cat(sprintf("Directory does not exist: %s\n", input_dir_abs))
    quit(status = 1)
  }
  
  cat(sprintf("Directory contents: %s\n", 
              paste(list.files(input_dir_abs), collapse = ", ")))
  
  # Find all H5 files
  data_files <- list.files(
    input_dir_abs, 
    pattern = "\\.h5$", 
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  # Exclude files ending with _sim
  data_files <- data_files[!grepl("_sim\\.h5$", data_files, ignore.case = TRUE)]
  
  cat(sprintf("Found %d H5 files to process\n", length(data_files)))
  
  if (length(data_files) == 0) {
    cat("No H5 files found.\n")
    quit(status = 1)
  }
  
  # Process each file
  timing_records <- list()
  
  for (idx in seq_along(data_files)) {
    file_path <- data_files[idx]
    cat(sprintf("\n=== Processing file %d/%d: %s ===\n", 
                idx, length(data_files), basename(file_path)))
    
    record <- process_data_file(file_path)
    timing_records[[idx]] <- record
  }
  
  # Combine timing records
  timing_df <- do.call(rbind, timing_records)
  
  # Save timing results
  output_dir <- file.path(input_dir_abs, 'oscope')
  csv_path <- file.path(output_dir, 'runtimes.csv')
  write.csv(timing_df, csv_path, row.names = FALSE)
  
  cat(sprintf("\n=== COMPLETED ===\n"))
  cat(sprintf("Runtimes saved to: %s\n", csv_path))
  cat("\nSummary:\n")
  print(timing_df)
}

proj_root <- here::here()
setwd(proj_root)
dyngen_dir  <- file.path(proj_root, "data/dyngen")

main(dyngen_dir)


