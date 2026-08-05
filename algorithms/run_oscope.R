library(Oscope)
library(hdf5r)

source(file.path(proj_root, "common/hdfrw.R"))

apply_oscope <- function(gene_expr_matrix, oscope_config) {

  KM_maxK           <- oscope_config$KM_maxK
  KM_quan           <- oscope_config$KM_quan
  CalcMV_MeanCutLow <- oscope_config$CalcMV_MeanCutLow
  FlagCluster_qt    <- oscope_config$FlagCluster_qt
  FlagCluster_thre  <- oscope_config$FlagCluster_thre

  cat(sprintf("Config: KM_maxK=%d, KM_quan=%.2f, MeanCutLow=%.3f, FlagCluster_qt=%.2f, FlagCluster_thre=%.3f\n",
              KM_maxK, KM_quan, CalcMV_MeanCutLow, FlagCluster_qt, FlagCluster_thre))

  gene_expr_matrix <- as.matrix(gene_expr_matrix)
  all_genes <- rownames(gene_expr_matrix)
  cat(sprintf("Input matrix: %d genes x %d samples\n",
              nrow(gene_expr_matrix), ncol(gene_expr_matrix)))

  # Cross-sample library-size normalization, always applied (no opt-out): MedianNorm()
  # reproduces DESeq's median-of-ratios size factors, and GetNormalizedMat() applies them.
  # This is Oscope's own recommended normalization (its vignette section 2.2), operating
  # directly on raw counts with no gene-length correction - appropriate for our raw UMI
  # input, the same way Cyclum/scPrisma are fed raw counts and normalize internally.
  #
  # Caveat: Oscope (Leng et al. 2015, Nature Methods) was developed and validated on the
  # Fluidigm C1 platform - full-length, plate-based, non-UMI, predating droplet protocols.
  # Correct normalization doesn't resolve this: Oscope's sine-fitting, its CalcMV_MeanCutLow
  # dropout threshold, and its K-medoids clustering were never validated against droplet
  # UMI data's much sparser, more zero-inflated profile. Treat Oscope's results on our
  # datasets with more caution than Cyclum's/scPrisma's, whose own papers explicitly
  # address droplet/UMI applicability.
  Sizes      <- MedianNorm(gene_expr_matrix)
  DataNorm   <- GetNormalizedMat(gene_expr_matrix, Sizes)
  MV         <- CalcMV(Data = DataNorm, Sizes = NULL, NormData = TRUE, MeanCutLow = CalcMV_MeanCutLow)
  DataSubset <- DataNorm[MV$GeneToUse,]

  # # Select high mean / high variance genes for sine model input
  # MV <- CalcMV(Data = gene_expr_matrix, Sizes = NULL, NormData = TRUE, MeanCutLow = CalcMV_MeanCutLow)
  # cat(sprintf("CalcMV: %d genes suggested (GeneToUse) out of %d total\n",
  #             length(MV$GeneToUse), length(all_genes)))
  # 
  # if (length(MV$GeneToUse) == 0) {
  #   warning("No genes passed mean-variance filtering")
  #   return(list(gene_classification = NULL, SineRes = NULL))
  # }
  # cat(sprintf("  %d / %d genes passed MV filter\n", length(MV$GeneToUse), nrow(DataNorm)))
  
  # Rescale to [-1, 1] for sine model
  DataInput <- NormForSine(DataSubset)
  # DataInput_full <- as.matrix(DataNormScaled[MV$GeneToUse, , drop = FALSE])
  # cat(sprintf("After MV, before complete.cases: %d genes\n", nrow(DataInput_full)))
  # 
  # DataInput <- DataInput_full[complete.cases(DataInput_full), , drop = FALSE]
  # cat(sprintf("After complete.cases: %d genes (dropped %d)\n",
  #             nrow(DataInput), nrow(DataInput_full) - nrow(DataInput)))
  # 
  # if (nrow(DataInput) < 2) {
  #   warning(sprintf("Too few genes after filtering: %d", nrow(DataInput)))
  #   return(list(gene_classification = NULL, SineRes = NULL))
  # }
  # cat(sprintf("  %d genes passed all filters for sine model\n", nrow(DataInput)))
  
  # Paired-sine model
  cat("Running OscopeSine...\n")
  SineRes <- OscopeSine(DataInput)
  cat("OscopeSine completed.\n")
  
  # K-medoids clustering
  cat(sprintf("Running OscopeKM (maxK=%d, quan=%.2f)...\n", KM_maxK, KM_quan))
  KMRes <- tryCatch(
    OscopeKM(SineRes, maxK = KM_maxK, quan = KM_quan),
    error = function(e) {
      cat(sprintf("OscopeKM failed (likely too few gene pairs): %s\n", e$message))
      return(NULL)
    }
  )
  cat(sprintf("OscopeKM returned %d cluster candidate sets (length(KMRes))\n",
              length(KMRes)))
  
  clustered_oscillating_genes <- character(0)
  clusters  <- list()
  KMResUse  <- list()
  
  if (is.null(KMRes) || length(KMRes) == 0) {
    warning("No clusters found from OscopeKM")
  } else {
    # Flag clusters with small within-cluster phase shift
    cat(sprintf("Running FlagCluster (qt=%.2f, thre=%.3f)...\n",
                FlagCluster_qt, FlagCluster_thre))
    ToRM <- FlagCluster(SineRes, KMRes, DataInput,
                        qt = FlagCluster_qt, thre = FlagCluster_thre)
    
    cat(sprintf("  Clusters before flagging: %d\n", length(KMRes)))
    cat("  FlagID from FlagCluster: ", paste(ToRM$FlagID, collapse = ", "), "\n")
    
    for (k in seq_along(KMRes)) {
      flagged <- k %in% ToRM$FlagID
      cat(sprintf("  Cluster %d: %d genes [%s]\n",
                  k, length(KMRes[[k]]), if (flagged) "FLAGGED" else "kept"))
    }
    
    if (length(ToRM$FlagID) > 0) {
      cat(sprintf("  Removing %d flagged cluster(s)\n", length(ToRM$FlagID)))
      KMResUse <- KMRes[-ToRM$FlagID]
    } else {
      KMResUse <- KMRes
    }
    
    cat(sprintf("  Clusters after flagging: %d\n", length(KMResUse)))
    
    if (length(KMResUse) > 0) {
      clustered_oscillating_genes <- unique(unlist(KMResUse))
      clusters <- KMResUse
      cat(sprintf("  Found %d oscillating genes in %d cluster(s)\n",
                  length(clustered_oscillating_genes), length(KMResUse)))
    } else {
      warning("No valid clusters after flagging")
    }
  }
  
  # Build per-gene classification table
  n_osc <- sum(all_genes %in% clustered_oscillating_genes)
  cat(sprintf("Final oscillating genes: %d\n", n_osc))
  
  gene_classification <- data.frame(
    gene           = all_genes,
    is_oscillating = all_genes %in% clustered_oscillating_genes,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  return(list(
    gene_classification = gene_classification,
    SineRes = SineRes
  ))
}


# ===== READ DATA FILE =====
read_data_file <- function(file_path) {
  matrix <- hdf2mat(file_path)
  cat("H5 file contents:\n")
  return(matrix)
}

# ===== PROCESS SINGLE FILE =====
process_data_file <- function(file_path, oscope_config, run_number = NULL) {
  cat(sprintf("Processing file: %s\n", file_path))
  
  tryCatch({
    df <- read_data_file(file_path)
    cat(sprintf("Matrix shape: %d x %d (genes x samples)\n", nrow(df), ncol(df)))
    
    start_time <- Sys.time()
    result     <- apply_oscope(df, oscope_config)
    elapsed    <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    filepath   <- dirname(file_path)
    fname      <- tools::file_path_sans_ext(basename(file_path))
    output_dir <- file.path(filepath, 'oscope')
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Create run number suffix if provided (inserted before tool name)
    if (!is.null(run_number)) {
      run_suffix <- sprintf("_r%d", run_number)
    } else {
      run_suffix <- ""
    }
    
    gene_classification <- result$gene_classification
    if (!is.null(gene_classification)) {
      csv_path <- file.path(output_dir, sprintf('%s%s.csv', fname, run_suffix))
      
      gene_df <- gene_classification[, c("gene", "is_oscillating")]
      colnames(gene_df) <- c("symbol", "score")
      gene_df$score <- as.integer(gene_df$score)  # TRUE/FALSE -> 1/0
      gene_df <- gene_df[order(-gene_df$score), ]
      
      write.csv(gene_df, csv_path, row.names = FALSE)
      cat(sprintf("Gene classification saved to: %s\n", csv_path))
      
      n_oscillating <- sum(gene_df$score, na.rm = TRUE)
      n_total       <- nrow(gene_df)
      cat(sprintf("  --> %d / %d genes are oscillating (%.1f%%)\n",
                  n_oscillating, n_total, 100 * n_oscillating / n_total))
    }
    
    return(data.frame(
      file            = fname,
      runtime_seconds = elapsed,
      n_genes         = nrow(df),
      n_samples       = ncol(df),
      n_oscillating   = if (!is.null(gene_classification)) sum(gene_classification$is_oscillating, na.rm = TRUE) else 0,
      error           = NA_character_,
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    cat(sprintf("Error processing file %s: %s\n", file_path, e$message))
    return(data.frame(
      file            = basename(file_path),
      runtime_seconds = NA_real_,
      n_genes         = NA_integer_,
      n_samples       = NA_integer_,
      n_oscillating   = NA_integer_,
      error           = e$message,
      stringsAsFactors = FALSE
    ))
  })
}

#' Run Oscope on every .h5 file in a directory
#'
#' Unlike Cyclum/scPrisma (separate Python processes invoked via reticulate), Oscope
#' is an R/Bioconductor package, so this runs in-process rather than shelling out.
#'
#' @param input_dir Directory containing input .h5 count files (non-`_sim.h5` files only)
#' @param oscope_config List with elements KM_maxK, KM_quan, CalcMV_MeanCutLow,
#'   FlagCluster_qt, FlagCluster_thre. Cross-sample library-size normalization is
#'   always applied internally (not configurable) - see apply_oscope() for why.
#' @param run_number Optional integer appended to output filenames (for grid search / repeat runs)
#' @return Invisibly NULL; per-file gene classifications are written to
#'   <input_dir>/oscope/, and a combined runtimes.csv summarising all files
# ===== MAIN FUNCTION =====
run_oscope <- function(input_dir, oscope_config, run_number = NULL) {
  input_dir_abs <- normalizePath(input_dir, mustWork = FALSE)
  
  cat(sprintf("Current working directory: %s\n", getwd()))
  
  if (!dir.exists(input_dir_abs)) {
    stop("Directory does not exist: ", input_dir_abs)
  }

  cat(sprintf("Directory contents: %s\n",
              paste(list.files(input_dir_abs), collapse = ", ")))

  data_files <- list.files(
    input_dir_abs,
    pattern = "\\.h5$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  data_files <- data_files[!grepl("_sim\\.h5$", data_files, ignore.case = TRUE)]

  cat(sprintf("Found %d H5 files to process\n", length(data_files)))
  if (length(data_files) == 0) {
    stop("No H5 files found in: ", input_dir_abs)
  }
  
  # Print run_number info if provided
  if (!is.null(run_number)) {
    cat(sprintf("Run number: %d\n", run_number))
  }
  
  timing_records <- list()
  for (idx in seq_along(data_files)) {
    file_path <- data_files[idx]
    cat(sprintf("\n=== Processing file %d/%d: %s ===\n",
                idx, length(data_files), basename(file_path)))
    timing_records[[idx]] <- process_data_file(file_path, oscope_config, run_number)
  }
  
  timing_df <- do.call(rbind, timing_records)
  
  output_dir <- file.path(input_dir_abs, 'oscope')
  csv_filename <- 'runtimes.csv'
  csv_path <- file.path(output_dir, csv_filename)
  write.csv(timing_df, csv_path, row.names = FALSE)
  
  cat(sprintf("\n=== COMPLETED ===\n"))
  cat(sprintf("Runtimes saved to: %s\n", csv_path))
  cat("\nSummary:\n")
  print(timing_df)
}

