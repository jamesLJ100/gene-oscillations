library(here)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(readr)
library(msigdbr)
library(BiocParallel)
library(ggplot2)
library(dorothea)
library(TFEA.ChIP)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

register(SerialParam())  # force single threaded - otherwise fgsea hangs...windows things
proj_root <- here::here()
setwd(proj_root)

source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "common/gsea_functions.R"))

hct116_dir      <- file.path(proj_root, "myc/data/GSM4286760")
input_data_file <- list.files(hct116_dir, pattern = "\\.h5$", full.names = TRUE)[1]
fname           <- tools::file_path_sans_ext(basename(input_data_file))

# Run the full pipeline for each algorithm. Outputs go to a per-algorithm
# subdirectory (myc/figures/<algorithm>/) so the two runs never overwrite
# each other, matching the convention used for nfkb/segmentation_clock.
for (algorithm in c("cyclum", "scPrisma")) {

  cat("\n############################################\n")
  cat("### Algorithm:", algorithm, "\n")
  cat("############################################\n")

  figures_dir <- file.path(proj_root, "myc/figures", algorithm)
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

  results_dir <- file.path(proj_root, "myc/results", algorithm)

  # Results live in myc/results/<algorithm>/ (separate from the raw/preprocessed
  # data), not nested inside hct116_dir the way utils.R's get_model_scores()
  # assumes - so load directly here instead of using that shared helper.
  if (algorithm == "cyclum") {
    weight_file <- file.path(results_dir, paste0(fname, ".h5"))
    results_df  <- if (file.exists(weight_file)) get_cyclum_scores(input_data_file, weight_file) else NULL
  } else if (algorithm == "scPrisma") {
    results_file <- file.path(results_dir, paste0(fname, ".csv"))
    results_df   <- if (file.exists(results_file)) read.csv(results_file) else NULL
  } else {
    stop("Unknown algorithm: ", algorithm)
  }
  if (is.null(results_df)) {
    cat("No", algorithm, "results found for", fname, "in", results_dir, "- skipping.\n")
    next
  }

  geneList <- build_ranked_gene_list(results_df, org_db = org.Hs.eg.db)

  pdf(file.path(figures_dir, paste0(fname, "_score_distribution.pdf")))
  hist(geneList, breaks = 50, main = paste(algorithm, "score distribution -", fname),
       xlab = paste(algorithm, "magnitude score"))
  dev.off()

  run_pathway_analyses(
    geneList    = geneList,
    species     = "human",
    figures_dir = figures_dir,
    fname       = fname,
    algorithm   = algorithm
  )
}
