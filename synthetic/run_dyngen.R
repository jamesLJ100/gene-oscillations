library(dyngen)
library(tidyverse)
library(here)
library(hdf5r)

## define project root
proj_root <- here::here()
setwd(proj_root)

source(file.path(proj_root, "common/hdfrw.R"))

# Real gene regulatory network / real scRNA-seq count distribution used to ground each
# simulation (rather than an arbitrary/unrealistic one). Not checked into the repo (data/ is
# gitignored) - downloaded automatically here on first use, then cached locally.
realnet_path   <- here::here("synthetic", "data", "realnet.rds")
realcount_path <- here::here("synthetic", "data", "realcount.rds")

if (!file.exists(realnet_path)) {
  cat("Downloading realnet.rds (one-time)...\n")
  download.file(
    "https://github.com/dynverse/dyngen/raw/data_files/regulatorycircuits_26_neuron-associated_cells_cancer.rds",
    destfile = realnet_path, mode = "wb"
  )
}
if (!file.exists(realcount_path)) {
  cat("Downloading realcount.rds (one-time)...\n")
  download.file(
    "https://github.com/dynverse/dyngen/raw/data_files/zenodo_1443566_real_gold_psc-astrocyte-maturation-neuron_sloan.rds",
    destfile = realcount_path, mode = "wb"
  )
}

realnet_matrix <- readRDS(realnet_path)
realcount_data <- readRDS(realcount_path)

# definition
backbone_cycle_simple <- function() {
  module_info <- tribble(
    ~module_id, ~basal, ~burn, ~independence,
    "A", 1, TRUE, 1,
    # "X", 1, TRUE, 1,
    "B", 0, TRUE, 1,
    "C", 0, TRUE, 1,
    "D", 0, TRUE, 1
    # "Y", 0, FALSE, 1,
    # "Z", 1, TRUE, 1
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~hill,
    "A", "B", 1L, 1, 2,
    "B", "C", 1L, 1, 2,
    "C", "D", 1L, 1, 2,
    "D", "B", -1L, 100, 2
    # "X", "Y", 1L, 10, 2,
    # "Y", "Z", -1L, 10, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    # "s0", "s1", "+A,+X,+B,+C,+Z", TRUE, TRUE, 100,
    # "s1", "s2", "+D,-C,+Y", FALSE, FALSE, 100,
    # "s2", "s3", "+C,-A,-B", FALSE, FALSE, 100,
    # "s3", "s1", "+A,+B,-D", FALSE, FALSE, 100
    "s0", "s1", "+A,+B,+C", TRUE, TRUE, 100,
    "s1", "s2", "+D,-C", FALSE, FALSE, 100,
    "s2", "s3", "+C,-B", FALSE, FALSE, 100,
    "s3", "s1", "+B,-D", FALSE, FALSE, 100
  )
  
  backbone(module_info, module_network, expression_patterns)
}

make_config <- function(n_cells, n_genes) {
  list(
    num_cells    = n_cells,
    num_tfs      = 4L,
    num_targets  = n_genes,
    num_hks      = 0,
    census_interval = 2,
    tau_seconds     = 30 / 3600,
    weight_bw       = 0.1,
    distance_metric = "euclidean",
    min_tfs_per_module = 1L,
    max_in_degree      = 1L,
    damping            = 0.01,
    target_resampling  = Inf,
    num_simulations = 50
  )
}

run_simulation <- function(model_config) {
  model <- initialise_model(
    backbone = backbone_cycle_simple(),
    num_cells = model_config$num_cells,
    num_tfs = model_config$num_tfs,
    num_targets = model_config$num_targets,
    num_hks = model_config$num_hks,
    distance_metric = model_config$distance_metric,
    tf_network_params = tf_network_default(
      min_tfs_per_module = model_config$min_tfs_per_module,
      sample_num_regulators = function() 1,
      weighted_sampling = FALSE
    ),
    feature_network_params = feature_network_default(
      realnet = realnet_matrix,
      damping = model_config$damping,
      target_resampling = model_config$target_resampling,
      max_in_degree = model_config$max_in_degree
    ),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(
      census_interval = model_config$census_interval,
      ssa_algorithm = ssa_etl(tau = model_config$tau_seconds),
      experiment_params = bind_rows(
        simulation_type_wild_type(num_simulations = model_config$num_simulations),
        simulation_type_knockdown(num_simulations = 0)
      )
    ),
    experiment_params = experiment_snapshot(
      realcount = realcount_data,
      map_reference_cpm = TRUE,
      map_reference_ls = TRUE,
      weight_bw = model_config$weight_bw
    ),
    verbose = TRUE,
    download_cache_dir = getOption("dyngen_download_cache_dir"),
    num_cores = getOption("Ncpus") %||% 1L,
    id = NULL
  )
  
  out <- dyngen::generate_dataset(model, format = "dyno")
  model <- out$model
  # dyngen's dyno-format `dataset$expression` is log2(counts + 1) - already log-transformed,
  # meant for trajectory-inference tools, not for feeding into Cyclum/scPrisma (which each do
  # their own CPM/log-normalization internally and expect raw counts as input, same as the
  # real UMI datasets elsewhere in this repo). `dataset$counts` (= counts_mrna + counts_premrna)
  # is the raw, non-negative-integer simulated molecule count matrix - use that instead.
  counts <- t(as.matrix(out$dataset$counts))

  list(
    model = model,
    counts = counts
  )
}

model_config <- make_config(200, 200)
sim          <- run_simulation(model_config)

fname_base <- sprintf("test_%dgenes_%dcells", nrow(sim$counts), ncol(sim$counts))

subdir   <- file.path(proj_root, "synthetic", "data")
h5_file  <- here::here(subdir, paste0(fname_base, ".h5"))
sim_file <- here::here(subdir, paste0(fname_base, "_sim.rds"))
dir.create(dirname(h5_file), showWarnings = FALSE, recursive = TRUE)

mat2hdf(sim$counts, h5_file)
saveRDS(sim, sim_file)

cat("Saved to:", h5_file, "\n")
cat("Dimensions:", nrow(sim$counts), "genes ×", ncol(sim$counts), "cells\n")
cat("Row names (genes):", length(rownames(sim$counts)), "saved\n")
cat("Column names (cells):", length(colnames(sim$counts)), "saved\n")

plot_simulation_expression(sim$model, 1, what = "mol_mrna")
