library(dyngen)
library(tidyverse)
library(here)
library(hdf5r)

## define project root
proj_root <- here::here()
setwd(proj_root)

source(file.path(proj_root, "hdfrw.R"))

##Todo add explicit downloads rather than mystery files 
#dyngen::list_realnetworks()
realnet_matrix <- readRDS(here::here("data", "realnet.rds"))
realcount_data <- readRDS(here::here("data", "realcount.rds"))

# definition
backbone_cycle_simple <- function() {
  module_info <- tribble(
    ~module_id, ~basal, ~burn, ~independence,
    "A", 1, TRUE, 1,
    "B", 0, FALSE, 1,
    "C", 0, FALSE, 1,
    "D", 0, FALSE, 1
  )
  
  module_network <- tribble(
    ~from, ~to, ~effect, ~strength, ~hill,
    "A", "B", 1L, 1, 2,
    "B", "C", 1L, 1, 2,
    "C", "D", 1L, 1, 2,
    "D", "B", -1L, 100, 2
  )
  
  expression_patterns <- tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "sBurn", "s1", "+A", TRUE, TRUE, 100,
    "s1", "s2", "+B,-A", FALSE, FALSE, 100,
    "s2", "s3", "+C", FALSE, FALSE, 100,
    "s3", "s1", "+D", FALSE, FALSE, 500
  )
  
  backbone(module_info, module_network, expression_patterns)
}

make_config <- function(n_cells, n_genes) {
  list(
    num_cells    = n_cells,
    num_tfs      = 4L,
    num_targets  = n_genes/2,
    num_hks      = n_genes/2,
    census_interval = 2,
    tau_seconds     = 30 / 3600,
    weight_bw       = 0.1,
    distance_metric = "euclidean",
    min_tfs_per_module = 1L,
    max_in_degree      = 1L,
    damping            = 0.01,
    target_resampling  = Inf
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
      ssa_algorithm = ssa_etl(tau = model_config$tau_seconds)
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
  expression <- t(as.matrix(out$dataset$expression))

  list(
    model = model,
    expression = expression
  )
}
