library(testthat)
library(here)
library(tidyverse)
library(readr)

proj_root <- here::here()
setwd(proj_root)

source(here::here("dyngen_utils.R"))
test_that("propagate_module_assignments matches ground truth", {
  sim <- readRDS(here::here("tests", "data", "c1000g504_1_sim.rds"))
  
  feature_net <- sim$model$feature_network
  
  tf_modules <- feature_net %>%
    filter(!is.na(from_module)) %>%
    select(gene = from, module = from_module) %>%
    distinct()
  
  target_edges <- feature_net %>%
    filter(grepl("^Target", to)) %>%
    select(from, to)
  
  ground_truth <- read_csv(here::here("tests", "data", "ground_truth_modules.csv"), 
                           show_col_types = FALSE) %>%
    rename(gene = gene_name, module = module_assignment)
  
  result <- propagate_module_assignments(tf_modules, target_edges)
  
  for (i in seq_len(nrow(ground_truth))) {
    gene <- ground_truth$gene[i]
    expected_module <- ground_truth$module[i]
    
    expect_true(gene %in% names(result),
                label = sprintf("Gene '%s' should be present in result", gene))
    
    expect_equal(result[[gene]], expected_module,
                 label = sprintf("Gene '%s' should be in module '%s', got '%s'", 
                                 gene, expected_module, result[[gene]]))
  }
  
  missing <- ground_truth$gene[!(ground_truth$gene %in% names(result))]
  expect_equal(length(missing), 0L,
               label = sprintf("Unresolved genes: %s", paste(missing, collapse = ", ")))
})