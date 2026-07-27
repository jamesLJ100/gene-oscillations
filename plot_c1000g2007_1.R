library(dyngen)
library(tidyverse)

proj_root <- here::here()

fname    <- "c1000g5004_1"
sim_file <- file.path(proj_root, "data/dyngen_new", paste0(fname, "_sim.rds"))

if (!file.exists(sim_file)) {
  stop(sprintf("Sim file not found: %s", sim_file))
}

sim <- readRDS(sim_file)

p <- plot_feature_network(sim$model)
print(p)
