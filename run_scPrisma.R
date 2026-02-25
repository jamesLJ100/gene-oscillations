library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)

input_dir <- file.path(proj_root, "data/dyngen")

system2(
  command = py_config()$python,
  args    = c("-u", "python/run_scPrisma.py", input_dir)
)