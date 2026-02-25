library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("cyclum_env", required = TRUE)

input_dir <- file.path(proj_root, "data/GSM4286760")

system2(
  command = py_config()$python,
  args    = c("-u", "python/run_cyclum.py", input_dir)
)