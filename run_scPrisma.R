library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)

# input_dir <- file.path(proj_root, "data/dyngen")
# 
# system2(
#   command = py_config()$python,
#   args    = c("-u", "python/run_scPrisma.py", input_dir)
# )


#Run on HCT116 data

hct116_dir <- file.path(proj_root, "data/GSM4286760")

python_exe <- py_config()$python

reg_str <- 0.1
iter <- 100

args <- c(
  "-u", 
  "python/run_scPrisma.py",
  hct116_dir,
  "--regularisation_strength", as.character(reg_str),
  "--iternum", as.character(iter)
)

# Run the Python script with real-time output
exit_code <- system2(
  command = python_exe,
  args = args,
  stdout = "",  
  stderr = ""   
)


#Run on mme95
# mme95_dir <- file.path(proj_root, "data/mme95/tpm")
# 
# python_exe <- py_config()$python
# 
# reg_str <- 0.1
# iter <- 100
# run_counter <- 0
# 
# args <- c(
#   "-u", 
#   "python/run_scPrisma.py",
#   mme95_dir,
#   "--regularisation_strength", as.character(reg_str),
#   "--iternum", as.character(iter),
#   "--run_number", as.character(run_counter)
# )
# 
# # Run the Python script with real-time output
# exit_code <- system2(
#   command = python_exe,
#   args = args,
#   stdout = "",  
#   stderr = ""   
# )

