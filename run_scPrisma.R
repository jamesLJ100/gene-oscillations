rm(list=ls())
library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)

use_condaenv("scPrisma_env", required = TRUE)

dyngen_dir <- file.path(proj_root, "data/dyngen_new")


python_exe <- py_config()$python

reg_str <- 0.1
iter <- 100

args <- c(
  "-u",
  "python/run_scPrisma.py",
  dyngen_dir,
  "--regularisation_strength", as.character(reg_str),
  "--iternum", as.character(iter)
)

exit_code <- system2(
  command = python_exe,
  args = args,
  stdout = "",  
  stderr = ""   
)


#Run on HCT116 data

# hct116_dir <- file.path(proj_root, "data/GSM4286760")
# 
# python_exe <- py_config()$python
# 
# reg_str <- 0.1
# iter <- 100
# 
# args <- c(
#   "-u", 
#   "python/run_scPrisma.py",
#   hct116_dir,
#   "--regularisation_strength", as.character(reg_str),
#   "--iternum", as.character(iter)
# )
# 
# # Run the Python script with real-time output
# exit_code <- system2(
#   command = python_exe,
#   args = args,
#   stdout = "",  
#   stderr = ""   
# )


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

# #Run on NFKB
# nfkb_dir <- file.path(proj_root, "data/GSE162992/processed")
# 
# python_exe <- py_config()$python
# 
# reg_str <- 0.2
# iter <- 100
# 
# args <- c(
#   "-u",
#   "python/run_scPrisma.py",
#   nfkb_dir,
#   "--regularisation_strength", as.character(reg_str),
#   "--iternum", as.character(iter)
# )
# 
# # Run the Python script with real-time output
# exit_code <- system2(
#   command = python_exe,
#   args = args,
#   stdout = "",  
#   stderr = ""   
# )

