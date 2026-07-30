library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
use_condaenv("scPrisma_env", required = TRUE)

# Define grid search parameters
regularisation_strength <- c(0, 0.1, 0.2, 0.3, 0.4)

iternum <- 100

# Input directory
input_dir <- file.path(proj_root, "synthetic/data/dyngen/gridsearch")

# Get Python executable
python_exe <- py_config()$python

# Total combinations
total_runs <- length(regularisation_strength)
cat(sprintf("Starting grid search with %d total combinations\n", total_runs))

# Initialize data frame to store run information
run_log <- data.frame(
  run_counter = integer(),
  reg_strength = numeric(),
  stringsAsFactors = FALSE
)

# Counter for progress tracking
run_counter <- 0

# Grid search loop
for (reg_str in regularisation_strength) {
  run_counter <- run_counter + 1
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Run %d/%d\n", run_counter, total_runs))
  cat(sprintf("Regularisation strength: %.1f\n", reg_str))
  cat(sprintf("========================================\n\n"))
  
  # Build command arguments
  args <- c(
    "-u", 
    "python/run_scPrisma.py",
    input_dir,
    "--regularisation_strength", as.character(reg_str),
    "--iternum", as.character(iternum),
    "--run_number", as.character(run_counter)
  )
  
  # Run the Python script with real-time output
  exit_code <- system2(
    command = python_exe,
    args = args,
    stdout = "",  
    stderr = ""   
  )
  
  # Check if the run was successful
  if (exit_code != 0) {
    cat(sprintf("WARNING: Run %d exited with code %d\n", run_counter, exit_code))
  }
  
  # Log this run
  run_log <- rbind(run_log, data.frame(
    run_counter = run_counter,
    reg_strength = reg_str,
    stringsAsFactors = FALSE
  ))

}

cat("\n========================================\n")
cat("Grid search completed!\n")
cat(sprintf("Total runs: %d\n", run_counter))
cat("========================================\n")

# Save run log to CSV
output_csv <- file.path(proj_root, "synthetic/data/dyngen/gridsearch/scPrisma/grid_search_log.csv")
write.csv(run_log, output_csv, row.names = FALSE)
cat(sprintf("\nGrid search log saved to: %s\n", output_csv))

# Print summary
cat("\nParameter combinations tested:\n")
print(run_log)