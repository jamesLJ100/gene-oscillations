library(reticulate)
library(here)

proj_root <- here::here()
setwd(proj_root)
use_condaenv("cyclum_env", required = TRUE)

#TODO nonlinear_reg & check for others.

# Define grid search parameters
encoder_widths <- list(
  c(30, 20),
  c(40, 30),
  c(50, 50),
  c(50, 40, 30),
  c(60, 40, 20)
)

epochs_list <- c(300, 500, 700)

learning_rates <- c(1e-4, 2e-4, 5e-4, 1e-3)

# Input directory
input_dir <- file.path(proj_root, "synthetic/data/dyngen/tuning")

# Get Python executable
python_exe <- py_config()$python

# Total combinations
total_runs <- length(encoder_widths) * length(epochs_list) * length(learning_rates)
cat(sprintf("Starting grid search with %d total combinations\n", total_runs))
cat(sprintf("Encoder widths: %d variations\n", length(encoder_widths)))
cat(sprintf("Epochs: %d variations\n", length(epochs_list)))
cat(sprintf("Learning rates: %d variations\n\n", length(learning_rates)))

# Initialize data frame to store run information
run_log <- data.frame(
  run_counter = integer(),
  encoder_width = character(),
  epochs = integer(),
  learning_rate = numeric(),
  stringsAsFactors = FALSE
)

# Counter for progress tracking
run_counter <- 0

# Grid search loop
for (encoder_width in encoder_widths) {
  for (epochs in epochs_list) {
    for (learning_rate in learning_rates) {
      run_counter <- run_counter + 1
      
      # Format encoder width for display and storage
      ew_str <- paste(encoder_width, collapse=", ")
      
      cat(sprintf("\n========================================\n"))
      cat(sprintf("Run %d/%d\n", run_counter, total_runs))
      cat(sprintf("Encoder width: [%s]\n", ew_str))
      cat(sprintf("Epochs: %d\n", epochs))
      cat(sprintf("Learning rate: %.0e\n", learning_rate))
      cat(sprintf("========================================\n\n"))
      
      # Build command arguments
      args <- c(
        "-u", 
        "python/run_cyclum.py",
        input_dir,
        "--encoder_width", encoder_width,
        "--epochs", as.character(epochs),
        "--learning_rate", sprintf("%.0e", learning_rate),
        "--run_number", run_counter
      )
      
      # Run the Python script with real-time output
      exit_code <- system2(
        command = python_exe,
        args = args,
        stdout = "",  
        stderr = ""   
      )
      
      # Log this run
      run_log <- rbind(run_log, data.frame(
        run_counter = run_counter,
        encoder_width = ew_str,
        epochs = epochs,
        learning_rate = learning_rate,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n========================================\n")
cat("Grid search completed!\n")
cat(sprintf("Total runs: %d\n", run_counter))
cat("========================================\n")

# Save run log to CSV
output_csv <- file.path(proj_root, "synthetic/data/dyngen/tuning/cyclum/grid_search_log.csv")
write.csv(run_log, output_csv, row.names = FALSE)
cat(sprintf("\nGrid search log saved to: %s\n", output_csv))