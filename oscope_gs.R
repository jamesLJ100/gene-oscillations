library(here)

proj_root <- here::here()
setwd(proj_root)

# Source the run_oscope.R script
source(file.path(proj_root, "run_oscope.R"))

# Define grid search parameters
#KM_quan_values <- c(0.95, 0.75, 0.5, 0.25, 0.05)
KM_quan_values <- c(0.05)
#FlagCluster_qt_values <- c(0.90, 0.95, 0.85)
FlagCluster_qt_values <- c(0.90)
#FlagCluster_thre_values <- c(pi/5, pi/4, pi/3)
FlagCluster_thre_values <- c(pi/4)

# Fixed parameters (not part of grid search)
fixed_KM_maxK <- 10
fixed_CalcMV_MeanCutLow <- 0.1

# Input directory
input_dir <- file.path(proj_root, "data/dyngen/tuning")

# Total combinations
total_runs <- length(KM_quan_values) * length(FlagCluster_qt_values) * length(FlagCluster_thre_values)
cat(sprintf("Starting Oscope grid search with %d total combinations\n", total_runs))

# Initialize data frame to store run information
run_log <- data.frame(
  run_counter = integer(),
  KM_quan = numeric(),
  FlagCluster_qt = numeric(),
  FlagCluster_thre = numeric(),
  KM_maxK = integer(),
  CalcMV_MeanCutLow = numeric(),
  stringsAsFactors = FALSE
)

# Counter for progress tracking
run_counter <- 0

# Grid search loop
for (km_quan in KM_quan_values) {
  for (flag_qt in FlagCluster_qt_values) {
    for (flag_thre in FlagCluster_thre_values) {
      run_counter <- run_counter + 1
      
      cat(sprintf("\n========================================\n"))
      cat(sprintf("Run %d/%d\n", run_counter, total_runs))
      cat(sprintf("KM_quan: %.2f\n", km_quan))
      cat(sprintf("FlagCluster_qt: %.2f\n", flag_qt))
      cat(sprintf("FlagCluster_thre: %.4f (%.2f rad)\n", flag_thre, flag_thre))
      cat(sprintf("KM_maxK: %d (fixed)\n", fixed_KM_maxK))
      cat(sprintf("CalcMV_MeanCutLow: %.2f (fixed)\n", fixed_CalcMV_MeanCutLow))
      cat(sprintf("========================================\n\n"))
      
      # Build oscope_config for this run
      oscope_config <- list(
        KM_maxK           = fixed_KM_maxK,
        KM_quan           = km_quan,
        CalcMV_MeanCutLow = fixed_CalcMV_MeanCutLow,
        FlagCluster_qt    = flag_qt,
        FlagCluster_thre  = flag_thre,
        Normalise=TRUE
      )
      
      # Run Oscope with this configuration
      tryCatch({
        run_oscope(input_dir, oscope_config, run_counter)
        cat(sprintf("Run %d completed successfully\n", run_counter))
      }, error = function(e) {
        cat(sprintf("ERROR in Run %d: %s\n", run_counter, e$message))
      })
      
      # Log this run
      run_log <- rbind(run_log, data.frame(
        run_counter = run_counter,
        KM_quan = km_quan,
        FlagCluster_qt = flag_qt,
        FlagCluster_thre = flag_thre,
        KM_maxK = fixed_KM_maxK,
        CalcMV_MeanCutLow = fixed_CalcMV_MeanCutLow,
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
output_csv <- file.path(proj_root, "data/dyngen/tuning/oscope/grid_search_log.csv")
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(run_log, output_csv, row.names = FALSE)
cat(sprintf("\nGrid search log saved to: %s\n", output_csv))

# Print summary
cat("\nParameter combinations tested:\n")
print(run_log)