# CMEE 2024 HPC exercises - Cluster run script
# This script runs neutral model simulations on the HPC cluster

# -------------------------------------------
# 1. Clear the environment and close the graphics device
rm(list = ls())  # Clear all variables
graphics.off()   # Close all graphics windows

# -------------------------------------------
# 2. Source the main code to ensure all functions are available
source("tz124_HPC_2024_main.R")

# -------------------------------------------
# 3. Get the task number iter
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))  # 从Get iter from HPC
if (is.na(iter) | iter == "") {  
  #stop("Error: PBS_ARRAY_INDEX is not set or is invalid.")  # Encountered an error and terminated the execution
}

# ensure iter is integer
iter <- as.integer(iter)
print(paste("✅ Running simulation with iter =", iter))  #

# -------------------------------------------
# 4. set seed 
set.seed(iter)

# -------------------------------------------
# 5. choose community size
if (iter >= 1 & iter <= 25) {
  size <- 500
} else if (iter >= 26 & iter <= 50) {
  size <- 1000
} else if (iter >= 51 & iter <= 75) {
  size <- 2500
} else {
  size <- 5000
}

# -------------------------------------------
# 6. Set the speciation rate 
speciation_rate <- 0.1  

# -------------------------------------------
# 7. Set the output file name to ensure that each iter output file is unique
output_file_name <- paste0("simulation_output_", iter, ".rda")

# -------------------------------------------
# 8. run simulation
neutral_cluster_run(
  speciation_rate = speciation_rate,
  size = size,
  wall_time = 690,  # run 690 min（11.5 h）
  interval_rich = 1,
  interval_oct = size / 10,
  burn_in_generations = 8 * size,
  output_file_name = output_file_name
)

# -------------------------------------------

if (file.exists(output_file_name)) {
  print(paste("✅ Simulation completed successfully. Output saved as", output_file_name))
} else {
  print("❌ Error: .rda file was NOT created!")
}

# 10. end
print("🎯 HPC job finished.")
