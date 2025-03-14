#!/bin/bash

# Enable error handling
set -e

# Print the start message
echo "Starting local execution of all R scripts..."

# Run the main simulation script
echo "Running main simulation script: tz124_HPC_2024_main.R"
Rscript tz124_HPC_2024_main.R

# Run the demographic cluster simulation
echo "Running demographic cluster simulation: tz124_HPC_2024_demographic_cluster.R"
Rscript tz124_HPC_2024_demographic_cluster.R

# Run the neutral cluster simulation
echo "Running neutral cluster simulation: tz124_HPC_2024_neutral_cluster.R"
Rscript tz124_HPC_2024_neutral_cluster.R

# Print completion message
echo "✅ All simulations have completed successfully!"
