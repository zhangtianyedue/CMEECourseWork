#!/bin/bash

# Print start message
echo "Starting MiniProject..."

# Define directories
PROJECT_DIR="$(dirname "$0")"
RESULTS_DIR="$PROJECT_DIR/results"
LOG_FILE="$RESULTS_DIR/run_log.txt"

# Ensure the results directory exists
mkdir -p "$RESULTS_DIR"

# Log start time
echo "Run started at: $(date)" > "$LOG_FILE"

# Step 1: Check for dependencies (Python & R)
echo "Checking dependencies..."
if ! command -v python3 &> /dev/null; then
    echo "Error: Python3 is not installed!" | tee -a "$LOG_FILE"
    exit 1
fi
if ! command -v Rscript &> /dev/null; then
    echo "Error: R is not installed!" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Python and R found! Proceeding..." | tee -a "$LOG_FILE"

# Step 2: Install required R packages
echo "Checking and installing required R packages..."
Rscript -e 'required_packages <- c("ggplot2", "dplyr", "minpack.lm", "gridExtra", "grid", "tibble", "rstan", "loo", "bayesplot");
           missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package")];
           if(length(missing_packages)) install.packages(missing_packages, repos="https://cloud.r-project.org/")' | tee -a "$LOG_FILE"

echo "R packages check completed." | tee -a "$LOG_FILE"

# Step 3: Install required Python packages
echo "Checking and installing required Python packages..."
PYTHON_PACKAGES=("matplotlib" "numpy" "pandas" "scipy" "seaborn" "stan")
for pkg in "${PYTHON_PACKAGES[@]}"; do
    if ! python3 -c "import $pkg" &> /dev/null; then
        echo "Installing Python package: $pkg" | tee -a "$LOG_FILE"
        python3 -m pip install --user $pkg
    fi
done
echo "Python packages check completed." | tee -a "$LOG_FILE"

# Step 4: Run data preprocessing (must be run first)
echo "Running data preprocessing..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/miniprojectdatacleaning.R" || { echo "Data cleaning failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 5: Run global model fitting
echo "Running global model fitting..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/globlefitting5model.R" || { echo "Global model fitting failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 6: Run time-segmented model fitting
echo "Running time-segmented model fitting..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/time_sigmentfitting5model.R" || { echo "Time-segmented model fitting failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 7: Run F-complexity test
echo "Running F-complexity test..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/Fcomplexicitytest.R" || { echo "F-complexity test failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 8: Run Bayesian analysis for Gompertz
echo "Running Bayesian analysis for Gompertz model..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/GOMBayesian.R" || { echo "Gompertz Bayesian analysis failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 9: Run Bayesian analysis for Logistic
echo "Running Bayesian analysis for Logistic model..." | tee -a "$LOG_FILE"
Rscript "$PROJECT_DIR/LogisticBayesian.R" || { echo "Logistic Bayesian analysis failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 10: Run plotting for time-segmented models
echo "Generating plots for time-segmented models..." | tee -a "$LOG_FILE"
python3 "$PROJECT_DIR/plot_timesigment.py" || { echo "Plotting time-segmented models failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 11: Run plotting for global models
echo "Generating plots for global models..." | tee -a "$LOG_FILE"
python3 "$PROJECT_DIR/plotgloble.py" || { echo "Plotting global models failed!" | tee -a "$LOG_FILE"; exit 1; }

# Step 12: Completion message
echo "MiniProject completed successfully!" | tee -a "$LOG_FILE"
echo "Results saved in $RESULTS_DIR"

echo "Run completed at: $(date)" >> "$LOG_FILE"
exit 0