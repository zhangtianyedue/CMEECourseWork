# 📌 Mini Project Code Documentation

## 🚀 Overview
This directory contains all the scripts necessary to run the **Mini Project** for microbial growth model selection and Bayesian inference. The workflow integrates **R and Python** to preprocess data, fit models, perform Bayesian analysis, and generate visualizations.

## ⚠️ Important Notice
**To ensure smooth execution, always run scripts from within the `code/` directory.**  
Since relative paths are used in the scripts, changing directories before running commands may lead to errors.

---

## 📜 Workflow Summary
The analysis follows a structured pipeline:

1. **Check Dependencies:** Verify that Python and R are installed, and required packages are available.
2. **Preprocess Data:** Clean and format the input dataset for model fitting.
3. **Global Model Fitting:** Fit various growth models (Logistic, Gompertz, Richards, etc.).
4. **Segmented Model Fitting:** Fit time-segmented models for better resolution.
5. **Complexity Testing:** Perform F-tests to compare model complexities.
6. **Bayesian Inference:** Apply Bayesian analysis to Logistic and Gompertz models.
7. **Generate Plots:** Visualize global and segmented model results.
8. **Store Results:** Outputs are saved in the `results/` directory.

---

## 📂 File Structure
code/ │── Run Miniproject.sh # Master script to automate the pipeline │── miniprojectdatacleaning.R # Preprocesses and cleans dataset │── globlefitting5model.R # Fits multiple global models │── time_sigmentfitting5model.R # Fits time-segmented models │── Fcomplexicitytest.R # Performs complexity analysis │── GOMBayesian.R # Bayesian analysis for Gompertz model │── LogisticBayesian.R # Bayesian analysis for Logistic model │── plot_timesigment.py # Plots results for segmented models │── plotgloble.py # Plots results for global models │── main.tex # LaTeX source file for documentation │── references.bib # Bibliography file │── supplementary.tex # Additional LaTeX materials


---

## 🛠️ How to Run
To execute the full pipeline, **ensure you are inside the `code/` directory** and then run:

```bash
bash Run\ Miniproject.sh

# Step 1: Data preprocessing
Rscript miniprojectdatacleaning.R

# Step 2: Global model fitting
Rscript globlefitting5model.R

# Step 3: Time-segmented model fitting
Rscript time_sigmentfitting5model.R

# Step 4: Complexity test
Rscript Fcomplexicitytest.R

# Step 5: Bayesian analysis
Rscript GOMBayesian.R
Rscript LogisticBayesian.R

# Step 6: Generate plots
python3 plot_timesigment.py
python3 plotgloble.py

📊 Output
All results, including fitted models, statistical outputs, and plots, are saved in the ../results/ directory.

Processed data: Cleaned_LogisticGrowthData.csv
Model fitting results: *.csv
Bayesian analysis results: *.rds
Plots: *.png



