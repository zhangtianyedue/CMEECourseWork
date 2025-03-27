# 📌 Mini Project Code Documentation

## 🚀 Overview

This directory contains all the scripts necessary to run the **Mini Project** for **microbial growth model selection and Bayesian inference**. The workflow integrates **R and Python** to preprocess data, fit models, perform Bayesian analysis (including MCMC inference using Stan), and generate visualizations.

### **⚠️ Important Notice**
- Always run scripts from within the `code/` directory to ensure smooth execution.
- **Changing directories before running commands may lead to errors** since relative paths are used.
- The provided `Run Miniproject.sh` script **has been extensively tested on a local machine** and runs without any issues when executed properly.

---

## 📜 **Workflow Summary**
This project follows a structured pipeline:

1. **Check Dependencies**: Ensures Python and R (and required packages) are installed.
2. **Preprocess Data**: Cleans and formats the dataset for model fitting.
3. **Global Model Fitting**: Fits multiple microbial growth models (Logistic, Gompertz, Richards, etc.).
4. **Segmented Model Fitting**: Fits time-segmented models for better resolution.
5. **Complexity Testing**: Performs F-tests to compare model complexities.
6. **Bayesian Inference**: Applies MCMC Bayesian analysis to **Logistic and Gompertz models** using `Stan`.
7. **Generate Plots**: Produces visualizations for global and segmented model results.
8. **Store Results**: Saves output in the `results/` directory.

---

## 📂 **File Structure**
MiniProject/ │── code/ │ │── Run Miniproject.sh # Master script to automate the pipeline │ │── miniprojectdatacleaning.R # Preprocesses and cleans dataset │ │── globlefitting5model.R # Fits multiple global models │ │── time_sigmentfitting5model.R # Fits time-segmented models │ │── Fcomplexicitytest.R # Performs complexity analysis │ │── GOMBayesian.R # Bayesian analysis for Gompertz model │ │── LogisticBayesian.R # Bayesian analysis for Logistic model │ │── plot_timesigment.py # Plots results for segmented models │ │── plotgloble.py # Plots results for global models │ │── bayesian_logistic.stan # Stan model for Bayesian Logistic inference │ │── bayesian_gompertz.stan # Stan model for Bayesian Gompertz inference │ │── main.tex # LaTeX source file for documentation │ │── references.bib # Bibliography file │ │── supplementary.tex # Additional LaTeX materials │ │── results/ # Processed output files │ │── README.md # Project documentation


---

## 🛠️ **How to Run**
To execute the full pipeline, navigate to the `code/` directory and run:
```bash
bash Run\ Miniproject.sh

Step 1: Data Preprocessing
Rscript miniprojectdatacleaning.R

Step 2: Global Model Fitting
Rscript globlefitting5model.R

Step 3: Time-Segmented Model Fitting
Rscript time_sigmentfitting5model.R

Step 4: Complexity Test
Rscript Fcomplexicitytest.R

Step 5: Bayesian Analysis (MCMC)
Rscript GOMBayesian.R
Rscript LogisticBayesian.R

Step 6: Generate Plots
python3 plot_timesigment.py
python3 plotgloble.py


📊 Output
All results, including fitted models, statistical outputs, and plots, are saved in the results/ directory.

File Type	Description
Processed Data	Cleaned_LogisticGrowthData.csv
Model Fitting Results	*.csv
Bayesian Analysis Results	*.rds
MCMC Stan Files	bayesian_logistic.stan, bayesian_gompertz.stan
Plots	*.png
🚀 Additional Notes
The Bayesian analysis relies on Stan models (bayesian_logistic.stan, bayesian_gompertz.stan), which must be present in the directory for successful execution.
The Run Miniproject.sh script has been thoroughly tested on a local machine and runs without issues.
If encountering package-related errors, refer to the dependency check inside Run Miniproject.sh which automatically installs missing Python and R packages.
👨‍💻 Author & Contact
Name: Tianye Zhang
Email: tz124@imperial.ac.uk
For any issues or troubleshooting, feel free to reach out.