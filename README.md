
# Final Version of Assignment

This archive contains the **final version of the coursework**, integrating feedback from the first three weeks. All issues and suggestions from the feedback have been carefully addressed. This version represents the culmination of all efforts, showcasing a polished and complete submission.

## Folder Structure

Below is the detailed file structure of the archive, organized by week and category:

### Week 1
- **Folder Structure**:
  - **code/**: Scripts for basic analysis and data manipulation.
  - **data/**: Datasets including `testcsv.csv` and `sequences.csv`.
  - **results/**: Placeholder for output files.
  - **sandbox/**: Reserved for experimental or temporary files.

### Week 2
- **Folder Structure**:
  - **code/**: Scripts for more advanced operations, such as data cleaning and statistical modeling.
  - **data/**: Includes datasets like `bodymass.csv`.
  - **results/**: Placeholder for generated results.
  - **sandbox/**: Reserved for drafts and exploratory tasks.

### Week 3
- **Folder Structure**:
  - **code/**: Contains R and Python scripts, including:
    - `align_seqs_fasta.py`: A Python script for sequence alignment.
    - `oaks_debugme.py`: Debugging script for tree data.
    - Other R scripts for data visualization and control flow tasks.
  - **data/**: Contains datasets like `PoundHillData.csv`, `.fasta` files, and metadata files.
  - **results/**: Placeholder for analysis outputs.
  - **sandbox/**: Reserved for exploratory work.

### Week 4
- **Folder Structure**:
  - **code/**: R scripts for trend analysis and LaTeX files for documentation, including:
    - `Florida.R`: Analysis script for Florida temperature trends.
    - `PP_Regress.R`: Regression analysis script.
  - **data/**: Includes datasets such as `KeyWestAnnualMeanTemperature.RData` and supporting images (`1.png`, `2.png`).
  - **results/**: Placeholder for final outputs.
  - **sandbox/**: Reserved for drafts and experimental tasks.

### Other Folders
- **results/**: Consolidated results and placeholders for outputs across all weeks.
- **data/**: Shared datasets like `TestOaksData.csv` and `bodymass.csv`.
- **sandbox/**: Reserved for experimental or temporary files.

## Complete File List

```plaintext
归档/
├── week1/
│   ├── code/
│   │   ├── script1.R
│   │   ├── script2.R
│   ├── data/
│   │   ├── testcsv.csv
│   │   ├── sequences.csv
│   ├── results/
│   │   ├── .gitkeep
│   ├── sandbox/
│       ├── .gitkeep
├── week2/
│   ├── code/
│   │   ├── script1.R
│   │   ├── script2.R
│   ├── data/
│   │   ├── bodymass.csv
│   ├── results/
│   │   ├── .gitkeep
│   ├── sandbox/
│       ├── .gitkeep
├── week3/
│   ├── code/
│   │   ├── align_seqs_fasta.py
│   │   ├── oaks_debugme.py
│   │   ├── Girko.R
│   ├── data/
│   │   ├── PoundHillData.csv
│   │   ├── PoundHillMetaData.csv
│   │   ├── 407228412.fasta
│   ├── results/
│   │   ├── .gitkeep
│   ├── sandbox/
│       ├── .gitkeep
├── week4/
│   ├── code/
│   │   ├── Florida.R
│   │   ├── PP_Regress.R
│   ├── data/
│   │   ├── KeyWestAnnualMeanTemperature.RData
│   │   ├── 1.png
│   │   ├── 2.png
│   ├── results/
│   │   ├── .gitkeep
│   ├── sandbox/
│       ├── .gitkeep
├── results/
│   ├── README.md
├── data/
│   ├── bodymass.csv
│   ├── sequences.csv
│   ├── TestOaksData.csv
├── sandbox/
│   ├── .gitkeep
```

## Key Features and Improvements
1. **Feedback Integration**: Adjusted code, data handling, and documentation based on feedback from Weeks 1–3.
2. **Comprehensive Documentation**: Each week includes a `README.md` file with specific instructions.
3. **Enhanced Modularity**: Scripts and data are well-organized to ensure ease of use.

## How to Use
1. Navigate to the relevant `weekX/` folder for specific tasks.
2. Review the `README.md` file in each folder for guidance on running scripts.
3. Install required dependencies for R or Python using the provided environment files or package lists.

## 🚀 Latest Additions (March 14, 2025)
On **March 14, 2025, at 02:31 AM**, I added two major components:
- **Mini Project**: A structured analysis pipeline for microbial growth model selection and Bayesian inference.
- **HPC Work**: A high-performance computing (HPC) workflow for large-scale simulations and data analysis.

These components integrate **R, Python, and HPC resources**, enabling **automated model selection, Bayesian inference, and large-scale simulations**.

---

## 📂 Project Structure

### 🔬 Mini Project (`miniproject/`)
miniproject/ │── code/ │ ├── Run Miniproject.sh │ ├── miniprojectdatacleaning.R │ ├── globlefitting5model.R │ ├── time_sigmentfitting5model.R │ ├── Fcomplexicitytest.R │ ├── GOMBayesian.R │ ├── LogisticBayesian.R │ ├── plot_timesigment.py │ ├── plotgloble.py │ ├── main.tex │ ├── references.bib │ └── supplementary.tex │ │── data/ │ ├── LogisticGrowthData.csv │ │── results/ │ ├── *.csv │ ├── *.png │ └── *.rds │ │── sandbox/ │ └── .gitkeep


### 🖥️ HPC Work (`hpcwork/`)
hpcwork/ │── simulation_results/ # Contains multiple .rda simulation files │── demographic_cluster/ # Stores HPC output logs (.o & .e files) │── scripts/ # R scripts for running and analyzing simulations │── plots/ # Generated result plots │── processed_data/ # Processed .rda results for further analysis


---

## 🔧 How to Use
1. **Mini Project**
   - Run all analysis scripts in the `code/` folder.
   - Outputs are stored in the `results/` folder.

2. **HPC Work**
   - Includes HPC job submission scripts and result processing.
   - Large-scale simulation results are stored in `simulation_results/` and `processed_data/`.

---

### ✨ Notes
- The **Mini Project** automates model selection using **AIC, BIC, WAIC, and Bayes Factor (BF)**.
- The **HPC Work** enables **parallelized simulations and demographic clustering analysis**.
- All key results and figures are saved in the `results/` and `plots/` directories.
