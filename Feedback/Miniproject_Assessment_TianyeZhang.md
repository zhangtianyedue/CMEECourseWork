# Miniproject Feedback and Assessment

## Report

**"Guidelines" below refers to the MQB report [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) [here](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html) (which were provided to the students in advance).**

**Title:** “Comparing Growth Models for Microbial Data: Gompertz and Logistic Outperform Richards and Linear Models”

- **Introduction (15%)**  
  - **Score:** 12/15  
  - Introduces Richards vs. logistic/Gompertz. Could emphasize the direct research aim or gap more explicitly.

- **Methods (15%)**  
  - **Score:** 12/15  
  - Data filtering approach is clear, mentions possible Bayesian or segmented approaches but snippet is brief on parameter details. The [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) recommends more thorough replicability details.

- **Results (20%)**  
  - **Score:** 15/20  
  - States Richards can overfit compared to Gompertz/Logistic. More numeric breakdown (e.g., how many subsets each model “won”) would follow [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) better.

- **Tables/Figures (10%)**  
  - **Score:** 6/10  
  - Mentions formula tables, but not deeply integrated in text. The [MQB Miniproject report guidelines](https://mulquabio.github.io/MQB/notebooks/Appendix-MiniProj.html#the-report) encourage explicit references to each table or figure.

- **Discussion (20%)**  
  - **Score:** 15/20  
  - Explains how the flexible Richards can lead to overfitting, while Gompertz/Logistic remain robust. Future directions or real-world usage can be expanded.

- **Style/Structure (20%)**  
  - **Score:** 15/20  
  - Good flow overall, some minor style oversights. Could cross-reference results in the discussion more thoroughly.

**Summary:** Very good presentation on how different models compare, highlighting Richards's potential to overfit. Would benefit from more explicit numeric detail in the results.

**Report Score:** 75

---

## Computing

### Project Structure & Workflow

**Strengths**

*The `Run Miniproject.sh` script sequences data cleaning, global and segmented fitting, complexity tests, Bayesian inference, and plotting across R and Python scripts, giving a clear end-to-end pipeline.

* Responsibilities are well-separated across modular scripts (but maybe too many!)

* The master shell script checks for Python and R availability, installs missing packages, and logs each step, good for reliability

**Suggestions**

1. **Shell script enhancements:**

   * Use a portable shebang (`#!/usr/bin/env bash`) and strict mode (`set -euo pipefail`) to exit on errors/undefined vars .
   * Rename `Run Miniproject.sh` to avoid spaces (e.g. `run_miniproject.sh`) and consistently quote file paths.
   * Parameterize data/result directories (flags or environment variables) for flexible reuse.
   * Consolidate all R steps into the driver rather than manual lines in the README.

2. **Reproducible environments:**

   *  Adopt **renv**; commit `renv.lock`, and call `renv::restore()` at script start or in the driver, avoiding ad-hoc `install.packages()` .
   * Provide a Python `requirements.txt` and use `pip install -r requirements.txt` rather than per-package installs.
   * Move `.stan` files to `code/stan/` and reference consistently in R scripts.

---

### README File

**Strengths**

* Clearly outlines each pipeline step, file structure, and dependencies .
* Includes author contact and notes on dataset origin.

**Suggestions**

1. Add explicit commands for venv (Python) and `renv` (R):

   ```bash
   # Python
   python3 -m venv venv && source venv/bin/activate
   pip install -r requirements.txt
   # R
   Rscript -e "renv::restore()"
   ```
2. Embed a tree diagram and provide one-line descriptions of each script's inputs/outputs for quick orientation.
3. Include a `LICENSE` file and attribution for `LogisticGrowthData.csv`.
4. Ensure code fences match actual script names (escape spaces) and remove duplicated or outdated sections.

---

## Script-Level Feedback

#### `miniprojectdatacleaning.R`

* Wrap steps (NA removal, type conversion, log transform, ID creation) into a `prepare_data()` function to improving testability.
* After each filter (e.g. `Time>=1`, `PopBio>=0`), print retained row and ID counts to document data loss.
* Use `dplyr::group_indices()` to generate numeric IDs or ensure string concatenation avoids ambiguous separators.

####  `globlefitting5model.R`

* Abstract common operations—initial guess generation, AICc/RSS computation, model evaluation—into dedicated functions.
* Instead of `for` + `rbind()`, use `purrr::map_dfr()` over `split(data, Unique_ID)` for efficiency and readability .
* Supply `lower`/`upper` in `nlsLM()` for all non-linear fits, and consider `nls.multstart::nls_multstart()` to streamline multi-start sampling and capture convergence diagnostics.
* The sliding-window max-slope heuristic is clever; extract it into a reusable function with configurable `window_size`.

####  `time_sigmentfitting5model.R`

* Share model definitions and helper functions with the global fitting script (e.g. via a common R file sourced at top).
* Unify `compute_rss()` and `compute_aicc()` definitions across scripts to avoid drift in formula implementations .

#### `Fcomplexicitytest.R`

* Replace manual `group_by()` + `summarise()` + `filter()` with a single pipeline that computes RSS for each model and then F-statistics, using `broom::glance()` for model summaries .
* Verify that `df1` and `df2` in `compute_f_stat()` reflect the correct difference in parameter counts and sample sizes per ID.
* Write both `stat_results` and `summary_stats` (including p-value distributions) to CSV in `results/` for downstream reporting.

### Bayesian Inference

* Include `Gompertz_Bayes.stan` and `Logistic_Bayes.stan` in a `code/stan/` directory, and reference them explicitly in R scripts via `here::here('code','stan','Gompertz_Bayes.stan')` and similarly for the logistic model.
* Move `stan()` calls into functions `run_bayesian_model(stan_file, data, init, control)` to avoid duplicated code.
* Pin Stan version and consider using **cmdstanr** for faster compilation, caching compiled models, and robust diagnostics.
* Use **tidybayes** to extract and visualize posterior intervals; combine plots into a multipanel figure rather than separate PNGs.
* Use `loo_compare()` for WAIC/LOO rather than manual Bayes factor approximations for consistent inference.
* Vectorize posterior predictive generation using matrix operations rather than `for` loops over samples to improve speed and clarity.

---

## Summary

Great job going above and beyond the call of duty! 

Hopefully then above feedback is helpful!

### **Score: 73**

---

## Overall Score: (75 + 73)/2 = 74