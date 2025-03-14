# HPC Stochastic Demographic and Neutral Model Simulations

This repository contains R scripts for running **stochastic demographic** and **neutral model simulations** on an **HPC cluster** as part of the CMEE 2024 HPC exercises. The scripts are designed to efficiently simulate population dynamics and analyze biodiversity patterns using parallel computing.

---

## **📌 Project Overview**
This project involves:
- **Stochastic Demographic Model**: Simulating population growth and extinction under different initial conditions.
- **Neutral Model Simulations**: Simulating biodiversity dynamics using a neutral theory framework.
- **HPC Parallel Execution**: Running simulations on multiple nodes using job submission scripts.

---

## **📂 Files and Directory Structure**
### **1️⃣ Core R Scripts**
| File | Description |
|------|------------|
| `tz124_HPC_2024_main.R` | Contains core functions for stochastic simulations and neutral model analysis. |
| `tz124_HPC_2024_demographic_cluster.R` | Runs demographic model simulations on the HPC cluster. |
| `tz124_HPC_2024_neutral_cluster.R` | Runs neutral model simulations on the HPC cluster. |

### **2️⃣ HPC Job Submission Scripts**
| File | Description |
|------|------------|
| `submit_hpc_job_neutral.sh` | Submits neutral model simulation jobs to the HPC cluster. |
| `run_job.sh` | Submits demographic model simulation jobs to the HPC cluster. |

### **3️⃣ Output Files**
| File | Description |
|------|------------|
| `simulation_results_*.rda` | Stores simulation output for demographic models. |
| `simulation_output_*.rda` | Stores simulation output for neutral model runs. |
| `population_size_data.csv` | Stores formatted simulation data for visualization. |
| `Challenge_A.png` | Visualization of stochastic population trends. |

---

## **🚀 Running Simulations on the HPC**
### **1️⃣ Submitting a Job**
To run simulations, submit jobs using the provided shell scripts:

#### **Run the Demographic Model**
```bash
qsub run_job.sh

#### **Run the neutral Model**
```bash
qsub submit_hpc_job_neutral.sh

