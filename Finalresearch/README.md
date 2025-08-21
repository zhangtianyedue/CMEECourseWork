# Spatial Structure Effects on Shannon Diversity Index Research

## Project Overview

This project investigates the effects of different spatial structures (1D, 2D, 3D) on Shannon diversity indices in ecological communities. Through numerical simulations, it analyzes how spatial diffusion, distance decay, and other factors influence species abundance distributions and community stability.

## Important Notice

**⚠️ CRITICAL: All code files must be in the same folder as `param.py` when running!**

## File Structure

```
Finalresearch/
├── param.py                           # Parameter module (required)
├── compare_space_structure_shannon.py # Baseline ecological simulation and Shannon diversity indices
├── frequency_boxplot_analysis.py      # Frequency distribution and boxplot analysis
├── analyze_distribution.py            # Distribution feature analysis
├── dimensionality_comparison.py       # Spatial dimension comparison analysis
└── distance_decay.py                  # Distance decay effect analysis
```

## Execution Order

**Follow this exact sequence when running the code:**

1. **compare_space_structure_shannon.py** - Generate baseline data and Shannon diversity indices first
2. **frequency_boxplot_analysis.py** - Frequency distribution and boxplot analysis (depends on step 1)
3. **analyze_distribution.py** - Distribution feature analysis (depends on step 1)
4. **dimensionality_comparison.py** - Spatial dimension comparison
5. **distance_decay.py** - Distance decay effect analysis

## File Function Details

### 1. param.py
**Function:** Core parameter module providing fundamental functions for ecological models
- `modular_uptake()`: Modular resource uptake function
- `modular_leakage()`: Modular leakage function  
- `generate_l_tensor()`: Generate leakage tensor

**Output Files:** None
**Dependencies:** None

### 2. compare_space_structure_shannon.py
**Function:** Generate baseline ecological simulation data and Shannon diversity indices
- Run ecological community simulations with different spatial structures
- Calculate alpha, beta, and gamma Shannon diversity indices
- Provide fundamental data for subsequent analyses

**Output Files:** None (provides data for other scripts)
**Dependencies:** `param.py`

### 3. frequency_boxplot_analysis.py
**Function:** Create combined analysis of frequency distributions and boxplots for Shannon diversity indices
- Generate frequency distribution histograms + kernel density estimation
- Create boxplots + scatter plots
- Provide statistical summaries

**Output Files:**
- `combined_frequency_boxplot.png` - Combined analysis plot

**Dependencies:** `param.py` + `compare_space_structure_shannon.py`

### 4. analyze_distribution.py
**Function:** In-depth analysis of Shannon diversity index distribution characteristics
- Normality tests (Shapiro-Wilk, D'Agostino K², Anderson-Darling)
- Distribution visualization (histograms, kernel density estimation, theoretical normal distribution)
- QQ plot analysis
- Statistical feature analysis (mean, standard deviation, skewness, kurtosis)

**Output Files:**
- `distribution_analysis.png` - Distribution analysis plot
- `qq_plots.png` - QQ plots

**Dependencies:** `param.py` + `compare_space_structure_shannon.py`

### 5. dimensionality_comparison.py
**Function:** Compare effects of 1D, 2D, 3D spatial structures on species abundance distributions
- Dimensionless parameter processing
- Multi-dimensional spatial simulation (with/without diffusion)
- Rank-abundance distribution analysis
- Spatial configuration comparison

**Output Files:**
- `species_abundance_boxplots.png` - Species abundance boxplots
- `median_abundance_trends.png` - Median abundance trend plots

**Dependencies:** `param.py`

### 6. distance_decay.py
**Function:** Analyze distance decay effects on community diversity
- Distance decay matrix generation
- Multiple diversity index calculations (Shannon, Simpson, Berger-Parker, etc.)
- Distribution fitting (Log-Normal vs Log-Series)
- Diffusion coefficient sensitivity analysis

**Output Files:**
- `rank_abundance_distributions.png` - Complete Rank-Abundance distributions
- `shannon_diversity.png` - Shannon diversity vs diffusion coefficient
- `simpson_diversity.png` - Simpson diversity vs diffusion coefficient
- `berger_parker_dominance.png` - Berger-Parker dominance vs diffusion coefficient
- `gini_coefficient.png` - Gini coefficient vs diffusion coefficient
- `shannon_evenness.png` - Shannon evenness vs diffusion coefficient
- `all_diversity_indices.png` - Combined plot of all diversity indices
- Multiple distribution fitting plots

**Dependencies:** `param.py`

## Environment Requirements

### Python Version
- Python 3.7+

### Required Dependencies
```bash
pip install numpy matplotlib seaborn pandas scipy
```

### Dependency Versions
- numpy >= 1.19.0
- matplotlib >= 3.3.0
- seaborn >= 0.11.0
- pandas >= 1.2.0
- scipy >= 1.6.0

## Installation and Execution

1. **Clone or download project files**
2. **Ensure all files are in the same directory**
3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
4. **Run code in sequence:**
   ```bash
   python compare_space_structure_shannon.py
   python frequency_boxplot_analysis.py
   python analyze_distribution.py
   python dimensionality_comparison.py
   python distance_decay.py
   ```

## Output Results

### Chart Files
- Distribution analysis plots
- Frequency distribution and boxplot plots
- Spatial dimension comparison plots
- Distance decay effect plots
- Diversity index analysis plots
- Distribution fitting plots

### Console Output
- Statistical summaries
- Normality test results
- Diversity index matrices
- Distribution fitting results
- Slope analysis results

## Important Notes

1. **File Dependencies:** `param.py` must be in the same directory as all code files
2. **Execution Order:** Follow the specified sequence strictly to avoid dependency errors
3. **Memory Requirements:** Some analyses may require significant memory, recommend 8GB+
4. **Computation Time:** Complete analysis may take several minutes to hours, depending on computer performance

## Technical Features

- **Dimensionless Processing:** All parameters are dimensionless for improved numerical stability
- **Multi-dimensional Support:** Supports 1D, 2D, 3D spatial structure simulation
- **Statistical Testing:** Multiple normality test methods
- **Distribution Fitting:** Log-Normal and Log-Series distribution fitting
- **Rich Visualization:** Multiple chart types with high-resolution output support

## Troubleshooting

### Common Errors
1. **ModuleNotFoundError: No module named 'param'**
   - Solution: Ensure `param.py` is in the current directory

2. **FileNotFoundError: compare_space_structure_shannon.py**
   - Solution: Ensure all dependency files exist

3. **Memory Insufficient Error**
   - Solution: Reduce `n_simulations` parameter or close other programs

### Performance Optimization
- Reduce `n_simulations` parameter value
- Adjust `t_span` and `t_eval` parameters
- Use faster numerical integration methods

## Contact Information

For questions or suggestions, please contact the project maintainer.

---

**Last Updated:** 2024
**Version:** 1.0 