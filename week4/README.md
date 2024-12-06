My CMEE Coursework Repository

***WEEK4***
First and foremost, please note that I have used relative paths in this coding assignment. Therefore, to ensure the code runs correctly, make sure to set your working directory to `week4/code`. Otherwise, the files stored in the `data` folder will not be accessible.

## README

This folder contains four R scripts designed for different analytical tasks. Below is a brief description of each script, including its purpose, required packages, and usage.

### 1. **Florida.R**
- **Purpose**: This script analyzes temperature data from Florida, focusing on trends and correlations over time. It includes permutation testing to assess the statistical significance of observed correlations.
- **Required Packages**: None explicitly, but the script relies on standard base R functions for statistical analysis and plotting.
- **Key Features**:
  - Computes the actual correlation between year and temperature.
  - Performs permutation tests to generate a null distribution of random correlations.
  - Visualizes the results using a histogram and marks the observed correlation on the plot.
- **Usage**:
  Ensure the dataset `KeyWestAnnualMeanTemperature.RData` is stored in the `data` folder and set your working directory to `week4/code`.

---

### 2. **PP_Regress.R**
- **Purpose**: This script explores the relationship between predator and prey masses using regression analysis. It includes visualizations and calculates regression statistics for different feeding interaction types and predator life stages.
- **Required Packages**:
  - `ggplot2`: For visualizing predator-prey relationships.
  - `dplyr`: For data manipulation.
- **Key Features**:
  - Automatically checks for and installs missing packages.
  - Converts prey mass from milligrams to grams for consistency.
  - Creates a PDF visualization of predator-prey mass relationships with regression lines.
  - Saves regression results (slope, intercept, R-squared, F-statistic, p-value) as a CSV file.
- **Usage**:
  Ensure the `EcolArchives-E089-51-D1.csv` dataset is in the `data` folder, and set the working directory to `week4/code`.

---

### 3. **TAutoCorr.R**
- **Purpose**: This script performs time-series analysis to explore temporal autocorrelation in temperature or related data.
- **Required Packages**:
  - `ggplot2`: For visualizing temporal trends and autocorrelation results.
  - `dplyr`: For data manipulation.
- **Key Features**:
  - Computes temporal autocorrelation metrics.
  - Visualizes trends and autocorrelation with ggplot2.
  - Automatically installs missing packages if necessary.
- **Usage**:
  Place the relevant time-series dataset in the `data` folder, and set the working directory to `week4/code`.

---

### 4. **TreeHeight.R**
- **Purpose**: This script calculates the height of trees based on the angle of elevation and the distance from the tree base using trigonometric principles.
- **Required Packages**: None.
- **Key Features**:
  - Uses the formula: `height = distance * tan(angle in radians)` to compute tree heights.
  - Allows for quick and accurate height calculations for field data.
  - Includes example usage for testing.
- **Usage**:
  This script requires no additional data files. Simply run the script to calculate tree heights for given angles and distances.


### Latex routine Warning 

Please note that within the `code` folder, I have included a LaTeX report related to the `Florida.R` script. The report provides detailed analysis and visualizations generated from the script. 

The images in the report are referenced using **relative paths**, ensuring consistency with the folder structure. Therefore, it is crucial to set your working directory to `week4/code` when compiling the LaTeX document. Failing to do so may result in the images not being displayed correctly in the compiled report.

Ensure you have the necessary LaTeX setup to compile the document and verify that all required image files are present in the appropriate directories.