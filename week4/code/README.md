***WEEK4***
First and foremost, please note that I have used relative paths in this coding assignment. Therefore, to ensure the code runs correctly, make sure to set your working directory to `week4/code`. Otherwise, the files stored in the `data` folder will not be accessible.


In the file `PP_Regress.R`, I utilized packages such as **ggplot2** and **dplyr** for data visualization and manipulation. To ensure the smooth execution of the code, I included a section at the very beginning of the script that checks whether these required packages are installed on your system. If any of the necessary packages are missing, the script will automatically install them before proceeding. This approach ensures that all dependencies are resolved, allowing the subsequent code to run without issues. 

Please note that within the `code` folder, I have included a LaTeX report related to the `Florida.R` script. The report provides detailed analysis and visualizations generated from the script. 

The images in the report are referenced using **relative paths**, ensuring consistency with the folder structure. Therefore, it is crucial to set your working directory to `week4/code` when compiling the LaTeX document. Failing to do so may result in the images not being displayed correctly in the compiled report.

Ensure you have the necessary LaTeX setup to compile the document and verify that all required image files are present in the appropriate directories.

### ！！！！！！！！！！！！！: Data Unit Standardization in `PP_Regress.R`

In the `PP_Regress.R` file, I aimed to fit the data to match the given plot as required in the assignment. Using the program I wrote, I successfully produced the same result as shown in the provided example. However, I believe there is a slight issue with the assignment instructions. Specifically, the units of measurement were not standardized before performing the fitting process. 

To ensure data consistency, I argue that all units (e.g., milligrams and grams) should first be converted into a common unit, such as grams, before proceeding with the analysis. This can be achieved by adding the following lines to the code:

```R
# Standardize units to grams
MyDF <- MyDF %>%
  mutate(Prey.mass = ifelse(Prey.mass.unit == "mg", Prey.mass / 1000, Prey.mass),
         Prey.mass.unit = "g") # Update the unit to grams for all entries
```

When this adjustment is applied, the resulting plot—specifically, the third graph from the top—shows fewer points on the right-hand side. This discrepancy occurs because data originally in milligrams have much larger numerical values compared to grams, leading to this visible effect. 

However, if the above code block is commented out, the final output matches the expected results exactly as required in the assignment. Despite this, I maintain that it is good practice to standardize units during data processing to avoid potential biases or inconsistencies. In this case, the impact is minimal because the number of milligram entries is much smaller than those in grams, so the discrepancy does not significantly affect the overall results.

By ensuring proper unit standardization, we uphold data integrity and reduce the risk of misleading results in future analyses.