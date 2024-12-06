# Clear the workspace by removing all objects
rm(list = ls())

# Load the dataset containing annual mean temperatures in Key West
load("../data/KeyWestAnnualMeanTemperature.RData") 

# List all objects in the environment to confirm data is loaded
ls()

# Check the class and structure of the dataset
class(ats) # Confirm the type of `ats` (e.g., data.frame)
head(ats) # Display the first few rows of the dataset
plot(ats) # Plot the dataset to visualize trends

# Calculate the actual correlation between Year and Temp
actual_corr <- cor(ats$Year, ats$Temp)

# Initialize parameters for permutation testing
set.seed(123)  # Set the random seed for reproducibility
n_permutations <- 10000 # Define the number of permutations
random_corrs <- numeric(n_permutations)  # Create a numeric vector to store random correlations

# Permutation test to calculate random correlations
for (i in 1:n_permutations) {
  shuffled_temp <- sample(ats$Temp)  # Shuffle the temperature values randomly
  random_corrs[i] <- cor(ats$Year, shuffled_temp)  # Compute the correlation with shuffled data
}

# Plot the histogram of random correlations
hist(random_corrs, 
     main = "Random Correlation Coefficients", # Set the title
     xlab = "Correlation") # Set the x-axis label

# Calculate the p-value for the two-tailed test
p_value <- mean(abs(random_corrs) >= abs(actual_corr))  # Proportion of random correlations greater than or equal to the actual correlation
print(abs(random_corrs) >= abs(actual_corr)) # Logical vector showing where the condition is met
print(paste("P-value:", p_value)) # Print the calculated p-value

# Save the histogram of random correlations to a PDF file
pdf("../results/random_correlation_histogram.pdf", width = 8, height = 6)

# Plot the histogram of random correlations
hist(random_corrs, 
     main = "Random Correlation Coefficients", # Set the title
     xlab = "Correlation", # Set the x-axis label
     xlim = range(c(random_corrs, actual_corr, -actual_corr)) + c(-0.1, 0.1),  # Extend the range for better visualization
     col = "lightblue", # Set the color of the bars
     breaks = 30)  # Set the number of bins

# Add vertical lines for the actual correlation and its symmetric counterpart
abline(v = actual_corr, col = "red", lwd = 2, lty = 2)  # Add a red dashed line for the actual correlation
abline(v = -actual_corr, col = "red", lwd = 2, lty = 2)  # Add a red dashed line for the symmetric correlation

# Add a legend to the plot
legend("topright",  # Place the legend in the top-right corner
       legend = c("Actual Correlation"), # Label the legend
       col = c("red"), # Set the color of the legend line
       lty = c(2), # Set the line type for the legend
       lwd = c(2), # Set the line width for the legend
       box.lty = 0)  # Remove the legend box border

# Close the PDF device to save the plot
dev.off()


