# Clear workspace
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(ggplot2)
# Load weather data
load("../data/KeyWestAnnualMeanTemperature.RData")

# Plot temperatures over time on line graph
pdf("../results/Temperature-Group.pdf")
plot(ats$Year, 
     ats$Temp, 
     xlab = "Year",
     ylab = "Temperature (°C)",
     type = "l",
     main = "Annual Mean Temperature in Key West, Florida (1901-2000)")
dev.off()

# Create two vectors of temperatures, one with the first row deleted to align for comparison
Temp_t0 <- ats$Temp[2:100]
Temp_t1 <- ats$Temp[1:99]

# Calculate the correlation between successive years
CorCoeff <- cor(Temp_t0, Temp_t1)
cat("Correlation between successive years is", CorCoeff, "\n")

# Create a matrix of 10,000 random permutations of temperature column
set.seed(123)
Temps1 <- replicate(10000, sample(ats$Temp, replace = FALSE))

# For each permutation, realign as before and calculate correlation
RdmCors <- apply(Temps1, 2, function(x) cor(x[2:100], x[1:99]))

# Generate histogram comparing p-values
pdf("../results/Temperature-Group_corr.pdf")
hist(RdmCors, 
     xlim = c(-0.4, 0.4),
     xlab = "Correlation Coefficients",
     main = "Histogram of Random Correlations")
abline(v = CorCoeff, col = "red", lwd = 2, lty = 2)
text(CorCoeff - 0.004, 1500, paste("Correlation coefficient \nfor successive years: ", round(CorCoeff, 3), sep = ""), pos = 2, col = "red", cex = 0.8)
dev.off()

# Combine vectors into a data frame

temp_data <- data.frame(Temp_t1 = Temp_t1, Temp_t0 = Temp_t0)

# Plot temperature correlation between successive and previous years using ggplot2
pdf("../results/Temperature_Correlation.pdf", width = 7, height = 5)
ggplot(temp_data, aes(x = Temp_t1, y = Temp_t0)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", col = "blue") +
  labs(x = "Previous year (degrees)", y = "Successive years (degrees)",
       title = "Temperature Correlation Between Successive and Previous Years") +
  annotate("text", x = min(Temp_t1) + 0.1, y = max(Temp_t0) - 0.2,
           label = paste("Observed correlation coefficient:", round(CorCoeff, 3)),
           color = "red", size = 4, hjust = 0) +
  theme_minimal()
dev.off()


# Calculate estimated p-value (fraction of correlation coefficients more extreme than CorCoeff)
p_estimate <- mean(abs(RdmCors) >= abs(CorCoeff))

# Print p_estimate to screen
cat("P-value estimate is", p_estimate, "\n")