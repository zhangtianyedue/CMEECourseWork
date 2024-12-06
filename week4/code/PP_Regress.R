# Define the list of required packages
required_packages <- c("ggplot2", "dplyr")

# Check if each package is installed, and install it if it's missing
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE) # Install the package
    library(pkg, character.only = TRUE) # Load the package
  } else {
    library(pkg, character.only = TRUE) # Load the package if already installed
  }
}
# Clear the workspace
rm(list = ls()) 
# Remove all objects in the current R session to start fresh

# Import necessary libraries
library(ggplot2) # For data visualization
library(dplyr) # For data manipulation

# Load the data
MyDF <- read.csv("../data/EcolArchives-E089-51-D1.csv") 
# Load the dataset from a CSV file

# Convert masses from mg to g
MyDF <- MyDF %>%
  mutate(
    Prey.mass = ifelse(Prey.mass.unit == "mg", Prey.mass / 1000, Prey.mass), 
    # Convert Prey mass to grams if the unit is mg
    Prey.mass.unit = "g" 
    # Update the unit column to "g" for consistency
  )

# Create a PDF to save the plot
pdf("../results/PP_Regress.pdf", width = 9.5, height = 12) 
# Open a PDF graphics device to save the plot

# Plot predator and prey mass by feeding type and predator lifestage
p <- ggplot(MyDF, aes(
  x = Prey.mass, 
  y = Predator.mass, 
  color = Predator.lifestage
)) +
  geom_point(shape = 3) + 
  # Add points to the plot, using "+" to layer graphical elements
  geom_smooth(
    method = "lm", 
    formula = y ~ x, 
    fullrange = TRUE, 
    na.rm = TRUE
  ) + 
  # Add linear regression lines for each group
  scale_x_log10() + 
  # Log-transform the x-axis
  scale_y_log10() + 
  # Log-transform the y-axis
  xlab("Prey mass in grams") + 
  # Label the x-axis
  ylab("Predator mass in grams") + 
  # Label the y-axis
  facet_grid(Type.of.feeding.interaction ~ .) + 
  # Create separate plots for each feeding interaction type
  theme_bw() + 
  # Use a black-and-white theme
  theme(
    legend.position = "bottom", 
    # Position the legend at the bottom
    panel.border = element_rect(colour = "grey", fill = NA), 
    # Add a grey border to the panels
    legend.title = element_text(size = 9, face = "bold")
    # Format the legend title
  ) +
  guides(colour = guide_legend(nrow = 1)) 
  # Display the legend in a single row

print(p) # Render the plot

# Close the PDF device to save the plot
graphics.off()

# Check and remove duplicate rows in the data
MyDF <- MyDF %>%
  distinct() 
  # Remove completely duplicate rows from the dataframe

# Identify duplicates based on specific columns
duplicates <- MyDF %>%
  group_by(Record.number) %>%
  filter(n() > 1) 
  # Identify rows with duplicate Record.number values

# Handle duplicates if found
if (nrow(duplicates) > 0) {
  cat("以下记录是重复的，将被移除：\n")
  # Print a message if duplicates are found
  print(duplicates) 
  # Display the duplicate rows
  
  # Keep only the first occurrence of each Record.number
  MyDF <- MyDF %>%
    distinct(Record.number, .keep_all = TRUE)
} else {
  cat("没有发现重复数据。\n") 
  # Print a message if no duplicates are found
}

# Calculate regression results for predator-prey relationships
LM <- MyDF %>%
  dplyr::select(
    Record.number, 
    Predator.mass, 
    Prey.mass, 
    Predator.lifestage, 
    Type.of.feeding.interaction
  ) %>%
  # Select relevant columns for analysis
  group_by(Type.of.feeding.interaction, Predator.lifestage) %>%
  # Group the data by feeding interaction type and predator lifestage
  do(mod = lm(Predator.mass ~ Prey.mass, data = .)) %>%
  # Fit linear models for each group
  mutate(
    Regression.slope = summary(mod)$coefficients[2, 1], 
    # Extract the slope of the regression line
    Regression.intercept = summary(mod)$coefficients[1, 1], 
    # Extract the intercept of the regression line
    R.squared = summary(mod)$adj.r.squared, 
    # Extract the adjusted R-squared value
    Fstatistic = summary(mod)$fstatistic[1], 
    # Extract the F-statistic
    P.value = summary(mod)$coefficients[2, 4]
    # Extract the p-value for the slope
  ) %>%
  dplyr::select(-mod) 
  # Remove the model object to keep only the summary data

# Save the regression results to a CSV file
write.csv(LM, "../results/PP_Regress_Results1.csv", row.names = FALSE)
# Write the regression summary to a CSV file without row names
