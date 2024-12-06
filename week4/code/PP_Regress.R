# Clear the workspace to ensure no leftover variables from previous runs
rm(list = ls())

# Import necessary libraries for data manipulation and visualization
library(ggplot2) # For creating elegant data visualizations
library(dplyr)   # For efficient data manipulation

# Load the dataset containing predator-prey interaction data
MyDF <- read.csv("../data/EcolArchives-E089-51-D1.csv")

# Convert prey mass from milligrams (mg) to grams (g) where applicable
# This ensures consistency in the unit of measurement across the dataset
#MyDF <- MyDF %>%
  #mutate(Prey.mass = ifelse(Prey.mass.unit == "mg", Prey.mass / 1000, Prey.mass),
         #Prey.mass.unit = "g") # Update the unit to grams for all entries

# Create a PDF file to save the resulting plot for documentation and sharing purposes
pdf("../results/PP_Regress.pdf", width = 9, height = 12) # Set dimensions for better readability

# Generate a scatter plot to visualize predator and prey mass by feeding type and predator lifestage
# Logarithmic scales are used for better visualization of mass relationships
p <- ggplot(MyDF, aes(x = Prey.mass, y = Predator.mass, color = Predator.lifestage)) +
  geom_point(shape = 3) + # Add scatter points to the plot
  geom_smooth(method = "lm", formula = y ~ x, fullrange = TRUE, na.rm = TRUE) + # Add linear regression lines
  scale_x_log10() + # Logarithmic scale for x-axis (Prey mass)
  scale_y_log10() + # Logarithmic scale for y-axis (Predator mass)
  xlab("Prey mass in grams") + # Label for x-axis
  ylab("Predator mass in grams") + # Label for y-axis
  facet_grid(Type.of.feeding.interaction ~ .) + # Facet by feeding type for detailed breakdown
  theme_bw() + # Use a clean, black-and-white theme
  theme(legend.position = "bottom", # Place legend below the plot
        panel.border = element_rect(colour = "grey", fill = NA), # Add grey borders around panels
        legend.title = element_text(size = 9, face = "bold")) + # Bold legend title for emphasis
  guides(colour = guide_legend(nrow = 1)) # Format legend to display in a single row

# Render the plot and save it to the PDF file
print(p)

# Close the PDF device to finalize the file
graphics.off()

# Calculate regression results corresponding to the linear models fitted in the figure
LM <- MyDF %>%
  # Remove specific subsets with only 2 examples that might skew results
  # These subsets contain identical species of prey and predator
  filter(Record.number != "30914" & Record.number != "30929") %>%
  # Select relevant columns needed for regression analysis
  dplyr::select(Record.number, Predator.mass, Prey.mass, Predator.lifestage, Type.of.feeding.interaction) %>%
  # Group data by feeding type and predator lifestage for individual regression analyses
  group_by(Type.of.feeding.interaction, Predator.lifestage) %>%
  # Perform linear regression for each group and store models temporarily
  do(mod = lm(Predator.mass ~ Prey.mass, data = .)) %>%
  # Extract key regression statistics and add them as new columns in the dataframe
  mutate(
    Regression.slope = summary(mod)$coefficients[2, 1],       # Slope of the regression line
    Regression.intercept = summary(mod)$coefficients[1, 1],  # Intercept of the regression line
    R.squared = summary(mod)$adj.r.squared,                  # Adjusted R-squared value for model fit
    Fstatistic = summary(mod)$fstatistic[1],                 # F-statistic for the model
    P.value = summary(mod)$coefficients[2, 4]                # P-value for the slope coefficient
  ) %>%
  # Remove the temporary model object to clean up the dataframe
  dplyr::select(-mod)

# Save the regression results to a CSV file for further analysis or reporting
write.csv(LM, "../results/PP_Regress_Results.csv", row.names = FALSE)
