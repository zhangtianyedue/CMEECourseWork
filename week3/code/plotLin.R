# Load necessary library
library(ggplot2) # Load ggplot2 for data visualization

# Generate a sequence of x values from 0 to 100 with a step of 0.1
x <- seq(0, 100, by = 0.1)

# Generate y values based on a linear relationship with added random noise
y <- -4. + 0.25 * x + # Linear equation: y = -4 + 0.25 * x
  rnorm(length(x), mean = 0., sd = 2.5) # Add random noise with mean 0 and standard deviation 2.5

# Put the x and y values into a dataframe for analysis
my_data <- data.frame(x = x, y = y)

# Perform a linear regression on y as a function of x and summarize the results
my_lm <- summary(lm(y ~ x, data = my_data)) 

# Plot the data
p <- ggplot(my_data, aes(
  x = x, y = y, # Define x and y variables
  colour = abs(my_lm$residual) # Colour points based on the absolute residuals from the regression
)) +
  geom_point() + # Add points to the plot
  scale_colour_gradient(low = "black", high = "red") + # Use a colour gradient from black to red
  theme(legend.position = "none") + # Remove the legend from the plot
  scale_x_continuous(expression(alpha^2 * pi / beta * sqrt(Theta))) 
  # Add a mathematical expression as the x-axis label

# Add the regression line to the plot
p <- p + geom_abline(
  intercept = my_lm$coefficients[1][1], # Use the intercept from the linear regression
  slope = my_lm$coefficients[2][1], # Use the slope from the linear regression
  colour = "red" # Set the line colour to red
)

# Add a mathematical annotation to the plot
p <- p + geom_text(
  aes(x = 60, y = 0, # Specify the position of the text on the plot
      label = "sqrt(alpha) * 2* pi"), # Add a LaTeX-style label
  parse = TRUE, size = 6, # Parse the label as a mathematical expression and set text size
  colour = "blue" # Set the text colour to blue
)

# Display the plot
p

# Save the plot as a PDF in the "../results/" directory
pdf("../results/MyLinReg.pdf", width = 8, height = 6)  # Open a PDF device with specified dimensions
print(p)  # Render the plot and save it to the file
dev.off()  # Close the PDF device to finalize the file
