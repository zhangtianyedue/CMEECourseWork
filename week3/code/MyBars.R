# Load necessary library
library(ggplot2) # Load ggplot2 for data visualization

# Load the dataset from a text file located in the "../data/" directory
a <- read.table("../data/Results.txt", header = TRUE) # Read the data with headers
head(a) # Display the first few rows of the dataset to inspect its structure

# Add a new column "ymin" initialized to zeros
a$ymin <- rep(0, dim(a)[1]) # Set the ymin column to 0 for all rows

# Initialize the ggplot object
p <- ggplot(a)

# Add the first set of line ranges to the plot
p <- p + geom_linerange(
  data = a, 
  aes(x = x, ymin = ymin, ymax = y1, size = 0.5), # Define x, ymin, and ymax for the range
  colour = "#E69F00", # Set the color of the lines
  alpha = 1/2, # Adjust transparency
  show.legend = FALSE # Hide the legend
)

# Add the second set of line ranges to the plot
p <- p + geom_linerange(
  data = a, 
  aes(x = x, ymin = ymin, ymax = y2, size = 0.5), # Define x, ymin, and ymax for the second range
  colour = "#56B4E9", # Set a different color
  alpha = 1/2, # Adjust transparency
  show.legend = FALSE # Hide the legend
)

# Add the third set of line ranges to the plot
p <- p + geom_linerange(
  data = a, 
  aes(x = x, ymin = ymin, ymax = y3, size = 0.5), # Define x, ymin, and ymax for the third range
  colour = "#D55E00", # Set a different color
  alpha = 1/2, # Adjust transparency
  show.legend = FALSE # Hide the legend
)

# Add labels to the plot
p <- p + geom_text(
  data = a, 
  aes(x = x, y = -500, label = Label) # Annotate the plot with labels from the "Label" column
)

# Customize axis labels, ticks, and remove legend
p <- p + scale_x_continuous(
  "My x axis", # Set x-axis label
  breaks = seq(3, 5, by = 0.05) # Define x-axis tick marks
) + 
  scale_y_continuous("My y axis") + # Set y-axis label
  theme_bw() + # Use a black-and-white theme for better printing
  theme(legend.position = "none") # Remove the legend from the plot

# Display the plot
p

# Save the plot to the "../results/" directory as MyBars.pdf
pdf("../results/MyBars.pdf", width = 8, height = 6)  # Open a PDF device with specified dimensions
print(p)  # Render and save the plot
dev.off()  # Close the PDF device to finalize the file
