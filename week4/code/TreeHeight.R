# Function to calculate tree height given angle and distance
TreeHeight <- function(degrees, distance) {
  radians <- degrees * pi / 180
  height <- distance * tan(radians)
  return(height) # Return height, no need for print statements in batch processing
}

# Set working directory to the script's location (optional)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Use if you want relative paths

# Load the input data (relative path from script location)
input_file <- "../data/trees.csv" # Adjust based on your folder structure
trees <- read.csv(input_file, header = TRUE)

# Calculate tree heights for each row
trees$Tree.Height.m <- TreeHeight(trees$Angle.degrees, trees$Distance.m)

# Output the modified data to a new CSV file in the results folder
output_file <- "../results/TreeHts.csv"
write.csv(trees, file = output_file, row.names = FALSE)

# Confirm success
print("Tree heights have been successfully calculated and saved to TreeHts.csv!")

