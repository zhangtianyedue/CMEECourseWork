library(ggplot2) # Load the ggplot2 library for data visualization

# Define the build_ellipse function
build_ellipse <- function(hradius, vradius) {
  # This function creates a dataframe representing an ellipse
  # hradius: horizontal radius of the ellipse
  # vradius: vertical radius of the ellipse
  
  npoints <- 250 # Number of points used to construct the ellipse
  a <- seq(0, 2 * pi, length = npoints + 1) # Generate angles from 0 to 2π
  x <- hradius * cos(a) # Compute x-coordinates of the ellipse
  y <- vradius * sin(a) # Compute y-coordinates of the ellipse
  return(data.frame(x = x, y = y)) # Return a dataframe of the ellipse coordinates
}

# Assign the size of the square matrix
N <- 250 # Number of rows and columns in the matrix

# Build the matrix with random values from a normal distribution
M <- matrix(rnorm(N * N), N, N) # Create an N x N matrix with random values

# Compute the eigenvalues of the matrix
eigvals <- eigen(M)$values # Extract eigenvalues from the matrix

# Create a dataframe for eigenvalues
eigDF <- data.frame("Real" = Re(eigvals), "Imaginary" = Im(eigvals)) 
# Real parts and imaginary parts of eigenvalues are stored as separate columns

# Define the radius of the Girko circle
my_radius <- sqrt(N) # Theoretical radius of the Girko circle

# Create the ellipse dataframe using the build_ellipse function
ellDF <- build_ellipse(my_radius, my_radius) 

# Rename the columns of the ellipse dataframe for consistency
names(ellDF) <- c("Real", "Imaginary") 

# Create the ggplot object
p <- ggplot(eigDF, aes(x = Real, y = Imaginary)) + # Plot eigenvalues
  geom_point(shape = 3) +  # Add points for eigenvalues
  geom_hline(aes(yintercept = 0)) +  # Add a horizontal line at y = 0
  geom_vline(aes(xintercept = 0)) +  # Add a vertical line at x = 0
  geom_polygon(data = ellDF, aes(x = Real, y = Imaginary, alpha = 1/20, fill = "red")) + 
  # Add the Girko circle as a polygon with transparency and fill color
  theme(legend.position = "none")  # Remove the legend for clarity

# Save the plot as a PDF file
pdf("../results/Girko.pdf", width = 8, height = 6)  # Open a PDF device with specified dimensions
print(p)  # Render the plot to the PDF
dev.off()  # Close the PDF device to save the file

