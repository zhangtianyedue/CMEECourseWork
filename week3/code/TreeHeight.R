# Define the TreeHeight function to calculate the height of a tree
TreeHeight <- function(degrees, distance) {
    # Convert the angle from degrees to radians
    radians <- degrees * pi / 180 # Formula: radians = degrees × π / 180
    
    # Calculate the height of the tree using the trigonometric formula
    # height = distance × tan(radians)
    height <- distance * tan(radians)
    
    # Print the calculated height to the console
    print(paste("Tree height is:", height))
  
    # Return the calculated height as the function output
    return(height)
}

# Example usage: Calculate the height of a tree
# with an angle of elevation of 37 degrees and a distance of 40 meters
TreeHeight(37, 40)
