Exponential <- function(N0 = 1, r = 1, generations = 10) {
  # Simulates exponential growth over a given number of generations
  # Arguments:
  #   N0: Initial population size (default is 1)
  #   r: Growth rate (default is 1)
  #   generations: Number of generations to simulate (default is 10)
  # Returns:
  #   A vector of population sizes for each generation
  
  N <- rep(NA, generations)    # Create a vector of length `generations`, initialized with NA values
  
  N[1] <- N0                   # Set the initial population size as the first element of the vector
  for (t in 2:generations) {   # Loop through generations from 2 to the specified number
    N[t] <- N[t-1] * exp(r)    # Calculate the population size for generation `t` using exponential growth formula
    browser()                  # Enter debug mode at this point for each generation
  }
  return(N)                    # Return the vector of population sizes
}

# Plot the exponential growth simulation
plot(Exponential(),            # Call the `Exponential` function with default arguments
     type = "l",               # Set plot type to "line"
     main = "Exponential growth") # Add a title to the plot
