Ricker <- function(N0 = 1, r = 1, K = 10, generations = 50) {
  # Simulates population dynamics using the Ricker model
  # Arguments:
  #   N0: Initial population size (default is 1)
  #   r: Intrinsic growth rate (default is 1)
  #   K: Carrying capacity of the environment (default is 10)
  #   generations: Number of generations to simulate (default is 50)
  # Returns:
  #   A vector of population sizes over the specified number of generations
  
  N <- rep(NA, generations)  # Create a vector of length `generations`, initialized with NA
  
  N[1] <- N0  # Set the initial population size as the first value in the vector
  for (t in 2:generations) {
    # Calculate the population size for generation t based on the Ricker model equation
    N[t] <- N[t-1] * exp(r * (1.0 - (N[t-1] / K)))
    # The formula models population growth with density dependence
  }
  return(N)  # Return the vector of population sizes
}

# Plot the results of the Ricker model simulation
plot(
  Ricker(generations = 10),  # Call the Ricker function for 10 generations
  type = "l"                 # Use a line plot to visualize population changes over time
)


