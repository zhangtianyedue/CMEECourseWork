# Runs the stochastic Ricker equation with Gaussian fluctuations

rm(list = ls()) # Clear the workspace to remove any existing objects

# Function to simulate population dynamics using the stochastic Ricker model
stochrick <- function(p0 = runif(1000, 0.5, 1.5), r = 1.2, K = 1, sigma = 0.2, numyears = 100) {
  
  # Initialize an empty matrix to store population sizes
  N <- matrix(NA, numyears, length(p0)) # Rows represent years, columns represent populations
  
  # Set the initial population sizes
  N[1, ] <- p0
  
  # Loop through each population
  for (pop in 1:length(p0)) { 
    
    # For each population, loop through the years
    for (yr in 2:numyears) {
      
      # Calculate the population size for the current year using the Ricker model
      # Add a Gaussian random fluctuation (with mean = 0 and standard deviation = sigma)
      N[yr, pop] <- N[yr - 1, pop] * exp(r * (1 - N[yr - 1, pop] / K) + rnorm(1, 0, sigma))
    }
  }
  return(N) # Return the matrix of population sizes
}

# Vectorized implementation of the stochastic Ricker model for improved performance
stochrickvect <- function(p0 = runif(1000, 0.5, 1.5), r = 1.2, K = 1, sigma = 0.2, numyears = 100) {
  
  # Initialize an empty matrix to store population sizes
  N <- matrix(NA, nrow = numyears, ncol = length(p0))
  
  # Set the initial population sizes
  N[1, ] <- p0
  
  # Pre-generate a matrix of random fluctuations for all years and populations
  random_fluctuations <- matrix(
    rnorm((numyears - 1) * length(p0), mean = 0, sd = sigma), 
    nrow = numyears - 1, ncol = length(p0)
  )
  
  # Use a vectorized loop to calculate population dynamics year by year
  for (yr in 2:numyears) {
    # Calculate population sizes for all populations in the current year
    N[yr, ] <- N[yr - 1, ] * exp(r * (1 - N[yr - 1, ] / K) + random_fluctuations[yr - 1, ])
  }
  
  return(N) # Return the matrix of population sizes
}

# Measure and print the runtime of the vectorized implementation
print("Vectorized Stochastic Ricker takes:")
print(system.time(res2 <- stochrickvect())) # Test the performance of the vectorized function

# Measure and print the runtime of the original loop-based implementation
print("Original Stochastic Ricker takes:")
print(system.time(res1 <- stochrick())) # Test the performance of the original function

# Verify that the results from both implementations are identical
print("Results are identical:")
print(all.equal(res1, res2)) # Check if the two result matrices are equal


