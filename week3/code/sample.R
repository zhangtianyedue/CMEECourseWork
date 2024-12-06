######### Functions ##########

## A function to take a sample of size n from a population "popn" and return its mean
myexperiment <- function(popn, n) {
  # Randomly sample `n` elements from the population `popn` without replacement
  pop_sample <- sample(popn, n, replace = FALSE)
  return(mean(pop_sample)) # Return the mean of the sampled elements
}

## Calculate means using a FOR loop on a vector without preallocation
loopy_sample1 <- function(popn, n, num) {
  result1 <- vector() # Initialize an empty vector
  for(i in 1:num) { # Loop through `num` iterations
    # Append the mean of a random sample to the vector
    result1 <- c(result1, myexperiment(popn, n))
  }
  return(result1) # Return the resulting vector of means
}

## To run "num" iterations of the experiment using a FOR loop on a vector with preallocation
loopy_sample2 <- function(popn, n, num) {
  result2 <- vector(, num) # Preallocate a vector of size `num`
  for(i in 1:num) { # Loop through `num` iterations
    result2[i] <- myexperiment(popn, n) # Assign the mean of a random sample to the ith position
  }
  return(result2) # Return the resulting vector of means
}

## To run "num" iterations of the experiment using a FOR loop on a list with preallocation
loopy_sample3 <- function(popn, n, num) {
  result3 <- vector("list", num) # Preallocate a list of size `num`
  for(i in 1:num) { # Loop through `num` iterations
    result3[[i]] <- myexperiment(popn, n) # Assign the mean of a random sample to the ith position in the list
  }
  return(result3) # Return the resulting list of means
}

## To run "num" iterations of the experiment using vectorization with lapply
lapply_sample <- function(popn, n, num) {
  # Use lapply to run the experiment `num` times
  # This applies the anonymous function `function(i) myexperiment(popn, n)` to each iteration
  result4 <- lapply(1:num, function(i) myexperiment(popn, n))
  return(result4) # Return the resulting list of means
}

## To run "num" iterations of the experiment using vectorization with sapply
sapply_sample <- function(popn, n, num) {
  # Use sapply to run the experiment `num` times
  # This is similar to lapply but returns a vector instead of a list
  result5 <- sapply(1:num, function(i) myexperiment(popn, n))
  return(result5) # Return the resulting vector of means
}

# Set a random seed for reproducibility
set.seed(12345)

# Generate a population of 10,000 random numbers from a normal distribution
popn <- rnorm(10000)

# Plot a histogram of the population distribution
hist(popn)

# Define the sample size for each experiment and the number of iterations
n <- 100 # Number of elements in each sample
num <- 10000 # Number of experiments to run

# Test and time each method
print("Using loops without preallocation on a vector took:" )
print(system.time(loopy_sample1(popn, n, num))) # Time the method without preallocation

print("Using loops with preallocation on a vector took:" )
print(system.time(loopy_sample2(popn, n, num))) # Time the method with preallocation on a vector

print("Using loops with preallocation on a list took:" )
print(system.time(loopy_sample3(popn, n, num))) # Time the method with preallocation on a list

print("Using the vectorized sapply function (on a list) took:" )
print(system.time(sapply_sample(popn, n, num))) # Time the method using sapply

print("Using the vectorized lapply function (on a list) took:" )
print(system.time(lapply_sample(popn, n, num))) # Time the method using lapply

