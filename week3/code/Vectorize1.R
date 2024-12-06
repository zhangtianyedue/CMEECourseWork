# Create a 1000 x 1000 matrix filled with random numbers between 0 and 1
M <- matrix(runif(1000000), 1000, 1000)

# Function to calculate the sum of all elements in a matrix using nested loops
SumAllElements <- function(M) {
  Dimensions <- dim(M) # Get the dimensions of the matrix (number of rows and columns)
  Tot <- 0 # Initialize the total sum to 0
  for (i in 1:Dimensions[1]) { # Loop over each row
    for (j in 1:Dimensions[2]) { # Loop over each column
      Tot <- Tot + M[i, j] # Add the current element to the total
    }
  }
  return(Tot) # Return the total sum
}

# Measure and print the time taken to compute the sum using the loop-based function
print("Using loops, the time taken is:")
print(system.time(SumAllElements(M))) # Measure execution time for `SumAllElements`

# Measure and print the time taken to compute the sum using the built-in vectorized `sum()` function
print("Using the in-built vectorized function, the time taken is:")
print(system.time(sum(M))) # Measure execution time for `sum()`

# Function to demonstrate the inefficiency of not preallocating a vector
NoPreallocFun <- function(x) {
  a <- vector() # Create an empty vector
  for (i in 1:x) { # Loop from 1 to x
    a <- c(a, i) # Concatenate the current number to the vector
    print(a) # Print the current state of the vector
    print(object.size(a)) # Print the memory size of the vector
  }
}

# Measure the execution time of the function without preallocation
system.time(NoPreallocFun(10)) # Test the function with x = 10

# Function to demonstrate the efficiency of preallocating a vector
PreallocFun <- function(x) {
  a <- rep(NA, x) # Preallocate a vector of size x with `NA` values
  for (i in 1:x) { # Loop from 1 to x
    a[i] <- i # Assign the current number to the ith position in the vector
    print(a) # Print the current state of the vector
    print(object.size(a)) # Print the memory size of the vector
  }
}

# Measure the execution time of the function with preallocation
system.time(PreallocFun(10)) # Test the function with x = 10

