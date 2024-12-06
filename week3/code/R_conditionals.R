rm(list = ls()) 
# Clear the workspace by removing all objects in the current R environment.

# Function to check if an integer is even
is.even <- function(n = 2) {
  # Check if the remainder of n divided by 2 is zero
  if (n %% 2 == 0) {
    return(paste(n, 'is even!')) # Return a message indicating n is even
  } else {
    return(paste(n, 'is odd!')) # Return a message indicating n is odd
  }
}

# Test the is.even function with input 6
is.even(6) # Evaluate the function but do not print the result
print(is.even(6)) # Print the result of the function

# Function to check if a number is a power of 2
is.power2 <- function(n = 2) {
  # Check if the base-2 logarithm of n is an integer
  if (log2(n) %% 1 == 0) {
    return(paste(n, 'is a power of 2!')) # Return a message indicating n is a power of 2
  } else {
    return(paste(n, 'is not a power of 2!')) # Return a message indicating n is not a power of 2
  }
}

# Test the is.power2 function with input 4
is.power2(4) # Evaluate the function but do not print the result
print(is.power2(4)) # Print the result of the function

# Function to check if a number is prime
is.prime <- function(n) {
  # Handle special cases for n = 0 and n = 1
  if (n == 0) {
    return(paste(n, 'is a zero!')) # Return a message for n = 0
  } else if (n == 1) {
    return(paste(n, 'is just a unit!')) # Return a message for n = 1
  }
  
  # Create a sequence of integers from 2 to n-1
  ints <- 2:(n - 1)
  
  # Check if n is not divisible by any integers in the sequence
  if (all(n %% ints != 0)) {
    return(paste(n, 'is a prime!')) # Return a message indicating n is a prime number
  } else {
    return(paste(n, 'is a composite!')) # Return a message indicating n is a composite number
  }
}

# Test the is.prime function with input 3
is.prime(3) # Evaluate the function but do not print the result
print(is.prime(3)) # Print the result of the function
