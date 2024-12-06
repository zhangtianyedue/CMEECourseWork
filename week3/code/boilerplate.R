# A boilerplate R script

# Define a generic function that takes two arguments and prints their types
MyFunction <- function(Arg1, Arg2) {
  
  # Print the type (class) of Arg1
  print(paste("Argument", as.character(Arg1), "is a", class(Arg1))) 
  
  # Print the type (class) of Arg2
  print(paste("Argument", as.character(Arg2), "is a", class(Arg2))) 
  
  return (c(Arg1, Arg2)) # Return both arguments as a vector (optional but useful for testing)
}

# Test the function with numeric arguments
MyFunction(1, 2)

# Test the function with character arguments
MyFunction("Riki", "Tiki") 

# List all objects in the environment matching the pattern "MyFun*"
ls(pattern = "MyFun*") 

# Function to check if a number is even
is.even <- function(n = 2) {
  if (n %% 2 == 0) { # Check if n is divisible by 2
    return(paste(n, 'is even!')) # Return a message indicating n is even
  } else {
    return(paste(n, 'is odd!')) # Return a message indicating n is odd
  }
}

# Test the is.even function
is.even(6)

# Function to check if a number is a power of 2
is.power2 <- function(n = 2) {
  if (log2(n) %% 1 == 0) { # Check if log2(n) is an integer
    return(paste(n, 'is a power of 2!')) # Return a message indicating n is a power of 2
  } else {
    return(paste(n, 'is not a power of 2!')) # Return a message indicating n is not a power of 2
  }
}

# Test the is.power2 function
is.power2(4)

# Function to check if a number is prime
is.prime <- function(n) {
  if (n == 0) { # Check if n is 0
    return(paste(n, 'is a zero!')) # Return a message indicating n is zero
  } else if (n == 1) { # Check if n is 1
    return(paste(n, 'is just a unit!')) # Return a message indicating n is a unit
  }
  
  ints <- 2:(n-1) # Create a sequence of integers from 2 to n-1
  
  if (all(n %% ints != 0)) { # Check if n is not divisible by any integers in the sequence
    return(paste(n, 'is a prime!')) # Return a message indicating n is prime
  } else {
    return(paste(n, 'is a composite!')) # Return a message indicating n is composite
  }
}

# Test the is.prime function
is.prime(3)
