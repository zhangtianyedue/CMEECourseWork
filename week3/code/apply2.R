# Define a custom function SomeOperation
SomeOperation <- function(v) { 
  # This function checks if the sum of the input vector `v` is greater than 0.
  if (sum(v) > 0) { # If the sum of elements in `v` is positive
    return(v * 100) # Multiply each element of `v` by 100 and return
  } else { 
    return(v) # Otherwise, return `v` as is
  }
}

# Create a 10x10 matrix with random values from a normal distribution
M <- matrix(rnorm(100), 10, 10)

# Apply the SomeOperation function to each row of the matrix `M`
# `apply` function with margin = 1 applies SomeOperation to each row of the matrix
print(apply(M, 1, SomeOperation))

