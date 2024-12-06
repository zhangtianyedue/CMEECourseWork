# Function without preallocation
NoPreallocFun <- function(x) {
    # Initialize an empty vector
    a <- vector() # Creates an empty vector
    
    # Loop through numbers from 1 to x
    for (i in 1:x) {
        a <- c(a, i) # Concatenate the current number to the vector `a`
        print(a) # Print the current state of the vector
        print(object.size(a)) # Print the memory size of the vector
    }
}

# Measure the execution time of the NoPreallocFun function for input 10
system.time(NoPreallocFun(10)) 

# Function with preallocation
PreallocFun <- function(x) {
    # Initialize a preallocated vector with `NA` values
    a <- rep(NA, x) # Create a vector of length `x` with all elements set to `NA`
    
    # Loop through numbers from 1 to x
    for (i in 1:x) {
        a[i] <- i # Assign the current number to the ith position of the vector `a`
        print(a) # Print the current state of the vector
        print(object.size(a)) # Print the memory size of the vector
    }
}

# Measure the execution time of the PreallocFun function for input 10
system.time(PreallocFun(10))
