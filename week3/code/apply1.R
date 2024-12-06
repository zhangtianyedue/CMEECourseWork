## Build a random matrix
M <- matrix(rnorm(100), 10, 10) # Create a 10x10 matrix filled with random numbers from a normal distribution

## Take the mean of each row
RowMeans <- apply(M, 1, mean)  # Apply the 'mean' function to each row (margin = 1)
print(RowMeans)                # Print the calculated row means

## Now the variance
RowVars <- apply(M, 1, var)    # Apply the 'var' function to each row (margin = 1) to compute row variances
print(RowVars)                 # Print the calculated row variances

## By column
ColMeans <- apply(M, 2, mean)  # Apply the 'mean' function to each column (margin = 2)
print(ColMeans)                # Print the calculated column means

