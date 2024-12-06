a <- TRUE
# Check if the variable `a` is TRUE or FALSE
if (a == TRUE) { 
    print("a is TRUE") # Print a message if `a` is TRUE
} else {
    print("a is FALSE") # Print a message if `a` is FALSE
}

z <- runif(1) # Generate a single random number from a uniform distribution between 0 and 1
# Check if the random number `z` is less than or equal to 0.5
if (z <= 0.5) { 
    print("Less than a half") # Print a message if `z` is less than or equal to 0.5
}

# Loop through integers 1 to 10
for (i in 1:10) { 
  j <- i * i # Calculate the square of the current number
  print(paste(i, "squared is", j)) # Print the number and its square
}

# Loop through a vector of species names
for (species in c('Heliodoxa rubinoides', 
                  'Boissonneaua jardini', 
                  'Sula nebouxii')) {
  print(paste('The species is', species)) # Print the current species name
}

# Loop through a character vector `v1`
v1 <- c("a", "bc", "def") # Create a vector of strings
for (i in v1) { 
  print(i) # Print each string in the vector
}

# Initialize a variable `i` to 0
i <- 0 
# Loop until `i` reaches 10
while (i < 10) { 
  i <- i + 1 # Increment `i` by 1
  print(i^2) # Print the square of the current value of `i`
}
