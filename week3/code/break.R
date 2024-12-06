i <- 0 # Initialize i to 0 (starting value)

# Start an infinite while loop
while (i < Inf) { 
    # Check if i equals 10
    if (i == 10) { 
        break # Exit the while loop when i reaches 10
    } else { 
        # If i is not 10, print its current value
        cat("i equals ", i, "\n") 
        i <- i + 1 # Increment i by 1
    }
}
