# Loop through numbers from 1 to 10
for (i in 1:10) { 
  # Check if the current number `i` is even
  if ((i %% 2) == 0) 
    next # Skip the rest of the loop for this iteration and move to the next number
  
  # Print the current number if it is odd
  print(i) 
}
