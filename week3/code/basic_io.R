# A simple script to illustrate R input-output.  
# Run line by line and check inputs and outputs to understand what is happening  

# Read the data from a CSV file located in the '../data/' directory.
# The 'header = TRUE' option indicates that the first row contains column names.
MyData <- read.csv("../data/trees.csv", header = TRUE) 

# Write the entire data frame to a new CSV file in the '../results/' directory.
write.csv(MyData, "../results/MyData.csv") 

# Append the first row of the data frame to the same file.
# The 'append = TRUE' option ensures that this new data is added to the existing file instead of overwriting it.
write.table(MyData[1,], file = "../results/MyData.csv", append = TRUE) 

# Write the entire data frame to the CSV file again, including row names.
# The 'row.names = TRUE' option includes the row indices from the data frame as the first column in the output.
write.csv(MyData, "../results/MyData.csv", row.names = TRUE) 

# Write the entire data frame to the CSV file but exclude column names.
# The 'col.names = FALSE' option ensures that the column headers are not written to the file.
write.table(MyData, "../results/MyData.csv", col.names = FALSE) 


