################################################################
################## Wrangling the Pound Hill Dataset ############
################################################################

############# Load the dataset ###############
# Load the raw data as a matrix (header = FALSE because the dataset has no real headers)
MyData <- as.matrix(read.csv("../data/PoundHillData.csv", header = FALSE))
head(MyData) # Display the first few rows of the data

# Load the metadata with proper headers (header = TRUE because metadata has headers)
MyMetaData <- read.csv("../data/PoundHillMetaData.csv", header = TRUE, sep = ";")

############# Inspect the dataset ###############
# Check the first few rows of the dataset
head(MyData)

# Get the dimensions of the dataset (number of rows and columns)
dim(MyData)

# Examine the structure of the dataset (e.g., data types, dimensions)
str(MyData)

# Open an interactive viewer to inspect and edit the data (if necessary)
fix(MyData)
fix(MyMetaData)

############# Transpose ###############
# Transpose the dataset to switch rows and columns (species into columns, treatments into rows)
MyData <- t(MyData)
head(MyData) # Display the first few rows after transposing
dim(MyData) # Check the new dimensions of the dataset
colnames(MyData) # Check the column names (which are still missing)

############# Replace species absences with zeros ###############
# Replace empty values ("") in the dataset with zeros to indicate absence
MyData[MyData == ""] <- 0

############# Convert raw matrix to data frame ###############
# Convert the matrix (excluding the first row) into a data frame
# `stringsAsFactors = FALSE` ensures that strings remain as characters, not factors
TempData <- as.data.frame(MyData[-1,], stringsAsFactors = FALSE) 

# Assign column names to the data frame using the first row of the original matrix
colnames(TempData) <- MyData[1,]

head(TempData) # Display the first few rows of the new data frame

############# Convert from wide to long format ###############
# Load the `reshape2` package for reshaping the data
require(reshape2)

# Use the `melt` function to reshape the data from wide to long format
# `id` specifies the identifier columns, and `variable.name` and `value.name` name the new columns
MyWrangledData <- melt(TempData, id = c("Cultivation", "Block", "Plot", "Quadrat"), 
                       variable.name = "Species", value.name = "Count")

# Convert specific columns to factors
MyWrangledData[, "Cultivation"] <- as.factor(MyWrangledData[, "Cultivation"])
MyWrangledData[, "Block"] <- as.factor(MyWrangledData[, "Block"])
MyWrangledData[, "Plot"] <- as.factor(MyWrangledData[, "Plot"])
MyWrangledData[, "Quadrat"] <- as.factor(MyWrangledData[, "Quadrat"])

# Convert the "Count" column to integer
MyWrangledData[, "Count"] <- as.integer(MyWrangledData[, "Count"])

# Display the frequency of each level in the "Cultivation" column
table(MyWrangledData$Cultivation)

# Inspect the structure, first few rows, and dimensions of the wrangled data
str(MyWrangledData)
head(MyWrangledData)
dim(MyWrangledData)

############# Exploring the data (extend the script below) ###############
# Install the `tidyverse` package for data manipulation
install.packages("tidyverse")
library(tidyverse)

# List all `tidyverse` packages, including "tidyverse" itself
tidyverse_packages(include_self = TRUE)

# Convert the wrangled data to a tibble (a tidyverse-friendly data frame)
MyWrangledData <- dplyr::as_tibble(MyWrangledData)
MyWrangledData # Display the tibble

# Glimpse the structure of the tibble (similar to `str()`, but formatted better)
glimpse(MyWrangledData)

# Filter rows where "Count" is greater than 100 (similar to `subset()`)
filter(MyWrangledData, Count > 100)

# Select a specific range of rows (e.g., rows 10 to 15)
slice(MyWrangledData, 10:15)

# Summarize the average count for each species using `group_by` and `summarise`
MyWrangledData %>%
  group_by(Species) %>%
  summarise(avg = mean(Count))

# Perform the same summarization using the `aggregate` function
aggregate(MyWrangledData$Count, list(MyWrangledData$Species), FUN = mean)
