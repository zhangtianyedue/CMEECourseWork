################################################################
################## Wrangling the Pound Hill Dataset ############
################################################################

############# Load the dataset ###############
require(reshape2) # Load the reshape2 package for reshaping data
install.packages("tidyverse") # Install the tidyverse package for data wrangling
library(tidyverse) # Load the tidyverse package
library(tidyr) # Load tidyr for reshaping data
library(dplyr) # Load dplyr for data manipulation

# Set the working directory to the correct path
getwd() # Display the current working directory
setwd("/Users/tianyezhang/Desktop/CMEECourseWork/week4/code") # Change the working directory

# Read the raw data without headers, converting it into a matrix
MyData <- read_csv("../data/PoundHillData.csv", col_names = FALSE) %>% as.matrix() # Use pipe operator

# Read the metadata with headers
MyMetaData <- read_csv2("../data/PoundHillMetaData.csv") 

############# Inspect the dataset ###############
# Inspect the structure of MyData using glimpse (a better version of str())
glimpse(MyData)
# Check the dimensions of MyData (number of rows and columns)
dim(MyData)
# Display the first few rows of MyData
head(MyData)
# Open MyData and MyMetaData in an interactive viewer
fix(MyData)
fix(MyMetaData)

############# Transpose ###############
# Transpose the dataset to move species into columns and treatments into rows
MyData <- t(MyData)
# Check the first few rows of the transposed dataset
head(MyData)
# Check the new dimensions of MyData
dim(MyData)
# Display the column names of MyData
colnames(MyData)

############# Replace species absences with zeros ###############
# Replace empty values with zeros across the entire dataset using mutate_all
MyData <- as_tibble(MyData) %>%
  mutate_all(~ replace(., . == "", 0)) # Apply the replacement to all columns

############# Convert raw matrix to data frame ###############
# Remove the first row and set it as column names for the new data frame
TempData <- MyData[-1, ] %>% 
  as_tibble(.name_repair = "minimal") %>% # Convert to a tibble with minimal name repair
  set_names(MyData[1, ]) # Assign column names from the first row

# Check the first few rows of the converted data frame
head(TempData)

############# Reshape data using tidyr::pivot_longer() ###############
# Reshape the data from wide to long format
MyWrangledData <- TempData %>%
  pivot_longer(
    cols = -c(Cultivation, Block, Plot, Quadrat), # Exclude these columns
    names_to = "Species",                         # New column for species names
    values_to = "Count"                           # New column for counts
  )

# Convert specific columns to factors and Count column to integers
MyWrangledData <- MyWrangledData %>%
  mutate(
    Cultivation = as.factor(Cultivation), # Convert Cultivation to factor
    Block = as.factor(Block),             # Convert Block to factor
    Plot = as.factor(Plot),               # Convert Plot to factor
    Quadrat = as.factor(Quadrat),         # Convert Quadrat to factor
    Count = as.integer(Count)             # Convert Count to integer
  )

# Display the frequency table of the Cultivation column
table(MyWrangledData$Cultivation)

# Inspect the structure of the wrangled dataset
glimpse(MyWrangledData)

# Display the first few rows of the wrangled dataset
head(MyWrangledData)

# Check the dimensions of the wrangled dataset
dim(MyWrangledData)

############# Exploring the data ###############
# List all tidyverse packages including "tidyverse" itself
tidyverse_packages(include_self = TRUE)

# Convert MyWrangledData to a tibble for better data handling
MyWrangledData <- dplyr::as_tibble(MyWrangledData)

# Display the tibble
MyWrangledData

# Use glimpse to inspect the structure of the tibble
glimpse(MyWrangledData)

# Filter rows where Count is greater than 100
filter(MyWrangledData, Count > 100)

# Extract rows 10 to 15 from the dataset
slice(MyWrangledData, 10:15)

# Group by Species and calculate the average Count for each species
MyWrangledData %>%
  group_by(Species) %>%
  summarise(avg = mean(Count))

# Perform the same aggregation using the base R aggregate function
aggregate(MyWrangledData$Count, list(MyWrangledData$Species), FUN = mean)
