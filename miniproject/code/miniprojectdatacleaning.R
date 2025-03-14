# Load necessary libraries
library(readr)
library(dplyr)

# Load the dataset
df <- read.csv('../data/LogisticGrowthData.csv')

# Step 1: Remove NA values
data <- na.omit(df)

# Step 2: Ensure PopBio is numeric and remove only negative values
data$PopBio <- as.numeric(data$PopBio)  # Convert to numeric
data <- data[!is.na(data$PopBio) & data$PopBio >= 0, ]  # Remove negative values only

# Step 3: Apply log transformation without shift to preserve lag phase
data$LogN <- log(data$PopBio)  # Allow negative LogN values

# Step 4: Remove invalid values
data <- data[!is.na(data$LogN) & !is.nan(data$LogN) & !is.infinite(data$LogN), ]

# Step 5: Ensure Temp is always an integer
data$Temp <- as.integer(data$Temp)  # 强制转换为整数
data <- data[data$Temp > 0, ]  # 确保有效的温度值

# Check if Temp is truly an integer
print(unique(data$Temp))  # 确保这里不会出现 1.5

# Step 6: Create Unique_ID
data <- data %>%
  mutate(Unique_ID = paste(Species, Temp, Medium, Citation, sep = "_"))  # 直接用整数 Temp

# Step 7: Remove rows where Time is non-positive or too small
data <- data[data$Time >= 1, ]

# Step 8: Split data into subsets based on Unique_ID
subsets <- split(data, data$Unique_ID)

# Save the cleaned data
write.csv(data, file = '../results/Cleaned_LogisticGrowthData.csv', row.names = FALSE)

# # Check if Unique_ID is correct
print(unique(data$Unique_ID))   # Ensure Unique_ID structure is correct
# Count the number of unique IDs
num_unique_ids <- length(unique(data$Unique_ID))

# Print the result
print(paste("Total number of unique IDs:", num_unique_ids))
