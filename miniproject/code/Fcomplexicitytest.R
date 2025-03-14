# load data
library(dplyr)

# load the Optimized_Model_Comparison
results <- read.csv("../results/Optimized_Model_Comparison.csv")

# Filter out Unique_IDs that have RSS for **Gompertz, Logistic, and Richards models**.
valid_ids <- results %>%
  filter(Model %in% c("Gompertz", "Logistic", "Richards")) %>%
  group_by(Unique_ID) %>%
  summarise(
    has_Gompertz = any(Model == "Gompertz" & !is.na(RSS)),
    has_Logistic = any(Model == "Logistic" & !is.na(RSS)),
    has_Richards = any(Model == "Richards" & !is.na(RSS))
  ) %>%
  filter(has_Gompertz & has_Logistic & has_Richards) %>%  #Ensure that all three models have RSS
  pull(Unique_ID)  # "Extract only Unique_ID"

# "Filter the data, retaining only the Unique_IDs that meet the criteria."
filtered_results <- results %>%
  filter(Unique_ID %in% valid_ids)

# "Calculate the F-statistic."
compute_f_stat <- function(RSS1, df1, RSS2, df2) {
  return((RSS1 / df1) / (RSS2 / df2))
}

# "Calculate the F-statistic and p-value
stat_results <- filtered_results %>%
  group_by(Unique_ID) %>%
  summarise(
    RSS_Gompertz = mean(RSS[Model == "Gompertz"], na.rm = TRUE),
    RSS_Logistic = mean(RSS[Model == "Logistic"], na.rm = TRUE),
    RSS_Richards = mean(RSS[Model == "Richards"], na.rm = TRUE)
  ) %>%
  mutate(
    F_Gompertz_Richards = compute_f_stat(RSS_Gompertz, 3, RSS_Richards, 4),
    F_Logistic_Richards = compute_f_stat(RSS_Logistic, 3, RSS_Richards, 4),
    p_Gompertz_Richards = 1 - pf(F_Gompertz_Richards, df1 = 1, df2 = n() - 4),
    p_Logistic_Richards = 1 - pf(F_Logistic_Richards, df1 = 1, df2 = n() - 4)
  )

# Calculate descriptive statistics.
summary_stats <- stat_results %>%
  summarise(
    Mean_RSS_Gompertz = mean(RSS_Gompertz, na.rm = TRUE),
    Median_RSS_Gompertz = median(RSS_Gompertz, na.rm = TRUE),
    SD_RSS_Gompertz = sd(RSS_Gompertz, na.rm = TRUE),
    
    Mean_RSS_Logistic = mean(RSS_Logistic, na.rm = TRUE),
    Median_RSS_Logistic = median(RSS_Logistic, na.rm = TRUE),
    SD_RSS_Logistic = sd(RSS_Logistic, na.rm = TRUE),
    
    Mean_RSS_Richards = mean(RSS_Richards, na.rm = TRUE),
    Median_RSS_Richards = median(RSS_Richards, na.rm = TRUE),
    SD_RSS_Richards = sd(RSS_Richards, na.rm = TRUE),
    
    Mean_F_Gompertz = mean(F_Gompertz_Richards, na.rm = TRUE),
    Median_F_Gompertz = median(F_Gompertz_Richards, na.rm = TRUE),
    SD_F_Gompertz = sd(F_Gompertz_Richards, na.rm = TRUE),
    
    Mean_F_Logistic = mean(F_Logistic_Richards, na.rm = TRUE),
    Median_F_Logistic = median(F_Logistic_Richards, na.rm = TRUE),
    SD_F_Logistic = sd(F_Logistic_Richards, na.rm = TRUE)
  )

# **print statistical results**
print("==== Summary of Residual Sum of Squares (RSS) ====")
print(summary_stats[, c("Mean_RSS_Gompertz", "Median_RSS_Gompertz", "SD_RSS_Gompertz")])
print(summary_stats[, c("Mean_RSS_Logistic", "Median_RSS_Logistic", "SD_RSS_Logistic")])
print(summary_stats[, c("Mean_RSS_Richards", "Median_RSS_Richards", "SD_RSS_Richards")])

print("==== Summary of F Statistic (Gompertz vs. Richards) ====")
print(summary_stats[, c("Mean_F_Gompertz", "Median_F_Gompertz", "SD_F_Gompertz")])

print("==== Summary of F Statistic (Logistic vs. Richards) ====")
print(summary_stats[, c("Mean_F_Logistic", "Median_F_Logistic", "SD_F_Logistic")])

print("==== Individual Model Comparison Results ====")
print(stat_results)

# **save CSV**
write.csv(stat_results, "../results/Model_Statistical_Comparison.csv", row.names = FALSE)
