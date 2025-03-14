library(ggplot2)
library(dplyr)
library(minpack.lm)
library(gridExtra)
library(grid)
library(tibble)

# Load data
data <- read.csv("../results/Cleaned_LogisticGrowthData.csv")
data$Unique_ID <- as.factor(data$Unique_ID)

# Store final results
results <- data.frame(Unique_ID = character(), Model = character(), AIC = numeric(), BIC = numeric(), AICc = numeric(), RSS = numeric(), stringsAsFactors = FALSE)

# Compute AICc
compute_aicc <- function(model, data) {
  k <- length(coef(model))  # Number of parameters
  n <- nrow(data) # Sample size
  aic_val <- AIC(model)
  return(aic_val + (2 * k * (k + 1)) / (n - k - 1))
}

# Compute RSS
compute_rss <- function(model, data) {
  predicted_values <- tryCatch(predict(model, newdata = data), error = function(e) rep(NA, nrow(data)))
  residuals_vec <- data$PopBio - predicted_values
  if (is.numeric(residuals_vec) && !all(is.na(residuals_vec))) {
    return(sum(residuals_vec^2))
  }
  return(NA)
}

# Compute maximum slope
find_max_slope <- function(subset, window_size = 3) {
  max_slope <- -Inf
  if (window_size > nrow(subset)) return(0.1)
  for (i in 1:(nrow(subset) - window_size + 1)) {
    fit <- lm(PopBio ~ Time, data = subset[i:(i + window_size - 1), ])
    slope <- coef(fit)["Time"]
    if (!is.na(slope) && slope > max_slope) max_slope <- slope
  }
  return(ifelse(is.finite(max_slope) && max_slope > 0, max_slope, 0.1))
}

# Fit linear models
fit_linear_models <- function(subset) {
  models <- list()
  models$linear <- lm(PopBio ~ Time, data = subset)
  models$quadratic <- lm(PopBio ~ poly(Time, 2, raw = TRUE), data = subset)
  return(models)
}

# Define nonlinear models
gompertz_model <- function(t, log_r_max, log_K, log_N_0, log_t_lag) {
  r_max <- exp(log_r_max)
  K <- exp(log_K)
  N_0 <- exp(log_N_0)
  t_lag <- exp(log_t_lag)
  return(N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t) / ((K - N_0) * log(10)) + 1)))
}

logistic_model <- function(t, log_r_max, log_K, log_N_0) {
  r_max <- exp(log_r_max)
  K <- exp(log_K)
  N_0 <- exp(log_N_0)
  return(N_0 * K * exp(r_max * t) / (K + N_0 * (exp(r_max * t) - 1)))
}

richards_model <- function(t, log_r_max, log_K, log_N_0, log_v) {
  r_max <- exp(log_r_max)
  K <- exp(log_K)
  N_0 <- exp(log_N_0)
  v <- exp(log_v)
  return(K / ((1 + ((K/N_0 - 1) * exp(-r_max * t)))^(1/v)))
}

# Evaluate models
evaluate_models <- function(unique_id, subset) {
  model_results <- data.frame(Unique_ID = character(), Model = character(), AIC = numeric(), BIC = numeric(), AICc = numeric(), RSS = numeric(), stringsAsFactors = FALSE)
  
  # Fit linear models
  linear_models <- fit_linear_models(subset)
  for (name in names(linear_models)) {
    model <- linear_models[[name]]
    model_results <- rbind(model_results, data.frame(
      Unique_ID = unique_id, Model = name, AIC = AIC(model), BIC = BIC(model), AICc = compute_aicc(model, subset), RSS = compute_rss(model, subset)
    ))
  }
  
  # Compute initial values
  N_0_start <- log(min(subset$PopBio, na.rm = TRUE) * 0.9)
  K_start <- log(max(subset$PopBio, na.rm = TRUE) * 1.1)
  t_lag_start <- log(median(subset$Time, na.rm = TRUE))
  r_max_start <- log(find_max_slope(subset))
  v_start <- runif(1, 0.5, 2)
  
  # Avoid NaN or Inf
  if (!is.finite(N_0_start)) N_0_start <- log(1)
  if (!is.finite(K_start)) K_start <- log(max(subset$PopBio, na.rm = TRUE))
  if (!is.finite(t_lag_start)) t_lag_start <- log(min(subset$Time, na.rm = TRUE))
  if (!is.finite(r_max_start)) r_max_start <- log(0.1)
  
  # Fit nonlinear models
  tryCatch({
    fit_gompertz <- nlsLM(PopBio ~ gompertz_model(Time, log_r_max, log_K, log_N_0, log_t_lag), data = subset,
                          start = list(log_r_max = r_max_start, log_N_0 = N_0_start, log_K = K_start, log_t_lag = t_lag_start),
                          algorithm = "port", trace = TRUE)
    model_results <- rbind(model_results, data.frame(
      Unique_ID = unique_id, Model = "Gompertz", AIC = AIC(fit_gompertz), BIC = BIC(fit_gompertz), AICc = compute_aicc(fit_gompertz, subset), RSS = compute_rss(fit_gompertz, subset)
    ))
  }, error = function(e) { message("Gompertz failed for ", unique_id, ": ", e$message) })
  
  tryCatch({
    fit_logistic <- nlsLM(PopBio ~ logistic_model(Time, log_r_max, log_K, log_N_0), data = subset,
                          start = list(log_r_max = r_max_start, log_N_0 = N_0_start, log_K = K_start),
                          algorithm = "port", trace = TRUE)
    model_results <- rbind(model_results, data.frame(
      Unique_ID = unique_id, Model = "Logistic", AIC = AIC(fit_logistic), BIC = BIC(fit_logistic), AICc = compute_aicc(fit_logistic, subset), RSS = compute_rss(fit_logistic, subset)
    ))
  }, error = function(e) { message("Logistic failed for ", unique_id, ": ", e$message) })
  
  return(model_results)
}

# Iterate through the dataset
for (unique_id in unique(data$Unique_ID)) {
  subset <- subset(data, Unique_ID == unique_id)
  results <- rbind(results, evaluate_models(unique_id, subset))
}
# Compute the proportion of each model being the best
best_models <- results %>%
  group_by(Unique_ID) %>%
  filter(AIC == min(AIC)) %>%
  ungroup()

# Add Model_Type class
best_models <- best_models %>%
  mutate(Model_Type = case_when(
    Model %in% c("linear", "quadratic") ~ "Linear",
    Model %in% c("Gompertz", "Logistic", "Richards") ~ "Nonlinear"
  ))

# Compute the best model proportion 
model_count <- best_models %>%
  group_by(Model) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)
print(model_count)
# compute the proportion of linear VS non-linear model 
model_type_count <- best_models %>%
  group_by(Model_Type) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

print(model_type_count)
p1 <- ggplot(model_count, aes(x = Model, y = Percentage, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Best Model Proportion (Percentage)",
       x = "Model",
       y = "Percentage (%)") +
  theme_minimal()

ggsave("../results/Best_Model_Proportion_Percentage.png", plot = p1, width = 6, height = 4)
# save the results
write.csv(results, "../results/Model_Fitting_Results_LogTransformed.csv", row.names = FALSE)
