
library(ggplot2)
library(dplyr)
library(minpack.lm)
library(gridExtra)
library(grid)
library(tibble)

# read data
data <- read.csv("../results/Cleaned_LogisticGrowthData.csv")
data$Unique_ID <- as.factor(data$Unique_ID)

# Compute RSS
compute_rss <- function(model, data) {
  predicted_values <- tryCatch(predict(model, newdata = data), error = function(e) rep(NA, nrow(data)))
  residuals_vec <- data$PopBio - predicted_values
  if (is.numeric(residuals_vec) && !all(is.na(residuals_vec))) {
    return(sum(residuals_vec^2))
  }
  return(NA)
}

# Compute AICc
compute_aicc <- function(aic, n, k) {
  return(aic + (2 * k * (k + 1)) / (n - k - 1))
}

# define the linear model 
fit_linear_models <- function(subset) {
  models <- list()
  models$linear <- lm(PopBio ~ Time, data = subset)
  models$quadratic <- lm(PopBio ~ poly(Time, 2, raw = TRUE), data = subset)
  return(models)
}

# define the non_linear model （log version）
gompertz_model <- function(t, r_max, logK, logN_0, t_lag){
  K <- exp(logK)
  N_0 <- exp(logN_0)
  return(N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t) / ((K - N_0) * log(10)) + 1)))
}

logistic_model <- function(t, r_max, logK, logN_0){
  K <- exp(logK)
  N_0 <- exp(logN_0)
  return(N_0 * K * exp(r_max * t) / (K + N_0 * (exp(r_max * t) - 1)))
}

richards_model <- function(t, r_max, logK, logN_0, v) {
  K <- exp(logK)
  N_0 <- exp(logN_0)
  return(K / ((1 + ((K/N_0 - 1) * exp(-r_max * t)))^(1/v)))
}

# time segment（according to the median）
split_data_by_time <- function(subset) {
  median_time <- median(subset$Time, na.rm = TRUE)
  subset_1 <- subset %>% filter(Time <= median_time)
  subset_2 <- subset %>% filter(Time > median_time)
  return(list(subset_1, subset_2))
}

# Evaluate model 
evaluate_models <- function(unique_id, subset) {
  model_results <- data.frame(Unique_ID = character(), Model = character(), AIC = numeric(), BIC = numeric(), AICc = numeric(), RSS = numeric(), stringsAsFactors = FALSE)
  
  # linear model 
  linear_models <- fit_linear_models(subset)
  for (name in names(linear_models)) {
    model <- linear_models[[name]]
    k <- length(coef(model))
    n <- nrow(subset)
    aic_val <- AIC(model)
    aicc_val <- compute_aicc(aic_val, n, k)
    model_results <- rbind(model_results, data.frame(
      Unique_ID = unique_id, Model = name, AIC = aic_val, BIC = BIC(model), AICc = aicc_val, RSS = compute_rss(model, subset)
    ))
  }
  
  # Non-linear model ，seperate 2 parts
  subsets <- split_data_by_time(subset)
  
  for (sub in subsets) {
    if (nrow(sub) < 3) next  # Avoid fitting failures due to insufficient data.
    
    # compute the initial condition 
    N_0_start <- min(sub$PopBio, na.rm = TRUE) * 0.9
    K_start <- max(sub$PopBio, na.rm = TRUE) * 1.1
    t_lag_start <- median(sub$Time, na.rm = TRUE)
    r_max_start <- 0.1  # set an initial condition 
    v_start <- runif(1, 0.5, 2)
    
    # Avoid NaN or Inf
    if (!is.finite(N_0_start) || N_0_start <= 0) N_0_start <- 1
    if (!is.finite(K_start) || K_start <= 0) K_start <- max(sub$PopBio, na.rm = TRUE)
    if (!is.finite(t_lag_start) || t_lag_start < 0) t_lag_start <- min(sub$Time, na.rm = TRUE)
    
    # fitting non_linear model
    tryCatch({
      fit_gompertz <- nlsLM(PopBio ~ gompertz_model(Time, r_max, logK, logN_0, t_lag), data = sub,
                            start = list(r_max = r_max_start, logK = log(K_start), logN_0 = log(N_0_start), t_lag = t_lag_start),
                            algorithm = "port")
      k <- length(coef(fit_gompertz))
      n <- nrow(sub)
      aic_val <- AIC(fit_gompertz)
      aicc_val <- compute_aicc(aic_val, n, k)
      model_results <- rbind(model_results, data.frame(
        Unique_ID = unique_id, Model = "Gompertz", AIC = aic_val, BIC = BIC(fit_gompertz), AICc = aicc_val, RSS = compute_rss(fit_gompertz, sub)
      ))
    }, error = function(e) {})
    
    tryCatch({
      fit_logistic <- nlsLM(PopBio ~ logistic_model(Time, r_max, logK, logN_0), data = sub,
                            start = list(r_max = r_max_start, logK = log(K_start), logN_0 = log(N_0_start)),
                            algorithm = "port")
      k <- length(coef(fit_logistic))
      n <- nrow(sub)
      aic_val <- AIC(fit_logistic)
      aicc_val <- compute_aicc(aic_val, n, k)
      model_results <- rbind(model_results, data.frame(
        Unique_ID = unique_id, Model = "Logistic", AIC = aic_val, BIC = BIC(fit_logistic), AICc = aicc_val, RSS = compute_rss(fit_logistic, sub)
      ))
    }, error = function(e) {})
    
    tryCatch({
      fit_richards <- nlsLM(PopBio ~ richards_model(Time, r_max, logK, logN_0, v), data = sub,
                            start = list(r_max = r_max_start, logK = log(K_start), logN_0 = log(N_0_start), v = v_start),
                            algorithm = "port")
      k <- length(coef(fit_richards))
      n <- nrow(sub)
      aic_val <- AIC(fit_richards)
      aicc_val <- compute_aicc(aic_val, n, k)
      model_results <- rbind(model_results, data.frame(
        Unique_ID = unique_id, Model = "Richards", AIC = aic_val, BIC = BIC(fit_richards), AICc = aicc_val, RSS = compute_rss(fit_richards, sub)
      ))
    }, error = function(e) {})
  }
  
  return(model_results)
}

# iterate every Unique_ID
results <- do.call(rbind, lapply(unique(data$Unique_ID), function(id) evaluate_models(id, subset(data, Unique_ID == id))))
write.csv(results, "../results/Optimized_Model_Comparison.csv", row.names = FALSE)


# Compute the proportion of each model being the best
best_models <- results %>%
  group_by(Unique_ID) %>%
  filter(AICc == min(AICc)) %>%
  ungroup()

# Compute the proportion of the best model 
model_count <- best_models %>%
  group_by(Model) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)
print(model_count)
print(model_count)

# Visualize the proportion of model selection.
p1 <- ggplot(model_count, aes(x = Model, y = Percentage, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Best Model Proportion (Based on AICc)",
       x = "Model",
       y = "Percentage (%)") +
  theme_minimal()

ggsave("../results/Best_Model_Proportion_Percentage2.png", plot = p1, width = 6, height = 4)
print(p1)
