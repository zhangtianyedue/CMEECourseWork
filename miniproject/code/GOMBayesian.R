
library(rstan)
library(loo)       # compute WAIC
library(ggplot2)
library(bayesplot)
library(dplyr)

# load data
data <- read.csv("../results/Cleaned_LogisticGrowthData.csv")

# ensure the data structure
data <- data %>%
  filter(!is.na(PopBio) & !is.na(Time) & PopBio > 0)

# Pass data to the Stan model
stan_data <- list(
  N = nrow(data),
  Time = data$Time,
  PopBio = data$PopBio
)

# Run MCMC sampling.
fit_gompertz <- stan(
  file = "Gompertz_Bayes.stan",
  data = stan_data,
  iter = 2000, chains = 4, seed = 123,
  control = list(adapt_delta = 0.95, max_treedepth = 15)  # 增加稳定性
)

# save the model as rda file
saveRDS(fit_gompertz, "../results/bayes_gompertz.rds")

# Extract posterior samples.
posterior_samples <- as.matrix(fit_gompertz)

# ensure the correct colume name 
print(colnames(posterior_samples))

# Plot the posterior distribution of parameters.
png("../results/Posterior_r_max.png", width = 1000, height = 800)
mcmc_areas(as.array(posterior_samples[, "r_max", drop = FALSE])) +
  ggtitle("Posterior Distribution of r_max")
dev.off()

png("../results/Posterior_K.png", width = 1000, height = 800)
mcmc_areas(as.array(posterior_samples[, "K", drop = FALSE])) +
  ggtitle("Posterior Distribution of K")
dev.off()

png("../results/Posterior_t_lag.png", width = 1000, height = 800)
mcmc_areas(as.array(posterior_samples[, "t_lag", drop = FALSE])) +
  ggtitle("Posterior Distribution of t_lag")
dev.off()

# Output MCMC sampling diagnostics
print(fit_gompertz)

#Compute WAIC
log_lik_gompertz <- extract_log_lik(fit_gompertz, parameter_name = "log_lik")
waic_gompertz <- waic(log_lik_gompertz)
print(waic_gompertz)

# save WAIC results
write.csv(waic_gompertz$estimates, "../results/WAIC_Gompertz.csv", row.names = TRUE)

# Compute Bayes Factor (BF)
bayes_factor_gompertz <- sum(log_lik_gompertz)
print(paste("Bayes Factor for Gompertz:", bayes_factor_gompertz))

# save BF results
write.csv(data.frame(BF_Gompertz = bayes_factor_gompertz), "../results/BF_Gompertz.csv", row.names = FALSE)
