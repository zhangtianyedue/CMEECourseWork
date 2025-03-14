# Install required packages if not available
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos="https://cloud.r-project.org/")
}
library(dplyr)
library(rstan)
library(loo)
library(ggplot2)
library(bayesplot)
library(magrittr)  # Ensure pipeline operator is available

# Define results directory
results_dir <- "../results"

# Load data
data <- read.csv(file.path(results_dir, "Cleaned_LogisticGrowthData.csv"))

# Ensure the correct structure
data <- data %>%
  filter(!is.na(PopBio) & !is.na(Time) & PopBio > 0)

# Prepare data for Stan model
stan_data <- list(
  N = nrow(data),
  Time = data$Time,
  PopBio = data$PopBio
)

# **Load previously computed Gompertz model results**
gompertz_rds_path <- file.path(results_dir, "bayes_gompertz.rds")
if (file.exists(gompertz_rds_path)) {
  fit_gompertz <- readRDS(gompertz_rds_path)
} else {
  stop("Error: Gompertz model results not found! Expected file: ", gompertz_rds_path)
}

# **Run MCMC sampling for Logistic model**
fit_logistic <- stan(
  file = "Logistic_Bayes.stan",
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

# **Compute WAIC**
log_lik_logistic <- extract_log_lik(fit_logistic)
log_lik_gompertz <- extract_log_lik(fit_gompertz)

waic_logistic <- waic(log_lik_logistic)
waic_gompertz <- waic(log_lik_gompertz)

# **Compute Bayes Factor (BF)**
bf <- exp((waic_gompertz$estimates[3,1] - waic_logistic$estimates[3,1]) / 2)
print(paste("Bayes Factor (Gompertz vs Logistic):", bf))

# **Save Logistic model results**
logistic_rds_path <- file.path(results_dir, "bayes_model_logistic.rds")
saveRDS(fit_logistic, logistic_rds_path)

# **Extract posterior samples**
posterior_samples_logistic <- as.data.frame(as.matrix(fit_logistic))
posterior_samples_gompertz <- as.data.frame(as.matrix(fit_gompertz))

# **Save posterior samples**
posterior_samples_logistic_path <- file.path(results_dir, "posterior_samples_logistic.rds")
posterior_samples_gompertz_path <- file.path(results_dir, "posterior_samples_gompertz.rds")

saveRDS(posterior_samples_logistic, posterior_samples_logistic_path)
saveRDS(posterior_samples_gompertz, posterior_samples_gompertz_path)

# **Compute number of samples and time points**
n_samples <- nrow(posterior_samples_logistic)
n_time <- length(data$Time)

# **Initialize prediction matrices**
y_rep_logistic <- matrix(NA, nrow = n_samples, ncol = n_time)
y_rep_gompertz <- matrix(NA, nrow = n_samples, ncol = n_time)

print(dim(y_rep_logistic))
print(dim(y_rep_gompertz))

# **Generate predictions**
for (i in 1:min(n_samples, nrow(posterior_samples_gompertz))) {
  # Logistic prediction
  y_rep_logistic[i, ] <- posterior_samples_logistic[i, "K"] / 
    (1 + exp(-posterior_samples_logistic[i, "r_max"] * 
               (data$Time - posterior_samples_logistic[i, "N_0"])))
  
  # Gompertz prediction
  y_rep_gompertz[i, ] <- posterior_samples_gompertz[i, "K"] * exp(
    -exp(posterior_samples_gompertz[i, "r_max"] * exp(1) * 
           (posterior_samples_gompertz[i, "t_lag"] - data$Time) / 
           posterior_samples_gompertz[i, "K"] + 1))
}

# **Save prediction results**
saveRDS(y_rep_logistic, file.path(results_dir, "y_rep_logistic.rds"))
saveRDS(y_rep_gompertz, file.path(results_dir, "y_rep_gompertz.rds"))

# **Output WAIC and LOO results**
print(waic_logistic)
print(waic_gompertz)

loo_logistic <- loo(log_lik_logistic)
loo_gompertz <- loo(log_lik_gompertz)

print(loo_logistic)
print(loo_gompertz)

# **Extract 95% confidence intervals**
posterior_logistic_filtered <- posterior_samples_logistic[, c("r_max", "K", "N_0", "sigma")]
posterior_gompertz_filtered <- posterior_samples_gompertz[, c("r_max", "K", "t_lag", "sigma")]

logistic_summary <- summary(posterior_logistic_filtered, probs = c(0.025, 0.5, 0.975))
gompertz_summary <- summary(posterior_gompertz_filtered, probs = c(0.025, 0.5, 0.975))

print(logistic_summary)
print(gompertz_summary)

# **Plot posterior distributions**
mcmc_intervals(posterior_logistic_filtered, pars = c("r_max", "K", "N_0", "sigma"))
ggsave(file.path(results_dir, "logistic_posterior.png"), width=8, height=6)

mcmc_intervals(posterior_gompertz_filtered, pars = c("r_max", "K", "t_lag", "sigma"))
ggsave(file.path(results_dir, "gompertz_posterior.png"), width=8, height=6)

# **Save final results**
saveRDS(logistic_summary, file.path(results_dir, "logistic_summary.rds"))
saveRDS(gompertz_summary, file.path(results_dir, "gompertz_summary.rds"))

print("Logistic Bayesian analysis completed successfully!")
