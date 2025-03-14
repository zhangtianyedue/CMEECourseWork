data {
  int<lower=1> N;         // 数据点数量
  vector[N] Time;         // 时间
  vector[N] PopBio;       // 观测值（种群生物量）
}

parameters {
  real<lower=0> r_max;    // 最大生长速率
  real<lower=0> K;        // 生态位容量
  real<lower=0> N_0;      // 初始种群
  real<lower=0> sigma;    // 观测误差
}

model {
  vector[N] mu;

  // 先验分布（与 Gompertz 保持一致）
  r_max ~ lognormal(0, 1);
  K ~ lognormal(3, 1);
  N_0 ~ lognormal(1, 1);
  sigma ~ cauchy(0, 2);

  // **Logistic 模型**
  for (i in 1:N) {
    mu[i] = K / (1 + ((K - N_0) / N_0) * exp(-r_max * Time[i]));
  }

  // 观测数据 Log-Normal 似然
  PopBio ~ lognormal(log(mu), sigma);
}

generated quantities {
  vector[N] mu_pred;   // 重新计算 mu
  vector[N] log_lik;   // 对数似然

  for (i in 1:N) {
    mu_pred[i] = K / (1 + ((K - N_0) / N_0) * exp(-r_max * Time[i]));  // 重新计算 mu
    log_lik[i] = lognormal_lpdf(PopBio[i] | log(mu_pred[i]), sigma);   // 计算 log_lik
  }
}
