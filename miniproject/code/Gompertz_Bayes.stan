data {
  int<lower=1> N;         // 数据点数量
  vector[N] Time;         // 时间
  vector[N] PopBio;       // 观测值（种群生物量）
}

parameters {
  real<lower=0> r_max;    // 最大生长速率
  real<lower=0> K;        // 生态位容量
  real<lower=0> t_lag;    // 滞后期时间
  real<lower=0> sigma;    // 观测误差
}

model {
  vector[N] mu;           // 预测值

  // 先验分布（Log-Normal）
  r_max ~ lognormal(0, 1);
  K ~ lognormal(3, 1);      // K 期望值较大
  t_lag ~ lognormal(1, 1);
  sigma ~ cauchy(0, 2);     // 误差项先验

  // Gompertz 模型
  for (i in 1:N) {
    mu[i] = K * exp(-exp(r_max * exp(1) * (t_lag - Time[i]) / K + 1));
  }

  // 观测数据的对数正态似然函数
  PopBio ~ lognormal(log(mu), sigma);
}

generated quantities {
  vector[N] log_lik;  // 计算对数似然 log_lik 用于 WAIC

  for (i in 1:N) {
    log_lik[i] = lognormal_lpdf(PopBio[i] | log(K * exp(-exp(r_max * exp(1) * (t_lag - Time[i]) / K + 1))), sigma);
  }
}
