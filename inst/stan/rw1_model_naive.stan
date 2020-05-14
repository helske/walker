data {
  int<lower=0> k;
  int<lower=0> n;
  matrix[k, n] xreg;
  vector[n] y;
  vector[k] beta_mean;
  vector[k] beta_sd;
  vector[k + 1] sigma_mean;
  vector[k + 1] sigma_sd;

}
transformed data {
  vector[k] sigma_b_mean = sigma_mean[2:];
  vector[k] sigma_b_sd = sigma_sd[2:];
  real sigma_y_mean = sigma_mean[1];
  real sigma_y_sd = sigma_sd[1];
}

parameters {
  real<lower=0> sigma_b[k];
  real<lower=0> sigma_y;
  matrix[k, n] beta_raw;
}

transformed parameters {
  matrix[k, n] beta;
  
  vector[k] tmp;
  beta[, 1] = beta_mean + beta_sd .* beta_raw[, 1]; 
  for(t in 2:n) {
    tmp = beta[, t - 1];
    beta[, t] = tmp + to_vector(sigma_b) .* beta_raw[, t];  
  }
}

model {
  sigma_b ~ normal(sigma_b_mean, sigma_b_sd);
  sigma_y ~ normal(sigma_y_mean, sigma_y_sd);
  to_vector(beta_raw) ~ normal(0, 1);
  {
    row_vector[n] mu = columns_dot_product(xreg, beta);
    y ~ normal(mu, sigma_y);
  }
}
