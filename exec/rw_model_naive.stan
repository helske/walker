data {
  int<lower=0> k;
  int<lower=0> n;
  matrix[n, k] xreg;
  vector[n] y;
  vector[k] beta_mean;
  vector[k] beta_sd;
  vector[k + 1] sigma_mean;
  vector[k + 1] sigma_sd;

}
parameters {
  real<lower=0> sigma[1 + k];
  matrix[k, n] beta_raw;
}
transformed parameters {
  matrix[k, n] beta;
  vector[k] sigma_vec;
  
  vector[k] tmp;
  for(i in 1:k) sigma_vec[i] = sigma[1 + i];
  
  beta[, 1] = beta_mean + beta_sd .* beta_raw[, 1]; 
  for(t in 2:n) {
    tmp = beta[, t - 1];
    beta[, t] = tmp + sigma_vec .* beta_raw[, t];  
  }
}
model {
  target += normal_lpdf(sigma | sigma_mean, sigma_sd);
 
  for(t in 1:n) {
    target += normal_lpdf(beta_raw[, t] | 0, 1);
    target += normal_lpdf(y[t] | xreg[t, ] * beta[, t], sigma[1]);    
  }

}
