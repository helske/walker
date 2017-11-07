functions {
  #include "common_functions.stan"
}

data {
  int<lower=0> k_fixed;
  int<lower=0> k_rw1;
  int<lower=0> k_rw2;
  int<lower=0> m;
  int<lower=0> k;
  int<lower=0> n;
  matrix[n, k_fixed] xreg_fixed;
  matrix[k, n] xreg_rw;
  vector[n] y;
  real<lower=0> sigma_y_mean;
  real<lower=0> sigma_y_sd;
  
  real beta_fixed_mean;
  real beta_rw1_mean;
  real beta_rw2_mean;
  real<lower=0> beta_fixed_sd;
  real<lower=0> beta_rw1_sd;
  real<lower=0> beta_rw2_sd;
  
  real sigma_rw1_mean;
  real sigma_rw2_mean;
  real<lower=0> sigma_rw1_sd;
  real<lower=0> sigma_rw2_sd;
  
  real<lower=0> slope_mean;
  real<lower=0> slope_sd;
  real<lower=0> gamma[n];
}

transformed data {
  vector[m] a1;
  matrix[m, m] P1 = rep_matrix(0.0, m, m);
  matrix[m, m] Tt = diag_matrix(rep_vector(1.0, m));
  
  Tt[(k_rw1+1):k, (k+1):m] = diag_matrix(rep_vector(1.0, k_rw2));
  
  for(i in 1:k_rw1) {
    a1[i] = beta_rw1_mean;
    P1[i, i] = beta_rw1_sd^2;
  }
   for(i in (k_rw1 + 1):k) {
    a1[i] = beta_rw2_mean;
    P1[i, i] = beta_rw2_sd^2;
  }
   for(i in (k + 1):m) {
    a1[i] = slope_mean;
    P1[i, i] = slope_sd^2;
  }

}

parameters {
  vector[k_fixed] beta_fixed;
  real<lower=0> sigma_rw1[k_rw1];
  real<lower=0> sigma_rw2[k_rw2];
  real<lower=0> sigma_y;
}

transformed parameters {
  matrix[m, m] Rt = rep_matrix(0.0, m, m);
  vector[n] xbeta;
  vector[n] y_;
  if (k_fixed > 0) {
    xbeta = xreg_fixed * beta_fixed;
  } else {
    xbeta = rep_vector(0.0, n);
  }
   y_ = y - xbeta;

  for(i in 1:k_rw1) {
    Rt[i, i] = sigma_rw1[i]^2;
  }
   for(i in 1:k_rw2) {
    Rt[k + i, k + i] = sigma_rw2[i]^2;
  }
}

model {
  sigma_y ~ normal(sigma_y_mean, sigma_y_sd);
  beta_fixed ~ normal(beta_fixed_mean, beta_fixed_sd);
  sigma_rw1 ~ normal(sigma_rw1_mean, sigma_rw1_sd);
  sigma_rw2 ~ normal(sigma_rw2_mean, sigma_rw2_sd);

  target += gaussian_filter(y_, a1, P1, sigma_y^2, Tt, Rt, xreg_rw, gamma);
}

generated quantities{

  vector[n] y_rep;
  matrix[k, n] beta_rw;
  matrix[k_rw2, n] slope;
  vector[n] y_fit;
  
  // sample coefficients given sigma's (no conditioning on y)

  for(i in 1:k_rw1) {
     beta_rw[i, 1] = normal_rng(beta_rw1_mean, beta_rw1_sd);
  }
  for(i in 1:k_rw2) {
     beta_rw[k_rw1 + i, 1] = normal_rng(beta_rw2_mean, beta_rw2_sd);
     slope[i, 1] = normal_rng(slope_mean, slope_sd);
  }

  for (t in 1:(n - 1)) {
    for(i in 1:k_rw1) {
      beta_rw[i, t + 1] = normal_rng(beta_rw[i, t], sigma_rw1[i]);
    }
    for(i in 1:k_rw2) {
      beta_rw[k_rw1 + i, t+1] = beta_rw[k_rw1 + i, t] + slope[i, t];
      slope[i, t + 1] = normal_rng(slope[i, t], sigma_rw2[i]);
    }
  }
  // sample new observations given previously simulated beta
  for(t in 1:n) {
    y_rep[t] = normal_rng(dot_product(xreg_rw[, t], beta_rw[, t]), sigma_y);
  }
  // perform mean correction to obtain sample from the posterior
  {
    matrix[m, n] states = gaussian_smoother(y_ - y_rep, a1, P1,
      sigma_y^2, Tt, Rt, xreg_rw, gamma);
    beta_rw = beta_rw + states[1:k, 1:n];
    slope = slope + states[(k + 1):m, 1:n];
  }

  // replicated data from posterior predictive distribution
  for(t in 1:n) {
    y_fit[t] = xbeta[t] + dot_product(xreg_rw[, t], beta_rw[, t]);
    y_rep[t] = normal_rng(y_fit[t], sigma_y);
  }

}
