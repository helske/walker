functions {
// Functions for Kalman filter and smoother for dynamic regression model
// note that these functions are not fully optimised yet

// univariate Kalman filter for RW2 model, returns the log-likelihood
real gaussian_filter_rw2_lpdf(vector y, vector a1, vector P1, real Ht, matrix Rt, matrix xreg) {

  int n = rows(y);
  int k = rows(xreg);
  int m = 2 * k;
  
  real loglik = 0.0;
  
  vector[m] x = rep_vector(0, m);
  matrix[m, m] P = diag_matrix(P1);
  x[1:k] = a1;
  
  for (t in 1:n) {
   
    real F = quad_form(P[1:k, 1:k], xreg[, t]) + Ht;
    real v = y[t] - dot_product(xreg[, t], head(x, k));
    vector[m] K = P[1:m, 1:k] * xreg[, t] / F;
    vector[m] tmpv = K * v;
    matrix[m, m] tmpm;
    vector[m] x_new;
    matrix[m, m] P_new;
   
    x_new[1:k] = head(x, k) + tail(x, k) + head(tmpv, k) + tail(tmpv, k);
    x_new[(k+1):m] = tail(x, k) + tail(tmpv, k);
    x = x_new;
    tmpm = P - K * K' * F;
    P_new[1:k, 1:k] = tmpm[1:k, 1:k] + tmpm[(k+1):m, (k+1):m] + tmpm[1:k, (k+1):m] + tmpm[(k+1):m, 1:k];
    P_new[1:k, (k+1):m] = tmpm[1:k, (k+1):m] + tmpm[(k+1):m, (k+1):m];
    P_new[(k+1):m, 1:k] = tmpm[(k+1):m, 1:k] + tmpm[(k+1):m, (k+1):m]; //P_new[1:k, (k+1):m];
    P_new[(k+1):m, (k+1):m] = tmpm[(k+1):m, (k+1):m];
    P = P_new + Rt;
    loglik = loglik - 0.5 * (log(F) + v * v / F);
  }
   return loglik;
  }
  
matrix gaussian_smoother_rw2(vector y, vector a1, vector P1, real Ht, matrix Rt, matrix xreg) {

  int n = rows(y);
  int k = rows(xreg);
  int m = 2 * k;
  vector[m] x = rep_vector(0, m);
  
  matrix[m, m] P = diag_matrix(P1);
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;
  matrix[m, m] Tt = diag_matrix(rep_vector(1.0, m));
  Tt[1:k, (k+1):m] = diag_matrix(rep_vector(1.0, k));
  x[1:k] = a1;
 
  for (t in 1:n) {
    vector[m] tmpv;
    matrix[m, m] tmpm;
    vector[m] x_new;
    matrix[m, m] P_new;
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht;
    v[t] = y[t] - dot_product(xreg[, t], head(x, k));
    K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
    tmpv = K[, t] * v[t];
    x_new[1:k] = head(x, k) + tail(x, k) + head(tmpv, k) + tail(tmpv, k);
    x_new[(k+1):m] = tail(x, k) + tail(tmpv, k);
    x = x_new;
    tmpm = P - K[, t] * K[, t]' * F[t];
    P_new[1:k, 1:k] = tmpm[1:k, 1:k] + tmpm[(k+1):m, (k+1):m] + tmpm[1:k, (k+1):m] + tmpm[(k+1):m, 1:k];
    P_new[1:k, (k+1):m] = tmpm[1:k, (k+1):m] + tmpm[(k+1):m, (k+1):m];
    P_new[(k+1):m, 1:k] = tmpm[(k+1):m, 1:k] + tmpm[(k+1):m, (k+1):m]; //P_new[1:k, (k+1):m];
    P_new[(k+1):m, (k+1):m] = tmpm[(k+1):m, (k+1):m];
    P = P_new + Rt;
    //x = x + K[, t] * v[t];
    //P = P - K[, t] * K[, t]' * F[t] + Rt;
  }
  
  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[,t+1];
    vector[m] tmp2 = rep_vector(0.0, m);
    tmp2[1:k] = xreg[, t];
    r[ ,t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
  }
  
  tmpr = r[,1];
  r[1:k, 1] = a1 + P1[1:k] .* tmpr[1:k];
  r[(k+1):m, 1] = P1[(k+1):m] .* tmpr[(k+1):m];
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = Tt * tmp + Rt * tmp2;
  }
  return r[1:m, 1:n];
  }
}

data {
  int<lower=0> k;
  int<lower=0> n;
  matrix[k, n] xreg;
  vector[n] y;
  vector[k] beta_mean;
  vector[k] beta_sd;
  vector[k] slope_sd;
  vector[k + 1] sigma_mean;
  vector[k + 1] sigma_sd;
  int<lower=0> n_new;
  matrix[k, n_new] xreg_new;
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
}

transformed parameters {
  vector[2 * k] R_vector;
  vector[2 * k] P1_vector;
  R_vector[1:k] = rep_vector(0.0, k);
  for(i in 1:k) {
    R_vector[i+k] = sigma_b[i]^2;
    P1_vector[i] = beta_sd[i]^2;
  }
  P1_vector[(k+1):(2*k)] = slope_sd;
}

model {
  sigma_b ~ normal(sigma_b_mean, sigma_b_sd);
  sigma_y ~ normal(sigma_y_mean, sigma_y_sd);
  y ~ gaussian_filter_rw2(beta_mean, P1_vector, sigma_y^2, diag_matrix(R_vector), xreg);
}

generated quantities{

  vector[n] y_rep;
  matrix[k, n] beta;
  matrix[k, n] slope;
  vector[n_new] y_new;
  matrix[k, n_new] beta_new;
  matrix[k, n_new] slope_new;

  // sample coefficients given sigma's (no conditioning on y)
  for(i in 1:k) {
     beta[i, 1] = normal_rng(beta_mean[i], beta_sd[i]);
     slope[i, 1] = normal_rng(0, slope_sd[i]);
  }

  for (t in 1:(n - 1)) {
    for(i in 1:k) {
      beta[i, t+1] = beta[i, t] + slope[i, t];
      slope[i, t+1] = normal_rng(slope[i, t], sigma_b[i]);
    }
  }
  // sample new observations given previously simulated beta
  for(t in 1:n) {
    y_rep[t] = normal_rng(dot_product(xreg[, t], beta[1:k, t]), sigma_y);
  }
  // perform mean correction to obtain sample from the posterior
  {
    matrix[2 * k, n] states = gaussian_smoother_rw2(y - y_rep, beta_mean, P1_vector, sigma_y^2, diag_matrix(R_vector), xreg);
    beta = beta + states[1:k, 1:n];
    slope = slope + states[(k + 1):(2 * k), 1:n];
  }

  // replicated data from posterior predictive distribution
  for(t in 1:n) {
    y_rep[t] = dot_product(xreg[, t], beta[1:k, t]) + normal_rng(0, sigma_y);
  }

  // prediction
  if (n_new > 0) {
    for(i in 1:k) {
      beta_new[i, 1] = beta[i, n] + slope[i, n];
      slope_new[i, 1] = normal_rng(slope[i, n], sigma_b[i]);
    }
    for(t in 1:(n_new - 1)) {
      y_new[t] = dot_product(xreg_new[, t], beta_new[, t]) + normal_rng(0, sigma_y);
      for(i in 1:k) {
        beta_new[i, t+1] = beta_new[i, t] + slope_new[i, t];
        slope_new[i, t+1] = normal_rng(slope_new[i, t], sigma_b[i]);
      }
    }
    y_new[n_new] = dot_product(xreg_new[, n_new], beta_new[, n_new]) + normal_rng(0, sigma_y);
  }
}
