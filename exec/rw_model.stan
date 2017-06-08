functions {
// Functions for Kalman filter and smoother for dynamic regression model
// note that these functions are not fully optimised yet

// univariate Kalman filter, returns the log-likelihood
real gaussian_filter(vector y, vector a1, vector P1, real Ht, matrix Rt, matrix xreg) {

  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;

  vector[m] x = a1;
  matrix[m, m] P = diag_matrix(P1);
  for (t in 1:n) {
    real F = xreg[t,] * P * xreg[t,]' + Ht;
    if (F > 1.0e-8) {
      real v = y[t] - xreg[t,] * x;
      vector[m] K =  P * xreg[t, ]' / F;
      x = x + K * v;
      P = P - K * K' * F + Rt;
      loglik = loglik - 0.5 * (log(F) + v * v / F);
    } else {
      P = P + Rt;
    }
  }
   return loglik;
  }
  
matrix gaussian_smoother(vector y, vector a1, vector P1, real Ht, matrix Rt, matrix xreg) {

  int n = rows(y);
  int m = rows(a1);
  vector[m] x = a1;
  matrix[m, m] P = diag_matrix(P1);
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;

  for (t in 1:n) {
    F[t] = xreg[t,] * P * xreg[t,]' + Ht;
    if (F[t] > 1.0e-8) {
      v[t] = y[t] - xreg[t,] * x;
      K[, t] = P * xreg[t,]' / F[t];
      x = x + K[,t] * v[t];
      P = P - K[,t] * K[,t]' * F[t] + Rt;
    } else {
      P = P + Rt;
    }
  }
 
  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[,t+1];
    if(F[t] > 1.0e-8) {
      r[,t] =  xreg[t,]' * v[t] / F[t] + tmp - (K[,t] * xreg[t,])' * tmp;
    }
  }

  tmpr = r[,1];
  r[,1] = a1 + P1 .* tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = tmp + Rt * tmp2;
  }
  return r[1:m, 1:n];
  }
}

data {
  int<lower=0> k;
  int<lower=0> n;
  matrix[n, k] xreg;
  vector[n] y;
  vector[k] beta_mean;
  vector[k] beta_sd;
  vector[k + 1] sigma_mean;
  vector[k + 1] sigma_sd;
  int<lower=0> n_new;
  matrix[n_new, k] xreg_new;
}

parameters {
  real<lower=0> sigma[1 + k];
}
transformed parameters {
  vector[k] R_vector;
  vector[k] P1_vector;
  for(i in 1:k) {
    R_vector[i] = sigma[1 + i]^2;
    P1_vector[i] = beta_sd[i]^2;
  }
}
model {
  target += normal_lpdf(sigma | sigma_mean, sigma_sd);
  target += gaussian_filter(y, beta_mean, P1_vector, sigma[1]^2, diag_matrix(R_vector), xreg);
}

generated quantities{
    
  vector[n] y_rep;
  matrix[k, n] beta;
  vector[n_new] y_new;
  vector[k] beta_new;
  
  // sample coefficients given sigma's (no conditioning on y)  
  for(i in 1:k) {
     beta[i, 1] = normal_rng(beta_mean[i], beta_sd[i]);
  }
  for (t in 1:(n - 1)) {
    for(i in 1:k) {
      beta[i, t+1]  = normal_rng(beta[i, t], sigma[1 + i]);
    }
  }
  // sample new observations given previously simulated beta
  for(t in 1:n) {
    y_rep[t] = xreg[t,] * beta[1:k, t] + normal_rng(0, sigma[1]);
  }
  // perform mean correction to obtain sample from the posterior
  beta = beta + gaussian_smoother(y - y_rep, beta_mean, P1_vector, sigma[1]^2, diag_matrix(R_vector), xreg);
  
  // replicated data from posterior predictive distribution
  for(t in 1:n) {
    y_rep[t] = xreg[t,] * beta[1:k, t] + normal_rng(0, sigma[1]);
  }
  
  // prediction 
  if (n_new > 0) {
    for(i in 1:k) {
      beta_new[i] = normal_rng(beta[i, n], sigma[1 + i]);
    }
    for(t in 1:n_new) {
      y_new[t] = xreg_new[t,] * beta_new + normal_rng(0, sigma[1]);
      for(i in 1:k) {
        beta_new[i] = normal_rng(beta_new[i], sigma[1 + i]);
      } 
    }
  }
}
