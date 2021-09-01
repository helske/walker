functions {
  // Functions for Kalman filter and smoother for dynamic regression model
  // note that these functions are not fully optimised yet
  
  // univariate Kalman filter for RW1+RW2 model, returns the log-likelihood
  vector gaussian_filter(vector y, int[] y_miss, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg, vector gamma2_y) {
    
    int k = rows(xreg);
    int n = rows(y);
    int m = rows(a1);
    vector[n] loglik = rep_vector(0, n);
    
    vector[m] x = a1;
    matrix[m, m] P = P1;
    
    real log2pi = log(2*pi());
    
    for (t in 1:n) {
      real F = quad_form(P[1:k, 1:k], xreg[, t]) + gamma2_y[t] * Ht;
      
      if (y_miss[t] == 0) {
        real v = y[t] - dot_product(xreg[, t], head(x, k));
        vector[m] K = P[1:m, 1:k] * xreg[, t] / F;
        x = Tt * (x + K * v);
        P = quad_form_sym(P - K * K' * F, Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
        loglik[t] = -0.5 * (log2pi + log(F) + v * v / F);
      } else {
        x = Tt * x;
        P = quad_form_sym(P, Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
      }
    }
    return loglik;
    
  }
  
  matrix gaussian_smoother(vector y, int[] y_miss, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg,vector gamma2_y) {
    
    int k = rows(xreg);
    int n = rows(y);
    int m = rows(a1);
    real loglik = 0.0;
    vector[m] x = a1;
    matrix[m, m] P = P1;
    vector[n] v;
    vector[n] F;
    matrix[m, n] K;
    matrix[m, n + 1] r;
    vector[m] tmpr;
    
    for (t in 1:n) {
      
      F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + gamma2_y[t] * Ht;
      
      if (y_miss[t] == 0 && F[t] > 1.0e-12) {
        v[t] = y[t] - dot_product(xreg[, t], head(x, k));
        K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
        x = Tt * (x + K[,t] * v[t]);
        P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
        loglik -= 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
      } else {
        x = Tt * x;
        P = quad_form_sym(P, Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
      }
    }
    
    r[,n+1] = rep_vector(0.0, m);
    for (tt in 1:n) {
      int t = n + 1 - tt;
      vector[m] tmp = r[, t+1];
      if(y_miss[t] == 0 && F[t] > 1.0e-12) {
        vector[m] tmp2 = rep_vector(0.0, m);
        tmp2[1:k] = xreg[, t];
        r[ ,t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
      } else {
        r[,t] = Tt' * tmp;
      }
    }
    
    tmpr = r[,1];
    r[,1] = a1 + P1 * tmpr;
    for (t in 2:n) {
      vector[m] tmp = r[,t-1];
      vector[m] tmp2 = r[,t];
      r[,t] = Tt * tmp + Rt[, t] .* tmp2;
    }
    return r[1:m, 1:n];
  }
  
}

data {
  int<lower=0> k_fixed;
  int<lower=0> k_rw1;
  int<lower=0> k_rw2;
  int<lower=0> m;
  int<lower=0> k;
  int<lower=1> n;
  int<lower=1> n_lfo;
  matrix[n, k_fixed] xreg_fixed;
  matrix[k, n] xreg_rw;
  vector[n] y;
  int<lower=0> y_miss[n];
  real<lower=0> sigma_y_shape;
  real<lower=0> sigma_y_rate;
  
  real beta_fixed_mean;
  real<lower=0> beta_fixed_sd;
  real beta_rw1_mean;
  real<lower=0> beta_rw1_sd;
  real beta_rw2_mean;
  real<lower=0> beta_rw2_sd;
  
  real<lower=0> sigma_rw1_shape;
  real<lower=0> sigma_rw2_shape;
  real<lower=0> sigma_rw1_rate;
  real<lower=0> sigma_rw2_rate;
  
  real<lower=0> nu_mean;
  real<lower=0> nu_sd;
  vector[n] gamma_y;
  matrix[k_rw1, n] gamma_rw1;
  matrix[k_rw2, n] gamma_rw2;
}

transformed data {
  vector[m] a1;
  matrix[m, m] P1 = rep_matrix(0.0, m, m);
  matrix[m, m] Tt = diag_matrix(rep_vector(1.0, m));
  vector[n] gamma2_y = gamma_y .* gamma_y;
  
  if(k_rw2 > 0) {
    Tt[(k_rw1+1):k, (k+1):m] = diag_matrix(rep_vector(1.0, k_rw2));
  }
  for(i in 1:k_rw1) {
    a1[i] = beta_rw1_mean;
    P1[i, i] = beta_rw1_sd^2;
  }
  for(i in (k_rw1 + 1):k) {
    a1[i] = beta_rw2_mean;
    P1[i, i] = beta_rw2_sd^2;
  }
  for(i in (k + 1):m) {
    a1[i] = nu_mean;
    P1[i, i] = nu_sd^2;
  }
}

parameters {
  vector[k_fixed] beta_fixed;
  real<lower=0> sigma_rw1[k_rw1];
  real<lower=0> sigma_rw2[k_rw2];
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[n] log_lik;
  matrix[m, n] Rt = rep_matrix(0.0, m, n);
  vector[n] xbeta;
  vector[n] y_;
  
  if (k_fixed > 0) {
    xbeta = xreg_fixed * beta_fixed;
  } else {
    xbeta = rep_vector(0.0, n);
  }
  y_ = y - xbeta;
  for (t in 1:n) {
    for(i in 1:k_rw1) {
      Rt[i, t] = (gamma_rw1[i, t] * sigma_rw1[i])^2;
    }
    for(i in 1:k_rw2) {
      Rt[k + i, t] = (gamma_rw2[i, t] * sigma_rw2[i])^2;
    } 
  }
  log_lik = gaussian_filter(y_, y_miss, a1, P1, sigma_y^2, Tt, Rt, xreg_rw, gamma2_y);
}

model {
  beta_fixed ~ normal(beta_fixed_mean, beta_fixed_sd);
  sigma_y ~ gamma(sigma_y_shape, sigma_y_rate);
  sigma_rw1 ~ gamma(sigma_rw1_shape, sigma_rw1_rate);
  sigma_rw2 ~ gamma(sigma_rw2_shape, sigma_rw2_rate);
  
  target += sum(log_lik[1:n_lfo]);
}

generated quantities{
  vector[n * (n_lfo == n)] y_rep;
  matrix[k, n * (n_lfo == n)] beta_rw;
  matrix[k_rw2, n * (n_lfo == n)] nu;
  vector[n * (n_lfo == n)] y_fit;
  
  if(n_lfo == n) {
    // sample coefficients given sigma's (no conditioning on y)
    
    for(i in 1:k_rw1) {
      beta_rw[i, 1] = normal_rng(beta_rw1_mean, beta_rw1_sd);
    }
    for(i in 1:k_rw2) {
      beta_rw[k_rw1 + i, 1] = normal_rng(beta_rw2_mean, beta_rw2_sd);
      nu[i, 1] = normal_rng(nu_mean, nu_sd);
    }
    
    for (t in 1:(n - 1)) {
      for(i in 1:k_rw1) {
        beta_rw[i, t + 1] = normal_rng(beta_rw[i, t], gamma_rw1[i, t] * sigma_rw1[i]);
      }
      for(i in 1:k_rw2) {
        beta_rw[k_rw1 + i, t+1] = beta_rw[k_rw1 + i, t] + nu[i, t];
        nu[i, t + 1] = normal_rng(nu[i, t], gamma_rw2[i, t] * sigma_rw2[i]);
      }
    }
    // sample new observations given previously simulated beta
    for(t in 1:n) {
      y_rep[t] = normal_rng(dot_product(xreg_rw[, t], beta_rw[, t]), gamma_y[t] * sigma_y);
    }
    // perform mean correction to obtain sample from the posterior
    {
      matrix[m, n] states = gaussian_smoother(y_ - y_rep, y_miss, a1, P1,
      sigma_y^2, Tt, Rt, xreg_rw, gamma2_y);
      beta_rw += states[1:k, 1:n];
      if(k_rw2 > 0) nu += states[(k + 1):m, 1:n];
    }
    
    // replicated data from posterior predictive distribution
    for(t in 1:n) {
      y_fit[t] = xbeta[t] + dot_product(xreg_rw[, t], beta_rw[, t]);
      y_rep[t] = normal_rng(y_fit[t], gamma_y[t] * sigma_y);
    }
  }
}
