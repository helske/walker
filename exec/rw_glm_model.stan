functions {
// Functions for Kalman filter and smoother for dynamic regression model
// note that these functions are not fully optimised yet

// univariate Kalman filter, returns the log-likelihood
vector glm_approx_smoother(vector y, vector a1, vector P1, vector Ht, matrix Rt, 
  matrix xreg, int distribution, vector u, vector y_original) {

  int n = rows(y);
  int m = rows(a1);
  vector[m] x = a1;
  matrix[m, m] P = diag_matrix(P1);
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n + 1] r;
  vector[m] tmpr;
  vector[2] loglik = rep_vector(0, 2);
  
  // filter
  for (t in 1:n) {
    F[t] = quad_form(P, xreg[, t]) + Ht[t];
    v[t] = y[t] - dot_product(xreg[, t], x);
    K[, t] = P * xreg[, t] / F[t];
    x = x + K[, t] * v[t];
    P = P - K[, t] * K[, t]' * F[t] + Rt;
    loglik[1] = loglik[1] - 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
  }
  // smoother
  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[,t+1];
    r[,t] =  xreg[, t] * v[t] / F[t] + tmp - xreg[, t] * dot_product(K[,t], tmp);
  }

  tmpr = r[,1];
  r[,1] = a1 + P1 .* tmpr;
  for (t in 2:n) {
    vector[m] tmp = r[,t-1];
    vector[m] tmp2 = r[,t];
    r[,t] = tmp + Rt * tmp2;
  }

  // add a correction term
  // Poisson case, generalization to binomial etc straightforward, see for example Durbin and Koopman 2012
 // if (distribution == 1) {
    for(t in 1:n) {
      real xbeta = dot_product(xreg[,t], r[, t]);
      loglik[2] = loglik[2] + y_original[t] * xbeta - u[t] * exp(xbeta) +
          0.5 * (y[t] - xbeta)^2 / Ht[t];
    }
  //}  
  return loglik;
}
  
matrix gaussian_smoother(vector y, vector a1, vector P1, vector Ht, matrix Rt, matrix xreg) {

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
    F[t] = quad_form(P, xreg[, t]) + Ht[t];
    v[t] = y[t] - dot_product(xreg[, t], x);
    K[, t] = P * xreg[, t] / F[t];
    x = x + K[, t] * v[t];
    P = P - K[, t] * K[, t]' * F[t] + Rt;
  }
  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[,t+1];
    r[,t] =  xreg[, t] * v[t] / F[t] + tmp - xreg[, t] * dot_product(K[,t], tmp);
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
  matrix[k, n] xreg;
  vector[n] y;
  vector[k] beta_mean;
  vector[k] beta_sd;
  vector[k] sigma_mean;
  vector[k] sigma_sd;
  int<lower=0> n_new;
  matrix[k, n_new] xreg_new;
  vector[n] Ht;
  vector[n] y_original;
  vector[n] u;
  int distribution;
  int<lower=0> N;
}


parameters {
  real<lower=0> sigma_b[k];
}

transformed parameters {
  vector[k] R_vector;
  vector[k] P1_vector;
  vector[2] loglik;
  for(i in 1:k) {
    R_vector[i] = sigma_b[i]^2;
    P1_vector[i] = beta_sd[i]^2;
  }
  loglik = glm_approx_smoother(y, beta_mean, P1_vector, Ht, diag_matrix(R_vector), xreg, 
  distribution, u, y_original);
}

model {
  sigma_b ~ normal(sigma_mean, sigma_sd);
  target += sum(loglik);
}

generated quantities{
  
  real beta[k, n];
  real weights;
  
  {  
  vector[n] y_rep;
  matrix[k, n] beta_j;
  real beta_array[k, n, N];
  //vector[n_new] y_new; not used yet
  //vector[k] beta_new;
  
  vector[N] w; //importance sampling weights
  // This is the simplest but not most efficient way to sample multiple realizations
  // We can save a lot by running only one full Kalman smoother and then doing some 
  // tricks (see for example DK2002, and implementations in KFAS and bssm)
  for(j in 1:N) {
    // sample coefficients given sigma's (no conditioning on y)  
    for(i in 1:k) {
       beta_j[i, 1] = normal_rng(beta_mean[i], beta_sd[i]);
    }
    for (t in 1:(n - 1)) {
      for(i in 1:k) {
        beta_j[i, t+1] = normal_rng(beta_j[i, t], sigma_b[i]);
      }
    }
    // sample new observations given previously simulated beta
    for(t in 1:n) {
      y_rep[t] = normal_rng(dot_product(xreg[, t], beta_j[1:k, t]), sqrt(Ht[t]));
    }
    // perform mean correction to obtain sample from the posterior
    beta_j = beta_j + gaussian_smoother(y - y_rep, beta_mean, P1_vector, Ht, diag_matrix(R_vector), xreg);
    beta_array[1:k,1:n,j] = to_array_2d(beta_j);
    
    w[j] = -loglik[2];
    //if (distribution == 1) {
      for(t in 1:n) {
        real xbeta = dot_product(xreg[,t], beta_j[1:k,t]);
        w[j] = w[j] + y_original[t] * xbeta - u[t] * exp(xbeta) +
            0.5 * (y[t] - xbeta)^2 / Ht[t];
      }
    //}
  }
    w = exp(w);
    weights = mean(w);
    beta = beta_array[1:k,1:n,categorical_rng(w / sum(w))];
  
  }
  // These parts need bit more work (currently use only approximation)
  // The algorithm is in the KFAS paper
  
  // // replicated data from posterior predictive distribution
  // for(t in 1:n) {
  //   y_rep[t] = poisson_log_rng(dot_product(xreg[, t], beta[1:k, t]));
  // }
  // 
  // // prediction 
  // if (n_new > 0) {
  //   for(i in 1:k) {
  //     beta_new[i] = normal_rng(beta[i, n], sigma_b[i]);
  //   }
  //   for(t in 1:n_new) {
  //     y_new[t] = poisson_log_rng(dot_product(xreg_new[, t], beta_new));
  //     for(i in 1:k) {
  //       beta_new[i] = normal_rng(beta_new[i], sigma_b[i]);
  //     } 
  //   }
  // }
}
