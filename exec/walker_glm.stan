// Dynamic regression for non-Gaussian models
//
// The importance sampling type correction is introduced in:
// Vihola M, Helske J and Franks J (2016). “Importance sampling type
// correction of Markov chain Monte Carlo and exact approximations.”
// On ArXiv: https://arxiv.org/abs/1609.02541

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
  int<lower=0> y_miss[n];
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
  
  vector[n] Ht;
  vector[n] y_original;
  int<lower=0> u[n];
  int distribution;
  int<lower=0> N;
  matrix[k_rw1, n] gamma_rw1;
  matrix[k_rw2, n] gamma_rw2;
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
}

transformed parameters {
  matrix[n, m] Rt = rep_matrix(0.0, m, n);
  vector[n] xbeta;
  vector[n] y_;
  vector[2] loglik;
   
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
   
  loglik = glm_approx_loglik(y_, y_miss, a1, P1, Ht, 
    Tt, Rt, xreg_rw, distribution, u, y_original, xbeta);

}

model {
  beta_fixed ~ normal(beta_fixed_mean, beta_fixed_sd);
  sigma_rw1 ~ normal(sigma_rw1_mean, sigma_rw1_sd);
  sigma_rw2 ~ normal(sigma_rw2_mean, sigma_rw2_sd);
  target += sum(loglik);
}


generated quantities{
  
  matrix[k, n] beta_rw;
  matrix[k_rw2, n] slope;
  real weights;
  vector[n] y_fit;
  vector[n] y_rep;

  {

    vector[n] y_rep_j;
    matrix[k, n] beta_j;
    matrix[k_rw2, n] slope_j;
    real beta_array[k, n, N];
    real slope_array[k_rw2, n, N];
    vector[N] w = rep_vector(0.0, N); //importance sampling weights
    
    // This is the simplest but not most efficient way to sample multiple realizations
    // We could save a lot by running only one full Kalman smoother and then doing some
    // tricks (see for example Durbin and Koopman (2002), 
    // and implementations in KFAS (in Fortran) and in bssm (in C++))
    
    for(j in 1:N) {
  
      for(i in 1:k_rw1) {
        beta_j[i, 1] = normal_rng(beta_rw1_mean, beta_rw1_sd);
      }
      for(i in 1:k_rw2) {
        beta_j[k_rw1 + i, 1] = normal_rng(beta_rw2_mean, beta_rw2_sd);
        slope_j[i, 1] = normal_rng(slope_mean, slope_sd);
      }

      for (t in 1:(n - 1)) {
        for(i in 1:k_rw1) {
          beta_j[i, t + 1] = normal_rng(beta_j[i, t], gamma_rw1[i, t] * sigma_rw1[i]);
        }
        for(i in 1:k_rw2) {
          beta_j[k_rw1 + i, t+1] = beta_j[k_rw1 + i, t] + slope_j[i, t];
          slope_j[i, t + 1] = normal_rng(slope_j[i, t], gamma_rw2[i, t] * sigma_rw2[i]);
        }
      }
      // sample new observations given previously simulated beta
      for(t in 1:n) {
        y_rep_j[t] = normal_rng(dot_product(xreg_rw[, t], beta_j[, t]), sqrt(Ht[t]));
      }
      // perform mean correction to obtain sample from the posterior
      {
        matrix[m, n] states = glm_approx_smoother(y_ - y_rep_j, y_miss, a1, P1,
          Ht, Tt, Rt, xreg_rw);
        beta_j += states[1:k, 1:n];
        slope_j += states[(k + 1):m, 1:n];
      }
  
      beta_array[1:k,1:n,j] = to_array_2d(beta_j);
      slope_array[1:k_rw2,1:n,j] = to_array_2d(slope_j);

      w[j] = -loglik[2];
      if (distribution == 1) {
        for(t in 1:n) {
          real xbeta_tmp = xbeta[t] + dot_product(xreg_rw[,t], beta_j[1:k,t]);
          w[j] += y_original[t] * xbeta_tmp - u[t] * exp(xbeta_tmp) +
            0.5 * (y[t] - xbeta_tmp)^2 / Ht[t];
        }
      } else {
        for(t in 1:n) {
          real xbeta_tmp = xbeta[t] + dot_product(xreg_rw[,t], beta_j[1:k,t]);
          w[j] += y_original[t] * xbeta_tmp - u[t] * log1p(exp(xbeta_tmp)) +
            0.5 * (y[t] - xbeta_tmp)^2 / Ht[t];
        }
      }
    }
    
    {
      // store only one of the simulated states
      // we could store all as well but this is a compromise between
      // space and accuracy. Of course, we could compute the summary 
      // statistics of the states here already using all replications...
      
      // note that the results will be a weighted posterior sample
      int index;
      vector[N] expw = exp(w);
      weights = mean(expw);
      index = categorical_rng(expw / sum(expw));
      beta_rw = to_matrix(beta_array[, , index]);
      if (k_rw2 > 0) slope = to_matrix(slope_array[, , index]);
   
    // replicated data from posterior predictive distribution
   
      if (distribution == 1) {
        for(t in 1:n) {
           y_fit[t] = u[t] * exp(xbeta[t] + dot_product(xreg_rw[, t], beta_rw[, t]));
           y_rep[t] = poisson_rng(y_fit[t]);
        }
      } else {
        for(t in 1:n) {
          real tmp = exp(xbeta[t] + dot_product(xreg_rw[, t], beta_rw[, t]));
          y_fit[t] = tmp / (1.0 + tmp);
          y_rep[t] = binomial_rng(u[t], y_fit[t]);
        }
      }
    }
    
  }
}
