// Dynamic regression for non-Gaussian models
//
// The importance sampling type correction is introduced in:
// Vihola, M, Helske, J, Franks, J. 
// Importance sampling type estimators based on approximate 
// marginal Markov chain Monte Carlo. Scand J Statist. 2020; 
// 47: 1339– 1376. https://doi.org/10.1111/sjos.12492

functions {
  
  // univariate Kalman filter & smoother for non-gaussian model, 
  // returns the log-likelihood of the corresponding approximating Gaussian model
  // and a extra correction term
  matrix glm_approx_loglik(vector y, int[] y_miss, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg, int distribution, int[] u, 
  vector y_original, vector xbeta_fixed) {
    
    int k = rows(xreg);
    int n = rows(y);
    int m = rows(a1);
    matrix[n,2] loglik = rep_matrix(0, n, 2);
    vector[m] x = a1;
    matrix[m, m] P = P1;
    vector[n] v;
    vector[n] F;
    matrix[m, n] K;
    matrix[m, n+1] r;
    vector[m] tmpr;
    
    real log2pi = log(2*pi());
        
    for (t in 1:n) {
      
      F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht[t];
      
      if (y_miss[t] == 0) {
        v[t] = y[t] - dot_product(xreg[, t], head(x, k));
        K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
        x = Tt * (x + K[,t] * v[t]);
        P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
        loglik[t, 1] = -0.5 * (log2pi + log(F[t]) + v[t] * v[t] / F[t]);
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
      if(y_miss[t] == 0) {
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
    
    // add a correction term
    if (distribution == 1) {
      for(t in 1:n) {
        if (y_miss[t] == 0){
          real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
          loglik[t, 2] = y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
          u[t] * exp(xbeta_rw + xbeta_fixed[t]) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
        }
      }
    } else {
      for(t in 1:n) {
        if (y_miss[t] == 0){
          real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
          loglik[t, 2] = y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
          u[t] * log1p(exp(xbeta_rw + xbeta_fixed[t])) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
        }
      }
    }
    return loglik;
  }
  
  
  // univariate Kalman filter & smoother for non-gaussian model, 
  // returns the log-likelihood of the corresponding approximating Gaussian model
  // and a extra correction term
  matrix glm_approx_smoother(vector y, int[] y_miss, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg) {
    
    int k = rows(xreg);
    int n = rows(y);
    int m = rows(a1);
    vector[m] x = a1;
    matrix[m, m] P = P1;
    vector[n] v;
    vector[n] F;
    matrix[m, n] K;
    matrix[m, n+1] r;
    vector[m] tmpr;
    
    for (t in 1:n) {
      
      F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht[t];
      
      if (y_miss[t] == 0) {
        v[t] = y[t] - dot_product(xreg[, t], head(x, k));
        K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
        x = Tt * (x + K[,t] * v[t]);
        P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
      } else {
        x = Tt * x;
        P = quad_form_sym(P, Tt');
        for (i in 1:m) {
          P[i, i] += Rt[i, t];
        }
      }
    }
    
    r[, n + 1] = rep_vector(0.0, m);
    for (tt in 1:n) {
      int t = n + 1 - tt;
      vector[m] tmp = r[, t+1];
      if(y_miss[t] == 0) {
        vector[m] tmp2 = rep_vector(0.0, m);
        tmp2[1:k] = xreg[, t];
        r[, t] = tmp2 * v[t] / F[t] + (Tt - Tt * K[,t] * tmp2')' * tmp;
      } else {
        r[, t] = Tt' * tmp;
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
}

transformed parameters {
  matrix[m, n] Rt = rep_matrix(0.0, m, n);
  vector[n] xbeta;
  vector[n] y_;
  matrix[n, 2] loglik;
  
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
  sigma_rw1 ~ gamma(sigma_rw1_shape, sigma_rw1_rate);
  sigma_rw2 ~ gamma(sigma_rw2_shape, sigma_rw2_rate);
  // n_lfo is for easier computation of leave-future-out cross-validation
  // use only observations up to n_lfo (default n_lfo = n)
  // (not very optimal way to do this though, we running KF etc with y_1:n...)
  target += sum(loglik[1:n_lfo, 1:2]);
}


generated quantities{
  
  matrix[k, n * (n_lfo == n)] beta_rw;
  matrix[k_rw2, n * (n_lfo == n)] nu;
  real weights;
  vector[n * (n_lfo == n)] y_fit;
  vector[n * (n_lfo == n)] y_rep;
  // approximate log-likelihood, unbiased estimate is log_lik + mean(w)
  vector[n] log_lik = loglik[, 1] + loglik[, 2];
  
  if(n_lfo == n) {
    
    vector[n] y_rep_j;
    matrix[k, n] beta_j;
    matrix[k_rw2, n] nu_j;
    real beta_array[k, n, N];
    real nu_array[k_rw2, n, N];
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
        nu_j[i, 1] = normal_rng(nu_mean, nu_sd);
      }
      
      for (t in 1:(n - 1)) {
        for(i in 1:k_rw1) {
          beta_j[i, t + 1] = normal_rng(beta_j[i, t], gamma_rw1[i, t] * sigma_rw1[i]);
        }
        for(i in 1:k_rw2) {
          beta_j[k_rw1 + i, t+1] = beta_j[k_rw1 + i, t] + nu_j[i, t];
          nu_j[i, t + 1] = normal_rng(nu_j[i, t], gamma_rw2[i, t] * sigma_rw2[i]);
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
        if(k_rw2 > 0) nu_j += states[(k + 1):m, 1:n];
      }
      
      beta_array[1:k,1:n,j] = to_array_2d(beta_j);
      if(k_rw2 > 0) nu_array[1:k_rw2,1:n,j] = to_array_2d(nu_j);
      
      w[j] = -sum(loglik[,2]);
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
      if (k_rw2 > 0) nu = to_matrix(nu_array[, , index]);
      
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
