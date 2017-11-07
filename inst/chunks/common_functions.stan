// Functions for Kalman filter and smoother for dynamic regression model
// note that these functions are not fully optimised yet

// univariate Kalman filter for RW1+RW2 model, returns the log-likelihood
real gaussian_filter(vector y, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg, real[] gamma) {
  
  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;

  vector[m] x = a1;
  matrix[m, m] P = P1;

  for (t in 1:n) {

    real F = quad_form(P[1:k, 1:k], xreg[, t]) + Ht;
    
    if (F > 1.0e-12) { // protect against numerical issues
      real v = y[t] - dot_product(xreg[, t], head(x, k));
      vector[m] K = P[1:m, 1:k] * xreg[, t] / F;
      x = Tt * (x + K * v);
      P = quad_form_sym(P - K * K' * F, Tt') + gamma[t]^2 * Rt;
      loglik = loglik - 0.5 * (log(F) + v * v / F);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + gamma[t]^2 * Rt;
    }
  }
  return loglik;
  
}

matrix gaussian_smoother(vector y, vector a1, matrix P1, real Ht, 
  matrix Tt, matrix Rt, matrix xreg, real[] gamma) {

  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  real loglik = 0.0;
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;
  
  for (t in 1:n) {
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht;
    
    if (F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt') + gamma[t]^2 * Rt;
      loglik = loglik - 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + gamma[t]^2 * Rt;
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(F[t] > 1.0e-12) {
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
    r[,t] = Tt * tmp + gamma[t]^2 * Rt * tmp2;
  }
  return r[1:m, 1:n];
}



// univariate Kalman filter & smoother for non-gaussian model, 
// returns the log-likelihood of the corresponding approximating Gaussian model
// and a extra correction term
vector glm_approx_loglik(vector y, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg, int distribution, int[] u, 
  vector y_original, vector xbeta_fixed, real[] gamma) {

  int k = rows(xreg);
  int n = rows(y);
  int m = rows(a1);
  vector[2] loglik = rep_vector(0.0, 2);
  vector[m] x = a1;
  matrix[m, m] P = P1;
  vector[n] v;
  vector[n] F;
  matrix[m, n] K;
  matrix[m, n+1] r;
  vector[m] tmpr;
  
  for (t in 1:n) {
    
    F[t] = quad_form(P[1:k, 1:k], xreg[, t]) + Ht[t];
    
    if (F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt') + gamma[t]^2 * Rt;
      loglik[1] = loglik[1] - 0.5 * (log(F[t]) + v[t] * v[t] / F[t]);
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + gamma[t]^2 * Rt;
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(F[t] > 1.0e-12) {
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
    r[,t] = Tt * tmp + gamma[t]^2 * Rt * tmp2;
  }

  // add a correction term
  // Poisson case, generalization to binomial etc straightforward, see for example Durbin and Koopman 2012
  if (distribution == 1) {
    for(t in 1:n) {
      real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
      loglik[2] = loglik[2] + y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
      u[t] * exp(xbeta_rw + xbeta_fixed[t]) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
    }
  } else {
    for(t in 1:n) {
     real xbeta_rw = dot_product(xreg[,t], r[1:k, t]);
      loglik[2] = loglik[2] + y_original[t] * (xbeta_rw + xbeta_fixed[t]) - 
      u[t] * log1p(exp(xbeta_rw + xbeta_fixed[t])) +
          0.5 * (y[t] - xbeta_rw)^2 / Ht[t];
    }
  }
  return loglik;
}


// univariate Kalman filter & smoother for non-gaussian model, 
// returns the log-likelihood of the corresponding approximating Gaussian model
// and a extra correction term
matrix glm_approx_smoother(vector y, vector a1, matrix P1, vector Ht, 
  matrix Tt, matrix Rt, matrix xreg, real[] gamma) {

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
    
    if (F[t] > 1.0e-12) {
      v[t] = y[t] - dot_product(xreg[, t], head(x, k));
      K[, t] = P[1:m, 1:k] * xreg[, t] / F[t];
      x = Tt * (x + K[,t] * v[t]);
      P = quad_form_sym(P - K[,t] * K[,t]' * F[t], Tt') + gamma[t]^2 * Rt;
    } else {
      x = Tt * x;
      P = quad_form_sym(P, Tt') + gamma[t]^2 * Rt;
    }
  }

  r[,n+1] = rep_vector(0.0, m);
  for (tt in 1:n) {
    int t = n + 1 - tt;
    vector[m] tmp = r[, t+1];
    if(F[t] > 1.0e-12) {
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
    r[,t] = Tt * tmp + gamma[t]^2 * Rt * tmp2;
  }

  return r[1:m, 1:n];
}

