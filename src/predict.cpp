#include "predict.h"
#include <RcppArmadilloExtensions/sample.h>

// sigma_rw: k x N matrix of posterior samples of rw_sigmas
// sigma_y: vector of length N of posterior samples of sigma_y
// beta_fixed p x N matrix of posterior samples of fixed betas
// beta_rw: k x N matrix of posterior samples of beta_rws at time n
// xreg: k x n_new matrix of new covariates
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List predict_walker(const arma::mat& sigma_rw1, 
  const arma::mat& sigma_rw2, const arma::vec sigma_y,
  const arma::mat beta_fixed, const arma::mat& beta_rw, const arma::mat& slope,
  const arma::mat& xreg_fixed, const arma::mat& xreg_rw, 
  const arma::uword n, const arma::uword k, const arma::uword k_rw1, const arma::uword k_rw2,
  const bool response) {
  
  arma::uword k_rw = k_rw1 + k_rw2;
  arma::uword n_iter = sigma_y.n_elem;
  
  arma::cube beta_new(k_rw, n, n_iter);
  arma::cube slope_new(k_rw2, n, n_iter);
  arma::mat y(n, n_iter);
  
  // for each realization from posterior
  for (arma::uword i = 0; i < n_iter; i++) {
    
    // sample the states at first time point
    for (arma::uword j = 0; j < k_rw1; j++) {
      beta_new(j, 0, i) = R::rnorm(beta_rw(j, i), sigma_rw1(j, i));
    }
    for (arma::uword j = 0; j < k_rw2; j++) {
      beta_new(k_rw1 + j, 0, i) = beta_rw(k_rw1 + j, i) + slope(j, i);
      slope_new(j, 0, i) = R::rnorm(slope(j, i), sigma_rw2(j, i));
    }
    
    for(arma::uword t = 0; t < (n - 1); t++) {
      // sample observations
      y(t, i) = arma::dot(xreg_rw.col(t), beta_new.slice(i).col(t));
      if (response) {
        y(t, i) += R::rnorm(0, sigma_y(i)); 
      }
      // and states
      for (arma::uword j = 0; j < k_rw1; j++) {
        beta_new(j, t + 1, i) = R::rnorm(beta_new(j, t, i), sigma_rw1(j, i));
      }
      for (arma::uword j = 0; j < k_rw2; j++) {
        beta_new(k_rw1 + j, t + 1, i) = beta_new(k_rw1 + j, t, i) + 
          slope_new(j, t, i);
        slope_new(j, t + 1, i) = R::rnorm(slope_new(j, t, i), sigma_rw2(j, i));
      }
    }
    // and the observations at last time point
    y(n - 1, i) = arma::dot(xreg_rw.col(n - 1), beta_new.slice(i).col(n - 1));
    if (response) {
      y(n - 1, i) += R::rnorm(0, sigma_y(i)); 
    }
  }
  
  if (k > 0) {
    for (arma::uword i = 0; i < n_iter; i++) {
      y.col(i) += xreg_fixed * beta_fixed.col(i);
    }
  }
  return Rcpp::List::create(Rcpp::Named("y_new") = y, 
    Rcpp::Named("beta_new") = beta_new, Rcpp::Named("slope_new") = slope_new);
}

// as above, but for non-Gaussian model
// [[Rcpp::export]]
Rcpp::List predict_walker_glm(const arma::mat& sigma_rw1, 
  const arma::mat& sigma_rw2,
  const arma::mat beta_fixed, const arma::mat& beta_rw, const arma::mat& slope,
  const arma::mat& xreg_fixed, const arma::mat& xreg_rw, 
  const arma::vec& u, const int distribution, arma::vec weights, 
  const arma::uword n, const arma::uword k, const arma::uword k_rw1, const arma::uword k_rw2,
  const bool response) {
  
  arma::uword k_rw = k_rw1 + k_rw2;
  arma::uword n_iter = weights.n_elem;
  
  arma::cube beta_new(k_rw, n, n_iter);
  arma::cube slope_new(k_rw2, n, n_iter);
  arma::mat y(n, n_iter);
  
  // now our posterior is weighted, so we could again return weighted 
  // predictions.
  // but instead we sample with replacement using weights
  // this decreases the efficiency a bit but there results are simpler to interpret
  
  // as many samples as in our posterior sample
  arma::uvec seq = arma::linspace<arma::uvec>(0, n_iter - 1, n_iter);
  arma::uvec indices = Rcpp::RcppArmadillo::sample(seq, seq.n_elem, true, weights);
  
  for (arma::uword ii = 0; ii < n_iter; ii++) {
    
    unsigned int i = indices(ii); // just laziness
    
    for (arma::uword j = 0; j < k_rw1; j++) {
      beta_new(j, 0, i) = R::rnorm(beta_rw(j, i), sigma_rw1(j, i));
    }
    for (arma::uword j = 0; j < k_rw2; j++) {
      beta_new(k_rw1 + j, 0, i) = beta_rw(k_rw1 + j, i) + slope(j, i);
      slope_new(j, 0, i) = R::rnorm(slope(j, i), sigma_rw2(j, i));
    }
    
    for(arma::uword t = 0; t < n - 1; t++) {
      // linear predictor
      y(t, i) = arma::dot(xreg_rw.col(t), beta_new.slice(i).col(t));
      for (arma::uword j = 0; j < k_rw1; j++) {
        beta_new(j, t + 1, i) = R::rnorm(beta_new(j, t, i), sigma_rw1(j, i));
      }
      for (arma::uword j = 0; j < k_rw2; j++) {
        beta_new(k_rw1 + j, t + 1, i) = beta_new(k_rw1 + j, t, i) + 
          slope_new(j, t, i);
        slope_new(j, t + 1, i) = R::rnorm(slope_new(j, t, i), sigma_rw2(j, i));
      }
    }
    // linear predictor at last time point
    y(n - 1, i) =  arma::dot(xreg_rw.col(n - 1), beta_new.slice(i).col(n - 1));
  }
  
  if (k > 0) {
    for (arma::uword i = 0; i < n_iter; i++) {
      y.col(i) += xreg_fixed * beta_fixed.col(i);
    }
  }
  
  y = arma::exp(y);
  if (response) {
    if(distribution == 1) {
      for (arma::uword i = 0; i < n_iter; i++) {
        for(arma::uword t = 0; t < n; t++) {
          y(t, i) = R::rpois(u(t) * y(t, i));
        }
      }
    } else {
      for (arma::uword i = 0; i < n_iter; i++) {
        for(arma::uword t = 0; t < n; t++) {
          y(t, i) = R::rbinom(u(t),  y(t, i) / (1.0 + y(t, i)));
        }
      }
    }
  } else {
    if(distribution == 1) {
      for (arma::uword i = 0; i < n_iter; i++) {
        for(arma::uword t = 0; t < n; t++) {
          y(t, i) = u(t) * y(t, i);
        }
      }
    } else {
      for (arma::uword i = 0; i < n_iter; i++) {
        for(arma::uword t = 0; t < n; t++) {
          y(t, i) = y(t, i) / (1.0 + y(t, i));
        }
      }
    }
    
  }
  
  
  return Rcpp::List::create(Rcpp::Named("y_new") = y, 
    Rcpp::Named("beta_new") = beta_new, Rcpp::Named("slope_new") = slope_new);
}
