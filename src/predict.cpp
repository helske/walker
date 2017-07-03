#include "predict.h"


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
  const arma::mat& xreg_fixed, const arma::mat& xreg_rw) {

  arma::cube beta_new(xreg_rw.n_rows, xreg_rw.n_cols, sigma_y.n_elem);
  arma::cube slope_new(sigma_rw2.n_rows, xreg_rw.n_cols, sigma_y.n_elem);
  arma::mat y(xreg_rw.n_cols, sigma_y.n_elem);

  unsigned int k_rw1 = sigma_rw1.n_rows;
  unsigned int k_rw2 = sigma_rw2.n_rows;
  
  for (arma::uword i = 0; i <  sigma_y.n_elem; i++) {

    
    for (arma::uword j = 0; j < k_rw1; j++) {
      beta_new(j, 0, i) = R::rnorm(beta_rw(j, i), sigma_rw1(j, i));
    }
    for (arma::uword j = 0; j < k_rw2; j++) {
      beta_new(k_rw1 + j, 0, i) = beta_rw(k_rw1 + j, i) + slope(j, i);
      slope_new(j, 0, i) = R::rnorm(slope(j, i), sigma_rw2(j, i));
    }

    for(arma::uword t = 0; t < xreg_fixed.n_cols - 1; t++) {
      
      y(t, i) = arma::dot(xreg_fixed.col(t), beta_fixed.col(i)) + 
        arma::dot(xreg_rw.col(t), beta_new.slice(i).col(t)) + R::rnorm(0, sigma_y(i));
      
      for (arma::uword j = 0; j < k_rw1; j++) {
        beta_new(j, t + 1, i) = R::rnorm(beta_new(j, t, i), sigma_rw1(j, i));
      }
      for (arma::uword j = 0; j < k_rw2; j++) {
        beta_new(k_rw1 + j, t + 1, i) = beta_new(k_rw1 + j, t, i) + 
          slope_new(j, t, i);
        slope_new(j, t + 1, i) = R::rnorm(slope_new(j, t, i), sigma_rw2(j, i));
      }
    }
    y(xreg_fixed.n_cols - 1, i) =  arma::dot(xreg_fixed.col(xreg_fixed.n_cols - 1), beta_fixed.col(i)) + 
      arma::dot(xreg_rw.col(xreg_fixed.n_cols - 1), beta_new.slice(i).col(xreg_fixed.n_cols - 1));
  }

  return Rcpp::List::create(Rcpp::Named("y_new") = y, 
    Rcpp::Named("beta_new") = beta_new, Rcpp::Named("slope_new") = slope_new);
}
