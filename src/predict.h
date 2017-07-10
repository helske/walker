#ifndef PREDICT_WALKER_H
#define PREDICT_WALKER_H

#include <RcppArmadillo.h>
// sigma_rw: k x N matrix of posterior samples of rw_sigmas
// sigma_y: vector of length N of posterior samples of sigma_y
// beta_fixed p x N matrix of posterior samples of fixed betas
// beta_rw: k x N matrix of posterior samples of beta_rws at time n
// xreg: k x n_new matrix of new covariates
Rcpp::List predict_walker(const arma::mat& sigma_rw1, 
  const arma::mat& sigma_rw2, const arma::vec sigma_y,
  const arma::mat beta_fixed, const arma::mat& beta_rw, const arma::mat& slope,
  const arma::mat& xreg_fixed, const arma::mat& xreg_rw);

Rcpp::List predict_walker_glm(const arma::mat& sigma_rw1, 
  const arma::mat& sigma_rw2,
  const arma::mat beta_fixed, const arma::mat& beta_rw, const arma::mat& slope,
  const arma::mat& xreg_fixed, const arma::mat& xreg_rw, 
  const arma::vec& u, const int distribution, arma::vec weights);

#endif