#' Bayesian regression with random walk coefficients and instrumental variables
#' 
#' Function \code{walker} performs Bayesian inference of a linear 
#' regression model with time-varying, random walk regression coefficients, 
#' i.e. ordinary regression model where instead of constant coefficients the 
#' coefficients follow first or second order random walks. 
#' All Markov chain Monte Carlo computations are done using Hamiltonian 
#' Monte Carlo provided by Stan, using a state space representation of the model 
#' in order to marginalise over the coefficients for efficient sampling.
#' 
#' The \code{rw1} and \code{rw2} functions used in the formula define new formulas 
#' for the first and second order random walks. In addition, these functions 
#' need to be supplied with priors for initial coefficients and the 
#' standard devitions. For second order random walk model, these sigma priors 
#' correspond to the standard deviation of slope distrurbances. For \code{rw2}, 
#' also a prior for the initial slope needs to be defined. See examples.
#' 
#' @note Beware of overfitting and identifiability issues. In particular, 
#' be careful in not defining multiple intercept terms (only one should be present).
#' 
#' @import rstan Rcpp methods
#' @importFrom Rcpp loadModule evalCpp
#' @importFrom stats model.matrix model.response rnorm delete.response terms window ts end glm poisson
#' @rdname walker
#' @useDynLib walker, .registration = TRUE
#' @param formula An object of class \code{\link[stats]{formula}} with additional terms 
#' \code{rw1} and/or \code{rw2} e.g. \code{y ~ x1 + rw1(~ -1 + x2)}. See details.
#' @param data An optional data.frame or object coercible to such, as in \code{\link[stats]{lm}}.
#' @param beta_prior A length vector of length two which defines the 
#' prior mean and standard deviation of the Gaussian prior for time-invariant coefficients
#' @param sigma_y_prior A vector of length two, defining the truncated Gaussian prior for 
#' the observation level standard deviation. Not used in \code{walker_glm}. 
#' @param chains Number of Markov chains. Default is 4.
#' @param init Initial value specification, see \code{\link[rstan]{sampling}}. 
#' Note that compared to default in \code{rstan}, here the default is a to sample from the priors.
#' @param return_x_reg If \code{TRUE}, does not perform sampling, but instead returns the matrix of 
#' predictors after processing the \code{formula}.
#' @param ... Further arguments to \code{\link[rstan]{sampling}}.
#' @return A list containing the \code{stanfit} object, observations \code{y},
#'   and covariates \code{xreg} and \code{xreg_new}.
#' @seealso \code{\link{walker_glm}} for non-Gaussian models.
#' @export
#' @examples 
walker_iv <- function(formula, data, sigma_y_prior, beta_prior, init, chains,
  return_x_reg = FALSE, ...) {
  
  if (missing(data)) data <- environment(formula)
  # Modifying formula object, catching special functions
  mf <- mc <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf[[1L]] <- quote(stats::model.frame)
  mf$na.action <- as.name("na.pass")
  mf$drop.unused.levels <- TRUE
  specials <- c("rw1", "rw2")
  
  all_terms <- terms(formula, specials = specials, data = data)
  rws <- unlist(attr(all_terms, "specials"))
  if (length(rws) > 0) {
    if (length(attr(all_terms, "term.labels")) == length(rws)){
      all_terms <- terms(update.formula(all_terms, . ~ . + .emptyx.),
        specials = specials)
    }
    drops <- which(attr(all_terms, "term.labels") %in%
        rownames(attr(all_terms, "factors"))[rws])
    mf$formula <- formula(drop.terms(all_terms, drops, keep.response = TRUE))
    mf$formula <- update.formula(mf$formula, . ~ . - .emptyx., simplify = TRUE)
  }
  
  # build y and xreg
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  n <- length(y)
  xreg_fixed <- model.matrix(attr(mf, "terms"), mf)
  
  ## RWs
  
  vars <- attr(all_terms, "variables")
  if (!is.null(attr(all_terms, "specials")$rw1)) {
    comp <- vars[[1 + attr(all_terms, "specials")$rw1[1]]]
    rw1_out <- eval(comp, envir = data, enclos = parent.frame())
    # only intercept
    if (nrow(rw1_out$xreg) == 0) rw1_out$xreg <- matrix(1, n, 1)
  } else {
    rw1_out <- list(xreg = matrix(0, n, 0), 
      beta_prior = numeric(2), sigma_prior = numeric(2))
  }
  if (!is.null(attr(all_terms, "specials")$rw2)) {
    comp <- vars[[1 + attr(all_terms, "specials")$rw2[1]]]
    rw2_out <- eval(comp, envir = data, enclos = parent.frame())
    # only intercept
    if (nrow(rw2_out$xreg) == 0) rw2_out$xreg <- matrix(1, n, 1)
  } else {
    rw2_out <- list(xreg = matrix(0, n, 0), 
      beta_prior = numeric(2), sigma_prior = numeric(2), slope_prior = numeric(2))
  }
  
  xreg_rw <- cbind(rw1_out$xreg, rw2_out$xreg)
  
  k_fixed <- max(0, ncol(xreg_fixed))
  k_rw1 <- max(0, ncol(rw1_out$xreg))
  k_rw2 <- max(0, ncol(rw2_out$xreg))
  if (return_x_reg) return(list(xreg_fixed = xreg_fixed, xreg_rw = xreg_rw))
  
  if (any(is.na(xreg_fixed)) || any(is.na(xreg_rw))) stop("Missing values in covariates are not allowed.")
  if (any(is.na(y))) stop("Missing values in response are not (yet) allowed.")
  
  if(k_fixed > 0 && length(beta_prior) != 2) {
    stop("beta_prior should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of fixed coefficients. ")
  }
  if(length(sigma_y_prior) != 2) {
    stop("sigma_prior should be should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of the standard deviation of y. ")
  }
  
  stan_data <- list(
    k_fixed = k_fixed, 
    k_rw1 = k_rw1,
    k_rw2 = k_rw2,
    m = k_rw1 + 2 * k_rw2,
    k = k_rw1 + k_rw2,
    n = n, 
    xreg_fixed = xreg_fixed, 
    xreg_rw = t(xreg_rw), 
    y = y, 
    sigma_y_mean = sigma_y_prior[1],
    sigma_y_sd = sigma_y_prior[2],
    beta_fixed_mean = if (k_fixed > 0) beta_prior[1] else 0,
    beta_fixed_sd = if (k_fixed > 0) beta_prior[2] else 0,
    beta_rw1_mean = rw1_out$beta_prior[1],
    beta_rw1_sd = rw1_out$beta_prior[2],
    beta_rw2_mean = rw2_out$beta_prior[1],
    beta_rw2_sd = rw2_out$beta_prior[2],
    sigma_rw1_mean = rw1_out$sigma_prior[1],
    sigma_rw1_sd = rw1_out$sigma_prior[2],
    sigma_rw2_mean = rw2_out$sigma_prior[1],
    sigma_rw2_sd = rw2_out$sigma_prior[2],
    slope_mean = rw2_out$slope_prior[1],
    slope_sd = rw2_out$slope_prior[2]
  )
  
  if (missing(chains)) chains <- 4
  if (missing(init)) {
    init <- replicate(chains, 
      list(beta_fixed = 
          if (k_fixed > 0) {
            structure(rnorm(k_fixed, beta_prior[1], beta_prior[2] / 10), dim = k_fixed)
          } else {
            structure(numeric(0), dim = 0) 
          },
        sigma_y  = abs(rnorm(1, sigma_y_prior[1], sigma_y_prior[2] / 10)), 
        sigma_rw1 = 
          if (k_rw1 > 0) {
            structure(abs(rnorm(k_rw1, rw1_out$sigma_prior[1], rw1_out$sigma_prior[2] / 10)), dim = k_rw1) 
          } else {
            structure(numeric(0), dim = 0)
          }, 
        sigma_rw2 = 
          if (k_rw2 > 0) {
            structure(abs(rnorm(k_rw2, rw2_out$sigma_prior[1], rw2_out$sigma_prior[2] / 10)), dim = k_rw2) 
          } else {
            structure(numeric(0), dim = 0)
          }), 
      simplify = FALSE)
  }
  stanfit <- sampling(stanmodels$walker_lm,
    data = stan_data, chains = chains, init = init,
    pars = c("sigma_y", "sigma_rw1", "sigma_rw2", "beta_fixed", "beta_rw", 
      "slope", "y_fit", "y_rep"), ...)
  
  structure(list(stanfit = stanfit, y = y, xreg_fixed = xreg_fixed, 
    xreg_rw = xreg_rw, call = mc, distribution = "gaussian"), class = "walker_fit")
}
