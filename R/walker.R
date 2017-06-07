#' Fully Bayesian linear regression with time-varying coefficients
#' 
#' Function \code{walker} performs Bayesian inference of a linear 
#' regression model with time-varying regression coefficients, 
#' i.e. ordinary regression model where instead of constant coefficients the 
#' coefficients follow random walks. All computations are done using Hamiltonian 
#' Monte Carlo provided by Stan, using a state space representation of the model 
#' in order to marginalise over the coefficients for efficient sampling.
#' 
#' @importFrom rstan sampling
#' @rdname walker
#' @useDynLib walker, .registration = TRUE
#' @param formula An object of class \link{\code{"formula"}}. See \link{code{"lm"}} for details.
#' @param data An optional data.frame or object coercible to such, as in \link{code{"lm"}}.
#' @param beta_prior A matrix with \eqn{k} rows and 2 columns, where first columns defines the 
#' prior means of the Gaussian priors of the corresponding \eqn{k} regression coefficients, 
#' and the second column defines the the standard deviations of those prior distributions.
#' @param sigma_prior A matrix with \eqn{k + 1} rows and two colums with similar structure as 
#' \code{beta_prior}, with first row corresponding to the prior of the standard deviation of the 
#' observation level noise, and rest of the rows define the priors for the standard deviations of 
#' random walk noise terms. The prior distributions for all sigmas are 
#' Gaussians truncated to positive real axis.
#' @param naive If \code{TRUE}, use "standard" approach which samples the joint posterior 
#' \eqn{p(beta, sigma | y)}. If \code{FALSE} (the default), use marginalisation approach 
#' where we sample the marginal posterior \eqn{p(sigma | y)} and generate the samples of 
#' \eqn{p(beta | sigma, y)} using state space modelling techniques 
#' (namely simulation smoother by Durbin and Koopman (2002)). Both methods give asymptotically 
#' identical results, but the latter approach is computationally much more efficient.
#' @export
#' @examples 
#' \dontrun{
#' ## Comparing the approaches, note that with such a small data 
#' ## the differences aren't huge, but try same with n = 500 and/or more terms...
#' set.seed(123)
#' n <- 100
#' beta1 <- cumsum(c(0.5, rnorm(n - 1, 0, sd = 0.05)))
#' beta2 <- cumsum(c(-1, rnorm(n - 1, 0, sd = 0.15)))
#' x1 <- rnorm(n, 1)
#' x2 <- 0.25 * cos(1:n)
#' ts.plot(cbind(beta1 * x1, beta2 *x2), col = 1:2)
#' u <- cumsum(rnorm(n))
#' y <- rnorm(n, u + beta1 * x1 + beta2 * x2)
#' ts.plot(y)
#' lines(u + beta1 * x1 + beta2 * x2, col = 2)
#' kalman_walker <- walker(y ~ x1 + x2, iter = 2000, chains = 1, seed = 1,
#'   beta_prior = cbind(0, rep(2, 3)), sigma_prior = cbind(0, rep(2, 4)))
#' print(kalman_walker, pars = c("sigma", "beta[1,1]", "beta[1,100]", 
#'   "beta[2,1]", "beta[2,100]", "beta[3,1]", "beta[3,100]"))
#' sum(get_elapsed_time(kalman_walker))
#' 
#' naive_walker <- walker(y ~ x1 + x2, iter = 2000, chains = 1, seed = 1, 
#'   beta_prior = cbind(0, rep(2, 3)), sigma_prior = cbind(0, rep(2, 4)), 
#'   naive = TRUE)
#' print(naive_walker, pars = c("sigma", "beta[1,1]", "beta[1,100]", 
#'   "beta[2,1]", "beta[2,100]", "beta[3,1]", "beta[3,100]"))
#' # check rstan:::throw_sampler_warnings(kalman_walker) 
#' # (does not work automatically for single chain)
#' sum(get_elapsed_time(naive_walker))
#' 
#' ## Larger problem
#'
#' set.seed(123)
#' n <- 500
#' beta1 <- cumsum(c(1.5, rnorm(n - 1, 0, sd = 0.05)))
#' beta2 <- cumsum(c(-1, rnorm(n - 1, 0, sd = 0.5)))
#' beta3 <- cumsum(c(-1.5, rnorm(n - 1, 0, sd = 0.15)))
#' beta4 <- 2
#' x1 <- rnorm(n, 1)
#' x2 <- 0.25 * cos(1:n)
#' x3 <- 1:n/10
#' x4 <- runif(n, 1, 3)
#' ts.plot(cbind(beta1 * x1, beta2 * x2, beta3 * x3, beta4 * x4), col = 1:4)
#' u <- cumsum(rnorm(n))
#' signal <- u + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4
#' y <- rnorm(n, signal)
#' ts.plot(y)
#' lines(signal, col = 2)
#' kalman_walker <- walker(y ~ x1 + x2 + x3 + x4, iter = 2000, chains = 1, seed = 1,
#'   beta_prior = cbind(0, rep(2, 5)), sigma_prior = cbind(0, rep(2, 6)))
#' print(kalman_walker, pars = "sigma")
#' sum(get_elapsed_time(kalman_walker))
#' # need to increase adapt_delta in order to get rid of divergences
#' naive_walker <- walker(y ~ x1 + x2 + x3 + x4, iter = 2000, chains = 1, seed = 1,
#'   beta_prior = cbind(0, rep(2, 5)), sigma_prior = cbind(0, rep(2, 6)),
#'   naive = TRUE, control = list(adapt_delta = 0.9)) 
#' print(naive_walker, pars = c("sigma"))
#' # check rstan:::throw_sampler_warnings(naive_walker)
#' # (does not work automatically for single chain)
#' sum(get_elapsed_time(naive_walker))
#' }

walker <- function(formula, data, beta_prior, sigma_prior, init, chains,
  naive = FALSE, return_x_reg = FALSE, ...) {
  
  # build y and xreg
  
  if (missing(data)) {
    data <- environment(formula)
    period <- NULL
  } else {
    period <- tsp(data)[3]
    data <- as.data.frame(data)
  }
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$na.action <- as.name("na.pass")
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  n <- length(y)
  xreg <- model.matrix(attr(mf, "terms"), mf)
  k <- ncol(xreg)
  
  if(!identical(dim(beta_prior), c(k, 2L))) {
    stop("beta_prior should be k x 2 matrix containing columns of prior means and sds for each k coefficients. ")
  }
  if(!identical(dim(sigma_prior), c(k + 1L, 2L))) {
    stop("sigma_prior should be (k + 1) x 2 matrix containing columns of prior means and sds for each k + 1 standard deviations. ")
  }
  stan_data <- list(k = k, n = n, y = y, xreg = xreg,
    beta_mean = beta_prior[, 1], beta_sd = beta_prior[, 2], 
    sigma_mean = sigma_prior[, 1], sigma_sd = sigma_prior[, 2])
  
  if (missing(chains)) chains <- 4
  if (missing(init)) {
    init <- replicate(chains, 
      list(sigma = abs(rnorm(k + 1, sigma_prior[, 1], sigma_prior[, 2])), 
        beta = rnorm(k, beta_prior[, 1], beta_prior[, 2])), simplify = FALSE)
  }
  if (naive) {
    sampling(walker:::stanmodels$rw_model_naive,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma", "beta"), ...)
  } else {
  sampling(walker:::stanmodels$rw_model,
    data = stan_data, chains = chains, init = init,
    pars = c("sigma", "beta"), ...)
  }
}

