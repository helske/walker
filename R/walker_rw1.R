#' Comparison of naive and state space implementation of RW1 regression model
#' 
#' This function is the first iteration of the function \code{walker},
#' which supports only time-varying model where all coefficients ~ rw1.
#' This is kept as part of the package in order to compare "naive" and 
#' state space versions of the model in the vignette.
#' 
#' @export
#' @param formula An object of class \code{\link[stats]{formula}}. See \code{\link[stats]{lm}} for details.
#' @param data An optional data.frame or object coercible to such, as in \code{\link[stats]{lm}}.
#' @param beta_prior A matrix with \eqn{k} rows and 2 columns, where first columns defines the 
#' prior means of the Gaussian priors of the corresponding \eqn{k} regression coefficients, 
#' and the second column defines the the standard deviations of those prior distributions.
#' @param sigma_prior A matrix with \eqn{k + 1} rows and two colums with similar structure as 
#' \code{beta_prior}, with first row corresponding to the prior of the standard deviation of the 
#' observation level noise, and rest of the rows define the priors for the standard deviations of 
#' random walk noise terms. The prior distributions for all sigmas are 
#' Gaussians truncated to positive real axis. For non-Gaussian models, this should contain only k rows. 
#' For second order random walk model, these priors correspond to the slope level standard deviations.
#' @param naive Only used for \code{walker} function. 
#' If \code{TRUE}, use "standard" approach which samples the joint posterior 
#' \eqn{p(beta, sigma | y)}. If \code{FALSE} (the default), use marginalisation approach 
#' where we sample the marginal posterior \eqn{p(sigma | y)} and generate the samples of 
#' \eqn{p(beta | sigma, y)} using state space modelling techniques 
#' (namely simulation smoother by Durbin and Koopman (2002)). Both methods give asymptotically 
#' identical results, but the latter approach is computationally much more efficient.
#' @param return_x_reg If \code{TRUE}, does not perform sampling, but instead returns the matrix of 
#' predictors after processing the \code{formula}.
#' @param chains Number of Markov chains. Default is 4.
#' @param init Initial value specification, see \code{\link[rstan]{sampling}}. 
#' @param ... Additional arguments to \code{\link[rstan]{sampling}}. 
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
#' kalman_walker <- walker_rw1(y ~ -1 + 
#'   rw1(~ x1 + x2, beta_prior = c(0, 2), sigma_prior = c(0, 2)), 
#'   sigma_y_prior = c(0, 2), iter = 2000, chains = 1)
#' print(kalman_walker$stanfit, pars = c("sigma_y", "sigma_rw1"))
#' betas <- extract(kalman_walker$stanfit, "beta")[[1]]
#' ts.plot(cbind(u, beta1, beta2, apply(betas, 2, colMeans)), 
#'   col = 1:3, lty = rep(2:1, each = 3))
#' sum(get_elapsed_time(kalman_walker$stanfit))
#' naive_walker <- walker_rw1(y ~ x1 + x2, iter = 2000, chains = 1, 
#'   beta_prior = cbind(0, rep(2, 3)), sigma_prior = cbind(0, rep(2, 4)), 
#'   naive = TRUE)
#' print(naive_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' sum(get_elapsed_time(naive_walker$stanfit))
#' 
#' ## Larger problem, this takes some time with naive approach
#'
#' set.seed(123)
#' n <- 500
#' beta1 <- cumsum(c(1.5, rnorm(n - 1, 0, sd = 0.05)))
#' beta2 <- cumsum(c(-1, rnorm(n - 1, 0, sd = 0.5)))
#' beta3 <- cumsum(c(-1.5, rnorm(n - 1, 0, sd = 0.15)))
#' beta4 <- 2
#' x1 <- rnorm(n, 1)
#' x2 <- 0.25 * cos(1:n)
#' x3 <- runif(n, 1, 3)
#' ts.plot(cbind(beta1 * x1, beta2 * x2, beta3 * x3), col = 1:3)
#' a <- cumsum(rnorm(n))
#' signal <- a + beta1 * x1 + beta2 * x2 + beta3 * x3
#' y <- rnorm(n, signal)
#' ts.plot(y)
#' lines(signal, col = 2)
#' kalman_walker <- walker_rw1(y ~ x1 + x2 + x3, iter = 2000, chains = 1,
#'   beta_prior = cbind(0, rep(2, 4)), sigma_prior = cbind(0, rep(2, 5)))
#' print(kalman_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' betas <- extract(kalman_walker$stanfit, "beta")[[1]]
#' ts.plot(cbind(u, beta1, beta2, beta3, apply(betas, 2, colMeans)), 
#'   col = 1:4, lty = rep(2:1, each = 4))
#' sum(get_elapsed_time(kalman_walker$stanfit))
#' # need to increase adapt_delta in order to get rid of divergences
#' # and max_treedepth to get rid of related warnings
#' # and still we end up with low BFMI warning after hours of computation
#' naive_walker <- walker_rw1(y ~ x1 + x2 + x3, iter = 2000, chains = 1, 
#'   beta_prior = cbind(0, rep(2, 4)), sigma_prior = cbind(0, rep(2, 5)),
#'   naive = TRUE, control = list(adapt_delta = 0.9, max_treedepth = 15)) 
#' print(naive_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' sum(get_elapsed_time(naive_walker$stanfit))
#' }
#' 
walker_rw1 <- function(formula, data, beta_prior, sigma_prior, init, chains,
  naive = FALSE, return_x_reg = FALSE,  ...) {
  
  # build y and xreg
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  n <- length(y)
  xreg <- model.matrix(attr(mf, "terms"), mf)
  if (return_x_reg) return(xreg)
  k <- ncol(xreg)
  
    xreg_new <- matrix(0, 0, k)
    n_new <- 0L
  
  if (any(is.na(xreg))) stop("Missing values in covariates are not allowed.")
  if (any(is.na(y))) stop("Missing values in response are not (yet) allowed.")
  
  
  if(!identical(dim(beta_prior), c(k, 2L))) {
    stop("beta_prior should be k x 2 matrix containing columns of prior means and sds for each k coefficients. ")
  }
  if(!identical(dim(sigma_prior), c(k + 1L, 2L))) {
    stop("sigma_prior should be (k + 1) x 2 matrix containing columns of prior means and sds for each k + 1 standard deviations. ")
  }
  stan_data <- list(k = k, n = n, y = y, xreg = t(xreg), 
    n_new = n_new, xreg_new = t(xreg_new),
    beta_mean = structure(beta_prior[, 1], dim = k), 
    beta_sd = structure(beta_prior[, 2], dim = k),
    sigma_mean = sigma_prior[, 1], sigma_sd = sigma_prior[, 2])
  
  if (missing(chains)) chains <- 4
  if (missing(init)) {
    init <- replicate(chains, 
      list(sigma_y  = abs(rnorm(1, sigma_prior[1, 1], sigma_prior[1, 2])), 
        sigma_b = structure(abs(rnorm(k, sigma_prior[-1, 1], sigma_prior[-1, 2])), dim = k), 
        beta = structure(rnorm(k, beta_prior[, 1], beta_prior[, 2]), dim = k)), simplify = FALSE)
  }
  stanfit <- if (naive) {
    sampling(stanmodels$rw1_model_naive,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma_y", "sigma_b", "beta"), ...)
  } else {
    sampling(stanmodels$rw1_model,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma_y", "sigma_b", "beta", 
        "y_rep", if (n_new > 0) c("y_new", "beta_new")), ...)
  }
  structure(list(stanfit = stanfit, y = y, xreg = xreg, xreg_new = xreg_new), class = "walker_fit")
}

walker_rw1_glm <- function(formula, data, beta_prior, sigma_prior, init, chains,
  distribution = "poisson", initial_mode = "kfas", u, mc_sim = 50,
  return_x_reg = FALSE, ...) {
  
  distribution <- match.arg(distribution, choices = "poisson")
  
  # build y and xreg
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  n <- length(y)
  xreg <- model.matrix(attr(mf, "terms"), mf)
  if (return_x_reg) return(xreg)
  k <- ncol(xreg)
  
  if (any(is.na(xreg))) stop("Missing values in covariates are not allowed.")
  if (any(is.na(y))) stop("Missing values in response are not (yet) allowed.")
  
  xreg_new <- matrix(0, 0, k)
  n_new <- 0L
  
  if(!identical(dim(beta_prior), c(k, 2L))) {
    stop("beta_prior should be k x 2 matrix containing columns of prior means and sds for each k coefficients. ")
  }
  if(!identical(dim(sigma_prior), c(k, 2L))) {
    stop("sigma_prior should be k x 2 matrix containing columns of prior means and sds for each k standard deviations. ")
  }
  if (missing(u)) {
    u <- rep(1, n)
  }
  if(any(u <= 0)) stop("All values of 'u' must be positive. ")
  
  if (is.numeric(initial_mode)) {
    pseudo_H <- 1 / (u * exp(initial_mode))
    pseudo_y <- y * pseudo_H + initial_mode - 1
  } else {
    switch(initial_mode, 
      obs = {
        
        expmode <- y / u + 0.1
        pseudo_H <- 1 / (u * expmode)
        pseudo_y <- y * pseudo_H + log(expmode) - 1
        
      },
      glm = {
        
        if (missing(data)) {
          fit <- glm(formula, offset = log(u), family = poisson)
        } else {
          fit <- glm(formula, offset = log(u), data = data, family = poisson)
        }
        pseudo_H <- 1 / fit$fitted.values
        pseudo_y <- y * pseudo_H + fit$linear.predictors - log(u) - 1
        
      },
      kfas = {
        
        model <- SSModel(y ~ -1 +
            SSMregression(formula, Q = diag(NA, k), remove.intercept = FALSE, 
              a1 = beta_prior[,1],
              P1 = diag(beta_prior[,2]^2, k)),
          distribution = "poisson", u = u)
        fit <- fitSSM(model, inits = rep(-1, k), method = "BFGS")
        app <- approxSSM(fit$model)
        pseudo_H <- as.numeric(app$H)
        pseudo_y <- as.numeric(app$y)
        
      },
      stop("Argument 'initial_mode' should be either 'obs', 'glm', 'kfas', or a numeric vector.")
    )
  }
  
  
  stan_data <- list(k = k, n = n, y = pseudo_y, Ht = pseudo_H, 
    y_original = y, u = u, distribution = 1L, N = mc_sim, xreg = t(xreg), 
    n_new = n_new, xreg_new = t(xreg_new),
    beta_mean = structure(beta_prior[, 1], dim = k), 
    beta_sd = structure(beta_prior[, 2], dim = k),
    sigma_mean = structure(sigma_prior[, 1], dim = k), 
    sigma_sd = structure(sigma_prior[, 2], dim = k))
  
  if (missing(chains)) chains <- 4
  if (missing(init)) {
    init <- replicate(chains, 
      list(
        sigma_b = structure(abs(rnorm(k, sigma_prior[, 1], sigma_prior[, 2])), dim = k), 
        beta = structure(rnorm(k, beta_prior[, 1], beta_prior[, 2]), dim = k)), 
      simplify = FALSE)
  }
  
  sampling(stanmodels$rw_glm_model,
    data = stan_data, chains = chains, init = init,
    pars = c("sigma_b", "beta", "weights"), ...)
}
