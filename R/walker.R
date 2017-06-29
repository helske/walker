#' Bayesian linear regression with random walk coefficients
#' 
#' Function \code{walker} performs Bayesian inference of a linear 
#' regression model with time-varying, random walk regression coefficients, 
#' i.e. ordinary regression model where instead of constant coefficients the 
#' coefficients follow random walks. All computations are done using Hamiltonian 
#' Monte Carlo provided by Stan, using a state space representation of the model 
#' in order to marginalise over the coefficients for efficient sampling.
#' 
#' @note Beware of overfitting and identifiability issues. For example, one
#' might be tempted to use an intercept and time (1:n) as covariates, but if 
#' both coefficients are allowed to vary over time, this can lead to large uncertainty 
#' in posterior as these to components might not easily distinguishable.
#' 
#' @import rstan Rcpp methods
#' @importFrom stats model.matrix model.response rnorm delete.response terms window ts end glm poisson
#' @rdname walker
#' @useDynLib walker, .registration = TRUE
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
#' @param chains Number of Markov chains. Default is 4.
#' @param init Initial value specification, see \code{\link[rstan]{sampling}}. 
#' Note that compared to default in \code{rstan}, here the default is a to sample from the priors.
#' @param return_x_reg If \code{TRUE}, does not perform sampling, but instead returns the matrix of 
#' predictors after processing the \code{formula}.
#' @param return_y_rep If \code{TRUE} (default), \code{walker} also returns the samples from the 
#' posterior predictive distribution \eqn{p(y_rep | y)}. This argument is ignored if 
#' argument \code{naive} is \code{TRUE}.
#' @param newdata Optional data.frame containing covariates used for prediction. This argument is 
#' ignored if argument \code{naive} is \code{TRUE}.
#' @param ... Further arguments to \code{\link[rstan]{sampling}}.
#' @return A list containing the \code{stanfit} object, observations \code{y},
#'   and covariates \code{xreg} and \code{xreg_new}.
#' @seealso walker_rw2 walker_glm
#' @export
#' @examples 
#' y <- window(log10(UKgas), end = time(UKgas)[100])
#' cos_t <- cos(2 * pi * trend /4)
#' sin_t <- sin(2 * pi * trend /4)
#' dat <- data.frame(y, trend, cos_t, sin_t)
#' future <- length(y) + 1:8
#' new_data <- data.frame(cos_t = cos(2 * pi * future / 4), 
#'   sin_t = sin(2 * pi * future / 4))
#' fit <- walker(y ~ cos_t + sin_t, data = dat, chains = 1, iter = 500, 
#'   newdata = new_data, beta = cbind(0, rep(10, 3)), sigma = cbind(0, rep(10, 4)))
#' print(fit$stanfit, pars = c("sigma_y", "sigma_b"))
#' 
#' plot_coefs(fit)
#' # posterior predictive check:
#' pp_check(fit)
#' 
#' mean_fit <- matrix(summary(fit$stanfit, "beta")$summary[, "mean"], ncol = 4)
#' 
#' # still needs bit manual work..  
#' ts.plot(cbind(y, rowSums(mean_fit * cbind(1, as.matrix(dat[, -1])))),
#'   col = 1:2, lwd = 2:1)
#' intervals <- summary(fit$stanfit, pars = "y_new")$summary[, c("mean", "2.5%", "97.5%")]
#' ts.plot(log10(UKgas), ts(intervals, start = end(y) + c(0,1), frequency = 4),
#'   col = c(1, 2, 2, 2), lty = c(1, 1, 2, 2))
#' 
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
#' kalman_walker <- walker(y ~ x1 + x2, iter = 2000, chains = 1,
#'   beta_prior = cbind(0, rep(2, 3)), sigma_prior = cbind(0, rep(2, 4)))
#' print(kalman_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' betas <- extract(kalman_walker$stanfit, "beta")[[1]]
#' ts.plot(cbind(u, beta1, beta2, apply(betas, 2, colMeans)), 
#'   col = 1:3, lty = rep(2:1, each = 3))
#' sum(get_elapsed_time(kalman_walker$stanfit))
#' naive_walker <- walker(y ~ x1 + x2, iter = 2000, chains = 1, 
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
#' kalman_walker <- walker(y ~ x1 + x2 + x3, iter = 2000, chains = 1,
#'   beta_prior = cbind(0, rep(2, 4)), sigma_prior = cbind(0, rep(2, 5)))
#' print(kalman_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' betas <- extract(kalman_walker$stanfit, "beta")[[1]]
#' ts.plot(cbind(u, beta1, beta2, beta3, apply(betas, 2, colMeans)), 
#'   col = 1:4, lty = rep(2:1, each = 4))
#' sum(get_elapsed_time(kalman_walker$stanfit))
#' # need to increase adapt_delta in order to get rid of divergences
#' # and max_treedepth to get rid of related warnings
#' # and still we end up with low BFMI warning after hours of computation
#' naive_walker <- walker(y ~ x1 + x2 + x3, iter = 2000, chains = 1, 
#'   beta_prior = cbind(0, rep(2, 4)), sigma_prior = cbind(0, rep(2, 5)),
#'   naive = TRUE, control = list(adapt_delta = 0.9, max_treedepth = 15)) 
#' print(naive_walker$stanfit, pars = c("sigma_y", "sigma_b"))
#' sum(get_elapsed_time(naive_walker$stanfit))
#' }
#' 
walker <- function(formula, data, beta_prior, sigma_prior, init, chains, newdata,
  naive = FALSE, return_x_reg = FALSE, return_y_rep = TRUE, ...) {
  
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
  
  
  if (!missing(newdata)) {
    xreg_new <- model.matrix(delete.response(terms(mf)), data = newdata)
    n_new <- nrow(xreg_new)
  } else {
    xreg_new <- matrix(0, 0, k)
    n_new <- 0L
  }
  if (any(is.na(xreg)) || any(is.na(xreg_new))) stop("Missing values in covariates are not allowed.")
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
    sampling(stanmodels$rw_model_naive,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma_y", "sigma_b", "beta"), ...)
  } else {
    sampling(stanmodels$rw_model,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma_y", "sigma_b", "beta", 
        if (return_y_rep) "y_rep", if (n_new > 0) c("y_new", "beta_new")), ...)
  }
  structure(list(stanfit = stanfit, y = y, xreg = xreg, xreg_new = xreg_new), class = "walker_fit")
}
#' Bayesian linear regression with second order random walk coefficients
#' 
#' Function \code{walker_rw2} replaces the first order random walk coefficients 
#' with smooth second order random walks.
#' 
#' @inheritParams walker
#' @seealso walker walker_glm
#' @export
#' @example 
#' 
#' rw1 <- walker(Nile ~ 1, beta_prior = cbind(0, 1000), 
#'   sigma_prior = cbind(0, rep(1000, 2)), chains = 2)
#' rw2 <- walker_rw2(Nile ~ 1, beta_prior = cbind(0, 1000), 
#'   sigma_prior = cbind(0, rep(1000, 2)), chains = 2)
#'   
#' g_y <- geom_point(data = data.frame(y = Nile, x = 1:100), 
#'   aes(x, y, alpha = 0.5), inherit.aes = FALSE) 
#' g_rw1 <- plot_coefs(rw1) + g_y
#' g_rw2 <- plot_coefs(rw2) + g_y
#' if(require("gridExtra"))) {
#'   grid.arrange(g_rw1, g_rw2, ncol=2, top = "RW1 (left) versus RW2 (right)")
#' } else {
#'   g_rw1
#'   g_rw2
#' }
#' 
walker_rw2 <- function(formula, data, beta_prior, sigma_prior, init, chains, newdata,
  naive = FALSE, return_x_reg = FALSE, return_y_rep = TRUE, ...) {
  
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
  
  
  if (!missing(newdata)) {
    xreg_new <- model.matrix(delete.response(terms(mf)), data = newdata)
    n_new <- nrow(xreg_new)
  } else {
    xreg_new <- matrix(0, 0, k)
    n_new <- 0L
  }
  if (any(is.na(xreg)) || any(is.na(xreg_new))) stop("Missing values in covariates are not allowed.")
  if (any(is.na(y))) stop("Missing values in response are not (yet) allowed.")
  
  
  if(!identical(dim(beta_prior), c(k, 2L))) {
    stop("beta_prior should be k x 2 matrix containing columns of prior means and sds for each k coefficients. ")
  }
  if(!identical(dim(sigma_prior), c(k + 1L, 2L))) {
    stop("sigma_prior should be (k + 1) x 2 matrix containing columns of prior means and sds for each k + 1 standard deviations. ")
  }
  stan_data <- list(k = k, n = n, y = y, xreg = t(xreg), 
    n_new = n_new, xreg_new = t(xreg_new), 
    slope_sd = structure(rep(100.0, k), dim = k),
    beta_mean = structure(beta_prior[, 1], dim = k), 
    beta_sd = structure(beta_prior[, 2], dim = k),
    sigma_mean = sigma_prior[, 1], sigma_sd = sigma_prior[, 2])
  
  if (missing(chains)) chains <- 4
  if (missing(init)) {
    init <- replicate(chains, 
      list(sigma_y  = abs(rnorm(1, sigma_prior[1, 1], sigma_prior[1, 2])), 
        sigma_b = structure(abs(rnorm(k, sigma_prior[-1, 1], sigma_prior[-1, 2])), dim = k)), 
      simplify = FALSE)
  }
  stanfit <- 
    sampling(stanmodels$rw2_model,
      data = stan_data, chains = chains, init = init,
      pars = c("sigma_y", "sigma_b", "beta", "slope",
        if (return_y_rep) "y_rep", if (n_new > 0) c("y_new", "beta_new", "slope_new")), ...)
  
  structure(list(stanfit = stanfit, y = y, xreg = xreg, xreg_new = xreg_new), class = "walker_fit")
}


#' Bayesian generalized linear regression with time-varying coefficients
#' 
#' Function \code{walker_glm} is a generalization of \code{walker} for non-Gaussian 
#' models. Compared to \code{walker}, the returned samples are based on Gaussian approximation, 
#' which can be used for exact analysis by weighting the sample properly. These weights 
#' are also returned as a part of the \code{stanfit} (they are generated in the 
#' generated quantities block of Stan model). See details.
#' 
#' 
#' The underlying idea of \code{walker_glm} is based on 
#' Vihola M, Helske J and Franks J (2016), 
#' "Importance sampling type correction of Markov chain Monte Carlo and exact
#' approximations", which is available at ArXiv.
#' 
#' This function is not fully tested yet, so please file and issue and/or pull request 
#' on Github if you encounter problems. The reason there might be problems in some cases 
#' is the use of global approximation (i.e. start of the MCMC) instead of more accurate 
#' but slower local approximation (where model is approximated at each iteration). 
#' However for these restricted models global approximation should be sufficient, 
#' assuming the the initial estimate of the conditional mode of p(xbeta | y, sigma) not too 
#' far away from the true posterior. Therefore by default \code{walker_glm} first finds the 
#' maximum likelihood estimates of the standard deviation parameters 
#' (using \code{\link[KFAS]{KFAS}}) package, and 
#' constructs the approximation at that point, before running the Bayesian 
#' analysis.
#' 
#' @inheritParams walker
#' @importFrom KFAS SSModel SSMregression fitSSM approxSSM
#' @param distribution Currently only Poisson models are supported.
#' @param initial_mode The initial guess of the fitted values on log-scale. 
#' Defines the Gaussian approximation used in the MCMC.
#' Either \code{"obs"} (corresponds to log(y+0.1) in Poisson case), 
#' \code{"glm"} (mode is obtained from time-invariant GLM), \code{"mle"} 
#' (default; mode is obtained from maximum likelihood estimate of the model), 
#' or numeric vector (custom guess).
#' @param u For Poisson model, a vector of exposures i.e. E(y) = u*exp(x*beta). Defaults to 1.
#' @param mc_sim Number of samples used in importance sampling. Default is 50.
#' @return A list containing the \code{stanfit} object, observations \code{y},
#'   and covariates \code{xreg} and \code{xreg_new}.
#' @seealso Package \code{diagis} in CRAN, which provides functions for computing weighted 
#' summary statistics.
#' @export
#' @examples 
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n, 1, 1)
#' beta <- cumsum(c(1, rnorm(n - 1, sd = 0.1)))
#' 
#' level <- -1
#' u <- sample(1:10, size = n, replace = TRUE)
#' y <- rpois(n, u * exp(level + beta * x))
#' ts.plot(y)
#' 
#' out <- walker_glm(y ~ x, u = u, beta_prior = cbind(0, c(10, 10)), 
#'   sigma_prior = cbind(0, c(2, 2)))
#' print(out$stanfit, pars = "sigma_b") ## approximate results
#' library("diagis")
#' weighted_mean(extract(out$stanfit, pars = "sigma_b")$sigma_b, 
#'   extract(out$stanfit, pars = "weights")$weights)
walker_glm <- function(formula, data, beta_prior, sigma_prior, init, chains, newdata, 
  distribution = "poisson", initial_mode = "kfas", u, mc_sim = 50,
  return_x_reg = FALSE,  return_y_rep = TRUE,...) {
  
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
  
  if (!missing(newdata)) .NotYetUsed("newdata", error = TRUE)
  if (!missing(return_y_rep)) .NotYetUsed("return_y_rep", error = TRUE)
  # no predictions supported yet
  # if (!missing(newdata)) {
  #   xreg_new <- model.matrix(delete.response(terms(mf)), data = newdata)
  #   n_new <- nrow(xreg_new)
  # } else {
  xreg_new <- matrix(0, 0, k)
  n_new <- 0L
  #}
  
  if(!identical(dim(beta_prior), c(k, 2L))) {
    stop("beta_prior should be k x 2 matrix containing columns of prior means and sds for each k coefficients. ")
  }
  if(!identical(dim(sigma_prior), c(k, 2L))) {
    stop("sigma_prior should be k x 2 matrix containing columns of prior means and sds for each k standard deviations. ")
  }
  if (missing(u)) {
    u <- rep(1, n)
  }
  if(any(u) <= 0) stop("All values of 'u' must be positive. ")
  
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

