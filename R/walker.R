#' Bayesian regression with random walk coefficients
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
#' 
#' ## note very low number of iterations for the CRAN checks
#' 
#' rw1_fit <- walker(Nile ~ -1 + 
#'   rw1(~ 1, 
#'     beta_prior = c(1000, 100), 
#'     sigma_prior = c(0, 100)), 
#'   sigma_y_prior = c(0, 100), 
#'   iter = 50, chains = 1)
#'   
#' rw2_fit <- walker(Nile ~ -1 + 
#'   rw2(~ 1,
#'     beta_prior = c(1000, 100), 
#'     sigma_prior = c(0, 100), 
#'     slope_prior = c(0, 100)), 
#'   sigma_y_prior = c(0, 100), 
#'   iter = 50, chains = 1)
#'   
#' g_y <- geom_point(data = data.frame(y = Nile, x = time(Nile)), 
#'   aes(x, y, alpha = 0.5), inherit.aes = FALSE) 
#' g_rw1 <- plot_coefs(rw1_fit) + g_y
#' g_rw2 <- plot_coefs(rw2_fit) + g_y
#' if(require("gridExtra")) {
#'   grid.arrange(g_rw1, g_rw2, ncol=2, top = "RW1 (left) versus RW2 (right)")
#' } else {
#'   g_rw1
#'   g_rw2
#' }
#' \dontrun{
#' y <- window(log10(UKgas), end = time(UKgas)[100])
#' n <- 100
#' cos_t <- cos(2 * pi * 1:n / 4)
#' sin_t <- sin(2 * pi * 1:n / 4)
#' dat <- data.frame(y, cos_t, sin_t)
#' fit <- walker(y ~ -1 + 
#'   rw1(~ cos_t + sin_t, beta_prior = c(0, 10), sigma_prior = c(0, 2)), 
#'   sigma_y_prior = c(0, 10), data = dat, chains = 1, iter = 250)
#' print(fit$stanfit, pars = c("sigma_y", "sigma_rw1"))
#' 
#' plot_coefs(fit)
#' # posterior predictive check:
#' pp_check(fit)
#' 
#' newdata <- data.frame(
#'   cos_t = cos(2 * pi * 101:108 / 4), 
#'   sin_t = sin(2 * pi * 101:108 / 4))
#' pred <- predict(fit, newdata)
#' plot_predict(pred)
#' }
#' 
walker <- function(formula, data, sigma_y_prior, beta_prior, init, chains,
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

#' Bayesian generalized linear regression with time-varying coefficients
#' 
#' Function \code{walker_glm} is a generalization of \code{walker} for non-Gaussian 
#' models. Compared to \code{walker}, the returned samples are based on Gaussian approximation, 
#' which can then be used for exact-approximate analysis by weighting the sample properly. These weights 
#' are also returned as a part of the \code{stanfit} (they are generated in the 
#' generated quantities block of Stan model). Note that plotting functions \code{pp_check}, 
#' \code{plot_coefs}, and \code{plot_predict} resample the posterior based on weights 
#' before plotting, leading to "exact" analysis.
#' 
#' The underlying idea of \code{walker_glm} is based on 
#' Vihola M, Helske J and Franks J (2016), 
#' "Importance sampling type correction of Markov chain Monte Carlo and exact
#' approximations", which is available at ArXiv.
#' 
#' This function is not fully tested yet, so please file an issue and/or pull request 
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
#' @importFrom KFAS SSModel SSMcustom fitSSM approxSSM
#' @param distribution Either \code{"poisson"} or \code{"binomial"}.
#' @param initial_mode The initial guess of the fitted values on log-scale. 
#' Defines the Gaussian approximation used in the MCMC.
#' Either \code{"obs"} (corresponds to log(y+0.1) in Poisson case), 
#' \code{"glm"} (mode is obtained from time-invariant GLM), \code{"mle"} 
#' (default; mode is obtained from maximum likelihood estimate of the model), 
#' or numeric vector (custom guess).
#' @param u For Poisson model, a vector of exposures i.e. E(y) = u*exp(x*beta). 
#' For binomial, a vector containing the number of trials. Defaults 1.
#' @param mc_sim Number of samples used in importance sampling. Default is 50.
#' @return A list containing the \code{stanfit} object, observations \code{y},
#'   covariates \code{xreg_fixed}, and \code{xreg_rw}.
#' @seealso Package \code{diagis} in CRAN, which provides functions for computing weighted 
#' summary statistics.
#' @export
#' @examples 
#' 
#' ## note very low number of iterations for the CRAN checks
#' 
#' data("discoveries", package = "datasets")
#' out <- walker_glm(discoveries ~ -1 + 
#'   rw2(~ 1, beta_prior = c(0, 10), sigma_prior = c(0, 2), slope_prior = c(0, 2)), 
#'   distribution = "poisson", iter = 50, chains = 1, refresh = 0)
#' 
#' plot_fit(out)
#' 
#' \dontrun{
#' 
#' set.seed(1)
#' n <- 25
#' x <- rnorm(n, 1, 1)
#' beta <- cumsum(c(1, rnorm(n - 1, sd = 0.1)))
#' 
#' level <- -1
#' u <- sample(1:10, size = n, replace = TRUE)
#' y <- rpois(n, u * exp(level + beta * x))
#' ts.plot(y)
#' 
#' out <- walker_glm(y ~ -1 + rw1(~ x, beta_prior = c(0, 10), 
#'   sigma_prior = c(0, 10)), distribution = "poisson", 
#'   iter = 250, chains = 1, refresh = 0)
#' print(out$stanfit, pars = "sigma_rw1") ## approximate results
#' if (require("diagis")) {
#'   weighted_mean(extract(out$stanfit, pars = "sigma_rw1")$sigma_rw1, 
#'     extract(out$stanfit, pars = "weights")$weights)
#' }
#' plot_coefs(out)
#' pp_check(out)
#' 
#' }
#'              
walker_glm <- function(formula, data, beta_prior, init, chains,
  return_x_reg = FALSE, distribution ,
  initial_mode = "kfas", u, mc_sim = 50, ...) {
  
  distribution <- match.arg(distribution, choices = c("poisson", "binomial"))
  
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
 
  if (missing(u)) {
    u <- rep(1, n)
  }
  if(any(u <= 0)) stop("All values of 'u' must be positive. ")
  
  beta_fixed_mean = if (k_fixed > 0) beta_prior[1] else 0
  beta_fixed_sd = if (k_fixed > 0) beta_prior[2] else 0
  beta_rw1_mean = rw1_out$beta_prior[1]
  beta_rw1_sd = rw1_out$beta_prior[2]
  beta_rw2_mean = rw2_out$beta_prior[1]
  beta_rw2_sd = rw2_out$beta_prior[2]
  slope_mean = rw2_out$slope_prior[1]
  slope_sd = rw2_out$slope_prior[2]
  
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
        
        fit <- glm(y ~ ., data = data.frame(cbind(xreg_fixed, xreg_rw)), offset = log(u), family = poisson)
        pseudo_H <- 1 / fit$fitted.values
        pseudo_y <- y * pseudo_H + fit$linear.predictors - log(u) - 1
        
      },
      kfas = {
        m <- k_fixed + k_rw1 + 2 * k_rw2
        Zt <- array(0, dim = c(1, m, n))
        if (k_fixed > 0) {
          Zt[1, 1:k_fixed, ] <- t(xreg_fixed)
        }
        Zt[1, (k_fixed + 1):(k_fixed + k_rw1 + k_rw2),] <- t(xreg_rw)
        Tt <- Rt <- diag(m)
        if(k_rw2 > 0) {
          Tt[(k_fixed + k_rw1 + 1):(k_fixed + k_rw1 + k_rw2), 
            (k_fixed + k_rw1 + k_rw2 + 1):m] <- diag(k_rw2)
        }
        Qt <- diag(rep(c(0, NA, 0, NA), times = c(k_fixed, k_rw1, k_rw2, k_rw2)), m)
        a1 <- rep(c(beta_fixed_mean, beta_rw1_mean, beta_rw2_mean, slope_mean), 
          times = c(k_fixed, k_rw1, k_rw2, k_rw2))
        P1 <- diag(rep(c(beta_fixed_sd, beta_rw1_sd, beta_rw2_sd, slope_sd), 
          times = c(k_fixed, k_rw1, k_rw2, k_rw2)), m)
        P1inf <- diag(0, m)
        model <- SSModel(y ~ -1 + SSMcustom(Zt, Tt, Rt, Qt, a1, P1, P1inf),
          distribution = distribution, u = u)
        fit <- fitSSM(model, inits = rep(-1, k_rw1 + k_rw2), method = "BFGS")
        app <- approxSSM(fit$model)
        pseudo_H <- as.numeric(app$H)
        pseudo_y <- as.numeric(app$y)
        
      },
      stop("Argument 'initial_mode' should be either 'obs', 'glm', 'kfas', or a numeric vector.")
    )
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
    beta_fixed_mean = beta_fixed_mean,
    beta_fixed_sd = beta_fixed_sd,
    beta_rw1_mean = beta_rw1_mean,
    beta_rw1_sd = beta_rw1_sd,
    beta_rw2_mean = beta_rw2_mean,
    beta_rw2_sd = beta_rw2_sd,
    sigma_rw1_mean = rw1_out$sigma_prior[1],
    sigma_rw1_sd = rw1_out$sigma_prior[2],
    sigma_rw2_mean = rw2_out$sigma_prior[1],
    sigma_rw2_sd = rw2_out$sigma_prior[2],
    slope_mean = slope_mean,
    slope_sd = slope_sd,
    y = pseudo_y, 
    Ht = pseudo_H, 
    y_original = y, 
    u = as.integer(u), 
    distribution = pmatch(distribution, c("poisson", "binomial")), 
    N = mc_sim
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
  
  stanfit <- sampling(stanmodels$walker_glm,
    data = stan_data, chains = chains, init = init,
    pars = c("sigma_rw1", "sigma_rw2", "beta_fixed", "beta_rw", 
      "slope", "y_fit", "y_rep", "weights"), ...)
  
  structure(list(stanfit = stanfit, y = y, xreg_fixed = xreg_fixed, 
    xreg_rw = xreg_rw, u = u, distribution = distribution, call = mc), 
    class = "walker_fit")
}
