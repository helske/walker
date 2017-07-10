#' Construct a first-order random walk component 
#' 
#' Auxiliary function used inside of the formula of \code{walker}.
#' 
#' @export
#' @param formula Formula for RW1 part of the model. Only right-hand-side is used. 
#' @param data Optional data.frame.
#' @param beta_prior A length vector of length two which defines the 
#' prior mean and standard deviation of the Gaussian prior for coefficients at time 1.
#' @param sigma_prior A vector of length two, defining the truncated Gaussian prior for 
#' the coefficient level standard deviation. 
rw1 <- function(formula, data, beta_prior, sigma_prior) {
 
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- as.name("na.pass")
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  xreg <- model.matrix(attr(mf, "terms"), mf)
  
  if(length(beta_prior) != 2) {
    stop("beta_prior should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of coefficients. ")
  }
  if(length(sigma_prior) != 2) {
    stop("sigma_prior should be should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of standard deviations. ")
  }
  list(xreg = xreg, beta_prior = beta_prior, 
    sigma_prior = sigma_prior)
  
}
#' Construct a first-order random walk component 
#' 
#' Auxiliary function used inside of the formula of \code{walker}.
#' 
#' @export
#' @param formula Formula for RW1 part of the model. Only right-hand-side is used. 
#' @param data Optional data.frame.
#' @param beta_prior A vector of length two which defines the 
#' prior mean and standard deviation of the Gaussian prior for coefficients at time 1.
#' @param sigma_prior A vector of length two, defining the truncated Gaussian prior for 
#' the slope level standard deviation. 
#' @param slope_prior A vector of length two which defines the 
#' prior mean and standard deviation of the Gaussian prior for the slopes at time 1.
#' @export
rw2 <- function(formula, data, beta_prior, sigma_prior, slope_prior) {
  
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- as.name("na.pass")
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  xreg <- model.matrix(attr(mf, "terms"), mf)
  
  if(length(beta_prior) != 2) {
    stop("beta_prior should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of initial coefficients. ")
  }
  if(length(sigma_prior) != 2) {
    stop("sigma_prior should be should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of standard deviations. ")
  }
  if(length(slope_prior) != 2) {
    stop("slope_prior should be should be a vector of length two, defining the mean and standard deviation for the Gaussian prior of initial slope coeffients. ")
  }
  list(xreg = xreg, beta_prior = beta_prior, 
    sigma_prior = sigma_prior, slope_prior = slope_prior)
  
}