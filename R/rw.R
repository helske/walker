#' @export
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