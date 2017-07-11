#' Predictions for walker object
#' 
#' Given the new covariate data and output from \code{walker}, 
#' obtain samples from posterior predictive distribution.
#' 
#' @importFrom stats deltat tsp
#' @param object An output from \code{\link{walker}} or \code{\link{walker_glm}}.
#' @param newdata A data.frame containing covariates used for prediction.
#' @param ... Ignored.
#' @return A list containing samples from posterior predictive distribution.
#' @method predict walker_fit
#' @export
predict.walker_fit <- function(object, newdata, ...){
  
  y_name <- as.character(object$call$formula[[2]])
  
  if (!(y_name%in% names(newdata))) {
    newdata[[y_name]] <- rep(NA, nrow(newdata))
  }
  object$call$data <- newdata
  object$call$return_x_reg <- TRUE
  xregs <- eval(object$call)
  
  if (any(is.na(xregs$xreg_fixed)) || any(is.na(xregs$xreg_rw))) {
    stop("Missing values in covariates are not allowed.")
  }
  
  n <- length(object$y)
  
  beta_rw <- 
    extract(object$stanfit, pars = "beta_rw")$beta_rw[, , n, drop = FALSE]
  dim(beta_rw) <- dim(beta_rw)[1:2]
  n_iter <- nrow(beta_rw)
  
  slope <- 
    extract(object$stanfit, pars = "slope")$slope[, , n, drop = FALSE]
  if (is.null(slope)) {
    slope <- matrix(0, n_iter, 0)
  } else  dim(slope) <- dim(slope)[1:2]
  
  beta_fixed <- extract(object$stanfit, pars = "beta_fixed")$beta_fixed
  if (is.null(beta_fixed)) beta_fixed <- matrix(0, n_iter, 0)
  
  sigma_rw1 <- extract(object$stanfit, pars = "sigma_rw1")$sigma_rw1
  if (is.null(sigma_rw1)) sigma_rw1 <- matrix(0, n_iter, 0)
  
  sigma_rw2 <- extract(object$stanfit, pars = "sigma_rw2")$sigma_rw2
  if (is.null(sigma_rw2)) sigma_rw2 <- matrix(0, n_iter, 0)
  
  if (object$distribution != "gaussian") {
  pred <- predict_walker_glm(t(sigma_rw1), t(sigma_rw2),
    t(beta_fixed), t(beta_rw), t(slope), t(xregs$xreg_fixed), 
    t(xregs$xreg_rw), object$u, 
    pmatch(object$distribution, c("poisson", "binomial")), 
    extract(object$stanfit, pars = "weights")$weights)
  } else {
    pred <- predict_walker(t(sigma_rw1), t(sigma_rw2),
      extract(object$stanfit, pars = "sigma_y")$sigma_y,
      t(beta_fixed), t(beta_rw), t(slope), t(xregs$xreg_fixed), 
      t(xregs$xreg_rw))
  }
  pred$y <- object$y
 
  st <-  tsp(object$y)[2L]
  if (is.null(st)) st <- length(object$y)
  dimnames(pred$y_new) <- 
    list(time = seq(st + deltat(object$y), by = deltat(object$y), 
      length = nrow(pred$y_new)), iter = 1:ncol(pred$y_new))
  pred
}