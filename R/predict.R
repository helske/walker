#' Predictions for walker object
#' 
#' Given the new covariate data and output from `walker`, 
#' obtain samples from posterior predictive distribution for future time points.
#' 
#' @importFrom stats deltat tsp
#' @param object An output from [walker()] or [walker_glm()].
#' @param newdata A `data.frame` containing covariates used for prediction.
#' @param u For Poisson model, a vector of future exposures i.e. E(y) = u*exp(x*beta). 
#' For binomial, a vector containing the number of trials for future time points. Defaults 1.
#' @param type If `"response"` (default for Gaussian model), predictions are on the response level 
#' (e.g., number of successes for Binomial case, and for Gaussian case the observational 
#' level noise is added to the mean predictions).
#' If `"mean"` (default for non-Gaussian case), predict means (e.g., success probabilities in Binomial case).
#' If `"link"`, predictions for non-Gaussian models are returned before applying the inverse of the link-function.
#' @param ... Ignored.
#' @return A list containing samples from posterior predictive distribution.
#' @method predict walker_fit
#' @seealso [plot_predict()] for example.
#' @export
predict.walker_fit <- function(object, newdata, u, 
  type = ifelse(object$distribution == "gaussian", "response", "mean"), ...){
  
  type <- match.arg(type, c("response", "mean", "link"))
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
  
  nu <- 
    extract(object$stanfit, pars = "nu")$nu[, , n, drop = FALSE]
  if (is.null(nu)) {
    nu <- matrix(0, n_iter, 0)
  } else  dim(nu) <- dim(nu)[1:2]
  
  beta_fixed <- extract(object$stanfit, pars = "beta_fixed")$beta_fixed
  if (is.null(beta_fixed)) beta_fixed <- matrix(0, n_iter, 0)
  
  sigma_rw1 <- extract(object$stanfit, pars = "sigma_rw1")$sigma_rw1
  if (is.null(sigma_rw1)) sigma_rw1 <- matrix(0, n_iter, 0)
  
  sigma_rw2 <- extract(object$stanfit, pars = "sigma_rw2")$sigma_rw2
  if (is.null(sigma_rw2)) sigma_rw2 <- matrix(0, n_iter, 0)
  
  if (object$distribution != "gaussian") {
    if (missing(u)) {
      u <- rep(1, nrow(newdata))
    } else {
      if (length(u) != nrow(newdata)) {
        if(length(u) > 1) stop("Length of 'u' should be 1 or equal to the number of predicted time points. ")
        u <- rep(u, nrow(newdata))
      }
    }
    type_int <- pmatch(type, c("link", "response", "mean")) - 1L
    pred <- predict_walker_glm(t(sigma_rw1), t(sigma_rw2),
      t(beta_fixed), t(beta_rw), t(nu), xregs$xreg_fixed, 
      t(xregs$xreg_rw), u, 
      pmatch(object$distribution, c("poisson", "binomial")), 
      extract(object$stanfit, pars = "weights")$weights, 
      nrow(newdata), ncol(beta_fixed), ncol(sigma_rw1), ncol(sigma_rw2), 
      type_int)
    
    
    pred$mean <- colMeans(fitted(object, summary = FALSE))
    pred$u <- object$u
    
  } else {
    if (type == "link") type <- "mean"
    
    pred <- predict_walker(t(sigma_rw1), t(sigma_rw2),
      extract(object$stanfit, pars = "sigma_y")$sigma_y,
      t(beta_fixed), t(beta_rw), t(nu), xregs$xreg_fixed, 
      t(xregs$xreg_rw), nrow(newdata), ncol(beta_fixed), ncol(sigma_rw1), 
      ncol(sigma_rw2), type == "response")
    
    pred$mean <- colMeans(fitted(object, summary = FALSE))
  }

  
  st <-  tsp(object$y)[2L]
  if (is.null(st)) st <- length(object$y)
  s1 <- tsp(object$y)[1L]
  if (is.null(s1)) s1 <- 1
  d <- deltat(object$y)
  pred$y <- object$y
  pred$mean <- ts(pred$mean, start = s1, end = st, deltat = d)
  
  attr(pred, "type") <- type
  dimnames(pred$y_new) <- 
    list(time = seq(st + deltat(object$y), by = deltat(object$y), 
      length = nrow(pred$y_new)), iter = 1:ncol(pred$y_new))
  
  pred
}
