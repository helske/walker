#' Predictions for walker object
#' 
#' Given the new covariate data and output from \code{walker}, 
#' obtain samples from posterior predictive distribution for counterfactual case,
#' i.e. for past time points with different covariate values.
#' 
#' @importFrom stats deltat tsp rpois plogis rbinom
#' @param object An output from \code{\link{walker}} or \code{\link{walker_glm}}.
#' @param newdata A \code{data.frame} containing covariates used for prediction. 
#' Should have equal number of rows as the original data
#' @param u For Poisson model, a vector of exposures i.e. E(y) = u*exp(x*beta). 
#' For binomial, a vector containing the number of trials. Defaults 1.
#' @param summary If \code{TRUE} (default), return summary statistics. Otherwise returns samples.
#' @param type If \code{"response"} (default for Gaussian model), predictions are on the response level 
#' (e.g., number of successes for Binomial case, and for Gaussian case the observational 
#' level noise is added to the mean predictions).
#' If \code{"mean"} (default for non-Gaussian case), predict means (e.g., success probabilities in Binomial case).
#' If \code{"link"}, predictions for non-Gaussian models are returned before applying the inverse of the link-function.
#' @return If \code{summary=TRUE}, time series containing summary statistics of predicted values. 
#' Otherwise a matrix of samples from predictive distribution.
#' @export
#' @examples 
#' \dontrun{
#' set.seed(1)
#' n <- 50
#' x1 <- rnorm(n, 0, 1)
#' x2 <- rnorm(n, 1, 0.5)
#' x3 <- rnorm(n)
#' beta1 <- cumsum(c(1, rnorm(n - 1, sd = 0.1)))
#' beta2 <- cumsum(c(0, rnorm(n - 1, sd = 0.1)))
#' beta3 <- -1
#' u <- sample(1:10, size = n, replace = TRUE)
#' y <- rbinom(n, u, plogis(beta3 * x3 + beta1 * x1 + beta2 * x2))
#' 
#' d <- data.frame(y, x1, x2, x3)
#' out <- walker_glm(y ~ x3 + rw1(~ -1 + x1 + x2, beta = c(0, 2), 
#'   sigma = c(2, 10)), distribution = "binomial", beta = c(0, 2), 
#'   u = u, data = d,
#'   iter = 2000, chains = 1, refresh = 0)
#' 
#' # what if our covariates were constant?
#' newdata <- data.frame(x1 = rep(0.4, n), x2 = 1, x3 = -0.1)
#' 
#' fitted <- fitted(out)
#' pred <- predict_counterfactual(out, newdata, type = "mean")
#' 
#' ts.plot(cbind(fitted[, c(1, 3, 5)], pred[, c(1, 3, 5)]), 
#'   col = rep(1:2, each = 3), lty = c(1, 2, 2))
#'}
predict_counterfactual <- function(object, newdata, u, summary = TRUE,
  type = ifelse(object$distribution == "gaussian", "response", "mean")){
  
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
  if (n != nrow(newdata)) stop("Number of rows in 'newdata' should match with the original data. ")
  
  beta_fixed <- extract(object$stanfit, pars = "beta_fixed")$beta_fixed
  beta_rw <- 
    extract(object$stanfit, pars = "beta_rw")$beta_rw
  n_iter <- nrow(beta_rw)
  
  y_new <- matrix(0, n, n_iter)
  if (!is.null(beta_fixed)) {
    beta_fixed <- t(beta_fixed)
    xregs$xreg_fixed <- t(xregs$xreg_fixed)
    for (i in 1:n_iter) {
      for(t in 1:n) {
        y_new[t, i] <- beta_fixed[, i] %*% xregs$xreg_fixed[, t]
      }
    }
  }
  
  for (i in 1:n_iter) {
    for (t in 1:n) {
      y_new[t, i] <- y_new[t, i] + beta_rw[i, , t] %*% xregs$xreg_rw[t, ]
    }
  }
  
  if (object$distribution != "gaussian") {
    if (missing(u)) {
      u <- rep(1, nrow(newdata))
    } else {
      if (length(u) != nrow(newdata)) {
        if(length(u) > 1) stop("Length of 'u' should be 1 or equal to the number of predicted time points. ")
        u <- rep(u, nrow(newdata))
      }
    }
    if (type != "link") {
      if (object$distribution == "poisson") {
        for (i in 1:n_iter) {
          y_new[, i] <- u * exp(y_new[, i])
        }
      } else {
        y_new <- plogis(y_new)
      }
      if (type == "response") {
        if (object$distribution == "poisson") {
          for (i in 1:n_iter) {
            y_new[, i] <- rpois(n, y_new[,i])
          }
        } else {
          for (i in 1:n_iter) {
            y_new[, i] <- rbinom(n, u, y_new[, i])
          }
        }
      }
    }
  }
  if (summary) {
    y_new <- t(apply(y_new, 1, function(x) {
      q <- quantile(x, c(0.025, 0.5, 0.975))
      c(mean = mean(x), sd = sd(x), q)
    }))
    rownames(y_new) <- time(object$y)
  }
  y_new
}
