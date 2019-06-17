#' Extract Fitted Values of Walker Fit
#'
#' Returns fitted values (posterior means) from output of \code{walker} or \code{walker_glm}.
#' 
#' @export
#' @importFrom stats fitted sd
#' @name fitted.walker_fit
#' @param object Output of \code{walker} or \code{walker_glm}.
#' @param summary If \code{TRUE} (default), return summary statistics. Otherwise returns samples.
#' @param ... Ignored.
#' @return Time series containing fitted values.
fitted.walker_fit <- function(object, summary = TRUE, ...) {
  
  y_fit <- extract(object$stanfit, pars = "y_fit", permuted = TRUE)$y_fit
  if (object$distribution != "gaussian") {
    y_fit <- y_fit[sample(1:nrow(y_fit), size = nrow(y_fit), replace = TRUE, 
      prob = extract(object$stanfit, pars = "weights", permuted = TRUE)$weights),  , drop = FALSE]
  }
  if (summary) {
    y_fit <- t(apply(y_fit, 2, function(x) {
      q <- quantile(x, c(0.025, 0.5, 0.975))
      c(mean = mean(x), sd = sd(x), q)
    }))
    rownames(y_fit) <- time(object$y)
    y_fit
  } else {
    y_fit
  }
}

#' Extract Coeffients of Walker Fit
#'
#' Returns the regression coeffients from output of \code{walker} or \code{walker_glm}.
#' 
#' @export
#' @importFrom stats coef
#' @name coef.walker_fit
#' @param object Output of \code{walker} or \code{walker_glm}.
#' @param summary If \code{TRUE} (default), return summary statistics. Otherwise returns samples.
#' @param transform Optional vectorized function for transforming the coefficients (for example exp).
#' @param ... Ignored.
#' @return Time series containing coeffients values.
coef.walker_fit <- function(object, summary = TRUE, transform = identity,  ...) {
  # N x k x n array
  coef_data <- transform(extract(object$stanfit, pars = "beta_rw", permuted = TRUE)$beta)
  
  if (object$distribution != "gaussian") {
    coef_data <- coef_data[sample(1:nrow(coef_data), size = nrow(coef_data), replace = TRUE, 
      prob = extract(object$stanfit, pars = "weights", permuted = TRUE)$weights), , , drop = FALSE]
  }
  dimnames(coef_data) <- 
    list(iter = 1:nrow(coef_data), 
      beta = colnames(object$xreg_rw), 
      time = as.numeric(time(object$y)))
  
  if (summary) {
    coef_data <- as.data.frame(as.table(coef_data))  
    names(coef_data)[4] <- "value"
    coef_data$time <- as.numeric(levels(coef_data$time))[coef_data$time]
    grouped <- group_by(coef_data, time, beta)
    summarise_(grouped, 
      .dots = list(
        mean = ~mean(value),
        sd = ~sd(value),
        "2.5%" = ~quantile(value, prob = 0.025), 
        "50%" = ~quantile(value, prob = 0.5),
        "97.5%" = ~quantile(value, prob = 0.975)))
  } else {
    coef_data  
  }
}