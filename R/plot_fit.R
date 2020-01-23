#' Plot the fitted values and sample quantiles for a walker object
#' 
#' @param object An output from \code{\link{walker}} or \code{\link{walker_glm}}.
#' @param level Level for intervals. Default is 0.05, leading to 90\% intervals.
#' @param alpha Transparency level for \code{geom_ribbon}.
#' @param ... Further arguments to \code{\link{ppc_ribbon}}.
#' @export
plot_fit <- function(object, level = 0.05, alpha = 0.33, ...){
  
  y_fit <- extract(object$stanfit, pars = "y_fit", permuted = TRUE)$y_fit
  if (object$distribution != "gaussian") {
    y_fit <- y_fit[sample(1:nrow(y_fit), size = nrow(y_fit), replace = TRUE, 
      prob = extract(object$stanfit, pars = "weights", permuted = TRUE)$weights),  , drop = FALSE]
  }
  noNA <- which(!is.na(object$y))
  ppc_ribbon(y = as.numeric(object$y[noNA]), 
    yrep = y_fit[,noNA],
    x = as.numeric(time(object$y))[noNA], ...) + 
    theme(legend.position = "none") + 
    scale_x_continuous(name = "time")
}