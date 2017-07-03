#' Plot the fitted values and sample quantiles for a walker object
#' 
#' @param object An output from \code{\link{predict.walker_fit}}.
#' @param level Level for intervals. Default is 0.05, leading to 90\% intervals.
#' @param alpha Transparency level for \code{geom_ribbon}.
#' @param ... Further arguments to \code{\link[bayesplot]{ppc_ribbon}}.
#' @export
plot_fit <- function(object, level = 0.05, alpha = 0.33, ...){
  
  ppc_ribbon(y = as.numeric(object$y), 
    yrep = extract(object$stanfit, pars = "y_fit", permuted = TRUE)$y_fit,
    x = as.numeric(time(object$y)), ...) + theme(legend.position = "none") + 
    scale_x_continuous(name = "time")
}