#' Posterior predictive check for walker object
#' 
#' Plots sample quantiles from posterior predictive sample. 
#' See \code{\link[bayesplot]{ppc_ribbon}} for details.
#' @importFrom bayesplot pp_check
#' @method pp_check walker_fit
#' @param Object An output from \code{\link{walker}}.
#' @param ... Further parameters to \code{\link[bayesplot]{ppc_ribbon}}.
#' @export
pp_check.walker_fit <- function(object, ...){
  
  ppc_ribbon(y = as.numeric(object$y), 
    yrep = extract(object$stanfit, pars = "y_rep", permuted = TRUE)$y_rep,
    x = as.numeric(time(object$y)),
    ...)
}