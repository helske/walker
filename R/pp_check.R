#' Posterior predictive check for walker object
#' 
#' Plots sample quantiles from posterior predictive sample. 
#' See \code{\link{ppc_ribbon}} for details.
#' @importFrom bayesplot pp_check
#' @method pp_check walker_fit
#' @param object An output from \code{\link{walker}}.
#' @param ... Further parameters to \code{\link{ppc_ribbon}}.
#' @export
pp_check.walker_fit <- function(object, ...){
  
  y_rep <- extract(object$stanfit, pars = "y_rep", permuted = TRUE)$y_rep
  if (object$distribution != "gaussian") {
    y_rep <- y_rep[sample(1:nrow(y_rep), size = nrow(y_rep), replace = TRUE, 
      prob = extract(object$stanfit, pars = "weights", permuted = TRUE)$weights), , drop = FALSE]
  }
  
  ppc_ribbon(y = as.numeric(object$y), 
    yrep = y_rep,
    x = as.numeric(time(object$y)),
    ...)
}