#' Posterior predictive check for walker object
#' 
#' Plots sample quantiles from posterior predictive sample. 
#' 
#' @details 
#' For other types of posterior predictive checks for example with `bayesplot`, 
#' you can extract the variable `yrep` from the output, see examples.
#' 
#' @importFrom bayesplot pp_check ppc_ribbon
#' @param object An output from [walker()].
#' @param ... Further parameters to [bayesplot::ppc_ribbon()].
#' @export
#' @examples 
#' \dontrun{
#' # Extracting the yrep variable for general use:
#' # extract yrep
#' y_rep <- extract(object$stanfit, pars = "y_rep", permuted = TRUE)$y_rep
#' 
#' # For non-gaussian model:
#' weights <- extract(object$stanfit, 
#' pars = "weights", permuted = TRUE)$weights
#' y_rep <- y_rep[sample(1:nrow(y_rep), 
#'   size = nrow(y_rep), replace = TRUE, prob = weights), , drop = FALSE]
#'}
#' 
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
