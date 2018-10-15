#' Print Summary of walker_fit Object
#' 
#' Prints the summary information of time-invariant model parameters.
#' 
#' @param x An output from \code{\link{walker}} or \code{\link{walker_glm}}.
#' @param ... Additional arguments to \code{\link{print.stanfit}}.
#' @method print walker_fit
#' @export
print.walker_fit <- function(x, ...) {
  pars <- setdiff(x$stanfit@sim$pars_oi, c("beta_rw", "slope", "y_fit", "y_rep"))
  print(x$stanfit, pars = pars, ...)
}