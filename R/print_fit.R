#' Print Summary of walker_fit Object
#' 
#' Prints the summary information of time-invariant model parameters. In case of non-Gaussian models, 
#' results based on approximate model are returned with a warning.
#' 
#' @param x An output from [walker()] or [walker_glm()].
#' @param ... Additional arguments to [rstan::print.stanfit()].
#' @method print walker_fit
#' @export
print.walker_fit <- function(x, ...) {
  pars <- setdiff(x$stanfit@sim$pars_oi, c("beta_rw", "nu", "y_fit", "y_rep", "log_lik"))
  if(x$distribution != "gaussian") warning("Results are based on approximate model, use summary method for exact results.")
  print(x$stanfit, pars = pars, ...)
}

#' Coerce Posterior Samples of walker Fit to a Data Frame
#' 
#' Creates a data.frame object from the output of walker fit.
#' 
#' @param x An output from [walker()] or [walker_glm()].
#' @param row.names `NULL` (default) or a character vector giving the row names 
#' for the data frame.
#' @param optional Ignored (part of generic `as.data.frame` signature).
#' @param type Either `tiv` (time-invariant parameters) or `tv` (time-varying coefficients).
#' @param ... Ignored.
#' @method as.data.frame walker_fit
#' @export
#' @examples 
#' \dontrun{
#'  as.data.frame(fit, "tiv") %>% 
#'  group_by(variable) %>%
#'  summarise(mean = mean(value),
#'            lwr = quantile(value, 0.05),
#'            upr = quantile(value, 0.95))
#' }
#'
as.data.frame.walker_fit <- function(x, row.names = NULL, optional = FALSE,  type, ...) {
  
  type <- match.arg(type, c("tiv", "tv"))
  
  if (type == "tiv") {
    pars <- setdiff(x$stanfit@sim$pars_oi, c("beta_rw", "nu", "y_fit", "y_rep", "lp__", "weights", "log_lik"))
    samples <- extract(x$stanfit, pars = pars, permuted = FALSE)
    n <- nrow(samples)
    k <- ncol(samples)
    d <- data.frame(iter = 1:n,
                    chain = rep(1:k, each = n),
                    value = c(samples), 
                    variable = rep(dimnames(samples)[[3]], each = n * k),
                    row.names = row.names)
    if (x$distribution != "gaussian") {
      d$weight <- c(extract(x$stanfit, pars = "weights", permuted = FALSE))
    }
  } else {
    pars <- intersect(x$stanfit@sim$pars_oi, c("beta_rw", "nu"))
    samples <- extract(x$stanfit, pars = pars, permuted = FALSE)
    n <- nrow(samples)
    k <- ncol(samples)
    d <- data.frame(iter = 1:n,
                    chain = rep(1:k, each = n),
                    time = rep(as.numeric(time(x$y)), each = n * k * ncol(x$xreg_rw)),
                    value = c(samples), 
                    variable = rep(paste0("beta_", colnames(x$xreg_rw)), 
                                   each = n * k),
                    row.names = row.names)
    if (x$distribution != "gaussian") {
      d$weight <- c(extract(x$stanfit, pars = "weights", permuted = FALSE))
    }
  }
  d
}


#' Summary of walker_fit Object
#' 
#' Return summary information of time-invariant model parameters.
#' 
#' @param object An output from [walker()] or [walker_glm()].
#' @param type Either `tiv` (time-invariant parameters, the default) or `tv` (time-varying coefficients).
#' @param ... Ignored.
#' @importFrom Hmisc wtd.mean wtd.var wtd.quantile
#' @importFrom coda spectrum0.ar
#' @method summary walker_fit
#' @export
summary.walker_fit <- function(object, type = "tiv", ...) {
  
  type <- match.arg(type, c("tiv", "tv"))
  
  if (type == "tiv") {
    pars <- setdiff(object$stanfit@sim$pars_oi, c("beta_rw", "nu", "y_fit", "y_rep", "lp__", "weights"))
  } else {
    pars <- intersect(object$stanfit@sim$pars_oi, c("beta_rw", "nu"))
  }
  
  if (object$distribution == "gaussian") {
    
    d <- as.data.frame(summary(object$stanfit, pars = pars)$summary[, c("mean", "se_mean", "sd", "2.5%", "97.5%", "n_eff")])
    d$n_eff <- round(d$n_eff)
    
  } else {
    
    samples <- extract(object$stanfit, pars = pars, permuted = FALSE)
    w <- extract(object$stanfit, pars = "weights", permuted = FALSE)
 
    means <- apply(samples, 3, function(x) wtd.mean(x, c(w), normwt = TRUE))
    sds <- sqrt(apply(samples, 3, function(x) wtd.var(x, c(w), normwt = TRUE)))
    lwrs <- apply(samples, 3, function(x) wtd.quantile(x, c(w), 0.025, normwt = TRUE))
    uprs <- apply(samples, 3, function(x) wtd.quantile(x, c(w), 0.975, normwt = TRUE))
    
    ess <- numeric(dim(samples)[3])
    for (i in 1:seq_along(ncol(samples))) {
      c2_est <- mean(w[, i, 1])^2 * length(w[, i, 1])
      for (j in seq_along(ess)) {
        v <- spectrum0.ar(samples[,i,j] * w[, i, 1])$spec / c2_est
        ess[j] <- ess[j] + sds[j]^2 / v
      }
    }
    d <- data.frame(mean = means, se_mean = sds / sqrt(ess),
                    sd = sds, "2.5%" = lwrs, 
                    "97.5%" = uprs, n_eff = round(ess),
                    row.names = dimnames(samples)[[3]],
                    check.names = FALSE)
    
  }
  d
}
