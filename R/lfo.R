#' Leave-Future-Out Cross-Validation
#'
#' Estimates the leave-future-out (LFO) information criterion for \code{walker} and \code{walker_glm} models.
#' 
#' The LFO for non-Gaussian models is (currently) based on the corresponding Gaussian approximation and 
#' not the importance sampling corrected true posterior. 
#' 
#' @export
#' @importFrom loo psis pareto_k_values weights.importance_sampling
#' @param object Output of \code{walker} or \code{walker_glm}.
#' @param L Positive integer defining how many observations should be used for the initial fit.
#' @param exact If \code{TRUE}, computes exact 1-step predictions by re-estimating the model repeatedly. 
#' If \code{FALSE} (default), uses approximate method based on Bürkner, Gabry and Vehtari (2020).
#' @param verbose If \code{TRUE} (default), print the progress of the LFO computations to the console.
#' @param k_thres Threshold for the pareto k estimate triggering refit. Default is 0.7. 
#' @references Paul-Christian Bürkner, Jonah Gabry & Aki Vehtari (2020). 
#' Approximate leave-future-out cross-validation for Bayesian time series models, 
#' Journal of Statistical Computation and Simulation, 90:14, 2499-2523, DOI: 10.1080/00949655.2020.1783262.
#' @return List with components \code{ELPD} (Expected log predictive density), \code{ELPDs} (observation-specific ELPDs),
#' \code{ks} (Pareto k values in case of approximation was used), and \code{refits} (time points where model was re-estimated)
#' @examples 
#' \dontrun{
#' fit <- walker(Nile ~ -1 + 
#'   rw1(~ 1, 
#'     beta = c(1000, 100), 
#'     sigma = c(2, 0.001)), 
#'   sigma_y_prior = c(2, 0.005), 
#'   iter = 2000, chains = 1)
#'  
#' fit_lfo <- lfo(fit, L = 20, exact = FALSE)
#' fit_lfo$ELPD
#' }
lfo <- function(object, L, exact = FALSE, verbose = TRUE, k_thres = 0.7) {
  
  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
  }
  sum_log_ratios <- function(loglik, ids = NULL) {
    if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
    rowSums(loglik)
  }
  
  if (is.null(object$data)) {
    stop("Data used to fit the model is missing. Please rerun the model with argument `return_data = TRUE`.")
  }
  d <- object$data
  
  n_samples <- sum(object$stanfit@sim$n_save - object$stanfit@sim$warmup2)
  
  # use posterior means as initial values to make refitting more robust
  samples <- extract(object$stanfit)
  
  if (object$distribution == "gaussian") {
    # initial values
    init <- replicate(object$stanfit@sim$chains, list(
      beta_fixed = if(!is.null(samples$beta_fixed)) 
        array(colMeans(samples$beta_fixed), dim = ncol(samples$beta_fixed)),
      sigma_rw1 = if(!is.null(samples$sigma_rw1)) 
        array(colMeans(samples$sigma_rw1), dim = ncol(samples$sigma_rw1)),
      sigma_rw2 = if(!is.null(samples$sigma_rw2)) 
        array(colMeans(samples$sigma_rw2), dim = ncol(samples$sigma_rw2)),
      sigma_y = mean(samples$sigma_y)), simplify = FALSE)
    
    if (exact) {
      ks <- NULL
      refits <-  L:(d$n - 1)
      ll <- matrix(NA, n_samples, d$n - L)
      for(i in L:(d$n - 1)) {
        if (verbose) print(paste0("Estimating model with ", i, " observations.")) 
        # increase number of observations used for estimating the parameters
        d$n_lfo <- i
        f <- stan(fit = object$stanfit, data = d, 
          init = init, 
          chains = object$stanfit@sim$chains,
          iter = object$stanfit@sim$iter,
          warmup = object$stanfit@sim$warmup,
          thin = object$stanfit@sim$thin,
          pars = "logLik",
          refresh = 0
        )
        ll[, i - L + 1] <- extract(f, "logLik")$logLik[, i + 1]
      }
      elpds <- apply(ll, 2, log_mean_exp)
      elpd <- sum(elpds)
    } else {
      # Based on the Bürkner et al.: https://mc-stan.org/loo/articles/loo2-lfo.html
      elpds <- numeric(d$n - L)
      d$n_lfo <- L
      if (verbose) print(paste0("Estimating model with ", L, " observations.")) 
      f <- stan(fit = object$stanfit, data = d, 
        init = init, 
        chains = object$stanfit@sim$chains,
        iter = object$stanfit@sim$iter,
        warmup = object$stanfit@sim$warmup,
        thin = object$stanfit@sim$thin,
        pars = "logLik",
        refresh = 0
      )
      elpds[1] <- log_mean_exp(extract(f, "logLik")$logLik[, L + 1])
      i_refit <- L
      refits <- L
      ks <- numeric(d$n - L - 1)
      
      for (i in (L + 1):(d$n - 1)) {
        loglik <- extract(f, "logLik")$logLik
        logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
        psis_obj <- suppressWarnings(loo::psis(logratio))
        k <- loo::pareto_k_values(psis_obj)
        ks[i] <- k
        if (k > k_thres) {
          if (verbose) print(paste0("Estimating model with ", i, " observations.")) 
          # refit the model based on the first i observations
          i_refit <- i
          refits <- c(refits, i)
          d$n_lfo <- i
          f <- stan(fit = object$stanfit, data = d, 
            init = init, 
            chains = object$stanfit@sim$chains,
            iter = object$stanfit@sim$iter,
            warmup = object$stanfit@sim$warmup,
            thin = object$stanfit@sim$thin,
            pars = "logLik",
            refresh = 0
          )
          loglik <- extract(f, "logLik")$logLik
          elpds[i - L + 1] <- log_mean_exp(loglik[, i + 1])
        } else {
          lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)[, 1]
          elpds[i - L + 1] <- log_sum_exp(lw + loglik[, i + 1])
        }
      }
      elpd <- sum(elpds)
    }
  } else {
    warning("LFO for non-Gaussian models is based on the approximating Gaussian model.")
    
    init <- replicate(object$stanfit@sim$chains, list(
      beta_fixed = if(!is.null(samples$beta_fixed)) 
        array(colMeans(samples$beta_fixed), dim = ncol(samples$beta_fixed)),
      sigma_rw1 = if(!is.null(samples$sigma_rw1)) 
        array(colMeans(samples$sigma_rw1), dim = ncol(samples$sigma_rw1)),
      sigma_rw2 = if(!is.null(samples$sigma_rw2)) 
        array(colMeans(samples$sigma_rw2), dim = ncol(samples$sigma_rw2))), 
      simplify = FALSE)
    
    if (exact) {
      ks <- NULL
      refits <-  L:(d$n - 1)
      ll <- matrix(NA, n_samples, d$n - L)
      for(i in L:(d$n - 1)) {
        if (verbose) print(paste0("Estimating model with ", i, " observations.")) 
        d$n_lfo <- i
        f <- stan(fit = object$stanfit, data = d, 
          init = init, 
          chains = object$stanfit@sim$chains,
          iter = object$stanfit@sim$iter,
          warmup = object$stanfit@sim$warmup,
          thin = object$stanfit@sim$thin,
          pars = "logLik",
          refresh = 0
        )
        ll[, i - L + 1] <- extract(f, "logLik")$logLik[, i + 1]
      }
      elpds <- apply(ll, 2, log_mean_exp)
      elpd <- sum(elpds)
      
    } else {
      # Based on the Bürkner et al.: https://mc-stan.org/loo/articles/loo2-lfo.html
      elpds <- numeric(d$n - L)
      d$n_lfo <- L
      f <- stan(fit = object$stanfit, data = d, 
        init = init, 
        chains = object$stanfit@sim$chains,
        iter = object$stanfit@sim$iter,
        warmup = object$stanfit@sim$warmup,
        thin = object$stanfit@sim$thin,
        pars = "logLik",
        refresh = 0
      )
      elpds[1] <- log_mean_exp(extract(f, "logLik")$logLik[, L + 1])
      i_refit <- L
      refits <- L
      ks <- numeric(d$n - L - 1)
      
      for (i in (L + 1):(d$n - 1)) {
        loglik <- extract(f, "logLik")$logLik
        logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
        psis_obj <- suppressWarnings(loo::psis(logratio))
        k <- loo::pareto_k_values(psis_obj)
        ks[i] <- k
        if (k > k_thres) {
          if (verbose) print(paste0("Estimating model with ", i, " observations.")) 
          # refit the model based on the first i observations
          i_refit <- i
          refits <- c(refits, i)
          d$n_lfo <- i
          f <- stan(fit = object$stanfit, data = d, 
            init = init, 
            chains = object$stanfit@sim$chains,
            iter = object$stanfit@sim$iter,
            warmup = object$stanfit@sim$warmup,
            thin = object$stanfit@sim$thin,
            pars = "logLik",
            refresh = 0
          )
          loglik <- extract(f, "logLik")$logLik
          elpds[i - L + 1] <- log_mean_exp(loglik[, i + 1])
        } else {
          lw <- loo::weights.importance_sampling(psis_obj, normalize = TRUE)[, 1]
          elpds[i - L + 1] <- log_sum_exp(lw + loglik[, i + 1])
        }
      }
      elpd <- sum(elpds)
    }
  }
  list(ELPD = elpd, ELPDs = elpds, ks = ks, refits = refits)
}
