[![R-CMD-check](https://github.com/helske/walker/workflows/R-CMD-check/badge.svg)](https://github.com/helske/walker/actions)
[![cran version](https://www.r-pkg.org/badges/version/walker)](https://cran.r-project.org/package=walker)
[![downloads](https://cranlogs.r-pkg.org/badges/walker)](https://cranlogs.r-pkg.org/badges/walker)

walker: Bayesian Generalized Linear Models with Time-Varying Coefficients
==========================================================================

The R package walker provides a method for fully Bayesian generalized linear regression where the 
regression coefficients are allowed to vary over time as a first or second order integrated random walk. 

The Markov chain Monte Carlo (MCMC) algorithm uses Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.
For non-Gaussian models the MCMC targets approximate marginal posterior based on Gaussian approximation, which is then corrected using importance sampling as in [Vihola, Helske, Franks (2020)](https://onlinelibrary.wiley.com/doi/10.1111/sjos.12492).

See the corresponding paper in [softwareX](https://www.sciencedirect.com/science/article/pii/S235271102200022X) for short introduction, and the package [vignette](https://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) and [documentation manual](https://cran.r-project.org/package=walker/walker.pdf) for details and further examples.

You can download the development version of `walker` from Github using the [`devtools`](https://cran.r-project.org/package=devtools) package:

```R
devtools::install_github("helske/walker")
```

NEWS
---------------------------------------------

### 3.3.2022, version 1.0.4

* Added an example of counterfactual predictions to the vignette.
* Added citation info for the softwareX paper.

### 24.9.2021, version 1.0.3-1

* Added a flag for stanc3 compatibility, thanks to Andrew Johnson. Also added an 
  import for RcppParallel.

### 16.8.2021

* Internal changes to make `walker` compatible with upcoming `StanHeaders`.

### 6.4.2021

* Changed the name of the `logLik` variable to `log_lik` so it is compatible with `loo`.

### 27.1.2021

* Fixed some issues in the vignette which resulted CRAN warnings.

### 25.1.2021

* For linear-Gaussian models the stanfit object now returns partial log-likelihood terms
  p(y_t | y_1,...,y_t-1,theta) which can be used for leave-future-out cross-validation (see function `lfo`). 
* New function `lfo` for estimating the leave-future-out information criterion.
* Priors for the standard deviation parameters are now Gamma instead of truncated normal, which helps to avoid (rare) problems where sampler wonders close to degenerate case of having all variances near zero. There are also default prior Gamma(2, 0.0001) for these parameters now.
* Fixed some issues in the vignette added a reference to the walker paper.
  
### 3.11.2020

* stanfit object of walker output now contains also variable `logLik`.
  For non-Gaussian models this is the approximate log-likelihood, the
  unbiased estimate is then `logLik + mean(w)`, where `w` are the returned weights.

### 19.10.2020

* Predict method now allow predictions on link scale.
* Added argument for plot_predict for controlling the drawing of past observations.
* Fix out-of-sample predictions for non-Gaussian models.
* New function: `predict_counterfactual` which can be used to predict the past assuming new 
  values for the covariates.

### 13.8.2020

* Proper export of `pp_check` for `bayesplot`, fixed some minor technical issues.

### 19.5.2020

* Added default values for `row.names` and `optional` for `as.data.frame` function.

### 12.5.2020

* Added as.data.frame function for `walker` and `walker_glm` output.
* Added a `summary` method.
* The print method now correctly warns about approximate results in case of non-Gaussian model.
* Changed arguments `*_prior` to more concise versions (e.g. `sigma_prior` is now just `sigma`). 
* Changed the name of the slope terms to `nu` as in vignette formulas. 
* Updated to rstantools 2.0.0 package structure and removed dependency on soft-depracated functions of `dplyr`.

### 23.1.2020

* Removed check for missing values in function `walker` which threw an error even though missing values in responses have been in principle supported since 2018...

### 20.9.2019

* Switched from GPL2+ to GPL3 in order to be compatible with future Stan versions.

### 04.03.2019

* Added methods fitted and coef for extracting the posterior means and and regression coefficents from the 
  walker_fit object.
* Fixed issue with Makevars and clang4 per request by CRAN.
* Added option to predict on mean-scale, e.g, probabilities instead of 0/1 in Bernoulli case.
* Fixed a bug in the Gaussian predictions, last time point was missing the observational level noise.

### 25.02.2019

* Issue with upcoming staged installation in CRAN fixed by Tomas Kalibera.

### 14.02.2019

* Dimension bug in GLM case fixed.

### 8.11.2018

* Fixed StanHeaders search in Makevars. 

### 22.10.2018

* Pull request by Ben Goodrich for fixing the issue with clang4. New version on it's way to CRAN.

### 15.10.2018
* Missing values in response variable are now supported.
* Added gamma variables to models which can be used to damp the variance of the random walks. 
* Tidied some Stan codes in order to reduce deep copying.
* Moved stan codes under `src`.
* Increased the iteration counts in examples in order to pass CRAN tests.
<
