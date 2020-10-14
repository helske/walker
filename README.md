[![cran version](https://www.r-pkg.org/badges/version/walker)](https://cran.r-project.org/package=walker)
[![downloads](https://cranlogs.r-pkg.org/badges/walker)](https://cranlogs.r-pkg.org/badges/walker)

walker: Bayesian Generalized Linear Models with Time-Varying Coefficients
==========================================================================

The R package walker provides a method for fully Bayesian generalized linear regression where the 
regression coefficients are allowed to vary over time as a first or second order integrated random walk. 

The Markov chain Monte Carlo (MCMC) algorithm uses Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.
For non-Gaussian models the MCMC targets approximate marginal posterior based on Gaussian approximation, which is then corrected using importance sampling as in [Vihola, Helske, Franks (2020)](https://arxiv.org/abs/1609.02541v6).

See the package [vignette](https://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) and [documentation manual](https://cran.r-project.org/package=walker/walker.pdf) for details and examples.

You can download the development version of `walker` from Github using the [`devtools`](https://cran.r-project.org/package=devtools) package:

```R
devtools::install_github("helske/walker")
```

NEWS
---------------------------------------------
### 14.10.2020

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
