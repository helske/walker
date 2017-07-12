[![Build Status](https://travis-ci.org/helske/walker.png?branch=master)](https://travis-ci.org/helske/walker)
[![cran version](http://www.r-pkg.org/badges/version/walker)](http://cran.r-project.org/package=walker)
[![downloads](http://cranlogs.r-pkg.org/badges/walker)](http://cranlogs.r-pkg.org/badges/walker)

# walker: Efficient Baysian dynamic linear regression models with Stan/R

Walker provides a method for fully Bayesian generalized linear regression where the 
regression coefficients are allowed to vary over "time" as a first or second order integrated random walk. 

All computations are done using Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.

See the package [vignette](http://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) for details and an examples.

