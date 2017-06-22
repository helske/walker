[![cran version](http://www.r-pkg.org/badges/version/walker)](http://cran.r-project.org/package=walker)

# walker: Efficient Baysian dynamic linear regression models with Stan/R

Walker provides a method for fully Bayesian linear regression where the 
regression coefficients are allowed to vary over "time", either as independent random walks. 

*Update: walker now supports also Poisson regression with time-varying coefficients!*

All computations are done using Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.

See the package [vignette](http://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) for details and an example.

## The package is fully functional as is, but more work is needed for better automatic handling of the output, visualization tools, testing, etc (see issues). Pull requests are very welcome.
