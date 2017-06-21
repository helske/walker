[![cran version](http://www.r-pkg.org/badges/version/walker)](http://cran.r-project.org/package=walker)

# walker: Efficient Baysian dynamic linear regression models with Stan/R

Walker provides a method for fully Bayesian linear regression where the 
regression coefficients are allowed to vary over "time", either as independent random walks. 

All computations are done using Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.

See the package [vignette](http://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) for details and an example.

## It is possible to extend walker to time-varying GLMs as well. Unfortunately I do not currently have time for implement it, but if someone wants to do it please contact me.
