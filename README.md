[![Build Status](https://travis-ci.org/helske/walker.png?branch=master)](https://travis-ci.org/helske/walker)

# walker: Efficient Baysian dynamic linear regression models with Stan/R

Walker provides a method for fully Bayesian linear regression where the 
regression coeffients are allowed to vary over "time", either as independent random walks. 

All computations are done using Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for accurate and efficient sampling.

See the package [vignette](http://htmlpreview.github.io/?https://github.com/helske/walker/blob/master/walker_html/walker.html) for details and an example.
