# walker: Baysian dynamic linear regression models with Stan/R

Walker provides a method for fully Bayesian linear regression where the 
regression coeffients are allowed to vary over "time", either as independent random walks. 

All computations are done using Hamiltonian Monte Carlo provided by Stan, 
using a state space representation of the model in order to marginalise over the coefficients for efficient sampling.
