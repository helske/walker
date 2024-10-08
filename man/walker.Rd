% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/walker.R
\name{walker}
\alias{walker}
\title{Bayesian regression with random walk coefficients}
\usage{
walker(
  formula,
  data,
  sigma_y_prior = c(2, 0.01),
  beta,
  init,
  chains,
  return_x_reg = FALSE,
  gamma_y = NULL,
  return_data = TRUE,
  ...
)
}
\arguments{
\item{formula}{An object of class \code{{formula}} with additional terms
\code{rw1} and/or \code{rw2} e.g. \code{y ~ x1 + rw1(~ -1 + x2)}. See details.}

\item{data}{An optional data.frame or object coercible to such, as in \code{{lm}}.}

\item{sigma_y_prior}{A vector of length two, defining the a Gamma prior for
the observation level standard deviation with first element corresponding to the shape parameter and
second to rate parameter. Default is Gamma(2, 0.0001). Not used in \code{walker_glm}.}

\item{beta}{A length vector of length two which defines the
prior mean and standard deviation of the Gaussian prior for time-invariant coefficients}

\item{init}{Initial value specification, see \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.
Note that compared to default in \code{rstan}, here the default is a to sample from the priors.}

\item{chains}{Number of Markov chains. Default is 4.}

\item{return_x_reg}{If \code{TRUE}, does not perform sampling, but instead returns the matrix of
predictors after processing the \code{formula}.}

\item{gamma_y}{An optional vector defining known non-negative weights for the standard
deviation of the observational level noise at each time point.
More specifically, the observational level standard deviation sigma_t is
defined as \eqn{\sigma_t = gamma_t * \sigma_y} (in default case
\eqn{\sigma_t = sigma_y})}

\item{return_data}{if \code{TRUE}, returns data input to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.
This is needed for \code{lfo}.}

\item{...}{Further arguments to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.}
}
\value{
A list containing the \code{stanfit} object, observations \code{y},
and covariates \code{xreg} and \code{xreg_new}.
}
\description{
Function \code{walker} performs Bayesian inference of a linear
regression model with time-varying, random walk regression coefficients,
i.e. ordinary regression model where instead of constant coefficients the
coefficients follow first or second order random walks.
All Markov chain Monte Carlo computations are done using Hamiltonian
Monte Carlo provided by Stan, using a state space representation of the model
in order to marginalise over the coefficients for efficient sampling.
}
\details{
The \code{rw1} and \code{rw2} functions used in the formula define new formulas
for the first and second order random walks. In addition, these functions
need to be supplied with priors for initial coefficients and the
standard deviations. For second order random walk model, these sigma priors
correspond to the standard deviation of slope disturbances. For \code{rw2},
also a prior for the initial slope nu needs to be defined. See examples.
}
\note{
Beware of overfitting and identifiability issues. In particular,
be careful in not defining multiple intercept terms
(only one should be present).
By default \code{rw1} and \code{rw2} calls add their own time-varying
intercepts, so you should use \code{0} or \code{-1} to remove some of them
(or the time-invariant intercept in the fixed-part of the formula).
}
\examples{

\dontrun{
set.seed(1)
x <- rnorm(10)
y <- x + rnorm(10)

# different intercept definitions:

# both fixed intercept and time-varying level,
# can be unidentifiable without strong priors:
fit1 <- walker(y ~ rw1(~ x, beta = c(0, 1)), 
  beta = c(0, 1), chains = 1, iter = 1000, init = 0) 

# only time-varying level, using 0 or -1 removes intercept:
fit2 <- walker(y ~ 0 + rw1(~ x, beta = c(0, 1)), chains = 1, iter = 1000,
  init = 0)

# time-varying level, no covariates:
fit3 <- walker(y ~ 0 + rw1(~ 1, beta = c(0, 1)), chains = 1, iter = 1000)

# fixed intercept no time-varying level:
fit4 <- walker(y ~ rw1(~ 0 + x, beta = c(0, 1)), 
  beta = c(0, 1), chains = 1, iter = 1000) 

# only time-varying effect of x:
fit5 <- walker(y ~ 0 + rw1(~ 0 + x, beta = c(0, 1)), chains = 1, iter = 1000) 
}

\dontrun{

rw1_fit <- walker(Nile ~ -1 + 
  rw1(~ 1, 
    beta = c(1000, 100), 
    sigma = c(2, 0.001)), 
  sigma_y_prior = c(2, 0.005), 
  iter = 2000, chains = 1)
  
rw2_fit <- walker(Nile ~ -1 + 
  rw2(~ 1,
    beta = c(1000, 100), 
    sigma = c(2, 0.001), 
    nu = c(0, 100)), 
  sigma_y_prior = c(2, 0.005), 
  iter = 2000, chains = 1)
  
g_y <- geom_point(data = data.frame(y = Nile, x = time(Nile)), 
  aes(x, y, alpha = 0.5), inherit.aes = FALSE) 
g_rw1 <- plot_coefs(rw1_fit) + g_y
g_rw2 <- plot_coefs(rw2_fit) + g_y
if(require("gridExtra")) {
  grid.arrange(g_rw1, g_rw2, ncol=2, top = "RW1 (left) versus RW2 (right)")
} else {
  g_rw1
  g_rw2
}

y <- window(log10(UKgas), end = time(UKgas)[100])
n <- 100
cos_t <- cos(2 * pi * 1:n / 4)
sin_t <- sin(2 * pi * 1:n / 4)
dat <- data.frame(y, cos_t, sin_t)
fit <- walker(y ~ -1 + 
  rw1(~ cos_t + sin_t, beta = c(0, 10), sigma = c(2, 1)), 
  sigma_y_prior = c(2, 10), data = dat, chains = 1, iter = 2000)
print(fit$stanfit, pars = c("sigma_y", "sigma_rw1"))

plot_coefs(fit)
# posterior predictive check:
pp_check(fit)

newdata <- data.frame(
  cos_t = cos(2 * pi * 101:108 / 4), 
  sin_t = sin(2 * pi * 101:108 / 4))
pred <- predict(fit, newdata)
plot_predict(pred)

# example on scalability
set.seed(1)
n <- 2^12
beta1 <- cumsum(c(0.5, rnorm(n - 1, 0, sd = 0.05)))
beta2 <- cumsum(c(-1, rnorm(n - 1, 0, sd = 0.15)))
x1 <- rnorm(n, mean = 2)
x2 <- cos(1:n)
rw <- cumsum(rnorm(n, 0, 0.5))
signal <- rw + beta1 * x1 + beta2 * x2
y <- rnorm(n, signal, 0.5)

d <- data.frame(y, x1, x2)

n <- 2^(6:12)
times <- numeric(length(n))
for(i in seq_along(n)) {
  times[i] <- sum(get_elapsed_time(
    walker(y ~ 0 + rw1(~ x1 + x2, 
      beta = c(0, 10)), 
      data = d[1:n[i],],
      chains = 1, seed = 1, refresh = 0)$stanfit))
}
plot(log2(n), log2(times))
}
}
\seealso{
\code{\link[=walker_glm]{walker_glm()}} for non-Gaussian models.
}
