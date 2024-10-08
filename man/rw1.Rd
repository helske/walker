% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rw.R
\name{rw1}
\alias{rw1}
\title{Construct a first-order random walk component}
\usage{
rw1(formula, data, beta, sigma = c(2, 1e-04), gamma = NULL)
}
\arguments{
\item{formula}{Formula for RW1 part of the model. Only right-hand-side is used.}

\item{data}{Optional data.frame.}

\item{beta}{A length vector of length two which defines the
prior mean and standard deviation of the Gaussian prior for coefficients at time 1.}

\item{sigma}{A vector of length two, defining the Gamma prior for
the coefficient level standard deviation. First element corresponds to the shape parameter and
second to the rate parameter. Default is Gamma(2, 0.0001).}

\item{gamma}{An optional k times n matrix defining a known non-negative weights of the
random walk noises, where k is the number of coefficients and n is the
number of time points. Then, the standard deviation of the random walk noise
for each coefficient is of form gamma_t * sigma (instead of just sigma).}
}
\description{
Auxiliary function used inside of the formula of \code{walker}.
}
