% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.walker_fit}
\alias{predict.walker_fit}
\title{Predictions for walker object}
\usage{
\method{predict}{walker_fit}(
  object,
  newdata,
  u,
  type = ifelse(object$distribution == "gaussian", "response", "mean"),
  ...
)
}
\arguments{
\item{object}{An output from \code{\link[=walker]{walker()}} or \code{\link[=walker_glm]{walker_glm()}}.}

\item{newdata}{A \code{data.frame} containing covariates used for prediction.}

\item{u}{For Poisson model, a vector of future exposures i.e. E(y) = u\emph{exp(x}beta).
For binomial, a vector containing the number of trials for future time points. Defaults 1.}

\item{type}{If \code{"response"} (default for Gaussian model), predictions are on the response level
(e.g., number of successes for Binomial case, and for Gaussian case the observational
level noise is added to the mean predictions).
If \code{"mean"} (default for non-Gaussian case), predict means (e.g., success probabilities in Binomial case).
If \code{"link"}, predictions for non-Gaussian models are returned before applying the inverse of the link-function.}

\item{...}{Ignored.}
}
\value{
A list containing samples from posterior predictive distribution.
}
\description{
Given the new covariate data and output from \code{walker},
obtain samples from posterior predictive distribution for future time points.
}
\seealso{
\code{\link[=plot_predict]{plot_predict()}} for example.
}
