% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGeneric.R
\name{coef_uncertainty}
\alias{coef_uncertainty}
\title{Uncertainty quantification of the fixed effects}
\usage{
coef_uncertainty(
  object,
  percentiles = c(5, 95),
  nsim = 400,
  random_effects = FALSE
)
}
\arguments{
\item{object}{object of class \code{SRE} returned from the constructor \code{SRE()} containing all the parameters and information on the SRE model}

\item{percentiles}{(applicable only if \code{method} = "TMB") a vector of scalars in (0, 100) specifying the desired percentiles of the posterior predictive distribution; if \code{NULL}, no percentiles are computed}

\item{nsim}{number of Monte Carlo samples used to compute the confidence intervals}

\item{random_effects}{logical; if set to true, confidence intervals will also be provided for the random effects random effects \ifelse{html}{\out{<i> <b>&gamma;</b></i>}}{\eqn{\boldsymbol{\gamma}}{\gamma}} (see `?SRE` for details on these random effects)}
}
\description{
Compute confidence intervals for the fixed effects (upper  
and lower bound specifed by percentiles; default 90\% confidence central interval)
}
