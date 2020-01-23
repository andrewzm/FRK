#' Link and inverse-link function generators.
#'
#' Create the link function \eqn{g}  and inverse-link function \eqn{\psi},
#' which link the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y} such that \eqn{g(\mu) = Y} and
#' \eqn{\mu = \psi(Y)}, respectively.
#'
#' For families lacking a "known constant" parameter,
#' \code{.link_fn} and \code{.inv_link_fn} generate the functions linking
#' the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y}. However, for families with a "known constant"
#' parameter (such as the "binomial" and "negative-binomial"), \code{.link_fn} and
#' \code{.inv_link_fn} generate the functions linking
#' the probability parameter \eqn{p} to the latent
#' geostatistical process \eqn{Y}; then, the conditional mean \eqn{mu} is constructed from \eqn{p}
#' using the known form of the mean for that distribution.
#'
#' @describeIn LINK_FNS create the link function \eqn{g} which
#' links the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y}, such that \eqn{g(\mu) = Y}.
#'
#' @param link A string indicating the assumed link function.
#' @return The link function or inverse-link function.
#' @examples
#' g <- .link_fn("log")
#' g(20)
.link_fn <- function(link = c("log", "identity", "logit", "probit", "cloglog", "inverse", "inverse-squared", "square-root")) {
  
  ## Argument check
  link <- tolower(link)   # make input lower case
  link <- match.arg(link) # check that specified link is valid
  
  if (link == "log")             g <- function(mu) log(mu)
  else if (link == "identity")   g <- function(mu) mu
  else if (link == "logit")      g <- function(mu) log(mu/(1-mu))
  else if (link == "probit")     g <- function(mu) qnorm(mu)
  else if (link == "cloglog")    g <- function(mu) log(-log(1 - mu))
  else if (link == "inverse") g <- function(mu) 1/mu
  else if (link == "inverse-squared") g <- function(mu) 1/(mu^2)
  else if (link == "square-root")g <- function (mu) sqrt(mu)
  
  return(g)
}



#' @describeIn LINK_FNS create the inverse link function \eqn{\psi}
#' which links the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y}, such that \eqn{\mu = \psi(Y)}.
#'
#' @examples
#' psi <- .inv_link_fn("log")
#' psi(3)
.inv_link_fn <- function(link = c("log", "identity", "logit", "probit", "cloglog", "inverse", "inverse-squared", "square-root")) {
  
  ## Argument check
  link <- tolower(link)
  link <- match.arg(link)
  
  
  if (link == "log")             psi <- function(Y) exp(Y)
  else if (link == "identity")   psi <- function(Y) Y
  else if (link == "logit")      psi <- function(Y) 1/(1 + exp(-Y))
  else if (link == "probit")     psi <- function(Y) pnorm(Y)
  else if (link == "cloglog")    psi <- function(Y) 1 - exp(-exp(Y))
  else if (link == "inverse") psi <- function(Y) 1/Y
  else if (link == "inverse-squared") psi <- function(Y) 1/(sqrt(Y))
  else if (link == "square-root")psi <- function(Y) Y^2
  
  return(psi)
}


#' Variance by rows.
#'
#' Computes the row-wise variance of a matrix.
#'
#' @param X An array like object (with a dim-attribute).
#' @return A vector containing the row-wise variances of the matrix \code{X}.
.rowVars <- function(X, ...) {
  rowSums((X - rowMeans(X, ...))^2, ...)/(dim(X)[2] - 1)
}

## FIX: Add documentation to this function
.concat_percentiles_to_df <- function (samples, df, name) {
  temp <- t(apply(samples, 1, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))
  colnames(temp) <- paste(name, sep = "_", c("percentile_05", "percentile_25", "percentile_50", "percentile_75", "percentile_95"))
  df <- cbind(df, temp)
  return(df)
}
