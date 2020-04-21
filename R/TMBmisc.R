#' Link and inverse-link function generator.
#'
#' Create the link function \eqn{g}  and inverse-link function \eqn{\psi},
#' which link the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y} such that \eqn{g(\mu) = Y} and
#' \eqn{\mu = \psi(Y)}, respectively.
#'
#' For families lacking a "known constant" parameter,
#' \code{.link_fn} generates the functions linking
#' the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y}. However, for families with a "known constant"
#' parameter (such as the "binomial" and "negative-binomial"), \code{.link_fn} generate the functions linking
#' the probability parameter \eqn{p} to the latent
#' geostatistical process \eqn{Y}, and then the parameter \eqn{p} to the conditional mean \eqn{\mu}.
#'
#' @param kind A character string indicating which kind of link function is desired. 
#' Valid values are 
#' \describe{
#'   \item{\code{"Y_to_mu"}}{Provides the function \eqn{\psi} such that \eqn{\mu = \psi(Y)}.}
#'   \item{\code{"mu_to_Y"}}{Provides the function \eqn{g} such that \eqn{g(\mu) = Y}.}
#'   \item{\code{"Y_to_prob"}}{Provides the function \eqn{\zeta} such that \eqn{p = \zeta(Y)}.}
#'   \item{\code{"prob_to_Y"}}{Provides the function \eqn{h} such that \eqn{h(p) = Y}.}
#'   \item{\code{"prob_to_mu"}}{Provides the function \eqn{\chi} such that \eqn{\mu = \chi(p)}.}
#'   \item{\code{"mu_to_prob"}}{Provides the function \eqn{f} such that \eqn{f(\mu) = p}.}
#' }
#' Note that the latter four values are relevant only to the binomial and negative-binomial distributions with logit, probit, or cloglog link functions. 
#' @param link A character string indicating the assumed link function. \emph{Not} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}. 
#' @param response A character string indicating the assumed response distribution. \emph{Only} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}.
#' @return A function.
#' @examples
#' g <- .link_fn(kind = "mu_to_Y", link = "log")
#' g(20)
#' psi <- .link_fn(kind = "Y_to_mu", link = "log")
#' psi(3)
#' chi <- .link_fn(kind = "prob_to_mu, response = "negative-binomial")
.link_fn <- function (kind, link, response) {
  
  if (kind == "Y_to_mu") {
    if (link == "log")             psi <- function(Y) exp(Y)
    if (link == "identity")        psi <- function(Y) Y
    if (link == "logit")           psi <- function(Y) 1/(1 + exp(-Y))
    if (link == "probit")          psi <- function(Y) pnorm(Y)
    if (link == "cloglog")         psi <- function(Y) 1 - exp(-exp(Y))
    if (link == "inverse")         psi <- function(Y) 1/Y
    if (link == "inverse-squared") psi <- function(Y) 1/(sqrt(Y))
    if (link == "square-root")     psi <- function(Y) Y^2
    return(psi)
    
  } else if (kind == "mu_to_Y") {
    if (link == "log")             g <- function(mu) log(mu)
    if (link == "identity")        g <- function(mu) mu
    if (link == "logit")           g <- function(mu) log( mu /(1 - mu))
    if (link == "probit")          g <- function(mu) qnorm(mu)
    if (link == "cloglog")         g <- function(mu) log(-log(1 - mu))
    if (link == "inverse")         g <- function(mu) 1/mu
    if (link == "inverse-squared") g <- function(mu) 1/(mu^2)
    if (link == "square-root")     g <- function (mu) sqrt(mu)
    return(g)
    
  } else if (kind == "Y_to_prob") {
    if (link == "logit")      zeta <- function(Y) 1/(1 + exp(-Y))
    if (link == "probit")     zeta <- function(Y) pnorm(Y)
    if (link == "cloglog")    zeta <- function(Y) 1 - exp(-exp(Y))
    return(zeta)
    
  } else if (kind == "prob_to_Y") {
    if (link == "logit")      h <- function(p) log(p /(1 - p))
    if (link == "probit")     h <- function(p) qnorm(p)
    if (link == "cloglog")    h <- function(p) log(-log(1 - p))
    return(h)
    
  } else if (kind == "prob_to_mu") {
    if (response == "binomial")          chi <- function(p, k) k * p
    if (response == "negative-binomial") chi <- function(p, k) k * (1 / p - 1)
    return(chi)
    
  } else if (kind == "mu_to_prob") {
    if (response == "binomial")          f <- function(mu, k) mu / k
    if (response == "negative-binomial") f <- function(mu, k) k / (k + mu)
    return(f)
  }
  
  stop("Invalid arguments.")
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


#' Computation and concatenation of percentiles to a dataframe.
#'
#' Given a matrix of MC samples \code{X} (where rows correspond to locations and 
#' columns correspond to samples), this function computes the percentiles at 
#' each location and appends the result to \code{df} (which must have the same 
#' number of rows as \code{X}). Note that we use percentiles rather than quantiles
#' because we including a "dot" in the dataframe column name may cause issues. 
#'
#' @param X A matrix of Monte Carlo samples (rows are observation, columns are 
#' samples).
#' @param df The prediction dataframe, with the same number of rows as \code{X}.
#' @param name The name of the quantity of interest. The names of the percentile
#' columns are "name_percentile_percents". In \code{FRK} it will be Y, mu, or Z. 
#' @param percents A vector containing the desired percentiles which will be 
#' included in the prediction dataframe. If \code{NULL}, 
#' no percentiles are computed.
#' @return The dataframe \code{df} with appended percentiles.
.concat_percentiles_to_df <- function (X, df, name, 
                                       percents = c(5, 25, 50, 75, 95)) {
  if (is.null(percents)) return(df)
  temp           <- t(apply(X, 1, quantile, percents/100))
  colnames(temp) <- paste(name, "percentile", as.character(percents), sep = "_")
  df             <- cbind(df, temp)
  return(df)
}


#' Prediction interval width. 
#'
#' This function is used to compute prediction interval width
#' given \code{X}, a matrix of Monte Carlo samples (wherein rows correpsond to locations, columns samples)
#' 
#' @param X A matrix of Monte Carlo samples (rows are observation, columns are 
#' samples).
#' @param l The lower probability.
#' @param u The upper probability.
#' @return The prediction interval width at each location.
.intervalWidth <- function(X, l = 0.025, u = 0.975) {
  Q <- apply(X, 1, function(x) quantile(x, c(l, u)))
  return(Q[2, ] - Q[1, ])
}