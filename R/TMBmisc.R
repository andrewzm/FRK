#' Argument checker for non-Gaussian FRK.
#'
#' Checks that the arguments relevant to FRKTMB functions are valid.
#' This function does not make any alterations to its arguments,
#' it simply calls \code{stop()} if an invalid argument is detected.
#'
#' @inheritParams .TMB_prep
#' @inheritParams .Y_pred
#' @inheritParams .MC_sampler
.check_arg_FRKTMB <- function(response = c("gaussian", "poisson", "bernoulli", "gamma",
                                           "inverse-gaussian", "negative-binomial", "binomial"),
                              link = c("identity", "log", "square-root", "logit", "probit", "cloglog", "inverse", "inverse-squared"),
                              K_type = c("precision", "block-exponential"),
                              taper,
                              k_Z,
                              k_BAU) {
  
  ## ??? problem with using match.arg() in these functions
  ## NEED TO MAKE THIS CHECK FOR EXACT MATCHES, BECAUSE
  ## match.arg(K_type)
  ## WOULD ALLOW "prec" as a valid input for K_type.
  ## This is fine if we can change the K_type of the SRE object,
  ## But the way this function is written this is not possible.
  ## The alternative is to accept M as the argument,
  ## then have M <- .check_arg_FRKTMB(M).
  
  ## Check K_type
  if (!missing(K_type)) {
    match.arg(K_type)
  }
  
  
  
  ## Check that a valid response-link combination has been chosen
  if (!missing(response) & !missing(link)) {
    
    match.arg(response)
    match.arg(link)
    
    if (response == "gaussian" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
        response == "poisson" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
        response == "gamma" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
        response == "inverse-gaussian" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
        response == "negative-binomial" & !(link %in% c("log", "square-root", "logit", "probit", "cloglog")) ||
        response == "binomial" & !(link %in% c("logit", "probit", "cloglog")) ||
        response == "bernoulli" & !(link %in% c("logit", "probit", "cloglog"))) {
      stop("Invalid response-link combination selected")
    }
    
  }
  
  
  ## Check taper
  if(!missing(taper)) {
    if (!is.numeric(taper) | taper <= 0){
      stop("taper must be positive.")
    }
  }
  
  
  ## Check k_Z
  if (!missing(k_Z)) {
    if (response %in% c("binomial", "negative-binomial")) {
      
      if (is.null(k_Z)) {
        stop("For binomial or negative-binomial data, the known constant parameter k must be provided for each observation.")
      } else if (!is.numeric(k_Z) | any(k_Z <= 0) | any(k_Z != round (k_Z))) {
        stop("The known constant parameter k must contain only positive integers.")
      }
      
    }
  }
  
  ## Check k_BAU
  if (!missing(k_BAU)) {
    if (response %in% c("binomial", "negative-binomial")) {
      
      if (is.null(k_BAU)){
        stop("For binomial or negative-binomial response, the known constant parameter k must be provided for each BAU.")
      } else if (!is.numeric(k_BAU) | any(k_BAU <= 0) | any(k_BAU != round (k_BAU))) {
        stop("The known constant parameter k must contain only positive integers.")
      } else if (length(k_BAU) != nrow(M@S0)) {
        stop("length(k) not equal to the number of BAUs." )
      }
      
    }
  }
  
}




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