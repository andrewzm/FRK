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



#' Computation and concatenation of percentiles to a dataframe.
#'
#' Computes the percentiles or HPD interval bounds at each prediction location and appends 
#' the result to \code{data}. Note that we use percentiles rather than quantiles
#' because we including a "dot" (corresponding to the decimal place) in the 
#' dataframe column name may cause issues. 
#'
#' @inheritParams .prediction_interval
#' @param data The dataframe we will append percentiles to; the number of rows (if \code{X} is a matrix) or the length (if \code{X} is a list) of \code{X} must equal the number of rows in \code{data}
#' @param name The name of the quantity of interest. The names of the percentile columns are "name_percentile". In \code{FRK} it will be Y, mu, prob, or Z
#' @return The dataframe \code{data} with appended percentiles
.concat_percentiles_to_df <- function (data, X, name, percentiles, cred_mass) {
  if (is.null(percentiles) & is.null(cred_mass)) 
    return(data)
  
  
  if (!is.null(percentiles)) {
    Q           <- .prediction_interval(X, interval_type = "central", percentiles = percentiles, cred_mass = cred_mass)
    colnames(Q) <- paste(name, "percentile", as.character(percentiles), sep = "_")
    data        <- cbind(data, Q)
  }
  
  if(!is.null(cred_mass)) {
    Q           <- .prediction_interval(X, interval_type = "HPD", percentiles = percentiles, cred_mass = cred_mass)
    colnames(Q) <- paste(colnames(Q), "HPD_bound", name, sep = "_")
    data        <- cbind(data, Q)
  }
  
  return(data)
}


#' Prediction intervals. 
#'
#' Given a matrix of Monte Carlo samples \code{X} (wherein rows correspond to 
#' locations, columns correspond to samples), this function computes symmetric 
#' or HPD prediction intervals, and returns the prediction interval lower and 
#' upper bound
#' 
#' @param X a matrix (wherein rows correspond to prediction locations and columns are samples) of Monte Carlo samples
#' @param interval_type string indicating whether a \code{"central"} or \code{"HPD"} interval is desired
#' @param percentiles a vector of scalars in [0, 100] specifying the desired percentiles; if \code{percentiles = NULL}, no percentiles are computed 
#' @param cred_mass a scalar [0, 1] specifying the mass within the credible interval; if \code{cred_mass = NULL}, no HPD credible interval is computed
#' @return The prediction interval at each location (or width at each location if \code{width = TRUE})
.prediction_interval <- function(X, interval_type, percentiles, cred_mass) {

  X <- as.matrix(X)
  if (interval_type == "central") Q <- t(apply(X, 1, quantile, percentiles / 100))
  if (interval_type == "HPD")     Q <- t(HDInterval::hdi(t(X), cred_mass = cred_mass)) # Note that hdi() expects locations in columns, samples in rows
  
  return(Q)
}


## Since we will use ggplot2 we will first convert our objects to data frames.
## The fortify command in ggplot2 was found to not work well with the sp polys.
## The following function converts a SpatialPolygonsDataFrame to a data.frame.
## sp_polys: Object of class SpatialPolygonsDataFrame
## vars: vector of characters identifying the field names to retain
.SpatialPolygonsDataFrame_to_df <- function (sp_polys, vars = names(sp_polys)) {
  ## Extract the names of the polygons
  polynames <- as.character(row.names(sp_polys))
  
  ## For each polygon do the following
  list_polys <- lapply(1:length(sp_polys), function(i) {
    ## Extract coordinates
    coords <- sp_polys@polygons[[i]]@Polygons[[1]]@coords
    row.names(coords) <- NULL
    coords <- data.frame(coords)
    
    ## Create data frame with coordinates and with "id" as the polygon name
    poldf <- cbind(coords, id = polynames[i], stringsAsFactors = FALSE)
    rownames(poldf) <- NULL
    poldf
  })
  
  ## Bind all the returned data frames together
  df_polys <- bind_rows(list_polys)
  
  ## Make sure the data frame and SpatialPolygonsDataFrame have the same "key" id var
  df_polys$id <- as.character(df_polys$id)
  sp_polys$id <- row.names(sp_polys)
  
  cnames <- coordnames(sp_polys)                     # coordinate names
  vars_no_coords <- vars[which(!vars %in% cnames)]   # all selected var names
  
  ## Now bind the variables to the data frame
  if (length(vars_no_coords) > 0)
    df_polys <- left_join(df_polys, sp_polys@data[c("id",
                                                    vars_no_coords)], by = "id")
  df_polys
}
