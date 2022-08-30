.response_name <- function(object) all.vars(object@f)[1]

#' @rdname SRE
#' @export
simulate <- function(object, newdata = NULL, nsim = 400, conditional_fs = FALSE, ...) {
  
  if (is.null(newdata)) {
    if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the fitted values will only be computed for the first data set.")
    newdata <- object@data[[1]]
  }
  
  if (object@method == "EM") {
    pred <- predict(object, newdata = newdata, ...) # TODO need conditional_fs for EM side. Once implemented, add conditional_fs = conditional_fs to this call to predict()
    Zhat       <- pred$mu
    sigma2e    <- object@Ve[1, 1]
    sdZ        <- sqrt(pred$sd^2 + sigma2e)
    MC_samples <- t(sapply(seq_along(Zhat), function(i) rnorm(nsim, Zhat[i], sdZ[i]))) 
  } else {
    pred <- predict(object, newdata = newdata, type = "response", percentiles = NULL, 
                    conditional_fs = conditional_fs, nsim = nsim, ...)
    MC_samples <- pred$MC$Z_samples
  }
  
  return(MC_samples)
}


# TODO marginal simulation (i.e., unconditionally on the random effects) would 
# be straightforward too, and we can use predict() to avoid code repetition. 
# Perhaps the best way to implement it would be by adding an argument, 
# conditional_basis, which controls if we are to simulate conditionally on the 
# fitted basis-function coefficients. If both conditional_fs = FALSE and 
# conditional_basis = FALSE, then we have marginal simulation.
#   S <- object@S 
#   Q <- object@Khat_inv 
#   L <- t(chol(Q))
#   r <- nrow(L)
#   eta <- matrix(rnorm(r * nsim), nrow = r, ncol = nsim)


#' @rdname SRE
#' @export
setMethod("fitted", signature="SRE", function(object, ...) {
  
  if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the fitted values will only be computed for the first data set.")
  newdata <- object@data[[1]]
  
  if (object@method == "EM") {
    Zhat <- predict(object, newdata = newdata, ...)$mu
  } else {
    Zhat <- predict(object, newdata = newdata, type = "mean", percentiles = NULL, ...)$newdata$p_mu
  }
  
  return(Zhat)
})


#' @rdname SRE
#' @export
setMethod("residuals", signature="SRE", function(object, type = "pearson") { 

  # Extract the data
  if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the residuals will only be computed for the first data set.")
  newdata <- object@data[[1]]
  Z       <- newdata@data[, .response_name(object)]
  
  # Fitted values. Use simulate() so that we can estimate the variance when 
  # computing pearson residuals. 
  Zsim <- simulate(object, newdata, conditional_fs = TRUE)
  Zhat <- rowMeans(Zsim)
  
  # Compute the "response" residuals
  res <- Zhat - Z
  
  if (type == "pearson") {
    Vhat <- apply(Zsim, 1, var)
    res <- res/sqrt(Vhat)
  } 
  
  return(res)
})

#' @rdname SRE
#' @export
setMethod("AIC", signature="SRE", function(object, k = 2) {
  
  ## Initialise the number of covariance parameters to zero
  p <- 0
  
  ## Fine-scale variance parameters
  if (object@fs_by_spatial_BAU) {
    ## If we are estimating a unique fine-scale variance parameter at each 
    ## spatial BAU, we will have as many fine-scale variance parameters as 
    ## spatial BAUs
    p <- p + dim(object@BAUs)[1]
  } else {
    ## Otherwise, we just have a single fine-scale variance parameter
    p <- p + 1
  }
  
  ## Dispersion parameter
  ## The following distributions do not have a dispersion parameter
  if (!(object@response %in% c("poisson", "binomial", "negative-binomial"))) {
    p <- p + 1
  }
  
  ## Basis-function covariance parameters
  if (is(object@basis,"TensorP_Basis")) {
    # We model the temporal basis-function coefficients as an AR1 process, which 
    # has a marginal variance parameter and a correlation parameter.
    temporal_cov_parameters <- 2 
    spatial_basis <- object@basis@Basis1
  } else {
    temporal_cov_parameters <- 0
    spatial_basis <- object@basis
  }

  R <- nres(spatial_basis)
  if (object@K_type == "block-exponential") {
    spatial_cov_parameters <- 2 * R   
  } else if (object@K_type == "precision") {
    
    if (object@basis@regular) {
      spatial_cov_parameters <- 2 * R 
    } else {
      spatial_cov_parameters <- 3 * R 
    }
    
  } else if (object@K_type == "unstructured") {
    r <- table(spatial_basis@df$res) # number of spatial basis functions at each resolution
    spatial_cov_parameters <- sum(r*(r+1)/2)
  }
  
  p <- p + spatial_cov_parameters + temporal_cov_parameters 
  
  ## Number of fixed effects
  q <- length(coef(object))
  
  ## The criterion
  k * (p + q) - 2 * loglik(object)
})



#' @rdname SRE
#' @export
setMethod("BIC", signature="SRE", function(object) AIC(object, k = log(nobs(object))))






# #' @title Compute model diagnostics
# #' @description Takes an object of class \code{SRE} and computes several diagnostics useful for model assessment. The default diagnostics are the root-mean-squared error (RMSE), the mean-absolute error (MAE), the continuous ranked probability score (CRPS), and the empirical coverage and interval scores resulting from a prediction interval with a nominal coverage of (1-\code{alpha}).
# #' @param object object of class \code{SRE}.
# #' @param newdata an object of class \code{Spatial*DataFrame} or \code{STFDF} which contains the response variable and any covariates in the model. If \code{NULL} (default), the diagnostics are computed with respect to the data used for model fitting.
# #' @param diagnostic_fns a named list of user-specified diagnostic functions (i.e., default diagnostics will also be computed). Each function should take as input a list and return a scalar value. The list that is passed into each function has elements:
# #'  \itemize{
# #'  \item{\code{Z}}{The data,}
# #'  \item{\code{Zhat}}{The predictions,}
# #'  \item{\code{Zhat_lower}}{The lower bound of the prediction interval,}
# #'  \item{\code{Zhat_upper}}{The upper bound of the prediction interval,}
# #'  \item{\code{MC_samples}}{A matrix of Monte Carlo samples, where rows correspond to observations in \code{Z}.}
# #'  }
# #' Note that a given function does not need to make use of all of the above elements (e.g., the RMSPE would simply use \code{Z} and \code{Zhat}).
# #' @param nominal_coverage the nominal coverage of the prediction intervals, which should be strictly between 0 (0% coverage) and 1 (100% coverage).
# #' @export
# #' @examples
# #' # See example in the help file for FRK
# setGeneric("diagnostics", function(object, newdata = NULL, diagnostic_fns = NULL, nominal_coverage = 0.9)
#   standardGeneric("diagnostics"))

# diagnostics <- function(simulatedResponse, observedResponse, diagnostic_fns = NULL, nominal_coverage = 0.9) {
#   
#   if (nominal_coverage <= 0 | nominal_coverage >= 1) stop("nominal_coverage should be between 0 and 1")
#   
#   # Deprecation coercion
#   Z <- observedResponse
#   MC_samples <- simulatedResponse
#   
#   # Define alpha from the nominal coverage
#   alpha <- 1 - nominal_coverage
#   quantiles <- c(alpha/2, 1 - alpha/2)   # symmetric prediction intervals
#   
#   Zhat       <- rowMeans(MC_samples)
#   Zhat_lower <- apply(MC_samples, 1, quantile, alpha/2)
#   Zhat_upper <- apply(MC_samples, 1, quantile, 1 - alpha/2)
#   
#   # Check that the number of predictions match the number of observed data
#   if (length(Z) != length(Zhat)) {
#     stop("The number of observations, observedResponse, does not match the number of predictions, rowMeans(simulatedResponse).")
#   }
#   
#   # Compute default diagnostics
#   df <- data.frame(
#     RMSE          = sqrt(mean((Zhat - Z)^2)),
#     MAE           = mean(abs(Zhat - Z)),
#     CRPS          = mean(crps_sample(y = Z, dat = MC_samples)),
#     Coverage      = mean(Zhat_lower < Z & Z < Zhat_upper),
#     intervalScore = mean(.intervalScore(Z = Z, l = Zhat_lower, u = Zhat_upper, a = alpha))
#   )
#   
#   if (!is.null(diagnostic_fns)) {
#     if (!is.list(diagnostic_fns)) stop("diagnostic_fns should be a named list")
#     if (is.null(names(diagnostic_fns))) stop("diagnostic_fns should be a named list")
#     if (!all(sapply(diagnostic_fns, class) == "function")) stop("Each element of diagnostic_fns should be a function")
#     
#     l  <- list(Z = Z, Zhat = Zhat, Zhat_lower = Zhat_lower, Zhat_upper = Zhat_upper, MC_samples = MC_samples)
#     x  <- lapply(diagnostic_fns, function(f) f(l)) %>% as.data.frame
#     df <- cbind(df, x)
#   }
#   
#   return(df)
# }
# 
# .intervalScore <- function(Z, l, u, a = 0.1) {
#   (u - l) + (2 / a) * (l - Z) * (Z < l) + (2 / a) * (Z - u) * (Z > u) 
# }

# diagnostics(simulatedResponse, observedResponse, nominal_coverage = 0.9)
# diagnostics(simulatedResponse, observedResponse, diagnostic_fns = list(MAD = function(l) mad(l$Zhat - l$Z)), nominal_coverage = 0.9)

