#' Fitting stage of non-Gaussian FRK.
#'
#' Performs model fitting using \code{TMB}.
#' Prepares an object of class \code{SRE} for the prediction stage
#' (performed internally with \code{\link{.FRKTMB_pred}}).
#'
#' @param M An object of class \code{SRE}.
#' The \code{SRE} object slots particularly relevant to \code{FRKTMB_fit} are
#' \describe{
#'   \item{\code{K_type}}{A string indicating the desired formulation of the prior variance/precision matrix of eta.}
#'   \item{\code{response}}{A string indicating the assumed distribution of the response variable.}
#'   \item{\code{link}}{A string indicating the desired link function.}
#'   \item{\code{taper}}{A positve numeric indicating the strength of the covariance tapering (only applicable if \code{K_type == "covariance"}).}
#' }
#' Furthermore in the case of binomial or negative-binomial response variables,
#' the \code{data} slot of \code{M} must contain a column named \code{k}
#' which contains the 'known-constant' parameters for each observation
#' (the number of trials for binomial data, or the target number of successes for negative-binomial data).
#' @param optimiser the optimising function used for model fitting when \code{method = 'TMB'} (default is \code{nlminb}). Users may pass in a function object or a string corresponding to a named function. Optional parameters may be passed to \code{optimiser} via \code{...}. The only requirement of \code{optimiser} is that the first three arguments correspond to the initial parameters, the obejctive function, and the gradient, respectively (note that this may be achieved by rearranging the order of arguments before passing into \code{optimiser}) 
#' @param ... other parameters passed on to \code{auto_basis} and \code{auto_BAUs} when calling \code{FRK}, or the user specified \code{optimiser} function when calling \code{FRK} or \code{SRE.fit}
#' @return This function updates the following slots of \code{M}:
#' \describe{
#'   \item{Q_eta_xi}{An estimate of the joint precision matrix of all random effects in the model (the random weights \eqn{\eta} and observed fine-scale variation \eqn{\xi_O}).}
#'   \item{mu_eta}{Posterior expectation of the random weights \eqn{\eta}.}
#'   \item{mu_xi_O}{Posterior expectation of the observed-fine scale variation \eqn{\xi_O}.}
#'   \item{alphahat}{Estimate of the fixed-effects.}
#'   \item{sigma2fshat}{Estimate of the variance parameter of the fine-scale variation.}
#'   \item{phi}{Estimate of the dispersion parameter (only for applicable response distributions).}
#'   \item{log_likelihood}{The log-likelihood of the model evaluated at the final parameter estimates. Can be obtained by calling loglik(M).}
#' }
.FRKTMB_fit <- function(M, optimiser, ...) {

  ## Data and parameter preparation for TMB
  data_params_init <- .TMB_prep(M)

  ## TMB model compilation
  obj <- MakeADFun(data = data_params_init$data,
                   parameters = data_params_init$parameters,
                   random = c("eta", "xi_O"),
                   DLL = "FRK")
  
  ## View the sparsity pattern.
  ## Note that this can be done before model fitting.
  # tmp <- obj$env$spHess(random = TRUE)
  # image(tmp)
  # length(tmp@x) # number of non-zeros

  ## the following means we want to print every parameter passed to obj$fn.
  obj$env$tracepar <- TRUE
  
  # ---- Model fitting ----
  
  ## The optimiser should have arguments: start, objective, gradient. 
  ## The remaining arguments can be whatever.
  fit <- optimiser(obj$par, obj$fn, obj$gr, ...)
    
  ## Log-likeihood
  log_likelihood <- -obj$fn() # could also use -fit$objective
  
  ## Add the C(Z, phi) terms for one-parameter exponential families
  ## (We do this here rather than in the C++ template because it only needs to 
  ## be done once, as C(Z, phi) depends only on Z for one-parameter exponential 
  ## families)
  Z <- data_params_init$data$Z
  k_Z <- data_params_init$data$k_Z
  if(M@response == "poisson") {
    cZphi = -lfactorial(Z)
  } else if (M@response == "negative-binomial") {
    cZphi   = lfactorial(Z + k_Z - 1.0) - lfactorial(Z) - lfactorial(k_Z - 1.0)
  } else if (M@response == "binomial") {
    cZphi   = lfactorial(k_Z) - lfactorial(Z) - lfactorial(k_Z - Z)
  } 
  
  ## Extract parameter and random effect estimates
  par <- obj$env$last.par.best
  estimates <- split(par, names(par)) # convert to named list object


  # ---- Joint precision/covariance matrix of random effects ----

  ## TMB treats all parameters (fixed effects, variance components, 
  ## and random effects) as random quantities, and so the joint precision 
  ## matrix obtained using sdreport(obj, getJointPrecision = TRUE) contains 
  ## the precision matrix for fixed and random effects. 
  ## However, we assume the regression parameters and variance components are 
  ## fixed effects, NOT random quantities, and so they should not have a 
  ## randomness associated to them. We overcome this by considering the fixed 
  ## effects and parameters as random during the fitting process, and then 
  ## post-fitting we condition on theta = \hat{theta}, the ML estimate.
  ## By conditioning on \hat{theta}, we can consider the precision matrix of  
  ## the random effects in isolation.

  ## Precision matrix of the random effects only
  Q <- obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
  
  ## Update the slots of M
  ## Convert to Matrix as these SRE slots require class "Matrix"
  M@alphahat <- as(estimates$alpha, "Matrix")
  M@mu_eta   <- as(estimates$eta, "Matrix")
  M@mu_xi_O  <- as(estimates$xi_O, "Matrix")
  
  M@sigma2fshat <- exp(estimates$logsigma2fs)
  M@Q_eta_xi <- Q          
  M@phi <- unname(exp(estimates$logphi))
  
  ## log-likelihood evaluated at optimal estimates. After fitting, can be obtained from loglik(M).
  M@log_likelihood <- log_likelihood 

  ## Posterior variance and precision matrix of eta random effects
  r <- nbasis(M)
  M@Q_eta <- Q[1:r, 1:r] # Don't need it, but this matrix is easily obtained, so provide anyway
  M@S_eta <- Matrix()    # Don't need this matrix so set it to NA. It is also not easily obtained because to obtain it we need to invert the joint precision matrix Q_eta_xi, which may be very large.
  
  ## For TMB, we do not need to compute these matrices; return an empty matrix.
  M@Khat_inv <- Matrix()
  M@Khat <- Matrix()
  
  return(M)
}


#' TMB data and parameter preparation.
#'
#' Prepares data and parameter/fixed-effects/random-effects for input to \code{TMB}.
#'
#' @param M An object of class \code{SRE}.
#' @return A list object containing:
#' \describe{
#'   \item{data}{The data.}
#'   \item{parameters}{The initialised parameters/fixed-effects/random-effects.}
#' }
.TMB_prep <- function (M) {

  k_Z  <- as.vector(M@k)
  nres <- max(M@basis@df$res)      # Number of resolutions
  r    <- M@basis@n                # Number of basis functions in total (space x time).
  N    <- nrow(M@S0)               # Number of BAUs
  Z    <- as.vector(M@Z)           # Binned data
  m    <- length(Z)                # Number of data points AFTER BINNING
  X    <- as.matrix(M@X)          

  ## Common to all
  data    <- list(Z = Z, X = X, S = M@S,
                  K_type = M@K_type,
                  response = M@response, 
                  link = M@link,
                  k_Z = k_Z, 
                  temporal = as.integer(is(M@basis,"TensorP_Basis")))
  

  ## The location of the stored spatial basis functions depend on whether the 
  ## data is spatio-temporal or not. Here, we also define the number of temporal
  ## basis function locations (equal to 1 if in a spatial setting).
  if (data$temporal == TRUE) {
    spatial_dist_matrix <- M@D_basis[[1]]
    spatial_basis <- M@basis@Basis1  
    data$r_t <- M@basis@Basis2@n
  } else {
    spatial_dist_matrix <- M@D_basis
    spatial_basis <- M@basis 
    data$r_t <- 1
  }
  
  data$r_si <- as.vector(table(spatial_basis@df$res))

  ## Data which depend on K_type:
  ## Provide dummy data, and only change the ones we actually need.
  data$beta <- data$nnz <- data$row_indices <- data$col_indices <- -1       
  data$x <- data$n_r <- data$n_c  <- -1 

  if (M@K_type == "block-exponential") {
    tmp         <- .cov_tap(spatial_dist_matrix, taper = M@taper)
    data$beta   <- tmp$beta 
    R            <- as(tmp$D_tap, "dgTMatrix")
    data$nnz         <- tmp$nnz 
    data$row_indices <- R@i
    data$col_indices <- R@j
    data$x           <- R@x
    
  } else if (M@K_type == "neighbour") {
    tmp <- .sparse_Q_block_diag(spatial_basis@df, kappa = 0, rho = 1)
    R <- as(tmp$Q, "dgTMatrix")
    data$nnz         <- tmp$nnz
    data$row_indices <- R@i
    data$col_indices <- R@j
    data$x           <- R@x
    
  } else if (M@K_type == "separable") {
    for (i in unique(spatial_basis@df$res)) {
      tmp <- spatial_basis@df[spatial_basis@df$res == i, ]
      data$n_r[i] <- length(unique(tmp$loc1))
      data$n_c[i] <- length(unique(tmp$loc2))
    }
  } 
  
  
  # ---- Parameter and random effect initialisations. ----

  ## Create the relevant link functions. When a probability parameter is 
  ## present in a model and the link-function is appropriate for modelling 
  ## probabilities (i.e., the link function maps to [0, 1]), we may use 
  ## hierarchical linking to first link the probability parameter to the 
  ## Gaussian Y-scale, and then the probability parameter to the conditional 
  ## mean at the data scale. In other situations, we simply map from Y to the mean.
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    h     <- .link_fn(kind = "prob_to_Y", link = M@link)
    f     <- .link_fn(kind = "mu_to_prob", response = M@response)
  } else {
    g     <- .link_fn(kind = "mu_to_Y", link = M@link) 
  }

  ## Create altered data to avoid the problems of applying g() to Z.
  ## This altered data is used only during the initialisation stage.
  Z0 <- Z
  if (M@response == "gaussian" & M@link %in% c("log", "square-root")) {
    Z0[Z <= 0] <- 0.1      
  } else if (M@response %in% c("poisson", "gamma", "inverse-gaussian") & M@link == "log") {
    Z0[Z == 0] <- 0.1    
  } else if (M@response == "negative-binomial" & M@link %in% c("logit", "probit", "cloglog", "log")) {
    Z0[Z == 0]   <- 0.1
  } else if (M@response == "binomial" & M@link %in% c("logit", "probit", "cloglog")) {
    Z0[Z == 0]   <- 0.1
    Z0[Z == k_Z] <- k_Z[Z == k_Z] - 0.1
  } else if (M@response == "bernoulli" & M@link %in% c("logit", "probit", "cloglog")) {
    Z0 <- Z + 0.05 * (Z == 0) - 0.05 * (Z == 1)
  } else if (M@link %in% c("inverse-squared", "inverse")) {
    Z0 <- Z + 0.05 * (Z == 0)
  } 

  ## Transformed data: convert from data scale to Gaussian Y-scale.
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    p0 <- f(Z0, k_Z) # Convert from response (mean) scale to probability scale
    Z0_t <- h(p0)    # Convert from probability scale to Gauussian Y scale
  } else if (M@response == "negative-binomial" & M@link %in% c("log", "square-root")) {
    Z0_t <- g(Z0 / k_Z) # Convert from response (mean) scale to Gaussian Y-scale (accounting for k_Z)
  } else {
    Z0_t <- g(Z0) # Convert from response (mean) scale to Gaussian Y-scale
  }

  ## Parameter and random effect initialisations.
  parameters <- list()
  tmp       <- as.data.frame(as.matrix(M@X))
  tmp$Z0_t  <- Z0_t
  parameters$alpha         <- coef(lm(update(formula(M@f), Z0_t ~ .), tmp)) 
  parameters$logsigma2fs  <- log(var(Z0_t))

  
  ## Dispersion parameter depends on response; some require it to be 1. 
  if (M@response %in% c("poisson", "bernoulli", "binomial", "negative-binomial")) {
    parameters$logphi <- log(1)
  } else if (M@response == "gaussian") {
    parameters$logphi <- log(M@Ve[1, 1]) # FIXME: may have to change this from simply getting the first element of Ve
  } else {
    ## we may not estimate sigma2e when the data is non-Gaussian, so use the 
    ## variance of the data as our estimate of the dispersion parameter.
    parameters$logphi <- log(var(Z0))
  }
  
  parameters$logsigma2      <- log(exp(parameters$logsigma2fs) * (0.1)^(0:(nres - 1)))
  parameters$logtau         <- log((1 / 3)^(1:nres))
  parameters$logsigma2_t    <- log(5)  
  parameters$frho_t         <- 0

  
  
  if (M@K_type == "separable") {
    ## Separability means we have twice as many spatial basis function variance 
    ## components. So, just replicate the already defined parameters. 
    parameters$logsigma2 <- rep(parameters$logsigma2, 2)
    parameters$logtau <- rep(parameters$logtau, 2)
    parameters$logdelta <- log(1)
  } else if (M@K_type == "precision_exp") {
    ## Precision exp requires one extra parameter
    parameters$logdelta <- rnorm(length(data$r_si))
  } else {
    parameters$logdelta <- log(1)
  }

  ## Initialise the latent random effects
  ## FIXME: Should the following be in terms of M@basis@df? Perhaps we could just initialise at one time point, and use the same initialisation for all times e.g. with spatial_basis@df and rep.
  ## FIXME: Add a check if nres == 4. If so, then it is wise to avoid a solve() here. 
  ## FIXME: perhaps change solve() to chol2inv(chol()). However, this requires positive definite matrix.
  ## FIXME: maybe change this intialisation to be only in terms of Qinit. Currently leave as is because the theory I wrote is in terms of Kinit and the block-exponential.
  tmp  <- exp(parameters$logphi) + exp(parameters$logsigma2fs) # note that in the Gaussian case, we have phi = sigma2e
  if (M@K_type == "block-exponential") {
    KInit <- .K_tap_matrix(M@D_basis,
                           beta = data$beta,
                           sigma2 =  exp(parameters$logsigma2),
                           tau = exp(parameters$logtau))
    parameters$eta  <- as.vector((1 / tmp) * solve(Matrix::t(M@S) %*% M@S / tmp + solve(KInit) ) %*% (Matrix::t(M@S)  %*% (Z0_t - X %*% parameters$alpha)))
    parameters$xi_O <- as.vector((exp(parameters$logsigma2fs) / tmp) * (Z0_t - X %*% parameters$alpha - M@S %*% parameters$eta))
  } else {
    kappa <- exp(-parameters$logsigma2)
    rho <- exp(parameters$logtau)
    QInit <- .sparse_Q_block_diag(M@basis@df, kappa = 0.5, rho = 1)$Q
    if (nres >= 4) { # Avoid full inverse if we have four resolutions (in which case r approx 10,000)
      
      ## Matrix we need to invert:
      mat <- Matrix::t(M@S) %*% M@S / tmp + QInit
      
      ## Permuted Cholesky factor
      mat_L <- sparseinv::cholPermute(Q = mat) 
      
      ## Sparse-inverse-subset (a proxy for the inverse)
      mat_inv <- sparseinv::Takahashi_Davis(Q = mat,
                                            cholQp = mat_L$Qpermchol,
                                            P = mat_L$P)
  
      parameters$eta  <- as.vector((1 / tmp) * mat_inv %*% Matrix::t(M@S)  %*% (Z0_t - X %*% parameters$alpha))
      parameters$xi_O <- as.vector((exp(parameters$logsigma2fs) / tmp) * (Z0_t - X %*% parameters$alpha - M@S %*% parameters$eta))  
    } else {
      parameters$eta  <- as.vector((1 / tmp) * solve(Matrix::t(M@S) %*% M@S / tmp + QInit ) %*% (Matrix::t(M@S)  %*% (Z0_t - X %*% parameters$alpha)))
      parameters$xi_O <- as.vector((exp(parameters$logsigma2fs) / tmp) * (Z0_t - X %*% parameters$alpha - M@S %*% parameters$eta))  
    }
    
  }
    
  return(list(data = data, parameters = parameters))
}




#' Covariance tapering based on distances.
#'
#' Computes the covariance tapering parameters \eqn{\beta} (which are dependent
#' on resolution), number of non-zeros in each block of the tapered
#' covariance matrix K_tap, and the tapered distance matrix whereby some distances
#' have been set to zero post tapering (although the remaining non-zero distances 
#' are unchanged).
#'
#' \code{taper} determines how strong the covariance tapering is; the ith taper parameter
#' \eqn{\beta_i} is equal to \code{taper[i]} * \code{minDist[i]}, where
#' \code{minDist[i]} is the minimum distance between basis functions at the
#' \code{i}th resolution.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param taper The strength of the taper (either a vector or a single number).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{A vector of taper parameters.}
#'   \item{nnz}{A vector containing the number of non-zeros in each block of the 
#'   tapered prior covariance matrix, K_tap.}
#'   \item{D_tap}{A sparse block-diagonal matrix containing the distances with 
#'   some distances set to zero post tapering. }
#' }
#' @seealso \code{\link{.K_matrix}}
.cov_tap <- function(D_matrices, taper = 8){

  ri <- sapply(D_matrices, nrow) # No. of basis functions at each res
  nres <- length(ri)

  ## Minimum distance between neighbouring basis functions.
  ## (add a large number to the diagonal, which would otherwise be 0)
  ## FIXME: This approach is flawed. What if we have a very skinny rectangle for the domain of interest?
  minDist <- vector()
  for(i in 1:nres) minDist[i] <- min(D_matrices[[i]] + 10^8 * diag(ri[i]))

  ## Taper parameters
  beta   <- taper * minDist
  
  ## Construct D matrix with elements set to zero after tapering
  ## FIXME: I could also apply the taper here, whereby we compute the spherical 
  ## taper and multiply the distances by the taper. This would have the advtange 
  ## that we no longer need to construct the taper within C++. However, it is 
  ## probably simpler to leave as is.
  D_tap <- list()
  nnz <- c()
  for (i in 1:nres) {
    indices    <- D_matrices[[i]] < beta[i]
    D_tap[[i]] <- as(D_matrices[[i]] * indices, "sparseMatrix")
    
    # Add explicit zeros
    D_tap[[i]] <- D_tap[[i]] + sparseMatrix(i = 1:ri[i], j = 1:ri[i], x = 0) 
    
    # Number of non-zeros in tapered covariance matrix at each resolution
    nnz[i] <- length(D_tap[[i]]@x) # Note that this approach DOES count explicit zeros
  }

  D_tap <- Matrix::bdiag(D_tap)

  return(list("beta" = beta, "nnz" = nnz, "D_tap" = D_tap))
}



#' K matrix, the un-tapered prior covariance matrix.
#'
#' Construct the prior covariance matrix of the random effects (NOT-TAPERED).
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param sigma2 The variance components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @param tau The correlation components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @return The block-diagonal (class \code{dgCmatrix}) prior covariance matrix, K.
.K_matrix <- function(D_matrices, sigma2, tau){

  exp_cov_function <- function(dist_matrix, sigma2, tau) {
    sigma2 * exp(- dist_matrix / tau)
  }

  if (class(D_matrices) == "list") {

    ## Construct the individual blocks of the K-matrix
    K_matrices <- mapply(exp_cov_function,
                         dist_matrix = D_matrices,
                         sigma2 = sigma2, tau = tau,
                         SIMPLIFY = FALSE)

    ## Construct the full K-matrix
    K <- Matrix::bdiag(K_matrices)

  } else {
    K <- exp_cov_function(D_matrices, sigma2 = sigma2, tau  = tau)
  }

  return(K)
}



#' \eqn{K_\beta} matrix, the taper matrix.
#'
#' Constructs \eqn{K_\beta}, the taper matrix, using the Spherical taper.
#'
#' Denoting by \eqn{d(s,s^*)} the distance between two spatial locations,
#' the Spherical taper is given by:
#' \deqn{C_\beta(s, s^*) = (1-\frac{d(s,s^*)}{\beta})^2_{+}  (1+\frac{d(s,s^*)}{2\beta}), }
#' where \eqn{x_+ = max(x, 0)}.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param beta The taper parameters (which vary based on resolution).
#' @return A taper matrix, which is block-diagonal and of class \code{dgCmatrix}.
.K_beta_matrix <- function(D_matrices, beta){

  spherical_taper <- function(matrix, beta) {
    apply(matrix, c(1,2), function(h){(1 + h / (2 * beta)) * (max(1 - h / beta, 0))^2})
  }

  if (class(D_matrices) == "list") {
    K_beta <- mapply(spherical_taper, matrix = D_matrices, beta = beta, SIMPLIFY = FALSE)
  } else {
    K_beta <- spherical_taper(matrix = D_matrices, beta = beta)
  }

  K_beta <- Matrix::bdiag(K_beta)

  return(K_beta)
}



#' K_tap, the tapered prior covariance matrix.
#'
#' Construct K_tap, the \emph{tapered} prior covariance matrix, using the Spherical taper.
#'
#' @inheritParams .K_matrix
#' @inheritParams .K_beta_matrix
#' @return The \emph{tapered} prior covariance matrix, which is
#' block-diagonal and of class \code{dgCmatrix}.
.K_tap_matrix <- function(D_matrices, beta, sigma2, tau){

  K <- .K_matrix(D_matrices, sigma2, tau)
  K_beta <- .K_beta_matrix(D_matrices, beta)

  K_tap <- K * K_beta

  return(K_tap)
}



#' Test equality within a given tolerance.
#'
#' Tests equality between objects \code{x} and \code{y} (within a given tolerance, \code{tol}).
#' The primary purpose of this function is to avoid deeming objects to be unequal if they differ
#' by some very tiny amount due to floating point inaccuracies.
#' Of particular note is that the function can accept a matrix argument for \code{x} and a single numeric
#' for \code{y}, and the output will be a matrix with elements 1 and 0 if elements of \code{x} are equal to
#' \code{y} or not, respectively; i.e., it does elementwise comparisons.
#'
#' @param x \code{R} object.
#' @param y \code{R} object we wish to compare to \code{x}.
#' @param tol Tolerance.
#' @return If \code{x} and \code{y} are single numbers, then the function
#' returns 1 if \code{x} and \code{y} are equal (within \code{tol}), and 0
#' otherwise. However matrices may also be passed, in which case the function
#' returns a matrix of equal size with elementwise comparisons.
## Examples:
## Note that .equal_within_tol is not exported, so it cannot be used in the 
## usual roxygen @examples format. We could use FRK:::.equal_within_tol, 
## but this produces a warning. For unexported functions, it is best to leave
## examples as regular comments.
# .equal_within_tol(2, 2 + 0.00000000000001)
# .equal_within_tol(2, 2 + 0.1)
# 
# A <- matrix(1:4, ncol = 2, byrow = TRUE)
# B <- matrix(c(1:3, 5), ncol = 2, byrow = TRUE)
# .equal_within_tol(A, B)
# .equal_within_tol(A, 3)
.equal_within_tol <- function(x, y, tol = 1e-8) {
  return(1 * (abs(x - y) < tol))
}


## Determine if range of vector is FP 0.
## (thise function tests whether all elements of a vector are equal, with a tolerance) 
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

## Determine if the basis functions are in a regular rectangular grid
## Inspired by: https://stackoverflow.com/a/32916893
.test_regular_rect_grid <- function(x, y) {
  
  a <- sort(unique(x))
  b <- sort(unique(y))
  
  if(
    .zero_range(diff(a)) &
    .zero_range(diff(b)) &
    (length(a) * length(b)) == length(x)
  ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}


#' Neighbour matrix.
#'
#' Creates a matrix \eqn{A} with elements \eqn{A_{i, j}} equal to 1 if basis
#' functions i and j (i not equal to j) are first order neighbours, 1/2 if they are second order neighbours, and so on.
#' Neighbours with a larger order than specified by \code{order} have 0 in this matrix.
#' The diagonal elements of \eqn{A} (i.e. \eqn{A_{i, i}}) indicate the row sums (if 
#' the \code{order == 1}, then it is the totalnumber of first order neighbours associated with basis function i).
#'
#' This function is only designed for basis functions
#' at \emph{one resolution}. It also assumes the basis functions are in a
#' regularly-spaced lattice; the shape of the lattice is not important, however
#' there \emph{must} be constant spacing between basis functions in a given
#' direction (horizontal and vertical spacing can be different).
#'
#' @seealso \code{\link{.sparse_Q}}
#'
#' @param df A dataframe containing spatial coordinates.
#' @param loc1 A string indicating the name of the column storing the x-coordinate.
#' @param loc2 A string indicating the name of the column storing the y-coordinate.
#' @param order If order == 1, only first order neighbours are considered. If order == 2, second order neighbours are also considered, and so on.
#' @param diag_neighbours Indicates whether to consider the diagonal neighbours. If FALSE (the default), only the horizontal and vertical neighbours are considered.
#' @return A "neighbour" matrix with element (i, j), for i not equal to j,
#' equal to 1/l if basis functions i and j are lth order neighbours (provided \code{l <= order}), and 0 otherwise.
#' Diagonal elements indicate the row sums.
.neighbour_matrix <- function(df, loc1 = "loc1", loc2 = "loc2", order = 1, diag_neighbours = FALSE) {
  
  A <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  
  ## absolute difference in x- and y-coordinate
  abs_diff_x <- abs(outer(df[, loc1], df[, loc1], "-"))
  abs_diff_y <- abs(outer(df[, loc2], df[, loc2], "-"))
  
  ## Vectors containing all x and y distances. 
  ## Note that the first elements in each vector is 0. 
  ## Note also that we only use 1 row from the abs_diff matrices (this helps 
  ## to prevents problems with unique() and floating point accuracy and is a 
  ## bit faster)
  x_distances <- sort(unique(abs_diff_x[1, ]))
  y_distances <- sort(unique(abs_diff_y[1, ]))
  
  for (i in 1:order) { ## Order will usually be 1
    x_min <- x_distances[i + 1] # x distance we are focused on for this order
    y_min <- y_distances[i + 1] # y distance we are focused on for this order
    
    ## Find the "horizontal" neighbours
    ## This produces a matrix with (i, j)th entry equal to 1 if i and j are 
    ## horizontal neighbours, and zero otherwise.
    ## Similarly, find the "vertical" neighbours.
    horizontal_neighbours <- .equal_within_tol(abs_diff_x, x_min) & .equal_within_tol(abs_diff_y, 0) 
    vertical_neighbours   <- .equal_within_tol(abs_diff_y, y_min) & .equal_within_tol(abs_diff_x, 0)
    
    ## Consider the diagonal neighbours, if specified.
    if (diag_neighbours == TRUE) {
      ## Two conditions which must be met for basis functions to be diagonal neighours
      diagonal_neighbours <- .equal_within_tol(abs_diff_y, y_min) & .equal_within_tol(abs_diff_x, x_min)
      ## Update A
      A <- A + 1/i * diagonal_neighbours
    }
    
    ## Update neighbour matrix (with zeros along the diagonal)
    ## We weight the neighbours by their order. 
    A <- A + 1/i * vertical_neighbours + 1/i * horizontal_neighbours
  }

  ## Add the sums of each row to the diagonal (required for use in the 
  ## precision matrix computation later)
  diag(A) <- rowSums(A)
  
  return(A)
}



#' Sparse precision matrix.
#'
#' Creates a sparse precision matrix \eqn{Q} with off diagonal elements equal
#' to -1 if the basis functions are neighbours, and zero otherwise.
#' The diagonal elements are equal to the number of neighbours for that basis
#' function, plus some amount given by \code{kappa}.
#'
#' @param A A "neighbour" matrix with element (i, j), for i not equal to j,
#' equal to 1 if basis functions i and j are neighbours, and 0 otherwise,
#' diagonal elements indicating the number of neighbours for that basis function.
#' @param kappa Quantity to add to the diagonal elements. This must be positive if Q is to be positive definite.
#' @param rho Quantity to multiply matrix by. This must be positive if Q is to be positive definite.
#' @return A sparse precision matrix of class \code{dgCMatrix}.
#' @seealso \code{\link{.neighbour_matrix}}, \code{\link{.sparse_Q_block_diag}}
.sparse_Q <- function(A, kappa, rho) {

  Q <- -A
  diag(Q) <- diag(A) + kappa
  Q <- rho * Q
  Q <- as(Q, "sparseMatrix")

  return(Q)
}



#' Block-diagonal sparse precision matrix.
#'
#' Creates a block-diagonal sparse precision matrix, where the blocks are created
#' using \code{sparse_Q}.
#'
#' @inheritParams .sparse_Q
#' @inheritParams .neighbour_matrix
#' @param df A dataframe containing the spatial coordinates (named "loc1" and "loc2") and a blocking column (named "res").
#' @return A list containing the sparse block-diagonal precision matrix (Q) of class "dgCMatrix", and the number of non-zero elements (nnz) at each resolution.
#' @seealso \code{\link{.sparse_Q}}, \code{\link{.neighbour_matrix}}
.sparse_Q_block_diag <- function(df, kappa, rho, order = 1, diag_neighbours = FALSE) {

  
  nres <- max(df$res)
  if (length(kappa) == 1) kappa <- rep(kappa, nres)
  if (length(rho) == 1) rho <- rep(rho, nres)

  ## Construct the blocks
  Q_matrices  <- list()
  nnz <- c()
  for (i in 1:nres) {
    A_i <- .neighbour_matrix(df[df$res == i, ], order = order, diag_neighbours = diag_neighbours)
    Q_matrices[[i]] <- .sparse_Q(A = A_i,
                                 kappa = kappa[i],
                                 rho = rho[i])
    nnz[i] <- Matrix::nnzero(Q_matrices[[i]]) # note that nnzero does not count explicit zeros
  }

  ## Block diagonal
  Q <- Matrix::bdiag(Q_matrices)

  return(list(Q = Q, nnz = nnz))
}
