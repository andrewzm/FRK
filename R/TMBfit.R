#' Fitting stage of non-Gaussian FRK.
#'
#' Performs model fitting using \code{TMB}.
#' Prepares an object of class \code{SRE} for the prediction stage
#' (which is performed via the function \code{\link{FRKTMB_pred}}).
#'
#' @param M An object of class \code{SRE}.
#' The slots of an \code{SRE} object particularly relevant to \code{FRKTMB_fit} are
#' \describe{
#'   \item{\code{K_type}}{A string indicating the desired formulation of the eta prior.
#'   Permitted values are "precision" and "block-exponential".}
#'   \item{\code{response}}{A string indicating the assumed distribution of the response variable.}
#'   \item{\code{link}}{A string indicating the assumed link function.}
#'   \item{\code{taper}}{A positve numeric indicating the strength of the covariance tapering
#'   (only applicable if \code{K_type == "covariance"}).}
#' }
#' Furthermore in the case of binomial or negative-binomial response variables,
#' the \code{data} slot of \code{M} must contain a column named \code{k}
#' which contains the 'known-constant' parameters for each observation
#' (the number of trials for binomial data, or the target number of successes
#' for negative-binomial data).
#' @return This function updates the following slots of \code{M}:
#' \describe{
#'   \item{Q}{An estimate of the joint precision matrix of all random effects in the model.}
#'   \item{estimates}{A named list containing estimates of all fixed effects and parameters, and predictions of all random effects.}
#'   \item{log_likelihood}{The log-likelihood evaluated at the optimal estimates.}
#' }
#' This function also makes the slots \code{K_type}, \code{response}, and \code{link} lower case.
#' @seealso \code{\link{FRKTMB_pred}}
#' @export
#' @useDynLib FRK
.FRKTMB_fit <- function(M) {

  ## Strings that must be lower-case
  M@K_type    <- tolower(M@K_type)
  M@response  <- tolower(M@response)
  M@link      <- tolower(M@link)

  # ------ Data preparation ------

  data_params_init <- .TMB_prep(M = M,
                                K_type = M@K_type,
                                response = M@response,
                                link = M@link,
                                taper = M@taper,
                                k_Z = M@data$k_Z)


  # ------- Model Compilation -------

  obj <- MakeADFun(data = data_params_init$data,
                   parameters = data_params_init$parameters,
                   random = c("eta", "xi_O"),
                   DLL = "FRK")


  # ------ Fitting and Parameter Estimates/Random Effect Predictions ------

  ## Fit the model
  fit <- nlminb(obj$par, obj$fn, obj$gr)

  ## Log-likeihood
  log_likelihood <- -obj$fn() # could also use - fit$objective

  ## Extract parameter and random effect estimates
  par <- obj$env$last.par
  estimates <- split(par, names(par)) # convert to named list object


  # ------ Joint precision/covariance matrix of random effects ------

  ## The precision matrix provided by sdreport() is the joint precision matrix
  ## of the fixed parameters AND the random effects.
  ## However, we only want the block corresponding to the random effects.
  ## Note, that the ordering of the precision matrix follows the original
  ## order supplied to TMB i.e. unique(names(obj$env$par)).

  report <- sdreport(obj, getJointPrecision = TRUE, skip.delta.method = TRUE)

  ## Joint precision matrix of fixed AND random effects
  Q <- report$jointPrecision

  # ## For cases in which the dispersion parameter does not change, remove it from Q.
  # if (response %in% c("gaussian", "poisson", "bernoulli", "binomial", "negative-binomial")) {
  #   idx <- which(rownames(Q) == "logphi")
  #   Q <- Q[-idx, -idx]
  # }

  ## The parameters and fixed effects should not be included as random variables,
  ## so we remove them from the precision matrix BEFORE inverting.
  m <- length(M@Z)
  r <- M@basis@n
  npar <- nrow(Q) - (r + m) # number of parameters/fixed effects
  Q    <- Q[-(1:npar), -(1:npar)]


  ## Update the slots of M
  ## NOTE: I had to convert to sparseMatrix as the SRE slot required value 
  M@sigma2fshat <- unname(exp(estimates$logsigma2xi))
  M@alphahat <- as(estimates$beta, "sparseMatrix")
  M@mu_eta <- as(estimates$eta, "sparseMatrix")
  M@mu_xi_O <- as(estimates$xi_O, "sparseMatrix")
  M@Q_eta_xi <- Q                   
  
  ## Not sure if we need to provide kappa and rho (the parameters of the prior precision matrix)?
  ## Also not sure if we need to provide phi (the dispersion paramter), as it is not used beyond this point.
  ## The log-likelihood should be obtained from loglik() function.
  M@log_likelihood <- log_likelihood # log-likelihood evaluated at optimal estimates

  return(M)
}


#' TMB data and parameter preparation.
#'
#' Prepares data and parameter/fixed-effects/random-effects for input to \code{TMB}.
#'
#' @inheritParams .cov_tap
#' @inheritParams .link_fn
#' @param M An object of class \code{SRE}.
#' @param K_type A string indicating whether to use a prior covariance matrix
#' or prior precision matrix for the random effect weights eta.
#' @param response A string indicating the assumed distribution of the response variable.
#' @param k_Z An integer vector. Relevant only to the binomial and negative-binomial cases.
#' For negative-binomial data, the ith element of \code{k_Z} indicates the number of failures until the experiment was stopped for the ith osbervation.
#' For binomial data, the ith element of \code{k_Z} indicates the number of trials for the ith observation.
#' @return A list object containing:
#' \describe{
#'   \item{data}{The data.}
#'   \item{parameters}{The initialised parameters/fixed-effects/random-effects.}
#' }
.TMB_prep <- function (M,
                       K_type = c("precision", "block-exponential"),
                       response = c("gaussian", "poisson", "bernoulli", "gamma",
                                    "inverse-gaussian", "negative-binomial", "binomial"),
                       link = c("log", "identity", "logit", "probit", "cloglog", "reciprocal", "reciprocal-squared"),
                       taper,
                       k_Z) {

  if(is.null(k_Z)){
    k_Z <- -2 # some arbitrary number to keep TMB happy
  }

  nres        <- max(M@basis@df$res)      # Number of resolutions
  r           <- M@basis@n                # Number of basis functions in total.
  N           <- nrow(M@S0)               # Number of BAUs

  Z   <- as.vector(M@Z)
  m   <- length(Z)        # Number of data points AFTER BINNING
  X   <- as.matrix(M@X)
  S   <- as(M@S, "sparseMatrix")

  ## Common to all
  data    <- list(Z = Z, X = X, S = S,
                  ri = as.vector(table(M@basis@df$res)), # Size of each "block": number of basis functions for each resolution
                  sigma2e  = M@Ve[1, 1],
                  response = response, link = link,
                  k_Z = k_Z)

  ## Data which depend on K_type
  if (K_type == "block-exponential") {

    TaperValues  <- .cov_tap(M@D_basis, taper = taper)
    data$D       <- Matrix::bdiag(M@D_basis)
    data$alpha   <- TaperValues$alpha
    data$nnz_tap <- TaperValues$nnz_tap

  } else if (K_type == "precision") {

    temp <- .sparse_Q_block_diag(M@basis@df, kappa = 0, rho = 1)
    R <- as(temp$Q, "dgTMatrix")
    data$row_indices = R@i
    data$col_indices = R@j
    data$x = R@x
    data$nnz = temp$nnz

  }

  ## Parameter and random effect initialisations.
  g <- .link_fn(link)

  ## Create altered data to avoids the problems of applying g() to Z.
  Z0 <- Z
  if (response == "gaussian" & link %in% c("log", "square-root")) {
    Z0[Z <= 0] <- 1      # set non-positive elements of Z to 1
  } else if (response %in% c("poisson", "gamma", "inverse-gaussian") & link == "log") {
    Z0[Z == 0] <- 0.1    # add a small quantity to zero elements
  } else if (response == "negative-binomial" & link %in% c("logit", "probit", "cloglog")) {
    Z0[Z == 0]   <- 0.1
  } else if (response == "negative-binomial" & link == "log") {
    Z0[Z == 0]   <- 0.1
  } else if (response == "binomial" & link %in% c("logit", "probit", "cloglog")) {
    Z0[Z == 0]   <- 0.1
    Z0[Z == k_Z] <- k_Z[Z == k_Z] - 0.1
    Z0 <- Z0/k_Z
  } else if (response == "bernoulli" & link %in% c("logit", "probit", "cloglog")) {
    Z0 <- Z + 0.05 * (Z == 0) - 0.05 * (Z == 1)
  }

  ## Transformed data
  if (response == "binomial" & link %in% c("logit", "probit", "cloglog")) {
    p0   <- Z0 / k_Z
    Z0_t <- g(p0)
  } else if (response == "negative-binomial" & link %in% c("logit", "probit", "cloglog")) {
    p0   <- k_Z / (k_Z + Z0)
    Z0_t <- g(p0)
  } else if (response == "negative-binomial" & link == "log") {
    Z0_t <- g(Z0 / k_Z)
  } else {
    Z0_t <- g(Z0)
  }

  ## Parameter and random effect initialisations.
  parameters <- list()
  parameters$beta         <- coef(lm(Z0_t ~ 1))
  parameters$logsigma2xi  <- log(var(Z0_t))
  parameters$logphi       <- log(data$sigma2e)

  if (K_type == "block-exponential") {
    parameters$logsigma2        <- log(exp(parameters$logsigma2xi) * (0.1)^(0:(nres - 1)))
    parameters$logtau           <- log((1 / 3)^(1:nres))

    ## Tapered prior covariance matrix (required for eta initialisation)
    KInit <- .K_tap_matrix(M@D_basis,
                           alpha = data$alpha,
                           sigma2 =  exp(parameters$logsigma2),
                           tau = exp(parameters$logtau))

  } else if (K_type == "precision") {
    parameters$logrho           <- log(exp(parameters$logsigma2xi) * (0.1)^(0:(nres - 1)))
    parameters$logkappa         <- log((1 / 3)^(1:nres))

    ## Prior covariance matrix (required for eta initialisation)
    KInit <- solve(as(R, "sparseMatrix") + Matrix::sparseMatrix(i = 1:r, j = 1:r, x = 0.5))

  }

  ## Initialise the random effects
  ## THIS SOLVE IS OF AN rxr MATRIX: FOR NRES = 4, THIS WILL TAKE A LONG TIME.
  ## FOR AN INITIALISATION, IT WOULD BE NICE TO FIND SOMETHING THAT DOESN'T INVOLVE
  ## A MATRIX INVERSE.
  temp            <- data$sigma2e+exp(parameters$logsigma2xi)
  parameters$eta  <- as.vector((1 / temp) * solve(Matrix::t(S) %*% S / temp + solve(KInit) ) %*% (Matrix::t(S)  %*% (Z0_t - X %*% parameters$beta)))
  parameters$xi_O <- as.vector((exp(parameters$logsigma2xi) / temp) * (Z0_t - X %*% parameters$beta - S %*% parameters$eta))


  return(list(data = data, parameters = parameters))
}




#' Covariance tapering based on distances.
#'
#' Computes the covariance tapering parameters \eqn{\alpha} (which are dependent
#' on resolution) and the number of non-zeros in each block of the tapered
#' covariance matrix K_tap.
#'
#' \code{taper} determines how strong the covariance tapering is; the ith taper parameter
#' \eqn{\alpha_i} is equal to \code{taper[i]} * \code{minDist[i]}, where
#' \code{minDist[i]} is the minimum distance between basis functions at the
#' \code{i}th resolution.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param taper The strength of the taper (either a vector or a single number).
#'
#' @return A list containing:
#' \describe{
#'   \item{alpha}{A vector of taper parameters.}
#'   \item{nnz_tap}{A vector containing the number of non-zeros in each block.
#' of the tapered prior covariance matrix, K_tap.}
#' }
#' @seealso \code{\link{.K_matrix}}
.cov_tap <- function(D_matrices, taper = 8){

  ri <- sapply(D_matrices, nrow) # No. of basis functions at each res
  nres <- length(ri)

  ## Minimum distance between neighbouring basis functions.
  ## (add a large number to the diagonal, which would otherwise be 0)
  minDist <- vector()
  for(i in 1:nres) minDist[i] <- min(D_matrices[[i]] + 10^5 * diag(ri[i]))

  ## Taper parameters
  alpha   <- taper * minDist

  ## Number of non-zeros in tapered covariance matrix
  nnz_tap <- 0
  for (i in 1:nres) nnz_tap <- nnz_tap + sum(D_matrices[[i]] < alpha[i])

  return(list("alpha" = alpha, "nnz_tap" = nnz_tap))
}



#' K matrix, the un-tapered prior covariance matrix.
#'
#' Construct the prior covariance matrix of the random effects (NOT-TAPERED).
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution
#' (see \code{dist_matrices}).
#' @param sigma2 The variance components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @param tau The correlation components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @return The block-diagonal (class \code{dgCmatrix}) prior covariance matrix, K.
#' @seealso \code{\link{dist_matrices}}
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



#' \eqn{K_\alpha} matrix, the taper matrix.
#'
#' Constructs \eqn{K_\alpha}, the taper matrix, using the Spherical taper.
#'
#' Denoting by \eqn{d(s,s^*)} the distance between two spatial locations,
#' the Spherical taper is given by:
#' \deqn{C_\alpha(s, s^*) = (1-\frac{d(s,s^*)}{\alpha})^2_{+}  (1+\frac{d(s,s^*)}{2\alpha}), }
#' where \eqn{x_+ = max(x, 0)}.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution
#' (see \code{dist_matrices}).
#' @param alpha The taper parameters (which vary based on resolution).
#' @return A taper matrix, which is block-diagonal and of class \code{dgCmatrix}.
.K_alpha_matrix <- function(D_matrices, alpha){

  spherical_taper <- function(matrix, alpha) {
    apply(matrix, c(1,2), function(h){(1 + h / (2 * alpha)) * (max(1 - h / alpha, 0))^2})
  }

  if (class(D_matrices) == "list") {
    K_alpha <- mapply(spherical_taper, matrix = D_matrices, alpha = alpha, SIMPLIFY = FALSE)
  } else {
    K_alpha <- spherical_taper(matrix = D_matrices, alpha = alpha)
  }

  K_alpha <- Matrix::bdiag(K_alpha)

  return(K_alpha)
}



#' K_tap, the tapered prior covariance matrix.
#'
#' Construct K_tap, the \emph{tapered} prior covariance matrix, using the Spherical taper.
#'
#' @inheritParams .K_matrix
#' @inheritParams .K_alpha_matrix
#' @return The \emph{tapered} prior covariance matrix, which is
#' block-diagonal and of class \code{dgCmatrix}.
.K_tap_matrix <- function(D_matrices, alpha, sigma2, tau){

  K <- .K_matrix(D_matrices, sigma2, tau)
  K_alpha <- .K_alpha_matrix(D_matrices, alpha)

  K_tap <- K * K_alpha

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
#' @return If \code{x} and \code{y} are single numbers, then the function
#' returns 1 if \code{x} and \code{y} are equal (within \code{tol}), and 0
#' otherwise. However matrices may also be passed, in which case the function
#' returns a matrix of equal size with elementwise comparisons.
#' @examples
#' .equal_within_tol(2, 2 + 0.00000000000001)
#' .equal_within_tol(2, 2 + 0.1)
#'
#' A <- matrix(1:4, ncol = 2, byrow = TRUE)
#' B <- matrix(c(1:3, 5), ncol = 2, byrow = TRUE)
#' .equal_within_tol(A, B)
#' .equal_within_tol(A, 3)
.equal_within_tol <- function(x, y, tol = 1e-8) {
  return(1 * (abs(x - y) < tol))
}



#' Neighbour matrix.
#'
#' Creates a matrix \eqn{A} with elements \eqn{A_{i, j}} equal to 1 if basis
#' functions i and j (i not equal to j) are neighbours and 0 otherwise.
#' The diagonal elements of \eqn{A} (i.e. \eqn{A_{i, i}}) indicate the total
#' number of neighbours associated with basis function i.
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
#' @return A "neighbour" matrix with element (i, j), for i not equal to j,
#' equal to 1 if basis functions i and j are neighbours, and 0 otherwise.
#' Diagonal elements indicate number of neighbours for that basis function.
.neighbour_matrix <- function(df, loc1 = "loc1", loc2 = "loc2") {

  ## difference in x-coordinate
  diff_x <- outer(df[, loc1], df[, loc1], "-")
  abs_diff_x <- abs(diff_x)

  ## difference in y-coordinate
  diff_y <- outer(df[, loc2], df[, loc2], "-")
  abs_diff_y <- abs(diff_y)

  ## Minimum x and y distances
  x_min <- sort(unique(abs_diff_x[1, ]))[2]
  y_min <- sort(unique(abs_diff_y[1, ]))[2]

  ## 1. Find the "vertical" neighbours
  ## Two conditions which must be met for basis functions to be neighours
  condition1 <- .equal_within_tol(abs_diff_y, y_min)
  condition2 <- .equal_within_tol(diff_x, 0)
  vertical_neighbours <- condition1 & condition2

  ## 2. Find the "horizontal" neighbours
  ## Two conditions which must be met for basis functions to be neighours
  condition1 <- .equal_within_tol(abs_diff_x, x_min)
  condition2 <- .equal_within_tol(diff_y, 0)
  horizontal_neighbours <- condition1 & condition2

  ## Full neighbour matrix (with zeros along the diagonal)
  A <- vertical_neighbours + horizontal_neighbours

  ## number of neighbours for each basis function
  nn <- rowSums(A)

  ## Add the number of neighbours to the diagonal
  diag(A) <- nn

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
#' @return A list containing the sparse block-diagonal precision matrix (Q) of class "dgCMatrix", and the number of non-zero elements (nnz) at each resolution.
#' @seealso \code{\link{.sparse_Q}}, \code{\link{.neighbour_matrix}}
.sparse_Q_block_diag <- function(df, kappa, rho) {

  nres <- max(df$res)
  if (length(kappa) == 1) kappa <- rep(kappa, nres)
  if (length(rho) == 1) rho <- rep(rho, nres)

  ## Construct the blocks
  Q_matrices  <- list()
  nnz <- c()
  for (i in 1:nres) {
    A_i <- .neighbour_matrix(df[df$res == i, ])
    Q_matrices[[i]] <- .sparse_Q(A = A_i,
                                 kappa = kappa[i],
                                 rho = rho[i])
    nnz[i] <- Matrix::nnzero(Q_matrices[[i]])
  }

  ## Block diagonal
  Q <- Matrix::bdiag(Q_matrices)

  return(list(Q = Q, nnz = nnz))
}
