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
#' @param optimiser the optimising function used for model fitting when \code{method = 'TMB'} (default is \code{nlminb}). Users may pass in a function object or a string corresponding to a named function. Optional parameters may be passed to \code{optimiser} via \code{...}. The only requirement of \code{optimiser} is that the first three arguments correspond to the initial parameters, the objective function, and the gradient, respectively (note that this may be achieved by rearranging the order of the arguments before passing into \code{optimiser}) 
#' @param known_sigma2fs known value of the fine-scale variance. If \code{NULL} (the default), the fine-scale variance \eqn{\sigma^2_\xi} is estimated as usual. If \code{known_sigma2fs} is not \code{NULL}, the fine-scale variance is fixed to the supplied value; this may be a scalar, or vector of length equal to the number of spatial BAUs (if fs_by_spatial_BAU = TRUE)
#' @param ... other parameters passed on to \code{auto_basis} and \code{auto_BAUs} when calling \code{FRK}, or the user specified \code{optimiser} function when calling \code{FRK} or \code{SRE.fit}
#' @return This function updates the following slots of \code{M}:
#' \describe{
#'   \item{Q_posterior}{An estimate of the joint precision matrix of all random effects in the model (the random weights \eqn{\eta} and fine-scale variation \eqn{\xi}).}
#'   \item{mu_eta}{Posterior expectation of the random weights \eqn{\eta}.}
#'   \item{mu_xi}{Posterior expectation of the fine-scale variation \eqn{\xi}.}
#'   \item{alphahat}{Estimate of the fixed-effects.}
#'   \item{sigma2fshat}{Estimate of the variance parameter of the fine-scale variation.}
#'   \item{phi}{Estimate of the dispersion parameter (only for applicable response distributions).}
#'   \item{log_likelihood}{The log-likelihood of the model evaluated at the final parameter estimates. Can be obtained by calling loglik(M).}
#' }
.FRKTMB_fit <- function(M, optimiser, known_sigma2fs, ...) {
  
  ## Parameter and data preparation for TMB
  parameters <- .TMB_initialise(M)
  data <- .TMB_data_prep(M, sigma2fs_hat = exp(parameters$logsigma2fs))
  
  ## If we are estimating a unique fine-scale variance at each spatial BAU, 
  ## simply replicate sigma2fs ns times. 
  ns <- dim(M@BAUs)[1]
  if (M@fs_by_spatial_BAU) {
    data$sigma2fs_hat <- rep(data$sigma2fs_hat, ns)
    parameters$logsigma2fs <- rep(parameters$logsigma2fs, ns)
  }

  ## Fix sigma2fs to the known value provided by the user (if provided). 
  if (!is.null(known_sigma2fs)) {
    data$fix_sigma2fs <- 1
    data$sigma2fs_hat <- known_sigma2fs
    parameters$logsigma2fs <- log(known_sigma2fs) 
  }
   
  ## Don't want to pass in variance components that are "too small"
  parameters$logsigma2 <- pmax(parameters$logsigma2, -3)
  parameters$logsigma2_t <- pmax(parameters$logsigma2_t, -3)
  parameters$logtau <- pmax(parameters$logtau, -3)
  
  ## TMB model compilation
  obj <- MakeADFun(data = data,
                   parameters = parameters,
                   random = c("random_effects"),
                   DLL = "FRK")

  ## The following means we want to print every parameter passed to obj$fn.
  obj$env$tracepar <- TRUE
  
  # ---- Model fitting ----
  
  ## The optimiser should have arguments: start, objective, gradient. 
  ## The remaining arguments can be whatever.
  fit <- optimiser(obj$par, obj$fn, obj$gr, ...)
  
  ## Log-likeihood (negative of the negative-log-likelihood)
  M@log_likelihood <- -obj$fn() # could also use -fit$objective
  
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

  ## Number of parameters and fixed effects:
  p <- length(obj$par)
  s <- length(estimates$random_effects)
  
  ## We need to use sdreport() if we wish to use the uncertainty of the
  ## parameters and fixed effects, as opposed to using obj$env$spHess(par = obj$env$last.par.best, random = TRUE) 
  ## We will only retain the uncertainty in the fixed effects
  ## (i.e., in alpha), and not the parameters.
  Q_posterior <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  retain_idx  <- rownames(Q_posterior) %in% c("alpha", "random_effects") 
  Q_posterior <- Q_posterior[retain_idx, retain_idx]

  ## Update the slots of M
  ## Convert to Matrix as these SRE slots require class "Matrix"
  r  <- nbasis(M)
  mstar <- ncol(M@C_O)
  M@alphahat <- as(estimates$alpha, "Matrix")
  M@mu_eta   <- as(estimates$random_effects[1:r], "Matrix")
  if (M@include_fs) {
    M@mu_xi  <- as(estimates$random_effects[(r+1):(r + mstar)], "Matrix")
  } else {
    M@mu_xi  <- as(rep(0, mstar), "Matrix")
  }
  
  M@sigma2fshat <- unname(exp(estimates$logsigma2fs))
  M@Q_posterior <- Q_posterior
  M@phi <- unname(exp(estimates$logphi))
  
  return(M)
}



## Initalise the fixed effects, random effects, and parameters for method = 'TMB'
.TMB_initialise <- function(M) {   
  
  # ---- Set-up ----
  
  l <- list() # list of initial values
  
  nres    <- max(M@basis@df$res)  # Number of basis function resolutions
  X_O     <- M@X_O
  S_O     <- M@S_O
  C_O     <- M@C_O
  Z       <- as.vector(M@Z)         
  k_Z     <- M@k_Z       
  k_BAU_O <- M@k_BAU_O
  r       <- M@basis@n   
  m <- nrow(C_O)
  mstar <- ncol(C_O)
  
  ## Create the relevant link functions. When a probability parameter is 
  ## present in a model and the link-function is appropriate for modelling 
  ## probabilities (i.e., the link function maps to [0, 1]), we may use 
  ## hierarchical linking to first link the probability parameter to the 
  ## Gaussian Y-scale, and then the probability parameter to the conditional 
  ## mean at the data scale. In other situations, we simply map from Y to the mean.
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    f     <- .link_fn(kind = "prob_to_Y", link = M@link)
    h     <- .link_fn(kind = "mu_to_prob", response = M@response)
  } else {
    g     <- .link_fn(kind = "mu_to_Y", link = M@link) 
  }
  
  ## Create altered data to avoid the problems of applying g() to Z.
  ## This altered data is used only during the initialisation stage.
  Z0 <- Z
  if (M@link %in% c("log", "square-root")) {
    Z0[Z <= 0] <- 0.1      
  } else if (M@response == "negative-binomial" & M@link %in% c("logit", "probit", "cloglog")) {
    Z0[Z == 0]   <- 0.1
  } else if (M@response == "binomial" & M@link %in% c("logit", "probit", "cloglog")) {
    Z0 <- Z + 0.1 * (Z == 0) - 0.1 * (Z == k_Z)
  } else if (M@response == "bernoulli" & M@link %in% c("logit", "probit", "cloglog")) {
    Z0 <- Z + 0.05 * (Z == 0) - 0.05 * (Z == 1)
  } else if (M@link %in% c("inverse-squared", "inverse")) {
    Z0[Z == 0] <- 0.05
  } 
  
  
  # ---- Estimate mu_Z, mu_O, and Y_O ----
  
  ## First, we estime mu_Z with the (adjusted) data
  mu_Z <- Z0
  
  ## I don't do this anymore, because including the size parameter makes the
  ## elements of mu_O extremely small. 
  # ## Excluding k from the C matrix makes sense when we are using C to 
  # ## aggregate the mean process, but NOT when we are trying to go the 
  # ## other way, because we need to account for k when splitting up mu_Z into mu.
  # if(M@response %in% c("binomial", "negative-binomial")) {
  #   ## if any observations are associated with multiple BAUs, we want to account 
  #   ## for the size parameter. Note that we cannot normalise the rows of this 
  #   ## matrix afterwards, because this could lead to some mu_O[i] > k_BAU_O[i] 
  #   if (any(as(C_O, "dgTMatrix")@j > 1))
  #     C_O@x <- k_BAU_O[as(C_O, "dgTMatrix")@j + 1]
  # }
  
  
  ## Construct mu_O from mu_Z
  if (M@response %in% c("binomial", "negative-binomial")) {
    ## Need to split mu_Z up based on the size parameter associated with each BAU
    mu_O <- numeric(mstar)
    for (x in 1:m) { # for each observation
      ## Find the BAU indices associated with mu_Z[x]
      idx <- which((C_O@i+1) == x)
      for (j in idx) # fill in the corresponding mu_O
        mu_O[j] <- (k_BAU_O[j] / sum(k_BAU_O[idx])) * mu_Z[x]
    }
    ## Sanity check: any(mu_O > k_BAU_O)
  }
  
  else {
    ## Use Moores-Penrose inverse to "solve" C_O mu_O = mu_Z. 
    ## NB: if C_O is m x m*, then Cp is m* x m.
    
    if (suppressWarnings(rankMatrix(C_O, method = 'qr') == min(dim(C_O)))) { 
      ## if C_O is full rank, simple algebraic forms for Moores-Penrose inverse exist
      if (m == mstar) {
        Cp <- solve(C_O)
      } else if (m < mstar) {
        Cp <- t(C_O) %*% solve(C_O %*% t(C_O)) 
      } else if (m > mstar) {
        Cp <- solve(t(C_O) %*% C_O) %*% t(C_O)  
      } 
    } else {
      ## Use rank-deficient computation of Moores-Penrose inverse 
      Cp <- VCA::MPinv(C_O) 
      Cp <- as(Cp, "dgCMatrix") # NB: not sure if VCA::MPinv() takes advantage of sparsity.
      ## Sparsity of C_O and Cp: nnzero(C_O) / (m * mstar), nnzero(Cp) / (m * mstar)
    }
    mu_O <- as.vector(Cp %*% mu_Z)
    ## Sanity check: max(abs(C_O %*% mu_O - mu_Z)) 
    ## sum((C_O %*% mu_O - mu_Z)^2)
    
    ## Development: it may be better to use the QR factorisation so solve the
    ## system C_O mu_O = mu_Z.
    #system.time(tmp <- qr.solve(C_O, mu_Z)) # this takes a long time (~60s)
  }
  
  ## For some link functions, mu_0 = 0 causes NaNs; set these to a small positive value.
  ## The size parameter being 0 also causes NaNs.
  k_BAU_O[k_BAU_O == 0] <- 1
  mu_O <- mu_O + 0.05 * (mu_O == 0) - 0.05 * (mu_O == k_BAU_O)
  
  
  ## Transformed data: convert from data scale to Gaussian Y-scale.
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    Y_O <- f(h(mu_O, k_BAU_O)) 
  } else if (M@response == "negative-binomial" & M@link %in% c("log", "square-root")) {
    Y_O <- g(mu_O / k_BAU_O) 
  } else {
    Y_O <- g(mu_O)
  } 
  
  
  
  # ---- Parameter and random effect initialisations ----
  
  ## Now that we have a crude estimate of Y_O, we may initialise the variance
  ## components and random effects

  ## i. Fixed effects alpha (OLS solution)
  l$alpha <- solve(t(X_O) %*% X_O) %*% t(X_O) %*% Y_O # OLS solution
  
  ## ii. Variance components
  ## Dispersion parameter depends on response; some require it to be 1. 
  if (M@response %in% c("poisson", "bernoulli", "binomial", "negative-binomial")) {
    l$phi <- 1
  } else if (M@response == "gaussian") {
    l$phi <- mean(diag(M@Ve))
  } else {
    ## Use the variance of the data as our estimate of the dispersion parameter.
    ## This will almost certainly be an overestimate, as the mean-variance 
    ## relationship is not considered.
    l$phi <- var(Z0)
  }
  
  l$sigma2  <- var(as.vector(Y_O)) * (0.1)^(0:(nres - 1))
  l$tau     <- (1 / 3)^(1:nres)
  if (M@K_type != "block-exponential") {
    l$sigma2   <- 1 / exp(l$sigma2)
    l$tau      <- 1 / exp(l$tau)
  }
  ## variance components of the basis function random weights
  ## FIXME: Not sure what to do for time yet. Also haven't really thought about when K_type = 'separable'. 
  l$sigma2_t    <- 1
  l$rho_t      <- 0.1
  if (M@K_type == "separable") {
    ## Separability means we have twice as many spatial basis function variance
    ## components. So, just replicate the already defined parameters.
    l$sigma2 <- rep(l$sigma2, 2)
    l$tau <- rep(l$tau, 2)
  }
  

  for (iteration_dummy in 1:5) {
    
    ## iii. Basis function random-effects 
    regularising_weight <- if (!is.null(l$sigma2fs)) l$sigma2fs else l$sigma2[1] 
    
    QInit <- .sparse_Q_block_diag(M@basis@df, 
                                  kappa = exp(l$sigma2), 
                                  rho = exp(l$tau))$Q
    
    ## Matrix we need to invert
    mat <- Matrix::t(S_O) %*% S_O / regularising_weight + QInit 
    
    ## Avoid full inverse if we have too many basis functions
    if (r > 4000) { 
      mat_L <- sparseinv::cholPermute(Q = mat) # Permuted Cholesky factor
      ## Sparse-inverse-subset (a proxy for the inverse)
      mat_inv <- sparseinv::Takahashi_Davis(Q = mat, cholQp = mat_L$Qpermchol, P = mat_L$P)
    } else {
      mat_inv <- solve(mat) 
    }
    
    ## MAP estimate of eta
    l$eta  <- (1 / regularising_weight) * mat_inv %*% Matrix::t(S_O)  %*% (Y_O - X_O %*% l$alpha)
    
    ## iv. Observed fine-scale random effects xi_O
    l$xi_O <- Y_O - X_O %*% l$alpha - S_O %*% l$eta
    
    ## v. Fine-scale variance
    l$sigma2fs <- var(as.vector(l$xi_O)) 
  }

  ## Return list of parameter initialisations for TMB
  transform_minus_one_to_one_inverse <- function(x) -0.5 * log(2 / (x + 1) - 1)
  return(list(
    alpha = as.vector(l$alpha),
    logphi = log(l$phi),
    logsigma2 = log(l$sigma2),
    logtau = log(l$tau),
    logsigma2_t = log(l$sigma2_t),  
    frho_t = transform_minus_one_to_one_inverse(l$rho_t),
    logsigma2fs = log(l$sigma2fs),
    random_effects = c(as.vector(l$eta), if(M@include_fs) as.vector(l$xi_O))
  ))
}


.TMB_data_prep <- function (M, sigma2fs_hat) {

  obsidx <- observed_BAUs(M)       # Indices of observed BAUs
  
  ## measurement error variance (NB: this is a vector)
  sigma2e <- if (M@response == "gaussian") diag(M@Ve) else -1
  
  ## Data passed to TMB which is common to all
  data <- list(Z = as.vector(M@Z),  # Binned data
               X_O = M@X_O, S_O = M@S_O, C_O = M@C_O,
               K_type = M@K_type, response = M@response, link = M@link,
               k_BAU_O = M@k_BAU_O, k_Z = M@k_Z,         
               temporal = as.integer(is(M@basis,"TensorP_Basis")), 
               fs_by_spatial_BAU = M@fs_by_spatial_BAU, sigma2e = sigma2e, 
               BAUs_fs = M@BAUs$fs[obsidx])

  ## Define the number of temporal basis function (r_t), and number of spatial BAUs (ns).
  ns <- dim(M@BAUs)[1]
  if (data$temporal) {
    spatial_dist_matrix <- M@D_basis[[1]]
    spatial_basis <- M@basis@Basis1  
    data$r_t <- M@basis@Basis2@n
  } else {
    spatial_dist_matrix <- M@D_basis
    spatial_basis <- M@basis 
    data$r_t <- 1
  }
  
  data$spatial_BAU_id <-  (obsidx - 1) %% ns
  data$r_si <- as.vector(table(spatial_basis@df$res))

  ## Data which depend on K_type: provide dummy data (can't provide nothing)
  ## and only change the ones we actually need.
  data$beta <- data$nnz <- data$row_indices <- data$col_indices <- 
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
    ## Compute number of basis functions in each row (n_r) and each column (n_c)
    for (i in unique(spatial_basis@df$res)) {
      tmp <- spatial_basis@df[spatial_basis@df$res == i, ]
      data$n_r[i] <- length(unique(tmp$loc1))
      data$n_c[i] <- length(unique(tmp$loc2))
    }
  } 

  ## Create a data entry of sigma2fs_hat (one that will stay constant if we are 
  ## not estimating sigma2fs within TMB)
  data$sigma2fs_hat <- sigma2fs_hat
  ## Only estimate sigma2fs if all the observations are associated with exactly 
  ## one BAU; otherwise, we must fix sigma2fs, or else TMB will explode.
  data$fix_sigma2fs <- as.integer( !all(tabulate(M@Cmat@i + 1) == 1) )
  data$include_fs   <- as.integer(M@include_fs)

  return(data)
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
## (tests whether all elements of a vector are equal, with a tolerance) 
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

## Determine if the basis functions are in a regular rectangular grid
## Inspired by: https://stackoverflow.com/a/32916893
.test_regular_grid <- function(x, y, rectangular) {
  
  a <- sort(unique(x))
  b <- sort(unique(y))
  
  condition <- .zero_range(diff(a)) && .zero_range(diff(b))
  
  if (rectangular)
    condition <- condition & ((length(a) * length(b)) == length(x))
  
  return(condition)
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



## Exctracts the covariate design matrix of the fixed effects at the BAU level.
.extract_BAU_X_matrix <- function (formula, BAUs) {
  ## Retrieve the dependent variable name
  depname <- all.vars(formula)[1]
  
  ## Set the dependent variable in BAUs to something just so that 
  ## .extract.from.formula doesn't throw an error. We are not returning M, 
  ## so we don't need to worry about NULL-ing afterwards.
  BAUs[[depname]] <- 0.1
  
  ## Extract covariates from BAUs
  L <- .extract.from.formula(formula, data = BAUs)
  X <- as(L$X,"Matrix")
  
  return(X)
}