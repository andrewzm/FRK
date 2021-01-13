#' Prediction stage of non-Gaussian FRK.
#'
#' @inheritParams .Y_var
#' @inheritParams .MC_sampler
#' @inheritParams .concat_percentiles_to_df
#' @param newdata object of class \code{SpatialPoylgons} indicating the regions over which prediction will be carried out. The BAUs are used if this option is not specified
#' @param CP Polygon prediction matrix
#' @param predict_BAUs Logical indicating whether or not we are predicting over the BAUs
#' @param pred_time vector of time indices at which prediction will be carried out. All time points are used if this option is not specified
#' @param type A character string (possibly vector) indicating the quantities for which predictions and prediction uncertainty is desired. If \code{"link"} is in \code{type}, the latent \eqn{Y} process is included; If \code{"mean"} is in \code{type}, the conditional mean \eqn{\mu} is included (and the probability parameter if applicable); If \code{"response"} is in \code{type}, the response variable \eqn{Z} is included. Note that any combination of these character strings can be provided. For example, if \code{type = c("link", "response")}, then predictions of the latent \eqn{Y} process and the response variable \eqn{Z} are provided
#' @param kriging A string indicating the kind of kriging (\code{"simple"} or \code{"universal"})
#' @return A list object containing:
#' \describe{
#'   \item{newdata}{An object of class \code{newdata}, with predictions and prediction uncertainty at each prediction location of the latent \eqn{Y} process, the conditional mean of the data \eqn{\mu}, the probability of success parameter \eqn{\pi} (if applicable), and the response variable \eqn{Z}}
#'   \item{MC}{A list with each element being an \code{N * n_MC} matrix of Monte Carlo samples of the quantities specified by \code{type} (some combination of \eqn{Y}, \eqn{\mu}, \eqn{p} (if applicable), and \eqn{Z}) at each prediction location}
#' }
#' Note that for all link functions other than the log- and identity-link functions, the predictions and prediction uncertainty of \eqn{\mu} contained in \code{newdata} are computed using the Monte Carlo samples contained in \code{MC}.
#' When the log- or identity-link functions are used, the expectation and variance of the \eqn{\mu} may be computed exactly.
.FRKTMB_pred <- function(M, newdata, CP, predict_BAUs, pred_time, type, n_MC, 
                         obs_fs, k, percentiles, cred_mass, kriging) {
  
  
  ## FIXME: predict over only the observed BAUs needed for newdata locations
  
  # ## If the user has specified which polygons he wants we can remove the ones we don't need
  # ## We only need those BAUs that are influenced by observations and prediction locations
  # ## For ST, use all BAUs, as it gets complicated
  # if(!predict_BAUs & is(newdata,"Spatial")) {
  #   
  #   ## The needed BAUs are the nonzero column CP
  #   needed_BAUs <- as(CP,"dgTMatrix")@j
  # 
  #   ## Filter the BAUs and the matrices
  #   ## (Note that we do not update the SRE object so this is safe to do)
  #   M@BAUs <- M@BAUs[needed_BAUs, ]
  #   CP <- CP[, needed_BAUs]
  #   # CZ <- CZ[, needed_BAUs]
  #   M@S0 <- M@S0[needed_BAUs, ]
  # }

  # ---- Create objects needed thoughout the function ----
  
  ## Id of observed BAUs
  obsidx <- observed_BAUs(M) # FIXME: Use this everywhere, so could just make it into a slot of SRE object
  
  ## The covariate design matrix, X (at the BAU level i.e. for all BAUs)
  X <- as(.extract_BAU_X_matrix(formula = M@f, BAUs = M@BAUs), "matrix")

  
  # ---- Compute the Cholesky factor of the permuted precision matrix ----
  
  ## Number of fixed and random effects
  p <- length(M@alphahat)
  mstar <- length(obsidx)
  r <- ncol(M@S0)
  s <- r + mstar * M@include_fs
  
  ## Permuted Cholesky factor. If we are doing universal kriging, keep the joint precision 
  ## matrix of the fixed and random effects. Otherwise, if we are doing simple kriging, 
  ## use only the random effect block of the precision matrix.
  if (kriging == "universal") {
    Q_joint <- M@Q_eta_xi
  } else if (kriging == "simple") {
    Q_joint <- M@Q_eta_xi[-(1:p), -(1:p)]
  }
  Q_L <- sparseinv::cholPermute(Q = Q_joint)
  
  
  # ------ Latent process Y prediction and Uncertainty ------
  
  ## Note that this does NOT depend on the response of Z.
    
  ## Posterior expectation E(Y|Z) at each prediction location (i.e., at each BAU).
  ## large- and medium-scale variation terms:
  ## FIXME: change from vector and leave as a matrix form. Then we can remove other as.vector() calls too
  # p_Y <- as.vector(X %*% M@alphahat + M@S0 %*% M@mu_eta)  # see Equation (2.1.2)

  ## Add posterior estimate of xi_O at observed BAUs
  # p_Y[obsidx]   <-  p_Y[obsidx] + as.vector(M@mu_xi)
  
  ## Posterior variance of Y at each prediction location.
  ## Note that MSPE(E(Y|Z), Y) is approximated by var(Y|Z).
  
  
  # Analytic <- TRUE
  # MSPE_Y  <- .Y_var(M = M, Q_joint = Q_joint, Q_L = Q_L, obsidx = obsidx, X = X, kriging = kriging)
  
  
  
  # ------ Monte Carlo sampling ------
  
  ## Generate Monte Carlo samples at all BAUs
  MC <- .MC_sampler(M = M, X = X, type = type, obs_fs = obs_fs, 
                    n_MC = n_MC, k = k, Q_L = Q_L, obsidx = obsidx, 
                    predict_BAUs = predict_BAUs, CP = CP, kriging = kriging)
  
  ## We do not allow aggregation of the Y-process when predicting over arbitrary polygons
  if(!predict_BAUs)
    MC$Y_samples <- NULL
  
  ## Remove other quantities if the user is not interested in them
  if(!("link" %in% type)) 
    MC$Y_samples <- NULL
  
  if(!("mean" %in% type)) 
    MC$mu_samples <- MC$prob_samples <- NULL
  

  # ------ Create Prediction data ------
  
  ## Produce prediction and RMSPE matrices. 
  ## The columns are the quantity of interest (Y, mu, prob, or Z), and the rows are prediction locations.
  predictions <- sapply(MC, rowMeans)
  RMSPE <- sapply(MC, apply, 1, sd)

  ## If we are predicting over BAUs, newdata is NULL, so set it to the BAUs.
  ## Note that the BAUs must be a Spatial*DataFrame, so coercion is unnecessary, but we will do it 
  ## to be safe.
  if (predict_BAUs) {
    newdata <- M@BAUs 
    newdata$dummy_blah123 <- rep(1, length(newdata))
    newdata$dummy_blah123 <- NULL
  }
    
  ## Now update newdata with the predictions, RMSPE, and percentiles. 
  ## (See https://datascience.stackexchange.com/a/8924 for a description of what this gsub is doing.)
  QOI <- gsub("_.*", "", names(MC)) # Quantities Of Interest
  
  ## Predictions and RMSPE
  colnames(predictions) <- paste0("p_", QOI)
  colnames(RMSPE) <- paste0("RMSPE_", QOI)
  newdata@data <- cbind(newdata@data, predictions, RMSPE)
  
  ## Use p_Y (computed with estimates from TMB) as the predictor for Y, rather than the noisy MC estimates
  # if("link" %in% type & predict_BAUs) {
  #   newdata$p_Y <- p_Y
  #   newdata$RMSPE_Y <- sqrt(MSPE_Y)
  # }
  
  ## Percentiles and HPD interval bounds
  newdata@data <- .concat_percentiles_to_df(data = newdata@data, MC = MC, 
                                            cred_mass = cred_mass, percentiles = percentiles)


  # ## Previous code: The main difference is that here we use p_Y and MSPE_Y 
  # ## rather than the Y_samples in MC, and that for some link functions (namely
  # ## the identity and lok link functions), we can analytically compute the mean. 
  # ## It would be good to retain these optimisations, but I will come back to this later.
  # ## I think this would probably need to go in the MC sampling function.. maybe not though.
  # if (predict_BAUs) {
  # 
  #   ## Conditional mean of the data, mu
  #   ## If a log- or identity-link function is used, then expectations and variance
  #   ## of the conditional mean may be evaluated analytically.
  #   ## Otherwise, use the Monte Carlo simulations for prediction.
  #   if (M@link == "log" & M@response != "negative-binomial") {
  #     p_mu         <- exp(p_Y + MSPE_Y / 2)
  #     RMSPE_mu     <- sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  #   } else if (M@link == "log" & M@response == "negative-binomial") {
  #     p_mu         <- k * exp(p_Y + MSPE_Y / 2)
  #     RMSPE_mu     <- k * sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  #   } else if (M@link == "identity") {
  #     p_mu <- p_Y
  #     RMSPE_mu <- sqrt(MSPE_Y)
  #   } else {
  #     ## USE THE SAMPLES.
  #   }
  # 
  # }
  
  # ---- Add a column indicating the time point in space-time setting ----
  
  ## FIXME: how does pred_time come into play? I think pred_time gives the time 
  ## indices to predict over. Perhaps we need to simply subset the final 
  # predictions based on pred_time
  if (is(M@basis,"TensorP_Basis")) {
    newdata$t <- M@BAUs@data$t
  }
  
  ## It is convenient to have the spatial coordinates in the @data slot of the
  ## return newdata object. Only add those coordinates not already in the data.
  tmp <- which(!(colnames(coordinates(newdata)) %in% names(newdata@data)))
  if (length(tmp))
    newdata@data <- cbind(newdata@data, coordinates(newdata)[, tmp]) 

    
  ## Return the predictions, and the MC samples at either the BAUs (if we are 
  ## predicting over BAUs) or over the user specified arbitrary polygons.
  return(list(newdata = newdata, MC = MC))
}



#' Posterior variance of the latent Y process.
#'
#' Computes the variance of the latent process \eqn{Y} at every BAU. Note that MSPE(E(Y|Z), Y) is approximated by var(Y|Z).
#'
#' To compute the prediction uncertainty of Y we require the joint
#' covariance matrix of the random effects \eqn{(\eta', \xi_O')'}. \code{TMB}
#' provides an approximation of the joint \emph{precision} matrix of
#' \eqn{(\eta', \xi_O')'}, which we must invert to obtain the approximate
#' covariance matrix. However, due to the potentially very large number of
#' random effects (the number of observations \eqn{m} is not restricted),
#' in practice we compute the sparse-inverse-subset of the joint precision matrix,
#' \emph{not} the true joint covariance matrix.
#'
#' Note that as we are using E(\eqn{Y|Z}) to predict \eqn{Y}, the posterior variance acts as an approximation of the mean-squared prediction error (see pg. 72 of Honours thesis).
#' 
#' @param M An object of class SRE.
#' @param Q_L A list containing the Cholesky factor of the permuted precision matrix (stored as \code{Q$Qpermchol}) and the associated permutationmatrix (stored as \code{Q_L$P}).
#' @param obsidx Vector containing the observed locations.
#' @return A vector of the posterior variance of Y at every BAU. 
.Y_var <- function(M, Q_joint, Q_L, obsidx, X, kriging){
  
  r <- ncol(M@S0)
  mstar <- length(obsidx)
  
  ## Number of fixed and random effects
  p <- length(M@alphahat)
  s <- r + mstar * M@include_fs
  
  ## number of spatial and temporal BAUs
  if (is(M@basis,"TensorP_Basis")) {
    ns <- length(M@BAUs@sp)
    nt <- length(unique(M@BAUs@endTime))
  } else {
    ns <- length(M@BAUs)
  }
  
  
  # ---- Sparse-inverse-subset of Q (acting as a proxy for the true covariance matrix) ----

  
  if (r + mstar < 4000) {
    Sigma <- chol2inv(chol(M@Q_eta_xi))
  } else {
    ## Sparse-inverse-subset of random effects (eta and xi_O)
    ## (a proxy for the covariance matrix)
    Sigma <- sparseinv::Takahashi_Davis(Q = Q_joint,
                                        cholQp = Q_L$Qpermchol,
                                        P = Q_L$P)
  }


  if (kriging == "universal") {
    Sigma_alpha <- Sigma[1:p, 1:p, drop = FALSE]
    Sigma_random <- Sigma[-(1:p), -(1:p), drop = FALSE]
    Cov_alpha_eta <- Sigma[1:p, (p+1):(p+r), drop = FALSE]
    Cov_alpha_xi <- Sigma[1:p, (p + r +1):(p+r+mstar), drop = FALSE]
  } else if (kriging == "simple") {
    Sigma_random <- Sigma
  }
  
  Sigma_eta <- Sigma_random[1:r, 1:r]
  if (M@include_fs) {
    Sigma_xi    <- Sigma_random[(r + 1):(r + mstar), (r + 1):(r + mstar)]
    Cov_eta_xi  <- Sigma_random[1:r, (r + 1):(r + mstar)] # Covariances between xi_O and eta
  }

  
  # ----- Uncertainty: Posterior variance of Y at each BAU ------
  
  ## To extract the variances of eta|Z, we need diag(S0 %*% Sigma_eta %*% t(S0)).
  ## Also, to extract the covariance terms, we need: diag(S %*% COV_{eta, xi}).
  ## This in very inefficient to do directly, it much better to use the identity:
  ##      diag(AB) = (A*B')1
  

  ## Add common terms for both observed and unobserved locations:
  vY <- as.vector( (M@S0 %*% Sigma_eta * M@S0) %*% rep(1, r) )
  
  if(kriging == "universal") {
    
    ## Variance due to alpha
    vY <- vY + as.vector( (X %*% Sigma_alpha * X) %*% rep(1, p) )
    
    ## Covariance terms between alpha and eta
    cov_alpha_eta <- as.vector( (X %*% Cov_alpha_eta * M@S0) %*% rep(1, r) )
    vY <- vY + 2 * cov_alpha_eta
  }

  
  if (M@include_fs) {
    
    ## UNOBSERVED locations
    
    ## simply add the estimate of sigma2fs to the variance.
    ## If we have a unique fine-scale variance at each spatial BAU (spatio-temporal 
    ## case only), add the sigma2fs associated with that BAU.
    if (M@fs_by_spatial_BAU) {
      unobsidx <- unobserved_BAUs(M)
      spatial_BAU_id <- ((unobsidx - 1) %% ns) + 1
      vY[unobsidx] <- vY[unobsidx] + M@sigma2fshat[spatial_BAU_id]
    } else {
      vY[-obsidx] <- vY[-obsidx] + M@sigma2fshat
    }
    
    ## OBSERVED locations
    
    ## add both var(xi_O|Z) and cov(xi_O, eta | Z)
    vY[obsidx] <- vY[obsidx] + diag(Sigma_xi) + 2 * (M@S_O * t(Cov_eta_xi)) %*% rep(1, r)

    ## Add covariance between alpha and xi_O if we are using universal kriging
    if(kriging == "universal") 
      vY[obsidx] <- vY[obsidx] + 2 * (M@X_O * t(Cov_alpha_xi)) %*% rep(1, p)
    
  }

  
  # ---- Output ----
  
  ## Return variance of Y
  return(vY)
}



#' Monte Carlo sampling of the conditional mean of the data (a function of the
#' latent process Y).
#'
#' Computes a Monte Carlo sample of \eqn{Y}, the conditional mean of the data
#' \eqn{\mu = \psi(Y)} (which is a deterministic function of Y), the response variable \eqn{Z}, and, for response-link
#' combinations to which it is applicable, the probability of success parameter
#' p. It does so for every BAU location. 
#' 
#' For negative-binomial and binomial data, the \code{BAUs} slot of the \code{SRE} object must contain a field \code{k}, which is the known constant parameter for each BAU.
#' For negative-binomial data, the ith element of \code{k} indicates the number of failures until the experiment is stopped at the ith BAU.
#' For binomial data, the ith element of \code{k} indicates the number of trials at the ith BAU.
#'
#'
#' @param M An object of class \code{SRE}
#' @param X The design matrix of the covariates at the BAU level (often simply an Nx1 column vector of 1's)
#' @param type A character string (possibly vector) indicating the quantities which are the focus of inference. Note: unlike in the predict() function, \emph{all} computed quantities are returned. That is, the latent \eqn{Y} process samples are always provided; If \code{"mean"} \emph{OR} \code{"response"} is in \code{type}, then the samples of \eqn{Y}, the conditonal mean \eqn{\mu}, and the probability parameter (if applicable) are provided. If \code{"response"} is in \code{type}, the response variable \eqn{Z} samples, and the samples of all other quantities are provided
#' @param n_MC A postive integer indicating the number of MC samples at each location
#' @param obs_fs Logical indicating whether the fine-scale variation is included in the latent Y process. If \code{obs_fs = FALSE} (the default), then the fine-scale variation term \eqn{\xi} is included in the latent \eqn{Y} process. If \code{obs_fs = TRUE}, then the the fine-scale variation terms \eqn{\xi} are removed from the latent Y process; \emph{however}, they are re-introduced for computation of the conditonal mean \eqn{\mu} and response variable \eqn{Z}
#' @param k vector of size parameters at each BAU (applicable only for binomial and negative-binomial data)
#' @param Q_L A list containing the Cholesky factor of the permuted precision matrix (stored as \code{Q$Qpermchol}) and the associated permutationmatrix (stored as \code{Q_L$P})
#' @param obsidx A vector containing the indices of observed locations
#' @return A list containing Monte Carlo samples of various quantites of interest. The list elements are (N x n_MC) matrices, whereby the ith row of each matrix corresponds to \code{n_MC} samples of the given quantity at the ith BAU. The available quantities are:
#' \describe{
#'   \item{Y_samples}{Samples of the latent, Gaussian scale Y process}
#'   \item{mu_samples}{Samples of the conditional mean of the data}
#'   \item{prob_samples}{Samples of the probability of success parameter (only for the relevant response distributions)}
#'   \item{Z_samples}{Samples of the response variable}
#' }
.MC_sampler <- function(M, X, type, n_MC, obs_fs, k, Q_L, obsidx, predict_BAUs, CP, kriging){
  
  MC <- list()              # object we will return holding the quantities of interest (Y, mu, Z)
  N   <- nrow(M@S0)   
  mstar <- length(obsidx)
  r   <- ncol(M@S0)   # Total number of basis functions
  
  ## Number of fixed and random effects
  p <- length(M@alphahat)
  s <- r + mstar * M@include_fs
  
  ## number of spatial and temporal BAUs
  if (is(M@basis,"TensorP_Basis")) {
    ns <- length(M@BAUs@sp)
    nt <- length(unique(M@BAUs@endTime))
  } else {
    ns <- length(M@BAUs)
  }
  
  
  # ---- Generate samples from (eta', xi_O')' ----
  
  ## Must generate samples jointly, as eta and xi_O are correlated.
  
  if (M@include_fs) {
    ## Construct the mean vector of (eta', xi')',
    ## then make an (r + m*) x n_MC matrix whose columns are the mean vector of (eta', xi')'.
    ## Finally, generate (r + m*) x n_MC samples from Gau(0, 1) distribution.
    if (kriging == "universal") {
      mu_eta_xi_O <- c(as.numeric(M@alphahat), as.numeric(M@mu_eta), as.numeric(M@mu_xi))
    } else {
      mu_eta_xi_O <- c(as.numeric(M@mu_eta), as.numeric(M@mu_xi))
    }
    
    mu_eta_xi_O_Matrix  <- matrix(rep(mu_eta_xi_O, times = n_MC), ncol = n_MC)
    
    if (kriging == "universal") {
      z <- matrix(rnorm((p + r + mstar) * n_MC), nrow = p + r + mstar, ncol = n_MC)
    } else {
      z <- matrix(rnorm((r + mstar) * n_MC), nrow = r + mstar, ncol = n_MC)
    }
    
    
    ## Compute the Cholesky factor of  Q (the joint precision matrix of (eta', xi')').
    ## Then, to generate samples from (eta, xi_O), 
    ## use eta_xi = L^{-T} z + mu = U^{-1} z + mu, 
    ## where U upper cholesky factor of Q, so that Q = U'U.
    U <- Matrix::t(Q_L$Qpermchol) # upper Cholesky factor of permuted joint precision matrix M@Q_eta_xi
    # x <- backsolve(U, z)          # x ~ Gau(0, A), where A is the permuted precision matrix i.e. A = P'QP
    x <- solve(U, z)          
    y <- Q_L$P %*% x              # y ~ Gau(0, Q^{-1})
    eta_xi_O  <- as.matrix(y + mu_eta_xi_O_Matrix) # add the mean to y
  
  
    ## Separate the eta and xi samples
    if (kriging == "universal") {
      alpha <- eta_xi_O[1:p, , drop = FALSE]
      eta   <- eta_xi_O[(p + 1):(p + r), ]
      xi_O  <- eta_xi_O[(p + r + 1):(p + r + mstar), ]
    } else {
      eta   <- eta_xi_O[1:r, ]
      xi_O  <- eta_xi_O[(r + 1):(r + mstar), ]
    }
    
  } else {
    
    
    ## Construct the mean vector of eta,
    ## then make an r x n_MC matrix whose columns are the mean vector of eta.
    ## Finally, generate r x n_MC samples from Gau(0, 1) distribution.
    mu_eta         <- as.numeric(M@mu_eta)
    mu_eta_Matrix  <- matrix(rep(mu_eta, times = n_MC), ncol = n_MC)
    z <- matrix(rnorm(r * n_MC), nrow = r, ncol = n_MC)
    
    ## Compute the Cholesky factor of  Q, where Q is the precision matrix of eta.
    ## Then, to generate samples from eta, use eta = L^{-T} z + mu = U^{-1} z + mu, 
    ## where U upper cholesky factor of Q, so that Q = U'U.
    U <- Matrix::t(Q_L$Qpermchol) # upper Cholesky factor of permuted joint precision matrix M@Q_eta_xi
    # x <- backsolve(U, z)          # x ~ Gau(0, A), where A is the permuted precision matrix i.e. A = P'QP
    x <- solve(U, z)
    y <- Q_L$P %*% x              # y ~ Gau(0, Q^{-1})
    eta  <- as.matrix(y + mu_eta_Matrix) # add the mean to y
  }
  

  ## We now have two matrices, eta and xi_O:
  ## row i of eta corresponds to n_MC MC samples of eta_i,
  ## row i of xi_O corresponds to n_MC MC samples of the fine-scale variation at the ith observed location.
  
  # ---- Generate samples from xi_U ----
  
  ## This is straightforward as each element of xi_U is independent of
  ## all other random effects in the model.
  ## All we have to do is make an (N-m) x n_MC matrix of draws from the
  ## Gaussian distribution with mean zero and variance equal to the fine-scale variance.
  if (M@include_fs) {
    
    if (M@fs_by_spatial_BAU) {
      unobsidx <- unobserved_BAUs(M)
      spatial_BAU_id <- ((unobsidx - 1) %% ns) + 1
      sigma2fs_U <- M@sigma2fshat[spatial_BAU_id]
    } else {
      sigma2fs_U <- M@sigma2fshat
    }
    
    xi_U <- rnorm((N - mstar) * n_MC, mean = 0, sd = sqrt(sigma2fs_U)) %>% 
      matrix(nrow = N - mstar, ncol = n_MC)

    xi_samples <- rbind(xi_O, xi_U)
  }
  
  
  
  
  # ---- Construct samples from the latent process Y ----
  
  ## We break the latent process down as: Y = Y_smooth + xi, 
  ## so that we may separate the fine-scale variation. 

  ## Split the covariate design matrix based on observed and unobserved samples
  X_U <- X[-obsidx, ]   # Unobserved fixed effect 'design' matrix
  S_U <- M@S0[-obsidx, ] # Unobserved random effect 'design' matrix

  ## Simulate the smooth Y-process (excluding fs variation) over the observed and unobserved BAUs
  if (kriging == "universal") {
    Y_smooth_O <- M@X_O %*% alpha + M@S_O %*% eta
    Y_smooth_U <- X_U %*% alpha + S_U %*% eta
  } else {
    Y_smooth_O <- M@X_O %*% M@alphahat + M@S_O %*% eta
    Y_smooth_U <- X_U %*% M@alphahat + S_U %*% eta
  }
  
  ## Combine samples
  Y_smooth_samples  <- rbind(Y_smooth_O, Y_smooth_U)
  
  ## Use permutation matrix to get the correct (original) ordering
  unobsidx         <- unobserved_BAUs(M)   # Unobserved BAUs indices
  ids              <- c(obsidx, unobsidx)  # All indices (observed and unobserved)
  P                <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  Y_smooth_samples <- Matrix::t(P) %*% Y_smooth_samples

  # ## Sanity check:
  # N = 10
  # obsidx <- sample(1:N, 6, replace = F)
  # unobsidx <- (1:N)[-obsidx]
  # ids              <- c(obsidx, unobsidx)  # All indices (observed and unobserved)
  # P                <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  # n_MC <- 5
  # (Y_smooth_samples <- rep(c(obsidx, unobsidx), each = n_MC) %>% matrix(nrow = N, byrow = T))
  # Matrix::t(P) %*% Y_smooth_samples
  
  ## Construct the samples from the latent process Y 
  if (M@include_fs) {
    xi_samples  <- Matrix::t(P) %*% xi_samples
    Y_samples   <- Y_smooth_samples + xi_samples
  } else {
    Y_samples <- Y_smooth_samples
  }
    
  Y_samples        <- as.matrix(Y_samples)
  Y_smooth_samples <- as.matrix(Y_smooth_samples)

  ## Outputted Y value depend on obs_fs
  if (obs_fs) {
    MC$Y_samples <- Y_smooth_samples
  } else 
    MC$Y_samples <- Y_samples
    
  
  ## If Y is the ONLY quantity of interest, exit the function.
  if (!("mean" %in% type) & !("response" %in% type)) return(MC) 
  
  
  # ---- Apply inverse-link function to the samples to obtain conditional mean ----
  
  ## Past this point we must have xi in the Y process (the model breaks down otherwise).
  ## In the case of type == "all", we simply export Y_smooth_samples as the samples of Y.

  ## For families with a known constant parameter (binomial, negative-binomial),
  ## zeta() maps the Gaussian scale Y process to the probability parameter p.
  ## Then, we map p to the conditional mean mu via chi().
  ## For all other families, psi() maps Y directly to mu.
  ## The exception is negative-binomial with a log or square-root link, 
  ## in which case we map directly from Y to mu.
  
  ## Note that for all cases other than type == "link", we need to compute the conditional mean samples.
  
  ## Create the relevant link functions.
  ## FIXME: change these names to ginv, hinv, finv
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    zeta    <- .link_fn("Y_to_prob", link = M@link)
    chi     <- .link_fn("prob_to_mu", response = M@response)
  } else {
    psi     <- .link_fn("Y_to_mu", link = M@link) 
  }
  
  ## Create the mu samples (and prob parameter if applicable)
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    prob_samples <- zeta(Y = Y_samples)
    mu_samples   <- chi(p = prob_samples, k = k)
  } else if (M@response == "negative-binomial" & M@link %in% c("log", "square-root")) {
    mu_samples   <- k * psi(Y_samples)
    f            <- .link_fn(kind = "mu_to_prob", response = M@response)
    prob_samples <- f(mu = mu_samples, k = k)
  } else {
    mu_samples <- psi(Y_samples)
  }



  
  # ---- Predicting over arbitrary polygons ----
  
  
  ## FIXME: Should add some checks for this. The user should only provide "mean" or "response" in \code{type}
  ## if they wish to predict over arbitrary polygons other than the BAUs.
  if (!predict_BAUs) 
    mu_samples <- as.matrix(CP %*% mu_samples)

    

  
  ## Output the mean samples. If probability parameter was computed, and 
  ## we are predicting over the BAUs, also output.
  MC$mu_samples <- mu_samples
  if (exists("prob_samples") & predict_BAUs) 
    MC$prob_samples <- prob_samples
  
  
  ## If the response is not a quanitity of interest, exit the function
  if (!("response" %in% type)) return(MC)
 
  
  # ---- Sample the response variable, Z ----

  n <- nrow(CP) * n_MC
  
  if (M@response == "poisson") {
    Z_samples <- rpois(n, lambda = c(t(mu_samples)))
  } else if (M@response == "gaussian") {
    sigma2e <- M@Ve[1, 1] # measurement error standard deviation
    Z_samples <- rnorm(n, mean = c(t(mu_samples)), sd = sqrt(sigma2e))
  } else if (M@response == "bernoulli") {
    Z_samples <- rbinom(n, size = 1, prob = c(t(mu_samples)))
  } else if (M@response == "gamma") {
    theta <- 1 / c(t(mu_samples)) # canonical parameter
    alpha <- 1/M@phi                 # shape parameter
    beta  <- theta * alpha           # rate parameter (1/scale)
    Z_samples <- rgamma(n, shape = alpha, rate = beta)
    Z_samples <- statmod::rinvgauss(n, mean = c(t(mu_samples)), dispersion = M@phi)
  } else if (M@response == "negative-binomial") {
    k_vec <- rep(k, each = n_MC)
    Z_samples <- rnbinom(n, size = k_vec, mu = c(t(mu_samples)))
  } else if (M@response == "binomial") {
    k_vec <- rep(k, each = n_MC)
    theta <- log((c(t(mu_samples))/k_vec) / (1 - (c(t(mu_samples))/k_vec)))
    p <- 1 / (1 + exp(-theta))
    ## NAs will occur if k = 0. Fortunately, if k = 0, we know Z will be 0. 
    ## Hence, simply replace the NA occurences in p with 0.
    p[is.na(p)] <- 0
    Z_samples <- rbinom(n, size = k_vec, prob = p)
  }
  
  ## Convert from a long vector to an n_pred_locs x n_MC matrix
  Z_samples <- matrix(Z_samples, ncol = n_MC, byrow = TRUE)
  
  ## Add Z_samples to list object
  MC$Z_samples <- Z_samples

  return(MC)
}
