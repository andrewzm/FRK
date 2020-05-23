#' Prediction stage of non-Gaussian FRK.
#'
#' @inheritParams .Y_var
#' @inheritParams .MC_sampler
#' @inheritParams .concat_percentiles_to_df
#' @param type A character string (possibly vector) indicating the quantities for which predictions and prediction uncertainty is desired.
#' If \code{"link"} is in \code{type}, the latent \eqn{Y} process is included; 
#' If \code{"mean"} is in \code{type}, the conditional mean \eqn{\mu} is included (and the probability parameter if applicable);
#' If \code{"response"} is in \code{type}, the response variable \eqn{Z} is included. 
#' Any combination of these character strings can be provided. For example, if \code{type = c("link", "response")}, then predictions of the latent \eqn{Y} process and the response variable \eqn{Z} are provided.
#' @return A list object containing:
#' \describe{
#'   \item{newdata}{A dataframe with predictions and prediction uncertainty at each prediction location of the latent \eqn{Y} process, the conditional mean of the data \eqn{\mu}, the probability of success parameter \eqn{\pi} (if applicable), and the response variable \eqn{Z}. The dataframe also contains percentiles of each term.}
#'   \item{MC}{A list with each element being an \code{N * n_MC} matrix of Monte Carlo samples of the quantities specified by \code{type} (some combination of \eqn{Y}, \eqn{\mu}, \eqn{p} (if applicable), and \eqn{Z}) at each prediction location.}
#'   \item{pred_time}{A list of timings for each step in the prediction stage.}
#' }
#' Note that for all link functions other than the log-link and identity-link, the predictions and prediction uncertainty of \eqn{\mu} contained in \code{newdata} are computed using the Monte Carlo samples contained in \code{MC}.
#' When the log- or identity-link functions are used the expectation and variance of the \eqn{\mu} may be computed exactly.
.FRKTMB_pred <- function(M, type = "mean", n_MC = 400, obs_fs = FALSE, 
                         k = NULL, 
                         percents = c(5, 25, 50, 75, 95)) {
  
  pred_time <- list()
  
  # ---- Create objects needed thoughout the function ----
  
  ## Id of observed BAUs
  obsidx <- apply(M@Cmat, 1, function(x) which(x == 1))
  
  #### The covariate design matrix, X (at the BAU level i.e. for all BAUs)
  
  ## Retrieve the dependent variable name
  depname <- all.vars(M@f)[1]
  
  ## Set the dependent variable in BAUs to something just so that .extract.from.formula doesn't
  ## throw an error.. we will NULL it shortly after
  M@BAUs[[depname]] <- 0.1
  
  ## Extract covariates from BAUs
  L <- .extract.from.formula(M@f, data = M@BAUs)
  X <- as(L$X,"Matrix")
  M@BAUs[[depname]] <- NULL
  
  rm(depname, L)
  
  # ---- Computed the Cholesky factor of the permuted precision matrix ----
  
  ## Permuted Cholesky factor
  pred_time$compute_cholesky_factor <- system.time(
    Q_L <- sparseinv::cholPermute(Q = M@Q_eta_xi)
  )
  

  
  # ------ Latent process Y prediction and Uncertainty ------
  
  ## Note that this does NOT depend on the response of Z.
  
  pred_time$Y_pred_and_uncertainty <- system.time({
    
  #### Posterior expectation E(Y|Z) at each prediction location.
  ## large- and medium-scale variation terms:
  p_Y <- as.vector(X %*% M@alphahat + M@S0 %*% M@mu_eta)

  ## Add posterior estimate of xi_O at observed BAUs
  p_Y[obsidx]   <-  p_Y[obsidx] + as.vector(M@mu_xi_O)
  
  #### Posterior variance of Y at each prediction location.
  ## Note that MSPE(E(Y|Z), Y) is approximated by var(Y|Z).

  MSPE_Y  <- .Y_var(M = M, Q_L = Q_L, obsidx = obsidx) 
  })
  
  # ------ Conditional mean (mu) prediction and uncertainty ------
  
  ## Compute Monte Carlo samples of conditional mean at each location
  pred_time$MC_sample_total <- system.time(
    MC <- .MC_sampler(M = M, X = X, type = type, obs_fs = obs_fs, 
                      n_MC = n_MC, k = k, Q_L = Q_L, obsidx = obsidx)
  ) 

  pred_time$MC_sample_backsolve <- MC$times$backsolve # time of backsolve specifically
  
  # ------ Create Prediction Dataframe ------
  
  ## This dataframe contains x and y, the spatial locations, of every BAU. 
  newdata <- as.data.frame(coordinates(M@BAUs)) 
  rownames(newdata) <- NULL
  
  ## Latent Y-process
  if ("link" %in% type) {
    newdata$p_Y     <- p_Y          
    newdata$RMSPE_Y <- sqrt(MSPE_Y) 
    newdata <- .concat_percentiles_to_df(X = MC$Y_samples, df = newdata, 
                                         name = "Y", percents = percents)  
  }
  
  ## If Y is the ONLY quantity of interest, exit the function.
  if (!("mean" %in% type) && !("response" %in% type)) return(list(newdata = newdata, MC = MC)) 
  
  ## Conditional mean of the data, mu
  ## If a log- or identity-link function is used, then expectations and variance
  ## of the conditional mean may be evaluated analytically.
  ## Otherwise, use the Monte Carlo simulations for prediction.
  if (M@link == "log" & M@response != "negative-binomial") {
    p_mu         <- exp(p_Y + MSPE_Y / 2)
    RMSPE_mu     <- sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  } else if (M@link == "log" & M@response == "negative-binomial") {
    p_mu         <- M@BAUs$k * exp(p_Y + MSPE_Y / 2)
    RMSPE_mu     <- M@BAUs$k * sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  } else if (M@link == "identity") {
    p_mu <- p_Y
    RMSPE_mu <- sqrt(MSPE_Y)
  } else {
    p_mu         <- rowMeans(MC$mu_samples)
    RMSPE_mu     <- sqrt(.rowVars(MC$mu_samples))
  }
  
  ## Output mu (and prob) predictions if it is requested
  if ("mean" %in% type) {
    newdata$p_mu <- p_mu
    newdata$RMSPE_mu <- RMSPE_mu
    newdata <- .concat_percentiles_to_df(X = MC$mu_samples, df = newdata, 
                                         name = "mu", percents = percents)
    
    ## For some response distributions, the probability of success parameter 
    ## was also computed (and is not equal to the conditonal mean, as is the 
    ## case for the Bernoulli distribution). 
    if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
      newdata$p_prob     <- rowMeans(MC$prob_samples)
      newdata$RMSPE_prob <- sqrt(.rowVars(MC$prob_samples))
      newdata <- .concat_percentiles_to_df(X = MC$prob_samples, df = newdata, 
                                           name = "prob", percents = percents)
    }
  }
  
  ## Response variable, Z
  if ("response" %in% type){
    newdata$p_Z <- p_mu 
    newdata$RMSPE_Z <- sqrt(.rowVars(MC$Z_samples))
    newdata <- .concat_percentiles_to_df(X = MC$Z_samples, df = newdata, 
                                         name = "Z", percents = percents) 
  }

  
  # ---- CHECK: Add a column indicating the time point in space-time setting ----
  
  if (is(M@basis,"TensorP_Basis")) {
    newdata$t <- M@BAUs@data$t
  }

  # ---- Return output ----

  ## Remove MC samples that were not asked for
  if (!("link" %in% type)) MC$Y_samples <- NULL
  if (!("mean" %in% type)) {MC$mu_samples <- NULL; MC$prob_samples <- NULL}
  
  return(list(newdata = newdata, MC = MC, pred_time = pred_time))
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
.Y_var <- function(M, Q_L, obsidx){
  
  r  <- ncol(M@S0)
  m  <- length(obsidx)
  
  # ---- Sparse-inverse-subset of Q (acting as a proxy for the true covariance matrix) ----

  ## Sparse-inverse-subset of fixed AND random effects
  ## (a proxy for the covariance matrix)
  Sigma <- sparseinv::Takahashi_Davis(Q = M@Q_eta_xi,
                                      cholQp = Q_L$Qpermchol,
                                      P = Q_L$P)
  
  # Sigma <- chol2inv(chol(M@Q_eta_xi))
  # warning("Removed sparse-inverse-subset and using full inverse: 
  #         we cannot do this for large m + r.")
  
  Sigma_eta   <- Sigma[1:r, 1:r]
  Sigma_xi    <- Sigma[(r + 1):(r + m), (r + 1):(r + m)]
  
  # Covariances between xi and eta
  Cov_eta_xi  <- Sigma[1:r, ( r+ 1):(r + m)]
  
  # ----- Uncertainty: Posterior variance of Y at each BAU ------
  
  ## To extract the variances of eta|Z, we need diag(S0 %*% Sigma_eta %*% t(S0)).
  ## Also, to extract the covariances term, we need: diag(S %*% COV_{eta, xi}).
  ## This in very inefficient to do directly, it much better to use the identity:
  ##      diag(AB) = (A*B')1
  
  ## Only one common term for both observed and unobserved locations:
  vY <- as.vector( (M@S0 %*% Sigma_eta * M@S0) %*% rep(1, r) )
  
  ## UNOBSERVED locations: simply add the estimate of sigma2xi to this quantity:
  vY[-obsidx] <- vY[-obsidx] + M@sigma2fshat
  
  ## OBSERVED location: add both var(xi_O|Z) and cov(xi_O, eta | Z)
  covar       <- (M@S * Matrix::t(Cov_eta_xi)) %*% rep(1, r)        # Covariance terms
  vY[obsidx]  <- vY[obsidx] + Matrix::diag(Sigma_xi) + 2 * covar
  
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
#' @param M An object of class \code{SRE}.
#' @param X The design matrix of the covariates at the BAU level (often simply an Nx1 column vector of 1's).
#' @param type A character string (possibly vector) indicating the quantities which are the focus of inference.
#' Note: unlike in the predict() function, \emph{all} computed quantities are returned. That is, 
#' the latent \eqn{Y} process samples are always provided; 
#' If \code{"mean"} \emph{OR} \code{"response"} is in \code{type}, then the samples of \eqn{Y}, the conditonal mean \eqn{\mu}, and the probability parameter (if applicable) are provided.
#' If \code{"response"} is in \code{type}, the response variable \eqn{Z} samples, and the samples of all other quantities are provided. 
#' @param n_MC A postive integer indicating the number of MC samples at each location.
#' @param obs_fs Logical indicating whether the fine-scale variation is included in the latent Y process. 
#' @param k vector of size parameters parameters at each BAU (applicable only for binomial and negative-binomial data).
#' If \code{obs_fs = FALSE} (the default), then the fine-scale variation term \eqn{\xi} is included in the latent \eqn{Y} process. 
#' If \code{obs_fs = TRUE}, then the the fine-scale variation terms \eqn{\xi} are removed from the latent Y process; \emph{however}, they are re-introduced for computation of the conditonal mean \eqn{\mu} and response variable \eqn{Z}. 
#' @param Q_L A list containing the Cholesky factor of the permuted precision matrix (stored as \code{Q$Qpermchol}) and the associated permutationmatrix (stored as \code{Q_L$P}).
#' @param obsidx A vector containing the indices of observed locations.
#' @return A list containing Monte Carlo samples of various quantites of interest. 
#' The list elements are (N x n_MC) matrices, whereby the ith row of each matrix corresponds to
#' n_MC samples of the given quantity at the ith BAU. The available quantities are:
#' \describe{
#'   \item{Y_samples}{Samples of the latent, Gaussian scale Y process.}
#'   \item{mu_samples}{Samples of the conditional mean of the data.}
#'   \item{prob_samples}{Samples of the probability of success parameter (only for the relevant response distributions).}
#'   \item{Z_samples}{Samples of the response variable.}
#' }
.MC_sampler <- function(M, X, type = "mean", n_MC = 400, obs_fs = FALSE, k = NULL, 
                        Q_L, obsidx){
  
  MC <- list() # object we will return 
  N   <- nrow(M@S0)
  m   <- length(M@Z)
  r   <- ncol(M@S0) # Total number of basis functions
  
  MC$times <- list() # object to track timings of function components
  
  # ---- Generate samples from (eta, xi_O) ----
  
  ## Must generate samples jointly, as eta and xi_O are correlated.

  
  ## Construct the mean vector of (eta, xi_O).
  ## Also make an (r+m) x n_MC matrix whose columns are the mean vector of (eta, xi_O).
  mu_eta_xi_O         <- c(as.numeric(M@mu_eta), as.numeric(M@mu_xi_O))
  mu_eta_xi_O_Matrix  <- matrix(rep(mu_eta_xi_O, times = n_MC), ncol = n_MC)
  
  ## Generate (r+m) x n_MC samples from Gau(0, 1) distribution
  z <- matrix(rnorm((r + m) * n_MC), nrow = r + m, ncol = n_MC)
  
  ## Compute the Cholesky factor of  Q (the joint precision matrix of (eta', xi_O')').
  ## Then, to generate samples from (eta, xi_O), 
  ## use eta_xi_O = L^{-T} z + mu = U^{-1} z + mu, 
  ## where U upper cholesky factor of Q, so that Q = U'U.
  # U           <- Matrix::chol(M@Q_eta_xi)
  # eta_xi_O    <- as.matrix(Matrix::solve(U, z) + mu_eta_xi_O_Matrix)
  
  ## Method 2: sparseinv package, Cholesky of permuted Q, then backsolve
  U <- Matrix::t(Q_L$Qpermchol) # upper Cholesky factor of permuted joint precision matrix M@Q_eta_xi
  MC$times$backsolve <- system.time(
    x <- backsolve(U, z)          # x ~ Gau(0, A), where A is the permuted precision matrix i.e. A = P'QP
  )
  y <- Q_L$P %*% x              # y ~ Gau(0, Q^{-1})
  eta_xi_O  <- as.matrix(y + mu_eta_xi_O_Matrix) # add the mean to y
  
  
  ## Separate the eta and xi_O samples
  eta     <- eta_xi_O[1:r, ]
  xi_O    <- eta_xi_O[(r + 1):(r + m), ]
  
  ## We now have two matrices, eta and xi_O:
  ## row i of eta corresponds to n_MC MC samples of eta_i,
  ## row i of xi_O corresponds to n_MC MC samples of the fine-scale variation at the ith observed location.
  
  
  # ---- Generate samples from xi_U ----
  
  ## This is straightforward as each element of xi_U is independent of
  ## all other random effects in the model.
  ## All we have to do is make an (N-m) x n_MC matrix of draws from the
  ## Gaussian distribution with mean zero and variance equal to the fine-scale variance.
  xi_U <- matrix(rnorm((N - m) * n_MC, mean = 0, sd = sqrt(M@sigma2fshat)), 
                 nrow = N - m, ncol = n_MC)

  
  # ---- Construct samples from latent process Y ----
  
  ## We break the latent process down as: Y = Y_smooth + xi, 
  ## so that we may separate the fine-scale variation. 
  
  ## Split the covariate design matrix based on observed and unobserved samples
  X_O <- X[obsidx, ]
  X_U <- X[-obsidx, ]
  
  ## Observed Samples
  Y_smooth_O <- X_O %*% M@alphahat + M@S %*% eta
  
  ## Unobserved Samples
  S_U          <- M@S0[-obsidx, ] # Unobserved random effect 'design' matrix
  Y_smooth_U  <- X_U %*% M@alphahat + S_U %*% eta
  
  ## Combine samples
  Y_smooth_samples  <- rbind(Y_smooth_O, Y_smooth_U)
  xi_samples        <- rbind(xi_O, xi_U) 
  
  ## Use permutation matrix to get the correct (original) ordering
  ## FIXME: could just use row indexing to avoid matrix multiplication here
  unobsidx         <- (1:N)[-obsidx]       # Unobserved BAUs indices
  ids              <- c(obsidx, unobsidx)  # All indices (observed and unobserved)
  P                <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  Y_smooth_samples <- Matrix::t(P) %*% Y_smooth_samples 
  xi_samples       <- Matrix::t(P) %*% xi_samples
  
  ## Construct the samples from the latent process Y 
  Y_samples <- Y_smooth_samples + xi_samples
  
  Y_samples        <- as.matrix(Y_samples)
  Y_smooth_samples <- as.matrix(Y_smooth_samples)


  ## Outputted Y value depend on obs_fs
  if (obs_fs == FALSE) MC$Y_samples <- Y_samples
  if (obs_fs == TRUE)  MC$Y_samples <- Y_smooth_samples
  
  ## If Y is the ONLY quantity of interest, exit the function.
  if (!("mean" %in% type) && !("response" %in% type)) return(MC) 
  
  
  # ---- Apply inverse-link function to the samples to obtain conditional mean ----
  
  ## Past this point we must have xi in the Y process (the model breaks down otherwise).
  ## In the case of type == "all", we simply export Y_smooth_samples as the samples of Y.

  ## For families with a known constant parameter (binomial, negative-binomial),
  ## zeta() maps the Gaussian scale Y process to the probability parameter p.
  ## Then, we map p to the conditional mean mu via chi().
  ## For all other families, psi() maps Y directly to mu.
  ## The exception is negative-binomial with a log or square-root link, 
  ## in which case we map directly from Y to mu.
  
  ## Note that for all cases other than type == "link", we need to compute the conditonal mean samples.
  
  ## Create the relevant link functions.
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

  ## Output the mean samples. If probability parameter p was computed, also output.
  MC$mu_samples <- mu_samples
  if (exists("prob_samples")) MC$prob_samples <- prob_samples

  
  ## If the response is not a quanitity of interest, exit the function
  if (!("response" %in% type)) return(MC)
 
  
  # ---- Sample the response variable, Z ----
  
  if (M@response == "poisson") {
    Z_samples <- rpois(n = N * n_MC, lambda = c(t(mu_samples)))
  } else if (M@response == "gaussian") {
    sigma_e <- sqrt(M@Ve[1, 1]) # measurement error variance
    Z_samples <- rnorm(n = N * n_MC, mean = c(t(mu_samples)), sd = sigma_e)
  } else if (M@response == "bernoulli") {
    Z_samples <- rbinom(n = N * n_MC, size = 1, prob = c(t(mu_samples)))
  } else if (M@response == "gamma") {
    theta <- 1 / c(t(mu_samples)) # canonical parameter
    alpha <- 1/M@phi                 # shape parameter
    beta  <- theta * alpha           # rate parameter (1/scale)
    Z_samples <- rgamma(n = N * n_MC, shape = alpha, rate = beta)
    Z_samples <- statmod::rinvgauss(n = N * n_MC, mean = c(t(mu_samples)), dispersion = M@phi)
  } else if (M@response == "negative-binomial") {
    k_vec <- rep(k, each = n_MC)
    Z_samples <- rnbinom(n = N * n_MC, size = k_vec, prob = c(t(prob_samples)))
  } else if (M@response == "binomial") {
    k_vec <- rep(k, each = n_MC)
    theta <- log((c(t(mu_samples))/k_vec) / (1 - (c(t(mu_samples))/k_vec)))
    p <- 1 / (1 + exp(-theta))
    Z_samples <- rbinom(n = N * n_MC, size = k_vec, prob = p)
  }
  
  ## Convert from a long vector to an N * n_MC matrix
  Z_samples <- matrix(Z_samples, ncol = n_MC, byrow = TRUE)
  
  ## Add Z_samples to list object
  MC$Z_samples <- Z_samples
  

  
  return(MC)
}
