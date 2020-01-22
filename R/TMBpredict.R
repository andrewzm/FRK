#' Prediction stage of non-Gaussian FRK.
#'
#' @inheritParams .Y_pred
#' @inheritParams .MC_sampler
#' @return A list object containing:
#' \describe{
#'   \item{xy}{a dataframe with predictions and prediction uncertainty of the latent Y process, the conditional mean of the data (\eqn{\mu}), and, if applicable, the probability of success parameter \eqn{p}, at each prediction location. The dataframe also contains quantiles of each term.}
#'   \item{MC}{a list containing several matrices, which contain Monte Carlo samples of Y, \eqn{\mu} and, if applicable, \eqn{p}, at each prediction location.}
#' }
#' Note that for all link functions other than the log-link, the predictions and prediction uncertainty contained in \code{xy} are computed using the Monte Carlo samples contained in \code{MC}.
#' When a log-link function is used the expectation and variance of the conditional mean \eqn{\mu} may be computed exactly.
.FRKTMB_pred <- function(M, n_MC, seed, obs_fs = FALSE, type = "mean") {
  
  
  #### Extract the covariate design matrix, X 
  
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
  
  # ------ Latent process (Y) prediction and Uncertainty ------
  
  ## Compute posterior expectation and variance of Y at each prediction location.
  ## This does NOT depend on the response of Z.
  temp   <- .Y_pred(M = M, X = X)
  
  p_Y     <- temp$EYgivenZ    # Prediction of Y given the data
  MSPE_Y  <- temp$varYgivenZ  # Conditonal variance of Y given Z
  
  
  # ------ Conditional mean (mu) prediction and uncertainty ------
  
  ## Compute Monte Carlo samples of conditional mean at each location
  MC <- .MC_sampler(M = M, 
                    X = X,
                    type = type,
                    obs_fs = obs_fs,
                    n_MC = n_MC,
                    seed = seed)

  
  # ------ Create Prediction Dataframe ------
  
  ## This dataframe contains x and y, the spatial locations. 
  ## It also contains the predictions, the RMSPE, and 5 quantiles for the specified quantities of interest. 

  xy <- as.data.frame(coordinates(M@BAUs)) # xy coordinates of every BAU
  rownames(xy) <- NULL
  
  ## Latent Y-process
  if ("link" %in% type) {
    xy$p_Y     <- p_Y          # Prediction
    xy$RMSPE_Y <- sqrt(MSPE_Y) # Uncertainty
    
    ## Compute the quantiles
    temp <- t(apply(MC$Y_samples, 1, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    colnames(temp) <- c("Y_percentile_05", "Y_percentile_25", "Y_percentile_50", "Y_percentile_75", "Y_percentile_95")
    xy <- cbind(xy, temp)
  }

  ## Conditional mean of the data, mu
  if ("mean" %in% type) {
    ## If a log-link function is used, then expectations and variance
    ## of the conditional mean may be evaluated analytically.
    ## Otherwise, use the Monte Carlo simulations for prediction.
    if (M@link == "log" & M@response != "negative-binomial") {
      xy$p_mu         <- exp(p_Y + MSPE_Y / 2)
      xy$RMSPE_mu     <- sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
    } else if (M@link == "log" & M@response == "negative-binomial") {
      xy$p_mu         <- M@BAUs$k * exp(p_Y + MSPE_Y / 2)
      xy$RMSPE_mu     <- M@BAUs$k * sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
    } else {
      xy$p_mu         <- rowMeans(MC$mu_samples)
      xy$RMSPE_mu     <- sqrt(.rowVars(MC$mu_samples))
    }
    
    ## Compute the quantiles
    temp <- t(apply(MC$mu_samples, 1, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    colnames(temp) <- c("mu_percentile_05", "mu_percentile_25", "mu_percentile_50", "mu_percentile_75", "mu_percentile_95")
    xy <- cbind(xy, temp)
    
    ## For some response and link combinations, the probability of success parameter was also computed
    ## (and is not equal to the conditonal mean, as is the case for the Bernoulli distribution)
    if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
      xy$p_pi     <- rowMeans(MC$p_samples)
      xy$RMSPE_pi <- sqrt(.rowVars(MC$p_samples))
    } else if (M@response == "negative-binomial" & M@link == "log") {
      ## IDEA: people may be interested in the probability of success parameter for negative-binomial with log-link.
      ## In this case, we can estimate it using the mean and the known formula for the mean in terms of the probability of success. 
    }
    
  }
  
  ## Response variable, Z
  if ("response" %in% type){
    
    ## Compute predictions and variance
    xy$p_Z_analytic <- xy$p_mu
    xy$p_Z_empirical <- rowMeans(MC$Z_samples)
    xy$RMSPE_Z <- sqrt(.rowVars(MC$Z_samples))
    
    ## Compute the quantiles
    temp <- t(apply(MC$Z_samples, 1, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95)))
    colnames(temp) <- c("Z_percentile_05", "Z_percentile_25", "Z_percentile_50", "Z_percentile_75", "Z_percentile_95")
    xy <- cbind(xy, temp)
  }


  # ---- Return output ----

  out <- list(xy = xy, MC = MC)
  
  ## change "xy" name to newdata
  
  return(out)
}



#' Prediction and uncertainty of the latent Y process.
#'
#' Computes the posterior mean and variance of the latent process Y at every BAU.
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
#' @param M An object of class SRE.
#' @param X The design matrix of the covariates at the BAU level (often simply an Nx1 column vector of 1's) 
#' @return A list object containing:
#' \describe{
#'   \item{EYgivenZ}{A vector of the posterior expectation of Y at every BAU.}
#'   \item{varYgivenZ}{A vector of the posterior variance of Y at every BAU.}
#' }
.Y_pred <- function(M, X){
  
  Q  <- M@Q_eta_xi
  S0 <- M@S0
  S  <- M@S
  obsidx <- apply(M@Cmat, 1, function(x) which(x == 1))
  r  <- ncol(S0)
  m  <- length(obsidx)
  
  # ---- Sparse-inverse-subset of Q (acting as a proxy for the true covariance matrix) ----
  
  ## Permuted Cholesky factor
  Q_L <- sparseinv:::cholPermute(Q = Q)
  
  ## Sparse-inverse-subset of fixed AND random effects
  ## (a proxy for the covariance matrix)
  Sigma <- sparseinv::Takahashi_Davis(Q = Q,
                                      cholQp = Q_L$Qpermchol,
                                      P = Q_L$P)

  Sigma_eta   <- Sigma[1:r, 1:r]
  Sigma_xi    <- Sigma[(r + 1):(r + m), (r + 1):(r + m)]
  
  # Covariances between xi and eta
  Cov_eta_xi  <- Sigma[1:r, ( r+ 1):(r + m)]
  
  # ----- Prediction: Posterior mean of Y at each BAU ------
  
  
  ## large- and medium-scale variation terms
  mY <-  as.vector(X %*% M@alphahat + S0 %*% M@mu_eta)
  
  ## Add posterior estimate of xi_O at observed BAUs
  mY[obsidx]   <-  mY[obsidx] + as.vector(M@mu_xi_O)
  
  
  # ----- Uncertainty: Posterior variance of Y at each BAU ------
  
  ## To extract the variances of eta|Z, we need diag(S0 %*% Sigma_eta %*% t(S0)).
  ## Also, to extract the covariances term, we need: diag(S %*% COV_{eta, xi}).
  ## This in very inefficient to do directly, it much better to use the identity:
  ##      diag(AB) = (A*B')1
  
  ## Only one common term for both observed and unobserved locations:
  vY <- as.vector( (S0 %*% Sigma_eta * S0) %*% rep(1, r) )
  
  ## UNOBSERVED locations: simply add the estimate of sigma2xi to this quantity:
  vY[-obsidx] <- vY[-obsidx] + M@sigma2fshat
  
  ## OBSERVED location: add both var(xi_O|Z) and cov(xi_O, eta | Z)
  covar       <- (S * Matrix::t(Cov_eta_xi)) %*% rep(1, r)        # Covariance terms
  vY[obsidx]  <- vY[obsidx] + Matrix::diag(Sigma_xi) + 2 * covar
  
  # ---- Output ----
  
  ## Return the posterior mean and variance of Y as a list object
  out <- list("EYgivenZ" = mY, "varYgivenZ" = vY)
  
  return(out)
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
#' If \code{"link"} is in \code{type}, the latent \eqn{Y} process samples are provided. 
#' If \code{"mean"} is in \code{type}, the conditonal mean \eqn{\mu} samples are provided (and the probability parameter if applicable).
#' If \code{"response"} is in \code{type}, the response variable \eqn{Z} samples are provided. 
#' For example, if \code{type = c("link", "response")}, then MC samples of the latent \eqn{Y} process and the response variable are provided.
#' @param n_MC A postive integer indicating the number of MC samples at each location.
#' @param obs_fs Logical indicating whether the fine-scale variation is included in the latent Y process. 
#' If \code{obs_fs == FALSE} (the default), then the fine-scale variation term \eqn{\xi} is included in the latent \eqn{Y} process. 
#' If \code{obs_fs == TRUE}, then the the fine-scale variation terms \eqn{\xi} are removed from the latent Y process; \emph{however}, they are re-introduced for computation of the conditonal mean, \eqn{mu}. 
#' @param seed A seed for reproducibility.
#' @return A list containing Monte Carlo samples of various quantites of interest. 
#' The list elements are (N x n_MC) matrices, whereby the ith row of each matrix corresponds to
#' n_MC samples of the given quantity at the ith BAU. The available quantities are:
#' \describe{
#'   \item{Y_samples}{Samples of the latent, Gaussian scale Y process.}
#'   \item{mu_samples}{Samples of the conditional mean of the data.}
#'   \item{Z_samples}{Samples of the response variable.}
#'   \item{p_samples}{Samples of the probability of success parameter at the ith BAU.}
#' }
.MC_sampler <- function(M, X, type = "mean", n_MC = 1600, obs_fs = FALSE, seed = NULL){
  
  MC <- list() # object we will return 
  
  k_BAU <- M@BAUs$k
  Q   <- M@Q_eta_xi
  S0  <- M@S0
  S   <- M@S
  N   <- nrow(S0)
  obsidx <- apply(M@Cmat, 1, function(x) which(x == 1))
  m   <- length(obsidx)
  r   <- ncol(S0)
  
  
  # ---- Generate samples from (eta, xi_O) ----
  
  ## Must generate samples jointly, as eta and xi_O are correlated.
  
  ## Construct the mean vector of (eta, xi_O).
  ## Also make an (r+m) x n_MC matrix whose columns are the mean vector of (eta, xi_O).
  mu_eta_xi_O         <- c(as.numeric(M@mu_eta), as.numeric(M@mu_xi_O))
  mu_eta_xi_O_Matrix  <- matrix(rep(mu_eta_xi_O, times = n_MC), ncol = n_MC)
  
  ## Generate (r+m) x n_MC samples from Gau(0, 1) distribution
  set.seed(seed)
  z <- matrix(rnorm((r + m) * n_MC), nrow = r + m, ncol = n_MC)
  
  ## Compute the Cholesky factor of  Q (the joint precision matrix of (eta', xi_O')').
  ## Then, to generate samples from (eta, xi_O), 
  ## use eta_xi_O = L^{-T} z + mu = U^{-1} z + mu, 
  ## where U upper cholesky factor of Q, so that Q = U'U.
  U           <- Matrix::chol(Q)
  eta_xi_O    <- as.matrix(Matrix::solve(U, z) + mu_eta_xi_O_Matrix)
  
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
  Y_smooth_O <- X_O %*% M@alphahat + S %*% eta
  
  ## Unobserved Samples
  S_U          <- S0[-obsidx, ] # Unobserved random effect 'design' matrix
  Y_smooth_U  <- X_U %*% M@alphahat + S_U %*% eta
  
  ## Combine samples
  Y_smooth_samples  <- rbind(Y_smooth_O, Y_smooth_U)
  xi_samples        <- rbind(xi_O, xi_U) 
  
  ## Use permutation matrix to get the correct (original) ordering
  unobsidx         <- (1:N)[-obsidx]       # Unobserved BAUs indices
  ids              <- c(obsidx, unobsidx)  # All indices (observed and unobserved)
  P                <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  Y_smooth_samples <- Matrix::t(P) %*% Y_smooth_samples 
  xi_samples       <- Matrix::t(P) %*% xi_samples
  
  ## Construct the samples from the latent process Y 
  Y_samples <- Y_smooth_samples + xi_samples
  
  ## Convert to matrix objects
  Y_samples        <- as.matrix(Y_samples)
  Y_smooth_samples <- as.matrix(Y_smooth_samples)


  ## Outputted Y value depend on obs_fs
  if ("link" %in% type) {
    if (obs_fs == FALSE) MC$Y_samples <- Y_samples
    if (obs_fs == TRUE)  MC$Y_samples <- Y_smooth_samples
  }
  
  ## If Y is the ONLY quantity of interest, exit the function.
  if (length(type) == 1) return(MC) 
  
  
  # ---- Apply inverse-link function to the samples to obtain conditional mean ----
  
  ## Past this point we must have xi in the Y process (the model breaks down otherwise).
  ## In the case of type == "all", we simply export Y_smooth_samples as the samples of Y.

  ## For families with a known constant parameter (binomial, negative-binomial),
  ## psi() maps the Gaussian scale Y process to the probability parameter p.
  ## Then, we map p to the conditional mean mu via psi_mu().
  ## For all other families, psi() maps Y directly to mu.
  
  ## Note that for all cases other than type == "link", we need to compute the conditonal mean samples.
  
  psi     <- .inv_link_fn(M@link)   # link function (either to p or directly to mu)
  
  #browser()
  ## FIX: I don't like how psi links to the mean or the prob parameter
  ## It would be better to always link to the mean, and for distributions that require it we can construct the probability parameter from the mean.
  if (M@response == "binomial" & M@link %in% c("logit", "probit", "cloglog")) {
    p_samples <- psi(Y_samples)
    mu_samples <- k_BAU * p_samples
  } else if (M@response == "negative-binomial" & M@link %in% c("logit", "probit", "cloglog")) {
    p_samples <- psi(Y_samples)
    mu_samples <- k_BAU * (1 / p_samples - 1)
  } else if (M@response == "negative-binomial" & M@link %in% c("log", "square-root")) {
    mu_samples <- k_BAU * psi(Y_samples)
  } else {
    mu_samples <- psi(Y_samples)
  }

  ## Output the mean samples if requested.
  ## If probability parameter p was computed, also output.
  if ("mean" %in% type) {
    MC$mu_samples <- mu_samples
    if (exists("p_samples")) MC$p_samples <- p_samples
  }
  
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
  } else if (M@response == "inverse-gaussian") {
    Z_samples <- statmod::rinvgauss(n = N * n_MC, mean = c(t(mu_samples)), dispersion = M@phi)
  } else if (M@response == "negative-binomial") {
    k_BAU_vec <- rep(k_BAU, each = n_MC)
    theta <- log(c(t(mu_samples)) / (k_BAU_vec + c(t(mu_samples))))
    p <- 1 - exp(theta)
    Z_samples <- rnbinom(n = N * n_MC, size = k_BAU_vec, prob = p)
  } else if (M@response == "binomial") {
    k_BAU_vec <- rep(k_BAU, each = n_MC)
    theta <- log((c(t(mu_samples))/k_BAU_vec) / (1 - (c(t(mu_samples))/k_BAU_vec)))
    p <- 1 / (1 + exp(-theta))
    Z_samples <- rbinom(n = N * n_MC, size = k_BAU_vec, prob = p)
  }
  
  ## Convert from a long vector to an N * n_MC matrix
  Z_samples <- matrix(Z_samples, ncol = n_MC, byrow = TRUE)
  
  ## Add Z_samples to list object
  MC$Z_samples <- Z_samples
  
  return(MC)
}