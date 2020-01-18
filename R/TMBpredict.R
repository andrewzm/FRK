#' Prediction stage of non-Gaussian FRK.
#'
#' 
#' @inheritParams .Y_pred
#' @inheritParams .MC_sampler
#' @return A list object containing:
#' \describe{
#'   \item{xy}{a dataframe with predictions and prediction uncertainty of the latent Y process, the conditional mean of the data (\eqn{\mu}), and, if applicable, the probability of success parameter \eqn{p}, at each prediction location.}
#'   \item{MC}{a list containing several matrices, which contain Monte Carlo samples of Y, \eqn{\mu} and, if applicable, \eqn{p}, at each prediction location.}
#' }
#' Note that for all link functions other than the log-link, the predictions and prediction uncertainty contained in \code{xy} are computed using the Monte Carlo samples contained in \code{MC}.
#' When a log-link function is used the expectation and variance of the conditional mean \eqn{\mu} may be computed exactly.
.FRKTMB_pred <- function(M, n_MC, seed, obs_fs = FALSE) {
  
  # ------ Latent process (Y) Prediction and Uncertainty ------
  
  ## Compute posterior expectation and variance of Y at each prediction location.
  ## This does NOT depend on the response of Z.
  temp   <- .Y_pred(M = M)
  
  p_Y     <- temp$EYgivenZ    # Prediction of Y given the data
  MSPE_Y  <- temp$varYgivenZ  # Conditonal variance of Y given Z
  
  
  # ------ Conditional mean (mu) prediction and uncertainty ------
  
  ## Compute Monte Carlo samples of conditional mean at each location
  MC <- .MC_sampler(M = M, 
                    n_MC = n_MC,
                    seed = seed,
                    k_BAU = M@BAUs$k)
  
  
  # ------ Create Prediction Dataframe ------
  
  ## Dataframe containing spatial location, predictions, and uncertainty.
  xy <- as.data.frame(coordinates(M@BAUs)) # xy coordinates of every BAU
  rownames(xy) <- NULL
  
  ## Latent Y-process
  xy$p_Y     <- p_Y          # Prediction
  xy$RMSPE_Y <- sqrt(MSPE_Y) # Uncertainty
  
  ## Conditional mean of the data, mu
  ## Compute the prediction (posterior mean) and uncertainty (posterior variance)
  
  ## If a log-link function is used, then expectations and variance
  ## of the conditional mean may be evaluated analytically.
  ## Otherwise, use the Monte Carlo simulations for prediction.
  if (M@link == "log" & M@response != "negative-binomial") {
    xy$p_mu         <- exp(p_Y + MSPE_Y / 2)
    xy$RMSPE_mu     <- sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  } else if (M@link == "log" & M@response == "negative-binomial") {
    xy$p_mu         <- k_BAU * exp(p_Y + MSPE_Y / 2)
    xy$RMSPE_mu     <- k_BAU * sqrt((exp(MSPE_Y) - 1) * exp(2 * p_Y + MSPE_Y))
  } else {
    xy$p_mu         <- rowMeans(MC$mu_samples)
    xy$RMSPE_mu     <- sqrt(.rowVars(MC$mu_samples))
  }
  
  ## For some response and link combinations, the probability of success parameter was also computed:
  if (M@response %in% c("binomial", "negative-binomial") & M@link %in% c("logit", "probit", "cloglog")) {
    xy$p_pi     <- rowMeans(MC$p_samples)
    xy$RMSPE_pi <- sqrt(.rowVars(MC$p_samples))
  }
  
  out <- list(xy = xy,
              MC = MC)
  
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
#' @param Q The joint precision matrix of the random effects \eqn{(\eta', \xi_O')'}.
#' @param estimates A list object containing the estimates of parameters and fixed-effects, and predictions of the random effects.
#' @param M An object of class SRE.
#' @return A list object containing:
#' \describe{
#'   \item{EYgivenZ}{A vector of the posterior expectation of Y at every BAU.}
#'   \item{varYgivenZ}{A vector of the posterior variance of Y at every BAU.}
#' }
.Y_pred <- function(M){
  
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
  mY           <-  as.vector(as.numeric(M@alphahat) + S0 %*% as.numeric(M@mu_eta))
  
  ## Add posterior estimate of xi_O at observed BAUs
  mY[obsidx]   <-  mY[obsidx] + as.numeric(M@mu_xi_O)
  
  
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
#' Computes a Monte Carlo sample of Y, the conditional mean of the data
#' \eqn{\mu = \psi(Y)} (which is a function of Y), and, for response-link
#' combinations to which it is applicable, the probability of success parameter
#' p. It does so for every BAU location.
#'
#'
#' @inheritParams .inv_link_fn
#' @param M An object of class \code{SRE}.
#' @param n_MC A postive integer indicating the number of MC samples at each location.
#' @param k_BAU An integer vector of length N (the number of BAUs). Relevant only to the binomial and negative-binomial cases.
#' @param response A string indicating the assumed distribution of the response variable.
#' For negative-binomial data, the ith element of \code{k_BAU} indicates the number of failures until the experiment is stopped at the ith BAU.
#' For binomial data, the ith element of \code{k_BAU} indicates the number of trials for the ith BAU.
#' @return A list containing:
#' 1. The Monte Carlo expectation of \eqn{\mu = \psi(Y)} at each BAU location.\cr
#' 2. The Monte Carlo variance of \eqn{\mu = \psi(Y)} at each BAU location.\cr
#' 3. The full Monte Carlo samples.
#' \describe{
#'   \item{Y_samples}{An (N x n_MC) matrix, whereby the ith row of the matrix correpsponds to
#'   n_MC samples of the latent Y process at the ith BAU.}
#'   \item{mu_samples}{An (N x n_MC) matrix, whereby the ith row of the matrix correpsponds to
#'   n_MC samples of the conditional mean of the data at the ith BAU.}
#'   \item{p_samples}{An (N x n_MC) matrix, whereby the ith row of the matrix correpsponds to
#'   n_MC samples of the probability of success parameter p at the ith BAU.}
#' }

.MC_sampler <- function(M, n_MC = 1600, seed = NULL, k_BAU, obs_fs = FALSE){
  
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
  
  if (obs_fs == FALSE) {
    xi_U <- matrix(rnorm((N - m) * n_MC, mean = 0, sd = sqrt(M@sigma2fshat)),
                         nrow = N - m, ncol = n_MC)
  } else if (obs_fs == TRUE) {
    xi_U <- 0
  }

  # ---- Construct samples from latent process Y ----
  
  ## Observed Samples
  Y_O <- as.matrix(as.numeric(M@alphahat) + S %*% eta + xi_O)
  
  ## Unobserved Samples
  SU      <- S0[-obsidx, ] # Unobserved random effect 'design' matrix
  Y_U     <- as.matrix(as.numeric(M@alphahat) + SU %*% eta + xi_U)
  
  ## Combine samples
  Y_samples <- rbind(Y_O, Y_U)
  
  ## Use permutation matrix to get the correct (original) ordering
  unobsidx        <- (1:N)[-obsidx]       # Unobserved BAUs indices
  ids             <- c(obsidx, unobsidx)  # All indices (observed and unobserved)
  P               <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  Y_samples       <- as.matrix(Matrix::t(P) %*% Y_samples)
  
  
  ## Y_samples is an N x n_MC matrix, whereby the ith row of the matrix
  ## correpsponds to n_MC samples of the latent Y process at the ith BAU.
  
  
  # ---- Apply inverse-link function to the samples ----
  
  ## For families with a known constant parameter (binomial, negative-binomial),
  ## psi() maps the Gaussian scale Y process to the probability parameter p.
  ## Then, we map p to the conditional mean mu via psi_mu().
  ## For all other families, psi() maps Y directly to mu.
  
  psi     <- .inv_link_fn(M@link)   # link function (either to p or directly to mu)
  
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
  
  ## Return a list containing the Monte Carlo samples of mu.
  out <- list(mu_samples = as.matrix(mu_samples),
              Y_samples = as.matrix(Y_samples))
  
  ## If probability parameter p was computed, also output:
  if (exists("p_samples")) out$p_samples = as.matrix(p_samples)
  
  
  return(out)
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
