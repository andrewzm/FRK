#' Deprecated: Please use \code{\link{predict}}
#'
#' @param ... (Deprecated)
#' @export
SRE.predict <- function(...) {
  stop("SRE.predict() is deprecated. Please use predict().")
}

#' @rdname SRE
#' @export
setMethod("predict", signature="SRE", function(object, newdata = NULL, obs_fs = FALSE,
                                               pred_time = NULL, covariances = FALSE, 
                                               nsim = 400, type = "mean", k = NULL, 
                                               percentiles = c(5, 95), 
                                               kriging = "simple", new_BAU_data = NULL) {
  
  ## Store the original newdata for comparison at the end of the function
  # newdata_original <- newdata
  
  ## Note that the following code must be done before .check_args3() because it
  ## alters the input, rather than simpling throwing a warning or error. 
  
  ## Need to add prediction and uncertainty at each location, and so newdata 
  ## must be a Spatial*DataFrame (not just a Spatial* object).
  newdata <- .Coerce_SpatialDataFrame(newdata)
  
  ## The user can either provide k in object@BAUs$k_BAU, or in the predict call.
  ## This is so that the user can change k without having to call SRE() and SRE.fit() again.
  ## The k supplied in predict() will take precedence over the k stored in object@BAUs$k.
  
  k_BAU <- object@BAUs$k_BAU # size parameter at the BAU level (could be NULL)
  N <- length(object@BAUs)   # number of BAUs
  
  if (object@response %in% c("binomial", "negative-binomial")) {
    if (is.null(k)) {
      if(is.null(k_BAU)) {
        k <- rep(1, N)
        cat("Assuming that the size parameter, k, is equal to 1 for all prediction locations: If you want to change this, use the argument 'k' in predict().\n")
      } else {
        if (any(is.na(k_BAU))) {
          
          ## If there's only a single non-NA value, it's safe to assume 
          ## That this value can be used for all BAUs (this typically happens
          ## in the case of Bernoulli data, where the non-NA value is 1).
          k_BAU_NA_removed <- k_BAU[!is.na(k_BAU)]
          if (length(unique(k_BAU_NA_removed)) == 1) {
            k <- rep(unique(k_BAU_NA_removed), N)
          } else {
            k <- rep(1, N)
            cat("Assuming that the size parameter, k, is equal to 1 for all prediction locations: If you want to change this, use the argument 'k' in predict.\n")
          }
        } else {
          k <- k_BAU
        }
      }
    }
  }
  
  ## Check the arguments are OK
  .check_args3(obs_fs = obs_fs, newdata = newdata, 
               pred_time = pred_time, covariances = covariances, 
               object = object, type = type, kriging = kriging,
               k = k, percentiles = percentiles)
  
  ## If newdata is SpatialPoints* object, then we predict over irregularly 
  ## spaced points. Do this by first predicting over the BAUs, and then 
  ## associating each prediction location with a BAU. 
  ## NB: Use a long name for the prediction points, exists() is used later to
  ## see if it is defined; using an obvious name such as "pred_points" (which 
  ## the user may have already defined in their global scope) may cause issues, 
  ## as global variables are visible from inside functions.
  if(is(newdata, "SpatialPoints") || is(newdata, "STI")) {
    pred_points_from_user <- newdata # save the prediction locations for use later
    newdata <- NULL # setting newdata to NULL means we predict over the BAUs
  }
  
  ## If the user does not specify time points to predict at when in space-time
  ## Then predict at every time point
  if(is.null(pred_time) & is(object@BAUs,"ST"))
    pred_time <- 1:length(object@BAUs@time)
  
  ## Construct the CP matrix (Polygon prediction matrix)
  CP <- .make_CP(newdata = newdata, object = object)
  
  ## We start by assuming that we will predict at BAUs
  predict_BAUs <- TRUE
  
  ## If even one polygon encompasses more than one BAU, then we need to
  ## predict over linear combinations of BAUs, and hence need to
  ## compute the full covariance matrix. Note this by setting
  ## predict_BAUs to be FALSE
  if (!is.null(newdata)) {
    tmp <- .as(CP, "dgTMatrix")
    if (!all(table(tmp@i) == 1))
      predict_BAUs <- FALSE
  }
  
  ## Replace BAU data if provided 
  if (!is.null(new_BAU_data)) { 
    if (!is.data.frame(new_BAU_data)) {
      stop("new_BAU_data must be a data.frame.")
    }
    if (nrow(new_BAU_data) != nrow(object@BAUs@data)) {
      stop("Number of rows in new_BAU_data does not match object@BAUs@data.")
    }
    if (!identical(names(new_BAU_data), names(object@BAUs@data))) {
      stop("Column names and their order in new_BAU_data must match object@BAUs@data.")
    }
    for (colname in names(new_BAU_data)) {
      if (class(new_BAU_data[[colname]])[1] != class(object@BAUs@data[[colname]])[1]) {
        stop(sprintf("Column '%s' in new_BAU_data does not have the same class as in object@BAUs@data.", colname))
      }
    }
    object@BAUs@data <- new_BAU_data
  }
  
  ## Call internal prediction functions depending on which method is used
  if (object@method == "EM") {
    pred_locs <- .SRE.predict_EM(object = object,              # Fitted SRE model
                                 obs_fs = obs_fs,             # Case 1 or Case 2?
                                 newdata = newdata,           # Prediction polygons
                                 pred_time = pred_time,       # Prediction time points
                                 covariances = covariances,   # Compute covariances?
                                 CP = CP,                     # Polygon prediction matrix
                                 predict_BAUs = predict_BAUs) # Are we predicting at BAUs?                   
    
  } else if (object@method == "TMB") {
    pred_locs <- .SRE.predict_TMB(object = object,               # Fitted SRE model
                                  newdata = newdata,           # Prediction polygons
                                  CP = CP,                     # Polygon prediction matrix
                                  predict_BAUs = predict_BAUs, # Are we predicting at BAUs?
                                  pred_time = pred_time,       # Prediction time points
                                  nsim = nsim,                 # Desired number of MC simulations
                                  obs_fs = obs_fs,             # Case 1 or Case 2?
                                  type = type,                 # Whether we are interested in the "link" (Y-scale), "mean", "response"
                                  k = k,                       # Size parameter
                                  percentiles = percentiles,   # Desired percentiles of MC samples
                                  kriging = kriging)          
  } 
  
  ## If the user provided irregular points (SpatialPoints* or STI*) to predict over, 
  ## associate each prediction location with a BAU, and use the prediction
  ## of this BAU as the prediction at the corresponding point. 
  if (exists("pred_points_from_user")) {
    
    ## The location of the BAU level predictions depends on 
    ## whether method = "TMB" or method = "EM".
    if(object@method == "TMB") {
      
      ## If method = "TMB", we have a list of MC samples which must
      ## also be altered. To do this, we select BAUs based on identifiers.
      ## We may not necessarily have identifiers, particularly if 
      ## the BAUs were not created with auto_BAUs(). Create some:
      BAU_UID <- .UIDs(object@BAUs)
      ## NB: BAU_name is created when calling map_data_to_BAUs(), but this isn't
      ## exactly what we want when in a spatio-temporal setting (it is not unique
      ## with respect to time).
      
      ## Add BAU_UID to prediction data so it gets included after map_data_to_BAUs()
      pred_locs$newdata@data$BAU_UID <- BAU_UID
      
      ## Map the data to the BAUs
      pred_locs$newdata <- map_data_to_BAUs(pred_points_from_user, pred_locs$newdata, 
                                            average_in_BAU = FALSE, silently = TRUE)
      
      ## For some reason, pred_locs$newdata@data$BAU_UID becomes a factor, 
      ## And this has a weird effect on the reordering below. Change to character 
      ## for the desired behaviour.
      pred_locs$newdata@data$BAU_UID <- as.character(pred_locs$newdata@data$BAU_UID)
      
      ## Add the BAU UIDs to the rows of the MC matrices, then subset based 
      ## on the BAU UIDs in the predictions.
      pred_locs$MC <- lapply(pred_locs$MC, function(X) { 
        rownames(X) <- BAU_UID
        X <- X[pred_locs$newdata$BAU_UID, ]
        rownames(X) <- NULL
        return(X)
      })
      
      ## Clean up returned prediction data
      pred_locs$newdata@data[, c("BAU_name", "BAU_UID")] <- NULL
      
    } else if (object@method == "EM") {
      ## If we don't have to reorder MC matrices, then our task is simple:
      pred_locs <- map_data_to_BAUs(pred_points_from_user, pred_locs, average_in_BAU = FALSE)
      ## Clean up returned prediction data
      pred_locs@data[, c("BAU_name")] <- NULL
    }
  }
  
  ## Return predictions
  return(pred_locs)
})

## If newdata is a SpatialPoints, SpatialPolygons, or SPatialPixels object, 
## this function coerces newdata to is Spatial*DataFrame variant (with an empty 
## data slot). If newdata is NULL, NULL is returned.
.Coerce_SpatialDataFrame <- function (newdata) {
  
  if(!is.null(newdata)) {
    ## Add dummy data (with a name that is extremely unlikely to have been used already)
    newdata$dummy_blah123 <- rep(1, length(newdata)) 
    ## NULL it so this dummy data is not returned to the user
    newdata$dummy_blah123 <- NULL 
  }
  
  return(newdata)
}


## Construct CP (the incidence matrix for prediction)
.make_CP <- function (newdata = NULL, object) {
  
  ## If the user has not specified polygons over which to predict, then CP is
  ## just the diagonal matrix and we predict over all the BAUs
  if(is.null(newdata)) {
    CP <- Diagonal(length(object@BAUs))
  } else {
    ## The user has specified arbitrary polygons
    ## Based on these polygons construct the C matrix
    newdata2 <- map_data_to_BAUs(newdata, object@BAUs,
                                 average_in_BAU = FALSE,
                                 sum_variables = NULL,
                                 est_error = FALSE)
    C_idx <- BuildC(newdata2, object@BAUs)
    
    CP <- sparseMatrix(i = C_idx$i_idx,
                       j = C_idx$j_idx,
                       x = C_idx$x_idx,
                       dims = c(length(newdata),
                                length(object@BAUs)))
    
    ## As in SRE(), make sure the polygons are averages (not sums) if requested
    if (object@normalise_wts)
      CP <- CP / rowSums(CP)
  }
  return(CP)
}


# ---- Prediction functions ----

.SRE.predict_EM <- function(object, obs_fs = FALSE, newdata = NULL, 
                            pred_time = NULL, covariances = FALSE,
                            CP, predict_BAUs) {
  
  ## Get BAUs from the SRE model
  BAUs <- object@BAUs
  
  ## If we have to compute too many covariances then stop and give error
  if(covariances & nrow(CP) > 4000)
    stop("Cannot compute covariances for so many prediction locations. Please reduce to less than 4000")
  
  ## Get the CZ matrix
  CZ <- object@Cmat
  
  ## If the user has specified which polygons he wants we can remove the ones we don't need
  ## We only need those BAUs that are influenced by observations and prediction locations
  ## For ST use all BAUs as it gets complicated
  if(is(newdata,"Spatial")) {
    
    ## The needed BAUs are the nonzero column indices of CZ and CP
    needed_BAUs <- union(.as(CP,"dgTMatrix")@j+1,
                         .as(CZ,"dgTMatrix")@j+1)
    
    ## Filter the BAUs and the matrices
    BAUs <- BAUs[needed_BAUs,]
    CP <- CP[,needed_BAUs]
    CZ <- CZ[,needed_BAUs]
    object@S0 <- object@S0[needed_BAUs,]
  }
  
  ## Extract covariates from BAUs
  X <- .extract_BAU_X_matrix(object@f, BAUs)
  
  ## Set variables to make code more concise
  S0 <- .as(object@S0,"dgCMatrix")   
  alpha <- object@alphahat
  K <- object@Khat
  sigma2fs <- object@sigma2fshat
  mu_eta <- object@mu_eta
  S_eta <- object@S_eta
  
  if(object@fs_model == "ind") {
    D <- sigma2fs*object@Vfs + object@Ve   # compute total variance (data)
    if(isDiagonal(D)) {
      D <- Diagonal(x=D@x)       # cast to diagonal
      Dchol <- sqrt(D)           # find the inverse of the variance-covariance matrix
      Dinv <- solve(D)           # if Diagonal the Cholesky and inverse are just sqrt and reciprocal
    } else {
      Dchol <- chol(D)           # otherwise they need to be computed in full
      Dinv <- chol2inv(Dchol)
    }
    if(sigma2fs > 0) {
         sig2_Vfs_pred <- Diagonal(x=sigma2fs*BAUs$fs)   # fine-scale variation including estimated factor
         Q <- solve(sig2_Vfs_pred)                       # precision of fine-scale variation
    }
  } else  {
    stop("Prediction for other models not yet implemented")
  }
  
  ## The prediction equations
  ## If !obs_fs, then we are in Case 2 (default). We have to cater for when the
  ## fine-scale variance is zero or non-zero
  
  ## Case 2 (fs variation in process)
  if(!obs_fs) {
    if(sigma2fs > 0) {   # fine-scale variance not zero
      
      ## The below equations implement Section 2.3
      LAMBDAinv <- bdiag(object@Khat_inv,Q)                # block diagonal precision matrix
      PI <- cbind(S0, .symDiagonal(n=length(BAUs)))    # PI = [S I]
      tC_Ve_C <- t(CZ) %*% solve(object@Ve) %*% CZ +       # summary matrix
        0*.symDiagonal(ncol(CZ))              # Ensure zeros on diagonal
      Qx <- t(PI) %*% tC_Ve_C %*% PI + LAMBDAinv       # conditional precision matrix
      chol_Qx <- cholPermute(.as(Qx,"dgCMatrix"))       # permute and do Cholesky
      ybar <- t(PI) %*%t(CZ) %*% solve(object@Ve) %*%      # Qx = ybar, see vignette
        (object@Z - CZ %*% X %*% alpha)
      x_mean <- cholsolve(Qx,ybar,perm=TRUE,           # solve Qx = ybar using permutations
                          cholQp = chol_Qx$Qpermchol, P = chol_Qx$P)
      
      if(predict_BAUs) {
        ## If we are predicting at BAUs then we only need a few covariance elements
        ## to compute the marginal variances. Find these elements
        Cov <- Takahashi_Davis(Qx,cholQp = chol_Qx$Qpermchol,P = chol_Qx$P)
        
        ## Compute the variance and std over the BAUs in batches
        BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = TRUE)
      } else {
        ## Do not compute covariance now, do it later
      }
      
      
    } else {
      ## If sigma2fs = 0 then the prediction is much simpler and all our
      ## predictions / uncertainty come from the random effects
      PI <- S0
      x_mean <- object@mu_eta  # conditional mean of eta
      Cov <- object@S_eta      # conditional covariance of eta
      
      ## Compute variances, this time indicating there is no fs variation in process
      BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = FALSE)
    }
    
    ## For both cases of sigma2fs,
    ## the conditional mean is simply given by fitted random effects + fitted fixed effects
    BAUs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)
    
  }
  
  ## Case 1 (fs variation in measurement equation)
  if(obs_fs) {
    ## All predictions and prediction uncertainties comes from our
    ## prediction of and uncertainty over eta
    x_mean <- object@mu_eta   # conditional mean
    Cov <- object@S_eta       # conditional covariance
    
    ## Compute variances, this time indicating there is no fs variation in process
    BAUs[["mu"]] <- as.numeric(X %*% alpha) + as.numeric(S0 %*% x_mean)
    BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = FALSE)
  }
  
  
  ## Now, if the user hasn't specified prediction polygons, and the user does not want covariances,
  ## our job is done and we just return the BAUs, possibly at selected time points
  if(predict_BAUs & !covariances) {
    BAUs[["sd"]] <- sqrt(BAUs[["var"]])  # compute the standard error
    if(is(newdata, "Spatial")) {
      # User had specified a specific set of BAUs. Return only these (spatial only for now)
      BAUs <- BAUs[row.names(newdata),]
    }
    if(!is.null(pred_time))
      BAUs <- BAUs[,pred_time]  # return only specified time points
    
    return(BAUs)
    
  } else {
    ## Otherwise we need to find the mean and variance of linear combinations of these BAUs
    
    ## The linear combination of the mean is easy
    newdata[["mu"]] <- as.numeric(CP %*% BAUs[["mu"]])
    
    ## If we have fs variation in the process layer we need to consider the
    ## fine-scale variation (PI = [S I]) when predicting over the polygons,
    ## otherwise we just need the variance over eta
    
    ## If there is no fine-scale variation then simply find linear combination
    if(obs_fs | sigma2fs == 0 )    {
      CPM <- CP %*% S0
      newdata[["var"]] <- diag2(CPM %*% Cov, t(CPM)) ## All Cov available
      if(covariances) {
        L <- t(chol(Cov))
        Covariances <- tcrossprod(CPM %*% L)
      }
    } else {
      ## Otherwise find full covariance matrix (including over fs-variation).
      ## This is a last-case resort and
      ## may crash the computer if there are several prediction polygons.
      ## However this shouldn't be the case if these polygons span multiple BAUs
      CPM <- CP %*% PI
      newdata[["var"]] <- diag2(CPM, cholsolve(Q=Qx,y=t(CPM),
                                               perm = TRUE,
                                               cholQp = chol_Qx$Qpermchol,
                                               P = chol_Qx$P))
      if(covariances) Covariances <- cholsolveAQinvAT(A = CPM,
                                                      Lp = chol_Qx$Qpermchol,
                                                      P = chol_Qx$P)
    }
    
    # Compute standard error
    newdata[["sd"]] <- sqrt(newdata[["var"]])
    
    ## Return the prediction polygons
    if(covariances) {
      newdata <- list(newdata = newdata,
                      Cov = Covariances)
    }
    return(newdata)
  }
}


.SRE.predict_TMB <- function(object, newdata, CP, predict_BAUs, pred_time, type, nsim, 
                             obs_fs, k, percentiles, kriging) {
  
  obsidx <- observed_BAUs(object)        # index of observed BAUs
  p      <- length(object@alphahat)       # number of fixed regression effects
  mstar  <- length(obsidx)                # number of observed BAUs
  r      <- nbasis(object)                # number of basis-function coefficients
  
  ## Extract covariates from BAUs
  X <- as(.extract_BAU_X_matrix(object@f, object@BAUs), "matrix") 
  
  ## If we are doing universal kriging, use the full joint precision 
  ## matrix of the fixed and random effects. If we are doing simple kriging, 
  ## use only the random effect block of the precision matrix.
  if (kriging == "universal") {
    Q_posterior <- object@Q_posterior
  } else if (kriging == "simple") {
    if (object@simple_kriging_fixed) {
      Q_posterior <- object@Q_posterior
    } else {
      Q_posterior <- object@Q_posterior[-(1:p), -(1:p)]
    }
  } 
  
  ## Compute the Cholesky factor of the permuted precision matrix.
  Q_L <- sparseinv::cholPermute(Q = Q_posterior)
  
  ## Generate Monte Carlo samples at all BAUs (or over arbitrary prediction 
  ## regions given in the argument newdata)
  MC <- .simulate(object = object, X = X, type = type, obs_fs = obs_fs, 
                  nsim = nsim, k = k, Q_L = Q_L,  
                  predict_BAUs = predict_BAUs, CP = CP, 
                  kriging = kriging, newdata = newdata)
  
  ## We do not allow aggregation of the Y-process when predicting over arbitrary polygons
  if(!predict_BAUs) MC$Y_samples <- NULL
  
  ## Remove other quantities if the user has not requested them
  if(!("link" %in% type)) MC$Y_samples <- NULL
  
  if(!("mean" %in% type)) MC$mu_samples <- MC$prob_samples <- NULL
  
  ## Produce prediction and RMSPE matrices. 
  ## The columns are the quantity of interest (Y, mu, prob, or Z), and the rows are prediction locations.
  predictions <- sapply(MC, rowMeans)
  RMSPE <- sapply(MC, apply, 1, sd)
  
  ## If we are predicting over BAUs, newdata is NULL, so set it to the BAUs.
  ## Note that BAUs must be a Spatial*DataFrame, so coercion isn't really necessary.
  if (predict_BAUs) newdata <- .Coerce_SpatialDataFrame(object@BAUs)
  
  ## Now update newdata with the predictions, RMSPE, and percentiles. 
  ## (See https://datascience.stackexchange.com/a/8924 for a description of what this gsub is doing.)
  QOI <- gsub("_.*", "", names(MC)) # Quantities Of Interest
  
  ## Predictions and RMSPE
  colnames(predictions) <- paste0("p_", QOI)
  colnames(RMSPE) <- paste0("RMSPE_", QOI)
  newdata@data <- cbind(newdata@data, predictions, RMSPE)
  
  ## Percentiles 
  newdata@data <- .concat_percentiles_to_df(newdata@data, MC = MC, percentiles = percentiles)
  
  ## subset based on pred_time
  if(!is.null(pred_time)) {
    newdata <- newdata[,pred_time]  # return only specified time points
    ## Also need to subset the Monte Carlo samples. Only do this if t is present 
    ## in the data; otherwise, don't worry about it, as the user can just do it manually.
    if (!is.null(newdata@data$t)) {
      idx <- which(newdata@data$t %in% pred_time)
      MC <- lapply(MC, function(X) X[idx, ])
    }
  }
  
  ## Return the predictions and uncertainty summary (in newdata) and the MC samples (in MC)
  return(list(newdata = newdata, MC = MC))
}


simulate_xi <- function(object, nsim, type) {

  ## number of spatial and temporal BAUs
  if (is(object@basis,"TensorP_Basis")) {
    ns <- length(object@BAUs@sp)
    nt <- length(unique(object@BAUs@endTime))
  } else {
    ns <- length(object@BAUs)
  }

  N     <- nrow(object@S0)
  mstar <- length(observed_BAUs(object))

  if (type == "observed") {
    f <- observed_BAUs
    nloc <- mstar
  } else {
    f <- unobserved_BAUs
    nloc <- N - mstar
  }


  if (object@fs_by_spatial_BAU) {
    idx <- f(object)
    spatial_BAU_id <- ((idx - 1) %% ns) + 1
    sigma2fs <- object@sigma2fshat[spatial_BAU_id]
  } else {
    sigma2fs <- object@sigma2fshat
  }

  xi <- rnorm(nsim * nloc, mean = 0, sd = sqrt(sigma2fs)) %>%
    matrix(nrow = nloc, ncol = nsim)

  return(xi)
}



.simulate <- function(object, X, type, nsim, obs_fs, k, Q_L, predict_BAUs, CP, kriging, newdata, fixed_effects_only = FALSE, fixed_effects_and_gamma = FALSE){
  
  ## Design matrices evaluated at observed BAUs only
  X_O <- .constructX_O(object) 
  S_O <- .constructS_O(object) 
  if (object@include_gamma) {
    G_O <- .constructG_O(object)
    G_O <- .dropzerocolumns(G_O) # Drop empty columns from G_O (some levels may be unobserved) 
    g_O <- sum(sapply(G_O, ncol))  # number of observed random effects
  }
  
  obsidx <- observed_BAUs(object)        # index of observed BAUs
  
  MC <- list()              
  N   <- nrow(object@S0)   
  mstar <- length(obsidx)
  r   <- nbasis(object)   
  
  ## Number of fixed and random effects
  p <- length(object@alphahat)
  
  # ---- Define alpha and generate samples from (eta', xi_O')' ----
  
  ## Construct the posterior mean of the regression (if kriging = "universal"),
  ## the basis-function random weights, the random effects gamma 
  ## (if object@include_gamma) and the fine-scale variation (if object@include_fs). 
  posterior_mean <- if (kriging == "universal") as.numeric(object@alphahat) else vector()
  posterior_mean <- c(posterior_mean, as.numeric(object@mu_eta))
  if (object@include_gamma) posterior_mean <- c(posterior_mean, as.numeric(object@mu_gamma))
  if (object@include_fs) posterior_mean <- c(posterior_mean, as.numeric(object@mu_xi))
  
  ## Generate samples from Gau(0, 1) distribution, and transform this 
  ## standard normal vector to one that has the desired posterior precision matrix.
  z <- matrix(rnorm(length(posterior_mean) * nsim), 
              nrow = length(posterior_mean), ncol = nsim)
  U <- Matrix::t(Q_L$Qpermchol) # upper Cholesky factor of permuted joint posterior precision matrix 
  x <- solve(U, z)              # x ~ Gau(0, A), where A is the permuted precision matrix i.e. A = P'QP
  y <- Q_L$P %*% x              # y ~ Gau(0, Q^-1)
  y <- as.matrix(y)
  
  ## Add the mean vector to the simulated random effects.
  samples <- y + matrix(rep(posterior_mean, times = nsim), ncol = nsim) 
  
  ## Separate the MC samples of each quantity 
  if (kriging == "universal") {
    alpha <- samples[1:p, , drop = FALSE]
    eta   <- samples[(p + 1):(p + r), ]
    if (object@include_gamma) {
      gamma_O <- samples[(p + r + 1):(p + r + g_O), ]
    }
  } else {
    eta <- samples[1:r, ]
    if (object@include_gamma) {
      gamma_O <- samples[(r + 1):(r + g_O), ]
    }
  }
  
  # Return the MC samples of the fixed effects if that's all we care about
  if (fixed_effects_only) return(alpha)
  
  ## Observed fine-scale variation
  if (object@include_fs) {
    if (kriging == "universal") {
      xi_O  <- samples[(p + r + 1):(p + r + mstar), ]
    } else {
      xi_O  <- samples[(r + 1):(r + mstar), ]
    }
  }
  
  ## We now have several matrices of Monte Carlo samples:
  ## alpha (if kriging = "universal"), 
  ## eta, 
  ## gamma_O,
  ## and xi_O (if include_fs = TRUE).
  ## row i of eta corresponds to nsim MC samples of eta_i,
  ## row i of xi_O corresponds to nsim MC samples of the fine-scale variation at the ith observed location.
  
  # ---- Generate samples from xi_U ----
  
  ## This is straightforward as each element of xi_U is independent of
  ## all other random effects in the model.
  ## All we have to do is make an (N-m) x nsim matrix of draws from the
  ## Gaussian distribution with mean zero and variance equal to the fine-scale variance.
  if (object@include_fs) {
    xi_U <- simulate_xi(object, nsim = nsim, type = "unobserved")
    xi_samples <- rbind(xi_O, xi_U)
  }
  
  # ---- Generate samples from gamma_U ----
  
  ## Now sample the unobserved random effects. This is easy as we simply need to 
  ## sample from a mean-zero Gaussian distribution with variance object@sigma2gamma
  
  if (object@include_gamma) {
    G0   <- object@G0
    B    <- length(G0)         # number of random-effect groups
    lb   <- sapply(G0, ncol)   # total number of levels for each random-effect group
    lb_O <- sapply(G_O, ncol)  # number of observed levels for each random-effect group
    lb_U <- lb - lb_O          # number of unobserved levels for each random-effect group
    
    gamma <- lapply(1:B, function(b) {
      ## Extract rows of gamma_O that correspond to random-effect group b
      cs  <- cumsum(lb_O) 
      idx <- (cs[b] - lb_O[b] + 1):cs[b] 
      gamma_O <- gamma_O[idx, ]
      
      ## Simulate unobserved rndom effects for group b  
      gamma_U <- matrix(rnorm(lb_U * nsim, sd = sqrt(object@sigma2gamma[b])), nrow = lb_U, ncol = nsim)
      
      ## Merge gamma_U and gamma_O. The ordering is very important, and needs to 
      ## reflect the way that G_O during model fitting was constructed. 
      G0 <- G0[[b]]
      G0 <- G0[obsidx, ]
      idx_U <- which(apply(G0, 2, sum) == 0)
      gamma <- matrix(NA, nrow = lb[b], ncol = nsim)
      # Now fill the rows of gamma, filling it with the rows of gamma_O or gamma_U
      counter_O = 1
      counter_U = 1
      for (i in 1:lb[b]) {
        if (i %in% idx_U) {
          gamma[i, ] <- gamma_U[counter_U, ]
          counter_U <- counter_U + 1
        } else {
          gamma[i, ] <- gamma_O[counter_O, ]
          counter_O <- counter_O + 1
        }
      }
      
      return(gamma)
    })
    gamma <- do.call(rbind, gamma)
  }
  
  # sanity check: apply(gamma, 1, var)
  if (fixed_effects_and_gamma) return(list(alpha = alpha, gamma = gamma))
  
  ## Split the design matrices based on observed and unobserved BAUs
  X_U <- X[-obsidx, ]         # Unobserved fixed effect 'design' matrix
  S_U <- object@S0[-obsidx, ] # Unobserved random effect 'design' matrix

  # ---- Construct samples from the latent process Y ----
  
  ## We break the latent process down as: Y = Y_smooth + xi, 
  ## so that we may separate the fine-scale variation. 
  
  ## Split the design matrices based on observed and unobserved BAUs
  X_U <- X[-obsidx, ]    # Unobserved fixed effect 'design' matrix
  S_U <- object@S0[-obsidx, ] # Unobserved random effect 'design' matrix
  
  ## Simulate the smooth Y-process (excluding fs variation) over the observed and unobserved BAUs
  if (kriging == "universal") {
    Y_smooth_O <- X_O %*% alpha + S_O %*% eta
    Y_smooth_U <- X_U %*% alpha + S_U %*% eta
  } else {
    Y_smooth_O <- X_O %*% object@alphahat + S_O %*% eta
    Y_smooth_U <- X_U %*% object@alphahat + S_U %*% eta
  }
  
  ## Combine samples
  Y_smooth_samples  <- rbind(Y_smooth_O, Y_smooth_U)
  
  ## Use permutation matrix to get the correct (original) ordering in terms of the BAUs
  unobsidx         <- unobserved_BAUs(object)   # Unobserved BAUs indices
  ids              <- c(obsidx, unobsidx)       # All indices (observed and unobserved)
  P                <- Matrix::sparseMatrix(i = 1:N, j = 1:N, x = 1)[ids, ]
  Y_smooth_samples <- Matrix::t(P) %*% Y_smooth_samples
  
  ## Add the random effects 
  if (object@include_gamma) {
    G <- do.call(cbind, G0)
    Y_smooth_samples = Y_smooth_samples + G %*% gamma
  }
  
  ## Construct the samples from the latent process Y 
  if (object@include_fs) {
    xi_samples  <- Matrix::t(P) %*% xi_samples
    Y_samples   <- Y_smooth_samples + xi_samples
  } else {
    Y_samples <- Y_smooth_samples
  }
  Y_samples        <- as.matrix(Y_samples)
  Y_smooth_samples <- as.matrix(Y_smooth_samples)
  
  ## Outputted Y value depend on obs_fs
  MC$Y_samples <- if (obs_fs) Y_smooth_samples else Y_samples
  
  ## If Y is the only quantity of interest, exit the function.
  if (!("mean" %in% type) & !("response" %in% type)) return(MC) 
  
  
  # ---- Apply inverse-link function to the samples to obtain conditional mean ----

  ## For all models other than the Gau-Gau model (i.e., Gaussian data with an 
  ## identity link function), the fine-scale variation process xi(.) is 
  ## henceforth included in the latent process Y(.).
  
  if (object@response == "gaussian" && object@link == "identity" && obs_fs) {
    mu_samples   <- Y_smooth_samples
  } else {
    tmp <- .inverselink(Y_samples, object, k = k)
    mu_samples   <- tmp$mu_samples
    prob_samples <- tmp$prob_samples
  }
  
  # ---- Predicting over arbitrary polygons ----
  
  if (!predict_BAUs) mu_samples <- as.matrix(CP %*% mu_samples)
  
  if (!predict_BAUs & object@response %in% c("binomial", "negative-binomial")) {
    k_P <- CP %*% k
    h   <- .link_fn("mu_to_prob", response = object@response)
    prob_samples <- as.matrix(h(mu_samples, k_P))
  }
  
  ## Output the mean samples. If probability parameter was computed, and 
  ## we are predicting over the BAUs, also output.
  MC$mu_samples <- mu_samples
  if (!is.null(prob_samples)) MC$prob_samples <- prob_samples
  
  ## If the response is not a quantity of interest, exit the function
  if (!("response" %in% type)) return(MC)
  
  
  # ---- Sample the response variable, Z ----
  
  if (object@response == "gaussian") {
    sigma2e <- .measurementerrorvariance(object, newdata)
    if (object@link == "identity" && obs_fs) {
      
      ## Add the fine-scale variation to mu(.) and possibly aggregate over 
      ## the prediction regions. Note that the mu_samples return in the list MC 
      ## is still the smooth version without fine-scale variation; we just add 
      ## it here so that it can be implicitly incorporated in the simulations of 
      ## the response variable. 
      mu_samples  <- Y_smooth_samples + xi_samples
      if (!predict_BAUs) mu_samples <- as.matrix(CP %*% mu_samples)
    }
  } 
  
  MC$Z_samples <- .sampleZ(mu_samples = mu_samples, object = object, obs_fs = obs_fs, sigma2e = sigma2e, k = k)
  
  return(MC)
}


.inverselink <- function(Y_samples, object, k) {
  
  ## For families with a known constant parameter (binomial, negative-binomial),
  ## finv() maps the Gaussian scale Y process to the probability parameter p.
  ## Then, we map p to the conditional mean mu via hinv().
  ## For all other families, ginv() maps Y directly to mu.
  ## The exception is negative-binomial with a log or square-root link, 
  ## in which case we map directly from Y to mu.
  
  ## Note that for all cases other than type == "link", we need to compute the conditional mean samples.
  
  ## Create the relevant link functions.
  if (object@response %in% c("binomial", "negative-binomial") & object@link %in% c("logit", "probit", "cloglog")) {
    finv  <- .link_fn("Y_to_prob", link = object@link)
  } else {
    ginv  <- .link_fn("Y_to_mu", link = object@link) 
  }
  
  if(object@response %in% c("binomial", "negative-binomial")) {
    hinv  <- .link_fn("prob_to_mu", response = object@response)
    h     <- .link_fn("mu_to_prob", response = object@response)
  }
  
  ## Create the mu samples (and prob parameter if applicable)
  if (object@response %in% c("binomial", "negative-binomial") & object@link %in% c("logit", "probit", "cloglog")) {
    prob_samples <- finv(Y = Y_samples)
    mu_samples   <- hinv(p = prob_samples, k = k)
  } else if (object@response == "negative-binomial" & object@link %in% c("log", "sqrt")) {
    mu_samples   <- k * ginv(Y_samples) 
    prob_samples <- h(mu = mu_samples, k = k)
  } else {
    mu_samples   <- ginv(Y_samples)
  }
  
  if (!exists("prob_samples")) prob_samples <- NULL
  
  return(list(mu_samples = mu_samples, prob_samples = prob_samples))
}



.sampleZ <- function(mu_samples, object, obs_fs, sigma2e, k) {
  
  n    <- length(mu_samples)
  nsim <- ncol(mu_samples)
  
  if (object@response == "poisson") {
    Z_samples <- rpois(n, lambda = c(t(mu_samples)))
  } else if (object@response == "gaussian") { 
    Z_samples <- rnorm(n, mean = c(t(mu_samples)), sd = sqrt(sigma2e))
  } else if (object@response == "gamma") {
    theta <- 1 / c(t(mu_samples))    # canonical parameter
    a     <- 1/object@phi            # shape parameter
    beta  <- theta * a               # rate parameter (1/scale)
    Z_samples <- rgamma(n, shape = a, rate = beta)
  } else if (object@response == "inverse-gaussian") {
    Z_samples <- statmod::rinvgauss(n, mean = c(t(mu_samples)), dispersion = object@phi) 
  } else if (object@response == "negative-binomial") {
    k_vec <- rep(k, each = nsim)
    Z_samples <- rnbinom(n, size = k_vec, mu = c(t(mu_samples)))
  } else if (object@response == "binomial") {
    k_vec <- rep(k, each = nsim)
    theta <- log((c(t(mu_samples))/k_vec) / (1 - (c(t(mu_samples))/k_vec)))
    p <- 1 / (1 + exp(-theta))
    ## NAs will occur if k = 0. Fortunately, if k = 0, we know Z will be 0. 
    ## Hence, simply replace the NA occurences in p with 0.
    p[is.na(p)] <- 0
    Z_samples <- rbinom(n, size = k_vec, prob = p)
  }
  
  ## Convert from a long vector to an n_pred_locs x nsim matrix
  Z_samples <- matrix(Z_samples, ncol = nsim, byrow = TRUE)
}

.measurementerrorvariance <- function(object, newdata) {
  
  if (!is.null(newdata$std)) {
    sigma2e <- (newdata$std)^2
  } else {
    
    # Ideally we would use all of Ve. However, Ve is of dimension m*, where m* 
    # is the number of observed BAUs (rather than m, the number of 
    # observations in newdata pre-binning). For now, just extract a single value 
    # from Ve. 
    
    if (!.zero_range(diag(object@Ve))) {
      sigma2e <- mean(diag(object@Ve))
      cat("Prediction of Gaussian-distributed data requires the measurement-error standard deviation; we are assuming that it is spatially invariant, and that it is the average of measurement-error standard deviation used during model fitting (which were supplied with the data in the field 'std', or estimated using variogram techniques).")
    } else {
      sigma2e <- object@Ve[1, 1]
    }
  }
  
  return(sigma2e)
}
