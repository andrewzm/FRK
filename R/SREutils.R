#' @rdname loglik
#' @export
setMethod("loglik", signature="SRE", function(object) {
  .Deprecated("logLik")
  return(logLik(object))
  
})

setMethod("binned_data", signature = "SRE", function(object) {
  # This will contain the binned data for each BAU, and will
  # feature many missing values in the case of point-referenced data.
  Cmat_dgT <- .as(object@Cmat, "dgTMatrix")
  obs_BAUs <- Cmat_dgT@j + 1 # BAUs associated with each observation
  data_idx <- Cmat_dgT@i + 1 # data index
  Z <- rep(NA, length(object@BAUs))
  Z[obs_BAUs] <- object@Z[data_idx] # TODO delete the following if the new method works without errors: Z[obs_BAUs] <- object@Z[, 1][data_idx]
  return(Z)
})

#' @rdname SRE
#' @export
setMethod("logLik", signature="SRE", function(object) {
  # This is structured this way so that extra models for fs-variation
  # can be implemented later
  if (length(object@log_likelihood) != 0) {
    return(object@log_likelihood)
  } else if (object@fs_model == "ind") {
    return(.loglik.ind(object))
  } else {
    stop("Currently only independent fine-scale model is implemented")
  }
})

#' @rdname SRE
#' @export
setMethod("nobs", signature="SRE", function(object, ...) length(object@Z))

#' @rdname observed_BAUs
#' @export
setMethod("observed_BAUs", signature(object = "SRE"), function (object) {
  return(object@obsidx)
})


## this function is defined so that we can call it in SRE(), before an object of 
## class SRE is defined
## Note that Cmat maps BAUs to the observations. The dimension of object@Cmat is
## (number of observations) * (number of BAUs).
.observed_BAUs_from_Cmat <- function(Cmat) { 
  return(unique(.as(Cmat, "dgTMatrix")@j) + 1)
}

#' @rdname observed_BAUs
#' @export
setMethod("unobserved_BAUs",signature(object = "SRE"), function (object) {
  ## ncol(object@Cmat) is the total number of BAUs:
  unobsidx <- (1:ncol(object@Cmat))[-object@obsidx]
  return(unobsidx)
})


## Print/Show SRE
print.SRE <- function(x,...) {
  print_list <- list(
    formula =  deparse(x@f),
    ndatasets = length(x@data),
    nbasis = x@basis@n,
    basis_class = class(x@basis)[1],
    nBAUs = length(x@BAUs),
    nobs = length(x@Z),
    mean_obsvar = mean(x@Ve@x), 
    fs_var = x@sigma2fshat,
    dimC = deparse(dim(x@Cmat)),
    dimS = deparse(dim(x@S)),
    ncovars = ncol(x@X), 
    response = x@response, 
    link = x@link, 
    method = x@method)
  
  cat("Formula:",print_list$formula,"\n")
  cat("Assumed response distribution:",print_list$response,"\n")
  cat("Specified link function:",print_list$link,"\n")
  cat("Method of model fitting:",print_list$method,"\n")
  cat("Number of datasets:",print_list$ndatasets,"\n")
  cat("Number of basis functions:",print_list$nbasis,"\n")
  cat("Class of basis functions:",print_list$basis_class,"\n")
  cat("Number of BAUs [extract using object@BAUs]: ",print_list$nBAUs,"\n")
  cat("Number of observations [extract using object@Z]: ",print_list$nobs,"\n")
  cat("Mean obs. variance at BAU level [extract using object@Ve]:",print_list$mean_obsvar,"\n")
  cat("Fine-scale variance proportionality constant [extract using object@sigma2fshat]:",print_list$fs_var,"\n")
  cat("Dimensions of C in Z = C*Y + e [extract using object@Cmat]: ",print_list$dimC,"\n")
  cat("Dimensions of S in Y = X*alpha + S*eta + delta [extract using object@S]: ",print_list$dimS,"\n")
  cat("Number of covariates:",print_list$ncovars,"\n\n")
}
setMethod("show",signature(object="SRE"),function(object) print(object))

## Summary SRE
summary.SRE <- function(object,...) {
  summ_list <- list(
    summ_Kdiag = summary(diag(object@Khat)),
    summ_mueta = summary(object@mu_eta[,1]),
    summ_vareta = summary(diag(object@S_eta)),
    summ_muxi = summary(object@mu_xi[, 1]),
    finescale_var = deparse(as.vector(object@sigma2fshat)),
    reg_coeff = deparse(as.vector(object@alphahat)) 
  )
  class(summ_list) <- "summary.SRE"
  summ_list
}
setMethod("summary",signature(object="SRE"),summary.SRE)

## Set method for retrieving info_fit
#' @rdname info_fit
#' @export
setMethod("info_fit", signature(object = "SRE"),
          function(object) {object@info_fit})


# Retrieve coefficients of SRE model
#' @rdname SRE
#' @export 
setMethod("coef",signature(object = "SRE"),function(object,...) {
  coeff <- as.numeric(object@alphahat)
  varnames <- all.vars(object@f)[-1]
  nms <- "Intercept"
  if(length(varnames) > 0) nms <- c(nms, varnames)
  names(coeff) <- nms
  return(coeff)
})


# Retrieve uncertainty quantification on the coefficients of SRE model
#' @rdname SRE
#' @export
setMethod("coef_uncertainty",signature(object = "SRE"),function(object, percentiles = c(5, 95),  nsim = 400) {
  
  if(object@simple_kriging_fixed) stop("Uncertainty quantification on the fixed effects is only possible if simple_kriging_fixed is set to FALSE in SRE.fit()")
  
  ## Compute the Cholesky factor of the permuted precision matrix.
  Q_L <- sparseinv::cholPermute(Q = object@Q_posterior)
  
  ## Generate sample of regression coefficients (most of these arguments can be 
  ## NULL because the function exits before they are used)
  alpha <- .simulate(
    object = object, nsim = nsim, X = NULL, G=NULL, Q_L = Q_L, fixed_effects_only = TRUE, kriging = "universal", 
    type = NULL, obs_fs = NULL, k = NULL,predict_BAUs = NULL, CP = NULL, newdata = NULL
  )
  
  x <- apply(alpha, 1, quantile, percentiles / 100)
  
  ## Add names to the percentiles
  varnames <- all.vars(object@f)[-1]
  nms <- "Intercept"
  if (length(varnames) > 0) nms <- c(nms, varnames)
  
  if (is.matrix(x)) {
    colnames(x) <- nms
  } else {
    names(x) <- nms
  }
  
  return(x)
})

## Print summary of SRE
print.summary.SRE <- function(x, ...) {
  cat("Summary of Var(eta) [extract using object@Khat]: \n")
  print(x$summ_Kdiag)
  cat("\n")
  cat("Summary of E(eta | Z) [extract using object@mu_eta]: \n")
  print(x$summ_mueta)
  cat("\n")
  cat("Summary of Var(eta | Z) [extract using object@S_eta]: \n")
  print(x$summ_vareta)
  cat("\n")    
  cat("Summary of E(xi | Z) [extract using object@mu_xi]: \n")
  print(x$summ_muxi)
  cat("\n")    
  cat("Fine-scale variance estimate [extract using object@sigma2fshat]:", x$finescale_var, "\n")
  cat("Regression coefficients [extract using object@alpha]:",x$reg_coeff,"\n")
  cat("For object properties use show().\n")
  invisible(x)
}



## The function below is used to facilitate the computation of multiple marginal variances
## by splitting up the problem into batches (prediction is an embarassingly parallel procedure)
.batch_compute_var <- function(S0,Cov,fs_in_process = TRUE) {
  
  ## Don't consider more than 1e4 elements at a time
  batch_size <- 1e4
  
  ## Create batching indices by simply dividing the n rows into batches of size 10000
  batching=cut(1:nrow(S0),
               breaks = seq(0,nrow(S0)+batch_size,
                            by=batch_size),labels=F)
  r <- ncol(S0) # number of columns
  
  ## At first this was parallelised, but the memory demand was becoming too large in many instances
  ## So for now parallelism is disabled. The following line can be uncommented if we wish to
  ## re-enable parallelism in the predictions
  
  #if(opts_FRK$get("parallel") > 1 & batch_size < nrow(X)) {
  if(0) {
    ## Export variables to the cluster
    clusterExport(opts_FRK$get("cl"),
                  c("batching","S0","Cov"),envir=environment())
    
    ## Compute the variances over max(10000) BAUs per batch
    var_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                          function(i) {
                            ## Find which BAUs this batch predicts over
                            idx = which(batching == i)
                            
                            ## Compute the marginal variance for these BAUs
                            if(fs_in_process) {
                              rowSums((S0[idx,] %*% Cov[1:r,1:r]) * S0[idx,]) +  # eta contribution
                                diag(Cov)[-(1:r)][idx] +                       # fs contribution
                                2*rowSums(Cov[r+idx,1:r] * S0[idx,])           # cross.cov between eta
                              # and fs variation
                            } else {
                              rowSums((S0[idx,] %*% Cov) * S0[idx,])             # just eta contribution
                            }
                          })
    clusterEvalQ(opts_FRK$get("cl"), {gc()})   # clear memory of cluster
    temp <- do.call(c,var_list)                # concatenate results
  } else {
    temp <- rep(0,nrow(S0))                    # initialise the vector of variances
    for(i in 1:max(unique(batching))) {        # for each batch
      
      ## Find which BAUs this batch predicts over
      idx = which(batching==i)
      
      ## If we have fs variation in the process
      if(fs_in_process)
        temp[idx] <- rowSums((S0[idx,] %*% Cov[1:r,1:r]) * S0[idx,]) +  # eta contribution
          diag(Cov)[-(1:r)][idx] +                                    # fs contribution
          2*rowSums(Cov[(r+idx),1:r] * S0[idx,])                      # cross-cov between eta
      else                                                                # and fs variation
        temp[idx] <- rowSums((S0[idx,] %*% Cov) * S0[idx,])             # otherwise we just have the eta
      # contribution
      
    }
  }
  
  ## Return all variances
  temp
}


## This function attempts to estimate the measurement error by fitting a variogram to the data
## and see where it crosses the y-axis. This captures the super-fine-scale variation that we
## characterise as measurement error. FRK then effectively fits a smooth variogram -- the difference
## between where this crosses the y-axis and the measurement error will be the fs-variation (estimated)
.est_obs_error <- function(sp_pts,variogram.formula,vgm_model = NULL,BAU_width = NULL) {
  
  ## Notify user (even if not verbose == TRUE)
  ## (after revision, now only notify user if verbose == TRUE)
  if(opts_FRK$get("verbose")) cat("Fitting variogram for estimating measurement error...\n")
  
  ## Basic checks
  if(!is(variogram.formula,"formula"))
    stop("variogram.formula needs to be of class formula")
  if(!is(sp_pts,"Spatial"))
    stop("sp_pts needs to be of class Spatial")
  if(!requireNamespace("gstat"))
    stop("gstat is required for variogram estimation. Please install gstat")
  
  ## Make sure we're not on sphere here, otherwise variogram fitting is too slow. Just remove
  ## CRS (this is only approximate anyway)
  if(!is.na(proj4string(sp_pts))) {
    sp_pts <- SpatialPointsDataFrame(coords = coordinates(sp_pts),
                                     data = sp_pts@data,proj4string = CRS())
  }
  
  
  ## If any of the coordinates are zero range (i.e., they contain only the same 
  ## value repeatedly), then we need to jitter the coordinates so that gstat can 
  ## estimate the measurement error using the variogram. 
  ## This situation typically only arises if we are trying to fit 1D spatial data.
  zero_range_dims <- which(apply(coordinates(sp_pts), 2, .zero_range))
  if (length(zero_range_dims) > 0) {
    new_coords <- coordinates(sp_pts)
    for (i in zero_range_dims) {
      new_coords[, i] <- jitter(new_coords[, i], amount = 1e-5)
    }
    
    ## An error occurs if we try to overwrite the original coordinates:
    # coordinates(sp_pts) <- new_coords 
    ## Just create a new Spatial* object:
    form <- as.formula(paste("~", paste(coordnames(sp_pts), collapse = "+")))
    sp_pts <- cbind(sp_pts@data, new_coords)
    coordinates(sp_pts) <- form
  }
  

  
  ## If we have many points (say > 50000) then subsample
  if(length(sp_pts) > 50000) {
    if(opts_FRK$get("verbose") > 0)
      cat("Selecting 50000 data points at random for estimating the measurement error variance\n")
    sp_pts_sub <- sp_pts[sample(1:length(sp_pts),50000),]
  } else sp_pts_sub <- sp_pts
  
  ## Find the maximum extent in each dimension
  coords_range <-  apply(coordinates(sp_pts_sub),2,
                         function(x) diff(range(x)))
  
  ## Find a maximum effective length scale by computing the "diagonal"
  diag_length <- sqrt(sum(coords_range^2))
  
  ## Compute the area of the domain
  area <- prod(coords_range)
  
  ## Consider the area that contains about 100 data points in it (we only want
  ## to study variogram points close to the origin)
  cutoff <- sqrt(area * 100 / length(sp_pts_sub))
  
  ## Remove covariates since we are only interested at behaviour close to iring
  ## and since the binning into BAUs hasn't (shouldn't have) occurred yet
  variogram.formula <- .formula_no_covars(variogram.formula)
  
  ## Extract  data values from data object
  L <- .extract.from.formula(variogram.formula,data=sp_pts_sub)
  
  ## Create a gstat object with this formula and data
  g <- gstat::gstat(formula=variogram.formula,data=sp_pts_sub)
  
  ## Compute the empirical variogram
  v <- gstat::variogram(g,cressie=T,            # Cressie's robust variogram estimate
                        cutoff = cutoff,        # maximum spatial separation distance
                        width = cutoff/10)      # width for semivariance estimates
  
  ## Fit the model. First, if the user did not supply any desired model, try to fit a linear model
  ## with initial conditions as given to vgm()
  if(is.null(vgm_model))
    vgm_model <-  gstat::vgm(psill = var(L$y)/2,
                             model = "Lin",
                             range = mean(v$dist),
                             nugget = var(L$y)/2)
  
  ## Try to fit the model.
  
  ## Try fitting using a linear model on just the first two to
  ## four points (cf. Kang and Cressie)
  OK <- FALSE                                    # init. OK
  num_points <- 2                                # start with 2 points
  while(!OK & num_points <= 4) {                 # while NOT OK
    ## fit only to first 1:num_points points
    linfit <- lm(gamma~dist,data=v[1:num_points,])
    
    ## ans. will be returned in vgm.fit$psill later on
    vgm.fit <- list(singular = 0)
    vgm.fit$psill <- coefficients(linfit)[1]   # extract psill from intercept
    if(vgm.fit$psill > 0) OK <- TRUE           # if variance > 0 OK = TRUE
    num_points <- num_points + 1               # increment num_points
  }
  
  
  ## If we have estimated a negative variance, then try to fit a linear variogram
  ## using more points
  if(vgm.fit$psill[1] <= 0) {
    
    vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))
    
    ## Check if the process of fitting generates a warning. If it did then OK == 0
    ## otherwise OK == 1 (this fits twice but since it's so quick it's not an issue)
    OK <- tryCatch({vgm.fit <- gstat::fit.variogram(v, model = vgm_model); 1},
                   warning=function(w) 0)
    ## If the reportef psill is less or equal to zero, or fit.variogram reported a singularity,
    ## or a Warning was thrown, then retry fitting using an exponential model
    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
      vgm_model <-  gstat::vgm(psill = var(L$y)/2,
                               model = "Exp",
                               range = mean(v$dist),
                               nugget = var(L$y)/2)
      OK <- tryCatch({vgm.fit = gstat::fit.variogram(v, model = vgm_model); OK <- 1},warning=function(w) 0)
      vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))
    }
    
    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
      ## Try with Gaussian, maybe process is very smooth or data has a large support
      vgm_model <-  gstat::vgm(var(L$y)/2, "Gau", mean(v$dist), var(L$y)/2)
      
      ## Try to fit the model.
      vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))
      
      ## Like above, we return OK = 0 if fit is still not good
      OK <- tryCatch({vgm.fit = gstat::fit.variogram(v, model = vgm_model); OK <- 1},warning=function(w) 0)
    }
    
    ## If we still have problems, then just take the first point of the empirical semivariogram and
    ## throw a warning that this estimate is probably not very good
    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
      vgm.fit$psill[1] <- v$gamma[1]
      warning("Estimate of measurement error is probably inaccurate.
                    Please consider setting it through the std variable
                    in the data object if known.")
    }
  }
  
  ## Return the sqrt of the psill as the measurement error
  sp_pts$std <- sqrt(vgm.fit$psill[1])
  sp_pts
}