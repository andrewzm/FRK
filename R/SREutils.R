# FRK: An R Software package for spatial and spatio-temporal prediction
# with large datasets.
# Copyright (c) 2017 University of Wollongong
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @name SRE
#' @title Construct SRE object, fit and predict
#' @description The Spatial Random Effects (SRE) model is the central object in FRK. The function \code{FRK} provides a wrapper for the construction and estimation of the SRE object from data, using the functions \code{SRE} (the object constructor) and \code{SRE.fit} (for fitting it to the data). Please see \code{\link{SRE-class}} for more details on the SRE object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame}, \code{SpatialPolygonsDataFrame}, \code{STIDF}, or  \code{STFDF}. If using space-time objects, the data frame must have another field, \code{t}, containing the time index of the data point. If the assumed response distribution is \code{"binomial"} or \code{"negative-binomial"}, the data frame must have another field, \code{k}, containing the known constant parameter \eqn{k} for each observation. 
#' @param basis object of class \code{Basis} (or \code{TensorP_Basis})
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame}, \code{STIDF}, or \code{STFDF}. The object's data frame must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality. If the function \code{FRK} is used directly, then BAUs are created automatically, but only coordinates can then be used as covariates
#' @param est_error flag indicating whether the measurement-error variance should be estimated from variogram techniques (only applicable for a Gaussian response). If this is set to 0, then \code{data} must contain a field \code{std}. Measurement-error estimation is currently not implemented for spatio-temporal datasets
#' @param include_fs flag indicating whether the fine-scale variation should be include in the model
#' @param average_in_BAU if \code{TRUE}, then multiple data points falling in the same BAU are averaged; the measurement error of the averaged data point is taken as the average of the individual measurement errors
#' @param sum_variables vector of strings indicating which variables are to be summed rather than averaged; only applicable if \code{average_in_BAU == TRUE}
#' @param normalise_wts if \code{TRUE}, the rows of the incidence matrices \eqn{C_Z} and \eqn{C_P} are normalised to sum to 1, so that the mapping represents a weighted average; if false, no normalisation of the weights occurs (i.e., the mapping corresponds to a weighted sum)
#' @param fs_by_spatial_BAU if \code{TRUE}, then each spatial BAU is associated with its own fine-scale variance parameter (only allowed if method = 'TMB'). Otherwise, a single fine-scale variance parameter is used
#' @param fs_model if "ind" then the fine-scale variation is independent at the BAU level. If "ICAR", then an ICAR model for the fine-scale variation is placed on the BAUs
#' @param vgm_model an object of class \code{variogramModel} from the package \code{gstat} constructed using the function \code{vgm}. This object contains the variogram model that will be fit to the data. The nugget is taken as the measurement error when \code{est_error = TRUE}. If unspecified, the variogram used is \code{gstat::vgm(1, "Lin", d, 1)}, where \code{d} is approximately one third of the maximum distance between any two data points
#' @param K_type the parameterisation used for the \code{K} matrix. If the EM algorithm is used for model fitting, \code{K_type} can be "unstructured" or "block-exponential". If TMB is used for model fitting, \code{K_type} can be "neighbour" or "block-exponential". The default is "block-exponential"
#' @param normalise_basis flag indicating whether to normalise the basis functions so that they reproduce a stochastic process with approximately constant variance spatially
#' @param SRE_model object returned from the constructor \code{SRE()} containing all the parameters and information on the SRE model
#' @param n_EM maximum number of iterations for the EM algorithm
#' @param tol convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently "EM" and "TMB" are supported
#' @param lambda ridge-regression regularisation parameter for when \code{K} is unstructured (0 by default). Can be a single number, or a vector (one parameter for each resolution)
#' @param print_lik flag indicating whether likelihood value should be printed or not after convergence of the EM estimation algorithm
# #' @param use_centroid flag indicating whether the basis functions are averaged over the BAU, or whether the basis functions are evaluated at the BAUs centroid in order to construct the matrix \eqn{S}. The flag can safely be set when the basis functions are approximately constant over the BAUs in order to reduce computational time
#' @param object object of class \code{SRE}
#' @param newdata object of class \code{SpatialPoylgons}, \code{SpatialPoints}, or \code{STI}, indicating the regions or points over which prediction will be carried out. The BAUs are used if this option is not specified. 
#' @param obs_fs flag indicating whether the fine-scale variation sits in the observation model (systematic error, Case 1) or in the process model (fine-scale process variation, Case 2, default)
#' @param pred_polys deprecated. Please use \code{newdata} instead
#' @param pred_time vector of time indices at which prediction will be carried out. All time points are used if this option is not specified
#' @param covariances logical variable indicating whether prediction covariances should be returned or not. If set to \code{TRUE}, a maximum of 4000 prediction locations or polygons are allowed.
#' @param response A character string indicating the assumed distribution of the response variable. It can be "gaussian", "poisson", "bernoulli", "gamma","inverse-gaussian", "negative-binomial", or "binomial".
#' @param link A character string indicating the desired link function. Can be "log", "identity", "logit", "probit", "cloglog", "reciprocal", or "reciprocal-squared". Note that only sensible link-function and response-distribution combinations are permitted. 
#' @param taper A positive numeric indicating the strength of the covariance tapering (only applicable if \code{K_type = "block-exponential"} and \code{TMB} is used to fit the data)
#' @inheritParams .FRKTMB_pred
#' @param optimiser the optimising function used for model fitting when \code{method = 'TMB'} (default is \code{nlminb}). Users may pass in a function object or a string corresponding to a named function. Optional parameters may be passed to \code{optimiser} via \code{...}. The only requirement of \code{optimiser} is that the first three arguments correspond to the initial parameters, the objective function, and the gradient, respectively (note that this may be achieved by rearranging the order of the arguments before passing into \code{optimiser}) 
#' @param known_sigma2fs known value of the fine-scale variance. If \code{NULL} (the default), the fine-scale variance \eqn{\sigma^2_\xi} is estimated as usual. If \code{known_sigma2fs} is not \code{NULL}, the fine-scale variance is fixed to the supplied value; this may be a scalar, or vector of length equal to the number of spatial BAUs (if fs_by_spatial_BAU = TRUE)
#' @param ... other parameters passed on to \code{auto_basis} and \code{auto_BAUs} when calling \code{FRK}, or the user specified \code{optimiser} function when calling \code{FRK} or \code{SRE.fit}
#' @details \code{SRE()} is the main function in the package: It constructs a spatial random effects model from the user-defined formula, data object, basis functions and a set of Basic Areal Units (BAUs). The function first takes each object in the list \code{data} and maps it to the BAUs -- this entails binning the point-referenced data into the BAUs (and averaging within the BAU) if \code{average_in_BAU = TRUE}, and finding which BAUs are influenced by the polygon datasets. Following this, the incidence matrix \code{Cmat} is constructed, which appears in the observation model \eqn{Z = CY + C\delta + e}, where \eqn{C} is the incidence matrix and \eqn{\delta} is systematic error at the BAU level.
#'
#' The SRE model for the hidden process is given by \eqn{Y = T\alpha + S\eta + \xi}, where \eqn{T} are the covariates at the BAU level, \eqn{\alpha} are the regression coefficients, \eqn{S} are the basis functions evaluated at the BAU level, \eqn{\eta} are the basis-function coefficients, and \eqn{\xi} is the fine scale variation (at the BAU level). The covariance matrix of \eqn{\xi} is diagonal, with its diagonal elements proportional to the field `fs' in the BAUs (typically set to one). The constant of proportionality is estimated in the EM algorithm. All required matrices (\eqn{S,T} etc.) are initialised using sensible defaults and returned as part of the object, please see \code{\link{SRE-class}} for more details.
#'
#'\code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown parameters, namely the covariance matrix \eqn{K}, the fine scale variance (\eqn{\sigma^2_{\xi}} or \eqn{\sigma^2_{\delta}}, depending on whether Case 1 or Case 2 is chosen; see the vignette) and the regression parameters \eqn{\alpha}. The only method currently implemented is the Expectation Maximisation (EM) algorithm, which the user configures through \code{n_EM} and \code{tol}. The log-likelihood (given in Section 2.2 of the vignette) is evaluated at each iteration at the current parameter estimate, and convergence is assumed to have been reached when this quantity stops changing by more than \code{tol}.
#'
#'The actual computations for the E-step and M-step are relatively straightforward. The E-step contains an inverse of an \eqn{r \times r} matrix, where \code{r} is the number of basis functions which should not exceed 2000. The M-step first updates the matrix \eqn{K}, which only depends on the sufficient statistics of the basis-function coefficients \eqn{\eta}. Then, the regression parameter \eqn{\alpha} is updated and a simple optimisation routine (a line search) is used to update the fine-scale variance \eqn{\sigma^2_{\delta}} or \eqn{\sigma^2_{\xi}}. If the fine-scale errors and measurement random errors are homoscedastic, then a closed-form solution is available for the update of \eqn{\sigma^2_{\xi}} or \eqn{\sigma^2_{\delta}}. Irrespectively, since the udpates of \eqn{\alpha}, and \eqn{\sigma^2_{\delta}} or \eqn{\sigma^2_{\xi}}, are dependent, these two updates are iterated until the change in \eqn{\sigma^2_{\cdot}} is no more than 0.1\%. Information on the fitting (convergence etc.) can be extracted using \code{info_fit(SRE_model)}.
#'
#'The function \code{FRK} acts as a wrapper for the functions \code{SRE} and \code{SRE.fit}. An added advantage of using \code{FRK} directly is that it automatically generates BAUs and basis functions based on the data. Hence \code{FRK} can be called using only a list of data objects and an \code{R} formula, although the \code{R} formula can only contain space or time as covariates when BAUs are not explicitly supplied with the covariate data.
#'
#'Once the parameters are fitted, the \code{SRE} object is passed onto the function \code{predict()} in order to carry out optimal predictions over the same BAUs used to construct the SRE model with \code{SRE()}. The first part of the prediction process is to construct the matrix \eqn{S} over the prediction polygons. This is made computationally efficient by treating the prediction over polygons as that of the prediction over a combination of BAUs. This will yield valid results only if the BAUs are relatively small. Once the matrix \eqn{S} is found, a standard Gaussian inversion (through conditioning) using the estimated parameters is used for prediction.
#'
#'\code{predict} returns the BAUs, which are of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame}, or \code{STFDF}, with two added attributes, \code{mu} and \code{var}. These can then be easily plotted using \code{spplot} or \code{ggplot2} (possibly in conjunction with \code{\link{SpatialPolygonsDataFrame_to_df}}) as shown in the package vignettes.
#' @seealso \code{\link{SRE-class}} for details on the SRE object internals, \code{\link{auto_basis}} for automatically constructing basis functions, and \code{\link{auto_BAUs}} for automatically constructing BAUs. See also the paper \url{https://arxiv.org/abs/1705.08105} for details on code operation.
#' @keywords spatial
#' @export
#' @examples
#' library(sp)
#'
#' ### Generate process and data
#' n <- 100
#' sim_process <- data.frame(x = seq(0.005,0.995,length=n))
#' sim_process$y <- 0
#' sim_process$proc <- sin(sim_process$x*10) + 0.3*rnorm(n)
#'
#' sim_data <- sim_process[sample(1:n,50),]
#' sim_data$z <- sim_data$proc + 0.1*rnorm(50)
#' sim_data$std <- 0.1
#' coordinates(sim_data) = ~x + y # change into an sp object
#' grid_BAUs <- auto_BAUs(manifold=real_line(),data=sim_data,
#'                        nonconvex_hull=FALSE,cellsize = c(0.01),type="grid")
#' grid_BAUs$fs = 1
#'
#' ### Set up SRE model
#' G <- auto_basis(manifold = real_line(),
#'                 data=sim_data,
#'                 nres = 2,
#'                 regular = 6,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#' f <- z ~ 1
#' S <- SRE(f,list(sim_data),G,
#'          grid_BAUs,
#'          est_error = FALSE)
#'
#' ### Fit with 5 EM iterations so as not to take too much time
#' S <- SRE.fit(S,n_EM = 5,tol = 0.01,print_lik=TRUE)
#'
#' ### Check fit info
#'
#'
#' ### Predict over BAUs
#' grid_BAUs <- predict(S)
#'
#' ### Plot
#' \dontrun{
#' library(ggplot2)
#' X <- slot(grid_BAUs,"data")
#' X <- subset(X, x >= 0 & x <= 1)
#'  g1 <- LinePlotTheme() +
#'     geom_line(data=X,aes(x,y=mu)) +
#'     geom_errorbar(data=X,aes(x=x,ymax = mu + 2*sqrt(var), ymin= mu - 2*sqrt(var))) +
#'     geom_point(data = data.frame(sim_data),aes(x=x,y=z),size=3) +
#'     geom_line(data=sim_process,aes(x=x,y=proc),col="red")
#'  print(g1)}
SRE <- function(f, data,basis,BAUs, est_error = TRUE, average_in_BAU = TRUE,
                sum_variables = NULL,
                normalise_wts = TRUE,
                fs_model = "ind", vgm_model = NULL, 
                K_type = c("block-exponential", "neighbour", "unstructured", "separable"), 
                normalise_basis = TRUE, 
                response = c("gaussian", "poisson", "bernoulli", "gamma",
                             "inverse-gaussian", "negative-binomial", "binomial"), 
                link = c("identity", "log", "square-root", "logit", "probit", "cloglog", "inverse", "inverse-squared"), 
                taper = 4, include_fs = TRUE, fs_by_spatial_BAU = FALSE,
                ...) {
  
  # stop()
    
    ## Strings that must be lower-case (this allows users to enter 
    ## response = "Gaussian", for example, without causing issues)
    response  <- tolower(response)
    link      <- tolower(link)
    K_type    <- tolower(K_type)
    
    ## Allow partial matching of string arguments
    K_type   <- match.arg(K_type)
    response <- match.arg(response)
    link     <- match.arg(link)
    
    ## The weights of the BAUs only really matter if the data are SpatialPolygons. 
    ## However, we still need a 1 for all other kinds, so we still want to set them to 1.
    ## only produce a warning if we have areal data.
    if (is.null(BAUs$wts)) {
        BAUs$wts <- 1
        if (any(sapply(data, function(x) is(x, "SpatialPolygons"))) &&
            !response %in% c("binomial", "negative-binomial")) # wts doesn't come into play for binomial or neg. binomial data, as it is forced to 1
            cat("SpatialPolygons were provided for the data support. No 'wts' field was found in the BAUs, so all BAUs are assumed to be of equal weight; if this is not the case, set the 'wts' field in the BAUs accordingly.\n")
    }
    
    ## When the response has a size parameter, restrict the incidence matrices 
    ## (Cz and Cp) to represent simple sums only. This behaviour is kept vague
    ## in the paper, but perhaps a warning should be thrown to notify users.
    ## Also enforce average_in_BAU = TRUE for simplicity.
    if (response %in% c("binomial", "negative-binomial")) {
      normalise_wts <- FALSE # Set to FALSE so that Cz represents an aggregation of the mean
      BAUs$wts <- 1
      average_in_BAU <- TRUE
    }
      
    ## Check that the arguments are OK
    .check_args1(f = f, data = data, basis = basis, BAUs = BAUs, est_error = est_error, 
                 response = response, link = link, taper = taper, K_type = K_type, 
                 fs_by_spatial_BAU = fs_by_spatial_BAU, normalise_wts = normalise_wts, 
                 sum_variables = sum_variables, average_in_BAU = average_in_BAU) 
    
    ## Extract the dependent variable from the formula
    av_var <- all.vars(f)[1]

    ## Number of data objects
    ndata <- length(data)

    ## Initialise list of matrices (We construct one for every data object then concatenate)
    S <- Ve <- Vfs <- X <- Z <- Cmat <- k_Z <- list()

    ## Number of spatial BAUs and basis functions
    ns <- dim(BAUs)[1] 
    n_basis_spatial <- if(is(basis,"TensorP_Basis")) nbasis(basis@Basis1) else nbasis(basis)
    
    ## Evaluate the basis functions over the BAUs. If we have fewer spatial BAUs 
    ## than basis functions, then we average the basis functions over the BAUs
    ## using Monte Carlo integration with 1000 samples per BAU. 
    ## Otherwise, evaluate the basis functions over the BAU centroids.
    S0 <- eval_basis(basis, if(ns < n_basis_spatial) BAUs else .polygons_to_points(BAUs))
    ## ^^THIS IS WHERE THE ERROR OCCURS^^
    
    ## Normalise basis functions for the prior process to have constant variance. This was seen to pay dividends in
    ## latticekrig, however we only do it once initially
    if(normalise_basis) {
        cat("Normalising basis function evaluations at BAU level ...\n")
        xx <- sqrt(rowSums((S0) * S0))                        # Find the standard deviation (assuming unit basis function weight)
        xx <- xx + 1*(xx == 0)                                # In the rare case all basis functions evaluate to zero don't do anything
        S0 <- S0 / (as.numeric(xx))                           # Normalise the S matrix
    }

    ## Find the distance matrix associated with the basis-function centroids
    D_basis <- BuildD(basis)

    ## For each data object
    for(i in 1:ndata) {

        ## If we are estimating measurement error
        ## (FIXME: I think that this measurement error term will be very dodgy when the
        ## link is not the identity function; perhaps we should add a transformation to
        ## the data.)
        if(est_error && response == "gaussian") {
            ## Algorithm for estimating measurement error in space-time objects still not implemented
            if(is(data[[i]],"ST"))
                stop("Estimation of error not yet implemented for spatio-temporal data")
            data[[i]]$std <- 0            # Set it to zero initially
            this_data <- data[[i]]        # Allocate current data object

            ## Now estimate the measurement error using variogram methods
            this_data <- .est_obs_error(this_data,variogram.formula=f,
                                        vgm_model = vgm_model)

            ## Allocate the measurement standard deviation from variogram analysis
            data[[i]]$std <- this_data$std
        } else if (est_error && response != "gaussian") 
            cat("Not estimating the measurement error variance since the response model is not Gaussian.\n")
      
        ## The next step is to allocate all data (both point and polygon referenced) to BAUs. 
        ## We can either average data points falling in the same 
        ## BAU (average_in_BAU == TRUE) or not (average_in_BAU == FALSE).
        cat("Binning data ...\n")
  
        ## sum_variables is a vector of variable names which are to be summed
        ## rather than averaged. Typically, we will wish to sum the data/size parameter 
        ## in a setting with a size parameter (binomial or negative-binomial). 
        if (response %in% c("binomial", "negative-binomial")) 
            sum_variables <- c(all.vars(f)[1], "k_Z", sum_variables)
        
        data_proc <- map_data_to_BAUs(data[[i]],       
                                      BAUs,           
                                      average_in_BAU = average_in_BAU, 
                                      sum_variables = sum_variables)  
        
        ## The mapping can fail if not all data are covered by BAUs. Throw an error message if this is the case
        if(any(is.na(data_proc@data[av_var])))
            stop("NAs found when mapping data to BAUs. Do you have NAs in your data?
                 If not, are you sure all your data are covered by BAUs?")

        ## Extract information from the data using the .extract.from.formula internal function
        L <- .extract.from.formula(f,data=data_proc)
        X[[i]] <- as(L$X,"Matrix")                # covariate information
        Z[[i]] <- Matrix(L$y)                     # data values
        Ve[[i]] <- Diagonal(x=data_proc$std^2)    # measurement-error variance
        
        ## Construct the incidence matrix mapping data to BAUs. This just returns 
        ## indices and values which then need to be assembled into a sparse matrix.
        C_idx <- BuildC(data_proc,BAUs)

        ## Construct the sparse incidence Matrix from above indices. This is the 
        ## matrix C_Z in the vignette
        Cmat[[i]] <- sparseMatrix(i=C_idx$i_idx,
                                  j=C_idx$j_idx,
                                  x=C_idx$x_idx,   
                                  dims=c(length(data_proc),  # ensure dimensions of C are good
                                         length(BAUs)))

        ## Every data should be affected by at least one BAU. If this is not the case throw an error message.
        if(any(rowSums(Cmat[[i]])==0))
            stop("I have found difficulty in associating the data with the BAUs.
                 If you have point-referenced data
                 then this could be because you have data outside BAUs. If you have
                 polygon data, then this could be because no BAUs centroids are
                 within the polygons. For polygon data, influence on a BAU is determined from
                 whether the BAU centroid falls within the polygon or not.")

        ## We ensure the polygon observations are a weighted average over the 
        ## BAUs. This just means dividing each row by its row sum, so that each 
        ## entry is between 0 and 1, and the row sums are all equal to 1. 
        if (normalise_wts) Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]]) 
        
        ## Only the independent model is allowed for now, future implementation will include CAR/ICAR (in development)
        if(fs_model == "ind") {
          ## FIXME: Originally, Cmat was constructed to have only 1's appear in the non-zero elements,
          ## To accomodate for non-equally weighted BAUs, Cmat now has non-zero elements
          ## If we want to retain the old functionality, just use:
          # tmp <- Cmat[[i]]; tmp@x <- rep(1, length(tmp@x))
          # Vfs[[i]] <- tcrossprod(tmp %*% Diagonal(x=sqrt(BAUs$fs)))
          Vfs[[i]] <- tcrossprod(Cmat[[i]] %*% Diagonal(x=sqrt(BAUs$fs)))
        } else stop("No other fs-model implemented yet")

        ## S0 is the matrix S in the vignette. Here S is the matrix SZ in the vignette.
        S[[i]] <- Cmat[[i]] %*% S0
        
        ## Construct k_Z the same way Z is constructed
        if(response %in% c("binomial", "negative-binomial")) 
          k_Z[[i]] <- Matrix(data_proc$k_Z) 
    }

    if(fs_model == "ind") {
        Qfs_BAUs <- Diagonal(x=1/BAUs$fs)
        Vfs_BAUs <- Diagonal(x=BAUs$fs)
    } else stop("No other fs-model implemented yet")

    ## Now concatenate the matrices obtained from all the observations together
    S <- do.call("rbind",S)
    X <- do.call("rbind",X)
    Cmat <- do.call("rbind",Cmat)
    Z <- do.call("rbind",Z)
    Ve <- do.call("bdiag",Ve)
    Vfs <- do.call("bdiag",Vfs)
    ## Some matrices evaluated at observed BAUs only:
    obsidx <- unique(as(Cmat, "dgTMatrix")@j) + 1 
    C_O <- Cmat[, obsidx, drop = FALSE] 
    S_O <- S0[obsidx, , drop = FALSE]
    S_O <- drop0(S_O) # For some reason, S0 results in explicit zeros when nres = 1.
    X_BAU <- as(.extract_BAU_X_matrix(f, BAUs), "matrix") # fixed-effect design matrix at BAU level
    X_O <- X_BAU[obsidx, , drop = FALSE]
    
    # Size parameter
    if(response %in% c("binomial", "negative-binomial")) {
      k_Z <- as.numeric(do.call("rbind", k_Z))
      
      ## Size parameter associated with observed BAUs.
      ## If any observations are associated with multiple BAUs, 
      ## we require the size parameter in a field of the BAUs;
      ## otherwise, we just use the observation size parameters.
      num_BAUs_each_data_support <- table(as(Cmat, "dgTMatrix")@i)
      if (!all(num_BAUs_each_data_support == 1)) {
        if (!("k_BAU" %in% names(BAUs))) stop("When dealing with binomial or negative-binomial data, and some data supports are associated with multiple BAUs (e.g., areal data), the size parameter must be provided in the BAUs objects, in a field named 'k_BAU'.") 
        k_BAU_O <- BAUs$k_BAU[obsidx]
      } else {
        ## Note that we have to re-order the observation size parameters, so that element i
        ## of k_Z is associated with the same BAU as element i of k_BAU_O;
        ## Cmat@i contains the index of the observation associated with BAU Cmat@j
        k_BAU_O <- k_Z[as(Cmat, "dgTMatrix")@i + 1]
      }
      
      if (any(is.na(k_BAU_O))) stop("The size parameter is required at all observed BAUs")
    } 
    
    ## If we are estimating a unique fine-scale variance at each spatial BAU, 
    ## simply replicate the initialisation of sigma2fs ns times. 
    ## Also check a few things that we couldn't check before computing the matrices.
    if (fs_by_spatial_BAU) {
      ## Check each spatial BAU is observed enough times for a unique fine-scale
      ## variance parameter to be associated with each spatial BAU
      fewest_obs <-  min(table(obsidx %% ns)) # count of observations from spatial BAU with fewest observations
      if(fewest_obs == 0) {
        stop("A unique fine-scale variance at each spatial BAU can only be fit (i.e., fs_by_spatial_BAU = TRUE) if all spatial BAUs are observed, which is not the case for the provided data and BAUs. Please set fs_by_spatial_BAU = FALSE.")
      } else if(fewest_obs < 10) {
        warning(paste0("The smallest number of observations associated with a spatial BAUs is: ", fewest_obs, 
                       ". As you have selected to fit a unique fine-scale variance at each spatial BAU (i.e., fs_by_spatial_BAU = TRUE), please consider if this is a sufficient number of observations."))
      }
      
      ## Throw a warning if the the number of spatial BAUs (and hence number 
      ## of fine-scale variance parameters) is very large
      if(ns > 500) warning(paste0("The number of spatial BAUs is relatively large (",ns,"). As you have chosen to fit a separate fine-scale variance parameter at each spatial BAU, there will be ",ns," fine-scale variance parameters to estimate, which may result in difficulties in model fitting. (However, it may not be an issue.)"))
    } 

    ## If the response should be an integer, round to be safe
    ## (otherwise factorials will not make sense).
    if (response %in% c("poisson", "binomial", "negative-binomial"))
       Z <- round(Z)
  
    ## Information on fitting
    info_fit <- list("SRE not fitted to data yet. Please use SRE.fit()")

    ## Initialise the fixed effects and parameters for the EM algorithm.
    ## Note that this is done here for backwards compatability, and that 
    ## the initalisations are not computationally demanding, so it is not a major
    ## issue that we always do them irrespective of which method is used.
    ## The initialisations for method = 'TMB' is more complicated, and is 
    ## performed in SRE.fit(); doing it there means we know the method, and
    ## we do not have to create extra slots to pass the initialised parameters 
    ## from SRE() to SRE.fit().
    l <- .EM_initialise(basis, Z, X, Ve, mstar = length(obsidx))
    
    ## Dummy values:
    if(!(response %in% c("binomial", "negative-binomial"))) k_BAU_O <- k_Z <- -1
    
    ## Construct the SRE object
    new("SRE",
        data=data,
        basis=basis,
        BAUs=BAUs,
        f = f,
        S = S,
        S0 = S0,
        S_O = S_O,
        D_basis = D_basis,
        Ve = Ve,
        Vfs = Vfs,
        Vfs_BAUs = Vfs_BAUs,
        Qfs_BAUs = Qfs_BAUs,
        Z = Z,
        Cmat = Cmat,
        C_O = C_O,
        X = X,
        X_O = X_O,
        mu_eta = l$mu_eta_init,
        S_eta = l$S_eta_init,
        Q_eta = l$Q_eta_init,
        K_type = K_type,
        Khat = l$K_init,
        Khat_inv = l$K_inv_init,
        alphahat = l$alphahat_init,
        sigma2fshat = l$sigma2fshat_init,
        fs_model = fs_model,
        info_fit = info_fit, 
        response = response, 
        link = link, 
        taper = taper, 
        mu_xi = l$mu_xi_init,
        k_Z = as.numeric(k_Z), 
        k_BAU_O = as.numeric(k_BAU_O), 
        include_fs = include_fs, 
        normalise_wts = normalise_wts, 
        fs_by_spatial_BAU = fs_by_spatial_BAU)
}

## Initalise the fixed effects and parameters for method = 'EM'
.EM_initialise <- function(basis, Z, X, Ve, mstar) {
  
  l <- list() # list of initial values
  nres <- max(basis@df$res)   # Number of resolutions
  
  ## Initialise the expectations and covariances from E-step to reasonable values
  l$mu_eta_init <- Matrix(0,nbasis(basis),1)
  l$mu_xi_init <- Matrix(0,mstar,1)
  l$S_eta_init <- Diagonal(x = rep(1,nbasis(basis)))
  l$Q_eta_init <- Diagonal(x = rep(1,nbasis(basis)))
  
  ## Start with reasonable parameter estimates (that will be updated in M-step)
  l$K_init = Diagonal(n=nbasis(basis),x = 1/(1/var(Z[,1])))
  l$K_inv_init = solve(l$K_init)
  
  if(!is.finite(determinant(t(X) %*% X)$modulus))
    stop("Matrix of covariates has columns that are linearly dependent. Please change formula or covariates.")
  l$alphahat_init <- solve(t(X) %*% X) %*% t(X) %*% Z
  l$sigma2fshat_init <- mean(diag(Ve)) / 4
  
  return(l)
}


#' @rdname SRE
#' @export
SRE.fit <- function(SRE_model, n_EM = 100L, tol = 0.01, method = c("EM", "TMB"),
                    lambda = 0, print_lik = FALSE, optimiser = nlminb, 
                    known_sigma2fs = NULL, ...) {

    method <- match.arg(method)
    optimiser <- match.fun(optimiser)
    
    if (!is.null(known_sigma2fs)) 
        SRE_model@sigma2fshat <- known_sigma2fs
    
    if (method == "TMB" & SRE_model@K_type == "block-exponential") {
      tmp <- readline(cat("You have selected method = 'TMB' and K_type = 'block-exponential'. Whilst this combination is allowed, it is significantly more computationally demanding than K_type = 'neighbour'. Please enter Y if you would like to continue with the block-exponential formulation, or N if you would like to change to the more efficient neighbour based sparse precision matrix formulation."))
      if (tmp != "Y" && tmp != "N") {
        stop("You did not enter Y or N.")
      } else if (tmp == "N") {
        SRE_model@K_type <- "neighbour"
      }
    }
    

    ## Check the arguments are OK
    ## (Note that we cannot pass the SRE model, because it complicates things 
    ## for the FRK() wrapper)
    .check_args2(n_EM = n_EM, tol = tol, lambda = lambda,
                 method = method, print_lik = print_lik, 
                 fs_by_spatial_BAU = SRE_model@fs_by_spatial_BAU, 
                 response = SRE_model@response, K_type = SRE_model@K_type, link = SRE_model@link, 
                 known_sigma2fs = known_sigma2fs,
                 BAUs = SRE_model@BAUs,
                 optimiser = optimiser, ...) # control parameters to optimiser() 
    
    ## Call internal fitting function with checked arguments
    SRE_model <- .SRE.fit(SRE_model = SRE_model, n_EM = n_EM, tol = tol, 
                          method = method, lambda = lambda, print_lik = print_lik, 
                          optimiser = optimiser, known_sigma2fs = known_sigma2fs, ...)
    return(SRE_model)
}

#' @rdname SRE
#' @export
SRE.predict <- function(SRE_model, obs_fs = FALSE, newdata = NULL, pred_polys = NULL,
                        pred_time = NULL, covariances = FALSE) {
    warning("SRE.predict is deprecated. Please use predict.")
    predict(SRE_model, obs_fs = obs_fs, newdata = newdata,
            pred_polys = pred_polys, pred_time = pred_time,
            covariances = covariances)

}

#' @rdname SRE
#' @export
setMethod("predict", signature="SRE", function(object, newdata = NULL, obs_fs = FALSE, pred_polys = NULL,
                                              pred_time = NULL, covariances = FALSE, 
                                              n_MC = 400, type = "mean", k = NULL, 
                                              percentiles = c(5, 95), 
                                              kriging = "simple") {


    SRE_model <- object
    ## Deprecation coercion
    if(!is.null(pred_polys))
        newdata <- pred_polys
    
    ## Need to add prediction and uncertainty at each location, and so newdata 
    ## must be a Spatial*DataFrame (not just a Spatial* object).
    newdata <- .Coerce_SpatialDataFrame(newdata)
    
    ## The user can either provide k in SRE_model@BAUs$k_BAU, or in the predict call.
    ## This is so that the user can change k without having to call SRE() and SRE.fit() again.
    ## The k supplied in predict() will take precedence over the k stored in SRE_model@BAUs$k.
    if (SRE_model@response %in% c("binomial", "negative-binomial")) 
        if (is.null(k)) {
            if(is.null(SRE_model@BAUs$k_BAU)) {
                k <- rep(1, length(SRE_model@BAUs))
                warning("k not provided for prediction: assuming k is equal to 1 for all prediction locations.")
            } else {
                k <- SRE_model@BAUs$k_BAU
            }
        }

    
    ## Check the arguments are OK
    .check_args3(obs_fs = obs_fs, newdata = newdata, pred_polys = pred_polys,
                 pred_time = pred_time, covariances = covariances, 
                 SRE_model = SRE_model, type = type, kriging = kriging,
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
    if(is.null(pred_time) & is(SRE_model@BAUs,"ST"))
        pred_time <- 1:length(SRE_model@BAUs@time)
    
    ## Construct the CP matrix (Polygon prediction matrix)
    CP <- .make_CP(newdata = newdata, SRE_model = SRE_model)
    
    ## We start by assuming that we will predict at BAUs
    predict_BAUs <- TRUE
    
    ## If even one polygon encompasses more than one BAU, then we need to
    ## predict over linear combinations of BAUs, and hence need to
    ## compute the full covariance matrix. Note this by setting
    ## predict_BAUs to be FALSE
    if (!is.null(newdata)) {
        tmp <- as(CP, "dgTMatrix")
        if (!all(table(tmp@i) == 1))
            predict_BAUs <- FALSE
    }

    ## Call internal prediction functions depending on which method is used
    if (SRE_model@method == "EM") {
        pred_locs <- .SRE.predict(Sm = SRE_model,              # Fitted SRE model
                                  obs_fs = obs_fs,             # Case 1 or Case 2?
                                  newdata = newdata,           # Prediction polygons
                                  pred_time = pred_time,       # Prediction time points
                                  covariances = covariances,   # Compute covariances?
                                  CP = CP,                     # Polygon prediction matrix
                                  predict_BAUs = predict_BAUs) # Are we predicting at BAUs?                   
        
    } else if (SRE_model@method == "TMB") {
        pred_locs <- .FRKTMB_pred(M = SRE_model,               # Fitted SRE model
                                  newdata = newdata,           # Prediction polygons
                                  CP = CP,                     # Polygon prediction matrix
                                  predict_BAUs = predict_BAUs, # Are we predicting at BAUs?
                                  pred_time = pred_time,       # Prediction time points
                                  n_MC = n_MC,                 # Desired number of MC simulations
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
      if(SRE_model@method == "TMB") {
        
        ## If method = "TMB", we have a list of MC samples which must
        ## also be altered. To do this, we select BAUs based on identifiers.
        ## We may not necessarily have identifiers, particularly if 
        ## the BAUs were not created with auto_BAUs(). Create some:
        BAU_UID <- .UIDs(S@BAUs)
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
        
      } else if (SRE_model@method == "EM") {
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


## unexported function to construct CP (the prediction polygon matrix)
## I made this into a function because it is common for the methods ("TMB" and
## "EM"), so it should be done before splitting. 
.make_CP <- function (newdata = NULL, SRE_model) {
    
    ## If the user has not specified polygons over which to predict, then CP is
    ## just the diagonal matrix and we predict over all the BAUs
    if(is.null(newdata)) {
        CP <- Diagonal(length(SRE_model@BAUs))
    } else {
        ## The user has specified arbitrary polygons
        ## Based on these polygons construct the C matrix
        newdata2 <- map_data_to_BAUs(newdata, SRE_model@BAUs,
                                     average_in_BAU = FALSE,
                                     sum_variables = NULL,
                                     est_error = FALSE)
        C_idx <- BuildC(newdata2, SRE_model@BAUs)
        
        CP <- sparseMatrix(i = C_idx$i_idx,
                           j = C_idx$j_idx,
                           x = C_idx$x_idx,
                           dims = c(length(newdata),
                                    length(SRE_model@BAUs)))
        
        ## As in SRE(), make sure the polygons are averages (not sums) if requested
        if (SRE_model@normalise_wts)
            CP <- CP / rowSums(CP)
    }
    return(CP)
}


#' @rdname SRE
#' @export
loglik <- function(SRE_model) {
    # This is structured this way so that extra models for fs-variation
    # can be implemented later
    if (length(SRE_model@log_likelihood) != 0) {
        return(SRE_model@log_likelihood)
    } else if (SRE_model@fs_model == "ind") {
        return(.loglik.ind(SRE_model))
    } else {
        stop("Currently only independent fine-scale model is implemented")
    }
}

## Print/Show SRE
print.SRE <- function(x,...) {
    print_list <- list(
        formula =  deparse(x@f),
        ndatasets = length(x@data),
        nbasis = x@basis@n,
        basis_class = class(x@basis)[1],
        nBAUs = length(x@BAUs),
        nobs = length(x@Z),
        mean_obsvar = mean(x@Ve@x), # FIXME: @
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
#' @aliases info_fit,SRE-method
setMethod("info_fit", signature(SRE_model = "SRE"),
          function(SRE_model) {SRE_model@info_fit})

# Retrieve coefficients of SRE model
#' @rdname coef
#' @aliases coef, SRE-method
setMethod("coef",signature(object = "SRE"),function(object,...) {
    coeff <- as.numeric(object@alphahat)
    varnames <- all.vars(object@f)[-1]
    nms <- "Intercept"
    if(length(varnames) > 0) {
        nms <- c(nms, varnames)
    }
    names(coeff) <- nms
    coeff
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



#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="SpatialPixelsDataFrame"),function(BAUs) {
    
    ## Number of dimensions of the BAUs
    n_coord <- dimensions(BAUs)
    
    ## Reverse the slots (must do this for each slot):
    ## First, deal with the @data slot.
    ## It is possible that there is additional data columns present, 
    ## so we need to cater for this.
    ## Desired coordinate order and their indices in the data:
    new_coord_order <- rev(coordnames(BAUs))
    coord_idx <- match(new_coord_order, names(BAUs@data), nomatch = 0)
    
    ## Indices of the remaining columns
    remaining_idx <- which(!names(BAUs@data) %in% new_coord_order)
    
    ## Subset the data
    BAUs@data <- BAUs@data[, c(coord_idx, remaining_idx)]
    
    ## Reverse the remaining slots
    BAUs@coords <- BAUs@coords[, new_coord_order]
    BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
    
    ## Note that we cannot subset BAUs@grid (S4 method). 
    ## Instead, we will manipulate the slots directly.
    tmp <- BAUs@grid
    tmp@cellcentre.offset <- tmp@cellcentre.offset[n_coord:1]
    tmp@cellsize <- tmp@cellsize[n_coord:1]
    tmp@cells.dim <- tmp@cells.dim[n_coord:1]
    ## NB: note sure how to access the coordinate names of tmp directly as of now (GridTopology
    ## objects have only an assignment method for coordnames, not a retrieval method). 
    ## For now I will just use the coordinate names of the BAUs.
    ## It doesn't seem to matter anyway, though, as it adjusts automatically.
    coordnames(tmp) <- coordnames(BAUs)[n_coord:1]
    
    BAUs@grid <- tmp
    
    return(BAUs)
})


#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="SpatialPointsDataFrame"),function(BAUs) {
    
    ## The documentation says that the @data slot may in fact contain the coordinates.
    ## So, we need to allow for this possibility.
    ## Desired coordinate order and their indices in the data:
    new_coord_order <- rev(coordnames(BAUs))
    coord_idx <- match(new_coord_order, names(BAUs@data), nomatch = 0)
    
    ## Indices of the remaining columns
    remaining_idx <- which(!names(BAUs@data) %in% new_coord_order)
    
    ## Subset the data
    BAUs@data <- BAUs@data[, c(coord_idx, remaining_idx), drop = FALSE]
    
    ## Reverse the order of the remaining slots:
    BAUs@coords <- BAUs@coords[, new_coord_order]
    BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
    
    return(BAUs)
})


#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="SpatialPoints"),function(BAUs) {
    
    ## The documentation does not explicitly say that coords can or cannot contain data, 
    ## but for safety I will assume it can. 
    ## For the spatial points object, I will simply reverse the coordinates matrix. 
    ## Hence, if we have (lon, lat, z), the reversed order will be (z, lat, lon).
    ## Desired coordinate order and their indices in the data:
    new_coord_order <- rev(coordnames(BAUs))
    
    ## Subset the data
    BAUs@coords <- BAUs@coords[, new_coord_order, drop = FALSE]
    
    ## Reverse the order of the remaining slot:
    BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
    
    return(BAUs)
})


#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="SpatialPolygonsDataFrame"),function(BAUs) {
    
    
    ## First, the @data slot:
    ## (Note that both coordinate and data columns may be present - we will just reverse all)
    new_coord_order <- 
    BAUs@data <- BAUs@data[, rev(colnames(BAUs@data))]
    
    ## Second, the @bbox:
    BAUs@bbox   <- BAUs@bbox[rev(row.names(BAUs@bbox)), ]
    
    ## Third, the @polygons
    ## FIXME: Have to do this for each polygons list (use lapply or something)
    tmp <- BAUs@polygons$`1` 
    class(BAUs@polygons$`1`)
    ## Swap the order of the label point ($labpt) slot:
    tmp@labpt <- rev(tmp@labpt)
    ## Swap the order of the coordinates ($coords) slot:
    coordinates(tmp@Polygons)
    coordinates(tmp)
    ## FIXME: not sure how to access the coordinates of the polygons.
    ## Don't want to spend any more time on this until Andrew says its worthwhile.
    
    return(BAUs)
})



#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="STFDF"),function(BAUs) {
    
    ## FIXME: should allow for SpatialPolygonsDF and SpatialPolygons too
    if(is(BAUs@sp, "SpatialPointsDataFrame") | 
       is(BAUs@sp, "SpatialPixelsDataFrame") | 
       is(BAUs@sp, "SpatialPoints"))
        BAUs@sp <- reverse_spatial_coords(BAUs@sp)
    else
        stop("The underlying spatial object should be of class 'SpatialPointsDataFrame' or 'SpatialPixelsDataFrame'.")
    
    return(BAUs)
})


#' @rdname reverse_spatial_coords
#' @export
setMethod("reverse_spatial_coords",signature(BAUs="STIDF"),function(BAUs) {
    
    ## FIXME: should allow for SpatialPolygonsDF and SpatialPolygons too
    if(is(BAUs@sp, "SpatialPointsDataFrame") | 
       is(BAUs@sp, "SpatialPixelsDataFrame") | 
       is(BAUs@sp, "SpatialPoints"))
        BAUs@sp <- reverse_spatial_coords(BAUs@sp)
    else
        stop("The underlying spatial object should be of class 'SpatialPointsDataFrame', 'SpatialPixelsDataFrame', or 'SpatialPoints'.")
    
    return(BAUs)
})


## Code common to the removal of spatial BAUs across methods
.remove_spatial_BAUs <- function(BAUs, rmidx, redefine_index = FALSE) {
    ntot <- nrow(BAUs@coords)
    if(!all(rmidx %in% 1:ntot))
        stop("Please ensure indices are numeric and within
             1 and the number of spatial BAUs.")
    
    BAUs_orig  <- BAUs
    BAUs <- BAUs[-rmidx, ]
    
    ## Check for indexing columns in @data slot of BAU object
    ## Note that the first check is for SpatialPoints and SpatialPolygons, each of which 
    ## do not have a @data slot.

    if ("data" %in% slotNames(class(BAUs_orig)) & redefine_index) {
        BAUs@data <- .redefine_indexing_variables(BAUs_orig@data, BAUs@data)
    }
    
    
    return(BAUs)
}


## Function to check if columns are "indexing variables".
## Returns the indices of these variables in the column names of df.
.redefine_indexing_variables <- function(data_orig, data_new) {
    
    ## Find the column indices of the indexing variables (if present)
    indexing_variables_idx <- apply(data_orig, 2, function(x) all(x == 1:nrow(data_orig))) %>%
        which()
    
    
    ## if indexing variables are present, redefine the corresponding column.
    if (length(indexing_variables_idx)) { 
        print(paste0("Possible indexing variables present, which will be redefined to maintain indexing following the removal of BAUs: ", names(indexing_variables_idx)))
        data_new[, indexing_variables_idx] <- 1:nrow(data_new)
    } 
    
    return(data_new)
}

#' @rdname remove_BAUs
#' @export
setMethod("remove_BAUs",signature(BAUs="SpatialPoints"),function(BAUs, rmidx, redefine_index = FALSE) {
    return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})


#' @rdname remove_BAUs
#' @export
setMethod("remove_BAUs",signature(BAUs="SpatialPointsDataFrame"),function(BAUs, rmidx, redefine_index = FALSE) {
    return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})


#' @rdname remove_BAUs
#' @export
setMethod("remove_BAUs",signature(BAUs="SpatialPixelsDataFrame"),function(BAUs, rmidx, redefine_index = FALSE) {
    return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})

#' @rdname remove_BAUs
#' @export
setMethod("remove_BAUs",signature(BAUs="STFDF"),function(BAUs, rmidx, redefine_index = FALSE) {
    
    BAUs_orig <- BAUs # Backup of BAUs for checks later
    
    ## Remove the spatial BAUs:
    BAUs@sp <- .remove_spatial_BAUs(BAUs@sp, rmidx)
    
    ## SUBSET @data to adjust for the new spatial BAUs:
    ## From the documentation, it says that @data is a data.frame containing
    ## measured values; space index cycling first, and time order preserved. 
    ## Hence, we should be able to remove the rows which correspond to rmidx
    ## by adding time to rmidx. 
    ## Also need to maintain space index cycling first. Achieve this by sorting the indcies.
    n_time     <- length(BAUs@time) # number of temporal frames
    n_spat     <- nrow(BAUs_orig@sp)      # ORIGINAL number of spatial BAUs
    increment   <- n_spat * 0:(n_time - 1) # Terms to add to each rmidx
    data_rmidx<- sort(c(outer(rmidx, increment, "+")))
    
    BAUs@data <- BAUs@data[-data_rmidx, ]
    
    ## Check for indexing columns in @data slot of BAU object
    ## Note that the first check is for SpatialPoints and SpatialPolygons, each of which 
    ## do not have a @data slot.
    if ("data" %in% slotNames(class(BAUs_orig)) & redefine_index)
        BAUs@data <- .redefine_indexing_variables(BAUs_orig@data, BAUs@data)
    
    return(BAUs)
})



#' @rdname remove_BAUs
#' @export
setMethod("remove_BAUs",signature(BAUs="STIDF"),function(BAUs, rmidx, redefine_index = FALSE) {
    
    BAUs_orig <- BAUs # Backup of BAUs for checks later
    
    ## Removal is straightforward for STIDF
    suppressWarnings(BAUs <- BAUs[-rmidx, ])
    
    ## Check for indexing columns in @data slot of BAU object
    ## Note that the first check is for SpatialPoints and SpatialPolygons, each of which 
    ## do not have a @data slot.
    if ("data" %in% slotNames(class(BAUs_orig)) & redefine_index)
        BAUs@data <- .redefine_indexing_variables(BAUs_orig@data, BAUs@data)
    
    return(BAUs)
})


#' @rdname observed_BAUs
#' @export
setMethod("observed_BAUs", signature(SRE_model = "SRE"), function (SRE_model) {
    
    ## Note that Cmat maps BAUs to the observations. The dimension of SRE_model@Cmat is
    ## (number of observations) * (number of BAUs).
    Cmat <- as(SRE_model@Cmat, "dgTMatrix")
    obsidx <- unique(Cmat@j) + 1 # unique() probably not necessary but want to be safe
    
    return(obsidx)
})


#' @rdname unobserved_BAUs
#' @export
setMethod("unobserved_BAUs",signature(SRE_model = "SRE"), function (SRE_model) {
    
    ## Id of observed BAUs:
    obsidx <- observed_BAUs(SRE_model)
    
    ## Id of unobserved BAUs (ncol(SRE_model@Cmat) is the total number of BAUs):
    unobsidx <- (1:ncol(SRE_model@Cmat))[-obsidx]
    
    return(unobsidx)
})


##################################
#### NOT EXPORTED ################
##################################

## Main prediction routine
.SRE.fit <- function(SRE_model, n_EM, tol, method, lambda, 
                     print_lik, optimiser = nlminb, known_sigma2fs, ...) {

    info_fit <- list()      # initialise info_fit
    
    if(method == "EM") {
        
        n <- nbasis(SRE_model)  # number of basis functions
        X <- SRE_model@X        # covariates

        info_fit$method <- "EM" # updated info_fit
        llk <- rep(0,n_EM)       # log-likelihood

        ## If user wishes to show progress show progress bar
        if(opts_FRK$get("progress"))
            pb <- utils::txtProgressBar(min = 0, max = n_EM, style = 3)

        ## For each EM iteration step
        for(i in 1:n_EM) {
            llk[i] <- loglik(SRE_model)                          # compute the log-lik
            SRE_model <- .SRE.Estep(SRE_model)                   # compute E-step
            SRE_model <- .SRE.Mstep(SRE_model, lambda = lambda)  # compute M-step
            if(opts_FRK$get("progress"))
                utils::setTxtProgressBar(pb, i)                  # update progress bar
            if(i>1)                                              # If we're not on first iteration
                if(abs(llk[i] - llk[i-1]) < tol) {                 # Compute change in log-lik
                    cat("Minimum tolerance reached\n")           # and stop if less than tol
                    break
                }
        }

        if(opts_FRK$get("progress")) close(pb)           # close progress bar
        info_fit$num_iterations <- i                     # update fit info


        ## If zero fine-scale variation detected just make sure user knows.
        ## This can be symptomatic of poor fitting
        if(SRE_model@sigma2fshat == 0) {
            info_fit$sigma2fshat_equal_0 <- 1
            if(opts_FRK$get("verbose") > 0)
                message("sigma2fs is being estimated to zero.
                        This might because of an incorrect binning
                        procedure or because too much measurement error
                        is being assumed (or because the latent
                        field is indeed that smooth, but unlikely).")
        } else {
            info_fit$sigma2fshat_equal_0 <- 0
        }

        ## If we have reached max. iterations, tell the user
        if(i == n_EM) {
            cat("Maximum EM iterations reached\n")
            info_fit$converged <- 0    # update info_fit
        } else {
            info_fit$converged <- 1   # update info_fit
        }

        ## Plot log-lik vs EM iteration plot
        info_fit$plot_lik <- list(x = 1:i, llk = llk[1:i],
                                  ylab = "log likelihood",
                                  xlab = "EM iteration")

        ## If user wants to see the log-lik vs EM iteration plot, plot it
        if(print_lik & !is.na(tol)) {
            plot(1:i, llk[1:i],
                 ylab = "log likelihood",
                 xlab = "EM iteration")

        }
    } else if (method == "TMB") {
        SRE_model <- .FRKTMB_fit(SRE_model, optimiser = optimiser, known_sigma2fs = known_sigma2fs, ...)
    } else {
        stop("No other estimation method implemented yet. Please use method = 'EM' or method = 'TMB'.")
    }

    ## Return fitted SRE model
    SRE_model@info_fit <- info_fit
    SRE_model@method <- method
    return(SRE_model)
    }

## E-Step
.SRE.Estep <- function(Sm) {
    # This is structured this way so that extra models for fs-variation
    # can be implemented later
    if(Sm@fs_model == "ind")
        Sm <- .SRE.Estep.ind(Sm)
    else stop("E-step only for independent fs-variation model currently implemented")
}

## Compute log-likelihood for independent fs-variation model
.loglik.ind <- function(Sm) {

    S <- Sm@S                               # basis-function matrix
    K <- Sm@Khat                            # random-effects cov. matrix
    chol_K <- chol(K)                       # its Cholesky
    Kinv <- chol2inv(chol_K)                # random-effects prec. matrix
    resid <- Sm@Z - Sm@X %*% Sm@alphahat    # residuals at fitted estimates
    N <- length(Sm@Z)                       # number of data points

    D <- Sm@sigma2fshat*Sm@Vfs + Sm@Ve      # total variance of data
    if(isDiagonal(D)) {                     # if this is diagonal
        D <- Diagonal(x = D@x)              # cast to diagonal matrix
        cholD <- sqrt(D)                    # just compute sqrt
        cholDinvT <- solve(cholD)           # and the inverse is just the reciprocal
    } else {
        cholD <- chol(D)                    # otherwise do the Cholesky
        cholDinvT <- t(solve(cholD))        # find the transposed (lower) inverse of the factor
    }

    ## Compute log-determinant. This is given by a formula in Section 2.2
    S_Dinv_S <-  crossprod(cholDinvT %*% S)
    R <- chol(Kinv + S_Dinv_S)
    log_det_SigmaZ <- logdet(R) +
        determinant(K,logarithm = TRUE)$modulus +
        logdet(cholD)  # this computes the log-determinant of a matrix from its Cholesky factor

    ## Alternatively: (slower but more direct)
    # Dinv <- chol2inv(chol(D))
    # SigmaZ_inv <- Dinv - Dinv %*% S %*% solve(Kinv + S_Dinv_S) %*% t(S) %*% Dinv
    # SigmaZ_inv2 <- Dinv - tcrossprod(Dinv %*% S %*% solve(R))

    ## Compute efficiently rDinv <- t(resid) %*% Dinv
    rDinv <- crossprod(cholDinvT %*% resid,cholDinvT)

    ## Compute the quadratic portion of the log-lik
    ## This is the same as quad_bit <- rDinv %*% resid - tcrossprod(rDinv %*% S %*% solve(R))
    ## but more efficient
    quad_bit <- crossprod(cholDinvT %*% resid) - tcrossprod(rDinv %*% S %*% solve(R))

    ## Now just add the bits together
    llik <- -0.5 * N * log(2*pi) -
        0.5 * log_det_SigmaZ -
        0.5 * quad_bit

    as.numeric(llik) # convert to numeric and return

}

## M-step
.SRE.Mstep <- function(Sm, lambda = 0) {
    # This is structured this way so that extra models for fs-variation
    # can be implemented later
    if(Sm@fs_model == "ind")
        Sm <- .SRE.Mstep.ind(Sm, lambda = lambda)
    else stop("M-step only for independent fs-variation model currently implemented")
}

## E-step for independent fs-variation model
.SRE.Estep.ind <- function(Sm) {
    alpha <- Sm@alphahat           # current regression coefficients estimates
    K <- Sm@Khat                   # current random effects covariance matrix estimate
    Kinv <- Sm@Khat_inv            # current random effects precision matrix estimate
    sigma2fs <- Sm@sigma2fshat     # current fs-variation factor estimate

    D <- sigma2fs*Sm@Vfs + Sm@Ve   # total variance-covariance of Z
    if(isDiagonal(D)) {            # if this is diagonal
        D <- Diagonal(x=D@x)       # cast to diagonal
        cholD <- sqrt(D)           # then the Cholesky is the sqrt
        cholDinv <- solve(cholD)   # the inverse Cholesky is the inverse-sqrt
        Dinv <- solve(D)           # the inverse is just reciprocal of diagonal elements
    } else {
        cholD <- Matrix::chol(D)   # if not diagonal then do Cholesky
        cholDinv <- solve(cholD)   # invery Cholesky factor
        Dinv <- chol2inv(cholD)    # find the inverse from the Cholesky factor
    }

    ## The below are simple Gaussian updating equations in Section 2.2 of the vignette
    Q_eta <- (crossprod(t(cholDinv) %*% Sm@S) + Kinv)
    S_eta <- chol2inv(chol(Q_eta))  # we can invert since we are low rank in FRK
    mu_eta <- (S_eta) %*%(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))


    ## Deprecated:
    # if(!is(Q_eta,"dsCMatrix")) Q_eta <- as(Q_eta,"dsCMatrix")
    # chol_Q_eta <- cholPermute(Q_eta)
    # mu_eta <- cholsolve(Q_eta,(t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha)),
    #                     perm=TRUE, cholQp = chol_Q_eta$Qpermchol,P = chol_Q_eta$P)
    # S_eta <- Matrix()

    ## Deprecated:
    # S_eta <- chol2inv(chol(crossprod(t(cholDinv) %*% Sm@S) + Kinv))
    # mu_eta <- S_eta %*% (t(Sm@S) %*% Dinv %*% (Sm@Z - Sm@X %*% alpha))

    Sm@mu_eta <- mu_eta  # update conditional mean
    Sm@S_eta <- S_eta    # update conditional covariance
    Sm@Q_eta <- Q_eta    # update conditional precision
    Sm                   # return SRE object
}

## M-step for the indepdent fine-scale variation model
.SRE.Mstep.ind <- function(Sm, lambda = 0) {

    mu_eta <- Sm@mu_eta              # current cond. mean of random effects
    S_eta <- Sm@S_eta                # current cond. cov. matrix of random effects
    alpha <- Sm@alphahat             # regression coefficients
    sigma2fs <- Sm@sigma2fshat       # fine-scale variance

    K <- .update_K(Sm,method=Sm@K_type,  # update the prior covariance matrix K
                   lambda = lambda)
    Khat_inv <- chol2inv(chol(K))        # compute the precision

    ## If the measurement and fs. variational covariance matricies
    ## are proportional to the identity then we have the
    ## special case of homoscedasticity
    if(all((a <- diag(Sm@Ve)) == a[1]) &
       all((b <- diag(Sm@Vfs)) == b[1]) &
       isDiagonal(Sm@Vfs))    {
        homoscedastic <- TRUE
    } else {
        homoscedastic <- FALSE
    }

    ## If the measurement and fs. variational covariance matricies
    ## are diagonal then we have another special case
    if(isDiagonal(Sm@Ve) & isDiagonal(Sm@Vfs))    {
        diagonal_mats <- TRUE
    } else {
        diagonal_mats <- FALSE
    }

    ## If we have some fine-scale variation terms
    if(!all(diag(Sm@Vfs) == 0))
        ## And we're not in the diagonal case (this is the most comp. intensive)
        if(!diagonal_mats) {
            ## We first need to create a function whose root is sigma2fshat
            ## See 2.2 of vignette for equation details
            J <- function(sigma2fs) {
                if(sigma2fs < 0) {
                    return(Inf)                              # cannot be less than 0
                } else {
                    D <- sigma2fs*Sm@Vfs + Sm@Ve             # total data variance-covariance
                    Dinv <- chol2inv(chol(D))                # it's inverse (this will be blocked so still sparse)
                    DinvV <- Dinv %*% Sm@Vfs                 # summary matrix
                    DinvVDinv <- Dinv %*% Sm@Vfs %*% Dinv    # summary matrix

                    alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%
                        t(Sm@X) %*% Dinv %*%                 # regression coefficients GLS estimates
                        (Sm@Z - Sm@S %*% mu_eta)
                    resid <- Sm@Z - Sm@X %*% alpha           # fitted residuals

                    Dinvr <- DinvVDinv %*% resid             # summary vector
                    DinvS <- DinvVDinv %*% Sm@S              # summary vector

                    tr1 <- tr(DinvV)                         # compute trace of first term

                    ## Compute trace of second term
                    tr2 <- sum(diag2(DinvS %*% (S_eta +  tcrossprod(mu_eta)),t(Sm@S))  -
                                   2*diag2(DinvS %*% mu_eta,t(resid)) +
                                   diag2(Dinvr,t(resid)))

                    ## return value of function
                    -(-0.5*tr1 +0.5*tr2)
                }
            }
        } else {

            ## If we have diagonal matrices then some simplifications are possible
            R_eta <- chol(S_eta + tcrossprod(mu_eta)) # Cholesky factor of E(eta eta^T)
            S_R_eta <- Sm@S %*% t(R_eta)              # summary matrix
            Omega_diag1 <- rowSums(S_R_eta^2)         # first part of diag(Omega) as in vignette
            J <- function(sigma2fs) {
                if(sigma2fs < 0) {
                    return(Inf)                       # sigma2fs >= 0
                } else {
                    D <- sigma2fs*Sm@Vfs + Sm@Ve    # total data variance. This must be diagonal
                    Dinv <- solve(D)                # just take reciprocal since D is defo. diagonal here
                    DinvV <- Dinv %*% Sm@Vfs        # summary matrix (diagonal)

                    alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%
                        t(Sm@X) %*% Dinv %*%          # regression coefficients GLS estimates
                        (Sm@Z - Sm@S %*% mu_eta)
                    resid <- Sm@Z - Sm@X %*% alpha    # fitted residuals
                    Omega_diag <- Omega_diag1 -       # other parts of diag(Omega)
                        2*diag2(Sm@S %*% mu_eta, t(resid)) +
                        diag2(resid,t(resid))
                    Omega_diag <- Diagonal(x=Omega_diag)  # compute a diagonal Omega matrix

                    ## Since DinvV and Dinv are diagonal and we only want the trace,
                    ## we only need the diagonal elements of Omega in the following
                    -(-0.5*tr(DinvV) +
                          0.5*tr(DinvV %*% Dinv %*% Omega_diag)
                    )
                }
            }
        }

    ## We need to find the root in J. For this we need to start uniroot with
    ## values on either side of sigma2fshat. The below implements
    ## a simple search algorithm for finding a good starting values.

    ## If we have some fine-scale variation terms
    if(!all(diag(Sm@Vfs) == 0))
        ## And we're not in the special homoscedastic case
        if(!homoscedastic) {
            amp_factor <- 10; OK <- 0  # initialise
            while(!OK) {
                amp_factor <- amp_factor * 10 # widen the interval

                ## If the signs are different, then we're OK, otherwise not
                if(!(sign(J(sigma2fs/amp_factor)) == sign(J(sigma2fs*amp_factor)))) OK <- 1

                ## If we have a really big amp_factor, it means we're not getting anywhere and
                ## sigma2fshat is probably tending to zero.
                if(amp_factor > 1e9) {
                    OK <- 1
                }
            }

            if(amp_factor > 1e9) {
                sigma2fs_new <- 0  # fix sigma2fshat to zero since we couldn't estimate it
            } else {
                ## Otherwise find the root of the equation with the sought initial conditions
                sigma2fs_new <- stats::uniroot(f = J,
                                               interval = c(sigma2fs/amp_factor,
                                                            sigma2fs*amp_factor))$root
            }
            D <- sigma2fs_new*Sm@Vfs + Sm@Ve  # total data variance-covariance
            if(isDiagonal(D)) {               # inverse of D (as above)
                D <- Diagonal(x=D@x)          # cast to Diagonal
                Dinv <- solve(D)
            } else {
                Dinv <- chol2inv(chol(D))
            }
            alpha <- solve(t(Sm@X) %*% Dinv %*% Sm@X) %*%      # alpha GLS estimate (as above)
                t(Sm@X) %*% Dinv %*% (Sm@Z - Sm@S %*% mu_eta)
        } else {
            ## Here we are in the homoscedastic (diagonal) case and we can
            ## solve for alpha independently of sigma2fshat
            alpha <- solve(t(Sm@X) %*% Sm@X) %*%      # alpha GLS estimate
                t(Sm@X) %*% (Sm@Z - Sm@S %*% mu_eta)
            resid <- Sm@Z - Sm@X %*% alpha            # residual
            Omega_diag <- Omega_diag1 -               # just compute Omega once
                2*diag2(Sm@S %*% mu_eta, t(resid)) +
                diag2(resid,t(resid))
            Omega_diag <- Diagonal(x=Omega_diag)

            ## Closed-form solution for sigma2fs (see vignette)
            sigma2fs_new <- 1/b[1]*(sum(Omega_diag)/length(Sm@Z) - a[1])
            if(sigma2fs_new < 0) {  # If we get less than zero because of numeric instability
                sigma2fs_new = 0    # just fix to zero
            }
        }

    ## If we do NOT have any fine-scale variation (e.g., estimated to zero in previous iteration)
    if(all(diag(Sm@Vfs) == 0)) {
        alpha <- solve(t(Sm@X) %*% solve(Sm@Ve) %*% Sm@X) %*% t(Sm@X) %*%  # just find GLS
            solve(Sm@Ve) %*% (Sm@Z - Sm@S %*% mu_eta)
        sigma2fs_new <- 0                                                  # and keep sigma2fs at zero
    }

    ## Update SRE model with estimated quantities
    Sm@Khat <- K
    Sm@Khat_inv <- Khat_inv
    Sm@alphahat <- alpha
    Sm@sigma2fshat <- sigma2fs_new

    ## Return SRE model
    Sm
}

## This routines updates the covariance matrix of the random effects
.update_K <- function(Sm,method="unstructured",
                      S_eta= NULL,mu_eta = NULL,
                      lambda = 0) {

    if (is.null(S_eta)) S_eta <- Sm@S_eta      # Conditional covariance matrix of random effects
    if (is.null(mu_eta)) mu_eta <- Sm@mu_eta   # Conditional mean of random effects

    if(method == "unstructured") {
        ## If K is unstructured, then the update is trivial, see vignette Section 2.2
        ## I allow for some regularisation through lambda should this be deemed required
        ## (This is useful for when we have lots of basis and few data points)
        K <- .regularise_K(Sm, lambda = lambda)
    } else if (method == "block-exponential") {
        ## If K is block exponential (blocked dby resolution) then
        ## we need to find the (i) precision, (ii) spatial length scale, and
        ## (iii) temporal length scale by resolution
        all_res <- count_res(Sm)               # number of resolutions
        eta2 <- lapply(1:nrow(all_res),function(i) {
            ## find which indices correspond to these basis functions
            idx <- which(data.frame(Sm@basis)$res == i)
            S_eta[idx,idx] +
                tcrossprod(mu_eta[idx])
        })

        ## (i) Find the precision associated with each resolution
        omega <- lapply(1:nrow(all_res),       # for each resolution
                        function(i) {
                            ## number of basis functions in i-th resolution
                            ni <- all_res[i,]$n

                            ## find which indices correspond to these basis functions
                            idx <- which(data.frame(Sm@basis)$res == i)

                            ## # find the current CORRELATION matrix associated with this resolution
                            Ki <- Sm@Khat[idx,idx]/Sm@Khat[idx[1],idx[1]]

                            ## Compute INVERSE CORRELATION matrix associated with this resolution
                            Ki_inv <- chol2inv(chol(Ki))

                            ## The precision is given by n / tr(Kinv %*% (S_eta + mu.mu'))
                            ni / sum(diag2(Ki_inv,eta2[[i]]))
                        })

        ## (ii,iii) Likelihood function for spatial/temporal length scales
        f_tau <- function(tau_i,i) {   # tau_i are the scales, i is the resolution
            if(any(tau_i <= 1e-10)) {  # do not let any of the taus be too small
                Inf
            } else {

                ## Find which bases are at this resolution
                idx <- which(data.frame(Sm@basis)$res == i)

                ## Since we're block exponential, the correlation matrix is simply
                ## computed from the distances using the appropriate decay parameters
                if(is(Sm@basis,"TensorP_Basis")) {
                    ## If we have a tensor basis then construct Ki using the Kronecker product
                    Ki1 <- exp(-Sm@D_basis$Basis2[[1]]/tau_i[2])  # temporal part
                    Ki2 <- exp(-Sm@D_basis$Basis1[[i]]/tau_i[1])  # spatial part
                    ## time runs slowest (and only one time resolution),  space runs fastest
                    Ki <- kronecker(Ki1,Ki2)

                    ## Compute the inverse correlation matrix
                    Qi1 <- chol2inv(chol(Ki1))
                    Qi2 <- chol2inv(chol(Ki2))
                    Ki_inv <- kronecker(Qi1,Qi2)

                    ## Compute log determinant
                    R1 <- chol(Qi1)
                    R2 <- chol(Qi2)
                    det_part <- 0.5*(nrow(R2)*logdet(R1) + nrow(R1)*logdet(R2))

                    ## Compute the log=likelihood. There doesn't seem to be a way to
                    ## simplify this using the Kronecker product
                    -as.numeric(det_part - omega[[i]]/2*sum(diag2(Ki_inv,eta2[[i]],symm=TRUE)))

                } else {
                    ## Just spatial, from distances between centroid
                    Ki <- exp(-Sm@D_basis[[i]]/tau_i)

                    ## Compute the inverse correlation matrix
                    Ki_inv <- chol2inv(chol(Ki))

                    ## Compute the log-likelihood
                    -as.numeric(0.5*determinant(Ki_inv)$modulus -
                                    omega[[i]]/2*sum(diag2(Ki_inv,eta2[[i]],symm=TRUE)))
                }

            }
        }

        ## (ii,iii) GRADIENT of the likelihood function for spatial/temporal length scales
        gr_f_tau <- function(tau_i,i) {
            idx <- which(Sm@basis@df$res == i)  # tau_i are the scales, i is the resolution

            if(is(Sm@basis,"TensorP_Basis")) {
                Ki1 <- exp(-Sm@D_basis$Basis2[[1]]/tau_i[2])  # temporal part
                Ki2 <- exp(-Sm@D_basis$Basis1[[i]]/tau_i[1])  # spatial part
                Ki <- kronecker(Ki1,Ki2)                      # Kronecker of the two

                ## Compute the inverse correlation matrix
                Qi1 <- chol2inv(chol(Ki1))
                Qi2 <- chol2inv(chol(Ki2))
                Ki_inv <- kronecker(Qi1,Qi2)

                ## d(X kron Y) = dX kron Y + X cron dY. Compute these below
                dKi <- kronecker(Ki1,(Sm@D_basis$Basis1[[i]]/(tau_i[1]^2))*Ki2)
                dKit <- kronecker((Sm@D_basis$Basis2[[1]]/(tau_i[2]^2))*Ki1,Ki2)

            } else {
                ## If only spatial then just compute derivative of exponential
                Ki <- exp(-Sm@D_basis[[i]]/tau_i)
                dKi <- (Sm@D_basis[[i]]/(tau_i^2))*exp(-Sm@D_basis[[i]]/tau_i)
                Ki_inv <- chol2inv(chol(Ki))  # inverse
            }

            ## derivative of log-likelihodd w.r.t tau_1 (spatial)
            tau_i1 <- -(-0.5*sum(diag2(dKi,Ki_inv)) +
                            0.5*omega[[i]]*sum(diag2(eta2[[i]]%*% Ki_inv,
                                                     dKi %*% Ki_inv)))
            tau_i1 <- as.numeric(tau_i1)

            if(length(tau_i) == 1) {  # Then we just have space
                return(tau_i1)
            } else {                  # We have time aswell

                ## derivative of log-likelihood w.r.t tau_2 (temporal)
                tau_i2 <-  -(-0.5*sum(diag2(dKit,Ki_inv)) +
                                 0.5*omega[[i]]*sum(diag2(eta2[[i]]%*% Ki_inv,
                                                          dKit %*% Ki_inv)))
                tau_i2 <- as.numeric(tau_i2)

                ## Return both derivatives
                return(c(tau_i1,tau_i2))
            }
        }

        ## Find the maximum spatial distance between centroids of all basis functions. This is used for initialisation
        max_l <- max(unlist(Sm@D_basis[[1]]))

        ## Below we actually estimate the parameters
        ## For each resolution
        tau <- lapply(1:nrow(all_res),
                      function(i) {
                          ## Find the basis functions for this resolution
                          idx <- which(Sm@basis@df$res == i)

                          ## Compute the correlation matrix
                          Ki <- Sm@Khat[idx,idx]/Sm@Khat[idx[1],idx[1]]

                          ## If we are in space-time
                          if(is(Sm@basis,"TensorP_Basis")) {
                              ## Extract previous estimate from current covariance matrix.
                              ## If zero (e.g., initial matrix is the identity), then pin to 1e-9
                              ## If only one basis function at this resolution then don't
                              ## attempt to estimate
                              par_init <- ifelse(all_res$n[i]>1,
                                                 max(-Sm@D_basis$Basis1[[i]][1,2]/log(Ki[1,2]),1e-9),
                                                 1e-9)


                              ## Same as above but for temporal (assume we have always more than one temporal basis function)
                              par_init[2] <- max(-Sm@D_basis$Basis2[[1]][1,2]/log(Ki[1,1+count_res(Sm@basis@Basis1)$n[i]]),1e-9) ## time

                              ## If we clamped the temporal length scale then set it initially to 1
                              if(par_init[2] == 1e-9) par_init[2] <- 1
                          } else {
                              ## As above but just for space
                              par_init <- ifelse(all_res$n[i]>1,
                                                 max(-Sm@D_basis[[i]][1,2]/log(Ki[1,2]),1e-9),
                                                 1e-9)
                          }

                          ## If we clamped the spatial length scale then set it initially to max(length) / 10
                          if(par_init[1] == 1e-9) par_init[1] <- max_l/10

                          ## Suppress warnings in case we hit max-iterations. If it hasn't converged we would be
                          ## in a GEM settings which is still OK
                          suppressWarnings(optim(par = par_init,
                                                 fn = f_tau,
                                                 gr = gr_f_tau,
                                                 i=i,control=list(maxit=100L,reltol=1e-4))$par)
                      })

        ## Reconstruct the K matrix based on above parameter estimates
        K <- lapply(1:nrow(all_res),
                    function(i) {
                        if(is(Sm@basis,"TensorP_Basis")) {
                            Ki <- kronecker(exp(-Sm@D_basis$Basis2[[1]]/tau[[i]][2]),
                                            exp(-Sm@D_basis$Basis1[[i]]/tau[[i]][1]))/omega[[i]]
                        } else {
                            Ki <- exp(-Sm@D_basis[[i]]/tau[[i]])/omega[[i]]
                        }

                    })

        ## Since we are in block diagonal mode we can just block-diagonalise across resolutions
        K <- do.call("bdiag",K)

        ## Now, if we have space AND time, block diagonalising by resolution is not correct
        ## as we have the following indices (res1t1....res1tN,res2t1,...,res2tN,...)
        ## This can be corrected by seeing how the indices were in the original data frame
        ## (which were correct by construction), and then permuting the K matrix using
        ## and internal function reverse_permute
        idx_all <- unlist(lapply(1:nrow(all_res),
                                 function(i) which(Sm@basis@df$res == i)))

        # reverse_permute rearranges the order of time/resolution when we have tensor products
        # When we don't have tensor product idx_all and 1:nrow(K) are the same so nothing changes
        K <- reverse_permute(K,idx_all)

        ## If user wants verbose output show estimates
        if( opts_FRK$get("verbose") > 0) {
            cat("  Estimates of omega: ",unlist(omega),"  ")
            cat("  Estimates of tau: ",unlist(tau),"  ")
        }

    }

    ## Return the estimated matrix
    K
}

## The function below regularises the K matrix when the K_type is "unstructured"
.regularise_K <- function(Sm,S_eta= NULL,mu_eta = NULL, lambda = 0) {

    if (is.null(S_eta)) S_eta <- Sm@S_eta      # extract from SRE model if not supplied
    if (is.null(mu_eta)) mu_eta <- Sm@mu_eta   # extract from SRE model if not supplied

    if(any(lambda > 0)) {  # if at least one lambda > 0

        ## If we have just one regulatisation parameter for all resolutions
        if(length(lambda) == 1) {
            reg_matrix <- lambda*Diagonal(nrow(S_eta)) # reg. matrix = lambda*I
        } else {
            ## If we have one regularisation parameter per resolution then the reg. matrix
            ## is diagonal but not proportional to the identity matrix
            ## We use the data frame returned by count_res which has the resolution number
            ## in the first column and the number of basis in the second column
            reg_matrix <- Diagonal(x = do.call("c",
                                               apply(count_res(Sm),1,
                                                     function(x) rep(lambda[x[1]],x[2]))))
        }

        ## Update K but this time regularising
        Q <- chol2inv(chol(S_eta + tcrossprod(mu_eta))) + reg_matrix
        K <- chol2inv(chol(Q))
    } else {
        ## If there is no regularisation then use the following simple update (see vignette for details)
        K <- S_eta + tcrossprod(mu_eta)
    }

    ## Return K
    K
}

.SRE.predict <- function(Sm, obs_fs = FALSE, newdata = NULL,
                         pred_time = NULL, covariances = FALSE, 
                         CP, predict_BAUs) {


    
    ## Get BAUs from the SRE model
    BAUs <- Sm@BAUs

    ## If we have to compute too many covariances then stop and give error
    if(covariances & nrow(CP) > 4000)
        stop("Cannot compute covariances for so many prediction locations. Please reduce
             to less than 4000")

    ## Get the CZ matrix
    CZ <- Sm@Cmat

    ## If the user has specified which polygons he wants we can remove the ones we don't need
    ## We only need those BAUs that are influenced by observations and prediction locations
    ## For ST use all BAUs as it gets complicated
    if(is(newdata,"Spatial")) {

        ## The needed BAUs are the nonzero column indices of CZ and CP
        needed_BAUs <- union(as(CP,"dgTMatrix")@j+1,
                             as(CZ,"dgTMatrix")@j+1)

        ## Filter the BAUs and the matrices
        BAUs <- BAUs[needed_BAUs,]
        CP <- CP[,needed_BAUs]
        CZ <- CZ[,needed_BAUs]
        Sm@S0 <- Sm@S0[needed_BAUs,]
    }

    # Deprecated:
    # if(is(BAUs,"ST")){
    #     needed_BAUs <- BAUs[,pred_time]$n
    #     BAUs <- BAUs[,pred_time]
    #     CP <- CP[,needed_BAUs]
    #     CZ <- CZ[,needed_BAUs]
    # }

    ## Retrieve the dependent variable name
    depname <- all.vars(Sm@f)[1]

    ## Set the dependent variable in BAUs to something just so that .extract.from.formula doesn't
    ## throw an error.. we will NULL it shortly after
    BAUs[[depname]] <- 0.1

    ## Extract covariates from BAUs
    L <- .extract.from.formula(Sm@f,data=BAUs)
    X = as(L$X,"Matrix")
    BAUs[[depname]] <- NULL

    ## Set variables to make code more concise
    S0 <- Sm@S0
    S0 <- as(S0,"dgCMatrix")   # ensure S0 is classified as sparse
    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat
    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta

    if(Sm@fs_model == "ind") {
        D <- sigma2fs*Sm@Vfs + Sm@Ve   # compute total variance (data)
        if(isDiagonal(D)) {
            D <- Diagonal(x=D@x)       # cast to diagonal
            Dchol <- sqrt(D)           # find the inverse of the variance-covariance matrix
            Dinv <- solve(D)           # if Diagonal the Cholesky and inverse are just sqrt and reciprocal
        } else {
            Dchol <- chol(D)           # otherwise they need to be computed in full
            Dinv <- chol2inv(Dchol)
        }
        sig2_Vfs_pred <- Diagonal(x=sigma2fs*BAUs$fs)   # fine-scale variation including estimated factor
        Q <- solve(sig2_Vfs_pred)                       # precision of fine-scale variation
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
            LAMBDAinv <- bdiag(Sm@Khat_inv,Q)                # block diagonal precision matrix
            PI <- cbind(S0, .symDiagonal(n=length(BAUs)))    # PI = [S I]
            tC_Ve_C <- t(CZ) %*% solve(Sm@Ve) %*% CZ +       # summary matrix
                0*.symDiagonal(ncol(CZ))              # Ensure zeros on diagonal
            Qx <- t(PI) %*% tC_Ve_C %*% PI + LAMBDAinv       # conditional precision matrix
            chol_Qx <- cholPermute(as(Qx,"dgCMatrix"))       # permute and do Cholesky
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*%      # Qx = ybar, see vignette
                (Sm@Z - CZ %*% X %*% alpha)
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
            x_mean <- Sm@mu_eta  # conditional mean of eta
            Cov <- Sm@S_eta      # conditional covariance of eta

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
        x_mean <- Sm@mu_eta   # conditional mean
        Cov <- Sm@S_eta       # conditional covariance

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
## between where this cross the y-axis and the measurement error will be the fs-variation (estimated)
.est_obs_error <- function(sp_pts,variogram.formula,vgm_model = NULL,BAU_width = NULL) {

    ## Notify user (even if not verbose == TRUE)
    cat("... Fitting variogram for estimating measurement error\n")

    ## Basic checks
    if(!is(variogram.formula,"formula"))
        stop("variogram.formula needs to be of class formula")
    if(!is(sp_pts,"Spatial"))
        stop("sp_pts needs to be of class Spatial")
    if(!requireNamespace("gstat"))
        stop("gstat is required for variogram estimation. Please install gstat")

    ## Make sure we're not on sphere here, otherwise variogram fitting is too slow. Just remove
    ## CRS (this is only approximate anyways)
    if(!is.na(proj4string(sp_pts))) {
        sp_pts <- SpatialPointsDataFrame(coords = coordinates(sp_pts),
                                         data = sp_pts@data,proj4string = CRS())
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

        ## If we still have problems, then just take the first point of the the empirical semivariogram and
        ## throw a warning that this estimate is probably not very good
        if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
            vgm.fit$psill[1] <- v$gamma[1]
            warning("Estimate of measurement error is probably inaccurate.
                    Please consider setting it through the std variable
                    in the data object if known.")
        }
        }
    cat("sigma2e estimate = ",vgm.fit$psill[1],"\n")

    ## Return the sqrt of the psill as the measurement error
    sp_pts$std <- sqrt(vgm.fit$psill[1])
    sp_pts

        }


## Checks arguments for the SRE() function. Code is self-explanatory
.check_args1 <- function(f, data,basis, BAUs, est_error, 
                         K_type, response, link, taper, fs_by_spatial_BAU, normalise_wts, 
                         sum_variables, average_in_BAU) {
  
    if(!is(f,"formula")) stop("f needs to be a formula.")
    if(!is(data,"list"))
        stop("Please supply a list of Spatial objects.")
    if(!all(sapply(data,function(x) is(x,"Spatial") | is(x,"ST"))))
        stop("All data list elements need to be of class Spatial or ST")
    if(!(is(BAUs,"SpatialPointsDataFrame") | is(BAUs,"SpatialPolygonsDataFrame") | is(BAUs,"SpatialPixelsDataFrame") | is(BAUs,"STFDF")))
        stop("BAUs should be a SpatialPolygonsDataFrame, SpatialPixelsDataFrame, or a STFDF object")
    if(is(BAUs,"STFDF")) if(!(is(BAUs@sp,"SpatialPointsDataFrame") | is(BAUs@sp,"SpatialPolygonsDataFrame") | is(BAUs@sp,"SpatialPixelsDataFrame")))
        stop("The spatial component of the BAUs should be a SpatialPolygonsDataFrame or SpatialPixelsDataFrame")
    if(is(BAUs,"STFDF") && nrow(BAUs@data) != nrow(BAUs@sp) * length(BAUs@time))
        stop("The number of rows in BAUs@data should be equal to the number of spatial BAUs, nrow(BAUs@sp), multiplied by the number of time indices, length(BAUs@time)")
  
    # if(is(BAUs,"SpatialPointsDataFrame"))
    #     stop("Implementation with Point BAUs is currently in progress")
    # if(is(BAUs,"STFDF")) if(is(BAUs@sp,"SpatialPointsDataFrame"))
    #     stop("Implementation with Point BAUs is currently in progress")

    if(!all(all.vars(f)[-1] %in% c(names(BAUs@data),coordnames(BAUs))))
        stop("All covariates need to be in the SpatialPolygons BAU object")
  
    if(any(sapply(data,function(x) any(names(x@data) %in% names(BAUs@data)))))
        stop("Please don't have overlapping variable names in data and BAUs. All covariates need to be in the BAUs")
    if(!all(sapply(data,function(x) all.vars(f)[1] %in% names(x@data))))
        stop("All data list elements to have values for the dependent variable")
    if(!all(sapply(data,function(x) identical(proj4string(x), proj4string(BAUs)))))
        stop("Please ensure all data items and BAUs have the same coordinate reference system")
    if(!(is(basis,"Basis") | is(basis,"TensorP_Basis")))
        stop("basis needs to be of class Basis  or TensorP_Basis (package FRK)")
    if(!("fs" %in% names(BAUs@data))) {
        stop("BAUs should contain a field 'fs' containing a basis
             vector for fine-scale variation. Do BAUs$fs <- 1 if you don't know what this is.")
    }
    if(!(all(BAUs$fs >= 0)))
        stop("fine-scale variation basis function needs to be nonnegative everywhere")
    if((is(manifold(basis),"sphere")) & !all((coordnames(BAUs) %in% c("lon","lat"))))
        stop("Since a sphere is being used, please ensure that
             all coordinates (including those of BAUs) are in (lon,lat)")
    if(response == "gaussian" & !est_error & !all(sapply(data,function(x) "std" %in% names(x@data))))
        stop("If the response is Gaussian and observational error is not going to be estimated,
             please supply a field 'std' in the data objects")

    
    if(!(K_type %in% c("block-exponential", "neighbour", "unstructured", "separable"))) 
        stop("Invalid K_type argument. Please select from 'block-exponential', 'neighbour', 'unstructured', or 'separable'")
    if(K_type == "unstructured" & response != "gaussian")
        stop("The unstructured covariance matrix (K_type = 'unstructured') is not implemented for non-Gaussian response (more specifically, when method = 'TMB')")
    if (K_type == "block-exponential" & response != "gaussian")
        warning("Using the block-exponential covariance matrix (K_type = 'block-exponential') is computationally inefficient with a non-Gaussian response (more specifically, when method = 'TMB'). For these situations, consider using K_type = 'neighbour' or K_type = 'separable'.")
    
  
  if (K_type == "separable" & is(basis,"TensorP_Basis"))
    stop("K_type = 'separable' is not yet implemented in a space-time setting.")
  
  ## If K_type == separable, the basis functions must be in a regular rectangular lattice
  if (K_type == "separable") 
    for (i in unique(basis@df$res)) {
      temp <- basis@df[basis@df$res == i, ]
      if (!.test_regular_grid(temp$loc1, temp$loc2, rectangular = TRUE) ) 
        stop("Basis functions must be arranged in a regular rectangular lattice when K_type = 'separable'.")
    }
  
  
  ## If K_type == "neighbour", we just need basis functions to be in a regular 
  ## lattice (does not need to be rectangular)
    if (K_type == "neighbour") 
      for (i in unique(basis@df$res)) {
        temp <- basis@df[basis@df$res == i, ]
        if (!.test_regular_grid(temp$loc1, temp$loc2, rectangular = FALSE) ) 
          stop("Basis functions must be arranged in a regular lattice when K_type = 'neighbour'.")
      }  
  
    ## Check that valid data model and link function have been chosen
    if (!(response %in% c("gaussian", "poisson", "bernoulli", "gamma", "inverse-gaussian", "negative-binomial", "binomial")))
            stop("Invalid response argument")
    if (!(link %in% c("identity", "log", "square-root", "logit", "probit", "cloglog", "inverse", "inverse-squared")))
            stop("Invalid link argument")
    ## Check that an appropriate response-link combination has been chosen
    if (response == "gaussian" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
            response == "poisson" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
            response == "gamma" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
            response == "inverse-gaussian" & !(link %in% c("identity", "inverse", "log", "inverse-squared", "square-root")) ||
            response == "negative-binomial" & !(link %in% c("log", "square-root", "logit", "probit", "cloglog")) ||
            response == "binomial" & !(link %in% c("logit", "probit", "cloglog")) ||
            response == "bernoulli" & !(link %in% c("logit", "probit", "cloglog"))) {
            stop("Invalid response-link combination selected. Please choose an appropriate link function for the specified response distribution.")
    }
    ## Provide a warning if a possibly problematic combination is chosen
    if (response == "gaussian" & link %in% c("log", "inverse-squared", "square-root") ||
        response == "poisson" & link %in% c("identity", "inverse", "inverse-squared") ||
        response == "gamma" & link %in% c("identity", "inverse", "inverse-squared") ||
        response == "inverse-gaussian" & link %in% c("inverse-squared")) {
        warning("Due to the implied range of the mean function, and the permitted support of the mean for the specified response, nonsensical results are possible with the chosen link function. Consider using a link function which ensures the mean is mapped to the correct support.")
    }
    
    ## Check taper (applicable to block-exponential when method = TMB)
    if (!(class(taper) %in% c("numeric", "integer"))) {
        stop("taper, the argument controlling the coveriance taper, must be numeric or integer.")
    } else if (taper <= 0) {
        stop("taper, the argument controlling the coveriance taper, must be positive.")
    }
    
    ## Check k_Z (size parameter for data)
    if (response %in% c("binomial", "negative-binomial")) {
        if (!all(sapply(data, function(l) "k_Z" %in% names(l)))) {
            stop("For binomial or negative-binomial data, the known constant size parameter must be provided for each observation. Please provide this in the data object, in a field called 'k_Z'.")
        } else if (!all(sapply(data, function(l) class(l$k_Z) %in% c("numeric", "integer")))) {
            stop("The known constant size parameter must contain only positive integers.")
        } else if (any(sapply(data, function(l) l$k_Z <= 0)) | 
                   !all(sapply(data, function(l) l$k_Z == round(l$k_Z)))) {
            stop("The known constant size parameter must contain only positive integers.")
        }
    } else {
      ## it is unlikely to occur, but here we ensure that k_Z and k_BAU 
      ## are not used for other variables: these terms are treated specially.
      if(any(sapply(data,function(x) any("k_Z" %in% names(x@data)))))
        stop("k_Z is a reserved keyword for the size parameter in a binomial or negative-binomial setting; please do not include it as a covariate name in the data objects.")
      if("k_BAU" %in% names(BAUs@data))
        stop("k_BAU is a reserved keyword for the size parameter in a binomial or negative-binomial setting; please do not include it as a covariate name in the BAUs object.")
    }
    
    ## If we wish to have unique fine-scale variance associated with each BAU, 
    ## we need to:
    ##  i) be in a spatio-temporal application
    ##  ii) ensure each spatial BAU is associated with a sufficient number of observations
    ## The second check requires the binned data and the indices of the observed 
    ## BAUs, so we need to check this condition later.
    if (fs_by_spatial_BAU & !is(BAUs, "STFDF")) 
        stop("A unique fine-scale variance can only be associated with each spatial BAU if the application is spatio-temporal (i.e., the BAUs are of class 'STFDF').
              Please either set fs_by_spatial_BAU to FALSE if you are not in a spatio-temporal application.")
  
    if(normalise_wts &
       response %in% c("poisson", "binomial", "bernoulli", "negative-binomial") & 
       any(sapply(data, function(x) is(x, "SpatialPolygons")))) {
        warning("You have specified a count data model with SpatialPolygons observations; consider setting normalise_wts = FALSE so that aggregation of the mean is a weighted sum rather than a weighted average.")
    }
    
    if ("n" %in% sum_variables) 
        stop("Summing the BAU indices will result in out of bounds errors and hence NAs in the incidence matrix C; please remove 'n' from sum_variables.") 
    
    if (!is.null(sum_variables) & !average_in_BAU) 
        warning("sum_variables is not considered when average_in_BAU = FALSE.")
    
}


## Checks arguments for the SRE.fit() function. Code is self-explanatory
.check_args2 <- function(n_EM, tol, lambda, method, print_lik, optimiser, 
                         response, K_type, link, fs_by_spatial_BAU, known_sigma2fs, 
                         BAUs, ...) {
  
    if(!is.numeric(n_EM)) stop("n_EM needs to be an integer")
    if(!(n_EM <- round(n_EM)) > 0) stop("n_EM needs to be greater than 0")
    if(!is.numeric(tol)) stop("tol needs to be a number greater than zero")
    if(!(tol > 0)) stop("tol needs to be a number greater than zero")
    if(!(is.logical(print_lik))) stop("print_lik needs to be a logical quantity")
    if(!(is.numeric(lambda))) stop("lambda needs to be a number")
    if(!(all(lambda >= 0))) stop("lambda needs to be greater or equal to zero")
    
    if(!(method %in% c("EM", "TMB"))) stop("Currently only the EM algorithm or TMB are implemented for parameter estimation.")
    
    if(method == "EM" & !(response == "gaussian")) stop("The EM algorithm is only available for response = 'gaussian'. Please use method = 'TMB' for all other assumed response distributions.")
    if(method == "EM" & !(link == "identity")) stop("The EM algorithm is only available for link = 'identity'. Please use method = 'TMB' for all other link functions.")
    if(method == "EM" & K_type == "neighbour") stop("The neighbour matrix formulation of the model is not implemented for method = 'EM'. Please choose K_type to be 'block-exponential' or 'unstructured'.")
    if(method == "EM" & K_type == "separable") stop("The separable spatial model is not implemented for method = 'EM'. Please choose K_type to be 'block-exponential' or 'unstructured'.")
    if(method == "TMB" & K_type == "unstructured") stop("The unstructured covariance matrix (K_type = 'unstructured') is not implemented for method = 'TMB'")
    if(method != "TMB" & fs_by_spatial_BAU) stop("fs_by_spatial_BAU can only be TRUE if method = 'TMB'. Please set method = 'TMB', or fs_by_spatial_BAU = FALSE.")
  
  
  ## Check known_sigma2fs 
  ns <- dim(BAUs)[1]
  if (!is.null(known_sigma2fs)) {
    if (any(known_sigma2fs < 0))
      stop("known_sigma2fs should contain only positive elements.")
    ## Check that the known_sigma2fs is of length 1, 
    ## or of length equal to the number of spatial BAUs
    if ((fs_by_spatial_BAU & length(known_sigma2fs) != ns) ||
        (!fs_by_spatial_BAU & length(known_sigma2fs) != 1))
      stop("The length of known_sigma2fs should be equal to the number of spatial BAUs when fs_by_spatial_BAU = TRUE, and 1 otherwise")
  }
  
  
    ## Check optional parameters to optimiser function are ok:
    l <- list(...)
    
    # The ellipsis argument in the wrapper FRK() allows users to pass in
    # optional parameters to auto_BAUs() and auto_basis(). This means that
    # we must check whether the parameters in ... match the formal parameters
    # of optimiser, auto_BAUs, and auto_basis.
    # This is not ideal, as the ... argument in SRE.fit() is present solely
    # for optimiser(), and so when SRE.fit() is called, we should really only
    # try to match with the parameters of optimiser.
    valid_param_names <- c(names(formals(auto_BAUs)), 
                           names(formals(auto_basis)), 
                           names(formals(optimiser)))

    ## If any of the formals have the same argument - omit the ellipsis "..." from this check
    optimiser_arg_in_other_functions <- names(formals(optimiser))[names(formals(optimiser)) != "..."] %in% c(names(formals(auto_BAUs)), names(formals(auto_basis)))
    if(any(optimiser_arg_in_other_functions)) {
        msg1 <- "The optimiser cannot have arguments which overlap with those of auto_basis() or auto_BAUs(). Offending arguments:"
        msg2 <- (names(formals(optimiser))[names(formals(optimiser)) != "..."])[which(optimiser_arg_in_other_functions)]
        stop(paste(c(msg1, msg2), collapse = " "))
    }
        
    if(!all(names(l) %in% valid_param_names)) {
        msg1 <- "Optional arguments for auto_BAUs(), auto_basis(), or optimiser function not matching declared arguments. It is also possible that you have misspelt an argument to SRE.fit() or FRK(). Offending arguments:"
        msg2 <- names(l)[!(names(l) %in% valid_param_names)]
        stop(paste(c(msg1, msg2), collapse = " "))
    }
}

## Checks arguments for the predict() function. Code is self-explanatory
.check_args3 <- function(obs_fs, newdata, pred_polys,
                         pred_time, covariances, SRE_model, type, 
                         k, percentiles, kriging, ...) {
  
    if(kriging != "simple" & SRE_model@method == "EM")
      stop("Universal kriging is only available when method = 'TMB'")
    
    if(!(obs_fs %in% 0:1)) stop("obs_fs needs to be logical")

    if(!(is(newdata,"Spatial") | is(newdata,"ST") | is.null(newdata)))
        stop("Predictions need to be over Spatial or ST objects")

    if(!is.null(newdata) & !is.null(pred_time))
        stop("Only one of newdata and pred_time can be not NULL")

    if(!is.null(pred_polys))
        warning("pred_polys is deprecated. Please use newdata instead")

    if(!(is.integer(pred_time) | is.null(pred_time))) stop("pred_time needs to be of class integer")
    if(!is.logical(covariances)) stop("covariances needs to be TRUE or FALSE")
    
    ## Quantities of interest
    if(!all(type %in% c("link", "mean", "response")))
        stop("type must be a vector containing combinations of 'link', 'mean', and 'response'")
    
    ## Check k (for predictions)
    if (SRE_model@response %in% c("binomial", "negative-binomial")) {
        if(length(k) == 1){
            warning("Single number k provided for all BAUs: assuming k is invariant over the whole spatial domain.")
        } else if (!(class(k) %in% c("numeric", "integer"))) {
            stop("k must contain only positive integers.")
        } else if (any(k < 0) | any(k != round (k))) {
            stop("k must contain only positive integers.")
        } else if (length(k) != nrow(SRE_model@S0)) {
            ## FIXME: If we are allowing the user to predict over arbitrary polygons specified
            ## by the argument "newdata", then k will not be equal to the number of BAUs. Perhaps
            ## we should check to see if newdata is null, and then test if k equals its length. 
            stop("length(k) must equal 1 or N (the number of BAUs)." )
        }      
    }
    
    ## Check requested percentiles 
    if (!is.null(percentiles) & !(class(percentiles) %in% c("numeric", "integer"))) 
        stop("percentiles must either be NULL or a numeric or integer vector with entries between 0 and 100.")
    else if (!is.null(percentiles)) {
        if (min(percentiles) < 0 | max(percentiles) > 100) 
            stop("percentiles must be a vector with entries between 0 and 100")   
    }
    
    # ## Check credible mass for HPD interval
    # if (!is.null(cred_mass)) {
    #     if (length(cred_mass) != 1)
    #         stop("cred_mass should be a single scalar, not a vector")
    #     if(any(cred_mass < 0 | cred_mass > 1))
    #         stop("Elements of cred_mass should be between 0 and 1")
    # }
}



## The following function is the internal prediction function
.SRE.predict.deprecated <- function(Sm, obs_fs = FALSE, pred_polys = NULL,
                                    pred_time = NULL, covariances = FALSE) {

    ## If the user does not specify time points to predict at when in space-time
    ## Then predict at every time point
    if(is.null(pred_time) & is(Sm@BAUs,"ST"))
        pred_time <- 1:length(time(Sm@BAUs))

    ## We start by assuming that we will predict at BAUs
    predict_BAUs <- TRUE

    ## Get BAUs from the SRE model
    BAUs <- Sm@BAUs

    ## If the user has not specified polygons over which to predict, then CP is
    ## just the diagonal matrix and we predict over all the BAUs
    if(is.null(pred_polys) & is.null(pred_time)) {
        CP <- Diagonal(length(BAUs))
    }
    if (!is.null(pred_polys)) {
        ## The user has maybe specified a subset of (could be all) the BAUs over which to predict.
        ## The following checks whether pred_polys is a subset of the BAUs through the row names
        pred_polys_are_BAUs <- all(row.names(pred_polys) %in% row.names(BAUs))

        ## If the user has specified a subset of BAUs
        if(pred_polys_are_BAUs) {
            ## See which BAUs the user has specified
            BAUs_idx <- match(row.names(pred_polys), row.names(BAUs))

            ## Construct an incidence matrix that picks out these BAUs
            CP <-  sparseMatrix(i = 1:length(pred_polys),
                                j = BAUs_idx,
                                x = 1,
                                dims = c(length(pred_polys),
                                         length(BAUs)))

        } else {
            ## The user has specified arbitrary polygons
            ## First try to coerce what the user supplied to Polygons (not pixels etc.)
            ## Recall that for now only Spatial pred_polys are allowed so the following is
            ## always valid
            pred_polys <- as(pred_polys,"SpatialPolygonsDataFrame")

            ## Based on these polygons construct the C matrix
            C_idx <- BuildC(pred_polys,BAUs)
            CP <- sparseMatrix(i=C_idx$i_idx,
                               j=C_idx$j_idx,
                               x=1,
                               dims=c(length(pred_polys),
                                      length(BAUs)))

            ## As in SRE(), make sure the polgons are averages (not sums)
            CP <- CP / rowSums(CP)

            ## If even one polygon encompasses more than one BAU, then we need to
            ## predict over linear combinations of BAUs, and hence need to
            ## compute the full covariance matrix. Note this by setting
            ## predict_BAUs <- FALSE
            if(!all(table(C_idx$i_idx) == 1))
                predict_BAUs <- FALSE   ## Need to compute full covariance matrix
        }
    }

    if(!is.null(pred_time)) {
        nspace <- dim(BAUs)[1]
        ntime <- length(pred_time)
        t1 <- ((pred_time[1] - 1)*nspace + 1)
        t2 <- (pred_time[ntime] *nspace)
        CP <- sparseMatrix(i = 1:(t2 - t1 + 1), j = t1:t2, x = 1,
                           dims = c((t2 -t1 + 1), length(BAUs)))
        pred_polys <- BAUs[,pred_time]
    }

    ## If we have to compute too many covariances then stop and give error
    if(covariances & nrow(CP) > 4000)
        stop("Cannot compute covariances for so many observations. Please reduce
             to less than 4000")

    ## Get the CZ matrix
    CZ <- Sm@Cmat

    ## If the user has specified which polygons he want we can remove the ones we don't need
    ## We only need those BAUs that are influenced by observations and prediction locations
    ## For space-time use all the BAUs
    if(!is.null(pred_polys)) {

        ## The needed BAUs are the nonzero column indices of CZ and CP
        if(!is.null(pred_time)) {
            needed_BAUs <- 1:length(BAUs)
        } else {
            needed_BAUs <- union(as(CP,"dgTMatrix")@j+1,
                                 as(CZ,"dgTMatrix")@j+1)
            ## Filter the BAUs and the matrices
            BAUs <- BAUs[needed_BAUs,]
            CP <- CP[,needed_BAUs]
            CZ <- CZ[,needed_BAUs]
            Sm@S0 <- Sm@S0[needed_BAUs,]
        }
    }

    # Deprecated:
    # if(is(BAUs,"ST")){
    #     needed_BAUs <- BAUs[,pred_time]$n
    #     BAUs <- BAUs[,pred_time]
    #     CP <- CP[,needed_BAUs]
    #     CZ <- CZ[,needed_BAUs]
    # }

    ## Retrieve the dependent variable name
    depname <- all.vars(Sm@f)[1]

    ## Set the dependent variable in BAUs to something just so that .extract.from.formula doesn't
    ## throw an error.. we will NULL it shortly after
    BAUs[[depname]] <- 0.1

    ## Extract covariates from BAUs
    L <- .extract.from.formula(Sm@f,data=BAUs)
    X = as(L$X,"Matrix")
    BAUs[[depname]] <- NULL

    ## Set variables to make code more concise
    S0 <- Sm@S0
    S0 <- as(S0,"dgCMatrix")   # ensure S0 is classified as sparse
    alpha <- Sm@alphahat
    K <- Sm@Khat
    sigma2fs <- Sm@sigma2fshat
    mu_eta <- Sm@mu_eta
    S_eta <- Sm@S_eta

    if(Sm@fs_model == "ind") {
        D <- sigma2fs*Sm@Vfs + Sm@Ve   # compute total variance (data)
        if(isDiagonal(D)) {
            D <- Diagonal(x=D@x)       # cast to diagonal
            Dchol <- sqrt(D)           # find the inverse of the variance-covariance matrix
            Dinv <- solve(D)           # if Diagonal the Cholesky and inverse are just sqrt and reciprocal
        } else {
            Dchol <- chol(D)           # otherwise they need to be computed in full
            Dinv <- chol2inv(Dchol)
        }
        sig2_Vfs_pred <- Diagonal(x=sigma2fs*BAUs$fs)   # fine-scale variation including estimated factor
        Q <- solve(sig2_Vfs_pred)                       # precision of fine-scale variation
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
            LAMBDAinv <- bdiag(Sm@Khat_inv,Q)                # block diagonal precision matrix
            PI <- cbind(S0, .symDiagonal(n=length(BAUs)))    # PI = [S I]
            tC_Ve_C <- t(CZ) %*% solve(Sm@Ve) %*% CZ +       # summary matrix
                0*.symDiagonal(ncol(CZ))              # Ensure zeros on diagonal
            Qx <- t(PI) %*% tC_Ve_C %*% PI + LAMBDAinv       # conditional precision matrix
            chol_Qx <- cholPermute(as(Qx,"dgCMatrix"))       # permute and do Cholesky
            ybar <- t(PI) %*%t(CZ) %*% solve(Sm@Ve) %*%      # Qx = ybar, see vignette
                (Sm@Z - CZ %*% X %*% alpha)
            x_mean <- cholsolve(Qx,ybar,perm=TRUE,           # solve Qx = ybar using permutations
                                cholQp = chol_Qx$Qpermchol, P = chol_Qx$P)

            if(predict_BAUs) {
                ## If we are predicting at BAUs then we only need a few covariance elements.
                ## Find these elements
                Cov <- Takahashi_Davis(Qx,cholQp = chol_Qx$Qpermchol,P = chol_Qx$P)

                ## Compute the variance and std over the BAUs in batches
                BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = TRUE)
            } else {
                ## Do not compute covariance now, do it later
                #  Cov <- cholsolve(Qx,Diagonal(nrow(Qx)),perm=TRUE,
                #                   cholQp = chol_Qx$Qpermchol, P = chol_Qx$P) # FULL
            }


        } else {
            ## If sigma2fs = 0 then the prediction is much simpler and all our
            ## predictions / uncertainty come from the random effects
            PI <- S0
            x_mean <- Sm@mu_eta  # conditional mean of eta
            Cov <- Sm@S_eta      # conditional covariance of eta

            ## Compute variances, this time indicating there is no fs variation in process
            BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = FALSE)
        }

        ## The conditional mean is simply given by fitted random effects + fitted fixed effects
        BAUs[["mu"]] <- as.numeric(X %*% alpha + PI %*% x_mean)

    }

    ## Case 1 (fs variation in measurement equation)
    if(obs_fs) {
        ## All predictions and prediction uncertainties comes from our prediction of and uncertainty over eta
        x_mean <- Sm@mu_eta   # conditional mean
        Cov <- Sm@S_eta       # conditional covariance

        ## Compute variances, this time indicating there is no fs variation in process
        BAUs[["mu"]] <- as.numeric(X %*% alpha) + as.numeric(S0 %*% x_mean)
        BAUs[["var"]] <- .batch_compute_var(S0,Cov,fs_in_process = FALSE)
    }


    ## Now, if the user hasn't specified prediction polygons, and the user does not want covariances,
    ## our job is done and we just return the BAUs, possibly at selected time points
    if(predict_BAUs & !covariances) {
        BAUs[["sd"]] <- sqrt(BAUs[["var"]])  # compute the standard error
        if(!is.null(pred_polys)) { # User had specified a specific set of BAUs. Return only these (spatial only for now)
            BAUs <- BAUs[row.names(pred_polys),]
        }
        if(!is.null(pred_time)) BAUs <- BAUs[,pred_time]  # return only specified time points
        BAUs

    } else {
        ## Otherwise we need to find the mean and variance of linear combinations of these BAUs

        ## The linear combination of the mean is easy
        pred_polys[["mu"]] <- as.numeric(CP %*% BAUs[["mu"]])

        ## If we have fs variation in the process layer we need to consider the
        ## fine-scale variation (PI = [S I]) when predicting over the polygons,
        ## otherwise we just need the variance over eta
        if(!obs_fs & sigma2fs > 0) CPM <- CP %*% PI else CPM <- CP %*% S0

        ## If there is no fine-scale variation then simply find linear combination
        if(sigma2fs == 0 | obs_fs)  {
            pred_polys[["var"]] <- diag2(CPM %*% Cov, t(CPM)) ## All Cov available
            if(covariances) {
                L <- t(chol(Cov))
                Covariances <- tcrossprod(CPM %*% L)
            }
        }

        ## Otherwise find full covariance matrix (including over fs-variation). This is a last-case resort and
        ## may crash the computer if there are several prediction polygons. However this shouldn't be the case
        ## if these polygons span multiple BAUs
        else {
            pred_polys[["var"]] <- diag2(CPM, cholsolve(Q=Qx,y=t(CPM),
                                                        perm = TRUE,
                                                        cholQp = chol_Qx$Qpermchol,
                                                        P = chol_Qx$P))
            if(covariances) Covariances <- cholsolveAQinvAT(A = CPM,
                                                            Lp = chol_Qx$Qpermchol,
                                                            P = chol_Qx$P)
        }

        # Compute standard error
        pred_polys[["sd"]] <- sqrt(pred_polys[["var"]])

        ## Return the prediction polygons
        if(covariances) {
            pred_polys <- list(pred_polys = pred_polys,
                               Cov = Covariances)
        }
        pred_polys
    }

}
