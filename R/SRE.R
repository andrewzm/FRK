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
#' 
#' @title Construct SRE object, fit and predict
#' @description The Spatial Random Effects (SRE) model is the central object in FRK. The function \code{FRK} provides a wrapper for the construction and estimation of the SRE object from data, using the functions \code{SRE} (the object constructor) and \code{SRE.fit} (for fitting it to the data). Please see \code{\link{SRE-class}} for more details on the SRE object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame}, \code{SpatialPolygonsDataFrame}, \code{STIDF}, or  \code{STFDF}. If using space-time objects, the data frame must have another field, \code{t}, containing the time index of the data point. If the assumed response distribution is "binomial" or "negative-binomial", the data frame must have another field, \code{k}, containing the known constant parameter \eqn{k} for each observation. 
#' @param basis object of class \code{Basis} (or \code{TensorP_Basis})
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame}, \code{STIDF}, or \code{STFDF}. The object's data frame must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality. If the function \code{FRK()} is used directly, then BAUs are created automatically, but only coordinates can then be used as covariates
#' @param est_error (applicable only if \code{response} = "gaussian") flag indicating whether the measurement-error variance should be estimated from variogram techniques. If this is set to 0, then \code{data} must contain a field \code{std}. Measurement-error estimation is currently not implemented for spatio-temporal datasets
#' @param include_fs (applicable only if \code{method} = "TMB") flag indicating whether the fine-scale variation should be included in the model
#' @param average_in_BAU if \code{TRUE}, then multiple data points falling in the same BAU are averaged; the measurement error of the averaged data point is taken as the average of the individual measurement errors
#' @param sum_variables if \code{average_in_BAU == TRUE}, the string \code{sum_variables} indicates which data variables (can be observations or covariates) are to be summed rather than averaged
#' @param normalise_wts if \code{TRUE}, the rows of the incidence matrices \eqn{C_Z} and \eqn{C_P} are normalised to sum to 1, so that the mapping represents a weighted average; if false, no normalisation of the weights occurs (i.e., the mapping corresponds to a weighted sum)
#' @param fs_by_spatial_BAU (applicable only in a spatio-temporal setting and if \code{method} = "TMB") if \code{TRUE}, then each spatial BAU is associated with its own fine-scale variance parameter; otherwise, a single fine-scale variance parameter is used
#' @param fs_model if "ind" then the fine-scale variation is independent at the BAU level. Only the independent model is allowed for now, future implementation will include CAR/ICAR (in development)
# #' If "ICAR", then an ICAR model for the fine-scale variation is placed on the BAUs
#' @param vgm_model (applicable only if \code{response} = "gaussian") an object of class \code{variogramModel} from the package \code{gstat} constructed using the function \code{vgm}. This object contains the variogram model that will be fit to the data. The nugget is taken as the measurement error when \code{est_error = TRUE}. If unspecified, the variogram used is \code{gstat::vgm(1, "Lin", d, 1)}, where \code{d} is approximately one third of the maximum distance between any two data points
#' @param K_type the parameterisation used for the basis-function covariance matrix, \code{K}. If \code{method} = "EM", \code{K_type} can be "unstructured" or "block-exponential". If \code{method} = "TMB", \code{K_type} can be "precision" or "block-exponential". The default is "block-exponential", however if \code{FRK()} is used and \code{method} = "TMB", for computational reasons \code{K_type} is set to "precision"
#' @param normalise_basis flag indicating whether to normalise the basis functions so that they reproduce a stochastic process with approximately constant variance spatially
#' @param object object of class \code{SRE} returned from the constructor \code{SRE()} containing all the parameters and information on the SRE model. Note that prior to v2.x, \code{loglik()} and \code{SRE.fit()} took the now-defunct argument \code{SRE_model} instead of \code{object}
#' @param n_EM (applicable only if \code{method} = "EM") maximum number of iterations for the EM algorithm
#' @param tol (applicable only if \code{method} = "EM") convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently "EM" and "TMB" are supported
#' @param lambda (applicable only if \code{K_type} = "unstructured") ridge-regression regularisation parameter (0 by default). Can be a single number, or a vector (one parameter for each resolution)
#' @param print_lik (applicable only if \code{method} = "EM") flag indicating whether to plot log-likelihood vs. iteration after convergence of the EM estimation algorithm
# #' @param use_centroid flag indicating whether the basis functions are averaged over the BAU, or whether the basis functions are evaluated at the BAUs centroid in order to construct the matrix \eqn{S}. The flag can safely be set when the basis functions are approximately constant over the BAUs in order to reduce computational time
#' @param newdata object of class \code{SpatialPoylgons}, \code{SpatialPoints}, or \code{STI}, indicating the regions or points over which prediction will be carried out. The BAUs are used if this option is not specified. 
#' @param obs_fs flag indicating whether the fine-scale variation sits in the observation model (systematic error; indicated by \code{obs_fs = TRUE}) or in the process model (process fine-scale variation; indicated by \code{obs_fs = FALSE}, default). For non-Gaussian data models, and/or non-identity link functions, if \code{obs_fs = TRUE}, then the fine-scale variation is removed from the latent process \eqn{Y}; however, they are re-introduced for computation of the conditonal mean \eqn{\mu} and response variable \eqn{Z}
#' @param pred_time vector of time indices at which prediction will be carried out. All time points are used if this option is not specified
#' @param covariances (applicable only for \code{method} = "EM") logical variable indicating whether prediction covariances should be returned or not. If set to \code{TRUE}, a maximum of 4000 prediction locations or polygons are allowed
#' @param response string indicating the assumed distribution of the response variable. It can be "gaussian", "poisson", "negative-binomial", "binomial", "gamma", or "inverse-gaussian". If \code{method} = "EM", only "gaussian" can be used
#' @param link  string indicating the desired link function. Can be "log", "identity", "logit", "probit", "cloglog", "reciprocal", or "reciprocal-squared". Note that only sensible link-function and response-distribution combinations are permitted. If \code{method} = "EM", only "identity" can be used
#' @param taper positive numeric indicating the strength of the covariance/partial-correlation tapering. Only applicable if \code{K_type} = "block-exponential", or if \code{K_type} = "precision" and the the basis-functions are irregular or the manifold is not the plane. If \code{taper} is \code{NULL} (default) and \code{method} = "EM", no tapering is applied; if \code{method} = "TMB", tapering must be applied (for computational reasons), and we set it to 3 if it is unspecified
#' @param optimiser (applicable only if \code{method} = "TMB") the optimising function used for model fitting when \code{method} = "TMB" (default is \code{nlminb}). Users may pass in a function object or a string corresponding to a named function. Optional parameters may be passed to \code{optimiser} via \code{...}. The only requirement of \code{optimiser} is that the first three arguments correspond to the initial parameters, the objective function, and the gradient, respectively (this may be achieved by simply constructing a wrapper function) 
#' @param known_sigma2fs known value of the fine-scale variance parameter. If \code{NULL} (the default), the fine-scale variance parameter is estimated as usual. If \code{known_sigma2fs} is not \code{NULL}, the fine-scale variance is fixed to the supplied value; this may be a scalar, or vector of length equal to the number of spatial BAUs (if \code{fs_by_spatial_BAU = TRUE})
#' @param kriging (applicable only if \code{method} = "TMB") string indicating the kind of kriging: "simple" ignores uncertainty due to estimation of the fixed effects, while "universal" accounts for this source of uncertainty
#' @param type (applicable only if \code{method} = "TMB") vector of strings indicating the quantities for which inference is desired. If "link" is in \code{type}, inference on the latent Gaussian process \eqn{Y(\cdot)} is included; if "mean" is in \code{type}, inference on the mean process \eqn{\mu(\cdot)} is included (and the probability process, \eqn{\pi(\cdot)},  if applicable); if "response" is in \code{type}, inference on the noisy data process \eqn{Z} is included. Only applicable if \code{method} = "TMB"
#' @param n_MC (applicable only if \code{method} = "TMB") a positive integer indicating the number of MC samples at each location
#' @param k (applicable only if \code{response} is "binomial" or "negative-binomial") vector of size parameters at each BAU
#' @param percentiles (applicable only if \code{method} = "TMB") a vector of scalars in (0, 100) specifying the desired percentiles of the posterior predictive distribution; if \code{NULL}, no percentiles are computed
#' @param ... other parameters passed on to \code{auto_basis()} and \code{auto_BAUs()} when calling \code{FRK()}, or the user specified function \code{optimiser()} when calling \code{FRK()} or \code{SRE.fit()}
#' @details 
#' \strong{Model description}
#' 
#' The hierarchical model implemented in \code{FRK} is a spatial generalised 
#' linear mixed model (GLMM), which may be summarised as
#' \deqn{Z_j \mid \mu_{Z,j}, \psi \sim EF(\mu_{Z, j}, \psi)}
#' \deqn{\mu_Z = C\mu}
#' \deqn{g(\mu) = Y}
#' \deqn{Y = T\alpha + S\eta + \xi}
#' \deqn{\alpha \mid \theta \sim N(0, K)}
#' \deqn{\xi \mid \sigma^2_\xi \sim N(0, \sigma^2_\xi V),}
#' where \eqn{Z_j} denotes a datum, \eqn{EF(\cdot)} denotes an exponential 
#' family member with mean parameter \eqn{\mu_{Z, j}}, \eqn{\mu} is the mean 
#' process evaluated over the BAUs, \eqn{g(\cdot)} is a link function that links
#' the mean process \eqn{\mu(\cdot)} to the latent Gaussian process \eqn{Y(\cdot)},  
#' \eqn{Y} is the latent Gaussian process evaluated over the BAUs, \eqn{T} are 
#' regression covariates at the BAU level associated with regression parameters
#' \eqn{\alpha}, \eqn{S} is a matrix of basis function evaluations over the BAUs, 
#' \eqn{\eta} are the random coefficients associated with the basis functions, and \eqn{\xi} is 
#' a vector containing fine-scale variation at the BAU level. The prior 
#' distribution of the basis-function coefficients, \eqn{\eta}, are formulated 
#' using either a covariance or precision matrix, depending on the argument 
#' \code{K_type}; the parameters of these matrices are estimated during model 
#' fitting. The covariance matrix of \eqn{\xi} is diagonal, with its 
#' diagonal elements proportional to the field `fs' in the 
#' BAUs (typically set to one). The constant of proportionality is estimated 
#' during model fitting. 
#' 
#' When the data is Gaussian, and an identity link function is used, the preceding 
#' model simplifies considerably: specifically,
#' \deqn{Z = CY + C\delta + e,} 
#' where \eqn{Z} is the data vector, \eqn{\delta} is systematic error at the 
#' BAU level, and \eqn{e} represents independent measurement error.  
#' 
#' \strong{Set-up}
#' 
#' \code{SRE()} 
#' constructs a spatial random effects model from the user-defined formula, data object (a list 
#' of spatially-referenced data), basis functions and a set of Basic Areal Units (BAUs). 
#' It first takes each object in the list \code{data} and maps it to the BAUs -- this 
#' entails binning point-referenced data into the BAUs (and averaging within the 
#' BAU if \code{average_in_BAU = TRUE}), and finding which BAUs are associated 
#' with observations. Following this, the incidence matrix, \eqn{C}, is 
#' constructed. 
#' All required matrices (\eqn{S}, \eqn{T}, \eqn{C}, etc.) 
#' are constructed within \code{SRE()} and returned as part of the \code{SRE} object. 
#' \code{SRE()} also intitialises the parameters and random effects using 
#' sensible defaults. Please see 
#' \code{\link{SRE-class}} for more details. 
#' The functions \code{observed_BAUs()} and \code{unobserved_BAUs()} return the 
#' indices of the observed and unobserved BAUs, respectively. 
#' 
#' \strong{Model fitting}
#'
#' \code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown
#' parameters, namely the covariance matrix \eqn{K}, the fine scale variance
#' (\eqn{\sigma^2_{\xi}} or \eqn{\sigma^2_{\delta}}, depending on whether Case 1
#' or Case 2 is chosen; see the vignette) and the regression parameters \eqn{\alpha}.
#' There are two methods of model fitting currently implemented, both of which 
#' implement maximum likelihood estimation (MLE).
#' \itemize{
#'  \item{MLE via the expectation maximisation
#'  (EM) algorithm. }{This method is implemented only
#'  for Gaussian data and an identity link function.
#'  The log-likelihood (given in Section 2.2 of the vignette) is evaluated at each
#' iteration at the current parameter estimate. Optimation continues until
#' convergence is reached (when the log-likelihood stops changing by more than
#' \code{tol}), or when the number of EM iterations reaches \code{n_EM}.
#' The actual computations for the E-step and M-step are relatively straightforward.
#' The E-step contains an inverse of an \eqn{r \times r} matrix, where \code{r}
#' is the number of basis functions which should not exceed 2000. The M-step
#' first updates the matrix \eqn{K}, which only depends on the sufficient
#' statistics of the basis-function coefficients \eqn{\eta}. Then, the regression
#' parameters \eqn{\alpha} are updated and a simple optimisation routine
#' (a line search) is used to update the fine-scale variance
#' \eqn{\sigma^2_{\delta}} or \eqn{\sigma^2_{\xi}}. If the fine-scale errors and
#' measurement random errors are homoscedastic, then a closed-form solution is
#' available for the update of \eqn{\sigma^2_{\xi}} or \eqn{\sigma^2_{\delta}}.
#' Irrespectively, since the updates of \eqn{\alpha}, and \eqn{\sigma^2_{\delta}}
#' or \eqn{\sigma^2_{\xi}}, are dependent, these two updates are iterated until
#' the change in \eqn{\sigma^2_{\cdot}} is no more than 0.1\%.}
#'  \item{MLE via \code{TMB}. }{This method is implemented for
#'  all available data models and link functions offered by \code{FRK}. Furthermore,
#'  this method faciliates the inclusion of many more basis function than possible
#'  with the EM algorithm (in excess of 10,000). \code{TMB} applies
#'  the Laplace approximation to integrate out the latent random effects from the
#'  complete-data likelihood. The resulting approximation of the marginal
#'  log-likelihood, and its derivatives with respect to the parameters, are then
#'  called from within \code{R} using the optimising function \code{optimiser}
#'  (default \code{nlminb()}).}
#' }
#' 
#' \code{info_fit()} extracts information on the fitting (convergence, etc.), 
#' \code{coef()} extracts the estimated regression regression coefficients, and 
#' \code{loglik()} returns the final log-likelihood. 
#' 
#' \emph{Wrapper for set-up and model fitting}
#'
#' The function \code{FRK()} acts as a wrapper for the functions \code{SRE()} and 
#' \code{SRE.fit()}. An added advantage of using \code{FRK()} directly is that it 
#' automatically generates BAUs and basis functions based on the data. Hence 
#' \code{FRK()} can be called using only a list of data objects and an \code{R} 
#' formula, although the \code{R} formula can only contain space or time as 
#' covariates when BAUs are not explicitly supplied with the covariate data.
#' 
#'
#' \strong{Prediction}
#'
#' Once the parameters are estimated, the \code{SRE} object is passed onto the 
#' function \code{predict()} in order to carry out optimal predictions over the 
#' same BAUs used to construct the SRE model with \code{SRE()}. The first part 
#' of the prediction process is to construct the matrix \eqn{S} over the 
#' prediction polygons. This is made computationally efficient by treating the 
#' prediction over polygons as that of the prediction over a combination of BAUs. 
#' This will yield valid results only if the BAUs are relatively small. Once the 
#' matrix \eqn{S} is found, a standard Gaussian inversion (through conditioning) 
#' using the estimated parameters is used for prediction.
#'
#' \code{predict} returns the BAUs (or an object specified in \code{newdata}), 
#' which are of class \code{SpatialPixelsDataFrame}, \code{SpatialPolygonsDataFrame}, 
#' or \code{STFDF}, with predictions and 
#' uncertainty quantification added. 
#' If \code{method} = "TMB", the returned object is a list, containing the 
#' previously described predictions, and a list of Monte Carlo samples. 
#' The predictions and uncertainties can be easily plotted using \code{\link{plot}}
#' or \code{spplot} from the package \code{sp}.
#' @seealso \code{\link{SRE-class}} for details on the SRE object internals, 
#' \code{\link{auto_basis}} for automatically constructing basis functions, and
#' \code{\link{auto_BAUs}} for automatically constructing BAUs. 
#' See also the paper \url{https://arxiv.org/abs/1705.08105} for details on code operation.
#' @export
#' @examples
#' library("sp")        
#' ## Generate process and data
#' m <- 250                                                   # Sample size
#' zdf <- data.frame(x = runif(m), y= runif(m))               # Generate random locs
#' zdf$Y <- 3 + sin(7 * zdf$x) + cos(9 * zdf$y)               # Latent process
#' zdf$z <- rnorm(m, mean = zdf$Y)                            # Simulate data
#' coordinates(zdf) = ~x+y                                    # Turn into sp object
#' 
#' ## Construct BAUs and basis functions 
#' BAUs <- auto_BAUs(manifold = plane(), data = zdf, 
#'                   nonconvex_hull = FALSE, cellsize = c(0.03, 0.03), type="grid") 
#' BAUs$fs <- 1 # scalar fine-scale covariance matrix
#' basis <- auto_basis(manifold =  plane(), data = zdf, nres = 2)
#' 
#' ## Fit the SRE model
#' S <- SRE(f = z ~ 1, list(zdf), basis = basis, BAUs = BAUs)
#' 
#' ## Compute observed and unobserved BAUs    
#' observed_BAUs(S)
#' unobserved_BAUs(S)   
#' 
#' ## Fit with 5 EM iterations so as not to take too much time
#' S <- SRE.fit(S,n_EM = 5, tol = 0.01, print_lik = TRUE)
#' 
#' ## Check fit info, final log-likelihood, and estimated regression coefficients
#' info_fit(S)
#' loglik(S)
#' coef(S)
#' 
#' ## Predict over BAUs
#' pred <- predict(S)
#' 
#' ## Plot
#' \dontrun{
#' plot_spatial_or_ST(zdf, "z")
#' plotlist <- plot(S, pred)
#' ggpubr::ggarrange(plotlist = plotlist)}
SRE <- function(f, data,basis,BAUs, est_error = TRUE, average_in_BAU = TRUE,
                sum_variables = NULL,
                normalise_wts = TRUE,
                fs_model = "ind", vgm_model = NULL, 
                K_type = c("block-exponential", "precision", "unstructured"), 
                normalise_basis = TRUE, 
                response = c("gaussian", "poisson", "bernoulli", "gamma",
                             "inverse-gaussian", "negative-binomial", "binomial"), 
                link = c("identity", "log", "square-root", "logit", "probit", "cloglog", "inverse", "inverse-squared"), 
                include_fs = TRUE, fs_by_spatial_BAU = FALSE,
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
  ## However, we still need a 1 for all other kinds; only produce a warning if 
  ## we have areal data.
  if (is.null(BAUs$wts)) {
    BAUs$wts <- 1
    if (any(sapply(data, function(x) is(x, "SpatialPolygons"))) &&
        !response %in% c("binomial", "negative-binomial")) # wts doesn't come into play for binomial or neg. binomial data, as it is forced to 1
      cat("SpatialPolygons were provided for the data support. No 'wts' field was found in the BAUs, so all BAUs are assumed to be of equal weight; if this is not the case, set the 'wts' field in the BAUs accordingly.\n")
  }
  
  ## When the response has a size parameter, restrict the incidence matrices 
  ## (Cz and Cp) to represent simple sums only. This behaviour is kept vague
  ## in the paper, but a note is provided to the user.
  ## Also enforce average_in_BAU = TRUE for simplicity.
  if (response %in% c("binomial", "negative-binomial")) {
    normalise_wts <- FALSE # Set to FALSE so that Cz represents an aggregation of the mean
    BAUs$wts <- 1
    average_in_BAU <- TRUE
    cat("You have selected a response distribution that has an associated size parameter. For simplicity, we enforce the elements of the incidence matrices (C_Z and C_P in the paper) to be 1, and average_in_BAU = TRUE.\n")
  }

  ## Check that the arguments are OK
  .check_args1(f = f, data = data, basis = basis, BAUs = BAUs, est_error = est_error, 
               response = response, link = link, K_type = K_type, 
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
    cat("Normalising basis function evaluations at BAU level...\n")
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
    cat("Binning data...\n")
    
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
  S    <- do.call("rbind",S)
  X    <- do.call("rbind",X)
  Cmat <- do.call("rbind",Cmat)
  Z    <- do.call("rbind",Z)
  Ve   <- do.call("bdiag",Ve)
  Vfs  <- do.call("bdiag",Vfs)
  
  ## Indices of observed BAUs
  obsidx <- .observed_BAUs_from_Cmat(Cmat)
  
  # Size parameter
  if(response %in% c("binomial", "negative-binomial")) {
    k_Z <- as.numeric(do.call("rbind", k_Z))
    
    ## Size parameter associated with observed BAUs.
    ## If any observations are associated with multiple BAUs, 
    ## we require the size parameter in a field of the BAUs;
    ## otherwise, we just use the observation size parameters.
    ## NB: THIS CODE ASSUMES THAT WE DO NOT HAVE OVERLAPPING DATA SUPPORTS 
    num_BAUs_each_data_support <- table(as(Cmat, "dgTMatrix")@i)
    if (!all(num_BAUs_each_data_support == 1)) {
      if (!("k_BAU" %in% names(BAUs))) 
        stop("When dealing with binomial or negative-binomial data, and some data supports are associated with multiple BAUs (e.g., areal data), the size parameter must be provided in the BAUs objects, in a field named 'k_BAU'.") 
      k_BAU_O <- BAUs$k_BAU[obsidx] 
    } else {
      ## Note that we have to re-order the observation size parameters, so that element i
      ## of k_Z is associated with the same BAU as element i of k_BAU_O;
      ## Cmat@i contains the index of the observation associated with BAU Cmat@j
      tmp <- as(Cmat, "dgTMatrix")
      k_BAU_O <- k_Z[tmp@i + 1]
      
      ## If all observations are associated with exactly one BAU, we can assign
      ## the data level size parameter (k_Z) to the BAU object. Note that it is 
      ## OK if non-observed BAUs do not have a size parameter. 
      BAUs@data[tmp@j + 1, "k_BAU"] <- k_BAU_O
    }

    if (any(is.na(k_BAU_O))) 
      stop("The size parameter is required at all observed BAUs")

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
    
    ## Throw a warning if the number of spatial BAUs (and hence number 
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
      D_basis = D_basis,
      Ve = Ve,
      Vfs = Vfs,
      Vfs_BAUs = Vfs_BAUs,
      Qfs_BAUs = Qfs_BAUs,
      Z = Z,
      Cmat = Cmat,
      X = X,
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
      mu_xi = l$mu_xi_init,
      k_Z = as.numeric(k_Z), 
      k_BAU_O = as.numeric(k_BAU_O), 
      include_fs = include_fs, 
      normalise_wts = normalise_wts, 
      fs_by_spatial_BAU = fs_by_spatial_BAU, 
      obsidx = obsidx)
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
  
  if (ncol(X) == 0) {
    stop("We need at least one covariate in the model")
  } else {
    l$alphahat_init <- solve(t(X) %*% X) %*% t(X) %*% Z 
  }
  l$sigma2fshat_init <- mean(diag(Ve)) / 4
  
  return(l)
}
