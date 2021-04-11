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

#' @rdname SRE
#' @export
FRK <- function(f,                     # formula (compulsory)
                data,                  # list of data objects (compulsory)
                basis = NULL,          # Basis object
                BAUs = NULL,           # BAUs
                est_error = TRUE,      # estimate measurement error
                average_in_BAU = TRUE, # average data into BAUs
                sum_variables = NULL,  # variables to sum rather than average
                normalise_wts = TRUE,
                fs_model = "ind",      # fine-scale variation component
                vgm_model = NULL,      # variogram model for error estimation
                K_type = c("block-exponential", "precision", "unstructured"), # type of K matrix
                n_EM = 100,            # max. no. of EM iterations
                tol = 0.01,            # tolerance at which EM is assumed to have converged
                method = c("EM", "TMB"),         # method for parameter estimation
                lambda = 0,            # regularisation parameter
                print_lik = FALSE,     # print log-likelihood at each iteration
                response = c("gaussian", "poisson", "bernoulli", "gamma",
                             "inverse-gaussian", "negative-binomial", "binomial"), 
                link = c("identity", "log", "square-root", "logit", "probit", "cloglog", "inverse", "inverse-squared"),
                optimiser = nlminb,    # Optimiser for fitting (applicable only if method = 'TMB')
                percentiles = c(5, 95),  # Desired percentiles of the quantitity of interest
                fs_by_spatial_BAU = FALSE,
                known_sigma2fs = NULL, 
                taper = 4, 
                ...)                   # other arguments for BAUs/basis-function construction, or 
{

    if(!is.list(data))          # Allow for user to supply data as a single object
        data <- list(data)      # If he/she does then put it into a list

    
    ## Strings that must be lower-case (this allows users to enter 
    ## response = "Gaussian", for example, without causing issues)
    response  <- tolower(response)
    link      <- tolower(link)
    K_type    <- tolower(K_type)

    ## Match arguments/function
    K_type    <- match.arg(K_type)
    response  <- match.arg(response)
    link      <- match.arg(link)
    method    <- match.arg(method)
    optimiser <- match.fun(optimiser)
    
    ## FRK() is intended to be used for a novice user; simply enforce 
    ## method = "TMB" for non-Gaussian data or when a link other than the 
    ## identity is used, and enforce K_type = "precision" when method = "TMB"
    if (response != "gaussian" || link != "identity") {
        if (method != "TMB") 
            cat("Setting method = 'TMB', as a non-Gaussian data model or a non-identity link function has been chosen.\n")
        method <- "TMB"
    }
    if (method == "TMB") {
        if (K_type != "precision") 
            cat("For computational efficiency, setting K_type = 'precision', because method = 'TMB'; if you really want to use the K_type = 'block-exponential' with method = 'TMB', use SRE() and SRE.fit().\n")
        K_type <- "precision"
    }
    
    .check_args_wrapper(f = f,              # check that the arguments are OK for SRE
                        data = data,
                        basis = basis,
                        BAUs = BAUs,
                        est_error = est_error, 
                        response = response)
    
    # check that the arguments are OK for SRE.fit
    .check_args2(n_EM = n_EM, tol = tol, method = method, print_lik = print_lik, 
                 response = response, link = link, K_type = K_type, lambda = lambda,
                 optimiser = optimiser, fs_by_spatial_BAU = fs_by_spatial_BAU, 
                 known_sigma2fs = known_sigma2fs, BAUs = BAUs, ...)                      

    ## if there is a measurement error declared in all datasets then
    ## don't estimate it
    if(all(sapply(data,function(x) !is.null(x@data$std)))) {
        cat("std already supplied with data -- not estimating the measurement error.
              If you wish to estimate measurement error then set the std field to NULL\n")
        est_error <- FALSE
    }

    ## Attempt to automatically find the manifold from the data
    manifold <- .choose_manifold_from_data(data[[1]])

    ## Automatic BAU construction. First find the dataset enclosing the largest area
    ## by finding the area of the enclosing bounding box for each dataset
    ## We will also use this for basis-function construction
    d_areas <- sapply(data, function(d)
        prod(apply(coordinates(d),2,function(x) diff(range(x)))))
    d <- which.max(d_areas)


    ## Now construct the BAUs around this dataset
    if(is.null(BAUs)) {

        cat("Constructing BAUs...\n")

        BAUs <- auto_BAUs(manifold = manifold, # Construct BAUs
                          data = data[[d]],    # Using the dataset with largest extent
                          ...)
        BAUs$fs <- 1                           #Default fine-scale variation at BAU level
    } else {
        cat("Assuming fine-scale variation is homoscedastic\n")
        if(is.null(BAUs$fs)) BAUs$fs <- 1      # If user supplied BAUs without fs field
                                               # then add on the default and inform user
    }

    if(is.null(basis)) {
        cat("Generating basis functions...\n")
        tot_data <- sum(sapply(data,length))         # Total number of data points available
        if(K_type == "unstructured") {               # If unstructured then limit the
            max_sp_basis <- min(tot_data^(0.5),2000) # amount of basis functions to be sqrt
        } else {                                     # of data points (or 2000), else
            max_sp_basis <- 2000                     # just limit to 2000
            if(is(manifold,"sphere")) {              # If we're on the sphere just hard
                max_basis <- NULL                    # code the default basis functions
                nres <- 3                            # to 3 ISEA3h resolutions which is OK
                isea3h_lo <- 2                       # in most applications
            }
        }

        ## If nres is provided, don't need to use the default max number of 
        ## basis functions. Only consider nres if method = "TMB"; we don't want 
        ## to allow nres >= 4 if method = "EM", as it would be too slow. 
        if(!is.null(list(...)$nres) && method == "TMB") max_sp_basis <- NULL

        if(!(grepl("ST",class(manifold)))) {         # If we are NOT in a space-time scenario

            G <- auto_basis(manifold =manifold,      # Automatically generate basis functions using
                            data=data[[d]],          # data with largest spatal extent
                            ...,
                            max_basis = max_sp_basis) # max. number of basis functions

        ## However if we ARE in a space-time setting
        } else {
            ## Fix the number of temporal knots to 10 (reasonable default)
            ntime <- 10

            ## Construct the SPATIAL basis functions
            spatial_manifold <- strsplit(class(manifold),"ST")[[1]][2] # Find the spatial manifold
            if(identical(spatial_manifold,"plane")) { # Manifold is either plane or sphere
                spatial_manifold <- plane()
            } else { spatial_manifold <- sphere() }
            max_sp_basis <- max_sp_basis/ntime        # Maximum number of spatial basis functions
                                                      # is then max_sp_basis / ntime

            ## Construct the spatial basis functions by projecting the space-time
            ## data onto the spatial domain
            G_spatial <- auto_basis(manifold = spatial_manifold,
                                    data=as(data[[d]],"Spatial"),
                                    ...,
                                    max_basis = max_sp_basis)

            ## Construct temporal basis functions
            ## The end time point is equivalent to the time point
            ## of the last spatial BAUs
            endTimePt <- ncol(BAUs)
            time_knots <- seq(0,endTimePt,length=ntime)      # equally space knots in time
            time_scales <- rep(1.7 * endTimePt/ntime,        # default bisquares with a decent
                               length(time_knots))           # amount of overlap

            G_temporal <- local_basis(manifold=real_line(),     # on R^1
                                      loc = matrix(time_knots), # locations
                                      scale = time_scales,      # scales
                                      type = "bisquare")

            ## Construct spatio-temporal basis functions
            G <- TensorP(G_spatial,G_temporal)         # take the tensor product

        }
    } else {
        G <- basis  # If user has provided basis functions, just use these
    }

    cat("Modelling using",nbasis(G),"basis functions\n")
    cat("Constructing SRE model...\n")


    
    S <- SRE(f = f,                            # formula
             data = data,                      # list of datasets
             basis = G,                        # basis functions
             BAUs = BAUs,                      # BAUs
             est_error = est_error,            # estimate measurement error?
             average_in_BAU = average_in_BAU,  # do not average data over BAUs
             normalise_wts = normalise_wts,
             sum_variables = sum_variables,  
             fs_model = fs_model,              # fs model (only "ind" for now)
             vgm_model = vgm_model,            # vgm model for error estimation
             K_type = K_type,                  # "block-exponential", "unstructured", "precision"
             response = response, 
             link = link, 
             fs_by_spatial_BAU = fs_by_spatial_BAU, 
             taper = taper)                  

    ## After constructing SRE model, fit it
    cat("Fitting SRE model...\n")
    S <- SRE.fit(SRE_model = S,       # SRE model
                 n_EM = n_EM,         # max. no. of EM iterations
                 tol = tol,           # tolerance at which EM is assumed to have converged
                 method = method,     # method ("EM" or "TMB")
                 lambda = lambda,     # regularisation parameter
                 print_lik = print_lik, # print log-likelihood at each iteration
                 optimiser = optimiser, 
                 known_sigma2fs = known_sigma2fs, 
                 ...) 

    ## Return fitted SRE model
    return(S)
}


## The function below checks the arguments for the function FRK. The code is self-explanatory
## This is similar, but slightly different to, .check_args1()
.check_args_wrapper <- function(f,data,basis,BAUs,est_error,response) {
    if(!is(f,"formula")) stop("f needs to be a formula.")
    if(!is(data,"list"))
        stop("Please supply a list of Spatial objects.")
    if(!all(sapply(data,function(x) is(x,"Spatial") | is(x,"ST"))))
        stop("All data list elements need to be of class Spatial or ST.")
    if(!all(sapply(data,function(x) all.vars(f)[1] %in% names(x@data))))
        stop("All data list elements to have values for the dependent variable.")
    if(!est_error && response == "gaussian" & !all(sapply(data,function(x) "std" %in% names(x@data))))
        stop("If observational error is not going to be estimated,
             please supply a field 'std' in the data objects.")
    if(!(is.null(BAUs))) {
        if(!(is(BAUs,"SpatialPolygonsDataFrame") | is(BAUs,"SpatialPixelsDataFrame") | is(BAUs,"STFDF")))
            stop("BAUs should be a SpatialPolygonsDataFrame, SpatialPixelsDataFrame, or a STFDF object")
        if(!all(sapply(data,function(x) identical(proj4string(x), proj4string(BAUs)))))
            stop("Please ensure all data items and BAUs have the same coordinate reference system")
        if(!(all(BAUs$fs >= 0)))
            stop("fine-scale variation basis function needs to be nonnegative everywhere")
        if(is(BAUs,"STFDF")) if(!(is(BAUs@sp,"SpatialPolygonsDataFrame") | is(BAUs@sp,"SpatialPixelsDataFrame")))
            stop("The spatial component of the BAUs should be a SpatialPolygonsDataFrame")
        if(any(sapply(data,function(x) any(names(x@data) %in% names(BAUs@data)))))
            stop("Please don't have overlapping variable names in data and BAUs. All covariates need to be in the BAUs.")
        if(!all(all.vars(f)[-1] %in% c(names(BAUs@data),coordnames(BAUs))))
            stop("All covariates need to be in the SpatialPolygons BAU object.")
    }

    if(!(is.null(basis))) {
        if(!(is(basis,"Basis") | is(basis,"TensorP_Basis")))
            stop("basis needs to be of class Basis  or TensorP_Basis (package FRK)")
    }
}

