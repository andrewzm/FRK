

## Checks arguments for the SRE() function. Code is self-explanatory
.check_args1 <- function(f, data,basis, BAUs, est_error, 
                         K_type, response, link, fs_by_spatial_BAU, normalise_wts, 
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
  if(!all(sapply(data,function(x) coordnames(x) == coordnames(BAUs))))
    stop("The coordinate names of data do not match those of the BAUs.")
  
  
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
  
  if(response == "gaussian" && !est_error && !all(sapply(data, function(x) x$std >= 0)))
    stop("If the response is Gaussian and observational error is not going to be estimated,
             the std field must contain only positive numbers")
  
  
  if(!(K_type %in% c("block-exponential", "precision", "unstructured")))
    stop("Invalid K_type argument. Please select from 'block-exponential', 'precision', 'unstructured', or 'separable'")
  if(K_type == "unstructured" & response != "gaussian")
    stop("The unstructured covariance matrix (K_type = 'unstructured') is not implemented for non-Gaussian response (more specifically, when method = 'TMB')")
  if (K_type == "block-exponential" & response != "gaussian")
    warning("Using the block-exponential covariance matrix (K_type = 'block-exponential') is computationally inefficient with a non-Gaussian response (more specifically, when method = 'TMB'). For these situations, consider using K_type = 'precision'.")
  
  # if (K_type == "separable" & is(basis,"TensorP_Basis"))
  # stop("K_type = 'separable' is not yet implemented in a space-time setting.")
  
  # ## If K_type == separable, the basis functions must be in a regular rectangular lattice
  # if (K_type == "separable") 
  #   for (i in unique(basis@df$res)) {
  #     temp <- basis@df[basis@df$res == i, ]
  #     if (!.test_regular_grid(temp$loc1, temp$loc2, rectangular = TRUE) ) 
  #       stop("Basis functions must be arranged in a regular rectangular lattice when K_type = 'separable'.")
  #   }
  
  # ## If K_type == "neighbour", we just need basis functions to be in a regular 
  # ## lattice (does not need to be rectangular)
  #   if (K_type == "neighbour") 
  #     for (i in unique(basis@df$res)) {
  #       temp <- basis@df[basis@df$res == i, ]
  #       if (!.test_regular_grid(temp$loc1, temp$loc2, rectangular = FALSE) ) 
  #         stop("Basis functions must be arranged in a regular lattice when K_type = 'precision'.")
  #     }  
  
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
                         BAUs, taper, ...) {
  
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
  if(method == "EM" & K_type == "precision") stop("The precision matrix formulation of the model is not implemented for method = 'EM'. Please choose K_type to be 'block-exponential' or 'unstructured'.")
  if(method == "EM" & K_type == "separable") stop("The separable spatial model is not implemented for method = 'EM'. Please choose K_type to be 'block-exponential' or 'unstructured'.")
  if(method == "TMB" & K_type == "unstructured") stop("The unstructured covariance matrix (K_type = 'unstructured') is not implemented for method = 'TMB'")
  if(method != "TMB" & fs_by_spatial_BAU) stop("fs_by_spatial_BAU can only be TRUE if method = 'TMB'. Please set method = 'TMB', or fs_by_spatial_BAU = FALSE.")
  
  if(print_lik && method == "TMB")
    cat("The likelihood at each iteration cannot be accessed because you have selected TMB for model fitting: print_lik will be ignored.")
  
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
  
  # SRE_model is present only because of a deprecation coercion 
  l$SRE_model <- NULL
  
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
  
  ## Check taper
  if (!is.null(taper)) {
    if (!(class(taper) %in% c("numeric", "integer"))) {
      stop("taper, the argument controlling the coveriance taper, must be numeric or integer.")
    } else if (taper <= 0) {
      stop("taper, the argument controlling the coveriance taper, must be positive.")
    }
  }

}

## Checks arguments for the predict() function. Code is self-explanatory
.check_args3 <- function(obs_fs, newdata, pred_polys,
                         pred_time, covariances, object, type, 
                         k, percentiles, kriging, ...) {
  if(!is.null(newdata)) {
    if(!all(coordnames(newdata) == coordnames(object@BAUs)))
      stop("The coordinate names of newdata do not match those of the BAUs.")
  }

  if(kriging != "simple" & object@method == "EM")
    stop("Universal kriging is only available when method = 'TMB'")
  
  if(!(obs_fs %in% 0:1)) stop("obs_fs needs to be logical")
  
  if(!(is(newdata,"Spatial") | is(newdata,"ST") | is.null(newdata)))
    stop("Predictions need to be over Spatial or ST objects")
  
  if(!is.null(newdata) & !is.null(pred_time))
    stop("Only one of newdata and pred_time can be not NULL")
  
  if(!(is.integer(pred_time) | is.null(pred_time))) stop("pred_time needs to be of class integer")
  if(!is.logical(covariances)) stop("covariances needs to be TRUE or FALSE")
  
  ## Quantities of interest
  if(!all(type %in% c("link", "mean", "response")))
    stop("type must be a vector containing combinations of 'link', 'mean', and 'response'")
  
  ## Check k (for predictions)
  if (object@response %in% c("binomial", "negative-binomial")) {
    if(length(k) == 1){
      warning("Single number k provided for all BAUs: assuming k is invariant over the whole spatial domain.")
    } else if (!(class(k) %in% c("numeric", "integer"))) {
      stop("k must contain only positive integers.")
    } else if (any(k < 0) | any(k != round (k))) {
      stop("k must contain only positive integers.")
    } else if (length(k) != nrow(object@S0)) {
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
}


