#' @rdname SRE
#' @export
FRK <- function(f,                     # formula (compulsory)
                data,                  # list of data objects (compulsory)
                basis = NULL,          # Basis object
                BAUs = NULL,           # BAUs
                est_error = TRUE,      # estimate measurement error
                average_in_BAU = TRUE, # average data into BAUs
                fs_model = "ind",      # fine-scale variation component
                vgm_model = NULL,      # variogram model for error estimation
                K_type = "block-exponential", # type of K matrix
                n_EM = 100,            # max. no. of EM iterations
                tol = 0.01,            # tolerance at which EM is assumed to have converged
                method = "EM",         # method for parameter estimation
                lambda = 0,            # regularisation parameter
                print_lik = FALSE,     # print log-likelihood at each iteration
                ...)                   # other arguments for BAUs/basis-function construction
{

    if(!is.list(data))          # Allow for user to supply data as a single object
        data <- list(data)      # If he/she does then put it into a list

    .check_args_wrapper(f = f,              # check that the arguments are OK for SRE
                        data = data,
                        basis = basis,
                        BAUs = BAUs,
                        est_error = est_error)
    .check_args2(...)                      # check that the arguments are OK for SRE.fit
    .check_args3(...)                      # check that the arguments are OK for prediction

    ## if there is a measurement error declared in all datasets then
    ## don't estimate it
    if(all(sapply(data,function(x) !is.null(x@data$std)))) {
        print("std already supplied with data -- not estimating the measurement error.
              If you wish to estimate measurement error then set the std field to NULL")
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

        print("Constructing BAUs...")

        BAUs <- auto_BAUs(manifold = manifold, # Construct BAUs
                          data = data[[d]],    # Using the dataset with largest extent
                          ...)
        BAUs$fs <- 1                           #Default fine-scale variation at BAU level
    } else {
        print("Assuming fine-scale variation is homoscedastic")
        if(is.null(BAUs$fs)) BAUs$fs <- 1      # If user supplied BAUs without fs field
                                               # then add on the default and inform user
    }

    if(is.null(basis)) {
        print("Generating basis functions...")
        tot_data <- sum(sapply(data,length))         # Total number of data points available
        if(K_type == "unstructured") {               # If unstructured then limit the
            max_sp_basis <- min(tot_data^(0.5),2000) # amount of basis functions to be sqrt
        } else {                                     # of data points (or 2000), else
            max_sp_basis <- 2000                     # just limit to 2000
            if(is(manifold,"sphere")) {              # If we're on the sphere just harc
                max_basis <- NULL                    # code the default basis functions
                nres <- 3                            # to 3 ISEA3h resolutions which is OK
                isea3h_lo <- 2                       # in most applications
            }
        }


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

    print(paste0("Modelling using ",nbasis(G)," basis functions"))
    print("Constructing SRE model...")

    S <- SRE(f = f,                            # formula
             data = data,                      # list of datasets
             basis = G,                        # basis functions
             BAUs = BAUs,                      # BAUs
             est_error=est_error,              # estimate measurement error?
             average_in_BAU = average_in_BAU,  # do not average data over BAUs
             fs_model = fs_model,              # fs model (only "ind" for now)
             vgm_model = vgm_model,            # vgm model for error estimation
             K_type = K_type)                  # "block-exponential" or "unstructured"

    ## After constructing SRE model, fit it
    print("Fitting SRE model...")
    S <- SRE.fit(SRE_model = S,       # SRE model
                 n_EM = n_EM,         # max. no. of EM iterations
                 tol = tol,           # tolerance at which EM is assumed to have converged
                 method = method,     # method (only "EM" for now)
                 lambda = lambda,     # regularisation parameter
                 print_lik=print_lik) # print log-likelihood at each iteration

    ## Return fitted SRE model
    S
}
