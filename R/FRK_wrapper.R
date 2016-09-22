#' @rdname SRE
#' @export
FRK <- function(f,                     # formula
                data,                  # list of data items
                basis = NULL,          # Basis object
                BAUs = NULL,           # BAUs (SpatialPolygonsDataFrame)
                K_type = "block-exponential", # type of K matrix
                lambda = 0,              # regularisation parameter
                fs_model = "ind",      # fine-scale variation component
                average_in_BAU = TRUE, # average data into BAUs
                est_error=TRUE,       # estimate measurement error
                n_EM = 100,             # max. no. of EM iterations
                tol = 0.01,            # tolerance at which EM is assumed to have converged
                method = "EM",         # method for parameter estimation
                print_lik=TRUE,        # print log-likelihood at each iteration
                cross_validate = 1L,
                ...)
{

    if(!is.list(data)) data <- list(data)
    .check_args_wrapper(f=f,data=data,basis=basis,BAUs = BAUs,est_error = est_error)
    .check_args2(...)
    .check_args3(...)

    if(all(sapply(data,function(x) !is.null(x@data$std)))) {
        print("std already supplied with data -- not estimating the measurement error")
        est_error <- FALSE
    }

    manifold <- .choose_manifold_from_data(data[[1]])
    if(is.null(BAUs)) {
        print("Constructing BAUs...")
        ## Find dataset enclosing largest box area on which to construct BAUs
        d <- which.max(sapply(data, function(d)
            prod(apply(coordinates(d),2,function(x) diff(range(x))))))
        BAUs <- auto_BAUs(manifold = manifold,
                          data = data[[d]],
                          ...)
        BAUs$fs <- 1   # fine-scale variation at BAU level
    }


    print("Generating basis functions...")
    tot_data <- length(data[[1]])
    if(K_type == "unstructured") {
        max_sp_basis <- min(tot_data^(0.5),2000)
    } else {
        max_sp_basis <- 2000
        if(is(manifold,"sphere")) {
            max_basis <- NULL
            nres <- 3
            isea3h_lo <- 2
        }
    }


    if(!(grepl("ST",class(manifold)))) {

        G <- auto_basis(manifold =manifold,
                        data=data[[1]],...,max_basis = max_sp_basis)
    } else {
        ## Construct spatial basis functions
        spatial_manifold <- strsplit(class(manifold),"ST")[[1]][2]
        if(identical(spatial_manifold,"plane")) spatial_manifold <- plane() else spatial_manifold <- sphere()
        max_sp_basis <- max_sp_basis/10

        G_spatial <- auto_basis(manifold = spatial_manifold,
                                data=as(data[[d]],"Spatial"),  # flatten data
                                ...,
                                max_basis = max_sp_basis)

        ## Construct temporal basis functions
        #ntime <- min(length(unique(time(data[[d]]))),10)
        ntime <- 10
        time_knots <- seq(0,ncol(BAUs)+1,length=ntime)
        G_temporal <- local_basis(manifold=real_line(),         # on R^1
                                  loc = matrix(time_knots),       # locations
                                  scale = rep(ncol(BAUs)/ntime/1.2,length(time_knots)))           # scales

        ## Construct spatio-temporal basis functions
        G <- TensorP(G_spatial,G_temporal)         # take the tensor product

    }

    print(paste0("Modelling using ",nbasis(G)," basis functions"))
    print("Constructing SRE model...")
    S <- SRE(f = f,                  # formula
             data = data,     # list of datasets
             BAUs = BAUs,       # BAUs
             basis = G,              # basis functions
             est_error=est_error,         # estimation measurement error
             average_in_BAU = average_in_BAU, # do not average data over BAUs
             fs_model = fs_model,
             K_type = K_type)

    print("Fitting SRE model...")
    S <- SRE.fit(SRE_model = S,    # SRE model
                 n_EM = n_EM,        # max. no. of EM iterations
                 tol = tol,       # tolerance at which EM is assumed to have converged
                 print_lik=print_lik,   # print log-likelihood at each iteration
                 cross_validate = cross_validate,
                 lambda = lambda)

    S
}
