#' @rdname SRE
#' @export
FRK <- function(f,                     # formula
                data,                  # list of data items
                basis = NULL,          # Basis object
                BAUs = NULL,           # BAUs (SpatialPolygonsDataFrame)
                est_error=TRUE,       # estimate measurement error
                average_in_BAU = TRUE, # average data into BAUs
                fs_model = "ind",      # fine-scale variation component
                SRE_model = NULL,      # SRE model
                n_EM = 20,             # max. no. of EM iterations
                tol = 0.01,            # tolerance at which EM is assumed to have converged
                method = "EM",         # method for parameter estimation
                print_lik=TRUE,        # print log-likelihood at each iteration
                use_centroid = TRUE,   # use BAU centroid to evaluate basis functions over BAUs
                obs_fs = TRUE,         # where to include fs variation?
                pred_polys = NULL,     # spatial prediction polygons
                pred_time = NULL,      # temporal prediction points
                K_type = "block-exponential", # type of K matrix
                lambda = 0,              # regularisation parameter
                cross_validate = 1L,
                ...)
                {

    .check_args_wrapper(f=f,data=data,basis=basis,BAUs = BAUs,est_error = est_error)
    .check_args2(...)
    .check_args3(...)

    if(!is.list(data)) data <- list(data)

    manifold <- .choose_manifold_from_data(data[[1]])
    if(is.null(BAUs)) {
        print("Constructing BAUs...")

        BAUs <- auto_BAUs(manifold = manifold,
                          data = data[[1]],
                           ...)
        BAUs$fs <- 1   # fine-scale variation at BAU level
    }

    print("Generating basis functions...")
    G <- auto_basis(manifold =manifold,
                    data=data[[1]],...)

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

    print("Predicting...")
    BAUs_pred <- SRE.predict(SRE_model = S,          # SRE model
                             use_centroid = use_centroid,     # use centroid as point reference
                             obs_fs = obs_fs,
                             pred_polys = pred_polys,
                             pred_time = pred_time)

    list(S=S,PredPolys = BAUs_pred)

}
