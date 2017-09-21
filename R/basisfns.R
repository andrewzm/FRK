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

#' @title Construct a set of local basis functions
#' @description Construct a set of local basis functions based on pre-specified location and scale parameters.
#' @param manifold object of class \code{manifold}, for example, \code{sphere}
#' @param loc a matrix of size \code{n} by \code{dimensions(manifold)} indicating centres of basis functions
#' @param scale vector of length \code{n} containing the scale parameters of the basis functions; see details
#' @param type either ``bisquare'', ``Gaussian'', ``exp'', or ``Matern32''
#' @details This functions lays out local basis functions in a domain of interest based on pre-specified location and scale parameters. If \code{type} is ``bisquare'', then
#'\deqn{\phi(u) = \left(1- \left(\frac{\| u \|}{R}\right)^2\right)^2 I(\|u\| < R),}
#' and \code{scale} is given by \eqn{R}, the range of support of the bisquare function. If \code{type} is ``Gaussian'', then
#' \deqn{\phi(u) = \exp\left(-\frac{\|u \|^2}{2\sigma^2}\right),}
#' and \code{scale} is given by \eqn{\sigma}, the standard deviation. If \code{type} is ``exp'', then
#'\deqn{\phi(u) = \exp\left(-\frac{\|u\|}{ \tau}\right),}
#' and \code{scale} is given by \eqn{\tau}, the e-folding length. If \code{type} is ``Matern32'', then
#'\deqn{\phi(u) = \left(1 + \frac{\sqrt{3}\|u\|}{\kappa}\right)\exp\left(-\frac{\sqrt{3}\| u \|}{\kappa}\right),}
#' and \code{scale} is given by \eqn{\kappa}, the function's scale.
#' @examples
#' library(ggplot2)
#' G <-  local_basis(manifold = real_line(),
#'                    loc=matrix(1:10,10,1),
#'                    scale=rep(2,10),
#'                    type="bisquare")
#' \dontrun{show_basis(G)}
#' @export
local_basis <- function(manifold=sphere(),          # default manifold is sphere
                        loc=matrix(c(1,0),nrow=1),  # one centroid at (1,0)
                        scale=1,                    # std = 1, and Gaussian RBF
                        type=c("bisquare","Gaussian","exp","Matern32")) {

    ## Basic checks
    if(!is.matrix(loc)) stop("loc needs to be a matrix")
    if(!(dimensions(manifold) == ncol(loc))) stop("number of columns in loc needs to be the
                                                  same as the number of manifold dimensions")
    if(!(length(scale) == nrow(loc))) stop("need to have as many scale parameters as centroids")
    type <- match.arg(type)

    n <- nrow(loc)                                # n is the number of centroids
    colnames(loc) <- c(outer("loc",1:ncol(loc),   # label dimensions as loc1,loc2,...,locn
                             FUN = paste0))

    fn <- pars <- list()  # initialise lists of functions and parameters
    for (i in 1:n) {
        ## Assign the functions as determined by user
        if(type=="bisquare") {
            fn[[i]] <-  .bisquare_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="Gaussian") {
            fn[[i]] <-  .GRBF_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="exp") {
            fn[[i]] <-  .exp_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        } else if (type=="Matern32") {
            fn[[i]] <-  .Matern32_wrapper(manifold,matrix(loc[i,],nrow=1),scale[i])
        }

        ## Save the parameters to the parameter list (loc and scale)
        pars[[i]] <- list(loc = matrix(loc[i,],nrow=1), scale=scale[i])
    }

    ## Create a data frame which summarises info about the functions. Set resolution = 1
    df <- data.frame(loc,scale,res=1)

    ## Create new basis function, using the manifold, n, functions, parameters list, and data frame.
    this_basis <- new("Basis", manifold=manifold,  n=n, fn=fn, pars=pars, df=df)
    return(this_basis)
}

#' @title Automatic basis-function placement
#' @description Generate automatically a set of local basis functions in the domain, and automatically prune in regions of sparse data.
#' @param manifold object of class \code{manifold}, for example, \code{sphere} or \code{plane}
#' @param data object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame} containing the data on which basis-function placement is based, or a list of these; see details
#' @param regular an integer indicating the number of regularly-placed basis functions at the first resolution. In two dimensions, this dictates the smallest number of basis functions in a row or column at the coarsest resolution. If \code{regular=0}, an irregular grid is used, one that is based on the triangulation of the domain with increased mesh density in areas of high data density; see details
#' @param nres the number of basis-function resolutions to use
#' @param prune a threshold parameter that dictates when a basis function is considered irrelevent or unidentifiable, and thus removed; see details
#' @param max_basis maximum number of basis functions. This overrides the parameter \code{nres}
#' @param subsamp the maximum amount of data points to consider when carrying out basis-function placement: these data objects are randomly sampled from the full dataset. Keep this number fairly high (on the order of 10^5), otherwise fine-resolution basis functions may be spuriously removed
#' @param type the type of basis functions to use; see details
#' @param isea3h_lo if \code{manifold = sphere()}, this argument dictates which ISEA3H resolution is the coarsest one that should be used for the first resolution
#' @param scale_aperture the aperture (in the case of the bisquare, but similar interpretation for other basis) width of the basis function is the minimum distance between all the basis function centroids multiplied by \code{scale_aperture}. Typically this ranges between 1 and 1.5 and is defaulted to 1 on the sphere and 1.25 on the other manifolds.
#' @param bndary a \code{matrix} containing points containing the boundary. If \code{regular == 0} this can be used to define a boundary in which irregularly-spaced basis functions are placed
#' @param verbose a logical variable indicating whether to output a summary of the basis functions created or not
#' @param ... unused
#' @details This function automatically places basis functions within the domain of interest. If the domain is a plane or the real line, then the object \code{data} is used to establish the domain boundary.
#'
#'
#'The argument \code{type} can be either ``Gaussian'', in which case
#'\deqn{\phi(u) = \exp\left(-\frac{\|u \|^2}{2\sigma^2}\right),}
#'``bisquare'', in which case
#'\deqn{\phi(u) = \left(1- \left(\frac{\| u \|}{R}\right)^2\right)^2 I(\|u\| < R),}
#'``exp'', in which case
#'\deqn{\phi(u) = \exp\left(-\frac{\|u\|}{ \tau}\right),}
#' or ``Matern32'', in which case
#'\deqn{\phi(u) = \left(1 + \frac{\sqrt{3}\|u\|}{\kappa}\right)\exp\left(-\frac{\sqrt{3}\| u \|}{\kappa}\right),}
#' where the parameters \eqn{\sigma, R, \tau} and \eqn{\kappa} are \code{scale} arguments.
#'
#' If the manifold is the real line, the basis functions are placed regularly inside the domain, and the number of basis functions at the coarsest resolution is dictated by the integer parameter \code{regular} which has to be greater than zero. On the real line, each subsequent resolution has twice as many basis functions. The scale of the basis function is set based on the minimum distance between the centre locations following placement. The scale is equal to the minimum distance if the type of basis function is Gaussian, exponential, or Matern32, and is equal to 1.5 times this value if the function is bisquare.
#'
#' If the manifold is a plane, and \code{regular > 0}, then basis functions are placed regularly within the bounding box of \code{data}, with the smallest number of basis functions in each row or column equal to the value of \code{regular} in the coarsest resolution (note, this is just the smallest number of basis functions). Subsequent resolutions have twice the number of basis functions in each row or column. If \code{regular = 0}, then the function \code{INLA::inla.nonconvex.hull} is used to construct a (non-convex) hull around the data. The buffer and smoothness of the hull is determined by the parameter \code{convex}. Once the domain boundary is found,  \code{INLA::inla.mesh.2d} is used to construct a triangular mesh such that the node vertices coincide with data locations, subject to some minimum and maximum triangular-side-length constraints. The result is a mesh that is dense in regions of high data density and not dense in regions of sparse data. Even basis functions are irregularly placed, the scale is taken to be a function of the minimum distance between basis function centres, as detailed above. This may be changed in a future revision of the package.
#'
#' If the manifold is the surface of a sphere, then basis functions are placed on the centroids of the discrete global grid (DGG), with the first basis resolution corresponding to the third resolution of the DGG (ISEA3H resolution 2, which yields 92 basis functions globally).  It is not recommended to go above \code{nres == 3} (ISEA3H resolutions 2--4) for the whole sphere; \code{nres=3} yields a total of 1176 basis functions. Up to ISEA3H resolution 6 is available with \code{FRK}; for finer resolutions; please install \code{dggrids} from \code{https://github.com/andrewzm/dggrids} using \code{devtools}.
#'
#' Basis functions that are not influenced by data points may hinder convergence of the EM algorithm when \code{K_type = ``unstructured''}, since the associated hidden states are, by and large, unidentifiable. We hence provide a means to automatically remove such basis functions through the parameter \code{prune}. The final set only contains basis functions for which the column sums in the associated matrix \eqn{S} (which, recall, is the value/average of the basis functions at/over the data points/polygons) is greater than \code{prune}. If \code{prune == 0}, no basis functions are removed from the original design.
#' @examples
#' \dontrun{
#' library(sp)
#' library(ggplot2)
#'
#' ### Create a synthetic dataset
#' set.seed(1)
#' d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
#'                 lat = runif(n=1000,min = -90, max = 90),
#'                 z = rnorm(5000))
#' coordinates(d) <- ~lon + lat
#' proj4string(d)=CRS("+proj=longlat +ellps=sphere")
#'
#' ### Now create basis functions over sphere
#' G <- auto_basis(manifold = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#'
#' ### Plot
#' \dontrun{show_basis(G,draw_world())}
#' }
#' @export
auto_basis <- function(manifold = plane(),
                       data,
                       regular = 1,
                       nres = 3,
                       prune = 0,
                       max_basis = NULL,
                       subsamp = 10000,
                       type = c("bisquare","Gaussian","exp","Matern32"),
                       isea3h_lo = 2,
                       bndary = NULL,
                       scale_aperture = ifelse(is(manifold,"sphere"),1,1.25),
                       verbose = 0L,
                       ...) {      # currently unused

    ## Basic checks
    type <- match.arg(type)
    if(!is(manifold,"manifold"))
        stop("manifold needs to be an object of class manifold")
    if(!is.numeric(prune) | prune < 0)
        stop("prune needs to be greater than zero")
    if(!is.numeric(subsamp) | subsamp < 0)
        stop("subsamp needs to be greater than zero")
    if((is(manifold,"sphere")  | is(manifold,"real_line")) & regular == 0)
        stop("Irregular basis only available on planes")
    if(!(is(isea3h_lo,"numeric")))
        stop("isea3h_lo needs to be an integer greater than 0")
    if(!(is.numeric(nres) | is.null(nres)))
        stop("nres needs to be greater than zero or NULL")

    ## If the user has specified a maximum number of basis functions then
    ## we need to add resolutions iteratively and stop when exceed the maximum
    ## of basis functions
    if(!is.null(max_basis)) {
        cat("...Automatically choosing number of functions...\n")
        tot_basis <- 0                  # start of with 0 basis functions
        tot_data <- length(data)        # number of data points
        nres <- 1                       # start off with one resolution (we have a minimum of one resolution)
        while(tot_basis <= max_basis) { # while we have less basis than the max
            nres <- nres + 1            # incremement the res number
            G <- .auto_basis(manifold = manifold,  # arguments as described above
                             data = data,
                             prune = prune,
                             regular = regular,
                             nres = nres,
                             subsamp = subsamp,
                             type = type,
                             isea3h_lo = isea3h_lo,
                             bndary = bndary,
                             scale_aperture = scale_aperture,
                             verbose = 0)  # force verbose to 0 for this procedure of finding the
            # number of basis functions
            tot_basis <- nbasis(G)      # record the number of basis functions
        }
        nres <- nres - 1  # nres resolutions was too much, deduct by one
    }

    ## Now call the local function with checked parameters
    .auto_basis(manifold=manifold,data=data,regular=regular,nres=nres,
                prune=prune,subsamp=subsamp,type=type,isea3h_lo = isea3h_lo,
                bndary=bndary, scale_aperture = scale_aperture, verbose=verbose)
}


.auto_basis <- function(manifold = plane(),
                        data,
                        regular=1,
                        nres=3,
                        prune=0,
                        subsamp=10000,
                        type = c("bisquare","Gaussian","exp","Matern32"),
                        isea3h_lo = 2,
                        bndary = NULL,
                        scale_aperture = ifelse(is(manifold,"sphere"),1,1.25),
                        verbose = 0L) {

    type <- match.arg(type)            # match type
    m <- manifold                      # abbreviate for convenience
    isea3h <- centroid <- res <- NULL  # (suppress warnings, these are loaded from data)
    coords <- coordinates(data)        # data coordinates or centroids

    ## Irregular basis function placement can only proceed with INLA. Throw an error
    ## if INLA is not installed
    if(is(m,"plane") & regular == 0 & is.null(bndary)) {
        if(!requireNamespace("INLA"))
            stop("For irregularly-placed basis-function generation INLA needs to be installed
                 for constructing basis function centres. Please install it
                 using install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\")")
    }

    ## Subsamp can be used to make the placement/pruning process more efficient with big data.
    ## If subsamp  < number of data points then sample the number of data randomly.
    if(nrow(coords)>subsamp) {
        coords <- coords[sample(nrow(coords),size=subsamp,replace=FALSE),]
    }

    ## Find the x and y extent of the (subsampled) data
    xrange <- range(coords[,1])
    yrange <- range(coords[,2])

    ## If we are on the plane and want an irregular function placement, then call INLA
    ## and find a nonconvex hull in which the basis functions can be enclosed. If
    ## a boundary is supplied as a matrix, convert it an inla.mesh.segment object.
    if(is(m,"plane") & regular == 0) {
        if(is.null(bndary)) {
            bndary_seg = INLA::inla.nonconvex.hull(coords,concave = 0)
        } else {
            if(!is(bndary,"matrix"))
                stop("bndary needs to be a matrix of points")
            bndary_seg <- INLA::inla.mesh.segment(bndary)
        }
    }

    ## If we want basis functions regularly placed on the plane then
    ## first calculate the aspect ratio (extent of y / extent of x)
    ## Then set the number of basis functions in y = regular if
    ## the aspect ratio < 1, or number of basis functions in x = regular
    ## if the aspect ratio >= 1. The number of basis functions in the other
    ## axis is then amplified by 1/aspect ratio, or aspect ratio, respectively
    if(is(m,"plane") & regular > 0) {
        asp_ratio <- diff(yrange) / diff(xrange)
        if(asp_ratio < 1) {
            ny <- regular
            nx <- ny / asp_ratio
        } else {
            nx <- regular
            ny <- nx * asp_ratio
        }
    }

    ## Initialise the variables loc and scale which we will fill in later (centroids and scales)
    loc <- scale <- NULL

    ## Initialise G which will contain a list of basis function objects by resolution (1 resolution per list item)
    G <- list()

    ## If we are on the sphere, then load the dggrids. This will be either
    ## from FRK if isea3h_lo + nres - 1 <= 6, or else from the dggrids package if
    ## sea3h_lo + nres - 1 > 6. See load_dggrids() for details
    if(is(m,"sphere")) {
        isea3h <- load_dggrids(res = isea3h_lo + nres - 1) %>%
            dplyr::filter(res >= isea3h_lo) # keep only the resolutions we need
    }


    ## Find possible grid points for basis locations in each resolution
    for(i in 1:nres) {

        ## If we are on the plane and we do not want a regular set
        if(is(m,"plane") & (regular == 0)) {
            ## Generate mesh using INLA and use the triangle vertices as notes
            ## The maximum and minimum triangle edges are a function of i
            ## The exact constants were found to be suitable by trial and error
            this_res_locs <- INLA::inla.mesh.2d(
                loc = matrix(apply(coords,2,mean),nrow=1),   # data locations
                boundary = list(bndary_seg),                 # boundary
                max.edge = max(diff(xrange),diff(yrange))/(2*2.5^(i-1)),
                cutoff = max(diff(xrange),diff(yrange))/(3*2.5^(i-1)))$loc[,1:2]

            ## Sometimes INLA returns overlapping points, therefore find the unique set of locations
            this_res_locs <- unique(this_res_locs)

            ## If we want a regular set of basis functions
        } else if(is(m,"plane") & (regular > 0)) {
            ## Make a grid that gets more dense as i increases
            xgrid <- seq(xrange[1], xrange[2], length = round(nx*(3^(i)))) # x coordinates of centroids
            ygrid <- seq(yrange[1], yrange[2], length = round(ny*(3^(i)))) # y coordinates of centroids
            this_res_locs <- xgrid %>%
                expand.grid(ygrid) %>%   # form the grid in long-table format
                as.matrix()              # convert to matrix
            ## If we are on the real line
        } else if(is(m,"real_line")) {
            ## Simply set the centroids equally spaced on the line
            this_res_locs <- matrix(seq(xrange[1],xrange[2],length=i*regular))
            ## If we are on the sphere
        } else if(is(m,"sphere")) {
            ## Simply take the centroids of the ISEA3H polygons at the appropriate resolutions to be the centroids
            this_res_locs <- filter(isea3h,centroid==1,res==(i + isea3h_lo - 1))[c("lon","lat")] %>%
                as.matrix() # and convert to matrix
        }

        ## Setting the scales and Refinement stage
        ##      Set scales: e.g., set to 1.5x the distance to nearest basis if bisquare
        ##      Refine/Prune: Remove basis which are not influenced by data and re-find the scales
        for(j in 1:2) { ## Need to go over this twice for refinement
            D <- FRK::distance(m,this_res_locs,this_res_locs)       # compute distance
            if(nrow(D) == 1) {                                      # if we only have one basis function
                this_res_scales <-max(diff(xrange),diff(yrange))/2  # set the "base" scale parameter to span most of the range
            } else {
                ## otherwise set the "base" scale parameter to be based
                ## on the distance to the nearest centroid (if not on sphere)
                diag(D) <- Inf
                this_res_scales <- apply(D,1,min)*scale_aperture
            }

            ## If we have more than one basis at this resolution (could be 0 because of pruning)
            if(nrow(D) >0)
                ## The following code can be used to verify that the following basis functions have a similar shape:
                ## Bisquare: R = 1.5,
                ## Gaussian: sigma = 0.7,
                ## Exp: tau = 0.7,
                ## Matern32: kappa = 0.7
                # f1 <- .bisquare_wrapper(real_line(),c = matrix(0),1.5)
                # f2 <- .GRBF_wrapper(real_line(),mu = matrix(0),0.7)
                # f3 <- .exp_wrapper(real_line(),c = matrix(0),tau = 0.7)
                # f4 <- .Matern32_wrapper(real_line(),c = matrix(0),kappa = 0.7)
                # x <- seq(-4,4,length=100)
                # plot(x,f1(x)); lines(x,f2(x)); lines(x,f3(x)); lines(x,f4(x))
            ## This explains the choice of scaling below
            this_res_basis <- local_basis(manifold = m,
                                          loc=this_res_locs,
                                          scale=ifelse(type=="bisquare",1.5,0.7)*this_res_scales,
                                          type=type)

            ## Now refine/prune these basis functions
            if(prune > 0 & nrow(D)>0) {
                ## Only in the first step
                if(j==1) {
                    ## The basis functions to remove are those which are not "considerably affected"
                    ## by observations. We determine this influence by evaluating the basis functions
                    ## at the data locations and summing over the function values. If the sum is less
                    ## then 'prune' we remove the basis function
                    rm_idx <- which(colSums(eval_basis(this_res_basis,coords)) < prune)

                    ## Throw a warning if all basis functions at a given resolution have been removed
                    if(length(rm_idx) == length(this_res_scales))
                        warning("prune is too large -- all functions at a resolution removed.
                                Consider also removing number of resolutions.")

                    ## If there are basis functions to remove,then remove them and reset scales in j==2
                    if(length(rm_idx) > 0)
                        this_res_locs <- this_res_locs[-rm_idx,,drop=FALSE]
                }
            } else {
                ## If we are not pruning just stop after j = 1 (no need to refind scales twice)
                break
            }
        }

        ## Print the number of basis functions at each resolution if verbose
        if(verbose)
            cat("Number of basis at resolution",i,"=",nrow(this_res_locs))

        ## Now that we actually have the centroids and scales we can construct the basis functions for this
        ## resolution. If all the basis functions have been removed at this resolution don't do anything
        if(nrow(D) > 0) {
            G[[i]] <-  local_basis(manifold = m,
                                   loc=this_res_locs,
                                   scale=ifelse(type=="bisquare",1.5,0.7)*this_res_scales,
                                   type=type)
            G[[i]]$res=i  # put resolution in the basis function data frame
        }
    }

    ## Finally concatenate all the basis functions together using the S4 method concat
    G_basis <- Reduce("concat",G)
    G_basis
}


#' @rdname TensorP
#' @aliases TensorP,Basis-Basis-method
setMethod("TensorP",signature(Basis1="Basis",Basis2="Basis"),function(Basis1,Basis2) {
    ## This function constructs the Tensor product of two sets of basis functions
    ## Currently there is the restriction that only one can be multiresolution
    ## (This is typically the case, as time usually has one resolution of basis functions in FRK)
    if(nres(Basis2) > 1)
        stop("Only one basis can be multiresolution (Basis 1)")

    ## Extract data frames and dimensions
    df1 <- data.frame(Basis1)
    df2 <- data.frame(Basis2)
    n1 <- dimensions(Basis1)
    n2 <- dimensions(Basis2)

    ## Create long data frame with all possible centroid combinations
    ## (Kronecker of space and time)
    expand.grid(df1[,1:n1,drop=FALSE],
                df2[,1:n2,drop=FALSE])

    ## We adopt a space-first approach, where the spatial index changes quicker than the
    ## temporal index. So we repeat the centroids of Basis1 for as many time points as
    ## we have, and we repeat the time points so they look like 1111111,222222,33333 etc.
    ## The resolution of the basis function is inherited from the spatial resolution
    ## NB: Only one resolution for Basis2 is assumed. This is checked aboce.
    df <- cbind(df1[rep(1:nrow(df1),times = nrow(df2)),1:n1,drop=FALSE],
                df2[rep(1:nrow(df2),each = nrow(df1)),1:n2,drop=FALSE],
                df1[rep(1:nrow(df1),times = nrow(df2)),"res",drop=FALSE])

    ## Change the names of the data frame to what is standard in this package
    names(df) <- c(paste0("loc",1:(n1 + n2)),"res")

    ## Create new Tensor Basis function from this information
    new("TensorP_Basis",
        Basis1=Basis1,
        Basis2=Basis2,
        n = nbasis(Basis1) * nbasis(Basis2),
        df = df)
})


#' @rdname eval_basis
#' @aliases eval_basis,Basis-matrix-method
setMethod("eval_basis",signature(basis="Basis",s="matrix"),
          function(basis,s){

              space_dim <- dimensions(basis)  # spatial dimensions
              n <- nrow(s)                    # number of points over which to evaluate basis functions

              batching=cut(1:n,                      # create batches of computation (in batches of 10000)
                           breaks = seq(0,           # break up into batches ot 10000
                                        n+10000,
                                        by=10000),
                           labels=F)                 # do not assign labels to batches

              if(opts_FRK$get("parallel") > 1L) {    # if we have a parallel backend when compute in parallel

                  clusterExport(opts_FRK$get("cl"),  # export variables to the cluster
                                c("batching","basis","s","space_dim"),
                                envir=environment())

                  ## Use parLapply to compute batches in parallel. The drop = FALSE in the end is required
                  ## for when we have one spatial dimension
                  pnt_eval_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                                             function(i) {
                                                 idx <- which(batching == i)     # see which spatial locations to evaluate on
                                                 return(.point_eval_fn(basis@fn, # evaluate on these locations
                                                                       s[idx,1:space_dim,drop = FALSE]))
                                             })
                  clusterEvalQ(opts_FRK$get("cl"), {gc()})  # clear data from the cluster workers
              } else  {
                  ## Same as above but using lapply
                  pnt_eval_list <- lapply(1:max(unique(batching)),
                                          function(i) {
                                              idx <- which(batching == i)
                                              return(.point_eval_fn(basis@fn,
                                                                    s[idx,1:space_dim,drop = FALSE]))
                                          })
              }

              ## Finally concatenate all the bits together using rBind
              do.call(rBind,pnt_eval_list)
          })

#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPointsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPointsDataFrame"),
          function(basis,s){

              ## Now just evaluate the basis functions at the coordinates of the SpatialPoints
              eval_basis(basis=basis,
                         s = coordinates(s))
          })


#' @rdname eval_basis
#' @aliases eval_basis,Basis-SpatialPolygonsDataFrame-method
setMethod("eval_basis",signature(basis="Basis",s="SpatialPolygonsDataFrame"),
          function(basis,s){

              ## Inform user this might take a while
              cat("Averaging over polygons...\n")


              ## If we have a parallel backend use parLapply
              if(opts_FRK$get("parallel") > 1L) {

                  ## Export variavles we need to the cluster
                  clusterExport(opts_FRK$get("cl"),
                                c("basis","s"),envir=environment())

                  ## Compute averaging over footprints in parallel
                  X <- parLapply(opts_FRK$get("cl"),1:length(s), function(i) {
                      samps <- .samps_in_polygon(basis,s,i)                # sample 1000 times (fixed) uniformly in polygon
                      colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)  # This is the averaging
                  })
                  clusterEvalQ(opts_FRK$get("cl"), {gc()})                 # clear the cluster memory

              } else  {
                  ## Same as above but serially
                  X <- lapply(1:length(s), function(i) {
                      samps <- .samps_in_polygon(basis,s,i)
                      colSums(.point_eval_fn(basis@fn,samps))/nrow(samps)
                  })
              }

              X <- Reduce("rBind",X)   # join the rows together
              as(X,"Matrix")           # coerce to Matrix if not already Matrix

          })

#' @rdname eval_basis
#' @aliases eval_basis,Basis-STIDF-method
setMethod("eval_basis",signature(basis="Basis",s="STIDF"),
          function(basis,s){
              ## Simply evaluate the basis functions at the spatial locations of the data
              ## (i.e., after projecting the data onto space)
              ## Note, tis assumes that the the Basis function is on a spatial manifold
              coords <- coordinates(s)
              if(!dimensions(basis) == ncol(coords))
                  stop("Basis functions need to be spatial")
              .point_eval_fn(basis@fn,coords)
          })

#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-matrix-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s="matrix"),
          function(basis,s){
              n1 <- dimensions(basis@Basis1) # spatial dimensions
              S1 <- eval_basis(basis@Basis1,s[,1:n1,drop=FALSE])    # evaluate at spatial locations
              S2 <- eval_basis(basis@Basis2,s[,-(1:n1),drop=FALSE]) # evaluate at temporal locations

              ## Order: Construct S matrix over space and time
              ## This is ordered as first space then time, therefore we take S2[,1]*S1 as our first block
              ## Then S2[,2]*S1 as our second block etc.
              XX <- lapply(1:ncol(S2),function(i)  (S2[,i] * S1))
              S <- quickcBind(XX)  # a quick cBind method using sparse matrix construction
              S
          })


#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-STIDF-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s = "STIDF"),function(basis,s){
    if(!("t" %in% names(s@data)))
        stop("FRK requires a column with a numeric value for time in the STIDF object.
             This can be, for example, in number of seconds or days from the first data point.")
    slocs <- coordinates(s)               # spatial coordinates
    tlocs <- matrix(s$t)                  # temporal coordinates
    eval_basis(basis,cbind(slocs,tlocs))  # evaluate basis functions over ST locations
})


#' @rdname eval_basis
#' @aliases eval_basis,TensorP_Basis-STFDF-method
setMethod("eval_basis",signature(basis="TensorP_Basis",s = "STFDF"),function(basis,s){
    if(!("t" %in% names(s@data)))
        stop("FRK requires a column with a numeric value for time in the STFDF object.
             This can be, for example, in number of seconds, hours or days from the
             first data point.")
    tlocs <- matrix(s$t)  # all time points
    nt <- length(time(s))  # number of unique time points

    S1 <- eval_basis(basis@Basis1,s[,1])                # evaluate over space (just take first time point)
    S1 <- do.call("rBind",lapply(1:nt,function(x) S1))  # now just repeat that for the nt time points
    S2 <- eval_basis(basis@Basis2,tlocs[,,drop=FALSE])  # evaluate over time (all time points)

    ## As in previous functions we compute S with space running fastest
    XX <- lapply(1:ncol(S2),function(i)  (S2[,i] * S1))
    S <- quickcBind(XX) # a quick cBind method using sparse matrix construction
    S
})

######################################################
########### HELPER FUNCTIONS #########################
######################################################

#' @rdname remove_basis
#' @aliases remove_basis,Basis-method
setMethod("remove_basis",signature(Basis="Basis"),function(Basis,rmidx) {
    ntot <- nbasis(Basis)
    if(!all(rmidx %in% 1:ntot))
        stop("Please make sure indices are numeric and within
             1 and the number of basis functions.")

    Basis_df <- data.frame(Basis)     # extract data frame
    Basis@fn <- Basis@fn[-rmidx]      # remove functions
    Basis@pars <- Basis@pars[-rmidx]  # remove parameters
    Basis@df <- Basis@df[-rmidx,]     # remove rows from data frame
    Basis@n <- nrow(Basis@df)         # reduce n as appropriate
    Basis                             # return basis object
})


#' @rdname local_basis
#' @export
radial_basis <- function(manifold=sphere(),loc=matrix(c(1,0),nrow=1),
                         scale=1,type=c("bisquare","Gaussian","exp","Matern32")) {
    stop("radial_basis is deprecated. Please use local_basis instead")
}

#' @rdname nres
#' @aliases nres_basis,Basis-method
setMethod("nres",signature(b="Basis"), # Returns number of resolutions for standard basis
          function(b){ length(unique(data.frame(b)$res))})

#' @rdname nres
#' @aliases nres_basis,Basis-method
setMethod("nres",signature(b="TensorP_Basis"),  # Returns number of resolutions for TP basis
          function(b){nres(b@Basis1) * nres(b@Basis2)})


#' @rdname nres
#' @aliases nres_SRE,SRE-method
setMethod("nres",signature(b="SRE"),  # Returns number of resolutions of basis in SRE model
          function(b){ nres(b@basis)})


#' @export
#' @rdname Basis_data.frame
setMethod( "$", "Basis",  # Allows easy access of fields in basis data frame
           function(x, name ){data.frame(x)[name][,1]} )

#' @export
#' @rdname Basis_data.frame
setMethod( "$<-", "Basis", # Allows easy assignment of fields in basis data frame
           function(x,name,value){
               df <- data.frame(x)
               df[name] <- value
               data.frame(x) <- df
               x
           })

#' @rdname Basis_data.frame
#' @aliases data.frame_Basis,Basis-method
#' @export
setMethod( "data.frame<-", "Basis", # Allows easy assignment of data.frame to Basis
           function(x, value){x@df <- value; x})

#' @rdname Basis_data.frame
#' @aliases data.frame_Basis,Basis-method
#' @export
setMethod( "data.frame<-", "TensorP_Basis", # Allows easy assignment of data.frame to TensorPbasis
           function(x, value){x@df <- value; x})

#' @title Basis-function data frame object
#' @description Tools for retrieving and manipulating the data frame within the Basis objects. Use the assignment \code{data.frame()<-} with care; no checks are made to make sure the data frame conforms with the object. Only use if you know what you're doing.
#' @param x the obect of class \code{Basis} we are assigning the new data to or retrieving data from
#' @param value the new data being assigned to the Basis object
#' @param name the field name to which values will be retrieved or assigned inside the Basis object's data frame
#' @param ... unused
#' @rdname Basis_data.frame
#' @examples
#' G <- local_basis()
#' df <- data.frame(G)
#' print(df$res)
#' df$res <- 2
#' data.frame(G) <- df
#' @export
as.data.frame.Basis = function(x,...) # Used to convert basis into its summary data frame
    x@df
setAs("Basis", "data.frame", function(from) as.data.frame.Basis(from))

#' @rdname Basis_data.frame
#' @export
as.data.frame.TensorP_Basis = function(x,...) # Used to convert basis into its summary data frame
    x@df
setAs("TensorP_Basis", "data.frame",
      function(from) as.data.frame.TensorP_Basis(from))

#' @rdname nbasis
#' @aliases nbasis,Basis_obj-method
setMethod("nbasis",signature(.Object="Basis_obj"), # Returns number of basis functions for Basis
          function(.Object) {return(.Object@n)})


#' @rdname nbasis
#' @aliases nbasis,SRE-method
setMethod("nbasis",signature(.Object="SRE"), # Returns number of basis functions for SRE model
          function(.Object) {return(nbasis(.Object@basis))})

#' @aliases count_res,Basis-method
setMethod("count_res",signature="Basis", # Returns count by resolution for Basis
          function(.Object) {
              res <- NULL # suppress bindings (it's in the data frame)
              count(data.frame(.Object),res)
          })

#' @aliases count_res,SRE-method
setMethod("count_res",signature="SRE", # Returns count by resolution for SRE model
          function(.Object) {
              count_res(.Object@basis)
          })

#' @aliases count_res,TensorP_Basis-method
setMethod("count_res",signature="TensorP_Basis", # Returns count by resolution for Tensor Basis
          function(.Object) {
              res <- NULL # suppress bindings (it's in the data frame)
              c1 <-  count(data.frame(.Object@Basis1),res) # count spatial by resolution
              c2 <-  count(data.frame(.Object@Basis2),res) # count temporal by resolution

              ## In the below we define resolutions with the spatial resolution moving fastest
              ## So for example if we have resolutions 1,2, and 3 for space and 1, and 2 for time,
              ## then we will have resolutions 1,2,3,4,5,6, where 1,2,3 correspond to space 1,2,3 and
              ## time 1, and 4,5,6 correspond to space 1,2,3 and time 2

              c_all <- NULL                  # initialise NULL
              max_res_c1 <- c1$res[1] - 1    # initialise at one less than coarsest resolution
              for( i in 1:nrow(c2)) {        # for each temporal resolution (should be only one)
                  new_res <- (max_res_c1 + 1):(max_res_c1 + nrow(c1)) # range of resolutions
                  temp_c1 <- c1              # spatial resolution count
                  temp_c1$res <- new_res     # but change resolution numbers
                  temp_c1$n <- temp_c1$n * c2$n[i] # we have n basis functions for each time point
                  c_all <- rbind(c_all,temp_c1) # append
                  max_res_c1 <- max(c_all$res)  # update maximum resolution count
              }
              c_all
          })

## Print/Show basis functions
print.Basis <- function(x,...) {
    cat("Number of basis functions:",nbasis(x),"\n")
    cat("Number of resolutions:",nres(x),"\n")
    cat("Type of manifold:",type(manifold(x)),"\n")
    cat("Dimension of manifold:",dimensions(manifold(x)),"\n")
    cat("First basis function:\n",deparse(x@fn[[1]]),"\n")
}
setMethod("show",signature(object="Basis"),function(object) print(object))

## Print/Show Tensor product basis functions
print.TensorP_Basis <- function(x,...) {
    cat("First set of basis functions\n")
    cat("----------------------------\n")
    print(x@Basis1)
    cat("\n\n")
    cat("Second set of basis functions\n")
    cat("-----------------------------\n")
    print(x@Basis2)
    cat("\n\nTotal number of basis functions:",nbasis(x))
}
setMethod("show",signature(object="TensorP_Basis"),function(object) print(object))


## Summary Basis
summary.Basis <- function(object,...) {
    summ <- summary(as.data.frame(object))
    class(summ) <- "summary.Basis"
    summ
}
setMethod("summary",signature(object="Basis"),summary.Basis)

print.summary.Basis <- function(x,...) {
    cat("Summary of Basis data frame:\n")
    print(as.table(x))
    cat("For object properties use show().\n")
    invisible(x)
}

######################################################
###########  FUNCTIONS NOT EXPORTED ##################
######################################################

## BuilD computes the distance matrices for the basis-function centroids
setMethod("BuildD",signature(G="Basis"),function(G){
    res <- NULL        # suppress bindings (it's in the data frame)
    nres <- nres(G)    # number of resoluations
    m <- manifold(G)   # basis manifold
    n <- dimensions(G) # manifld dimensions

    ## Since the first columns of df are ALWAYS the centroid coordinates
    ## we just take the first n columns when computing the distances
    D_basis = lapply(1:nres,function(i)  {
        x1 <- filter(data.frame(G),res == i)[,1:n] %>%
            as.matrix()
        distance(m,x1,x1)})
    D_basis
})

setMethod("BuildD",signature(G="TensorP_Basis"),function(G){
    nres1 <- nres(G@Basis1) # number of spatial resolutions
    nres2 <- nres(G@Basis2) # number of temporal resolutions
    stopifnot(nres2 == 1)   # only allow for one temporal dimensions
    D_basis <- list(Basis1 = BuildD(G@Basis1), # form distance functions for
                    Basis2 = BuildD(G@Basis2)) # space and time separately
    D_basis
})

## Takes a list of cuntions and evaluates each of these functions over
## the locations s (which here are definitely matrix or numeric)
.point_eval_fn <- function(flist,s) {
    x <- do.call("cbind",
                 sapply(flist,function(f) f(s),simplify=FALSE))
    as(x,"Matrix")
}

## Uniformly samples inside the polygon for carrying out the approximate integration
.samps_in_polygon <- function(basis,s,i) {
    nMC <- 1000                         # 1000 samples
    if(is(manifold(basis),"plane")) {   # If we are on the plane then use spsample
        samps <- coordinates(spsample(s[i,],n=nMC,type="random"))
        ## else sample on the sphere by first drawing a box around the polygon (on the sphere)
        ## sampling in the box, and then keeping only those samples inside the polygon
    } else if(is(basis@manifold,"sphere")){

        ## Find coordinates
        coords <- data.frame(slot(s@polygons[[i]]@Polygons[[1]],"coords"))

        ## Find lat/lon range
        rangelon <- range(coords$lon) + 180
        rangelat <- range(coords$lat) + 90

        ## Find limits in transformed space
        rangeu <- rangelon/360
        rangev <- (cos(rangelat*2*pi/360) + 1)/2

        ## Sample in transformed space
        samps <- cbind(runif(nMC,min=min(rangeu), max=max(rangeu)),
                       runif(nMC,min=min(rangev), max=max(rangev)))

        ## Re-transform back to spherical coordinates
        samps[,1] <- 360*samps[,1] - 180
        samps[,2] <- acos(2*samps[,2] -1) * 360 / (2*pi) - 90

        ## Find which points are in polygon (usually approx. 70%)
        pip <- over(SpatialPoints(samps),  # pip = points in polygon
                    SpatialPolygons(list(s@polygons[[i]]),1L))
        samps <- samps[which(pip == 1),]   # keep those samples for which pip == TRUE

    }
    samps
}

## Basic checks for the basis-function arguments
.check_basis_function_args <- function(manifold,loc,scale) {
    if(!is.matrix(loc))
        stop("Basis functions need to be evaluated using locations inside a matrix")
    if(!(dimensions(manifold) == ncol(loc)))
        stop("Incorrect number of columns for the argument loc")
    if(!is.numeric(scale))
        stop("The scale parmaeter needs to be numeric")
    if(!(scale > 0))
        stop("The scale parameter needs to be greater than zero")
}

# Gaussian Basis Function
.GRBF_wrapper <- function(manifold,mu,std) {
    .check_basis_function_args(manifold,mu,std)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold)) # Internal checking
        dist_sq <- distance(manifold,s,mu)^2
        exp(-0.5* dist_sq/(std^2) )
    }
}

# Bisquare Basis Function
.bisquare_wrapper <- function(manifold,c,R) {
    .check_basis_function_args(manifold,c,R)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold)) # Internal checking
        y <- distance(manifold,s,c)
        (1-(y/R)^2)^2 * (y < R)
    }
}

# Exponential Basis Function
.exp_wrapper <- function(manifold,c,tau) {
    .check_basis_function_args(manifold,c,tau)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold)) # Internal checking
        y <- distance(manifold,s,c)
        exp(-y/tau)
    }
}

# Matern32 Basis Function
.Matern32_wrapper <- function(manifold,c,kappa) {
    .check_basis_function_args(manifold,c,kappa)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold)) # Internal checking
        y <- distance(manifold,s,c)
        (1 + sqrt(3)*y/kappa)*exp(-sqrt(3)*y/kappa)
    }
}

#' @rdname concat
#' @aliases concat,Basis-method
#' @noRd
setMethod("concat",signature = "Basis",function(...) {
    l <- list(...)           # put arguments into list
    if(length(l) < 2)
        stop("Need more than one basis set to concatenate")
    if(!(length(unique(sapply(sapply(l,manifold),type))) == 1))
        stop("Basis need to be on the same manifold")
    G <- l[[1]]              # initialise first basis

    ## We are going to use internals (slots) since this is not exported
    ## This is more direct and safe in this case
    for (i in 2:length(l)) { # add on the other basis functions
        G@fn <- c(G@fn, l[[i]]@fn)        # append functions
        G@pars <- c(G@pars, l[[i]]@pars)  # append parameters
        G@df <- rbind(G@df, l[[i]]@df)    # append data frame
    }
    G@n <- length(G@fn)  # new number of basis functions
    G                    # return new basis
})


