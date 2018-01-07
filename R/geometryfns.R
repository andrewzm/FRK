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

#' @title manifold
#' @param .Object \code{manifold} object passed up from lower-level constructor
#' @description Manifold initialisation. This function should not be called directly as \code{manifold} is a virtual class.
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here if needed
    .Object
})

#' @title real line
#' @description Initialisation of the real-line (1D) manifold.
#' @param measure an object of class \code{measure}
#' @details A real line is initialised using a \code{measure} object. By default, the measure object (\code{measure}) describes the distance between two points as the absolute difference between the two coordinates.
#' @export
#' @examples
#' R <- real_line()
#' print(type(R))
#' print(sp::dimensions(R))
real_line <- function(measure=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(measure)==1L)   # measure needs to take coordinates of length 1
    new("real_line",measure=measure)     # init. object
}

## Constructor for real line
setMethod("initialize",signature="real_line",function(.Object,measure=Euclid_dist(dim=1L)) {
    .Object@type <- "real_line"  # set type
    .Object@measure <- measure   # set measure
    callNextMethod(.Object)})    # pass on to virtual manifold method


#' @title plane
#'
#' @description Initialisation of a 2D plane.
#'
#' @param measure an object of class \code{measure}
#'
#' @details A 2D plane is initialised using a \code{measure} object. By default, the measure object (\code{measure}) is the Euclidean distance in 2 dimensions, \link{Euclid_dist}.
#' @export
#' @examples
#' P <- plane()
#' print(type(P))
#' print(sp::dimensions(P))
plane <- function(measure=Euclid_dist(dim=2L)) {
    stopifnot(dimensions(measure)==2L) # measure needs to take coordinates of length 2
    new("plane",measure=measure)       # init. object
}

## Constructor for plane
setMethod("initialize",signature="plane",function(.Object,measure=Euclid_dist(dim=2L)) {
    .Object@type <- "plane"       # set type
    .Object@measure <- measure    # set measure
    callNextMethod(.Object)})     # pass on to virtual manifold method

#' @title plane in space-time
#' @description Initialisation of a 2D plane with a temporal dimension.
#' @param measure an object of class \code{measure}
#' @details A 2D plane with a time component added is initialised using a \code{measure} object. By default, the measure object (\code{measure}) is the Euclidean distance in 3 dimensions, \link{Euclid_dist}.
#' @export
#' @examples
#' P <- STplane()
#' print(type(P))
#' print(sp::dimensions(P))
STplane <- function(measure=Euclid_dist(dim=3L)) {
    stopifnot(dimensions(measure)==3L)   # measure needs to take coordinates of length 3
    new("STplane",measure=measure)       # init. object
}

## Constructor for STplane
setMethod("initialize",signature="STplane",function(.Object,measure=Euclid_dist(dim=3L)) {
    .Object@type <- "STplane"    # set type
    .Object@measure <- measure   # set measure
    callNextMethod(.Object)})    # pass on to virtual manifold method


#' @title sphere
#' @description Initialisation of the 2-sphere, S2.
#' @param radius radius of sphere
#' @details The 2D surface of a sphere is initialised using a \code{radius} parameter. The default value of the radius \code{R} is \code{R}=6371 km, Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}.
#' @export
#' @examples
#' S <- sphere()
#' print(sp::dimensions(S))
sphere <- function(radius=6371) {
    measure=gc_dist(R=radius)           # defuault measure is GC distance
    stopifnot(dimensions(measure)==2L)  # measure needs to take coordinates of length 2
    stopifnot(radius>0)                 # radius needs to be positive
    new("sphere",measure=measure,       # init. object
        radius=radius)
}

## Constructor for surface of sphere
setMethod("initialize",signature="sphere",function(.Object,radius=1,measure=gc_dist(R=radius)) {
    .Object@type <- "surface of sphere"    # set type
    .Object@measure <- measure             # set measure
    .Object@radius <- radius               # set radius of sphere
    callNextMethod(.Object)})              # pass on to virtual manifold method

#' @title Space-time sphere
#' @description Initialisation of a 2-sphere (S2) with a temporal dimension
#' @param radius radius of sphere
#' @details As with the spatial-only sphere, the sphere surface is initialised using a \code{radius} parameter. The default value of the radius \code{R} is \code{R}=6371, which is the Earth's radius in km, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}. By default Euclidean geometry is used to factor in the time component, so that dist((s1,t1),(s2,t2)) = sqrt(gc_dist(s1,s2)^2 + (t1 - t2)^2). Frequently this distance can be used since separate correlation length scales for space and time  are estimated in the EM algorithm (that effectively scale space and time separately).
#' @export
#' @examples
#' S <- STsphere()
#' print(sp::dimensions(S))
STsphere <- function(radius=6371) {
    measure=gc_dist_time(R=radius)       # the space-time distance function
    stopifnot(dimensions(measure)==3)    # measure needs to take coordinates of length 3
    stopifnot(radius>0)                  # radius needs to be positive
    new("STsphere",measure=measure,      # init. object
        radius=radius)
}

## Constructor for STsphere
setMethod("initialize",signature="STsphere",function(.Object,radius=6371,measure=gc_dist_time(R=radius)) {
    .Object@type <- "STsphere"   # set type
    .Object@measure <- measure   # set measure
    .Object@radius <- radius     # set radius
    callNextMethod(.Object)})    # pass on to virtual manifold method

#' @name distances
#' @aliases measure
#' @aliases Euclid_dist
#' @aliases gc_dist
#' @aliases gc_dist_time
#' @title Pre-configured distances
#' @description Useful objects of class \code{distance} included in package.
#' @param dist a function taking two arguments \code{x1,x2}
#' @param dim the dimension of the manifold (e.g., 2 for a plane)
#' @param R great-circle radius
#' @details Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points.  Currently the Euclidean distance and the great-circle distance are included with \code{FRK}.
#' @export
#' @examples
#' M1 <- measure(distR,2)
#' D <- distance(M1,matrix(rnorm(10),5,2))
measure <- function(dist,dim) {
    ## Basic checks
    if(!is.function(dist))
        stop("dist needs to be a function that accepts dim arguments")
    if(!(is.numeric(dim) | is.integer(dim)))
        stop("dim needs to be an integer, generally 1L, 2L or 3L")
    dim = as.integer(dim)              # coerce to integer if numeric
    new("measure",dist=dist,dim=dim)   # init. object
}

#' @rdname distances
#' @export
Euclid_dist <- function(dim=2L) {
    ## Euclidean distance
    stopifnot(is.integer(dim))   # dimension needs to be integer
    new("measure",               # init. measure object with the distR function
        dist=function(x1,x2) distR(x1,x2), dim=dim)
}

#' @rdname distances
#' @export
gc_dist <- function(R=NULL) {
    ## Great circle distance
    stopifnot(is.null(R) | R > 0)  # radius needs to be positive if specified
    new("measure",                 # init. measure object with the dist_sphere function
        dist=function(x1,x2=NULL)
            dist_sphere(x1,x2,R=R),dim=2L)
}

#' @rdname distances
#' @export
gc_dist_time <- function(R=NULL) {
    ## Great circle distance*time
    stopifnot(is.null(R) | R > 0)  # radius needs to be positive if specified
    new("measure",dist=function(x1,x2)  {
        spatdist <- dist_sphere(x1[,1:2,drop=FALSE],       # spatial distance
                                x2[,1:2,drop=FALSE],R=R)
        tdist <- distR(x1[,3],x2[,3])                      # temporal distance
        sqrt(spatdist^2 + tdist^2) } ,dim=3L)              # combination of the two
}

#' @name dist-matrix
#' @title Distance Matrix Computation from Two Matrices
#' @description This function extends \code{dist} to accept two arguments.
#' @param x1 matrix of size N1 x n
#' @param x2 matrix of size N2 x n
#' @details Computes the distances between the coordinates in \code{x1} and the coordinates in \code{x2}. The matrices \code{x1} and \code{x2} do not need to have the same number of rows, but need to have the same number of columns (dimensions).
#' @return Matrix of size N1 x N2
#' @export
#' @examples
#' A <- matrix(rnorm(50),5,10)
#' D <- distR(A,A[-3,])
distR <- function (x1, x2 = NULL)  {
    ## Try to coerce to matrix
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
    }

    ## If x2 is not specified set it equatl to x1
    if (is.null(x2)) {
        x2 <- x1
    }

    ## If it is specified, coerce it to matrix
    if (!is.matrix(x2)) {
        x2 <- as.matrix(x2)
    }

    ## Basic check
    if(!(ncol(x1) == ncol(x2)))
        stop("x1 and x2 have to have same number of columns")

    ## Compute the distance in C (distR_C is a wrapper)
    distR_C(x1,x2)
}

## Retrieve the border points
setMethod("coordinates",signature(obj="SpatialPolygons"),function(obj){
    coord_vals <- t(sapply(1:length(obj),                              # for each polygon
                           function(i)
                               obj@polygons[[i]]@Polygons[[1]]@labpt)) # retrieve border points

    ## Ensure the column names are the coordinate names
    colnames(coord_vals) <- colnames(obj@polygons[[1]]@Polygons[[1]]@coords)

    ## Return coordinates
    coord_vals
})

# Retrieve dimensions of measure from measure object
#' @aliases dimensions,measure-method
setMethod("dimensions",signature("measure"),function(obj){obj@dim})

# Retrieve dimensions of measure from manifold object
#' @aliases dimensions,manifold-method
setMethod("dimensions",signature("manifold"),function(obj){dimensions(obj@measure)})

# Retrieve dimensions of measure from Basis object
#' @aliases dimensions,Basis-method
setMethod("dimensions",signature("Basis"),function(obj){dimensions(obj@manifold)})

#' @rdname distance
#' @aliases distance,measure-method
setMethod("distance",signature("measure"),function(d,x1,x2=NULL){d@dist(x1,x2)})

# Compute distance on manifold
#' @rdname distance
#' @aliases distance,manifold-method
setMethod("distance",signature("manifold"),function(d,x1,x2=NULL){distance(d@measure,x1,x2)})

# Retrieve type of manifold
#' @rdname type
#' @aliases type,manifold-method
setMethod("type",signature(.Object="manifold"),function(.Object) {
    return(.Object@type)
})

# Retrieve manifold from Basis
#' @rdname manifold
#' @aliases manifold,Basis-method
setMethod("manifold",signature(.Object="Basis"),function(.Object) {
    return(.Object@manifold)
})

# Retrieve manifold from TensorP_Basis
#' @rdname manifold
#' @aliases manifold,TensorP_Basis-method
setMethod("manifold",signature(.Object="TensorP_Basis"),function(.Object) {
    return(list(manifold(.Object@Basis1),
                manifold(.Object@Basis2)))
})

# Retrieve coordnames from sp part of STFDF and add "t" as third coordinate
#' @aliases coordnames,STFDF-method
setMethod("coordnames",signature(x="STFDF"),function(x) {
    return(c(coordnames(x@sp),"t"))
})

# Retrieve coordnames from sp part of STIDF and add "t" as third coordinate
#' @aliases coordnames,STIDF-method
setMethod("coordnames",signature(x="STIDF"),function(x) {
    return(c(coordnames(x@sp),"t"))
})

#' @title Automatic BAU generation
#' @description This function calls the generic function \code{auto_BAU} (not exported) after a series of checks and is the easiest way to generate a set of Basic Areal Units (BAUs) on the manifold being used; see details.
#' @param manifold object of class \code{manifold}
#' @param type either ``grid'' or ``hex'', indicating whether gridded or hexagonal BAUs should be used
#' @param cellsize denotes size of gridcell when \code{type} = ``grid''. Needs to be of length 1 (square-grid case) or a vector of length \code{dimensions(manifold)} (rectangular-grid case)
#' @param isea3h_res resolution number of the isea3h DGGRID cells for when type is ``hex'' and manifold is the surface of a \code{sphere}
#' @param data object of class \code{SpatialPointsDataFrame}, \code{SpatialPolygonsDataFrame},  \code{STIDF}, or \code{STFDF}. Provision of \code{data} implies that the domain is bounded, and is thus necessary when the manifold is a \code{real_line, plane}, or \code{STplane}, but is not necessary when the manifold is the surface of a \code{sphere}
#' @param nonconvex_hull flag indicating whether to use \code{INLA} to generate a non-convex hull. Otherwise a convex hull is used
#' @param convex convex parameter used for smoothing an extended boundary when working on a bounded domain (that is, when the object \code{data} is supplied); see details
#' @param tunit temporal unit when requiring space-time BAUs. Can be either "secs", "mins", "hours" or "days"
#' @param xlims limits of the horizontal axis (overrides automatic selection)
#' @param ylims limits of the vertical axis (overrides automatic selection)
#' @param ... currently unused
#' @details \code{auto_BAUs} constructs a set of Basic Areal Units (BAUs) used both for data pre-processing and for prediction. As such, the BAUs need to be of sufficienly fine resolution so that inferences are not affected due to binning.
#'
#' Two types of BAUs are supported by \code{FRK}: ``hex'' (hexagonal) and ``grid'' (rectangular). In order to have a ``grid'' set of BAUs, the user should specify a cellsize of length one, or of length equal to the dimensions of the manifold, that is, of length 1 for \code{real_line} and of length 2 for the surface of a \code{sphere} and \code{plane}. When a ``hex'' set of BAUs is desired, the first element of \code{cellsize} is used to determine the side length by dividing this value by approximately 2. The argument \code{type} is ignored with \code{real_line} and ``hex'' is not available for this manifold.
#'
#'   If the object \code{data} is provided, then automatic domain selection may be carried out by employing the \code{INLA} function \code{inla.nonconvex.hull}, which finds a (non-convex) hull surrounding the data points (or centroids of the data polygons). This domain is extended and smoothed using the parameter \code{convex}. The parameter \code{convex} should be negative, and a larger absolute value for \code{convex} results in a larger domain with smoother boundaries (note that \code{INLA} was not available on CRAN at the time of writing).
#' @examples
#' ## First a 1D example
#' library(sp)
#' set.seed(1)
#' data <- data.frame(x = runif(10)*10, y = 0, z= runif(10)*10)
#' coordinates(data) <- ~x+y
#' Grid1D_df <- auto_BAUs(manifold = real_line(),
#'                        cellsize = 1,
#'                        data=data)
#' \dontrun{spplot(Grid1D_df)}
#'
#' ## Now a 2D example
#' data(meuse)
#' coordinates(meuse) = ~x+y # change into an sp object
#'  if(require(INLA)) {
#'     ## Grid BAUs
#'     GridPols_df <- auto_BAUs(manifold = plane(),
#'                              cellsize = 200,
#'                              type = "grid",
#'                              data = meuse,
#'                              convex=-0.05)
#'     \dontrun{plot(GridPols_df)}
#'
#'     ## Hex BAUs
#'     HexPols_df <- auto_BAUs(manifold = plane(),
#'                             cellsize = 200,
#'                             type = "hex",
#'                             data = meuse,
#'                             convex=-0.05)
#'     \dontrun{plot(HexPols_df)}
#' }
#' @export
auto_BAUs <- function(manifold, type=NULL,cellsize = NULL,
                      isea3h_res=NULL,data=NULL,nonconvex_hull=TRUE,
                      convex=-0.05,tunit=NULL,xlims=NULL,ylims=NULL,...) {

    ## Basic checks and setting of defaults
    if(!(is(data,"Spatial") | is(data,"ST") | is(data,"Date") | is.null(data)))
        stop("Data needs to be of class 'Spatial', 'ST', 'Date', or NULL")
    if(is(data,"Spatial") | is(data,"ST"))           # if data is not NULL
        if((class(coordnames(data)) == "NULL"))
            stop("data needs to have coordinate names")
    if(!is(manifold,"manifold"))
        stop("manifold needs to be of class 'manifold'")

    on_sphere <- grepl("sphere",type(manifold))

    if(is.null(type)) {
        type <- ifelse(on_sphere,"hex","grid")
    }

    ## If user has specified the ISEA3h resolution
    if(!is.null(isea3h_res)) {

        ## Check it's valid
        if(!is.numeric(isea3h_res) | is.integer(isea3h_res))
            stop("isea3h_res needs to be of type 'numeric' or 'integer'")
        if(!on_sphere)
            stop("The problem is not on the surface of sphere. Please set isea3h_res to NULL")

        ## Check it's not too big or too small
        if(!(isea3h_res >=0 & isea3h_res <= 9)) stop("isea3h_res needs to be between 0 and 9")

        ## Coerce type to hex if user wants to use the ISEA3h
        if(type=="grid") {
            type = "hex"
            message("Only hex BAUs possible when setting isea3h_res. Coercing type to 'hex'")
        }

        ## Ensure the resolution is an integer and assign to resl
        resl <- round(isea3h_res)
    } else {
        ## If user did not specify ISEA3h then just set resl to NULL
        resl <- NULL
    }

    ## If we are on the sphere set defaults for type and resolution if not already specified
    if(on_sphere) {
        if(is.null(type)) type <- "hex"
        if(is.null(isea3h_res)) isea3h_res <- 6
    }

    if(is.null(data) & !on_sphere)
        stop("Need to supply data for planar problems")

    ## If user has not supplied cellsize supply it. Note that cellsize is not relevant when not
    ## on sphere since BAUs are given by ISEA3h in this case
    if(is.null(cellsize)) {
        ## if we are on the plane (hence isea3h_res is not set) or we are on the sphere or user has
        ## specified type ``grid'' (even with sphere), then find the cellsize
        if(!on_sphere | type == "grid")
            cellsize <- .choose_BAU_cellsize_from_data(data)

        if(is(manifold,"real_line"))
            cellsize <- 1
    }

    ## If we are not on the sphere
    if(!on_sphere){

        ## Make cellsize have the same dimensions as the manifold if only one number specified
        if(length(cellsize) == 1) cellsize <- rep(cellsize,dimensions(manifold))

        ## If user has specified an incorrect number of cell edges length throw an error
        if(!length(cellsize) == dimensions(manifold))
            stop("cellsize needs to be of length equal to dimension of manifold")
    }

    ## Check xlims and ylims
    if(!is.null(xlims))
        if(!(length(xlims) == 2) & !is.numeric(xlims)) stop("xlims need to be numeric and  of length 2")
    if(!is.null(ylims))
        if(!(length(ylims) == 2) & !is.numeric(ylims)) stop("ylims need to be numeric and  of length 2")

    ## Check and set tunit if we are in a space-time setting
    if(grepl("ST",class(manifold)) & is.null(data) & is.null(tunit))
        stop("Need to specify tunit if data is not specified in ST case")
    if(grepl("ST",class(manifold)) & is.null(tunit))
        tunit <- .choose_BAU_tunit_from_data(data)

    ## Call the internal function with checked arguments
    auto_BAU(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=data,
             nonconvex_hull=nonconvex_hull,convex=convex,tunit=tunit,xlims=xlims,ylims=ylims)
}

## Automatically generate BAUs on the real line
setMethod("auto_BAU",signature(manifold="real_line"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,xlims=NULL,...) {

              if(is.null(d))
                  stop("Data must be supplied when generating BAUs on a plane")

              crs <- CRS(proj4string(d))     # CRS of data
              coords <- coordinates(d)       # coordinates of data

              if(is.null(xlims))               # if x limits are not specified
                  xrange <- range(coords[,1])    # range of coordinates
              else xrange <- xlims             # else just allocate

              drangex <- diff(xrange)        # range of data

              ## Make a SpatialPoints object. Set y = 0 so all points are on
              ## the x-axis
              xgrid <- data.frame(x=seq(xrange[1] - drangex*0.2,
                                        xrange[2] + drangex*0.2,
                                        by=cellsize[1]),
                                  y = 0)
              xy <- SpatialPoints(xgrid,proj4string = crs)

              ## Suppress warning of unknown y grid cell size when converting to
              ## SpatialPixels
              xy <- suppressWarnings(SpatialPixels(xy,proj4string = crs))

              ## Add UIDs
              row.names(xy) <- .UIDs(xy)

              ## Create SpatialPixelsDataFrame
              xy_df <- SpatialPixelsDataFrame(xy,
                                              data.frame(coordinates(xy),
                                                         row.names = row.names(xy)))
              return(xy_df)
          })

## if we have a Date object for BAU generation call a different function after first converting to POSIXct
setMethod("auto_BAU",signature(manifold="real_line",d="Date"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,
                   convex=-0.05,...) {
              d <- as.POSIXct(d)
              auto_BAU_time(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=d,convex=convex,...)
          })

## if we have a POSIXct object for BAU generation dispatch to a different function
setMethod("auto_BAU",signature(manifold="real_line",d="POSIXct"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,
                   convex=-0.05,...) {
              auto_BAU_time(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=d,convex=convex,...)
          })


## if we have a POSIXlt object for BAU generation dispatch to a different function after first converting to POSIXct
setMethod("auto_BAU",signature(manifold="real_line",d="POSIXlt"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,
                   convex=-0.05,...) {
              d <- as.POSIXct(d)
              auto_BAU_time(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=d,convex=convex,...)
          })

## if we have a xts object for BAU generation dispatch to a different function after first converting to POSIXct
setMethod("auto_BAU",signature(manifold="real_line",d="xts"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,
                   convex=-0.05,...) {
              d <- as.POSIXct(time(d))
              auto_BAU_time(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=d,convex=convex,...)
          })



## Construct the BAUs around some time data
auto_BAU_time <- function (manifold,type="grid",cellsize = 1,resl=resl,d=NULL,convex=-0.05,...) {

    ## Extract other user-supplied arguments
    l <- list(...)

    ## User needs to supply a time unit around which to construct the BAUs
    if(!"tunit" %in% names(l))
        stop("Need to supply argument tunit with value secs, mins, hours, days, months or years")

    ## Now we have the time unit and the time points
    tunit <- l$tunit
    tpoints <- d

    ## From which we can extract a range and the duration
    trange <- range(tpoints) # e.g., 1st January 2017, 4th January 2017
    dranget <- diff(trange)  # e.g., 4 days duration
    #dt <- as.difftime(cellsize, units = tunit) # time block size

    ## The time spacing
    tspacing <- paste(cellsize,tunit)      # e.g., paste(1,"days")
    tgrid <- seq(truncPOSIXt(trange[1],tunit), # create grid based on range and spacing by truncating
                 ceil(trange[2]+1,tunit),      # to this time unit (e.g., "days")
                 by=tspacing)                  # and making the interval equal to tunit

    ## Finally round to the time unit (probably not needed)
    tgrid <- suppressWarnings(roundPOSIXt(tgrid, tunit))

    ## NOTE: A suppresswarnings is needed since there seems to be a bug
    ## in trunc.POSIXt when using units faster than days and more than
    ## one element as input. This warning has no adverse affects as far
    ## as I can see.

    ## Ensure it's POSIXct, which is what FRK uses
    tgrid <- as.POSIXct(tgrid)

    attr(tgrid,"tzone") <- attr(d,"tzone") # Make sure the time zone is the same as in the data
    return(tgrid)                          # Return time BAUs
}

## Automatically generating BAUs on the plane
setMethod("auto_BAU",signature(manifold="plane"),
          function(manifold,type="grid",cellsize = c(1,1),resl=resl,d=NULL,
                   nonconvex_hull=TRUE,convex=-0.05,xlims=NULL,ylims=NULL,...) {

              ## To arrange BAUs in a nonconvex hull we need INLA to find the domain boundary
              if(nonconvex_hull)
                  if(!requireNamespace("INLA"))
                      stop("For creating a non-convex hull INLA needs to be installed. Please install it using
                           install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\"). Alternatively
                           please set nonconvex_hull=FALSE to use a simple convex hull.")

              if(is.null(d))
                  stop("Data must be supplied when generating BAUs on a plane")

              crs <- CRS(proj4string(d))   # CRS of data

              X1 <- X2 <- NULL             # Suppress bindings warning

              if(is(d,"SpatialPoints")){    # If data are spatial points
                  coords <- coordinates(d)  # extract coordinates
              } else if(is(d,"SpatialPolygons")){ # if polygons
                  ## get out all edges and concatenate into one big data frame
                  coords <- do.call("rbind",
                                    lapply(1:length(d),
                                           function(i) coordinates(d@polygons[[i]]@Polygons[[1]])))
              }
              coord_names <- coordnames(d) # extract coordinate names

              if(is.null(xlims))              # if xlims not specified
                  xrange <- range(coords[,1])  # find x-range of coordinates
              else xrange <- xlims            # else just allocate

              if(is.null(ylims))             # if ylims not specified
                  yrange <- range(coords[,2])  # y-range of coordinates
              else yrange = ylims            # else just allocate

              ## Increase convex until domain is contiguous and smooth
              ## (i.e., the distance betweeen successive points is small)
              ## This procedure helps avoid holes in the domain when using INLA
              OK <- 0                 # initialise
              while(!OK) {
                  ## Find the hull (convex or non-convex)
                  bndary_seg <- .find_hull(coords,
                                           nonconvex_hull=nonconvex_hull,
                                           convex=convex)

                  ## Find the distance between the boundary points
                  D <- as.matrix(dist(bndary_seg))

                  ## find the distribution of nearest-neighbour distances
                  distances <- unique(band(D,1,1)@x)[-1] # except the first element which is zero
                  if(nonconvex_hull) {
                      ## somtimes we get islands... the following is a simple check for islands. It's not
                      ## very robust and might need to be improved at a later stage (maybe using a simple persistent
                      ## homology algorithm?)
                      OK <- 0.5*sd(distances) < median(distances)

                      ## Update convex to make boundaries smoother
                      convex <- convex*2
                  } else OK <- 1
              }

              ## Consolidate bndary_seg into a SpatialPolygons object
              bndary_seg <- data.frame(x = bndary_seg[,1],
                                       y = bndary_seg[,2],
                                       id = 1)
              bndary_seg <- df_to_SpatialPolygons(bndary_seg,
                                                  keys="id",
                                                  coords=c("x","y"),
                                                  proj=crs)

              drangex <- diff(xrange)  # range of x
              drangey <- diff(yrange)  # range of y

              ## Create x and y grid with 20% buffer and selected cellsizes
              ## If the user has specified the limits do not do buffer
              bufferx <- ifelse(is.null(xlims),0.2,0) # x buffer
              buffery <- ifelse(is.null(ylims),0.2,0) # y buffer

              xgrid <- seq(xrange[1] - drangex*bufferx,
                           xrange[2] + drangex*bufferx,
                           by=cellsize[1])
              ygrid <- seq(yrange[1] - drangey*buffery,
                           yrange[2] + drangey*buffery,
                           by=cellsize[2])
              ## Make a SpatialPoints grid
              xy <- SpatialPoints(expand.grid(x=xgrid,y=ygrid),
                                  proj4string = crs)

              if(type == "hex") {
                  ## If user wants hexagons
                  HexPts <- spsample(xy,type="hexagonal",       # User spsample to generate hexagons on the grid
                                     cellsize = cellsize[1])    # of a certain size
                  idx <- which(!is.na(over(HexPts,bndary_seg))) # Find which hexagons are outside boundary
                  HexPols <- HexPoints2SpatialPolygons(HexPts[idx,])  # Convert hexagons to SpatialPolygons
                  coordnames(HexPols) <- coord_names            # assign coordinate names to both the points
                  coordnames(HexPts) <- coord_names             # and the polygons
                  row.names(HexPols) <- .UIDs(HexPols)           # Create UIDs

                  ## Now create a SpatialPolygonsDataFrame from the polygons
                  ## With the dataframe extracted from the information about the points
                  HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                                         data.frame(
                                                             coordinates(HexPts[idx,]),
                                                             row.names = row.names(HexPols)))
                  ## Return the hexagons
                  return(HexPols_df)
              } else if (type == "grid") {
                  coordnames(xy) <- coord_names               # assign coordinate names
                  xy <- SpatialPixels(xy,proj4string = crs)   # and convert to SpatialPixels

                  ## keep both the pixels inside boundary
                  idx1 <- which(!is.na(over(xy,bndary_seg)))

                  ## and the pixels on boundary
                  bndary_pts <- SpatialPoints(bndary_seg@polygons[[1]]@Polygons[[1]]@coords,
                                              proj4string = crs)
                  idx2 <- unique(over(bndary_pts,xy))

                  ## Double check no boundary points outside grid, otherwise
                  ## remove from the indices we will keep
                  if(any(is.na(idx2)))
                      idx2 <- idx2[-which(is.na(idx2))]

                  # Make sure to include all pixels that contain the data points
                  idx3 <- over(d,xy) # cannot contain NAs by definition of how xy was constructed

                  ## Now take the union of all the indices but only if xlims and ylims were not specified
                  if(is.null(xlims) & is.null(ylims))
                      xy <- xy[union(union(idx1,idx2),idx3),]

                  ## Add UIDs
                  row.names(xy) <- .UIDs(xy)

                  ## Finally we can form our SpatialPixelsDataFrame of all the pixels
                  ## we have left
                  xy_df <- SpatialPixelsDataFrame(xy,
                                                  data.frame(
                                                      coordinates(xy),
                                                      row.names = row.names(xy)))
                  ## Return the pixels
                  return(xy_df)
              }

          })

## Constructing BAUs on the surface of the sphere
setMethod("auto_BAU",signature(manifold="sphere"),
          function(manifold,type="grid",cellsize = c(1,1),resl=2,d=NULL,xlims=NULL,ylims=NULL,...) {

              ## For this function d (the data) may be NULL in which case the whole sphere is covered with BAUs
              if(is.null(d))                                  # set CRS if data not provided
                  prj <- CRS("+proj=longlat +ellps=sphere")
              else {
                  prj <- CRS(proj4string(d))                # extract CRS
                  coords <- data.frame(coordinates(d))      # extract coordinates

                  ## When modelling on the sphere, the CRS needs to be CRS("+proj=longlat +ellps=sphere")
                  if(!identical(prj,CRS("+proj=longlat +ellps=sphere")))
                      stop("If modelling on the sphere please set the CRS of
                           the data to CRS('+proj=longlat +ellps=sphere)")

                  ## When modelling on the sphere, the coordnames need to be (lon,lat)
                  if(!"lat" %in% names(coords) & "lon" %in% names(coords))
                      stop("The coordinate names when modelling on the sphere need to
                           be lat and lon")
              }

              ## If the user wants hexagonal BAUs
              if(type == "hex") {

                  ## Suppress bindings warnings
                  isea3h <- res <- lon <- centroid <- lat <- in_chull <- NULL

                  ## Load the discrete global grids at the desired resolution. This can be either
                  ## from the data in FRK or the dggrids package (depending on how fine the resolution is)
                  isea3h <- load_dggrids(res=resl)

                  ## Split the ISEA3H across the 180 degree boundary using process_isea3h
                  isea3h_res <- process_isea3h(isea3h,resl)

                  ## Now change the dggrid polygons into SpatialPolygons
                  isea3h_sp_pol <- df_to_SpatialPolygons(  # converts data frame to polygons
                      df=filter(isea3h_res,centroid==0),   # do not send in centroid as part of polygon
                      keys=c("id"),                        # ID of BAU
                      coords=c("lon","lat"),               # coordinate names
                      proj=prj)                            # projection

                  ## Create a data frame informing us on the BAUs
                  isea3h_df_info <- filter(isea3h_res,centroid==1)       # centroid od BAU
                  isea3h_df_info <- isea3h_df_info[c("id","lon","lat")]  # keep the ID, lon and lat
                  row.names(isea3h_df_info) <- row.names(isea3h_sp_pol)  # assign the names

                  ## Now attach the informative data frame to make SpatialPolygonsDataFrame
                  sphere_BAUs <- SpatialPolygonsDataFrame(isea3h_sp_pol,isea3h_df_info)

              }  else if (type == "grid") {
                  ## If the user wants a grid

                  ## If the user has specified limits assign xmin,xmax,ymin and ymin clamped to
                  ## the spherical coordinate limits
                  if(!is.null(xlims) & !is.null(ylims)) {
                      xmin <- max(xlims[1],-180)
                      xmax <- min(xlims[2],180)
                      ymin <- max(ylims[1],-90)
                      ymax <- min(ylims[2],90)

                  } else if(!is.null(d)) {
                      ## Else if there is data try to get limits from the data
                      xrange <- range(coords$lon)  # find x-range of coordinates
                      yrange <- range(coords$lat)  # y-range of coordinates

                      ## And how long/wide it is in a lon/lat sense
                      drangex <- diff(xrange)
                      drangey <- diff(yrange)

                      ## Formulate min/max lon and lats, clamping to the 180 lon boundary
                      ## and 90 lat boundary
                      xmin <- max(xrange[1] - drangex*0.2,-180)
                      xmax <- min(xrange[2] + drangex*0.2,180)
                      ymin <- max(yrange[1] - drangey*0.2,-90)
                      ymax <- min(yrange[2] + drangey*0.2,90)
                  } else {
                      ## If data is not supplied then just fill the whole sphere with BAUs
                      xmin <- -180
                      xmax <- 180
                      ymin <- -90
                      ymax <- 90
                  }

                  ## Create the lon/lat rectangular grid with cell centroids at the
                  ## boundaries
                  longrid <- seq(xmin + cellsize[1]/2,xmax - cellsize[1]/2,by=cellsize[1])
                  latgrid <- seq(ymin + cellsize[2]/2,ymax - cellsize[2]/2,by=cellsize[2])

                  ## Now create the lon-lat grid and convert to a GridTopology
                  lonlat <- expand.grid(lon=longrid,lat=latgrid)
                  lonlat <- points2grid(SpatialPoints(lonlat))
                  lonlat <- as.SpatialPolygons.GridTopology(lonlat, proj4string = prj)
                  row.names(lonlat) <- .UIDs(lonlat)

                  ## Ensure that the coordinate names are (lon,lat)
                  coordnames(lonlat) <- c("lon","lat")

                  ## Create a data frame informing us on the BAUs
                  lonlat_df_info <- data.frame(
                      lon = coordinates(lonlat)[,1],  # lon centroids
                      lat = coordinates(lonlat)[,2],  # lat centroids
                      row.names = row.names(lonlat))

                  ## Finally return the BAUs as SpatialPolygonsDataFrame
                  sphere_BAUs <- SpatialPolygonsDataFrame(lonlat,lonlat_df_info)
              }


              ## Now, if the user has supplied us with data, we should cut out BAUs "around" the data
              ## We do this by simply taking a convex hull around the data points
              if(!is.null(d)) {

                  ## Take convex hull of data. .find_hull is an FRK function which adds a bit of buffer
                  conv_hull <- .find_hull(d,                        # data
                                          nonconvex_hull = FALSE,   # convex hull
                                          convex=-0.01)             # buffer of 1%
                  conv_hull <- data.frame(conv_hull)                # .find_hull returns a matrix without column names
                  names(conv_hull) <- coordnames(d)                 # add columns names
                  conv_hull$id <- 1               # just one polygon; set id = 1
                  row.names(conv_hull) <- NULL    # remove row names
                  conv_hull_coords <- conv_hull   # save the coordinates; conv_hull will be a sp object next

                  ## Now make SpatialPolygons object from hull
                  conv_hull <- df_to_SpatialPolygons(conv_hull,                  # data frame
                                                     keys="id",                  # ID
                                                     coords=c("lon","lat"),      # coordinates
                                                     proj=prj)                   # CRS

                  ## Find which BAUs fall outside the hull
                  ## If we have a wide longitude extent then just filter by latitude
                  if(diff(range(coords$lon)) > 270)
                      sphere_BAUs$in_chull <- ifelse((sphere_BAUs$lat < max(conv_hull_coords[,"lat"])) &
                                                         (sphere_BAUs$lat > min(conv_hull_coords[,"lat"])),
                                                     1,NA)

                  ## Otherwise filter by convex hull
                  else  sphere_BAUs$in_chull <- over(sphere_BAUs,conv_hull)

                  ## Remove those BAUs
                  sphere_BAUs <- subset(sphere_BAUs,!is.na(in_chull))

                  ## Remove chull info
                  sphere_BAUs$in_chull <- NULL
              }

              ## Return the final BAUs
              sphere_BAUs

              })


## Constructing BAUs on the surface of the sphere x time
setMethod("auto_BAU",signature(manifold = c("STmanifold")),
          function(manifold,type="grid",cellsize = c(1,1,1),resl=resl,d=NULL,
                   nonconvex_hull=TRUE,convex=-0.05,xlims=NULL,ylims=NULL,...) {

              ## In this function user can opt to just supply a Date object, in which case
              ## the whole surface of the sphere is covered and the temporal part of the BAUs
              ## is extended so as to enclose the temporal span

              ## Now extract the spatial and temporal components from the dataset
              if(is(d,"ST")) {
                  space_part <- d@sp        # spatial
                  time_part <- .time.ST(d)  # temporal
              } else if (is(d,"Date")) {
                  space_part <- NULL        # spatial
                  time_part <- d            # temporal
              } else {
                  stop("Need to supply either a spatio-temporal dataset or an object of class
                       Date to construct BAUs.")
              }

              ## Currently we only have the plane and the sphere
              if(is(manifold,"STplane")) {
                  spat_manifold <- plane()
              } else if (is(manifold,"STsphere")) {
                  spat_manifold <- sphere()
              } else stop("Cannot recognise manifold")

              ## Set cellsize if not supplied. Time cellsize defaults to 1
              if(is.null(cellsize) & !is.null(space_part)) {
                  cellsize_spat <-  .choose_BAU_cellsize_from_data(space_part)
                  cellsize_temp <- 1
              } else {
                  cellsize_spat <- cellsize[1:2]
                  cellsize_temp <- cellsize[3]
              }


              ## Construct the spatial BAUs
              spatial_BAUs <- auto_BAU(manifold=spat_manifold,cellsize=cellsize_spat,
                                       resl=resl,type=type,d=space_part,nonconvex_hull=nonconvex_hull,
                                       convex=convex,xlims=xlims,ylims=ylims,...)

              ## Construct the temporal BAUs
              temporal_BAUs <- auto_BAU(manifold=real_line(), cellsize=cellsize_temp,
                                        resl=resl,type=type,d=time_part,convex=convex,...)

              ## Number of temporal and spatial BAUs
              nt <- length(temporal_BAUs)
              ns <- nrow(spatial_BAUs)

              ## Construct the info data frame on the ST-BAUs
              df_info <- data.frame(n = 1:(nt *ns),                     # BAU number
                                    t = rep(1:nt,each=ns))              # time index

              ## Construct an STFDF based on the spatial and temporal BAUs
              STBAUs <- STFDF(spatial_BAUs,
                              temporal_BAUs,
                              data = df_info)

              ## Return the ST BAUs
              return(STBAUs)

              })

#' @title Convert data frame to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object.
#' @param df data frame containing polygon information, see details
#' @param keys vector of variable names used to group rows belonging to the same polygon
#' @param coords vector of variable names identifying the coordinate columns
#' @param proj the projection of the \code{SpatialPolygons} object. Needs to be of class \code{CRS}
#' @details Each row in the data frame \code{df} contains both coordinates and labels (or keys) that identify to which polygon the coordinates belong. This function groups the data frame according to \code{keys} and forms a \code{SpatialPolygons} object from the coordinates in each group. It is important that all rings are closed, that is, that the last row of each group is identical to the first row. Since \code{keys} can be of length greater than one, we identify each polygon with a new key by forming an MD5 hash made out of the respective \code{keys} variables that in themselves are unique (and therefore the hashed key is also unique). For lon-lat coordinates use \code{proj = CRS("+proj=longlat +ellps=sphere")}.
#' @export
#' @examples
#' library(sp)
#' df <- data.frame(id = c(rep(1,4),rep(2,4)),
#'                  x = c(0,1,0,0,2,3,2,2),
#'                  y=c(0,0,1,0,0,1,1,0))
#' pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
#' \dontrun{plot(pols)}
df_to_SpatialPolygons <- function(df,keys,coords,proj) {

    ## Basic checks
    if(!is(df,"data.frame")) stop("df needs to be a data frame")
    if(!is(keys,"character")) stop("keys needs to be of class character")
    if(!is(coords,"character")) stop("coords needs to be of class character")
    if(!all(keys %in% names(df))) stop("All keys needs to be labels in data frame")
    if(!all(coords %in% names(df))) stop("All coordinate labels needs to be labels in data frame")
    if(!is(proj,"CRS")) stop("proj needs to be of class CRS")

    ## dfun takes a data frame with coordinates for 1 polygon, and makes one POLYGON object from it
    ## with a UID from the polygon key
    dfun <- function(d) {
        Polygons(list(Polygon(d[coords])),digest::digest(d[keys]))
    }

    ## Now apply dfun to all polygons in data frame
    df_poly <- plyr::dlply(df,keys,dfun)

    ## Frorm a SpatialPolygons object from all the returned Polygons
    Sr <- SpatialPolygons(df_poly,             # Polygons
                          1:length(df_poly),   # plotting order
                          proj4string=proj)    # CRS
}

#' @title SpatialPolygonsDataFrame to df
#' @description Convert \code{SpatialPolygonsDataFrame} object to data frame.
#' @param sp_polys object of class \code{SpatialPolygonsDataFrame}
#' @param vars variables to put into data frame (by default all of them)
#' @details This function is mainly used for plotting \code{SpatialPolygonsDataFrame} objects with \code{ggplot} rather than \code{spplot}. The coordinates of each polygon are extracted and concatenated into one long data frame. The attributes of each polygon are then attached to this data frame as variables that vary by polygon \code{id} (the rownames of the object).
#' @export
#' @examples
#' library(sp)
#' library(ggplot2)
#' opts_FRK$set("parallel",0L)
#' df <- data.frame(id = c(rep(1,4),rep(2,4)),
#'                  x = c(0,1,0,0,2,3,2,2),
#'                  y=c(0,0,1,0,0,1,1,0))
#' pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
#' polsdf <- SpatialPolygonsDataFrame(pols,data.frame(p = c(1,2),row.names=row.names(pols)))
#' df2 <- SpatialPolygonsDataFrame_to_df(polsdf)
#' \dontrun{ggplot(df2,aes(x=x,y=y,group=id)) + geom_polygon()}
SpatialPolygonsDataFrame_to_df <- function(sp_polys,vars = names(sp_polys)) {

    ## The names of the polygons is the same
    polynames <- as.character(row.names(sp_polys))


    ## Form a list of data frames, one for each polygon
    list_polys <- lapply(1:length(sp_polys),   # for each polygon
                         function(i) {
                             coords <- sp_polys@polygons[[i]]@Polygons[[1]]@coords # extract coordinates
                             row.names(coords) <- NULL                             # set row names to NULL
                             coords <- data.frame(coords)                          # convert to data frame
                             poldf <- cbind(coords,id=polynames[i],                # cbind the coordinates with the ID
                                            stringsAsFactors=FALSE)                # ID not factor

                             ## remove the rownames from the data frame
                             rownames(poldf) <- NULL

                             ## return the data frame
                             poldf })

    ## rbind the data frames for each polygon into one big data frame
    df_polys <- bind_rows(list_polys)

    ## merge other information from the sp_polys with the data frame (merge by polygon ID)
    df_polys$id <- as.character(df_polys$id)
    sp_polys$id <- row.names(sp_polys)
    cnames <- coordnames(sp_polys)
    vars_no_coords <- vars[which(!vars %in% cnames)]

    if(length(vars_no_coords) > 0)
        df_polys <- left_join(df_polys,
                              sp_polys@data[c("id",vars_no_coords)],by="id")

    ## Return df_polys
    df_polys
}


## Create very small square polygon BAUs around SpatialPoints
#' @rdname BAUs_from_points
#' @aliases BAUs_from_points,SpatialPoints-method
setMethod("BAUs_from_points",signature(obj = "SpatialPoints"),
          function(obj, offset = 1e-10) {

              sp_obj_pols <- NULL                 # Initialise polygons
              cnames <- coordnames(obj)           # coordinate names
              coords <- coordinates(obj)          # coordinates of SpatialPoints

              if(any(duplicated(coords)))
                  stop("Please remove any duplicated data locations from the object before proceeding.")

              ## Generate the Bottom Left, Bottom Right, Top Right, and Top Left, corners of the BAUs
              BL <- data.frame(X1 = coords[,1] - offset, X2 = coords[,2] - offset, id = 1:length(obj))
              BR <- data.frame(X1 = coords[,1] + offset, X2 = coords[,2] - offset, id = 1:length(obj))
              TR <- data.frame(X1 = coords[,1] + offset, X2 = coords[,2] + offset, id = 1:length(obj))
              TL <- data.frame(X1 = coords[,1] - offset, X2 = coords[,2] + offset, id = 1:length(obj))

              ## Interleave them appropriate so they form polygon paths and set names
              sp_obj_pols <- .interleave(BL,BR,TR,TL)
              names(sp_obj_pols) <- c(cnames,"id")

              ## Now create polygons from the above paths, and keep same projection
              sp_obj_pols <- df_to_SpatialPolygons(sp_obj_pols,coords=cnames,keys="id",
                                                   proj = CRS(proj4string(obj)))

              ## We assign the centroid of the BAU to the data object
              df_data <- as.data.frame(coords)

              ## If data points had other variables, add them aswell
              if(is(obj,"SpatialPointsDataFrame"))
                  df_data <- cbind(df_data,obj@data)

              ## Ensure the row names are the same and construct the SpatialPolygonsDataFrame BAUs
              row.names(df_data) <- row.names(sp_obj_pols)
              sp_obj_pols <- SpatialPolygonsDataFrame(sp_obj_pols,data = df_data)
          })

## Create very small square polygon BAUs around SpatialPoints
#' @rdname BAUs_from_points
#' @aliases BAUs_from_points,ST-method
setMethod("BAUs_from_points",signature(obj = "ST"),
          function(obj, offset = 1e-10) {
              cat("BAUs from points for space-time data not yet
                  implemented. Please contact the package maintainer.\n")
          })

## Print/Show manifold
print.manifold <- function(x,...) {
    cat("Type of manifold:",type(x),"\n")
    if(grepl("sphere",type(x)))
        cat("Radius of sphere:",x@radius,"\n")
    cat("Dimension of manifold:",dimensions(x),"\n")
    cat("Distance function:\n",deparse(x@measure@dist),"\n")
}
setMethod("show",signature(object="manifold"),function(object) print(object))


###########################################
########## Not exported ###################
###########################################

## Map the data to the BAUs. This is done after BAU construction
## data_sp: data (SpatialPoints object)
## sp_pols: BAUs (SpatialPolygonsDataFrame or SpatialPixelsDataFrame)
## average_in_BAU: flag indicating whether we want to average data/standard errors in BAUs
## Returns a SpatialPointsDataFrame with the points aligned at the BAU centroids
#' @aliases map_data_to_BAUs,Spatial-method
setMethod("map_data_to_BAUs",signature(data_sp="SpatialPoints"),
          function(data_sp,sp_pols,average_in_BAU = TRUE) {

              ## Suppress bindings warnings
              . <- BAU_name <- NULL

              ## Add BAU ID to the data frame of the SP object
              sp_pols$BAU_name <- as.character(row.names(sp_pols))

              ## Add coordinates to the @data aswell if not aleady there
              if(!(all(coordnames(sp_pols) %in% names(sp_pols@data))))
                  sp_pols@data <- cbind(sp_pols@data,coordinates(sp_pols))

              ## Time how long this takes
              timer <- system.time({

                  ## Find which fields in the data object are not already declared in the BAUs
                  ## These are the variables we will average over
                  diff_fields <- intersect(setdiff(names(data_sp),names(sp_pols)),names(data_sp))

                  ## Create a data frame just of these fields
                  data_df <- data_sp@data[diff_fields]

                  ## The following over returns a data frame equal in number of rows to data_sp
                  ## with the BAU info at the data location
                  data_over_sp <- .parallel_over(data_sp,sp_pols)

                  ## We now cbind the original data with data_over_sp
                  data_over_sp <- cbind(data_df,data_over_sp)

                  if(any(is.na(data_over_sp$BAU_name))) {  # data points at 180 boundary or outside BAUs -- remove
                      ii <- which(is.na((data_over_sp$BAU_name)))
                      data_sp <- data_sp[-ii,]
                      data_over_sp <- data_over_sp[-ii,]
                      warning("Removing data points that do not fall into any BAUs.
                              If you have simulated data, please ensure no simulated data fall on a
                              BAU boundary as these classify as not belonging to any BAU.")
                  }

                  ## We can have multiple data points falling the same BAU. If we wish to
                  ## average over the BAUs, we now apply the mean function to all  data falling
                  ## in the same BAU and convert to data frame. When the safe mean is asked to
                  ## take averages over quantities that are not numeric, it just returns the first
                  ## element of the vector (so, e.g., the below does not crash when averages over
                  ## BAU names are sought)
                  if(average_in_BAU)
                      Data_in_BAU <- group_by(data_over_sp,BAU_name) %>%  # group by BAU
                      summarise_all(.safe_mean) %>%             # apply safe mean to each column BAU
                      as.data.frame()                                     # convert to data frame
                  else Data_in_BAU <- data_over_sp                        # otherwise don't average
                  })                                                          # end timer


              ## We now create a new SpatialPointsDataFrame but this time the data
              ## is averaged over the BAUs, and we have at most one data point per BAU
              new_sp_pts <- SpatialPointsDataFrame(
                  coords=Data_in_BAU[coordnames(data_sp)],         # coordinates of summarised data
                  data=Data_in_BAU,                                # data frame
                  proj4string = CRS(proj4string(data_sp)))         # CRS of original data

              ## Report time taken to bin data
              cat("Binned data in",timer[3],"seconds\n")

              ## Return new matched data points
              new_sp_pts
              })

## Map the data to the BAUs. This is done after BAU construction
## data_sp: data (SpatialPolygons object)
## sp_pols: BAUs (SpatialPolygonsDataFrame or SpatialPixelsDataFrame)
## average_in_BAU: flag indicating whether we want to average data/standard errors in BAUs
#' @aliases map_data_to_BAUs,Spatial-method
setMethod("map_data_to_BAUs",signature(data_sp="SpatialPolygons"),
          function(data_sp,sp_pols,average_in_BAU = TRUE)
          {
              ## Suppress bindings warnings
              . <- BAU_name <- NULL

              ## SpatialPixels have equal area while SpatialPolygons need not.
              ## Currently we are not weighting by the different BAU area.
              ## Inform user of this
              if(!is(sp_pols,"SpatialPixels"))
                  message("BAUs are Polygons and not Pixels. Currently BAU of identical
                          area are being assumed when computing the incidence matrix
                          from observations having a large support.
                          Handling of different areas will be catered for in a future revision.
                          Please report this issue to the package maintainer.")

              ## Attach the ID of the data polygon to the data frame
              data_sp$id <- row.names(data_sp)

              ## Assume the BAUs are so small that it is sufficient to see whether the
              ## BAU centroid falls in the data polygon. To do this we first make
              ## A SpatialPointsDataFrame from the BAUs reflecting the BAU centroids
              BAU_as_points <- SpatialPointsDataFrame(coordinates(sp_pols),
                                                      sp_pols@data,
                                                      proj4string = CRS(proj4string(sp_pols)))

              ## Now see which centroids fall into the BAUs
              ## The following returns a data frame equal in number of rows to
              ## the data polygons, with all the BAU features averaged (hence if
              ## BAU_as_points$xx = c(1,2,3) for those BAUs inside the data polygon
              ## BAUs_aux_data$xx = 3.
              BAUs_aux_data <- .parallel_over(data_sp,BAU_as_points,fn=.safe_mean)

              ## Now include the ID in the table so we merge by it later
              BAUs_aux_data$id <- row.names(BAUs_aux_data)

              ## Do the merging
              updated_df <- left_join(data_sp@data,BAUs_aux_data,by="id")

              ## Make sure the rownames are OK
              row.names(updated_df) <- data_sp$id

              ## Allocated data frame to SpatialPolygons object
              data_sp@data <- updated_df

              ## Return Spatial object
              data_sp
          })

## Map the data to the BAUs. This is done after BAU construction
## data_sp: data (SpatialPixels object)
## sp_pols: BAUs (SpatialPolygonsDataFrame or SpatialPixelsDataFrame)
## average_in_BAU: flag indicating whether we want to average data/standard errors in BAUs
#' @aliases map_data_to_BAUs,Spatial-method
setMethod("map_data_to_BAUs",signature(data_sp="SpatialPixels"),
          function(data_sp,sp_pols,average_in_BAU = TRUE) {
              if(is(data_sp, "SpatialPixels")) {
                  data_sp <- as(data_sp, "SpatialPolygons")
              } else {
                  data_sp <- as(data_sp, "SpatialPolygonsDataFrame")
              }
              map_data_to_BAUs(data_sp, sp_pols, average_in_BAU = average_in_BAU)
          })

## Returns either a STIDF with the data at the BAU centroids (if data_sp is STIDF)
## Or else an STFDF with the original data shifted to the BAU time points and with BAU
## features averaged over the ST data polygons
setMethod("map_data_to_BAUs",signature(data_sp="ST"),
          function(data_sp,sp_pols,average_in_BAU = TRUE) {

              ## Initialise to no spatial field
              sp_fields <- NULL

              ## Coerce to STIDF if necessary and then project all the space-time data onto space
              data_all_spatial <- as(as(data_sp,"STIDF"),"Spatial")

              ## Now we require all dates to be POSIXct, therefore convert
              if(!all(class(data_all_spatial$time) == "POSIXct")) {
                  data_all_spatial$time <- as.POSIXct(data_all_spatial$time)
              }


              ## Bin every spatial frame separately. The following returns a list of Spatial objects
              ## that are either SpatialPoints or SpatialPolygons, depending on data_sp
              sp_fields <- lapply(seq_along(sp_pols@time),
                                  function(i) {
                                      ## See which initial time we are at
                                      t1 <- time(sp_pols)[i]

                                      ## If this is not the last (initial) time point
                                      if(i < last(sp_pols@time)) {
                                          ## Then mark the beginning of the next time interval as the end of this one
                                          t2 <- time(sp_pols)[i+1]

                                      } else {

                                          ## If we are the last time interval then lump all data into this interval
                                          t2 <- last(data_sp@endTime) + 1
                                      }

                                      ## Now we know which data to bin in space, those appearing between t1 and t2
                                      data_spatial <- subset(data_all_spatial, time >= t1 & time < t2)


                                      ## If there are not data points in this BAU interval then do nothing (return NULL).
                                      ## Otherwise
                                      if(nrow(data_spatial) > 0) {

                                          ## Remove time info from data now
                                          data_spatial$time <- NULL

                                          ## Extract the BAUs of this time interval
                                          BAU_spatial <- sp_pols[,i]

                                          ## For some reason, the above converts the time field to numeric.
                                          ## Replace with actual POSIXct object
                                          BAU_spatial$time <- time(sp_pols)[i]

                                          ## Map the now spatial data the now spatial BAUs
                                          map_data_to_BAUs(data_spatial,
                                                           BAU_spatial,
                                                           average_in_BAU = average_in_BAU)
                                      } else {
                                          NULL
                                      }})


              ## Next we are going to construct our data frame which will be part of the return ST object
              ## This is based on the spatial mapping done in the different time intervals

              ## Initialise
              time <- time_single <- sp <- n <- NULL


              ## For each BAU time point
              for(i in seq_along(sp_pols@time)) {

                  ## If there is some data in this BAU
                  if(!is.null(sp_fields[[i]])) {

                      ## Concatenate into a new data frame sp
                      sp <- rbind(sp,sp_fields[[i]]@data)

                      ## The time field is simply the current time * number of data points
                      n <- nrow(sp_fields[[i]])
                      this_time <- rep(time(sp_pols)[i],n)

                      ## If this is the first iteration then allocate this_time, otherwise
                      ## concatenate (cannot just use c() in both cases)
                      if(is.null(time)) time <- this_time else time <- c(time,this_time)
                      if(is.null(time_single)) time_single <- this_time[1] else time_single <- c(time_single,this_time[1])

                      ##F Finally this ensures the labels are from non-null field which could be first or last
                      coordlabels <- coordnames(sp_fields[[i]])
                  }
              }

              ## Now if original data was an STIDF, we are going to create a new STIDF with data centred at the BAU
              ## centroids (we can do this because BAUs are our smallest unit of consideration)
              if(is(data_sp,"STIDF")) {
                  coordinates(sp) <- coordlabels
                  sp@data <- cbind(sp@data,coordinates(sp))

                  STIDF(as(sp,"SpatialPoints"),
                        time,
                        data = sp@data)
              } else {
                  ## If data_sp is STFDF then we return the same data_sp, but with time shifted to those
                  ## of the BAUs and covariates averaged over the containing BAUs
                  STFDF(data_sp@sp,
                        time_single,
                        data = sp)
              }

          })

## The following three methods build the incidence C matrices. This can be
## C (over all the BAUs) CZ (at data locations, over data polygons) and CP (at pred. polys)
## For these functions to work, the data must have already been mapped to the BAUs using
## map_data_to_BAUs.

## The SpatialPoints object must be a SpatialPointsDataFrame, where the data frame
## contains a field "BAU_name" indicating which BAU the data point falls in. This would have
## been attached to the SpatialPointsDataFrame using map_data_to_BAUs.
setMethod("BuildC",signature(data="SpatialPoints"),
          function(data,BAUs) {
              if(!is(data,"SpatialPointsDataFrame"))
                  stop("In BuildC, the SpatialPoints must have a data frame
                       attached containing BAU_name. This is an internal function,
                       please contact package maintainer for assistance.")
              BAU_index <- data.frame(row.names=row.names(BAUs),  # names of BAUs
                                      n =1:length(BAUs))          # column number of C matrix
              i_idx <- 1:length(data)                             # row number (simply 1:ndata)
              j_idx <- BAU_index[data$BAU_name,]                  # column number reflects the BAU the data falls in
              list(i_idx=i_idx,j_idx=j_idx)                       # return the (i,j) indices of nonzeros
          })

## The BuildC method for when we have Polygon data. Note that in this case we haven't allocated
## the BAUs to the data yet
setMethod("BuildC",signature(data="SpatialPolygons"),
          function(data,BAUs) {
              data$id <- 1:length(data)                           # polygon number
              BAU_as_points <- SpatialPoints(coordinates(BAUs))   # convert BAUs to SpatialPoints
              i_idx <- j_idx <-  NULL                             # initialise
              for (i in 1L:length(data)) {                        # for each data point
                  this_poly <- SpatialPolygons(list(data@polygons[[i]]),1L) # extract polygon
                  overlap <- which(over(BAU_as_points,this_poly) == 1)      # see which BAUs are overlapped
                  i_idx <- c(i_idx,rep(i,length(overlap)))                  # the row index is the data number repeated
                  j_idx <- c(j_idx,as.numeric(overlap))                     # the column index is the BAU number
              }

              ## If no overlap was found it means the user generated the BAUs and that
              ## they don't overlap all the observations. This could also be because
              ## many datasets were provided and FRK did not correctly allocate a domain
              ## which covers all the data. In this case the user should construct
              ## the BAUs manually.
              if(any(is.na(j_idx))) stop("NAs when constructing observation from
                                         large support observations. Are you sure all
                                         observations are covered by BAUs?")

              list(i_idx=i_idx,j_idx=j_idx)  # return the (i,j) indices of nonzeros
          })

## If we have an STIDF the BAUs and the data both need to have a field "n"
## which can be used for mapping. Tese fields are automatically created
## by auto_BAUs and map_data_to_BAUs
setMethod("BuildC",signature(data="STIDF"),
          function(data,BAUs) {
              i_idx <- 1:length(data)       # the row index is simple 1:ndata
              j_idx <- BAUs$n[data@data$n]  # the column index is whatever BAU the data falls in
              list(i_idx=i_idx,j_idx=j_idx) # return the (i,j) indices of nonzeros
          })

setMethod("BuildC",signature(data="STFDF"),
          function(data,BAUs) {

              i_idx <- j_idx <- NULL   # suppress bindings warning
              count <- 0L              # initialise count

              ## Since data is STFDF, the C matrix for one time points can be found and then
              ## replicated. Without loss of generality, the one-time-point C matrix is found
              ## by mapping the first time point data with the first time point BAU
              C_one_time <- BuildC(data[,1],
                                   BAUs[,1])

              ## The first row and column indices are those returned by the spatial BuildC
              i <- C_one_time$i_idx
              j <- C_one_time$j_idx

              ## Now, for each time point we just replicate according to the time point
              ## We might need to skip many columns because no data falls into some BAUs
              ## Note that we have not catered for change of support in time. This is
              ## marked for future work
              for(k in seq_along(.time.ST(BAUs))) {
                  ## Find which time
                  overlap_time <- which(as.POSIXct(.time.ST(data)) == # Find which data time matches
                                            (.time.ST(BAUs)[k]))       # the BAU time. This will always
                  # work because we matched BAUs

                  ## If data is covering more than one time point throw error (currently we do not cater)
                  ## for temporal change of support, and all data is assumed to occupy just one temporal BAU
                  if(!length(overlap_time) == 1L)
                      stop("Something is wrong in binning polygon data into BAUs.
                           Note that currently we don't support temporal change of support.")

                  t_idx <- as.numeric(BAUs@time[k])            # find the appropriate time index
                  j_idx <- c(j_idx, (t_idx-1)*nrow(BAUs) + j)  # find the appropriate column indices and append
                  i_idx <- c(i_idx, count*nrow(data) + i)      # row indices are simply shifted by
                  # the amount of spatial locations in the data
                  count <- count + 1                           # increment count
              }
              list(i_idx=i_idx,j_idx=j_idx)                     # return the (i,j) indices of nonzeros
          })

## Does the over function in parallel
.parallel_over <- function(sp1,sp2,fn=NULL,batch_size = NULL) {

    if(!(opts_FRK$get("parallel") > 1)) {    # Either do serially
        over(sp1,sp2,fn=fn)
    } else {                                 # Or in parallel

        if(is.null(batch_size))              # if batch size not set, set such that
            # we get equal load balance
            batch_size <- ceil(length(sp1) / opts_FRK$get("parallel"))

        n1 <- length(sp1) # length of first object
        n2 <- length(sp2) # length of second object

        ## Break n1 into batches of size 1000
        batching=cut(1:n1,breaks = seq(0,n1+batch_size,by=batch_size),
                     labels=F)

        ## Export the objects to the cluster
        clusterExport(opts_FRK$get("cl"),
                      c("batching","sp1","sp2"),envir=environment())

        ## Do the over operation in parallel
        over_list <- parLapply(opts_FRK$get("cl"),1:max(unique(batching)),
                               function(i) {
                                   idx <- which(batching == i)  # subset the sp1 objects and
                                   over(sp1[idx,],sp2,fn=fn)    # do the over
                               })

        clusterEvalQ(opts_FRK$get("cl"), {gc()})               # clear the cluster memory

        if(is(over_list[[1]],"data.frame")) {                  # if the over have returned data frames
            over_res <- do.call(rbind,over_list)               # then concatenate them using rbind
        } else {
            over_res <- do.call(c,over_list)                   # otherwise concatenate using c (they are numbers)
        }
        over_res                                               # Return answer
    }
}

## Compute the great circle distance
dist_sphere <- function (x1, x2 = NULL, R = NULL)
{
    ## If R is null set to radius of Earth
    if (is.null(R)) R <- 6378.137

    ## If x2 is NULL set to x1
    if(is.null(x2)) x2 <- x1

    ## Convert lon/lat to radians
    x1 <- x1 * pi/180
    x2 <- x2 * pi/180

    # Formula from https://en.wikipedia.org/wiki/Great-circle_distance
    # d = r.acos(n1.n2) where n1 and n2 are the normals to the ellipsoid at the two positions
    n1 <- cbind(cos(x1[, 2]) * cos(x1[, 1]), cos(x1[, 2]) * sin(x1[, 1]), sin(x1[, 2]))
    n2 <- cbind(cos(x2[, 2]) * cos(x2[, 1]), cos(x2[, 2]) * sin(x2[, 1]), sin(x2[, 2]))
    delta <- sigma <- tcrossprod(n1,n2)

    ## Return gcdist
    return(R * acos(ifelse(abs(delta <- sigma) > 1,  # Clamp to one
                           sign(delta <- sigma),
                           delta <- sigma)))

}

## Loads the dggrids from either the FRK or the dggrids package
load_dggrids <- function (res = 3L){

    isea3h <- NA # suppress binding warning

    ## Basic check
    if(!is.numeric(res))
        stop("res needs to be an integer or vector of integers")

    ## We ship dggrids at res 6 or less with FRK. Finer resolutions are available with the dggrids package
    if(all(res <= 6L))  {
        data(isea3h, envir=environment(),package="FRK")  # load ISEA3h from FRK
    } else {
        if(!requireNamespace("dggrids")) {
            stop("Such fine DGGRID resolutions are not
                 shipped with the package FRK. For this
                 resolution please download and install the
                 package dggrids from https://github.com/andrewzm/dggrids")
        } else {
            data(isea3h,envir=environment(),package = "dggrids") # load ISEA3h from dggrids
        }
    }

    ## Return the ISEA3h data
    return(isea3h)
    }

## Extracts important information from the data based on the formula
.extract.from.formula <- function (formula, data)
{
    m = model.frame(terms(formula),            # create data frame based on terms
                    as(data, "data.frame"),    # after coercing data to data frame
                    na.action = na.fail)       # and do not accept NAs
    Y = model.extract(m, "response")           # Y is the response
    if (length(Y) == 0)                        # throw an error if there is no response variable in data
        stop("no response variable present in formula")
    Terms = attr(m, "terms")                   # extract the terms
    X = model.matrix(Terms, m)                 # and form the covariate matrix from these
    list(y = Y, X = X)                         # Return Y (data) and X (covariates)
}


process_isea3h <- function(isea3h,resl) {
    ## Splits the polygons at the 180 boundary
    ## Algorithm adapted from
    ## https://stat.ethz.ch/pipermail/r-sig-geo/2015-July/023168.html

    ## suppress bindings warning
    res <- lon <- probpoly <- centroid <- lat <- NULL

    ## We need rgeos to process these polygons
    if(!requireNamespace("rgeos"))
        stop("rgeos is required for processing hexagons on the sphere.
             Please install using install.packages().")

    isea3h_res <- filter(isea3h,res == resl) %>%
        arrange(id) %>%
        group_by(id) %>%
        data.frame() %>%
        group_by(id) %>%
        mutate(probpoly= diff(range(lon)) > 180)

    last_id <- max(isea3h_res$id)

    prob_polys <- filter(isea3h_res,
                         probpoly == TRUE)
    prob_polys2 <- prob_polys
    prob_polys2$lon <- prob_polys2$lon + 360*(prob_polys2$lon < 0)

    prob_polys2_sp <- df_to_SpatialPolygons(
        df=filter(prob_polys2,centroid==0),
        keys=c("id"),
        coords=c("lon","lat"),
        proj=CRS())

    line = SpatialLines(list(Lines(list(Line(cbind(lon=c(180,180),lat=c(-90,90)))),
                                   ID="line")))
    new_polys <- NULL
    for(i in 1:length(prob_polys2_sp)) {
        # Just ignore if cannot find intersection, might create some small gaps in sphere (?)
        lpi <- tryCatch(rgeos::gIntersection(prob_polys2_sp[i,], line),
                        error=function(e) {TRUE})

        if(!is(lpi,"logical")) {
            blpi <- rgeos::gBuffer(lpi, width = 0.000001)        # create a very thin polygon
            dpi <- rgeos::gDifference(prob_polys2_sp[i,], blpi)  # split using gDifference

            pol1 <- dpi@polygons[[1]]@Polygons[[1]]@coords
            pol2 <- dpi@polygons[[1]]@Polygons[[2]]@coords
            idx1 <- which((abs(pol1[,1] - 180)) < 0.00001)
            idx2 <- which((abs(pol2[,1] - 180)) < 0.00001)

            if(mean(pol1[,1]) > 180) {  # Polygon is to the right
                pol1[idx1,1] <- 180 + 0.00001
            } else {
                pol1[idx1,1] <- 180 - 0.00001
            }

            if(mean(pol2[,1]) > 180) {  # Polygon is to the right
                pol2[idx2,1] <- 180 + 0.00001
            } else {
                pol2[idx2,1] <- 180 - 0.00001
            }

            colnames(pol1) <- colnames(pol2) <- c("lon","lat")
            new_polys <- rbind(new_polys,
                               cbind(pol1,
                                     id = last_id  + (i-1)*2+1,
                                     res = resl,
                                     centroid = 0,
                                     probpoly = FALSE),
                               cbind(pol2,
                                     id = last_id  + i*2,
                                     res = resl,
                                     centroid = 0,
                                     probpoly = FALSE))
        }
    }
    new_polys <- as.data.frame(new_polys)

    centroids <- new_polys %>%
        group_by(id) %>%
        summarise(lon=mean(lon),lat=mean(lat),res=res[1],centroid=1,probpoly=0)
    new_polys <- rbind(new_polys,centroids)

    new_polys$lon <- new_polys$lon - 360*(new_polys$lon >= 180)

    ## Put polygon points back on boundary
    idx <- which(new_polys$lon > 179.999)
    new_polys$lon[idx] <- 180
    idx <- which(new_polys$lon < -179.999)
    new_polys$lon[idx] <- -180

    isea3h_res2 <- isea3h_res %>%
        data.frame() %>%
        filter(probpoly == FALSE) %>%
        rbind(new_polys)

    # ## Consolidate edges
    # isea3h_res2 <- isea3h_res2 %>%

    isea3h_res2

}


## Choose the manifold based nn data
.choose_manifold_from_data <- function(data) {

    ## Basic check
    if(!(is(data,"Spatial") | is(data,"ST")))
        stop("data needs to be Spatial or ST")

    if(is(data, "Spatial")) {             # if data is spatial
        p4 <- proj4string(data)           # extract proj4string
        manifold = plane()                # default to the plane
        if(!is.na(p4))                    # if there is a non-NA CRS
            if(grepl("longlat",p4))       # if longlat is in CRS
                manifold = sphere()       # then we're on the sphere

    } else {                              # if data is ST
        p4 <- proj4string(data@sp)        # extract proj4string
        manifold = STplane()              # default to the STplane
        if(!is.na(p4))                    # if there is a non-NA CRS
            if(grepl("longlat",p4))       # if longlat is in CRS
                manifold = STsphere()     # then we're on the STsphere
    }
    manifold                              # return the manifold
}

## Automatically choose the BAU cellsize from the data
.choose_BAU_cellsize_from_data <- function(data) {

    ## Basic check
    if(!(is(data,"Spatial") | is(data,"ST")))
        stop("data needs to be Spatial or ST")

    coords <- coordinates(data)             # extract coordinates
    xrange <- diff(range(coords[,1]))       # find range of x
    yrange <- diff(range(coords[,2]))       # find range of y

    cellsize <- c(xrange/100,yrange/100)   # plan for a 100 x 100 BAU grid
    if (is(data,"Spatial")) {
        cellsize                           # if there's no time we're done
    } else {
        c(cellsize,1)                       # otherwise add a cellsize of 1 time unit
    }
}

## Automatically choose the time unit from the data
.choose_BAU_tunit_from_data <- function(data) {

    ## Aim for no more than 40 BAUs in time
    trange <- range(.time.ST(data))  # find the range of time we have
    t1 <- trange[1]                  # initial time
    t2 <- trange[2]                  # final time

    ## We will try to choose between days/weeks/month/years
    tunits <- c("days","weeks","months","years")

    ## For each option
    for(i in seq_along(tunits)) {

        ## See how many units (e.g., days) we would need to cover the span
        l <- length(seq(t1,t2,by=tunits[i]))

        ## If we need more than 40 try again with coarser unit, otherwise stop
        if(l < 40) break
    }

    ## Return the chosen tunit
    tunits[i]
}

## Convert polygons to points (centroids)
.polygons_to_points <- function(polys) {

    ## Basic check
    if(!(is(polys,"STFDF") |  is(polys,"SpatialPixels")| is(polys,"SpatialPolygons")))
        stop("polys needs to be of class STFDF, SpatialPixels, or SpatialPolygons")

    ## If object is STFDF
    if(is(polys,"STFDF")) {
        if(!("t" %in% names(polys)))
            as.matrix(cbind(coordinates(polys),polys@data$t)) # return matrix in the form [x,y,t] if t is present
    } else {
        as.matrix(coordinates(polys))                     # return matrix in the form [x,y]
    }

}

## Find a hull (convex or nonconvex) around a set of points
.find_hull <- function(coords,nonconvex_hull=TRUE,convex = -0.05) {

    ## If we want a nonconvex hull we need to call INLA
    if(nonconvex_hull) {
        bndary_seg = INLA::inla.nonconvex.hull(coords,convex=convex)$loc

    } else {
        ## Otherwise we just find a convex hull
        chull_idx <- chull(coordinates(coords))         # find which points are on the hull
        conv_hull <- coordinates(coords)[chull_idx,]    # extract those points
        bound_box <- bbox(coords)                       # find the bounding box of the points
        centroid <- apply(bound_box,1,mean)             # and the centroid of the bounding box

        ## Now we expand the convex hull
        bndary_seg <- conv_hull                         # initialise hull
        delta <- max(apply(bound_box,1,                 # find the maximum extent (in both x and y)
                           function(x) diff(range(x))))

        ## Now take the hull and expand it in x and y by 5% of delta in the right direction
        bndary_seg[,1] <- conv_hull[,1] + sign(conv_hull[,1] - centroid[1])*delta*(-convex)
        bndary_seg[,2] <- conv_hull[,2] + sign(conv_hull[,2] - centroid[2])*delta*(-convex)

        ## We don't need any column names for this (to match what INLA gives)
        colnames(bndary_seg) <- NULL
    }

    ## Return the hull
    bndary_seg
}

## Return the time index of an ST object
.time.ST <- function (x, ...) {
    if(!is(x,"ST"))
        stop("x needs to be of class ST")
    index(x@time)
}

## Takes a spatial object and finds UIDs for it
.UIDs <- function(x) {
    n <- length(x)
    sapply(rnorm(n),function(x) digest::digest(x,algo="md5"))
}

## Computes the mean of a vector x if x is numeric or logical, otherwise just returns the first element
## This is useful as if we are averaging over several columns using summarise_each, then if one column
## is with characters it just returns the first element, while mean() would crash. Use with caution.
.safe_mean <- function(x) {
    if(is(x,"logical") | is(x,"numeric")) {
        mean(x)
    } else { x[1] }
}

## Takes a formula including covariates, and returns a formula with only the intercept term
.formula_no_covars <- function(f) {
    labs <- all.names(f)                                      # extract all sub-strings in formula
    dep_var <- all.vars(f)[1]                                 # extract dep. variable
    idx <- which(labs == dep_var)                             # first char is ~, second is the dep var
    if(idx > 2)  {                                            # if we have a transformation of dep var
        LHS <- paste0(paste0(labs[2:idx],collapse = "("),     # concatenate all transfrmations
                      paste0(rep(")",idx-2,collapse=""),      # and close the brackets
                             collapse=""))
    } else { LHS <- dep_var }                                 # otherwise it's just the dep var
    newf <- formula(paste0(LHS,"~1"))                         # now return formula without covariates
}

## Interleave data frames
.interleave <- function(...) {

    ## Extract data frames
    dfs <- list(...)

    ## Basic checks -- all data frames, same names and all same length
    stopifnot(all(sapply(dfs,function(x) is(x,"data.frame"))))
    stopifnot(all(sapply(dfs,function(x) names(x) == names(dfs[[1]]))))
    stopifnot(all(sapply(dfs,function(x) nrow(x) == nrow(dfs[[1]]))))

    ndfs <- length(dfs)   # number of data frames
    n <- nrow(dfs[[1]])   # number of rows

    stacked <- do.call("rbind",dfs)                       # stack all together
    idx <- rep(1:n, each = ndfs) + (0:(ndfs-1)) * n       # interleaving indices
    interleaved <- stacked[idx,]                          # re-order appropriately
    row.names(interleaved) <- 1:nrow(interleaved)         # reset row names
    interleaved                                           # return

}
