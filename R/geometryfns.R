#' @title manifold
#' @param .Object \code{manifold} object passed up from lower-level constructor
#' @description Manifold initialisation. This function should not be called directly as \code{manifold} is a virtual class.
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here
    .Object
})

#' @title Automatic BAU generation
#' @description This function calls the generic function \code{auto_BAU} (currently not exported) after a series of checks and is the easiest way to generate a set of Basic Areal Units (BAUs) on the manifold being used; see details.
#' @param manifold object of class \code{manifold}
#' @param type either ``grid'' or ``hex'', indicating whether gridded or hexagonal BAUs should be used
#' @param cellsize denotes size of gridcell when \code{type} = ``grid''. Needs to be of length 1 (isotropic-grid case) or a vector of length \code{dimensions(manifold)}
#' @param isea3h_res resolution number of the isea3h DGGRID cells for when type is ``hex'' and manifold is the surface of a \code{sphere}
#' @param data object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}. Provision of \code{data} implies that the domain is bounded, and is thus necessary when the manifold is a \code{real_line} or a \code{plane} but is not necessary when the manifold is the surface of a \code{sphere}
#' @param use_INLA flag indicating whether to use INLA to generate a non-convex hull. Otherwise a convex hull is used
#' @param convex convex parameter used for smoothing an extended boundary when working on a finite domain (that is, when the object \code{d} is supplied), see details.
#' @param tunit temporal unit when requiring space-time BAUs. Can be either "secs", "mins", "hours" or "days".
#' @param ... currently unused
#' @details \code{auto_BAUs} constructs a set of Basic Areal Units (BAUs) used both for data pre-processing and for prediction. As such, the BAUs need to be of sufficienly fine resolution so that data is not adversely affected.
#'
#' Two types of BAUs are supported by \code{FRK}: ``hex'' (hexagonal) and ``grid'' (rectangular). In order to have a ``grid'' set of BAUs, the user should specify a cellsize of length equal to the dimensions of the manifold, that is, of length 1 for \code{real_line} and 2 for the surface of a \code{sphere} and \code{plane}. When a ``hex'' set of BAUs is desired, the first element of \code{cellsize} is used to determine the side length by dividing this value by approximately 2. The argument \code{type} is ignored with \code{real_line} and ``hex'' is not available for this manifold.
#'
#'   If the object \code{data} is provided, then automatic domain selection is carried out by employing the \code{INLA} function \code{inla.nonconvex.hull}, which finds a (non-convex) hull surrounding the data points (or centroids of the data polygons). This domain is extended and smoothed using the \code{convex} parameter. The parameter \code{convex} should be negative, and a larger absolute value for \code{convex} results in a larger domain with smoother boundaries. Due to the dependency on hull construction, \code{INLA} needs to be installed in order to use this function unless BAUs on a sphere are desired (note that \code{INLA} was not available on CRAN at time of writing).
#' @examples
#' ## First a 1D example
#' library(sp)
#' data <- data.frame(x = runif(10)*10, y = 0, z= runif(10)*10)
#' coordinates(data) <- ~x+y
#' Grid1D_df <- auto_BAUs(manifold = real_line(),
#'                        cellsize = 1,
#'                        data=data)
#' spplot(Grid1D_df)
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
#'     plot(GridPols_df)
#'
#'     ## Hex BAUs
#'     HexPols_df <- auto_BAUs(manifold = plane(),
#'                             cellsize = 200,
#'                             type = "hex",
#'                             data = meuse,
#'                             convex=-0.05)
#'     plot(HexPols_df)
#' }
#' @export
auto_BAUs <- function(manifold, type="grid",cellsize = NULL,
                      isea3h_res=NULL,data=NULL,use_INLA=TRUE,
                      convex=-0.05,tunit=NULL,...) {


    if(!(is(data,"Spatial") | is(data,"ST") | is(data,"Date") | is.null(data)))
        stop("Data needs to be of class 'Spatial', 'ST', 'Date', or NULL")
    if(is(data,"Spatial") | is(data,"ST"))
        if((class(coordnames(data)) == "NULL"))
            stop("data needs to have coordinate names")
    if(!is(manifold,"manifold")) stop("manifold needs to be of class 'manifold'")
    if(is.null(type)) {
        if(grepl("longlat",proj4string(data[[1]]))) {
            BAU_type <- "hex"
            if(is.null(isea3h_res)) isea3h_res <- 6
        }
    }
    if(is.null(cellsize)) {
        if(is.null(data)) {
            if(!grepl("sphere",type(manifold))) stop("Need to supply data for planar problems")
        } else {
            if(!grepl("sphere",type(manifold))) cellsize <- .choose_BAU_cellsize_from_data(data)
        }
    }

    if(!is.null(isea3h_res)) {
        if(!is.numeric(isea3h_res) | is.integer(isea3h_res))
            stop("isea3h_res needs to be of type 'numeric' or 'integer'")
        if(!(isea3h_res >=0 & isea3h_res <= 9)) stop("isea3h_res needs to be between 0 and 9")
        if(type=="grid") {
            type = "hex"
            message("Only hex BAUs possible when setting isea3h_res. Coercing type to 'hex'")
        }
        resl <- round(isea3h_res)
    } else { resl <- NULL}

    if(!grepl("sphere",type(manifold))){
        if(length(cellsize) == 1) cellsize <- rep(cellsize,dimensions(manifold))
        if(!length(cellsize) == dimensions(manifold)) stop("cellsize needs to be of length equal to dimension of manifold")
    }


    if(grepl("ST",class(manifold)) & is.null(data) & is.null(tunit)) stop("Need to specify tunit if data is not specified in ST case")
    if(grepl("ST",class(manifold)) & is.null(tunit)) tunit  <-  .choose_BAU_tunit_from_data(data)

    auto_BAU(manifold=manifold,type=type,cellsize=cellsize,resl=resl,
             d=data,use_INLA=use_INLA,convex=convex,tunit=tunit)
}


setMethod("auto_BAU",signature(manifold="plane"),
          function(manifold,type="grid",cellsize = c(1,1),resl=resl,d=NULL,
                   use_INLA=TRUE,convex=-0.05,...) {

              if(use_INLA)
               if(!requireNamespace("INLA"))
                   stop("For creating a non-convex hull INLA needs to be installed. Please install it using
                        install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\"). Alternatively
                        please set use_INLA=FALSE to use a simple convex hull.")

              X1 <- X2 <- NULL # Suppress bindings warning
              if(is(d,"SpatialPoints")){
                  coords <- coordinates(d)
              } else if(is(d,"SpatialPolygons")){
                  ## get out all edges
                  coords <- do.call("rbind",
                                    lapply(1:length(d),
                                           function(i) coordinates(d@polygons[[i]]@Polygons[[1]])))
              }
              coord_names <- coordnames(d)


              xrange <- range(coords[,1])
              yrange <- range(coords[,2])

              ## Increase convex until domain is contiguous and smooth (distance betweeen successive points is small)
              OK <- 0
              while(!OK) {
                  bndary_seg <- .find_hull(coords,use_INLA=use_INLA,convex=convex)
                  D <- dist(bndary_seg) %>% as.matrix()
                  distances <- unique(band(D,1,1)@x)[-1]
                  OK <- 1
                  if(use_INLA) { # somtimes we get islands... check and redo
                      OK <- 0.5*sd(distances) < median(distances)
                      convex <- convex*2
                  }
              }

              bndary_seg <- bndary_seg %>%
                  data.frame() %>%
                  mutate(x=X1,y=X2,id = 1) %>%
                  select(-X1,-X2) %>%
                  df_to_SpatialPolygons(keys="id",coords=c("x","y"),proj=CRS())

              if(!is.null(d)) proj4string(bndary_seg) <- proj4string(d)
              drangex <- diff(xrange)
              drangey <- diff(yrange)
              xgrid <- seq(xrange[1] - drangex*0.2,xrange[2] + drangex*0.2,by=cellsize[1])
              ygrid <- seq(yrange[1] - drangey*0.2,yrange[2] + drangey*0.2,by=cellsize[2])

              xy <- expand.grid(x=xgrid,y=ygrid)  %>%
                  SpatialPoints()
              if(!is.null(d)) proj4string(xy) <- proj4string(d)

              if(type == "hex") {
                  HexPts <- spsample(xy,type="hexagonal",cellsize = cellsize[1])
                  idx <- which(!is.na(over(HexPts,bndary_seg)))
                  HexPols <- HexPoints2SpatialPolygons(HexPts[idx,])
                  coordnames(HexPols) <- coord_names
                  coordnames(HexPts) <- coord_names
                  HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                                         data.frame(
                                                             coordinates(HexPts[idx,]),
                                                             row.names = row.names(HexPols)))
                  return(HexPols_df)
              } else if (type == "grid") {
                  coordnames(xy) <- coord_names
                  # xy <- xy %>%
                  #   points2grid() %>%
                  #    as.SpatialPolygons.GridTopology2()
                  xy <- SpatialPixels(xy)

                  if(!all(coordnames(xy) == coord_names)) {
                      warning("Coordinate names different from (x,y).
                              Renaming polygons might take a while due to
                              structure of the function as.SpatialPolygons.GridTopology.")
                      coordnames(xy) <- coord_names
                  }

                  #keep pixels inside boundary
                  idx1 <- which(!is.na(over(xy,bndary_seg)))
                  # and pixels on boundary
                  bndary_pts <- SpatialPoints(bndary_seg@polygons[[1]]@Polygons[[1]]@coords)
                  if(!is.null(d)) proj4string(bndary_pts) <- proj4string(d)
                  idx2 <- unique(over(bndary_pts,xy))
                  if(any(is.na(idx2))) idx2 <- idx2[-which(is.na(idx2))]
                  # and pixels that contain any points
                  idx3 <- over(d,xy) ## cannot contain NAs
                  xy <- xy[union(union(idx1,idx2),idx3),]

                  xy_df <- SpatialPixelsDataFrame(xy,
                                                    data.frame(
                                                        coordinates(xy),
                                                        row.names = row.names(xy)))
                  return(xy_df)
              }

          })


setMethod("auto_BAU",signature(manifold="timeline"),
          function(manifold,type="grid",cellsize = c(1),resl=resl,d=NULL,
                   convex=-0.05,...) {

              l <- list(...)
              if(is.null(cellsize)) cellsize <- 1

              if(!"tunit" %in% names(l))
                  stop("Need to supply argument tunit with value secs, mins, hours, days, months or years")

              tunit <- l$tunit
              tpoints <- d

              if(is(tpoints,"Date"))
                  tpoints <- as.POSIXct(tpoints)

              trange <- range(tpoints)
              dranget <- diff(trange)

              tspacing <- paste(cellsize,tunit)
              tgrid <- seq(trunc(trange[1],tunit),
                           trunc(trange[2],tunit),
                           by=tspacing) %>%
                               ## Old code:
                               # switch(tunit,
                               #       secs    = cellsize,
                               #       mins    = cellsize*60,
                               #       hours   = cellsize*3600,
                               #       days    = cellsize*3600*24,
                               #       months  = "mon",
                               #       years   = "year")) %>%
                  round(tunit)

              attr(tgrid,"tzone") <- attr(d,"tzone")
              return(tgrid)
          })



setMethod("auto_BAU",signature(manifold = c("STmanifold")),
          function(manifold,type="grid",cellsize = c(1,1,1),resl=resl,d=NULL,
                   use_INLA=TRUE,convex=-0.05,...) {

              if(is(d,"ST")) {
                  space_part <- d@sp
                  time_part <- .time.ST(d)
              } else if (is(d,"Date")) {
                  space_part <- NULL
                  time_part <- d
              } else {
                  stop("Need to supply either a spatio-temporal dataset or
                       a timeline to construct BAUs.")

              }

              if(is(manifold,"STplane")) {
                  spat_manifold <- plane()
              } else {
                  spat_manifold <- sphere()
              }

              spatial_BAUs <- auto_BAU(manifold=spat_manifold,cellsize=cellsize[1:2],
                                       resl=resl,type=type,d=space_part,use_INLA=use_INLA,
                                       convex=convex,...)
              temporal_BAUs <- auto_BAU(manifold=timeline(), cellsize=cellsize[3],
                                        resl=resl,type=type,d=time_part,convex=convex,...)

              nt <- length(temporal_BAUs)
              ns <- nrow(spatial_BAUs)

              STBAUs <- STFDF(spatial_BAUs,
                              temporal_BAUs,
                              data = data.frame(n = 1:(nt *ns),
                                                time = rep(temporal_BAUs,each=ns),
                                                t = rep(1:nt,each=ns)))
              return(STBAUs)

          })


setMethod("auto_BAU",signature(manifold="real_line"),
          function(manifold,type="grid",cellsize = 1,resl=resl,d=NULL,...) {

    coords <- coordinates(d)
    xrange <- range(coords[,1])

    drangex <- diff(xrange)
    xgrid <- data.frame(x=seq(xrange[1] - drangex*0.2,xrange[2] + drangex*0.2,by=cellsize[1])) %>%
        cbind(y=0) %>%
        SpatialPoints()

    suppressWarnings(xy <- xgrid %>%
                         points2grid() %>%
                         as.SpatialPolygons.GridTopology())
    ## Suppress warning of unknown y grid cell size
    xy_df <- SpatialPolygonsDataFrame(xy,data.frame(coordinates(xy),
                                                    row.names = row.names(xy)))
    return(xy_df)
})


setMethod("auto_BAU",signature(manifold="sphere"),
          function(manifold,type="grid",cellsize = c(1,1),resl=2,d=NULL,...) {
    if(type == "hex") {
        isea3h <- res <- lon <- centroid <- lat <- in_chull <- NULL # Suppress bindings warnings


        isea3h <- load_dggrids(res=resl)

        isea3h_res <- process_isea3h(isea3h,resl)

        if(is.null(d)) prj <- CRS("+proj=longlat +ellps=sphere") else prj <-CRS(proj4string(d))
        isea3h_sp_pol <- df_to_SpatialPolygons(
            df=filter(isea3h_res,centroid==0),
            keys=c("id"),
            coords=c("lon","lat"),
            proj=prj)

        isea3h_sp_poldf <- SpatialPolygonsDataFrame(
            isea3h_sp_pol,
            cbind(data.frame(row.names=names(isea3h_sp_pol)),
                  (filter(isea3h_res,centroid==1) %>%
                       select(id,lon,lat))))
        sphere_BAUs <- isea3h_sp_poldf
    }  else if (type == "grid") {

        if(!is.null(d)) {
            coords <- coordinates(d) %>% data.frame()
            stopifnot("lat" %in% names(coords) &
                        "lon" %in% names(coords))
            xrange <- range(coords$lon)
            yrange <- range(coords$lat)
            drangex <- diff(xrange)
            drangey <- diff(yrange)
            xmin <- max(xrange[1] - drangex*0.2,-180)
            xmax <- min(xrange[2] + drangex*0.2,180)
            ymin <- max(yrange[1] - drangey*0.2,-90)
            ymax <- min(yrange[2] + drangey*0.2,90)
        } else {
            xmin <- -180
            xmax <- 180
            ymin <- -90
            ymax <- 90
        }


        longrid <- seq(xmin + cellsize[1]/2,xmax - cellsize[1]/2,by=cellsize[1])
        latgrid <- seq(ymin + cellsize[2]/2,ymax - cellsize[2]/2,by=cellsize[2])
        if(is.null(d)) prj <- CRS("+proj=longlat +ellps=sphere") else prj <-CRS(proj4string(d))
        lonlat <- expand.grid(lon=longrid,lat=latgrid) %>%
            SpatialPoints() %>%
            points2grid() %>%
            as.SpatialPolygons.GridTopology(proj4string = prj)

        coordnames(lonlat) <- c("lon","lat")

        sphere_BAUs <- lonlat %>%
            SpatialPolygonsDataFrame(data.frame(
                lon = coordinates(lonlat)[,1],
                lat = coordinates(lonlat)[,2],
                id = row.names(lonlat),
                row.names = row.names(lonlat)))
    }
    if(!is.null(d)) {
        sub_pols <- coordinates(d) %>%
            chull()
        conv_hull <- coordinates(d)[c(sub_pols,sub_pols[1]),] %>%
            data.frame(row.names=NULL) %>%
            mutate(id=1) %>%
            df_to_SpatialPolygons(keys="id",
                                  coords=c("lon","lat"),
                                  proj=CRS(proj4string(d)))
        sphere_BAUs$in_chull <- over(sphere_BAUs,conv_hull)
        sphere_BAUs <- subset(sphere_BAUs,!is.na(in_chull))
    }
    sphere_BAUs

})


#' @title sphere
#'
#' @description Initialisation of the 2-sphere, S2.
#'
#' @param radius radius of sphere
#'
#' @details The 2D surface of a sphere is initialised using a \code{radius} parameter. The default value of the radius \code{R} is \code{R}=6371 km, Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}.
#' @export
#' @examples
#' S <- sphere()
#' print(sp::dimensions(S))
sphere <- function(radius=6371) {
    measure=gc_dist(R=radius)
    stopifnot(dimensions(measure)==2L)
    stopifnot(radius>0)
    new("sphere",measure=measure,radius=radius)
}

setMethod("initialize",signature="sphere",function(.Object,radius=1,measure=gc_dist(R=radius)) {
    .Object@type <- "sphere"
    .Object@measure <- measure
    .Object@radius <- radius
    callNextMethod(.Object)})


#' @title Space-time sphere
#'
#' @description Initialisation of a 2-sphere (S2) with a temporal dimension
#'
#' @param radius radius of sphere
#'
#' @details As with the spatial-only sphere, the sphere surface is initialised using a \code{radius} parameter. The default value of the radius \code{R} is \code{R}=6371 km, the Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}.
#' @export
#' @examples
#' S <- STsphere()
#' print(sp::dimensions(S))
STsphere <- function(radius=6371) {
    measure=gc_dist_time(R=radius)
    stopifnot(dimensions(measure)==3)
    stopifnot(radius>0)
    new("STsphere",measure=measure,radius=radius)
}

setMethod("initialize",signature="STsphere",function(.Object,radius=6371,measure=gc_dist_time(R=radius)) {
    .Object@type <- "STsphere"
    .Object@measure <- measure
    .Object@radius <- radius
    callNextMethod(.Object)})


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
    stopifnot(dimensions(measure)==2L)
    new("plane",measure=measure)
}


setMethod("initialize",signature="plane",function(.Object,measure=Euclid_dist(dim=2L)) {
    .Object@type <- "plane"
    .Object@measure <- measure
    callNextMethod(.Object)})

#' @title plane in space-time
#'
#' @description Initialisation of a 2D plane with a temporal dimension.
#'
#' @param measure an object of class \code{measure}
#'
#' @details A 2D plane with a time component added is initialised using a \code{measure} object. By default, the measure object (\code{measure}) is the Euclidean distance in 3 dimensions, \link{Euclid_dist}.
#' @export
#' @examples
#' P <- STplane()
#' print(type(P))
#' print(sp::dimensions(P))
STplane <- function(measure=Euclid_dist(dim=3L)) {
    stopifnot(dimensions(measure)==3L)
    new("STplane",measure=measure)
}

setMethod("initialize",signature="STplane",function(.Object,measure=Euclid_dist(dim=3L)) {
    .Object@type <- "STplane"
    .Object@measure <- measure
    callNextMethod(.Object)})


#' @title timeline
#'
#' @description Initialisation of a timeline (real line).
#'
#' @param measure an object of class \code{measure}
#'
#' @details A time axis initialised using a \code{measure} object. By default, the measure object (\code{measure}) is the absolute difference.
#' @export
#' @examples
#' P <- timeline()
#' print(type(P))
#' print(sp::dimensions(P))
timeline <- function(measure=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(measure)==1L)
    new("timeline",measure=measure)
}

setMethod("initialize",signature="timeline",function(.Object,measure=Euclid_dist(dim=1L)) {
    .Object@type <- "timeline"
    .Object@measure <- measure
    callNextMethod(.Object)})


#' @title real line
#'
#' @description Initialisation of the real-line (1D) manifold.
#'
#' @param measure an object of class \code{measure}
#'
#' @details A real line is initialised using a \code{measure} object. By default, the measure object (\code{measure}) describes the distance between two points as the absolute difference between the two coordinates.
#' @export
#' @examples
#' R <- real_line()
#' print(type(R))
#' print(sp::dimensions(R))
real_line <- function(measure=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(measure)==1L)
    new("real_line",measure=measure)
}

setMethod("initialize",signature="real_line",function(.Object,measure=Euclid_dist(dim=1L)) {
    .Object@type <- "real_line"
    .Object@measure <- measure
    callNextMethod(.Object)})



#' @name distances
#' @aliases measure
#' @aliases Euclid_dist
#' @aliases gc_dist
#' @aliases gc_dist_time
#'
#' @title Pre-configured distances
#'
#' @description Useful objects of class \code{distance} included in package.
#'
#' @param dist a function taking two arguments \code{x1,x2}
#' @param dim the dimension of the manifold (e.g., 2 for a plane)
#' @param R great-circle radius
#' @details Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points.  Currently the Euclidean distance and the great-circle distance are included.
#' @export
#' @examples
#' M1 <- measure(distR,2)
#' D <- distance(M1,matrix(rnorm(10),5,2))
measure <- function(dist,dim) {
    if(!is.function(dist)) stop("dist needs to be a function that accepts dim arguments")
    if(!(is.numeric(dim) | is.integer(dim))) stop("dim needs to be an integer, generally 1L, 2L or 3L")
    dim = as.integer(dim)
    new("measure",dist=dist,dim=dim)

}

#' @rdname distances
#' @export
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  distR(x1,x2), dim=dim)
}

#' @rdname distances
#' @export
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2=NULL)  dist_sphere(x1,x2,R=R),dim=2L)
}

#' @rdname distances
#' @export
gc_dist_time <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  {
        spatdist <- dist_sphere(x1[,1:2,drop=FALSE],x2[,1:2,drop=FALSE],R=R)
        tdist <- distR(x1[,3],x2[,3])
        sqrt(spatdist^2 + tdist^2) } ,dim=3L)
}

#' @name dist-matrix
#' @title Distance Matrix Computation from Two Matrices
#'
#' @description This function extends \code{dist} to accept two arguments.
#'
#' @param x1 matrix of size N1 x n
#' @param x2 matrix of size N2 x n
#' @details Computes the distances between the coordinates in \code{x1} and the coordinates in \code{x2}. The matrices \code{x1} and \code{x2} do not need to have the same number of rows, but need to have the same number of columns (dimensions).
#' @return Matrix of size N1 x N2
#' @export
#' @examples
#' A <- matrix(rnorm(50),5,10)
#' D <- distR(A,A[-3,])
distR <- function (x1, x2 = NULL)  {
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
    }
    if (is.null(x2)) {
        x2 <- x1
    }
    if (!is.matrix(x2)) {
        x2 <- as.matrix(x2)
    }
    if(!(ncol(x1) == ncol(x2))) stop("x1 and x2 have to have same number of columns")
    distR_C(x1,x2)
}


#' @title Convert data frame to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object.
#' @param df data frame containing polygon information, see details
#' @param keys vector of variable names used to group rows belonging to the same polygon
#' @param  coords vector of variable names identifying the coordinate columns
#' @param proj the projection of the \code{SpatialPolygons} object. Needs to be of class \code{CRS}
#' @details Each row in the data frame \code{df} contains both coordinates and labels (or keys) that identify to which polygon the coordinates belong. This function groups the data frame according to \code{keys} and forms a \code{SpatialPolygons} object from the coordinates in each group. It is important that all rings are closed, that is, that the last row of each group is identical to the first row. Since \code{keys} can be of length greater than one, we identify each polygon with a new key by forming an MD5 hash made out of the respective \code{keys} variables that in themselves are unique (and therefore the hashed key is also unique). For lon-lat coordinates use \code{proj = CRS("+proj=longlat")}.
#' @export
#' @examples
#' library(sp)
#' df <- data.frame(id = c(rep(1,4),rep(2,4)),
#'                  x = c(0,1,0,0,2,3,2,2),
#'                  y=c(0,0,1,0,0,1,1,0))
#' pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
#' plot(pols)
df_to_SpatialPolygons <- function(df,keys,coords,proj) {
    if(!is(df,"data.frame")) stop("df needs to be a data frame")
    if(!is(keys,"character")) stop("keys needs to be of class character")
    if(!is(coords,"character")) stop("coords needs to be of class character")
    if(!all(keys %in% names(df))) stop("All keys needs to be labels in data frame")
    if(!all(coords %in% names(df))) stop("All coordinate labels needs to be labels in data frame")
    if(!is(proj,"CRS")) stop("proj needs to be of class CRS")

    dfun <- function(d) {
        Polygons(list(Polygon(d[coords])),digest::digest(d[keys]))
    }

    if(opts_FRK$get("parallel") > 1e10) { ## Do not enable, mostly overhead
        ## Deprecated to remove plyr:
        #doParallel::registerDoParallel(opts_FRK$get("parallel"))
        #df_poly <- plyr::dlply(df,keys,dfun,.parallel=TRUE)

        unique_keys <- unique(data.frame(df[keys]))[,1]

        clusterExport(opts_FRK$get("cl"),
                      c("df"),envir=environment())
        df_poly <- parLapply(opts_FRK$get("cl"),unique_keys,
                            function(key) {
                                df[df[keys]==key,] %>%
                                    data.frame() %>%
                                    dfun})
        clusterEvalQ(opts_FRK$get("cl"), {gc()})

        # df_poly <- mclapply(unique_keys,
        #                     function(key) {
        #                         df[df[keys]==key,] %>%
        #                             data.frame() %>%
        #                             dfun},
        #                     mc.cores = opts_FRK$get("parallel"))

    } else {
        df_poly <- plyr::dlply(df,keys,dfun)
    }

    ## Rhipe version (currently disabled)

    # df_poly <- rhwrapper(Ntot = nrow(df),
    #                      N = 4000,
    #                      f_expr = .rhdlply,
    #                      df=df,
    #                      keys=keys,
    #                      coords=coords,
    #                      dfun=parse(text = deparse(dfun)))

    Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=proj)
}

#' @title SpatialPolygonsDataFrame to df
#' @description Convert \code{SpatialPolygonsDataFrame} object to data frame.
#' @param sp_polys object of class \code{SpatialPolygonsDataFrame}
#' @param vars variables to put into data frame (by default all of them)
#' @details This function is mainly used for plotting \code{SpatialPolygonsDataFrame} objects with \code{ggplot} rather than \code{spplot}. The coordinates of each polygon are extracted and concatenated into one long data frame. The attributes of each polygon are then attached to this data frame as variables which vary by polygon \code{id}.  The returned \code{id} variable describes the polygon `id' and ranges from 1 to the number of polygons represented in the data frame.
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
#' ggplot(df2,aes(x=x,y=y,group=id)) + geom_polygon()
SpatialPolygonsDataFrame_to_df <- function(sp_polys,vars = names(sp_polys)) {
    #if(!("id" %in% names(sp_polys@data))) stop("sp_polys has to have an id columns in its data frame")
    #if("id" %in% vars) stop("vars should not contain the variable 'id' (this is implicitly assumed)")
    sp_polys$id <- 1 : length(sp_polys)
    polynames <- 1 : length(sp_polys)
    X <- data.frame(do.call("rbind",lapply(1:length(sp_polys),
                                           function(i) {
                                               poldf <- cbind(sp_polys@polygons[[i]]@Polygons[[1]]@coords,
                                                              id=polynames[i])
                                               rownames(poldf) <- NULL
                                               poldf }))) %>%
        left_join(sp_polys@data[c("id",vars)])
    X
}

.parallel_over <- function(sp1,sp2,fn=fn,batch_size = 1000) {
    n1 <- length(sp1)
    n2 <- length(sp2)
    batching=cut(1:n1,breaks = seq(0,n1+batch_size,by=batch_size),labels=F)

    clusterExport(opts_FRK$get("cl"),
                  c("batching","sp1","sp2"),envir=environment())
    over_list <- mclapply(1:max(unique(batching)),
                          function(i) {
                              idx <- which(batching == i)
                              over(sp1[idx,],sp2,fn=sum)
                          })
    clusterEvalQ(opts_FRK$get("cl"), {gc()})

    if(is(over_list[[1]],"data.frame")) {
        over_res <- do.call(rbind,over_list)
    } else {
        over_res <- do.call(c,over_list)
    }
    over_res
}

#' @aliases map_data_to_BAUs,Spatial-method
setMethod("map_data_to_BAUs",signature(data_sp="Spatial"),
          function(data_sp,sp_pols,av_var,average_in_BAU = TRUE)
          {
              ## Suprress bindings warnings
              BAU_name <- NULL
              . <- NULL
              if(is(data_sp,"SpatialPointsDataFrame")) {

                  if(average_in_BAU) {
                      Nobs <- NULL
                      data_sp$Nobs <- 1

                      ## Deprecated:: The below did the over the other way round which was very
                      ## inefficient if the BAUs could be represented as Pixels
                      # if(opts_FRK$get("parallel") > 1) {
                      #     browser()
                      #     timer <- system.time(Data_in_BAU <-
                      #          .parallel_over(sp_pols,data_sp[c(av_var,"Nobs","std")],fn=sum))
                      # } else {
                      #     timer <- system.time(Data_in_BAU <-
                      #             over(sp_pols,data_sp[c(av_var,"Nobs","std")],fn=sum))
                      # }
                      #
                      #
                      # ## Rhipe VERSION (Currently disabled)
                      # # print("Using RHIPE to find overlays")
                      # # timer <- system.time(
                      # # Data_in_BAU <- rhwrapper(Ntot = length(sp_pols),
                      # #                                  N = 4000,
                      # #                                  f_expr = .rhover,
                      # #                                  sp_pols = sp_pols,
                      # #                                  data_sp = data_sp,
                      # #                                  av_var=av_var)
                      # #     )
                      #
                      #
                      # sp_pols@data[av_var] <- Data_in_BAU[av_var]/Data_in_BAU$Nobs
                      # sp_pols@data["std"] <- Data_in_BAU["std"]/Data_in_BAU$Nobs
                      # sp_pols@data["Nobs"] <- Data_in_BAU$Nobs
                      # sp_pols@data["BAU_name"] <- as.character(row.names(sp_pols))
                      #
                      # new_sp_pts <- SpatialPointsDataFrame(
                      #     coords=sp_pols[coordnames(data_sp)]@data,
                      #     data=sp_pols@data,
                      #     proj4string = CRS(proj4string(data_sp)))
                      # ## If uncommented assumes uncorrelated observations
                      # #new_sp_pts$std <- sqrt(new_sp_pts$std^2 / new_sp_pts$Nobs)
                      # new_sp_pts <- subset(new_sp_pts,!is.na(Nobs))

                      sp_pols@data["BAU_name"] <- as.character(row.names(sp_pols))
                      safe_mean <- function(x) {
                          if(is(x,"logical") | is(x,"numeric")) {
                              mean(x)
                          } else { x[1] }
                      }

                      ## Add coordinates to @data
                      if(!(all(coordnames(sp_pols) %in% names(sp_pols@data)))) {
                          sp_pols@data <- cbind(sp_pols@data,coordinates(sp_pols))
                      }



                      timer <- system.time({
                              data_df <- data_sp@data[setdiff(names(data_sp),
                                                                  names(sp_pols)) %>%
                                                        intersect(names(data_sp))]

                              Data_in_BAU <- cbind(data_df,
                                                   over(data_sp[av_var],
                                                        sp_pols)) %>%
                                                group_by(BAU_name) %>%
                                                summarise_each(funs(safe_mean(.))) %>%
                                             as.data.frame()})

                      new_sp_pts <- SpatialPointsDataFrame(
                          coords=Data_in_BAU[coordnames(data_sp)],
                          data=Data_in_BAU,
                          proj4string = CRS(proj4string(data_sp)))
                  } else {
                      if(opts_FRK$get("parallel") > 1) {
                          timer <- system.time(Data_in_BAU <- .parallel_over(data_sp,
                                                                             as(sp_pols,"SpatialPolygons")))
                      } else {
                          timer <- system.time(Data_in_BAU <- over(data_sp,
                                                                   as(sp_pols,"SpatialPolygons")))
                      }

                      if(any(is.na(Data_in_BAU))) {  # data points at 180 boundary or outside BAUs -- remove

                          ii <- which(is.na((Data_in_BAU)))
                          data_sp <- data_sp[-ii,]
                          Data_in_BAU <- Data_in_BAU[-ii]
                          warning("Removing data points that do not fall into any BAUs.
                                  If you have simulated data, please ensure no simulated data fall on a
                                  BAU boundary as these classify as not belonging to any BAU.")

                      }
                      add_columns <- colnames(sp_pols@data[which(!colnames(sp_pols@data) %in%
                                                                     colnames(data_sp@data))])
                      new_sp_pts <- data_sp
                      new_sp_pts@data <- cbind(new_sp_pts@data,sp_pols@data[Data_in_BAU,][add_columns])
                      new_sp_pts@data["BAU_name"] <- as.character(row.names(sp_pols)[Data_in_BAU])
                      new_sp_pts@data["Nobs"] <- 1

                  }
                  print(paste0("Binned data in ",timer[3]," seconds"))
              } else {

                  if(!is(sp_pols,"SpatialPixels"))
                  warning("BAUs are Polygons and not Pixels. Currently BAU of identical
                           area are being assumed when computing the incidence matrix
                           from observations having a large support.
                           Handling of different areas will be catered for in a future revision.
                           Please report this issue to the package maintainer.")
                  #data_sp$id <- rownames(data_sp@data)
                  data_sp$id <- row.names(data_sp)
                  BAU_as_points <- SpatialPointsDataFrame(coordinates(sp_pols),sp_pols@data,
                                                          proj4string = CRS(proj4string(sp_pols)))
                  BAUs_aux_data <- over(data_sp,BAU_as_points)
                  ## The covariates are not averaged using this method... only the covariate in the last BAU
                  ## is recorded. In the following we find the average covariate over support
                  ## (needs to be done separately as my have overlapping observations)
                  for (i in 1L:length(data_sp)) {
                      #extract poly (a bit long-winded)
                      this_poly <- SpatialPolygons(list(data_sp@polygons[[i]]),1L,
                                                   proj4string=CRS(proj4string(data_sp)))
                      # find which points overlap observations
                      overlap <- which(over(BAU_as_points,this_poly) == 1)
                      # extract BAU data at these points
                      BAU_data <- BAU_as_points[names(overlap),1:ncol(BAU_as_points)]
                      this_attr <- data.frame(t(apply(BAU_data@data,2,mean))) # average over BAU data
                      # only columns not already in data so that they cannot be written over
                      BAUs_aux_data[data_sp[["id"]][i],] <- this_attr # assign to data
                  }

                  stopifnot(all(row.names(BAUs_aux_data) == row.names(data_sp)))
                  col_sel <- which(!(names(BAUs_aux_data) %in% names(data_sp)))

                  data_sp@data <- cbind(data_sp@data,BAUs_aux_data[,col_sel])
                  data_sp$Nobs <- 1
                  new_sp_pts <- data_sp
              }

              new_sp_pts
          })

setMethod("map_data_to_BAUs",signature(data_sp="ST"),
          function(data_sp,sp_pols,av_var,average_in_BAU = TRUE) {

              sp_fields <- NULL

              if(is(data_sp,"STIDF")) {
                  data_all_spatial <- as(data_sp,"Spatial")
              } else {
                  data_sp2 <- as(data_sp,"STIDF")
                  data_all_spatial <- as(data_sp2,"Spatial")
              }
              if(all(class(data_all_spatial$time) == "Date")) {
                  data_all_spatial$time <- as.POSIXlt( data_all_spatial$time)
              }

              ## Bin every spatial frame separately
              sp_fields <- lapply(seq_along(sp_pols@time),
                                  function(i) {
                                      t1 <- time(sp_pols)[i]
                                      if(i < last(sp_pols@time)) {
                                          t2 <- time(sp_pols)[i+1]
                                          trange <- paste0(format(t1),"::",format(t2-1))
                                      } else {
                                          trange <- paste0(format(t1),"::",format(last(data_sp@endTime)))
                                          t2 <- format(last(sp_pols@endTime))
                                      }
                                      data_spatial <- subset(data_all_spatial, time >= t1 &
                                                                                   time < t2)
                                      if(nrow(data_spatial) > 0) { ## If at least one point in this time period
                                          if(is(data_sp,"STFDF")) {
                                              ## Only subset polygons if there is at least one equivalent point
                                              data_spatial <- data_sp[,trange]
                                          }

                                          BAU_spatial <- sp_pols[,i]
                                          BAU_spatial@data <- filter(BAU_spatial@data,time == t1)
                                          BAU_spatial@data <- cbind(BAU_spatial@data,coordinates(BAU_spatial))
                                          map_data_to_BAUs(data_spatial,
                                                           BAU_spatial,
                                                           av_var=av_var,
                                                           average_in_BAU = average_in_BAU)
                                          } else {
                                              NULL
                                          }})

              if(is(data_sp,"STIDF")) {
                  ## Recast into a STIDF
                  time <- sp <- n <- NULL
                  for(i in seq_along(sp_pols@time)) {
                      if(!is.null(sp_fields[[i]])) {
                          sp <- rbind(sp,sp_fields[[i]]@data)
                          n <- nrow(sp_fields[[i]])
                          this_time <- rep(time(sp_pols)[i],n)
                          if(is.null(time)) time <- this_time else time <- c(time,this_time)
                          coordlabels <- coordnames(sp_fields[[i]]) # Ensures labels are from non-null field
                      }
                  }

                  coordinates(sp) <- coordlabels
                  sp@data <- cbind(sp@data,coordinates(sp))

                  STIDF(as(sp,"SpatialPoints"),
                        time,
                        data = sp@data)
              } else {
                  time <- sp <- n <- NULL
                  for(i in seq_along(sp_pols@time)) {
                      if(!is.null(sp_fields[[i]])) {
                          sp <- rbind(sp,sp_fields[[i]]@data)
                          n <- nrow(sp_fields[[i]])
                          this_time <- rep(time(sp_pols)[i],n)
                          if(is.null(time)) time <- this_time else time <- c(time,this_time)
                          coordlabels <- coordnames(sp_fields[[i]]) # Ensures labels are from non-null field
                      }
                  }
                  data_sp@data <- sp
                  data_sp

              }

          })

est_obs_error <- function(sp_pts,variogram.formula,vgm_model = NULL,BAU_width = NULL) {

    #stopifnot(is(variogram.formula,"formula"))
    stopifnot(is(sp_pts,"Spatial"))
    if(!("Nobs" %in% names(sp_pts)))
        stop("Nobs (number of observations in grid cell) needs to be a field of the Spatial object")
    if(!requireNamespace("gstat"))
        stop("gstat is required for variogram estimation. Please install gstat")
    if(!is.na(proj4string(sp_pts))) { ## Make sure we're not on sphere, otherwise variogram fitting is slow
        sp_pts <- SpatialPointsDataFrame(coords=coordinates(sp_pts),
                                         data = sp_pts@data,proj4string = CRS())
    }
    if(length(sp_pts) > 50000) {
        print("Selecting 50000 data points at random for estimating the measurement error variance")
        sp_pts_sub <- sp_pts[sample(1:length(sp_pts),50000),]
    } else {
        sp_pts_sub <- sp_pts
    }

    diag_length <- sqrt(sum(apply(coordinates(sp_pts_sub),2,
                                  function(x) diff(range(x)))^2))
    area <- prod(apply(coordinates(sp_pts_sub),2,
                       function(x) diff(range(x))))  ## full area
    cutoff <- sqrt(area * 100 / length(sp_pts_sub)) ## consider the area that contains about 100 data points in it

    print("... Fitting variogram for estimating measurement error")
    L <- .gstat.formula(variogram.formula,data=sp_pts_sub)
    g <- gstat::gstat(formula=variogram.formula,data=sp_pts_sub)
    v <- gstat::variogram(g,cressie=T,cutoff=cutoff,cutoff/10)

    if(is.null(vgm_model))
        vgm_model <-  gstat::vgm(var(L$y)/2, "Lin", mean(v$dist), var(L$y)/2)

    OK <- tryCatch({vgm.fit <- gstat::fit.variogram(v, model = vgm_model); 1},warning=function(w) 0)
    vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))

    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
        ## Try with line on first four points
        linfit <- lm(gamma~dist,data=v[1:4,])
        vgm.fit$psill[1] <- coefficients(linfit)[1]
    }

    if(vgm.fit$psill[1] <= 0) {
        ## Try with exponential
        vgm_model <-  gstat::vgm(var(L$y)/2, "Exp", mean(v$dist), var(L$y)/2)
        OK <- tryCatch({vgm.fit = gstat::fit.variogram(v, model = vgm_model); OK <- 1},warning=function(w) 0)
        vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))
    }

    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
        ## Try with Gaussian, maybe process is very smooth or data has a large support
        vgm_model <-  gstat::vgm(var(L$y)/2, "Gau", mean(v$dist), var(L$y)/2)
        OK <- tryCatch({vgm.fit = gstat::fit.variogram(v, model = vgm_model); OK <- 1},warning=function(w) 0)
        vgm.fit <- suppressWarnings(gstat::fit.variogram(v, model = vgm_model))
    }

    if(vgm.fit$psill[1] <= 0 | attributes(vgm.fit)$singular | !OK) {
        ## Just take the first point and print warning
        vgm.fit$psill[1] <- v$gamma[1]
        warning("Estimate of measurement error is probably inaccurate.
                 Please consider setting it through the std variable in the data object if known.")
    }

    #plot(v,vgm.fit)
    print(paste0("sigma2e estimate = ",vgm.fit$psill[1]))

    # sp_pts$std <- sqrt(vgm.fit$psill[1] / sp_pts$Nobs)
    sp_pts$std <- sqrt(vgm.fit$psill[1])  ## Assume observations have correlated error in BAU
    ##Observational error estimation could be improved. Currently a variogram is fitted to the data, and then the error variance of a single observation is assumed to be the partial sill. Then the variance of the averaged observations in the BAU is divided by Nobs. Currently there is no accounting for multiple data in the same grid box during variogram fitting as it's not straightforward with gstat

    sp_pts

}


setMethod("BuildC",signature(data="SpatialPolygons"),
          function(data,BAUs) {
              data$id <- 1:length(data)
              BAU_as_points <- SpatialPoints(coordinates(BAUs))
              i_idx <- j_idx <-  NULL
              for (i in 1L:length(data)) {
                  this_poly <- SpatialPolygons(list(data@polygons[[i]]),1L) #extract poly (a bit long-winded)
                  overlap <- which(over(BAU_as_points,this_poly) == 1) # find which points overlap observations
                  i_idx <- c(i_idx,rep(i,length(overlap)))
                  j_idx <- c(j_idx,as.numeric(overlap))
              }
              if(any(is.na(j_idx))) stop("NAs when constructing observation from
                                         large support observations. Are you sure all
                                         observations are covered by BAUs?")

              list(i_idx=i_idx,j_idx=j_idx)
          })

setMethod("BuildC",signature(data="SpatialPointsDataFrame"),
          function(data,BAUs) {
              BAU_index <- data.frame(row.names=row.names(BAUs),n =1:length(BAUs))
              i_idx <- 1:length(data)
              j_idx <- BAU_index[data$BAU_name,]
              ## Deprecated: #j_idx <- which(row.names(BAUs) %in% data$BAU_name)
              list(i_idx=i_idx,j_idx=j_idx)
          })

setMethod("BuildC",signature(data="STIDF"),
          function(data,BAUs) {
              i_idx <- 1:length(data)
              j_idx <- BAUs@data$n[data@data$n]
              ## Deprecated: #j_idx <- which(BAUs@data$n %in% data@data$n)
              list(i_idx=i_idx,j_idx=j_idx)
          })

setMethod("BuildC",signature(data="STFDF"),
          function(data,BAUs) {

              i_idx <- j_idx <- NULL
              count <- 0L
              C_one_time <- BuildC(data[,1],BAUs[,1])
              i <- C_one_time$i_idx
              j <- C_one_time$j_idx
              for(k in seq_along(.time.ST(BAUs))) {
                  overlap_time <- which(as.POSIXct(.time.ST(data)) ==
                                            (.time.ST(BAUs)[k]))
                  if(length(overlap_time) > 1L) stop("Something is wrong in binning polygon data into BAUs")
                  if(length(overlap_time) == 1) {
                      t_idx <- as.numeric(BAUs@time[k])
                      j_idx <- c(j_idx, (t_idx-1)*nrow(BAUs) + j)
                      i_idx <- c(i_idx, count*nrow(data) + i)
                      count <- count + 1
                  }
              }
              list(i_idx=i_idx,j_idx=j_idx)
          })


setMethod("coordinates",signature(obj="SpatialPolygons"),function(obj){
    coord_vals <- t(sapply(1:length(obj),function(i) obj@polygons[[i]]@Polygons[[1]]@labpt))
    colnames(coord_vals) <- colnames(obj@polygons[[1]]@Polygons[[1]]@coords)
    coord_vals
})


#' @aliases dimensions,measure-method
setMethod("dimensions",signature("measure"),function(obj){obj@dim})

#' @aliases dimensions,manifold-method
setMethod("dimensions",signature("manifold"),function(obj){dimensions(obj@measure)})

#' @aliases dimensions,Basis-method
setMethod("dimensions",signature("Basis"),function(obj){dimensions(obj@manifold)})


#' @rdname distance
#' @aliases distance,measure-method
setMethod("distance",signature("measure"),function(d,x1,x2=NULL){d@dist(x1,x2)})

#' @rdname distance
#' @aliases distance,manifold-method
setMethod("distance",signature("manifold"),function(d,x1,x2=NULL){distance(d@measure,x1,x2)})


#' @rdname type
#' @aliases type,manifold-method
setMethod("type",signature(.Object="manifold"),function(.Object) {
    return(.Object@type)
})

#' @rdname manifold
#' @aliases manifold,Basis-method
setMethod("manifold",signature(.Object="Basis"),function(.Object) {
    return(.Object@manifold)
})

#' @rdname manifold
#' @aliases manifold,TensorP_Basis-method
setMethod("manifold",signature(.Object="TensorP_Basis"),function(.Object) {
    return(list(manifold(.Object@Basis1),
                manifold(.Object@Basis2)))
})


#' @aliases coordnames,STFDF-method
setMethod("coordnames",signature(x="STFDF"),function(x) {
    return(coordnames(x@sp))
})

#' @aliases coordnames,STIDF-method
setMethod("coordnames",signature(x="STIDF"),function(x) {
    return(c(coordnames(x@sp),"t"))
})

dist_sphere <- function (x1, x2 = NULL, R = NULL)
    {
            # if R is null set to radius of Earth
            if (is.null(R)) R <- 6378.137
                if(is.null(x2)) x2 <- x1

                # Convert to radians
                x1 <- x1 * pi/180
                x2 <- x2 * pi/180

                # Formula from https://en.wikipedia.org/wiki/Great-circle_distance
                # d = r.acos(n1.n2) where n1 and n2 are the normals to the ellipsoid at the two positions
                n1 <- cbind(cos(x1[, 2]) * cos(x1[, 1]), cos(x1[, 2]) * sin(x1[, 1]), sin(x1[, 2]))
                n2 <- cbind(cos(x2[, 2]) * cos(x2[, 1]), cos(x2[, 2]) * sin(x2[, 1]), sin(x2[, 2]))
                delta <- sigma <- tcrossprod(n1,n2)

                ## Clamp to one
                return(R * acos(ifelse(abs(delta <- sigma) > 1,
                                           sign(delta <- sigma),
                                           delta <- sigma)))

}



load_dggrids <- function (res = 3L){
    if(!is.numeric(res))
        stop("res needs to be an integer or vector of integers")
    isea3h <- NA # suppress binding warning
    if(all(res <= 6L))  {
        data(isea3h, envir=environment(),package="FRK")
    } else {
        if(!requireNamespace("dggrids")) {
            stop("Such high DGGRID resolutions are not
                                       shipped with the package FRK. For this
                                       resolution please download and install the
                                       package dggrids from https://github.com/andrewzm/dggrids")
        } else {
            data(isea3h,envir=environment(),package = "dggrids")
        }
    }
    return(isea3h)
}

.extract.from.formula <- function (formula, data)
{
    m = model.frame(terms(formula), as(data, "data.frame"),
                           na.action = na.fail)
    Y = model.extract(m, "response")
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")
    grid = numeric(0)
    xlevels = .getXlevels(Terms, m)
    list(y = Y, locations = coordinates(data), X = X, call = call,
         has.intercept = has.intercept, grid = as.double(unlist(grid)),
         xlevels = xlevels)
}



process_isea3h <- function(isea3h,resl) {
    ## Splits the polygons at the 180 boundary
    ## Algorithm taken from
    ## https://stat.ethz.ch/pipermail/r-sig-geo/2015-July/023168.html

    res <- lon <- probpoly <- centroid <- lat <- NULL # suppress bindings warning

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


.prec_from_neighb <- function (neighb, intrinsic = 1, precinc = 1)
{
    num_v <- length(neighb)
    num_neighb <- lapply(neighb, length)
    if (intrinsic == 1) {
        i_list <- vector("list", num_v)
        for (k in 1:num_v) {
            i_list[[k]] <- rep(k, num_neighb[[k]])
        }
        i <- unlist(i_list)
        j <- unlist(neighb)
        z <- rep(-1, length(j))
        i <- c(i, 1:num_v)
        j <- c(j, 1:num_v)
        zdiag <- unlist(num_neighb)
        z <- c(z, zdiag)
    }
    if (intrinsic == 2) {
        i1 <- 1:num_v
        j1 <- 1:num_v
        z1 <- rep(0, num_v)
        for (k in 1:num_v) {
            z1[k] <- num_neighb[[k]]^2 + num_neighb[[k]]
        }
        count <- 1
        i2 <- rep(0, num_v * 10)
        j2 <- rep(0, num_v * 10)
        z2 <- rep(0, num_v * 10)
        for (k in 1:num_v) {
            for (l in neighb[[k]]) {
                i2[count] <- k
                j2[count] <- l
                z2[count] <- -(num_neighb[[k]] + num_neighb[[l]] -
                                   sum(duplicated(c(neighb[[l]], neighb[[k]]))))
                count <- count + 1
            }
        }
        i2 <- i2[1:count - 1]
        j2 <- j2[1:count - 1]
        z2 <- z2[1:count - 1]
        count <- 1
        i3 <- rep(0, num_v * 15)
        j3 <- rep(0, num_v * 15)
        z3 <- rep(0, num_v * 15)
        neighb2 <- vector("list", num_v)
        for (k in 1:num_v) {
            for (l in neighb[[k]]) {
                neighb2[[k]] <- c(neighb2[[k]], setdiff(neighb[[l]],
                                                        c(neighb[[k]], k)))
            }
            for (l in unique(neighb2[[k]])) {
                i3[count] <- k
                j3[count] <- l
                z3[count] <- sum(neighb2[[k]] == l)
                count <- count + 1
            }
        }
        i3 <- i3[1:count - 1]
        j3 <- j3[1:count - 1]
        z3 <- z3[1:count - 1]
        i <- c(i1, i2, i3)
        j <- c(j1, j2, j3)
        z <- c(z1, z2, z3)
    }
    z <- precinc * z
    Q <- sparseMatrix(i, j, x = z)
    return(Q)
}

## exactly as as.SpatialPolygons.GridTopology2 but with correct names to avoid name switch
as.SpatialPolygons.GridTopology2 <- function (grd, proj4string = CRS(as.character(NA)))
{
    coord_names <- names(grd@cellsize)
    grd_crds <- coordinates(grd)
    IDs <- IDvaluesGridTopology(grd)
    nPolygons <- nrow(grd_crds)
    cS <- grd@cellsize
    cS2 <- cS/2
    cS2x <- cS2[1]
    cS2y <- cS2[2]
    Srl <- vector(mode = "list", length = nPolygons)
    xi <- grd_crds[,1]
    yi <- grd_crds[,2]
    xall <- cbind(xi - cS2x, xi - cS2x, xi + cS2x, xi + cS2x, xi -
               cS2x)
    yall <- cbind(yi - cS2y, yi + cS2y, yi + cS2y, yi - cS2y, yi -
               cS2y)
    for (i in 1:nPolygons) {
        coords <- cbind(xall[i,],yall[i,])
        colnames(coords) <- coord_names
        Srl[[i]] <- Polygons(list(Polygon(coords = coords)),
                             ID = IDs[i])
        comment(Srl[[i]]) <- "0"
    }
    res <- SpatialPolygons(Srl, proj4string = proj4string)
    res
}

.choose_manifold_from_data <- function(data) {

    stopifnot(is(data,"Spatial") | is(data,"ST"))
    if(is(data, "Spatial")) {
        p4 <- proj4string(data)
        if(is.na(p4)) {
            manifold = plane()
        }  else {
            if(grepl("longlat",p4)) {
                manifold = sphere()
            } else {
                manifold = plane()
            }
        }
    } else {
        p4 <- proj4string(data@sp)
        if(is.na(p4)) {
            manifold = STplane()
        }  else {
            if(grepl(p4,"longlat")) {
                manifold = STsphere()
            } else {
                manifold = STplane()
            }
        }
    }
    manifold
}

.choose_BAU_cellsize_from_data <- function(data) {
    cellsize <- c(diff(range(coordinates(data)[,1]))/100,
                  diff(range(coordinates(data)[,2]))/100)
    if (is(data,"Spatial")) {
        cellsize
    } else {
       c(cellsize,1)
    }
}

.choose_BAU_tunit_from_data <- function(data) {
    # Aim for no more than 30 BAUs
    t1 <- range(.time.ST(data))[1]
    t2 <- range(.time.ST(data))[2]
    tunits <- c("days","weeks","months","years")
    for(i in seq_along(tunits)) {
        l <- length(seq(t1,t2,by=tunits[i]))
        if(l < 40) break
    }
    tunits[i]
}

.polygons_to_points <- function(polys) {
    stopifnot(is(polys,"STFDF") |  is(polys,"SpatialPixels")| is(polys,"SpatialPolygons"))
    if(is(polys,"STFDF")) {
        as.matrix(cbind(coordinates(polys),polys@data$t))
    } else {
        #as.matrix(polys[coordnames(polys)]@data)
        as.matrix(coordinates(polys))
    }

}

.find_hull <- function(coords,use_INLA=TRUE,convex = -0.05) {
    if(use_INLA) {
        bndary_seg = INLA::inla.nonconvex.hull(coords,convex=convex)$loc
    } else {
        conv_hull <- coordinates(coords)[chull(coordinates(coords)),]
        bound_box <- bbox(coords)
        centroid <- apply(bound_box,1,mean)
        # expand convex hull out
        bndary_seg <- conv_hull
        delta <- max(apply(bound_box,1,function(x) diff(range(x))))
        bndary_seg[,1] <- conv_hull[,1] + sign(conv_hull[,1] - centroid[1])*delta/5
        bndary_seg[,2] <- conv_hull[,2] + sign(conv_hull[,2] - centroid[2])*delta/5
        #bndary_seg[,1] <- centroid[1] + (conv_hull[,1] - centroid[1])*2
        #bndary_seg[,2] <- centroid[2] + (conv_hull[,2] - centroid[2])*2
        colnames(bndary_seg) <- NULL
    }
    bndary_seg
}

.time.ST <- function (x, ...)
    index(x@time)
