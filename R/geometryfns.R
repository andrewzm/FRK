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
#' @param convex convex parameter used for smoothing an extended boundary when working on a finite domain (that is, when the object \code{d} is supplied), see details.
#' @param tunit temporal unit when requiring space-time BAUs. Can be either "secs", "mins", "hours" or "days".
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
auto_BAUs <- function(manifold, type="grid",cellsize = rep(1,dimensions(manifold)),
                      isea3h_res=NULL,data=NULL,convex=-0.05,tunit="days") {
    if(!(is(data,"Spatial") | is(data,"ST") | is(data,"Date") | is.null(data)))
        stop("Data needs to be of class 'Spatial', 'ST', 'Date', or NULL")
    if(!is(manifold,"manifold")) stop("manifold needs to be of class 'manifold'")
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
    if(length(cellsize) == 1) cellsize <- rep(cellsize,dimensions(manifold))
    if(!length(cellsize) == dimensions(manifold)) stop("cellsize needs to be of length equal to dimension of manifold")

    auto_BAU(manifold=manifold,type=type,cellsize=cellsize,resl=resl,d=data,convex=convex,tunit=tunit)
}


setMethod("auto_BAU",signature(manifold="plane"),
          function(manifold,type="grid",cellsize = c(1,1),resl=resl,d=NULL,convex=-0.05,...) {

               if(!requireNamespace("INLA"))
                   stop("For automatic BAU generation INLA needs to be installed for
                        constructing non-convex hull. Please install it using
                        install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/stable\")")

              X1 <- X2 <- NULL # Suppress bindings warning

              if(is(d,"SpatialPoints")){
                coords <- coordinates(d)
              } else if(is(d,"SpatialPolygons")){
                 ## get out all edges
                 coords <- do.call("rbind",
                         lapply(1:length(d),
                                function(i) coordinates(d@polygons[[i]]@Polygons[[1]])))
              }
              xrange <- range(coords[,1])
              yrange <- range(coords[,2])
              bndary_seg = INLA::inla.nonconvex.hull(coords,convex=convex)$loc %>%
                  data.frame() %>%
                  mutate(x=X1,y=X2,id = 1) %>%
                  select(-X1,-X2) %>%
                  df_to_SpatialPolygons(keys="id",coords=c("x","y"),proj=CRS())

              drangex <- diff(xrange)
              drangey <- diff(yrange)
              xgrid <- seq(xrange[1] - drangex*1.2,xrange[2] + drangex*1.2,by=cellsize[1])
              ygrid <- seq(yrange[1] - drangey*1.2,yrange[2] + drangey*1.2,by=cellsize[2])

              xy <- expand.grid(x=xgrid,y=ygrid)  %>%
                  SpatialPoints()

              if(type == "hex") {
                  HexPts <- spsample(xy,type="hexagonal",cellsize = cellsize[1])
                  idx <- which(!is.na(over(HexPts,bndary_seg)))
                  HexPols <- HexPoints2SpatialPolygons(HexPts[idx,])
                  coordnames(HexPols) <- coordnames(d)
                  coordnames(HexPts) <- coordnames(d)
                  HexPols_df <- SpatialPolygonsDataFrame(HexPols,
                                                         data.frame(
                                                             coordinates(HexPts[idx,]),
                                                             row.names = row.names(HexPols)))
                  return(HexPols_df)
              } else if (type == "grid") {
                  xy <- xy %>%
                      points2grid() %>%
                      as.SpatialPolygons.GridTopology()
                  if(!all(coordnames(xy) == coordnames(d))) {
                      warning("Coordinate names different from (x,y).
                              Renaming polygons might take a while due to
                              structure of the function as.SpatialPolygons.GridTopology.")
                      coordnames(xy) <- coordnames(d)
                  }
                  idx <- which(!is.na(over(xy,bndary_seg)))
                  xy <- xy[idx,]
                  xy_df <- SpatialPolygonsDataFrame(xy,
                                                    data.frame(
                                                        coordinates(xy),
                                                        row.names = row.names(xy)))
                  return(xy_df)
              }

          })


setMethod("auto_BAU",signature(manifold="timeline"),
          function(manifold,type="grid",cellsize = c(1),resl=resl,d=NULL,convex=-0.05,...) {

              l <- list(...)

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
          function(manifold,type="grid",cellsize = c(1,1,1),resl=resl,d=NULL,convex=-0.05,...) {

              if(is(d,"ST")) {
                  space_part <- d@sp
                  time_part <- time(d)
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
                                       resl=resl,type=type,d=space_part,convex=convex,...)
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
    xgrid <- seq(xrange[1] - drangex*0.2,xrange[2] + drangex*0.2,by=cellsize[1]) %>%
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

        isea3h_sp_pol <- df_to_SpatialPolygons(
            df=filter(isea3h_res,centroid==0),
            keys=c("id"),
            coords=c("lon","lat"),
            proj=CRS("+proj=longlat +ellps=sphere"))

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
        lonlat <- expand.grid(lon=longrid,lat=latgrid) %>%
            SpatialPoints() %>%
            points2grid() %>%
            as.SpatialPolygons.GridTopology(proj4string = CRS("+proj=longlat +ellps=sphere"))

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
                                  proj=CRS("+proj=longlat +ellps=sphere"))
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
#' @details A 2D plane is initialised using a \code{measure} object. By default, the measure object (\code{measure}) is the Euclidean distance in 2 dimensions, \link{Euclid_dist}. Distances are computed using functions from the \code{fields} package.
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
#' @aliases Euclid_dist
#' @aliases gc_dist
#' @aliases gc_dist_time
#' @title Pre-configured distances
#'
#' @description Useful objects of class \code{distance} included in package.
#'
#' @param dim dimension of Euclidian space
#' @param R great-circle radius
#' @details Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points.  Currently the Euclidean distance and the great-circle distance are included. Distances are computed using functions extracted from the package \code{fields}.
#' @export
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  rdist(x1,x2), dim=dim)
}

#' @rdname distances
#' @export
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  rdist.earth(x1,x2,miles=F,R=R),dim=2L)
}

#' @rdname distances
#' @export
gc_dist_time <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  {
        spatdist <- rdist.earth(x1[,1:2,drop=FALSE],x2[,1:2,drop=FALSE],miles=F,R=R)
        tdist <- rdist(x1[,3],x2[,3])
        sqrt(spatdist^2 + tdist^2) } ,dim=3L)
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

    if(opts_FRK$get("parallel") > 1) {
        cl <- makeCluster(opts_FRK$get("parallel"))

        ## Deprecated to remove plyr:
        #doParallel::registerDoParallel(opts_FRK$get("parallel"))
        #df_poly <- plyr::dlply(df,keys,dfun,.parallel=TRUE)

        df_poly <- mclapply(unique(df[keys])[,1],
                            function(key) {
                                dfun(df[df[keys]==key,])},
                            mc.cores = opts_FRK$get("parallel"))
        stopCluster(cl)
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


#' @aliases map_data_to_BAUs,Spatial-method
setMethod("map_data_to_BAUs",signature(data_sp="Spatial"),
          function(data_sp,sp_pols,av_var,average_in_BAU = TRUE)
          {
              if(is(data_sp,"SpatialPointsDataFrame")) {

                  if(average_in_BAU) {
                      Nobs <- NULL
                      data_sp$Nobs <- 1

                      timer <- system.time(Data_in_BAU <- over(sp_pols,data_sp[c(av_var,"Nobs","std")],fn=sum))

                      ## Rhipe VERSION (Currently disabled)
                      # print("Using RHIPE to find overlays")
                      # timer <- system.time(
                      # Data_in_BAU <- rhwrapper(Ntot = length(sp_pols),
                      #                                  N = 4000,
                      #                                  f_expr = .rhover,
                      #                                  sp_pols = sp_pols,
                      #                                  data_sp = data_sp,
                      #                                  av_var=av_var)
                      #     )


                      sp_pols@data[av_var] <- Data_in_BAU[av_var]/Data_in_BAU$Nobs
                      sp_pols@data["std"] <- Data_in_BAU["std"]/Data_in_BAU$Nobs
                      sp_pols@data["Nobs"] <- Data_in_BAU$Nobs
                      sp_pols@data["BAU_name"] <- as.character(row.names(sp_pols))

                      new_sp_pts <- SpatialPointsDataFrame(coords=sp_pols[coordnames(data_sp)]@data,
                                                           data=sp_pols@data,
                                                           proj4string = CRS(proj4string(data_sp)))
                      #new_sp_pts$std <- sqrt(1 / new_sp_pts$Nobs)
                      new_sp_pts$std <- sqrt(new_sp_pts$std^2 / new_sp_pts$Nobs)
                      new_sp_pts <- subset(new_sp_pts,!is.na(Nobs))
                  } else {

                      timer <- system.time(Data_in_BAU <- over(data_sp,as(sp_pols,"SpatialPolygons")))
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
                  data_sp$id <- rownames(data_sp@data)
                  BAU_as_points <- SpatialPointsDataFrame(coordinates(sp_pols),sp_pols@data)
                  BAUs_aux_data <- over(data_sp,BAU_as_points)
                  ## The covariates are not averaged using this method... only the covariate in the last BAU
                  ## is recorded. In the following we find the average covariate over support
                  ## (needs to be done separately as my have overlapping observations)
                  for (i in 1L:length(data_sp)) {
                      this_poly <- SpatialPolygons(list(data_sp@polygons[[i]]),1L) #extract poly (a bit long-winded)
                      overlap <- which(over(BAU_as_points,this_poly) == 1) # find which points overlap observations
                      BAU_data <- BAU_as_points[names(overlap),1:ncol(BAU_as_points)] # extract BAU data at these points
                      this_attr <- data.frame(t(apply(BAU_data@data,2,mean))) # average over BAU data
                      BAUs_aux_data[data_sp[["id"]][i],] <- this_attr # assign to data
                  }

                  stopifnot(all(row.names(BAUs_aux_data) == row.names(data_sp)))
                  data_sp@data <- cbind(data_sp@data,BAUs_aux_data)
                  data_sp$Nobs <- 1
                  new_sp_pts <- data_sp
              }

              new_sp_pts
          })

setMethod("map_data_to_BAUs",signature(data_sp="ST"),
          function(data_sp,sp_pols,av_var,average_in_BAU = TRUE) {

              sp_fields <- NULL

              data_all_spatial <- as(data_sp,"Spatial")
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
                                      #data_spatial <- as(data_sp[,trange],"Spatial")
                                      data_spatial <- subset(data_all_spatial, time >= t1 &
                                                                               time < t2)
                                      BAU_spatial <- sp_pols[,i]
                                      BAU_spatial@data <- filter(BAU_spatial@data,time == t1)
                                      BAU_spatial@data <- cbind(BAU_spatial@data,coordinates(BAU_spatial))
                                      map_data_to_BAUs(data_spatial,
                                                       BAU_spatial,
                                                       av_var=av_var,
                                                       average_in_BAU = average_in_BAU)})


              if(is(data_sp@sp,"SpatialPolygons")) stop("ST not implemented for polygon observations yet.
                                                        Please contact package maintainer")

              ## Recast into a STIDF
              time <- sp <- n <- NULL
              for(i in seq_along(sp_pols@time)) {
                  sp <- rbind(sp,sp_fields[[i]]@data)
                  n <- nrow(sp_fields[[i]])
                  time[[i]] <- rep(time(sp_pols)[i],n)
              }
              coordinates(sp) <- coordnames(sp_fields[[1]])
              sp@data <- cbind(sp@data,coordinates(sp))
              time <- do.call("c",time)

              STIDF(as(sp,"SpatialPoints"),
                    time,
                    data = sp@data)

          })

est_obs_error <- function(sp_pts,variogram.formula) {

    #stopifnot(is(variogram.formula,"formula"))
    stopifnot(is(sp_pts,"Spatial"))
    if(!("Nobs" %in% names(sp_pts)))
        stop("Nobs (number of observations in grid cell) needs to be a field of the Spatial object")
    if(!requireNamespace("gstat"))
        stop("gstat is required for variogram estimation. Please install gstat")

    g <- gstat::gstat(formula=variogram.formula,data=sp_pts)
    v <- gstat::variogram(g,cressie=T)
    vgm.fit = gstat::fit.variogram(v, model = gstat::vgm(1, "Lin", mean(v$dist), 1))
    plot(v,vgm.fit)
    print(paste0("sigma2e estimate = ",vgm.fit$psill[1]))
    if(vgm.fit$psill[1] == 0)
        stop("Measurement error estimated to be zero. Please pre-specify measurement error. If
              unknown please specify a reasonable value in the field 'std' and set
              est_error = FALSE")
    sp_pts$std <- sqrt(vgm.fit$psill[1] / sp_pts$Nobs)

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

setMethod("coordinates",signature(obj="SpatialPolygons"),function(obj){
    coord_vals <- t(sapply(1:length(obj),function(i) obj@polygons[[i]]@Polygons[[1]]@labpt))
    colnames(coord_vals) <- colnames(obj@polygons[[1]]@Polygons[[1]]@coords)
    coord_vals
})


#' @aliases dimensions,measure-method
setMethod("dimensions",signature("measure"),function(obj){obj@dim})

#' @aliases dimensions,manifold-method
setMethod("dimensions",signature("manifold"),function(obj){dimensions(obj@measure)})

#' @rdname distance
#' @aliases distance,measure-method
setMethod("distance",signature("measure"),function(d,x1,x2){d@dist(x1,x2)})

#' @rdname distance
#' @aliases distance,manifold-method
setMethod("distance",signature("manifold"),function(d,x1,x2){distance(d@measure,x1,x2)})


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


## The functions rdist and rdist.earth were taken from the package fields
## fields is lincensed under GPL >=2
rdist <- function (x1, x2 = NULL, compact = FALSE)
{
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
    }
    if (is.null(x2)) {
        storage.mode(x1) <- "double"
        if (compact)
            return(dist(x1))
        else return(.Call("RdistC", x1, x1, PACKAGE = "FRK"))
    }
    else {
        if (!is.matrix(x2)) {
            x2 <- as.matrix(x2)
        }
        storage.mode(x1) <- "double"
        storage.mode(x2) <- "double"
        return(.Call("RdistC", x1, x2))
    }
}


rdist.earth <- function (x1, x2 = NULL, miles = TRUE, R = NULL)
{
    if (is.null(R)) {
        if (miles)
            R <- 3963.34
        else R <- 6378.388
    }
    coslat1 <- cos((x1[, 2] * pi)/180)
    sinlat1 <- sin((x1[, 2] * pi)/180)
    coslon1 <- cos((x1[, 1] * pi)/180)
    sinlon1 <- sin((x1[, 1] * pi)/180)
    if (is.null(x2)) {
        pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
            t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
        return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
    else {
        coslat2 <- cos((x2[, 2] * pi)/180)
        sinlat2 <- sin((x2[, 2] * pi)/180)
        coslon2 <- cos((x2[, 1] * pi)/180)
        sinlon2 <- sin((x2[, 1] * pi)/180)
        pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
            t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
        return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
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
                                       package dggrids from github/andrewzm")
        } else {
            data(isea3h,envir=environment(),package = "dggrids")
        }
    }
    return(isea3h)
}

.extract.from.formula <- function (formula, data)
{
    m = model.frame(terms(formula), as(data, "data.frame"), na.action = na.fail)
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

    isea3h_res2 <- isea3h_res %>%
        data.frame() %>%
        filter(probpoly == FALSE) %>%
        rbind(new_polys)

    isea3h_res2

}

