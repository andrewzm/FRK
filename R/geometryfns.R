#' @title manifold
#' @param .Object \code{manfold} object passed up from lower-level constructor
#' @description Manifold initialisation. This function should not be called directly as \code{manifold} is a virtual class.
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here
    .Object
})

#' @title Convert DGGRID gen. file to data frame
#'
#' @description Convert discrete global grid (DGGRID) .gen file to data frame
#' @param filename name of DGGRID output file of type .gen
#' @param res resolution of the DGGRID (this number is attached to the resulting data frame)
#'
#' @details DGGRID produces .gen files with three columns, \code{id}, \code{lon} and \code{lat}. Polygons each have a unique identifier and following each polygon there is a line "END" that needs to be filtered out. This function takes the .gen file and converts it into a data frame with fields \code{lon, lat, id} and \code{res}.
#'
#' @export
dggrid_gen_to_df <- function(filename,res) {
    if(!require(zoo)) stop("Please install package zoo") ## needed for na.locf
    lat <- lon <- isna <- NULL # Suppress bindings warning

    if(!is.character(filename)) stop("filename needs to of type 'character'")
    if(!require(zoo)) stop("zoo is needed for this function. Please install it separately")
    if(!is.numeric(res)) stop("res needs to be of type character")
    X <- read.table(filename,
                    sep=" ",fill=T,
                    header=F,
                    col.names = c("id","lon","lat")) %>%
        filter(!(id == "END")) %>%  ## Sometimes DGGRID doesn't put in first space, this mutate caters for this
        mutate(isna = is.na(lat),
               lat = ifelse(isna,lon,lat),
               lon = ifelse(isna,as.numeric(as.character(X$id)),lon))
    X$id[X$isna] <- NA
    X <- select(X,-isna) %>%
        mutate(res = res,
               id = as.numeric(as.character(id)),
               centroid = as.numeric(!is.na(id))) %>%
        transform(id = zoo::na.locf(id))
}

#' @title Automatic BAU generation
#' @description This function calls the generic function \code{auto_BAU} (currently not exported) after a series of checks and is the easiest way to generate a set of BAUs on the manifold being used.
#' @param manifold object of class \code{manifold}
#' @param res resolution number of the isea3h DGGRID cells for when type is ``hex'' and manifold is a \code{sphere}
#' @param cellsize denotes size of gridcell when \code{type} = ``grid''. Needs to be f length 1 (isotropic grid) or a vector of length \code{dimensions(manifold)}
#' @param type either ``hex'' or ``grid'', indicating whether gridded or hexagonal BAUs should be used
#' @param data object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}. Provision of \code{data} implies that the domain is bounded, and is thus necessary when the manifold is a \code{real_line} or a \code{plane} but is not necessary when the manifold is a \code{sphere}
#' @param convex convex parameter used for smoothing an extended boundary when working on a finite domain (that is, when the object \code{d} is supplied), see details.
#' @param tunit temporal unit when requiring space-time BAUs. Can be either "secs", "mins", "hours" or "days".
#' @details \code{auto_BAUs} constructs a set of basic aerial units (BAUs) used both for data pre-processing and for prediction. As such, the BAUs need to be of sufficienly fine resolution so that data is not adversely affected; a subset of BAUs (or a completely new set of polygons) may be used for prediction
#'
#' Two types of BAUs are supported by \code{FRK}: ``hex'' (hexagonal) and ``grid'' (rectangular). In order to have a ``grid'' set of BAUs, the user should specify a cellsize of length equal to the dimensions of the manifold, that is, of length 1 for \code{real_line} and 2 for \code{sphere} and \code{plane}. When a ``hex'' set of BAUs is desired, the first element of \code{cellsize} is used to determine the side length by dividing this value by approximately 2. If the manifold is a \code{sphere}, then \code{cellsize} is ignored and the resulting hexagonal grid is a discrete global grid with resolution \code{res}. The parameter type is ignored with \code{real_line} and ``hex'' is not available for this manifold.
#'
#'   If the object \code{data} is provided, then automatic domain selection is carried out by employing the \code{INLA} function \code{inla.nonconvex.hull}, which finds a (non-convex) hull surrounding the data points (or centroids of the data polygons). This domain is extended and smoothed using the \code{convex} parameter. The parameter \code{convex} should be negative, and a larger absolute value for \code{convex} results in a larger domain with smoother boundaries. Due to the dependency on hull construction, \code{INLA} needs to be installed in order to use this function unless BAUs on a sphere are desired (note that \code{INLA} is not available on CRAN at time of writing).
#'  @examples
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
#' if(require(INLA)) { # INLA is needed to find the non-convex hull
#'      GridPols_df <- auto_BAUs(manifold = plane(),
#'                              cellsize = 200,
#'                              type = "grid",
#'                              data = meuse,
#'                              convex=-0.05)
#'      plot(GridPols_df)
#'      HexPols_df <- auto_BAUs(manifold = plane(),
#'                             cellsize = 200,
#'                             type = "hex",
#'                             data = meuse,
#'                             convex=-0.05)
#'      plot(HexPols_df)
#'  }
#'  @export
auto_BAUs <- function(manifold,res=2,cellsize = rep(1,dimensions(manifold)), type="hex",data=NULL,convex=-0.05,tunit="days") {
    if(!(is(data,"Spatial") | is(data,"ST") | is.null(data)))
        stop("Data needs to be of class 'Spatial', 'ST', or NULL")
    if(!is(manifold,"manifold")) stop("manifold needs to be of class 'manifold'")
    if(!is.numeric(res) | is.integer(res)) stop("res needs to be of type 'numeric' or 'integer'")
    if(length(cellsize) == 1) cellsize <- rep(cellsize,dimensions(manifold))
    if(!length(cellsize) == dimensions(manifold)) stop("cellsize needs to be of length equal to dimension of manifold")
    if(!(res >=0 & res <= 9)) stop("res needs to be between 0 and 9")
    resl <- round(res)

    auto_BAU(manifold=manifold,resl=resl,cellsize=cellsize,type=type,d=data,convex=convex,tunit=tunit)

}


setMethod("auto_BAU",signature(manifold="plane"),
          function(manifold,cellsize = c(1,1),resl=resl,type="hex",d=NULL,convex=-0.05,...) {

              X1 <- X2 <- NULL # Suppress bindings warning

              coords <- coordinates(d)
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
                  coordnames(xy) <- coordnames(d)
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
          function(manifold,cellsize = c(1),resl=resl,type="grid",d=NULL,convex=-0.05,...) {

              l <- list(...)

              if(!"tunit" %in% names(l)) stop("Need to supply argument tunit with value secs, mins, hours, days, months or years")

              tunit <- l$tunit
              tpoints <- time(d)

              if(is(tpoints,"Date"))
                  tpoints <- as.POSIXct(tpoints)

              trange <- range(tpoints)
              dranget <- diff(trange)

              tspacing <- paste(cellsize,tunit)
              tgrid <- seq(trunc(trange[1],tunit),
                           trunc(trange[2],tunit),
                           by=tspacing) %>%
                               # switch(tunit,
                               #       secs    = cellsize,
                               #       mins    = cellsize*60,
                               #       hours   = cellsize*3600,
                               #       days    = cellsize*3600*24,
                               #       months  = "mon",
                               #       years   = "year")) %>%
                  round(tunit)

              attr(tgrid,"tzone") <- attr(d@time,"tzone")
              return(tgrid)
          })



setMethod("auto_BAU",signature(manifold = c("STmanifold")),
          function(manifold,cellsize = c(1,1,1),resl=resl,type="hex",d=NULL,convex=-0.05,...) {

              if(is(manifold,"STplane")) {
                  spat_manifold <- plane()
              } else {
                  spat_manifold <- sphere()
              }

              spatial_BAUs <- auto_BAU(manifold=spat_manifold,cellsize=cellsize[1:2],
                                       resl=resl,type=type,d=d@sp,convex=convex,...)
              temporal_BAUs <- auto_BAU(manifold=timeline(), cellsize=cellsize[3],
                                        resl=resl,type=type,d=d,convex=convex,...)

              nt <- length(temporal_BAUs)
              ns <- nrow(spatial_BAUs)
              STBAUs <- STFDF(spatial_BAUs,
                              temporal_BAUs,
                              data = data.frame(n = 1:(nt *ns),
                                                time = rep(temporal_BAUs,each=ns),
                                                t = rep(1:nt,each=ns)))
              return(STBAUs)

          })

# setMethod("auto_BAU",signature(manifold="STsphere"),
#           function(manifold,cellsize = c(1,1,1),resl=resl,type="hex",d=NULL,convex=-0.05,...) {
#
#               spatial_BAUs <- auto_BAU(manifold=sphere(),cellsize=cellsize[1:2],
#                                        resl=resl,type=type,d=STobj@sp,convex=convex,...)
#               temporal_BAUs <- auto_BAU(manifold=sphere(), cellsize=cellsize[3],
#                                         resl=resl,type=type,d=d,convex=convex,...)
#
#               nt <- length(temporal_BAUs)
#               ns <- nrow(spatial_BAUs)
#               STBAUs <- STFDF(spatial_BAUs,
#                               temporal_BAUs,
#                               data = data.frame(n = 1:(nt *ns),
#                                                 time = rep(temporal_BAUs,each=ns),
#                                                 t = rep(1:nt,each=ns)))
#
#               return(STBAUs)
#           })

setMethod("auto_BAU",signature(manifold="real_line"),function(manifold,cellsize = 1,resl=resl,type="grid",d=NULL,...) {

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


setMethod("auto_BAU",signature(manifold="sphere"),function(manifold,cellsize = c(1,1),resl=2,type="hex",d=NULL,...) {
    if(type == "hex") {
        isea3h <- res <- lon <- centroid <- lat <- in_chull <- NULL # Suppress bindings warnings
        data(isea3h,envir = environment())
        isea3h_res <- filter(isea3h,res == resl) %>%
            arrange(id) %>%
            group_by(id) %>%
            filter(diff(range(lon)) < 90) %>%
            data.frame()

        isea3h_sp_pol <- df_to_SpatialPolygons(
            df=filter(isea3h_res,centroid==0),
            keys=c("id"),
            coords=c("lon","lat"),
            proj=CRS("+proj=longlat"))

        isea3h_sp_poldf <- SpatialPolygonsDataFrame(
            isea3h_sp_pol,
            cbind(data.frame(row.names=names(isea3h_sp_pol)),
                  (filter(isea3h_res,centroid==1) %>%
                       select(id,lon,lat))))
        sphere_BAUs <- isea3h_sp_poldf
    }  else if (type == "grid") {
        if(!is.null(d)) {
            coords <- coordinates(d)
            xrange <- range(coords[,1])
            yrange <- range(coords[,2])
            drangex <- diff(xrange)
            drangey <- diff(yrange)
            xmin <- xrange[1] - drangex*1.2
            xmax <- xrange[2] + drangex*1.2
            ymin <- yrange[1] - drangey*1.2
            ymax <- yrange[1] + drangey*1.2
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
            as.SpatialPolygons.GridTopology(proj4string = CRS("+proj=longlat"))

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
                                  proj=CRS("+proj=longlat"))
        sphere_BAUs$in_chull <- over(sphere_BAUs,conv_hull)
        sphere_BAUs <- subset(sphere_BAUs,!is.na(in_chull))
    }
    sphere_BAUs

})


#' @title sphere
#'
#' @description Sphere initialisation.
#'
#' @param radius radius of sphere
#'
#' @details A sphere is initialised using a \code{radius} parameter. By default, the radius \code{R} is equal to \code{R}=6371 km, the Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}. Distances are computed using functions from the \code{fields} package.
#' @export
#' @examples
#' S <- sphere()
#' print(sp::dimensions(S))
sphere <- function(radius=6371) {
    metric=gc_dist(R=radius)
    stopifnot(dimensions(metric)==2L)
    stopifnot(radius>0)
    new("sphere",metric=metric,radius=radius)
}

setMethod("initialize",signature="sphere",function(.Object,radius=1,metric=gc_dist(R=radius)) {
    .Object@type <- "sphere"
    .Object@metric <- metric
    .Object@radius <- radius
    callNextMethod(.Object)})


#' @title Space-time sphere
#'
#' @description Domain for modellng space-time fields over a sphere.
#'
#' @param radius radius of sphere
#'
#' @details As with the spatial-only sphere, the sphere is initialised using a \code{radius} parameter. By default, the radius \code{R} is equal to \code{R}=6371 km, the Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}. Distances are computed using functions from the \code{fields} package.
#' @export
#' @examples
#' S <- STsphere()
#' print(sp::dimensions(S))
STsphere <- function(radius=6371) {
    metric=gc_dist_time(R=radius)
    stopifnot(dimensions(metric)==3)
    stopifnot(radius>0)
    new("STsphere",metric=metric,radius=radius)
}

setMethod("initialize",signature="STsphere",function(.Object,radius=6371,metric=gc_dist_time(R=radius)) {
    .Object@type <- "STsphere"
    .Object@metric <- metric
    .Object@radius <- radius
    callNextMethod(.Object)})


#' @title plane
#'
#' @description Plane initialisation.
#'
#' @param metric an object of class \code{measure}
#'
#' @details A 2D plane is initialised using a \code{measure} object. By default, the measure object (\code{metric}) is the Euclidean distance in 2 dimensions, \link{Euclid_dist}. Distances are computed using functions from the \code{fields} package.
#' @export
#' @examples
#' P <- plane()
#' print(type(P))
#' print(sp::dimensions(P))
plane <- function(metric=Euclid_dist(dim=2L)) {
    stopifnot(dimensions(metric)==2L)
    new("plane",metric=metric)
}


setMethod("initialize",signature="plane",function(.Object,metric=Euclid_dist(dim=2L)) {
    .Object@type <- "plane"
    .Object@metric <- metric
    callNextMethod(.Object)})

#' @title plane in space-time
#'
#' @description Space-time plane initialisation.
#'
#' @param metric an object of class \code{measure}
#'
#' @details A 2D plane with a time component added is initialised using a \code{measure} object. By default, the measure object (\code{metric}) is the Euclidean distance in 3 dimensions, \link{Euclid_dist}. Distances are computed using functions from the \code{fields} package.
#' @export
#' @examples
#' P <- STplane()
#' print(type(P))
#' warning("CHECK: This may not work for ST plane")
#' print(sp::dimensions(P))
STplane <- function(metric=Euclid_dist(dim=3L)) {
    stopifnot(dimensions(metric)==3L)
    new("STplane",metric=metric)
}

setMethod("initialize",signature="STplane",function(.Object,metric=Euclid_dist(dim=3L)) {
    .Object@type <- "STplane"
    .Object@metric <- metric
    callNextMethod(.Object)})


#' @title timeline
#'
#' @description Timeline initialisation.
#'
#' @param metric an object of class \code{measure}
#'
#' @details A time axis initialised using a \code{measure} object. By default, the measure object (\code{metric}) is the absolute differences.
#' @export
#' @examples
#' P <- timeline()
#' print(type(P))
#' warning("CHECK: This may not work for ST plane")
#' print(sp::dimensions(P))
timeline <- function(metric=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(metric)==1L)
    new("timeline",metric=metric)
}

setMethod("initialize",signature="timeline",function(.Object,metric=Euclid_dist(dim=1L)) {
    .Object@type <- "timeline"
    .Object@metric <- metric
    callNextMethod(.Object)})


#' @title real line
#'
#' @description Real-line initialisation.
#'
#' @param metric an object of class \code{measure}
#'
#' @details A real line is initialised using a \code{measure} object. By default, the measure object (\code{metric}) describes the distance between two points as the absolute difference between the two coordinates. Distances are computed using functions from the \code{fields} package.
#' @export
#' @examples
#' R <- real_line()
#' print(type(R))
#' print(sp::dimensions(R))
real_line <- function(metric=Euclid_dist(dim=1L)) {
    stopifnot(dimensions(metric)==1L)
    new("real_line",metric=metric)
}

setMethod("initialize",signature="real_line",function(.Object,metric=Euclid_dist(dim=1L)) {
    .Object@type <- "real_line"
    .Object@metric <- metric
    callNextMethod(.Object)})



#' @name distances
#' @aliases Euclid_dist
#' @aliases gc_dist
#' @aliases gc_dist_time
#' @title Pre-configured distances
#'
#' @description Distance objects included in package
#'
#' @param dim dimension of Eucledian space
#' @param R great-circle radius
#' @details Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points.  Currently the Euclidean distance and the great-circle distance are included. Distances are computed using functions from the package \code{fields}.
#' @export
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  fields::rdist(x1,x2), dim=dim)
}

#' @rdname distances
#' @export
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  fields::rdist.earth(x1,x2,miles=F,R=R),dim=2L)
}

#' @rdname distances
#' @export
gc_dist_time <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  {
        spatdist <- fields::rdist.earth(x1[,1:2,drop=FALSE],x2[,1:2,drop=FALSE],miles=F,R=R)
        tdist <- fields::rdist(x1[,3],x2[,3])
        sqrt(spatdist^2 + tdist^2) } ,dim=3L)
}

#' @title Convert data frame to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object
#' @param df data frame containing polygon information, see details
#' @param keys vector of variable names used to group rows belonging to the same polygon
#' @param  coords vector of variable names identifying the coordinate columns
#' @param proj the projection of the \code{SpatialPolygons} object. Needs to be of class \code{CRS}
#' @details Each row in the data frame \code{df} contains both coordinates and labels (or keys) that identify to which polygon the coordinates belong. This function groups the data frame according to \code{keys} and forms a \code{SpatialPolygons} object from the coordinates in each group. It is important that all rings are closed, that is, that the last row of each group is identical to the first row. Since the keys can be arbitrary, and can contain more than one element, we identify each polygon with a new key by forming an MD5 hash made out of the respective \code{keys} variables that in themselves are unique (and therefore the hashed key is also, for most practicaly purposes, unique). For lon-lat coordinates use \code{proj = CRS("+proj=longlat")}.
#' @export
#' @examples
#' library(sp)
#' opts_FRK$set("parallel",0L)
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
    if(0) {
        df_poly <- rhwrapper(Ntot = nrow(df),
                             N = 4000,
                             f_expr = .rhdlply,
                             df=df,
                             keys=keys,
                             coords=coords,
                             dfun=parse(text = deparse(dfun)))
    } else {
        if(opts_FRK$get("parallel") > 1) {
            cl <- makeCluster(opts_FRK$get("parallel"))
            doParallel::registerDoParallel(opts_FRK$get("parallel"))
            df_poly <- plyr::dlply(df,keys,dfun,.parallel=TRUE)
            stopCluster(cl)
        } else {
            df_poly <- plyr::dlply(df,keys,dfun)
        }
    }
    Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=proj)
}

#' @title SpatialPolygonsDataFrame to df
#' @description Convert \code{SpatialPolygonsDataFrame} object to data frame
#' @param sp_polys object of class \code{SpatialPolygonsDataFrame}
#' @param vars variables to put into data frame (by default all of them)
#' @details This function is mainly used for plotting \code{SpatialPolygonsDataFrame} objects with \code{ggplot} rather than \code{spplot}. The coordinates of each polygon are extracted, and concatenated into one long data frame. The attributes of each polygon are then attached to this data frame as variables which vary by polygon \code{id}.  The returned \code{id} variable describes the polygon `id' and varies from 1 to the number of polygons represented in the data frame.
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
          function(data_sp,sp_pols,av_var)
          {
              if(is(data_sp,"SpatialPointsDataFrame")) {
                  Nobs <- NULL
                  data_sp$Nobs <- 1

                  if(!(opts_FRK$get("Rhipe"))) {
                      timer <- system.time(Data_in_BAU <- over(sp_pols,data_sp[c(av_var,"Nobs","std")],fn=sum))
                  } else {
                      print("Using RHIPE to find overlays")
                      timer <- system.time(
                          Data_in_BAU <- rhwrapper(Ntot = length(sp_pols),
                                                   N = 4000,
                                                   f_expr = .rhover,
                                                   sp_pols = sp_pols,
                                                   data_sp = data_sp,
                                                   av_var=av_var)
                      )
                  }
                  print(paste0("Binned data in ",timer[3]," seconds"))

                  sp_pols@data[av_var] <- Data_in_BAU[av_var]/Data_in_BAU$Nobs
                  sp_pols@data["std"] <- Data_in_BAU["std"]/Data_in_BAU$Nobs
                  sp_pols@data["Nobs"] <- Data_in_BAU$Nobs

                  new_sp_pts <- SpatialPointsDataFrame(coords=sp_pols[coordnames(data_sp)]@data,
                                                       data=sp_pols@data,
                                                       proj4string = CRS(proj4string(data_sp)))
                  #new_sp_pts$std <- sqrt(1 / new_sp_pts$Nobs)
                  new_sp_pts$std <- sqrt(new_sp_pts$std^2 / new_sp_pts$Nobs)
                  new_sp_pts <- subset(new_sp_pts,!is.na(Nobs))

              } else {
                  BAUs_aux_data <- over(data_sp,SpatialPointsDataFrame(coordinates(sp_pols),sp_pols@data))
                  stopifnot(all(row.names(BAUs_aux_data) == row.names(data_sp)))
                  data_sp@data <- cbind(data_sp@data,BAUs_aux_data)
                  data_sp$Nobs <- 1
                  new_sp_pts <- data_sp
              }

              new_sp_pts
          })

setMethod("map_data_to_BAUs",signature(data_sp="ST"),
          function(data_sp,sp_pols,av_var) {

              sp_fields <- NULL

              ## Bin every spatial frame separately
              sp_fields <- lapply(seq_along(sp_pols@time),
                                  function(i) {
                                      t1 <- time(sp_pols)[i]
                                      if(i < last(sp_pols@time)) {
                                          t2 <- time(sp_pols)[i+1]
                                          trange <- paste0(format(t1),"::",format(t2-1))
                                      } else {
                                          trange <- paste0(format(t1),"::",format(last(data_sp@endTime)))
                                      }
                                      data_spatial <- as(data_sp[,trange],"Spatial")
                                      BAU_spatial <- sp_pols[,i]
                                      BAU_spatial@data <- filter(BAU_spatial@data,time == t1)
                                      BAU_spatial@data <- cbind(BAU_spatial@data,coordinates(BAU_spatial))
                                      map_data_to_BAUs(data_spatial,
                                                       BAU_spatial,
                                                       av_var=av_var)})


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
    if(!("Nobs" %in% names(sp_pts))) stop("Nobs (number of observations in grid cell) needs to be a field of the Spatial object")

    g <- gstat::gstat(formula=variogram.formula,data=sp_pts)
    v <- gstat::variogram(g,cressie=T)
    vgm.fit = gstat::fit.variogram(v, model = gstat::vgm(1, "Lin", mean(v$dist), 1))
    plot(v,vgm.fit)
    print(paste0("sigma2e estimate = ",vgm.fit$psill[1]))
    if(vgm.fit$psill[1] == 0)
        stop("Observational error estimated to be zero. Please consider using finer BAUs or do not attempt to estimate observation error")
    sp_pts$std <- sqrt(vgm.fit$psill[1] / sp_pts$Nobs)

    warning("Observational error estimation could be improved. Currently a variogram is fitted to the data, and then the error variance of a single observation is assumed to be the partial sill. Then the variance of the averaged observations in the BAU is divided by Nobs. Currently there is no accounting for multiple data in the same grid box during variogram fitting as it's not straightforward with gstat.")

    sp_pts

}


setMethod("BuildC",signature(data="SpatialPolygonsDataFrame"),
          function(data,BAUs) {
              data$id <- 1:length(data)
              overlap <- over(SpatialPoints(coordinates(BAUs)),data)
              i_idx <- as.numeric(na.exclude(overlap$id))
              j_idx <- which(!is.na(overlap$id))
              list(i_idx=i_idx,j_idx=j_idx)
          })

setMethod("BuildC",signature(data="SpatialPointsDataFrame"),
          function(data,BAUs) {
              i_idx <- 1:length(data)
              j_idx <- which(row.names(BAUs) %in% row.names(data))
              list(i_idx=i_idx,j_idx=j_idx)
          })

setMethod("BuildC",signature(data="STIDF"),
          function(data,BAUs) {
              i_idx <- 1:length(data)
              j_idx <- which(BAUs@data$n %in% data@data$n)
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
setMethod("dimensions",signature("manifold"),function(obj){dimensions(obj@metric)})

#' @rdname distance
#' @aliases distance,measure-method
setMethod("distance",signature("measure"),function(d,x1,x2){d@dist(x1,x2)})

#' @rdname distance
#' @aliases distance,manifold-method
setMethod("distance",signature("manifold"),function(d,x1,x2){distance(d@metric,x1,x2)})


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

setMethod("manifold",signature(.Object="TensorP_Basis"),function(.Object) {
    return(NA)
})


#' @aliases coordnames,STFDF-method
setMethod("coordnames",signature(x="STFDF"),function(x) {
    return(coordnames(x@sp))
})

#' @aliases coordnames,STIDF-method
setMethod("coordnames",signature(x="STIDF"),function(x) {
    return(c(coordnames(x@sp),"t"))
})


