#' @title manifold
#' @description Buffer manifold initialisation. This function should not be called directly.
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here
    .Object
})

#' @title Convert dggrid gen file to data frame
#'
#' @description Convert dggrid gen file to data frame
#' @export
dggrid_gen_to_df <- function(filename) {
    if(!is.character(filename)) stop("filename needs to of type 'character'")
    if(!require(zoo)) stop("zoo is needed for this function. Please install it separately")
    read.table(filename,
               sep=" ",fill=T,
               header=F,
               col.names = c("id","lon","lat")) %>%
        filter(!(id == "END")) %>%
        mutate(res = j,
               id = as.numeric(as.character(id)),
               centroid = as.numeric(!is.na(id))) %>%
        transform(id = na.locf(id))
}


#' @title Auatomatic BAU generation
#'
#' @description Autoatic BAU generation
#' @export
auto_BAUs <- function(manifold,res,data=NULL) {
    if(!(is(data,"Spatial") | is.null(data))) stop("Data needs to be of class 'Spatial' or NULL")
    if(!is(manifold,"manifold")) stop("manifold needs to be of class 'manifold'")
    if(!is.numeric(res) | is.integer(res)) stop("res needs to be of type 'numeric' or 'integer'")
    if(!(res >=0 & res <= 9)) stop("res needs to be between 0 and 9")
    resl <- round(res)

    if(is(manifold,"sphere")) {
        data(isea3h)
        isea3h_res <- filter(isea3h,res == resl) %>%
            arrange(id) %>%
            group_by(id) %>%
            filter(diff(range(lon)) < 90) %>%
            data.frame() %>%
            mutate(fs=1)

        isea3h_sp_pol <- df_to_SpatialPolygons(
            df=filter(isea3h_res,centroid==0),
            keys=c("id"),
            coords=c("lon","lat"),
            proj=CRS("+proj=longlat"))

        isea3h_sp_poldf <- SpatialPolygonsDataFrame(
            isea3h_sp_pol,
            cbind(data.frame(row.names=names(isea3h_sp_pol)),
                  (filter(isea3h_res,centroid==1) %>%
                       select(id,lon,lat,fs))))
        if(!is.null(data)) {
            sub_pols <- coordinates(data) %>%
                chull()
            conv_hull <- coordinates(data)[c(sub_pols,sub_pols[1]),] %>%
                data.frame() %>%
                mutate(id=1) %>%
                df_to_SpatialPolygons(keys="id",
                                      coords=c("lon","lat"),
                                      proj=CRS("+proj=longlat"))
            isea3h_sp_poldf$in_chull <- over(isea3h_sp_poldf,conv_hull)
            isea3h_sp_poldf <- subset(isea3h_sp_poldf,!is.na(in_chull))
        }
        isea3h_sp_poldf
    }

}

#' @title sphere
#'
#' @description Sphere initialisation.
#'
#' @param radius radius of sphere
#'
#' @details A sphere is initialised using a \code{radius} parameter. By default, the radius \code{R} is equal to \code{R}=6371 km, the Earth's radius, while the measure used to compute distances on the sphere is the great-circle distance on a sphere of radius \code{R}. Currently the user does not have the option to alter the distance function used on the sphere. Distances are computed using functions from the \code{fields} package.
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

#' @title plane
#'
#' @description Plane initialisation.
#'
#' @param metric an object of class \code{measure}
#'
#' @details A 2D plane is initialised using a \code{measure} object. By default, the measure object (\code{metric}) is the Euclidean distance in 2 dimensions. Distances are computed using functions from the \code{fields} package.
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



#' @title real line
#'
#' @description Real line initialisation.
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



#' @title Great-circle distance
#'
#' @description Great-circle distance object
#'
#' @param R radius of sphere
#' @description Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points. Distances are computed using functions from the package \code{fields}
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  fields::rdist.earth(x1,x2,miles=F,R=R),dim=2L)
}

#' @title Euclidean distance
#'
#' @description Euclidean distance object
#'
#' @param dim dimension of Eucledian space
#' @description Initialises an object of class \code{measure} which contains a function \code{dist} used for computing the distance between two points in Euclidean space.  Distances are computed using functions from the package \code{fields}.
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  fields::rdist(x1,x2), dim=dim)
}


#' @title df to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object
#' @param df data frame containing polygon information. See details
#' @param keys vector of variable names used to group rows pertaining to the same polygon
#' @param  coords vector of variable names identifying the coordinate columns
#' @param proj the projection of the SpatialPolygons object. Needs to be of class \code{CRS}
#' @details Each row in the data frame \code{df} contains both coordinates and labels (or keys) that identify to which polygon the coordinates belong. This function groups the data frame according to \code{keys} and forms a \code{SpatialPolygons} object from the coordinates in each group. It is important that all rings are closed, that is, that the last row of each group equals the first row. Since the keys can be arbitrary, we identify each polygon with a new key by forming an MD5 has of the respective \code{keys} variables. For lon-lat coordinates use \code{proj = CRS("+proj=longlat")}.
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
    if(0) {
        df_poly <- rhwrapper(Ntot = nrow(df),
                             N = 4000,
                             f_expr = .rhdlply,
                             df=df,
                             keys=keys,
                             coords=coords,
                             dfun=parse(text = deparse(dfun)))
    } else {
        if(defaults$parallel > 0) {
            cl <- makeCluster(defaults$parallel)
            doParallel::registerDoParallel(defaults$parallel)
            df_poly <- plyr::dlply(df,keys,dfun,.parallel=TRUE)
            stopCluster(cl)
        } else {
            df_poly <- plyr::dlply(df,keys,dfun)
        }
    }
    Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=proj)
}

#' @title SpatialPolygonsDataFrame to df
#' @description Convert SpatialPolygonsDataFrame object to data frame
#' @param sp_polys object of class SpatialPolygonsDataFrame
#' @param vars variables to put into data frame (by default all of them)
#' @details This function is mainly used for plotting SpatialPolygonsDataFrame objects with ggplot. The coordinates of each polygon are extracted, and concatenated into one long data frame. The attributes of each polygon are then attached to this data frame as variables which vary by polygon \code{id}.  The returned \code{id} variable describes the polygon 'id' and varies from 1 to the number of polygons represented in the data frame.
#' @export
#' @examples
#' library(sp)
#' library(ggplot2)
#' df <- data.frame(id = c(rep(1,4),rep(2,4)),
#'                  x = c(0,1,0,0,2,3,2,2),
#'                  y=c(0,0,1,0,0,1,1,0))
#' pols <- df_to_SpatialPolygons(df,"id",c("x","y"),CRS())
#' polsdf <- SpatialPolygonsDataFrame(pols,data.frame(p = c(1,2),row.names=row.names(pols)))
#' df2 <- SpatialPolygonsDataFrame_to_df(polsdf)
#' ggplot(df2,aes(x=x,y=y,group=id)) + geom_polygon()
SpatialPolygonsDataFrame_to_df <- function(sp_polys,vars = names(sp_polys)) {
    sp_polys$id <- 1 : length(sp_polys)
    polynames <- 1 : length(sp_polys)
    X <- data.frame(do.call("rbind",lapply(1:length(sp_polys),
                                           function(i)
                                               cbind(sp_polys@polygons[[i]]@Polygons[[1]]@coords,
                                                     id=polynames[i])))) %>%
    left_join(sp_polys@data[c("id",vars)])
    X
}

#' @title Bin data into BAUs
#' @description This is an internal function which bins data into BAUs or aggregates across BAUs if the data have a large footprint. If est_error == T, the observation error is estimated as in Katzfuss & Cressie (2011)
map_data_to_BAUs <- function(data_sp,sp_pols,av_var,variogram.formula=NULL,est_error=T) {
    if(is(data_sp,"SpatialPointsDataFrame")) {
        Nobs <- NULL
        data_sp$Nobs <- 1
        if(est_error) data_sp$std <- 0 ## Just set it to something, this will be overwritten later on

        if(!defaults$Rhipe) {
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
        new_sp_pts$std <- sqrt(1 / new_sp_pts$Nobs)
        new_sp_pts <- subset(new_sp_pts,!is.na(Nobs))

    } else {
        BAUs_aux_data <- over(data_sp,SpatialPointsDataFrame(coordinates(sp_pols),sp_pols@data))
        stopifnot(all(row.names(BAUs_aux_data) == row.names(data_sp)))
        data_sp@data <- cbind(data_sp@data,BAUs_aux_data)
        data_sp$Nobs <- 1
        new_sp_pts <- data_sp
    }

    if(is(variogram.formula,"formula") & (est_error==TRUE)) {
        g <- gstat::gstat(formula=variogram.formula,data=new_sp_pts)
        v <- gstat::variogram(g,cressie=T)
        warning("Not accounting for multiple data in the same grid box during variogram estimation. Need to see how to do this with gstat")
        vgm.fit = gstat::fit.variogram(v, model = gstat::vgm(1, "Lin", mean(v$dist), 1))
        plot(v,vgm.fit)
        print(paste0("sigma2e estimate = ",vgm.fit$psill[1]))
        if(vgm.fit$psill[1] == 0) stop("Observational error estimated to be zero. Please consider using finer BAUs")
        new_sp_pts$std <- sqrt(vgm.fit$psill[1] / new_sp_pts$Nobs)
    }
    new_sp_pts
}

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


