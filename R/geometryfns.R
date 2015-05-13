#' @title manifold
#'
#' @description n This function initialises a manifold
setMethod("initialize",signature="manifold",function(.Object) {
    ## General manifold checks can come in here
    .Object
})


#' @title sphere
#'
#' @description This function initialises a sphere
#' @export
sphere <- function(radius=1,metric=gc_dist(R=radius)) {
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
#' @description This function initialises a plane
#' @export
plane <- function(metric=Euclid_dist(dim=2L)) {
    stopifnot(dimensions(metric)==2L)
    new("plane",metric=metric)
}

setMethod("initialize",signature="plane",function(.Object,metric=Euclid_dist(dim=2L)) {
    .Object@type <- "plane"
    .Object@metric <- metric
    callNextMethod(.Object)})



#' @title real_line
#'
#' @description This function initialises a sphere
#' @export
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
#' @description This function initialises a great-circle distance measure
gc_dist <- function(R=NULL) {
    new("measure",dist=function(x1,x2)  fields::rdist.earth(x1,x2,miles=F,R=R),dim=2L)
}

#' @title Euclidean distance
#'
#' @description This function initialises a Euclidean distance measure on 2D
Euclid_dist <- function(dim=2L) {
    stopifnot(is.integer(dim))
    new("measure",dist=function(x1,x2)  fields::rdist(x1,x2), dim=dim)
}

#' @title domain
#'
#' @description This function initialises a boundary polygon
domain <- function(m=sphere(),bndary=matrix(c(0,1,0,1),2,2)) {
    stopifnot(nrow(bndary) >= 2)
    stopifnot(ncol(bndary) == dimensions(m))
    new("domain",m,bndary)
}

setMethod("initialize",signature="domain",function(.Object,m = sphere(), bndary=matrix(c(0,1,0,1),2,2)) {
    .Object@m <- m
    .Object@bndary <- bndary
    .Object
})



#' @title df to SpatialPolygons
#' @description Convert data frame to SpatialPolygons object
#' @export
df_to_SpatialPolygons <- function(df,keys,coords,proj) {
    df_poly <- plyr::dlply(df,keys,function(d) {
        Sr <- Polygons(list(Polygon(d[coords])),digest::digest(d[keys]))})
    Sr <- SpatialPolygons(df_poly,1:length(df_poly),proj4string=proj)
}

#' @title SpatialPolygons to df
#' @description Convert SpatialPolygons object to data frame
#' @export
SpatialPolygons_to_df <- function(sp_polys,vars) {
    sp_polys$id <- 1 : length(sp_polys)
    polynames <- 1 : length(sp_polys)
    X <- data.frame(do.call("rbind",lapply(1:length(sp_polys),
                                           function(i)
                                               cbind(sp_polys@polygons[[i]]@Polygons[[1]]@coords,
                                                     id=polynames[i])))) %>%
    left_join(sp_polys@data[c("id",vars)])
}


map_data_to_BAUs <- function(data_sp,sp_pols,av_var,variogram.formula=NULL) {
    if(is(data_sp,"SpatialPointsDataFrame")) {
        data_sp$Nobs <- 1
        Data_in_BAU <- over(sp_pols,data_sp[c(av_var,"Nobs")],fn=sum)
        sp_pols@data[av_var] <- Data_in_BAU[av_var]/Data_in_BAU$Nobs
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

    if(is(variogram.formula,"formula")) {
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

#' @aliases dimensions,domain-method
setMethod("dimensions",signature("domain"),function(obj){dimensions(obj@m)})

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


