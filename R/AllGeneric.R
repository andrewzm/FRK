#### GENERIC FUNCTIONS ######
#' @title Show basis functions
#' @description Generic plotting function for visualising the basis functions.
#' @param basis object of class \code{Basis}
#' @param g object of class \code{gg} (a \code{ggplot} object) over which to overlay the basis functions (optional)
#' @param ... not in use
#' @details The function \code{show_basis} adapts its behaviour to the manifold being used. With \code{real_line}, the 1D basis functions are plotted with colour distinguishing between the different resolutions. With \code{plane}, only local basis functions are supported (at present). Each basis function is shown as a circle with diameter equal to the \code{scale} parameter of the function. Linetype distinguishes the resolution. With \code{sphere}, the centres of the basis functions are shown as circles, with larger sizes corresponding to lower (i.e., coarser) resolutions. Space-time basis functions of subclass \code{TensorP_Basis} are be visualised by showing the spatial basis functions and the temporal basis functions in two separate plots.
#' @examples
#' library(ggplot2)
#' library(sp)
#' data(meuse)
#' coordinates(meuse) = ~x+y # change into an sp object
#' G <- auto_basis(manifold = plane(),data=meuse,nres = 2,regular=2,prune=0.1,type = "bisquare")
#' \dontrun{show_basis(G,ggplot()) + geom_point(data=data.frame(meuse),aes(x,y))}
#' @export
setGeneric("show_basis", function(basis,...) standardGeneric("show_basis"))

#' @title Retrieve manifold
#' @description Retrieve manifold from \code{FRK} object.
#' @param .Object \code{FRK} object
#' @examples
#' G <-  local_basis(manifold = plane(),
#'                    loc=matrix(0,1,2),
#'                    scale=0.2,
#'                    type="bisquare")
#' manifold(G)
#' @export
setGeneric("manifold", function(.Object) standardGeneric("manifold"))

#' @title Number of basis functions
#' @description Retrieve the number of basis functions from \code{Basis} or \code{SRE} object.
#' @param .Object object of class \code{Basis} or \code{SRE}
#' @examples
#' library(sp)
#' data(meuse)
#' coordinates(meuse) = ~x+y # change into an sp object
#' G <- auto_basis(manifold = plane(),
#'                 data=meuse,
#'                 nres = 2,
#'                 regular=1,
#'                 type = "Gaussian")
#' print(nbasis(G))
#' @export
setGeneric("nbasis", function(.Object) standardGeneric("nbasis"))

#' @title Type of manifold
#' @description Retrieve slot \code{type} from object
#' @param .Object object of class \code{Basis} or \code{manifold}
#' @examples
#' S <- sphere()
#' print(type(S))
#' @export
setGeneric("type", function(.Object) standardGeneric("type"))

#' @title Compute distance
#' @description Compute distance using object of class \code{measure} or \code{manifold}.
#' @param d object of class \code{measure} or \code{manifold}
#' @param x1 first coordinate
#' @param x2 second coordinate
#' @examples
#' distance(sphere(),matrix(0,1,2),matrix(10,1,2))
#' distance(plane(),matrix(0,1,2),matrix(10,1,2))
#' @export
setGeneric("distance", function(d,x1,x2=NULL) standardGeneric("distance"))

#' @title Evaluate basis functions
#' @description  Evaluate basis functions at points or average functions over polygons.
#' @param basis object of class \code{Basis}
#' @param s object of class \code{matrix}, \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame} containing the spatial locations/footprints
#' @details This function evaluates the basis functions at isolated points, or averages
#' the basis functions over polygons, for computing the matrix \eqn{S}. The latter
#' operation is carried out using Monte Carlo integration with 1000 samples per polygon. When
#' using space-time basis functions, the object MUST contain a field \code{t} containing a numeric
#' representation of the time, for example, containing the number of second, hours, or days since the first
#' data point.
#' @examples
#' library(sp)
#'
#' ### Create a synthetic dataset
#' set.seed(1)
#' d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
#'                 lat = runif(n=1000,min = -90, max = 90),
#'                 z = rnorm(5000))
#' coordinates(d) <- ~lon + lat
#' proj4string(d)=CRS("+proj=longlat")
#'
#' ### Now create basis functions on sphere
#' G <- auto_basis(manifold = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#'
#' ### Now evaluate basis functions at origin
#' S <- eval_basis(G,matrix(c(0,0),1,2))
#' @export
setGeneric("eval_basis", function(basis,s) standardGeneric("eval_basis"))

#' @title Tensor product of basis functions
#' @description Constructs a new set of basis by finding the tensor product of two sets of basis functions.
#' @param Basis1 first set of basis functions
#' @param Basis2 second set of basis functions
#' @export
#' @examples
#' library(spacetime)
#' library(sp)
#' library(dplyr)
#' sim_data <- data.frame(lon = runif(20,-180,180),
#'                        lat = runif(20,-90,90),
#'                        t = 1:20,
#'                        z = rnorm(20),
#'                        std = 0.1)
#' time <- as.POSIXct("2003-05-01",tz="") + 3600*24*(sim_data$t-1)
#' space <- sim_data[,c("lon","lat")]
#' coordinates(space) = ~lon+lat # change into an sp object
#' proj4string(space)=CRS("+proj=longlat +ellps=sphere")
#' STobj <- STIDF(space,time,data=sim_data)
#' G_spatial <- auto_basis(manifold = sphere(),
#'                         data=as(STobj,"Spatial"),
#'                         nres = 1,
#'                         type = "bisquare",
#'                         subsamp = 20000)
#' G_temporal <- local_basis(manifold=real_line(),loc = matrix(c(1,3)),scale = rep(1,2))
#' G <- TensorP(G_spatial,G_temporal)
#' # show_basis(G_spatial)
#' # show_basis(G_temporal)
setGeneric("TensorP", function(Basis1,Basis2) standardGeneric("TensorP"))

#' @title Return  the number of resolutions
#' @description Return the number of resolutions from a basis function object.
#' @param b object of class \code{Basis} or \code{SRE}
#' @examples
#' library(sp)
#' set.seed(1)
#' d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
#'                 lat = runif(n=1000,min = -90, max = 90),
#'                 z = rnorm(5000))
#' coordinates(d) <- ~lon + lat
#' proj4string(d)=CRS("+proj=longlat")
#'
#' ### Now create basis functions on sphere
#' G <- auto_basis(manifold = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#' nres(G)
#' @export
setGeneric("nres", function(b) standardGeneric("nres"))

#' @title Basis-function data frame object
#' @description Tools for retrieving and manipulating the data frame within the Basis objects. Use the assignment \code{data.frame()<-} with care; no checks are made to make sure the data frame conforms with the object. Only use if you know what you're doing
#' @param x the obect of class \code{Basis} we are assigning the new data to or retrieving data from
#' @param value the new data being assigned to the Basis object
#' @param name the field name to which values will be retrieved or assigned inside the Basis' data frame
#' @param ... unused
#' @rdname Basis_data.frame
#' @examples
#' G <- local_basis()
#' df <- data.frame(G)
#' print(df$res)
#' df$res <- 2
#' data.frame(G) <- df
#' @export
Basis_as_data.frame <- setAs("Basis", "data.frame",
                             function(from) as.data.frame.Basis(from))

#' @rdname Basis_data.frame
#' @export
TensorP_Basis_as_data.frame <- setAs("TensorP_Basis", "data.frame",
                                     function(from) as.data.frame.TensorP_Basis(from))

#' @rdname Basis_data.frame
setGeneric("data.frame<-", function(x, value) standardGeneric("data.frame<-"))


#' @title Creates pixels around points
#' @description Takes a SpatialPointsDataFrame and converts it into SpatialPolygonsDataFrame by constructing a tiny (within machine tolerance) BAU around each SpatialPoint.
#' @param obj object of class \code{SpatialPointsDataFrame}
#' @details This function is there to allow users to mimic standard geospatial analysis where BAUs are not used. Since \code{FRK} is build on the concept of a BAU, this function constructs tiny BAUs around the observation and prediction locations which can be subsequently passed on to the functions \code{SRE} and \code{FRK}. With \code{BAUs_from_points}, the user supplies both the data and prediction locations accompanied with covariates.
#' @export
#' @examples
#' library(sp)
#' opts_FRK$set("parallel",0L)
#' df <- data.frame(x = rnorm(10),
#'                  y = rnorm(10))
#' coordinates(df) <- ~x+y
#' BAUs <- BAUs_from_points(df)
setGeneric("BAUs_from_points", function(obj) standardGeneric("BAUs_from_points"))

########################
#### NOT EXPORTED ######
########################

#' @title Automatic BAU generation
#' @noRd
#' @description This generic function is called by \code{auto_BAUs} after a series of checks.
#' @param manifold object of class \code{manifold}
#' @param type either ``hex'' or ``grid'', indicating whether gridded or hexagonal BAUs should be used
#' @param cellsize denotes the length of the sides of the gridcell when \code{type} == ``grid''
#' @param resl resolution number of isea3h DGGRID cells for when \code{type} is ``hex'' and \code{manifold} is \code{sphere}
#' @param d data, that is, an object of class SpatialPointsDataFrame or SpatialPolygonsDataFrame. Provision of data implies that the domain is bounded (necessary with \code{real_line} and \code{plane} but not necessary with \code{sphere})
#' @param nonconvex_hull flag indicating whether \code{INLA} should be used to create a non-convex domain boundary
#' @param convex convex parameter for the \code{INLA} function \code{inla.nonconvex.hull} used for smoothing an extended boundary when working on a finite domain (that is, when the object \code{d} is supplied), see details
#' @param xlims limits of the horizontal axis (overrides automatic selection).
#' @param ylims limits of the vertical axis (overrides automatic selection).
#' @param ... currently unused
#' @details This generic function is not called directly. Please refer to \code{auto_BAUs} for more details.
setGeneric("auto_BAU", function(manifold,type,cellsize,resl,d,nonconvex_hull,convex,xlims,ylims,...) standardGeneric("auto_BAU"))

#' @title Bin data into BAUs
#' @description This is an internal function which bins data into BAUs or aggregates across BAUs if the data have a large footprint. If \code{est_error == TRUE}, the observation error is estimated using the variogram (see vignette for details).
#' @param data_sp object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param sp_pols object of class \code{SpatialPolygonsDataFrame} that contains the BAUs
#' @param variogram.formula formula used for detrending the data for variogram estimation of the observation error. Should be identical to that used for \code{SRE()}
#' @param est_error flag indicating whether variogram estimation of the observation error should be carried out or no. This can take a long time with large datasets
#' @param average_in_BAU flag indicating whether to summarise data that fall into a single BAU by simply taking an average of the data and the standard devitation of the data within each BAU (suitable for extremely large datasets)
#' @details This generic function is not called directly. It is called in the SRE function for binning data in BAUs
#' @noRd
setGeneric("map_data_to_BAUs", function(data_sp,sp_pols,variogram.formula=NULL,est_error=T,average_in_BAU = TRUE) standardGeneric("map_data_to_BAUs"))

#' @title Concatenation
#' @description Concatenates FRK objects of the same class together. This is primarily
#' used to join up \code{Basis} blocks together.
#' @param ... a series of \code{FRK} objects
#' @noRd
setGeneric("concat", function(...) standardGeneric("concat"))



#' @title Counts of basis per resolution
#' @description Returns a data frame with two columns, containing the resolution number and the number of basis functions at that resolution.
#' @param .Object an \code{FRK} object of class \code{Basis} or class \code{SRE}
#' @noRd
setGeneric("count_res", function(.Object) standardGeneric("count_res"))


#' @title Construct incidence matrix
#' @description Construct incidence matrix for Spatial and Spatiotemporal fields by mapping the data to the BAUs.
#' @param data object of class \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param BAUs object of class \code{SpatialPolygonsDataFrame} or \code{SpatialPixelsDataFrame}
#' @noRd
setGeneric("BuildC", function(data,BAUs) standardGeneric("BuildC"))

#' @title Construct distance matrices
#' @description Construct distance matrices for (local) basis functions by computing the distances between the respective centroids by resolution.
#' @param G an object of class \code{Basis}
#' @noRd
setGeneric("BuildD", function(G) standardGeneric("BuildD"))


