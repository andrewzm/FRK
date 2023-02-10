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

#### GENERIC FUNCTIONS ######
#' @title Show basis functions
#' @description Generic plotting function for visualising the basis functions.
#' @param basis object of class \code{Basis}
#' @param g object of class \code{gg} (a \code{ggplot} object) over which to overlay the basis functions (optional)
#' @param ... not in use
#' @details The function \code{show_basis} adapts its behaviour to the manifold being used. With \code{real_line}, the 1D basis functions are plotted with colour distinguishing between the different resolutions. With \code{plane}, only local basis functions are supported (at present). Each basis function is shown as a circle with diameter equal to the \code{scale} parameter of the function. Linetype distinguishes the resolution. With \code{sphere}, the centres of the basis functions are shown as circles, with larger sizes corresponding to coarser resolutions. Space-time basis functions of subclass \code{TensorP_Basis} are visualised by showing the spatial basis functions and the temporal basis functions in two separate plots.
#' @seealso  \code{\link{auto_basis}} for automatically constructing basis functions.
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
#' @seealso \code{\link{real_line}}, \code{\link{plane}}, \code{\link{sphere}}, \code{\link{STplane}} and \code{\link{STsphere}} for constructing manifolds.
#' @export
setGeneric("manifold", function(.Object) standardGeneric("manifold"))

#' @title Number of basis functions
#' @description Retrieve the number of basis functions from \code{Basis} or \code{SRE} object.
#' @param .Object object of class \code{Basis} or \code{SRE}
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions.
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
#' @seealso \code{\link{real_line}}, \code{\link{plane}}, \code{\link{sphere}}, \code{\link{STplane}} and \code{\link{STsphere}} for constructing manifolds.
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
#' @seealso \code{\link{real_line}}, \code{\link{plane}}, \code{\link{sphere}}, \code{\link{STplane}} and \code{\link{STsphere}} for constructing manifolds, and \code{\link{distances}} for the type of distances available.
#' @export
setGeneric("distance", function(d,x1,x2=NULL) standardGeneric("distance"))

#' @title Evaluate basis functions
#' @description  Evaluate basis functions at points or average functions over polygons.
#' @param basis object of class \code{Basis}
#' @param s object of class \code{matrix}, \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame} containing the spatial locations/footprints
#' @details This function evaluates the basis functions at isolated points, or averages
#' the basis functions over polygons, for computing the matrix \eqn{S}. The latter
#' operation is carried out using Monte Carlo integration with 1000 samples per polygon. When
#' using space-time basis functions, the object must contain a field \code{t} containing a numeric
#' representation of the time, for example, containing the number of seconds, hours, or days since the first
#' data point.
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions.
#' @examples
#' library(sp)
#'
#' ### Create a synthetic dataset
#' set.seed(1)
#' d <- data.frame(lon = runif(n=500,min = -179, max = 179),
#'                 lat = runif(n=500,min = -90, max = 90),
#'                 z = rnorm(500))
#' coordinates(d) <- ~lon + lat
#' slot(d, "proj4string") = CRS("+proj=longlat")
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
#' @description Constructs a new set of basis functions by finding the tensor product of two sets of basis functions.
#' @param Basis1 first set of basis functions
#' @param Basis2 second set of basis functions
#' @export
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions.
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
#' slot(space, "proj4string") = CRS("+proj=longlat +ellps=sphere")
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
#' d <- data.frame(lon = runif(n=500,min = -179, max = 179),
#'                 lat = runif(n=500,min = -90, max = 90),
#'                 z = rnorm(500))
#' coordinates(d) <- ~lon + lat
#' slot(d, "proj4string") = CRS("+proj=longlat")
#'
#' ### Now create basis functions on sphere
#' G <- auto_basis(manifold = sphere(),data=d,
#'                 nres = 2,prune=15,
#'                 type = "bisquare",
#'                 subsamp = 20000)
#' nres(G)
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions.
#' @export
setGeneric("nres", function(b) standardGeneric("nres"))


#' @rdname Basis_data.frame
setGeneric("data.frame<-", function(x, value) standardGeneric("data.frame<-"))

#' @title Creates pixels around points
#' @description Takes a SpatialPointsDataFrame and converts it into SpatialPolygonsDataFrame by constructing a tiny (within machine tolerance) BAU around each SpatialPoint.
#' @param obj object of class \code{SpatialPointsDataFrame}
#' @param offset edge size of the mini-BAU (default 1e-10)
#' @details This function allows users to mimic standard geospatial analysis where BAUs are not used. Since \code{FRK} is built on the concept of a BAU, this function constructs tiny BAUs around the observation and prediction locations that can be subsequently passed on to the functions \code{SRE} and \code{FRK}. With \code{BAUs_from_points}, the user supplies both the data and prediction locations accompanied with covariates.
#' @seealso \code{\link{auto_BAUs}} for automatically constructing generic BAUs.
#' @export
#' @examples
#' library(sp)
#' opts_FRK$set("parallel",0L)
#' df <- data.frame(x = rnorm(10),
#'                  y = rnorm(10))
#' coordinates(df) <- ~x+y
#' BAUs <- BAUs_from_points(df)
setGeneric("BAUs_from_points", function(obj,offset = 1e-10)
    standardGeneric("BAUs_from_points"))

#' @title Removes basis functions
#' @description Takes an object of class \code{Basis} and returns an object of class \code{Basis} with selected basis functions removed
#' @param Basis object of class \code{Basis}
#' @param rmidx indices of basis functions to remove. Or a \code{SpatialPolygons} object; basis functions overlapping this \code{SpatialPolygons} object will be \emph{retained}
#' @export
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions
#' @examples
#' library(sp)
#' df <- data.frame(x = rnorm(10),
#'                  y = rnorm(10))
#' coordinates(df) <- ~x+y
#' G <- auto_basis(plane(),df,nres=1)
#' data.frame(G) # Print info on basis
#' 
#' ## Removing basis functions by index
#' G_subset <- remove_basis(G, 1:(nbasis(G)-1))
#' data.frame(G_subset)
#' 
#' ## Removing basis functions using SpatialPolygons
#' x <- 1
#' poly <- Polygon(rbind(c(-x, -x), c(-x, x), c(x, x), c(x, -x), c(-x, -x)))
#' polys <- Polygons(list(poly), "1")
#' spatpolys <- SpatialPolygons(list(polys))
#' G_subset <- remove_basis(G, spatpolys)
#' data.frame(G_subset)
setGeneric("remove_basis", function(Basis,rmidx)
    standardGeneric("remove_basis"))


#' @title Combine basis functions
#' @description Takes a list of objects of class \code{Basis} and returns a 
#' single object of class \code{Basis}.   
#' @param Basis_list a list of objects of class \code{Basis}. Each element of the list is assumed to 
#' represent a single resolution of basis functions
#' @export
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions
#' @examples
#'## Construct two resolutions of basis functions using local_basis() 
#'Basis1 <- local_basis(manifold = real_line(), 
#'                      loc = matrix(seq(0, 1, length.out = 3), ncol = 1), 
#'                      scale = rep(0.4, 3))
#'
#'Basis2 <- local_basis(manifold = real_line(), 
#'                      loc = matrix(seq(0, 1, length.out = 6), ncol = 1), 
#'                      scale = rep(0.2, 6))
#'
#'## Combine basis-function resolutions into a single Basis object
#'combine_basis(list(Basis1, Basis2)) 
setGeneric("combine_basis", function(Basis_list)
    standardGeneric("combine_basis"))



#' @title Observed (or unobserved) BAUs
#' @description Computes the indices (a numeric vector) of the observed (or unobserved) BAUs
#' @param object object of class \code{SRE}
#' @seealso See \code{\link{FRK}} for more information on the SRE model and available fitting methods.
#' @examples 
#' # See example in the help file for FRK
#' @export
setGeneric("observed_BAUs", function(object)
    standardGeneric("observed_BAUs"))

#' @rdname observed_BAUs
#' @export
setGeneric("unobserved_BAUs", function(object)
    standardGeneric("unobserved_BAUs"))

#' @title Retrieve fit information for SRE model
#' @description Takes an object of class \code{SRE} and returns a list containing all the relevant information on parameter estimation
#' @param object object of class \code{SRE}
#' @seealso See \code{\link{FRK}} for more information on the SRE model and available fitting methods.
#' @examples
#' # See example in the help file for FRK
#' @export
setGeneric("info_fit", function(object)
    standardGeneric("info_fit"))

#' @title (Deprecated) Retrieve log-likelihood
#' @description This function is deprecated; please use \code{logLik}
#' @param object object of class \code{SRE}
#' @export
setGeneric("loglik", function(object)
    standardGeneric("loglik"))


#' @title Plot a Spatial*DataFrame or STFDF object
#' @description Takes an object of class \code{Spatial*DataFrame} or \code{STFDF}, and plots requested data columns using \code{ggplot2}
#' @param newdata an object of class \code{Spatial*DataFrame} or \code{STFDF}
#' @param column_names a vector of strings indicating the columns of the data to plot
#' @param map_layer (optional) a \code{ggplot} layer or object to add below the plotted layer, often a map
#' @param subset_time (optional) a vector of times to be included; applicable only for \code{STFDF} objects
#' @param palette the palette supplied to the argument \code{palette} of \code{scale_*_distiller()}. Alternatively, if \code{palette} = "nasa", a vibrant colour palette is created using \code{scale_*_gradientn()}
#' @param plot_over_world logical; if \code{TRUE}, \code{coord_map("mollweide")} and \code{\link{draw_world}} are used to plot over the world
#' @param labels_from_coordnames logical; if \code{TRUE}, the coordinate names of \code{newdata} (i.e., \code{coordnames(newdata)}) are used as the horizontal- and vertical-axis labels. Otherwise, generic names, s_1 and s_2, are used
#' @param ... optional arguments passed on to whatever geom is appropriate for the \code{Spatial*DataFrame} or \code{STFDF} object (\code{geom_point}, \code{geom_tile}, \code{geom_raster}, or \code{geom_polygon})
#' @return A list of \code{ggplot} objects corresponding to the provided \code{column_names}. This list can then be supplied to, for example, \code{ggpubr::ggarrange()}.
#' @seealso \code{\link{plot}}
#' @export
#' @examples 
#' ## See example in the help file for FRK
setGeneric("plot_spatial_or_ST", function(newdata, column_names,  map_layer=NULL, 
                                          subset_time=NULL, palette="Spectral", plot_over_world=FALSE, 
                                          labels_from_coordnames = TRUE, ...)
    standardGeneric("plot_spatial_or_ST"))


#' @title Uncertainty quantification of the fixed effects
#' @description Compute confidence intervals for the fixed effects (upper  
#' and lower bound specifed by percentiles; default 90\% confidence central interval)
#' @inheritParams SRE
#' @export
setGeneric("coef_uncertainty", function(object, percentiles = c(5, 95),  nsim = 400)
  standardGeneric("coef_uncertainty"))


## INHERITED GENERICS

#' @export
setGeneric("logLik", function(object)
  standardGeneric("logLik"))

#' @export
setGeneric("coef", function(object)
    standardGeneric("coef"))

#' @export
setGeneric("predict", function(object, ...)
    standardGeneric("predict"))

#' @export
setGeneric("fitted", function(object, ...)
  standardGeneric("fitted"))

#' @export
setGeneric("residuals", function(object, ...)
  standardGeneric("residuals"))

#' @export
setGeneric("nobs", function(object, ...)
  standardGeneric("nobs"))

#' @export
setGeneric("AIC", function(object, ..., k = 2)
  standardGeneric("AIC"))

#' @export
setGeneric("BIC", function(object, ...)
  standardGeneric("BIC"))

#' @title Plot predictions from FRK analysis 
#' @description This function acts as a wrapper around 
#' \code{\link{plot_spatial_or_ST}}. It plots the fields of the 
#' \code{Spatial*DataFrame} or \code{STFDF} object corresponding to 
#' prediction and prediction uncertainty quantification. It also uses the 
#' \code{@data} slot of \code{SRE} object to plot the training data set(s), 
#' and generates informative, latex-style legend labels for each of the plots. 
#' @param x object of class \code{SRE} 
#' @param y the \code{Spatial*DataFrame} or \code{STFDF} object resulting from the call \code{predict(x)}. 
#' Keep in mind that \code{predict()} returns a \code{list} when \code{method} = "TMB"; the element \code{$newdata} contains the required \code{Spatial}/\code{ST} object. 
#' If the list itself is passed, you will receive the error: "x" and "y" lengths differ.  
#' @param ... optional arguments passed on to \code{\link{plot_spatial_or_ST}}
#' @return A list of \code{ggplot} objects consisting of the observed data, predictions, and standard errors. This list can then be supplied to, for example, \code{ggpubr::ggarrange()}.
#' @export
#' @examples 
#' ## See example in the help file for SRE
setGeneric("plot", function(x, y, ...)
    standardGeneric("plot"))


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
#' @param sum_variables if \code{average_in_BAU == TRUE}, the string \code{sum_variables} indicates which data variables (can be observations or covariates) are to be summed rather than averaged
#' @details This generic function is not called directly. It is called in the SRE function for binning data in BAUs
#' @noRd
setGeneric("map_data_to_BAUs", function(data_sp,sp_pols,variogram.formula=NULL,est_error=T,average_in_BAU = TRUE, sum_variables = NULL, silently = FALSE) standardGeneric("map_data_to_BAUs"))

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

#' @title Retrieve the binned data. 
#' @description Takes an object of class \code{SRE} and returns the observed data binned into the BAUs. The binned data will feature many missing values in the case of point-referenced data.  The binned data is in the same order as the BAUs, facilitating plotting. 
#' @param object object of class \code{SRE}
#' @noRd
setGeneric("binned_data", function(object)
  standardGeneric("binned_data"))
