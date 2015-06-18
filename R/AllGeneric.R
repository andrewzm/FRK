#### GENERIC FUNCTIONS ######
#' @title Show basis functions
#' @description Generic plotting function for illustrating the basis functions.
#' @param g object of class \code{gg} (a \code{ggplot} object)
#' @param basis object of class \code{Basis}
#' @details The function \code{show_basis} adapts its behaviour to the manifold being used. With \code{real_line}, the one-dimensional basis functions are plotted with colour distinguishing between the different resolutions. With \code{plane}, only radial basis functions are supported (at present). Each basis function is shown as a circle with diameter equal to the \code{scale} parameter of the function. Linetype distinguishes the resolution. With \code{sphere}, the centres of the basis functions are shown as circles, with larger sizes corresponding to lower (i.e., coarser) resolutions.
#' @examples
#' library(ggplot2)
#' library(sp)
#' data(meuse)
#' coordinates(meuse) = ~x+y # change into an sp object
#' G <- auto_basis(m = plane(),data=meuse,nres = 2,regular=2,prune=10,type = "Gaussian")
#' show_basis(ggplot(),G) + geom_point(data=data.frame(meuse),aes(x,y))
#' @export
setGeneric("show_basis", function(g,basis) standardGeneric("show_basis"))

#' @title Retrieve manifold
#' @description Retrieve manifold from \code{FRK} object.
#' @param .Object \code{FRK} object
#' @examples
#' G <-  radial_basis(manifold = plane(),
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
#' G <- auto_basis(m = plane(),data=meuse,nres = 2,regular=2,prune=10,type = "Gaussian")
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
#' @description Compute distance using object of class \code{measure} or \code{manifold}
#' @param d object of class \code{measure} or \code{manifold}
#' @param x1 first coordinate
#' @param x2 second coordinate
#' @examples
#' distance(sphere(),matrix(0,1,2),matrix(10,1,2))
#' distance(plane(),matrix(0,1,2),matrix(10,1,2))
#' @export
setGeneric("distance", function(d,x1,x2) standardGeneric("distance"))

#' @title Evaluate basis functions
#' @description  Evaluate basis functions at points or average functions over polygons
#' @param basis object of class \code{Basis}
#' @param s object of class \code{matrix}, \code{SpatialPointsDataFrame} or \code{SpatialPolygonsDataFrame}
#' @param output either a "list" or "matrix", depending on desired output format.
#' @details This function evaluates the basis functions at isolated points, or averages the basis functions over polygons, for computing the matrix \eqn{S}. The latter operation is carried out using Monte Carlo integration with 1000 samples per polygon. This function is embarassingly parallelisable and is thus Hadoop-enabled (provided a Hadoop backend is available and \code{opts_FRK$get("Rhipe")} is TRUE).
#' @examples
#' library(sp)
#' data(meuse)
#' coordinates(meuse) = ~x+y # change into an sp object
#' G <- auto_basis(m = plane(),data=meuse,nres = 2,regular=2,prune=10,type = "Gaussian")
#' S <- eval_basis(G,matrix(c(180000,332000),1,2))
#' @export
setGeneric("eval_basis", function(basis,s,output="matrix") standardGeneric("eval_basis"))


#### NOT EXPORTED ####

#' @title Automatic BAU generation
#' @noRd
#' @description This generic function is called by \code{auto_BAUs} after a series of checks.
#' @param manifold object of class \code{manifold}
#' @param cellsize denotes size of gridcell when \code{type} == ``grid''
#' @param resl resolution number of isea3h DGGRID cells for when \code{type} is ``hex'' and \code{manifold} is a sphere
#' @param type either ``hex'' or ``grid'', indicating whether gridded or hexagonal BAUs should be used
#' @param d data, that is, an object of class SpatialPointsDataFrame or SpatialPolygonsDataFrame. Provision of data implies that the domain is bounded (necessary with \code{real_line} and \code{plane} but not necessary with \code{sphere})
#' @param convex convex parameter used for smoothing an extended boundary when working on a finite domain (that is, when the object \code{d} is supplied), see details
#' @param ... currently unused
#' @details This generic function is not called directly. Please refer to \code{auto_BAUs} for more details.
setGeneric("auto_BAU", function(manifold,cellsize,resl,type,d,convex,...) standardGeneric("auto_BAU"))


#' @title Concatenation
#' @description Concatenates FRK objects of the same class together. This is primarily used to join up \code{Basis} blocks together.
#' @param ... a series of \code{FRK} objects
#' @noRd
setGeneric("concat", function(...) standardGeneric("concat"))

