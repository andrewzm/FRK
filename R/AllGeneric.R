#### GENERIC FUNCTIONS ######
#' @title Show basis functions
#' @export
setGeneric("show_basis", function(g,basis,...) standardGeneric("show_basis"))

#' @title Automatic BAU generation
#' @export
setGeneric("auto_BAU", function(manifold,cellsize,res,type,data,...) standardGeneric("auto_BAU"))


#' @title Manifold
#' @export
setGeneric("manifold", function(.Object) standardGeneric("manifold"))

#' @title Number of basis functions
#' @export
setGeneric("nbasis", function(.Object) standardGeneric("nbasis"))

#' @title Type of manifold
#' @export
setGeneric("type", function(.Object) standardGeneric("type"))

#' @title Retrieve distance measure used
#' @export
setGeneric("distance", function(d,x1,x2) standardGeneric("distance"))

#' @title Evaluate basis functions at points or over polygons
#' @export
setGeneric("eval_basis", function(basis,s,output="list") standardGeneric("eval_basis"))

#' @title Concatenation
#' @description Concatenates MVST objects of the same class together. This is primarily used to join up \code{GMRF_basis} blocks and \code{Obs} blocks together.
#' @param ... a series of \code{MVST} objects
setGeneric("concat", function(...) standardGeneric("concat"))

