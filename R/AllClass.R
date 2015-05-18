#### CLASS DEFINITIONS ######

### Manifolds and measures
#'  @docType class
#'  @title measure
#'
#' @description Measure function
#'
#' @keywords Manifolds, spheres, planes
setClass("measure",representation(dist="function",dim="integer"),
         prototype(dist=fields::rdist,dim=2L))

#'  @docType class
#'  @title manifold
#'
#' @description Manifolds supported by the package.
#'
#' @keywords Manifolds, spheres, planes
setClass("manifold",representation(type="character", metric = "measure","VIRTUAL"))

#'  @docType class
#'  @title sphere
#' @description Sphere manifold
#'
#' @keywords Manifolds, spheres, planes
setClass("sphere",representation(radius="numeric"),contains="manifold")

#'  @docType class
#'  @title plane
#'
#' @description Plane
#'
#' @keywords Manifolds, spheres, planes
setClass("plane",contains="manifold")

#'  @docType class
#'  @title real_line
#'
#' @description 1D real line
#'
#' @keywords Manifolds, spheres, planes, lines
setClass("real_line",contains="manifold")



####  Basis functions ####
#'  @docType class
#'  @title Parent class for basis functions
#'
#' @description A basis function contains three slots, the first is a list of parameters (pars), the second is the number (n) of basis functions and the third is a list
#' of functions (fn). The pars field usually contains parameters which appear in fn.
#'
#' @keywords Basis functions
#' @rdname Basisclass
setClass("Basis", representation(manifold="manifold",n = "numeric",fn="list",pars="list", df="data.frame"))

#'  @docType class
#'  @title Gaussian radial basis functions basis
#'
#' @description GRBFBasis inherits from the virtual class Basis. A GRBF basis functions is initialsied using the function \code{initGRBFbasis}.
#'
#' @keywords Basis functions, GRBF
setClass("GRBFBasis",contains="Basis")

#'  @title Bisquare radial basis functions basis
#' @description bisquareBasis inherits from the virtual class Basis. A bisquare basis functions is initialsied using the function \code{initbisquarebasis}.
#'
#' @keywords Basis functions, GRBF
setClass("bisquareBasis",contains="Basis")


#'  @title SRE model
#' @description SRE model
#'
#' @keywords Spatial random effects, fixed rank kriging
setClass("SRE",representation(data="list",
                              basis="Basis",
                              BAUs="SpatialPolygonsDataFrame",
                              f = "formula",
                              S = "Matrix",
                              Ve = "Matrix",
                              Vfs = "Matrix",
                              Z = "Matrix",
                              Cmat = "Matrix",
                              X = "Matrix",
                              mu_eta = "Matrix",
                              S_eta = "Matrix",
                              Khat = "Matrix",
                              alphahat = "Matrix",
                              sigma2fshat = "numeric"))
