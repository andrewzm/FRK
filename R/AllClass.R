#### CLASS DEFINITIONS ######

### Manifolds and measures
#'  @docType class
#'  @title measure
#'
#' @description Measure class used for defining measures used to compute distances between points in objects constructed with the \code{FRK} package.
#' @details An object of class \code{measure} contains a distance function and a variable \code{dim} with the dimensions of the Riemannian manifold over which the distance is computed. By default, the distance function used is inherited from the package \code{fields}.
#' @keywords Manifolds, spheres, planes
setClass("measure",representation(dist="function",dim="integer"),
         prototype(dist=dist,dim=2L))

#'  @docType class
#'  @title manifold
#' @description The class \code{manifold} is virtual; other manifold classes inherit from this class.
#' @details A \code{manifold} object is characterised by a character variable \code{type} which contains a description of the manifold, and a variable \code{metric} of type \code{measure}. A typical measure is the Euclidean distance.
#'
#' \code{FRK} supports three manifolds, the real line (in one dimension), instantiated using \code{real_line()}, the 2D plane instantiated using \code{plane()}, and the S2 sphere, instantiated using \code{sphere()}. User-specific manifolds can also be specified, however helper functions that are manifold specific, such as \code{auto_BAU} and \code{auto_BAUs} and \code{auto_basis} only work with the three pre-configured manifolds. On the other hand, there one can always change the distance function used on the manifold to synthesise anisotropy or heterogeneity, see vignettes for an example.
#'
#' @keywords Manifolds, spheres, planes
setClass("manifold",representation(type="character", metric = "measure","VIRTUAL"))

#' @rdname manifold-class
#' @aliases STmanifold-class
setClass("STmanifold",contains="manifold")



#' @rdname manifold-class
#' @aliases sphere-class
setClass("sphere",representation(radius="numeric"),contains="manifold")

#' @rdname manifold-class
#' @aliases STsphere-class
setClass("STsphere",representation(radius="numeric"),contains="STmanifold")


#' @rdname manifold-class
#' @aliases plane-class
setClass("plane",contains="manifold")

#' @rdname manifold-class
#' @aliases STplane-class
setClass("STplane",contains="STmanifold")

#' @rdname manifold-class
#' @aliases real_line-class
setClass("real_line",contains="manifold")

#' @rdname manifold-class
#' @aliases timeline-class
setClass("timeline",contains="manifold")

####  Basis functions ####
#'  @docType class
#'  @title Basis functions
#'
#' @description An object of class \code{Basis} contains the basis functions used to construct the \eqn{S} matrix in fixed-rank kriging. It contains five slots, described below.
#' @slot manifold an object of class \code{manifold} that contains information on the manifold and the distance metric used on the manifold. See \code{\link{manifold-class}} for more details
#' @slot n  the number of basis functions in this set
#' @slot fn a list of length \code{n}, with each item the function of a specific basis function
#' @slot pars a list of parameters where the \eqn{i}-th item in the list contains the parameters of the \eqn{i}-th basis function, \code{fn[[i]]}
#' @slot df a data frame containing other attributes specific to each basis function (for example the geometric centre of the local basis function)
#' @details Basis functions are a critical component of the \code{FRK} package and the package is designed to work with user-defined specifications of these. For convenience however, several functions aid the user to construct specific sets for a given set of data points. Please see \code{\link{auto_basis}} for more details. The function \code{\link{local_basis}} helps the user construct a set of local basis functions from a collection of locations and scale parameters.
#' @keywords Basis functions
#' @rdname Basisclass
setClass("Basis_obj", representation(n = "numeric","VIRTUAL"))

#' @rdname Basisclass
setClass("Basis",contains="Basis_obj", representation(manifold="manifold",fn="list",pars="list", df="data.frame"))

#' @rdname Basisclass
setClass("TensorP_Basis", contains="Basis_obj", representation(Basis1="Basis",Basis2="Basis"))


#'  @title Spatial Random Effects class
#' @description Central object of package, containing the model and all other information required for estimation and prediction.
#' @details The spatial random effects (SRE) model is the model employed in fixed rank kriging, and the \code{SRE} object contains all information required for estimation and prediction from spatial data. Object slots contain both other objects (for example, an object of class \code{Basis}) and matrices derived from these objects (for example, the matrix \eqn{S}) in order to facilitate computations.
#'
#'@slot data the original data use to condition upon when training the model
#'@slot basis object of class \code{Basis} used to construct the matrix \eqn{S}
#'@slot BAUs object of class \code{SpatialPolygonsDataFrame} that contains the basic aerial units (BAUs) that are used to both (i) project the data onto a common discretisation if they are point-referenced and (ii) provide a BAU to data relationship if the data has a spatial footprint
#' @slot f formula used to define the SRE object. All covariates employed need to be specified in the object \code{BAUs}
#' @slot S matrix constructed by evaluating the basis functions at all BAUs affected by the data
#' @slot Ve observation variance-covariance matrix (typically diagonal)
#' @slot Vfs fine-scale variance-covariance matrix (typically digonal) up to a constant of proportionality estimated in the framework
#' @slot Z vector of observations (of class \code{Matirx})
#' @slot Cmat incidence matrix mapping the observations to the BAUs
#' @slot X matrix of covariates
#' @slot mu_eta updated expectation of random effects (estimated)
#' @slot S_eta updated covariance matrix of random effects (estimated)
#' @slot Khat prior covariance matrix of random effects (estimated)
#' @slot alphahat covariates weights (estimated)
#' @slot sigma2fshat fine-scale variation scaler (estimated)
#' @keywords Spatial random effects, fixed rank kriging
setClass("SRE",representation(data="list",
                              basis="Basis_obj",
                              BAUs="ANY",     # should be SpatialPolygonsDataFrame of STFDF
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
