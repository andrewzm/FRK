#### CLASS DEFINITIONS ######

### Manifolds and measures
#'  @docType class
#'  @title measure
#'
#' @description Measure class used for defining measures used to compute distances between points in objects constructed with the \code{FRK} package.
#' @details An object of class \code{measure} contains a distance function and a variable \code{dim} with the dimensions of the Riemannian manifold over which the distance is computed. By default, distance functions used are those extracted from the package \code{fields}.
#' @keywords Manifolds, spheres, planes
setClass("measure",representation(dist="function",dim="integer"),
         prototype(dist=dist,dim=2L))

#'  @docType class
#'  @title manifold
#' @description The class \code{manifold} is virtual; other manifold classes inherit from this class.
#' @details A \code{manifold} object is characterised by a character variable \code{type}, which contains a description of the manifold, and a variable \code{measure} of type \code{measure}. A typical measure is the Euclidean distance.
#'
#' \code{FRK} supports five manifolds; the real line (in one dimension), instantiated by using \code{real_line()}; the 2D plane, instantiated by using \code{plane()}; the 2D-sphere surface S2, instantiated by using \code{sphere()}; the R2 space-time manifold, instantiated by using \code{STplane()}, and the S2 space-time manifold, instantiated by using \code{STsphere()}. User-specific manifolds can also be specified, however helper functions that are manifold specific, such as \code{auto_BAU} and \code{auto_BAUs} and \code{auto_basis} only work with the pre-configured manifolds. Importantly, one can change the distance function used on the manifold to synthesise anisotropy or heterogeneity.
#'
#' @keywords Manifolds, spheres, planes
setClass("manifold",representation(type="character", measure = "measure","VIRTUAL"))

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
#' @docType class
#' @title Basis functions
#'
#' @description An object of class \code{Basis} contains the basis functions used to construct the matrix \eqn{S} in fixed-rank kriging. It contains five slots, described below.
#' @slot manifold an object of class \code{manifold} that contains information on the manifold and the distance measure used on the manifold. See \code{\link{manifold-class}} for more details
#' @slot n  the number of basis functions in this set
#' @slot fn a list of length \code{n}, with each item the function of a specific basis function
#' @slot pars a list of parameters where the \eqn{i}-th item in the list contains the parameters of the \eqn{i}-th basis function, \code{fn[[i]]}
#' @slot df a data frame containing other attributes specific to each basis function (for example the geometric centre of the local basis function)
#' @details Basis functions are a central component of \code{FRK}, and the package is designed to work with user-defined specifications of these. For convenience, however, several functions are available to aid the user to construct a basis set for a given set of data points. Please see \code{\link{auto_basis}} for more details. The function \code{\link{local_basis}} helps the user construct a set of local basis functions (e.g., bisquare functions) from a collection of locations and scale parameters.
#' @keywords Basis functions
#' @rdname Basisclass
setClass("Basis_obj", representation(n = "numeric","VIRTUAL"))

#' @rdname Basisclass
setClass("Basis",contains="Basis_obj", representation(manifold="manifold",fn="list",pars="list", df="data.frame"))

#' @rdname Basisclass
setClass("TensorP_Basis", contains="Basis_obj", representation(Basis1="Basis",Basis2="Basis"))


#' @title Spatial Random Effects class
#' @description This is the central class definition of the \code{FRK} package, containing the model and all other information required for estimation and prediction.
#' @details The spatial random effects (SRE) model is the model employed in Fixed Rank Kriging, and the \code{SRE} object contains all information required for estimation and prediction from spatial data. Object slots contain both other objects (for example, an object of class \code{Basis}) and matrices derived from these objects (for example, the matrix \eqn{S}) in order to facilitate computations.
#'
#'@slot data the original data from which the model's parameters are estimated
#'@slot basis object of class \code{Basis} used to construct the matrix \eqn{S}
#'@slot BAUs object of class \code{SpatialPolygonsDataFrame} that contains the Basic Areal Units (BAUs) that are used to both (i) project the data onto a common discretisation if they are point-referenced and (ii) provide a BAU-to-data relationship if the data has a spatial footprint
#' @slot f formula used to define the SRE object. All covariates employed need to be specified in the object \code{BAUs}
#' @slot S matrix constructed by evaluating the basis functions at all BAUs affected by the data
#' @slot Ve measurement-error variance-covariance matrix (typically diagonal)
#' @slot Vfs fine-scale variance-covariance matrix (typically diagonal) up to a constant of proportionality estimated in the framework
#' @slot Z vector of observations (of class \code{Matrix})
#' @slot Cmat incidence matrix mapping the observations to the BAUs
#' @slot X matrix of covariates
#' @slot mu_eta updated expectation of random effects (estimated)
#' @slot S_eta updated covariance matrix of random effects (estimated)
#' @slot Khat prior covariance matrix of random effects (estimated)
#' @slot Khat_inv prior covariance matrix of random effects (estimated)
#' @slot alphahat covariates weights (estimated)
#' @slot sigma2fshat fine-scale variation scaler (estimated)
#' @keywords Spatial random effects, fixed rank kriging
setClass("SRE",representation(data="list",
                              basis="Basis_obj",
                              BAUs="ANY",     # should be SpatialPolygonsDataFrame or STFDF
                              f = "formula",
                              S = "Matrix",
                              Ve = "Matrix",
                              Vfs = "Matrix",
                              Vfs_BAUs = "Matrix",
                              Qfs_BAUs = "Matrix",
                              Z = "Matrix",
                              Cmat = "Matrix",
                              X = "Matrix",
                              mu_eta = "Matrix",
                              mu_xi = "Matrix",
                              S_eta = "Matrix",
                              Khat = "Matrix",
                              alphahat = "Matrix",
                              Q_eta = "Matrix",
                              Khat_inv = "Matrix",
                              B_run = "Matrix",
                              v_run = "Matrix",
                              sigma2fshat = "numeric",
                              fs_model = "character",
                              regularise = "logical"))
