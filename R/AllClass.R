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

#### CLASS DEFINITIONS ######

### Manifolds and measures
#' @docType class
#' @title measure
#'
#' @description Measure class used for defining measures used to compute distances between points in objects constructed with the \code{FRK} package.
#' @details An object of class \code{measure} contains a distance function and a variable \code{dim} with the dimensions of the Riemannian manifold over which the distance is computed.
#' @seealso \code{\link{distance}} for computing a distance and \code{\link{distances}} for a list of implemented distance functions.
setClass("measure",representation(dist="function",dim="integer"),
         prototype(dist=dist,dim=2L))

#' @docType class
#' @title manifold
#' @description The class \code{manifold} is virtual; other manifold classes inherit from this class.
#' @details A \code{manifold} object is characterised by a character variable \code{type}, which contains a description of the manifold, and a variable \code{measure} of type \code{measure}. A typical measure is the Euclidean distance.
#'
#' \code{FRK} supports five manifolds; the real line (in one dimension), instantiated by using \code{real_line()}; the 2D plane, instantiated by using \code{plane()}; the 2D-sphere surface S2, instantiated by using \code{sphere()}; the R2 space-time manifold, instantiated by using \code{STplane()}, and the S2 space-time manifold, instantiated by using \code{STsphere()}. User-specific manifolds can also be specified, however helper functions that are manifold specific, such as \code{auto_BAUs} and \code{auto_basis}, only work with the pre-configured manifolds. Importantly, one can change the distance function used on the manifold to synthesise anisotropy or heterogeneity. See the vignette for one such example.
#' @seealso \code{\link{real_line}}, \code{\link{plane}}, \code{\link{sphere}}, \code{\link{STplane}} and \code{\link{STsphere}} for constructing manifolds.
setClass("manifold",representation(type="character", measure = "measure","VIRTUAL"))

#' @rdname manifold-class
#' @aliases STmanifold-class
setClass("STmanifold",contains="manifold")

#' @rdname manifold-class
#' @aliases real_line-class
setClass("real_line",contains="manifold")

#' @rdname manifold-class
#' @aliases plane-class
setClass("plane",contains="manifold")

#' @rdname manifold-class
#' @aliases STplane-class
setClass("STplane",contains="STmanifold")

#' @rdname manifold-class
#' @aliases sphere-class
setClass("sphere",representation(radius="numeric"),contains="manifold")

#' @rdname manifold-class
#' @aliases STsphere-class
setClass("STsphere",representation(radius="numeric"),contains="STmanifold")

####  Basis functions ####
#' @docType class
#' @title Basis functions
#'
#' @description An object of class \code{Basis} contains the basis functions used to construct the matrix \eqn{S} in FRK. 
#' @slot manifold an object of class \code{manifold} that contains information on the manifold and the distance measure used on the manifold. See \code{\link{manifold-class}} for more details
#' @slot n  the number of basis functions in this set
#' @slot fn a list of length \code{n}, with each item the function of a specific basis function
#' @slot pars a list of parameters where the \eqn{i}-th item in the list contains the parameters of the \eqn{i}-th basis function, \code{fn[[i]]}
#' @slot df a data frame containing other attributes specific to each basis function (for example the geometric centre of the local basis function)
#' @slot regular logical indicating if the basis functions (of each resolution) are in a regular grid
#' @details Basis functions are a central component of \code{FRK}, and the package is designed to work with user-defined specifications of these. For convenience, however, several functions are available to aid the user to construct a basis set for a given set of data points. Please see \code{\link{auto_basis}} for more details. The function \code{\link{local_basis}} helps the user construct a set of local basis functions (e.g., bisquare functions) from a collection of location and scale parameters.
#' @seealso \code{\link{auto_basis}} for automatically constructing basis functions and \code{\link{show_basis}} for visualising basis functions.
#' @rdname Basis-class
setClass("Basis_obj", representation(n = "numeric","VIRTUAL"))

#' @rdname Basis-class
setClass("Basis",contains="Basis_obj", representation(manifold="manifold",fn="list",pars="list", df="data.frame", regular = "numeric"))

#' @rdname Basis-class
setClass("TensorP_Basis", contains="Basis_obj", representation(Basis1="Basis",Basis2="Basis",n = "integer", df = "data.frame", regular="logical"))

#' @title Spatial Random Effects class
#' @description This is the central class definition of the \code{FRK} package, containing the model and all other information required for estimation and prediction.
#' @details The spatial random effects (SRE) model is the model employed in Fixed Rank Kriging, and the \code{SRE} object contains all information required for estimation and prediction from spatial data. Object slots contain both other objects (for example, an object of class \code{Basis}) and matrices derived from these objects (for example, the matrix \eqn{S}) in order to facilitate computations.
#'
#' @slot f formula used to define the SRE object. All covariates employed need to be specified in the object \code{BAUs}
#' @slot data the original data from which the model's parameters are estimated
#' @slot basis object of class \code{Basis} used to construct the matrix \eqn{S}
#' @slot BAUs object of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame} of \code{STFDF} that contains the Basic Areal Units (BAUs) that are used to both (i) project the data onto a common discretisation if they are point-referenced and (ii) provide a BAU-to-data relationship if the data has a spatial footprint
#' @slot S matrix constructed by evaluating the basis functions at all the data locations (of class \code{Matrix})
#' @slot S0 matrix constructed by evaluating the basis functions at all BAUs (of class \code{Matrix})
#' @slot D_basis list of distance-matrices of class \code{Matrix}, one for each basis-function resolution
#' @slot Ve measurement-error variance-covariance matrix (typically diagonal and of class \code{Matrix})
#' @slot Vfs fine-scale variance-covariance matrix at the data locations (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated using the EM algorithm
#' @slot Vfs_BAUs fine-scale variance-covariance matrix at the BAU centroids (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated using the EM algorithm
#' @slot Qfs_BAUs fine-scale precision matrix at the BAU centroids (typically diagonal and of class \code{Matrix}) up to a constant of proportionality estimated using the EM algorithm
#' @slot Z vector of observations (of class \code{Matrix})
#' @slot Cmat incidence matrix mapping the observations to the BAUs
#' @slot X design matrix of covariates at all the data locations
#' @slot G list of objects of class Matrix containing the design matrices for random effects at all the data locations
#' @slot G0 list of objects of class Matrix containing the design matrices for random effects at all BAUs
#' @slot K_type type of prior covariance matrix of random effects. Can be "block-exponential" (correlation between effects decays as a function of distance between the basis-function centroids), "unstructured" (all elements in \code{K} are unknown and need to be estimated), or "neighbour" (a sparse precision matrix is used, whereby only neighbouring basis functions have non-zero precision matrix elements).
#' @slot mu_eta updated expectation of the basis-function random effects (estimated)
#' @slot mu_gamma updated expectation of the random effects (estimated)
#' @slot S_eta updated covariance matrix of random effects (estimated)
#' @slot Q_eta updated precision matrix of random effects (estimated)
#' @slot Khat prior covariance matrix of random effects (estimated)
#' @slot Khat_inv prior precision matrix of random effects (estimated)
#' @slot alphahat fixed-effect regression coefficients (estimated)
#' @slot sigma2fshat fine-scale variation scaling (estimated)
#' @slot sigma2gamma random-effect variance parameters (estimated)
#' @slot fs_model type of fine-scale variation (independent or CAR-based). Currently only "ind" is permitted
#' @slot info_fit information on fitting (convergence etc.)
#' @slot response A character string indicating the assumed distribution of the response variable
#' @slot link A character string indicating the desired link function. Can be "log", "identity", "logit", "probit", "cloglog", "reciprocal", or "reciprocal-squared". Note that only sensible link-function and response-distribution combinations are permitted. 
#' @slot mu_xi updated expectation of the fine-scale random effects at all BAUs (estimated)
#' @slot Q_posterior updated joint precision matrix of the basis function random effects and observed fine-scale random effects (estimated)
#' @slot log_likelihood the log likelihood of the fitted model
#' @slot method the fitting procedure used to fit the SRE model
#' @slot phi the estimated dispersion parameter (assumed constant throughout the spatial domain)
#' @slot k_Z vector of known size parameters at the observation support level (only applicable to binomial and negative-binomial response distributions) 
#' @slot k_BAU vector of known size parameters at the observed BAUs (only applicable to binomial and negative-binomial response distributions) 
#' @slot include_fs flag indicating whether the fine-scale variation should be included in the model
#' @slot include_gamma flag indicating whether there are gamma random effects in the model 
#' @slot normalise_wts if \code{TRUE}, the rows of the incidence matrices \eqn{C_Z} and \eqn{C_P} are normalised to sum to 1, so that the mapping represents a weighted average; if false, no normalisation of the weights occurs (i.e., the mapping corresponds to a weighted sum)
#' @slot fs_by_spatial_BAU if \code{TRUE}, then each BAU is associated with its own fine-scale variance parameter
#' @slot obsidx indices of observed BAUs
#' @slot simple_kriging_fixed logical indicating whether one wishes to commit to simple kriging at the fitting stage: If \code{TRUE}, model fitting is faster, but the option to conduct universal kriging at the prediction stage is removed 
#' @seealso \code{\link{SRE}} for details on how to construct and fit SRE models.
#' @references
#' Zammit-Mangion, A. and Cressie, N. (2017). FRK: An R package for spatial and spatio-temporal prediction with large datasets. Journal of Statistical Software, 98(4), 1-48. doi:10.18637/jss.v098.i04.
setClass("SRE",representation(data="list",
                              basis="Basis_obj",
                              BAUs="ANY",     # should be SpatialPolygonsDataFrame, SpatialPixelsDataFrame or STFDF
                              f = "formula",
                              S = "Matrix",
                              S0 = "Matrix",
                              D_basis = "list",
                              Ve = "Matrix",
                              Vfs = "Matrix",
                              Vfs_BAUs = "Matrix",
                              Qfs_BAUs = "Matrix",
                              Z = "Matrix",
                              Cmat = "Matrix",
                              X = "Matrix",
                              G = "list",
                              G0 = "list",
                              mu_eta = "Matrix",
                              mu_gamma = "Matrix",
                              S_eta = "Matrix",
                              Q_eta = "Matrix",
                              K_type = "character",
                              Khat = "Matrix",
                              Khat_inv = "Matrix",
                              alphahat = "Matrix",
                              sigma2fshat = "numeric",
                              sigma2gamma = "numeric",
                              fs_model = "character",
                              info_fit = "list", 
                              response = "character", 
                              link = "character", 
                              mu_xi = "Matrix",
                              Q_posterior = "dsCMatrix",
                              log_likelihood = "numeric", 
                              method = "character",
                              phi = "numeric", 
                              k_Z = "numeric", 
                              k_BAU_O = "numeric", 
                              include_fs = "logical", 
                              include_gamma = "logical", 
                              normalise_wts = "logical", 
                              fs_by_spatial_BAU = "logical", 
                              obsidx = "numeric", 
                              simple_kriging_fixed = "logical"))
