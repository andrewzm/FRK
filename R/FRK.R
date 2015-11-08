#' Fixed Rank Kriging
#'
#' This package implements fixed rank kriging, a tool used for spatial modelling and prediction with large detasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of n basis functions, where n is typically much smaller than the number of data points (or polygons) m. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The projected field is typically the integral component of the spatial random effects (SRE) model, the central component of this package. The package FRK provides a means to model, fit and predict with the SRE with relative ease. The work is based on Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70, 209-226.
#' @name FRK-package
#' @docType package
#' @useDynLib FRK
#' @import methods
#' @import ggplot2
#' @import Matrix
#' @import sp
#' @import spacetime
#' @import parallel
#' @import gstat
#' @import dplyr
#' @import mapproj
#' @import foreach
#' @importFrom plyr ddply dlply rbind.fill
#' @importFrom PBSmapping clipPolys
#' @importFrom digest digest
#' @importFrom scales brewer_pal gradient_n_pal muted
#' @importFrom grid unit grid.newpage pushViewport viewport grid.layout
#' @importFrom fields rdist fields.rdist.near rdist.earth
#' @importFrom doParallel registerDoParallel
#' @name FRK
NULL
