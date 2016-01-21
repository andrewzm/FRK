#' Fixed Rank Kriging
#'
#' Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of n basis functions, where n is typically much smaller than the number of data points (or polygons) m. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The projected field is a key building block of the spatial random effects (SRE) model, on which this package is based. The package FRK provides  helper functions to model, fit, and predict using an SRE with relative ease. Reference: Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B, 70, 209-226..
#' @name FRK-package
#' @docType package
#' @useDynLib FRK
#' @import ggplot2
#' @import Matrix
#' @import sp
#' @import spacetime
#' @import parallel
#' @import dplyr
#' @importFrom plyr ddply dlply rbind.fill
#' @importFrom digest digest
#' @name FRK
NULL
