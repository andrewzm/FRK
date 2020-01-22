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

#' Fixed Rank Kriging
#'
#' Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of n basis functions, where dimension n is typically much smaller than the number of data points (or polygons) m. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The dimension-reduced field is a key building block of the Spatial Random Effects (SRE) model, upon which this package is based. The package FRK provides helper functions to model, fit, and predict using an SRE with relative ease. Reference: Cressie, N. and Johannesson, G. (2008) <DOI:10.1111/j.1467-9868.2007.00633.x>.
#' @name FRK-package
#' @docType package
#' @useDynLib FRK, .registration=TRUE
#' @import methods
#' @import ggplot2
#' @import Matrix
#' @import TMB
#' @import RcppEigen
#' @import Rcpp
#' @import sp
#' @import spacetime
#' @import parallel
#' @import dplyr
#' @import sparseinv
#' @importFrom Hmisc roundPOSIXt truncPOSIXt ceil
#' @importFrom plyr ddply dlply rbind.fill
#' @importFrom digest digest
#' @importFrom Rcpp cppFunction
#' @importFrom grDevices chull
#' @importFrom stats .getXlevels coefficients dist kmeans lm median model.extract model.frame model.matrix na.fail optim runif sd terms var time rnorm formula
#' @importFrom utils data
#' @importFrom statmod rinvgauss
NULL
