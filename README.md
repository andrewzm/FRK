[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](http://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/FRK?branch=master)


Fixed Rank Kriging
================

This package is operational but still in its early stages of development. It is based on the work by Cressie and Johannesson (2008) and carries out using Fixed Rank Kriging (FRK) using large datasets. A special feature of this package is that it allows for a Hadoop backend, allowing for the construction of massive L3 products.

Package: FRK

Type: Package

Title: Fixed Rank Kriging

Version: 0.1.0

Date: 2015-05-15

Author: Andrew Zammit-Mangion

Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>

VignetteBuilder: knitr

Description: This package implements fixed rank kriging, a tool used for spatial modelling and prediction with large detasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of n basis functions, where n is typically much smaller than the number of data points (or polygons) m. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The projected field is typically the integral component of the spatial random effects (SRE) model, the central component of this package. The package FRK provides a means to model, fit and predict with the SRE with relative ease. The work is based on Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70, 209-226.

BugReports: http://github.com/andrewzm/FRK/issues

Depends:
    R (>= 3.2)

Suggests:
    testthat,
    reshape2,
    shapefiles,
    knitr,
    zoo,
    INLA,
    splancs

Enhances:
    Rhipe

Imports:
    methods,
    Matrix,
    digest,
    gstat,
    ggplot2,
    dplyr,
    scales,
    grid,
    fields,
    plyr,
    PBSmapping,
    parallel,
    doParallel,
    sp,
    mapproj

License: GPL (>= 2)

NeedsCompilation: yes

LazyData: true


