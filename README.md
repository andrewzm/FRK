[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](http://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/FRK?branch=master)


Fixed Rank Kriging
================


Description
------------

Package: FRK

Type: Package

Title: Fixed Rank Kriging

Version: 0.1.0

Date: 2016-01-21

Author: Andrew Zammit-Mangion

Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>

Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of n basis functions, where *n* is typically much smaller than the number of data points (or polygons) *m*. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The projected field is a key building block of the spatial random effects (SRE) model, on which this package is based. The package FRK provides  helper functions to model, fit, and predict using an SRE with relative ease. Reference: Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B, 70, 209-226.


License: GPL (>= 2)

Installation 
------------

This package is operational but still in its early stages of development. It is based on the work by Cressie and Johannesson (2008) and carries out using Fixed Rank Kriging (FRK) using large datasets. To install, please isntall `devtools` and type

    install_github("andrewzm/FRK",dependencies=TRUE,build_vignettes=FALSE)


