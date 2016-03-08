[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](http://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/FRK?branch=master)


Fixed Rank Kriging
================

The package FRK is still under development. Both function interfaces and the underlying modelling approach are frequently being changed. This notice will be removed once a stable version is released.



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

`FRK` is still in its early stages of development. To install, first please install `INLA` from `http://www.r-inla.org/download`, then please load `devtools` and type

    install_github("andrewzm/FRK",dependencies=TRUE,build_vignettes=TRUE)

Installation will take a few minutes since the vignette is extensive and built from scratch. After installation see the  several examples provided in the vignette by typing

    vignette("FRK_intro")

Known Issues
------------

Currently `FRK` is not installing on OSX with `build_vignettes=TRUE` as it fails to find `texi2dvi`. Set `build_vignettes=FALSE` to ensure installation. Then download the `.Rnw` file in the `vignettes` folder and compile the pdf file separately in `RStudio` with `knitr`. 