[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](http://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/FRK?branch=master)


Fixed Rank Kriging
================

The package FRK is still under development. However, both function interfaces and the underlying modelling approach are now relatively stable. This notice will be removed once a stable version is released.

Description
------------

Package: FRK

Type: Package

Title: Fixed Rank Kriging

Version: 0.1.0

Date: 2016-04-19

Author: Andrew Zammit-Mangion

Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>

Description: Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of *n* basis functions, where *n* is typically much smaller than the number of data points (or polygons) *m*. The method naturally allows for non-stationary, anisotropic covariance functions and the use of observations with varying support (with known error variance). The projected field is a
    key building block of the Spatial Random Effects (SRE) model, on which this package is based. The package FRK provides helper functions to model, fit, and predict using an SRE with relative ease. Reference: Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B, 70, 209-226.

License: GPL (>= 2)

Installation 
------------

`FRK` is still under development but stable and ready for use with simple applications. To install, first please install `INLA` from `http://www.r-inla.org/download`, then please load `devtools` and type

    install_github("andrewzm/FRK",dependencies=TRUE,build_vignettes=TRUE)

A journel article, currently under preparation, will be uploaded as vignette once submitted. If you wish to see a draft of the paper please contact me on azm@uow.edu.au directly. Thank you.

Quick start
------------

```r
    library(sp)
    library(ggplot2)
    library(FRK)
    
    ## Setup
    set.seed(1)                                               # Fix seed
    zdf <- Z <- data.frame(x = runif(1000), y= runif(1000))   # Generate random locs
    zdf$z <- Z$z <- sin(8*Z$x) + cos(8*Z$y) + 0.5*rnorm(100)  # Simulate data
    coordinates(Z) = ~x+y                                     # Turn into sp object
    
    ## Run FRK
    S <- FRK(f = z~1,                                         # Formula to FRK
                list(Z),                                      # All datasets are supplied in list
                n_EM = 10)                                    # Max number of EM iterations
    Pred <- SRE.predict(SRE_model = S)                        # Prediction stage
    
    xy <- data.frame(coordinates(Pred))                       # Extract info from predictions
    xy$mu <- Pred$mu
    xy$se <- Pred$sd
    
    ## Plotting
    ggplot(zdf) + geom_point(aes(x,y,colour=z)) + 
      scale_colour_distiller(palette="Spectral") + theme_bw() + coord_fixed()
    ggplot(xy) + geom_raster(aes(x,y,fill=mu)) + 
      scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()
    ggplot(xy) + geom_tile(aes(x,y,fill=se)) + 
       geom_point(data=zdf,aes(x,y),pch=46) +
       scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()
```    

![alt tag](https://dl.dropboxusercontent.com/u/3028804/FRK/FRK_ex_data.png)
![alt tag](https://dl.dropboxusercontent.com/u/3028804/FRK/FRK_ex_mu.png)
![alt tag](https://dl.dropboxusercontent.com/u/3028804/FRK/FRK_ex_se.png)

[//]: # ( > ggsave(gdata,file="~/Dropbox/Public/FRK/FRK_ex_data.png",width=7,height=5) > ggsave(gmu,file="~/Dropbox/Public/FRK/FRK_ex_mu.png",width=7,height=5) > ggsave(gse,file="~/Dropbox/Public/FRK/FRK_ex_se.png",width=7,height=5) )

Known Issues
------------

[//]: # (Currently `FRK` is not installing on OSX with `build_vignettes=TRUE` as it fails to find `texi2dvi`. Set `build_vignettes=FALSE` to ensure installation. Then download the `.Rnw` file in the `vignettes` folder and compile the pdf file separately in `RStudio` with `knitr`. )


Demonstrations
--------------

The package `FRK` is currently being used to generate spatio-temporal animations of fields observed by satellite data. [Here](https://www.youtube.com/watch?v=_kPa8VoeSdM) we show a daily prediction of CO2 using data from the NASA OCO-2 between September 2014 and June 2016.

[![alt tag](https://img.youtube.com/vi/ENx4CIZdoQk/0.jpg)](https://www.youtube.com/watch?v=ENx4CIZdoQk)
