[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](https://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](https://codecov.io/github/andrewzm/FRK?branch=master)

Fixed Rank Kriging 
================

# pkgdown <img src="man/figures/FRK_logo.svg" alt="" width="100"/>

Installation 
------------

The package `FRK` is now at v2.1.0 and available on CRAN! To install, please type

```r
install.packages("FRK")
```

To install the most recent (development) version, first please install `INLA` from `https://www.r-inla.org/download`, then please load `devtools` and type

```r
install_github("andrewzm/FRK", dependencies = TRUE, build_vignettes = TRUE)
```

A document containing a description, details on the underlying maths and computations for the EM algorithm (only applicable for Gaussian data), as well as several examples, is available as a vignette titled "FRK_intro". Another vignette, "FRK_non-Gaussian", summarises the maths in a non-Gaussian setting, and contains a couple of examples using non-Gaussian data and the newly available plotting methods. To load these vignettes please type

```r
library("FRK")
vignette("FRK_intro")
vignette("FRK_non-Gaussian")
```

A paper with more details is available [here](https://www.jstatsoft.org/article/view/v098i04), and a draft paper detailing the approach in a non-Gaussian setting is available [here](https://arxiv.org/abs/2110.02507). If you use `FRK` in your work, please cite it using the information provided by `citation("FRK")`. 


Description
------------

Package: FRK

Type: Package

Title: Fixed Rank Kriging

Version: 2.0.1

Date: 2021-05-26

Author: Andrew Zammit-Mangion, Matthew Sainsbury-Dale

Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>

Description: Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach decomposes the field, and hence the covariance function, using a fixed set of *r* basis functions, where *r* is typically much smaller than the number of data points (or polygons) *m*. This low-rank basis-function representation facilitates the modelling of 'big' spatial/spatio-temporal data. The method naturally allows for non-stationary, anisotropic covariance functions. Discretisation of the spatial domain into so-called basic areal units (BAUs) facilitates the use of observations with varying support (i.e., both point-referenced and areal supports, potentially simultaneously), and prediction over arbitrary user-specified regions. The package `FRK` caters for Gaussian and non-Gaussian data, employing a spatial generalised linear mixed model (GLMM) framework  to cater for Poisson, binomial, negative-binomial, gamma, and inverse-Gaussian distributed data. `FRK` also supports inference over various manifolds, including the 2D plane and 3D sphere, and it provides helper functions to model, fit, predict, and plot with relative ease. Zammit-Mangion and Cressie (2021) describe `FRK` in a Gaussian setting, 
and detail it's use of basis functions and BAUs. 

* Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B, 70, 209–226.
* Zammit-Mangion A, Cressie N (2021). “FRK: an R package for spatial and spatio-temporal prediction with large datasets.” Journal of Statistical Software, In press.


License: GPL (>= 2)


Quick start
------------

### Gaussian data


```r
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")

## Setup
m <- 1000                                                  # Sample size
RNGversion("3.6.0"); set.seed(1)                           # Fix seed
zdf <- data.frame(x = runif(m), y= runif(m))               # Generate random locs
zdf$z <- sin(8 * zdf$x) + cos(8 * zdf$y) + 0.5 * rnorm(m)  # Simulate data
coordinates(zdf) = ~x+y                                    # Turn into sp object

## Run FRK
S <- FRK(f = z ~ 1,                         # Formula to FRK
         list(zdf),                         # All datasets are supplied in list
         n_EM = 10)                         # Max number of EM iterations
pred <- predict(S)                          # Prediction stage

## Plotting
plot_list <- plot(S, pred)
ggarrange(plotlist = plot_list, nrow = 1, legend = "top")

```

<!---
ggsave( 
  filename = "Gaussian_data_tmp.png", device = "png", 
  width = 10, height = 4,
  path = "~/Dropbox/FRK/man/figures/"
)
--->

![(Left) Gaussian data. (Centre) Predictions. (Right) Standard errors.](/man/figures/Gaussian_data.png?raw=true)

### Non-Gaussian data

Here we analyse simulated Poisson data. We signify a Poisson data model with a mean response that is modelled using the square-root link function by setting `response = "poisson"` and `link = "sqrt"` in `FRK()`. Other non-Gaussian response distributions available in `FRK` are the binomial, negative-binomial, gamma, and inverse-Gaussian distributions. 

```r
## Simulate Poisson data using the previous example's data to construct a mean 
zdf$z <- rpois(m, lambda = zdf$z^2)

## Run FRK
S <- FRK(f = z ~ 1, list(zdf),                          
         response = "poisson",         # Poisson data model
         link = "sqrt")                # square-root link function
pred <- predict(S)                            

## Plotting
plot_list <- plot(S, pred$newdata)
ggarrange(plot_list$z, plot_list$p_mu, plot_list$interval90_mu, 
          nrow = 1, legend = "top")
             
```    
<!---
ggsave( 
  filename = "Poisson_data.png", device = "png", 
  width = 10, height = 4,
  path = "~/Dropbox/FRK/man/figures/"
)
--->

![(Left) Poisson data. (Centre) Prediction of the mean response. (Right) Prediction interval width of the mean response.](/man/figures/Poisson_data.png?raw=true)


### Spatio-temporal data

We now analyse spatio-temporal data, using the NOAA dataset.

```r
## Setup
data("NOAA_df_1990")
Tmax <- subset(NOAA_df_1990, month %in% 7 & year == 1993)
Tmax <- within(Tmax, {time = as.Date(paste(year,month,day,sep="-"))})
STObj <- stConstruct(x = Tmax, space = c("lon","lat"), time = "time", interval = TRUE)

## BAUs: spatial BAUs are 1x1 pixels, temporal BAUs are 1 day intervals
BAUs <- auto_BAUs(manifold = STplane(), 
                       cellsize = c(1, 1, 1),    
                       data=STObj, tunit = "days")
BAUs$fs <- 1 # scalar fine-scale variance matrix, implicit in previous examples

## Basis functions
G <- auto_basis(manifold = STplane(), data = STObj, nres = 2, tunit = "days")

## Run FRK
STObj$std <- 2 # fix the measurement error variance
S <- FRK(f = z ~ 1 + lat, data = list(STObj), 
         basis = G, BAUs = BAUs, est_error = FALSE, method = "TMB")
pred <- predict(S, percentiles = NULL)

## Plotting: include only some times via the argument subset_time
plot_list <- plot(S, pred$newdata, subset_time = c(1, 7, 13, 19, 25, 31)) 
ggarrange(plotlist = plot_list, nrow = 1, legend = "top") 
```

<!---
## Apply a labeller so the facet shows day x rather than just x
facet_names <- paste0("day ", unique(pred$newdata$t))
names(facet_names) <- unique(pred$newdata$t)
plot_list <- lapply(
  plot_list, 
  function(gg) gg + facet_wrap(~t, labeller = as_labeller(facet_names)))
  
ggsave( 
  filename = "ST_data.png", device = "png", 
  width = 12.5, height = 3.8,
  path = "~/Dropbox/FRK/man/figures/"
)
--->

![(Left) Prediction of spatio-temporal process. (Right) Prediction interval width.](/man/figures/ST_data.png?raw=true)


[//]: # (Currently `FRK` is not installing on OSX with `build_vignettes=TRUE` as it fails to find `texi2dvi`. Set `build_vignettes=FALSE` to ensure installation. Then download the `.Rnw` file in the `vignettes` folder and compile the pdf file separately in `RStudio` with `knitr`. )


Demonstrations
--------------

The package `FRK` is currently being used to generate spatio-temporal animations of fields observed by satellite data. [Here](https://www.youtube.com/watch?v=_kPa8VoeSdM) we show a daily prediction of CO2 using data from the NASA OCO-2 between September 2014 and June 2016.

[![alt tag](https://img.youtube.com/vi/ENx4CIZdoQk/0.jpg)](https://www.youtube.com/watch?v=ENx4CIZdoQk)

Acknowledgements
--------------

Thanks to [Michael Bertolacci](https://mbertolacci.github.io/) for designing the FRK hex logo!
