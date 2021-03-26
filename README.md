[![Build Status](https://travis-ci.org/andrewzm/FRK.svg)](https://travis-ci.org/andrewzm/FRK)
[![codecov.io](http://codecov.io/github/andrewzm/FRK/coverage.svg?branch=master)](http://codecov.io/github/andrewzm/FRK?branch=master)

Fixed Rank Kriging 
================

<img src="./man/figures/FRK_logo.svg" width=100>



Installation 
------------

The package `FRK` is now at v1.0.0 and available on CRAN! To install, please type

```r
install.packages("FRK")
```

To install the most recent (development) version, first please install `INLA` from `http://www.r-inla.org/download`, then please load `devtools` and type

```r
install_github("andrewzm/FRK", dependencies = TRUE, build_vignettes = TRUE)
```

A document containing a description, details on the underlying maths and computations, as well as several examples, is available as a vignette. To load this vignette please type

```r
library("FRK")
vignette("FRK_intro")
```

A draft paper with more details, currently in press, is available from [here](https://arxiv.org/abs/1705.08105). A draft paper which details the approach to modelling non-Gaussian data is available **LINK TO FRK v2 PAPER**.

Description
------------

Package: FRK

Type: Package

Title: Fixed Rank Kriging

Version: 1.0.0

Date: 2020-03-25

Author: Andrew Zammit-Mangion, Matthew Sainsbury-Dale

Maintainer: Andrew Zammit-Mangion <andrewzm@gmail.com>

Description: Fixed Rank Kriging is a tool for spatial/spatio-temporal modelling and prediction with large datasets. The approach, discussed in Cressie and Johannesson (2008), decomposes the field, and hence the covariance function, using a fixed set of *r* basis functions, where *r* is typically much smaller than the number of data points (or polygons) *m*. The low-rank basis-function representation, combined with the use of sparse precision matrices, facilitates the modelling of big spatial/spatio-temporal data. The method naturally allows for non-stationary, anisotropic covariance functions, the use of observations with varying support, and any user-specified prediction regions. The package `FRK` employs a spatial generalised linear mixed model framework (Diggle et al., 1998) to cater for a range of non-Gaussian data models, including the Poisson, binomial, negative-binomial, gamma, and inverse-Gaussian distributions. In a non-Gaussian setting, `FRK` uses the package `TMB` (Kristensen et al., 2016) to make computationally efficient inferences through the use of Laplace approximations.  <!--The projected field is a key building block of the Spatial Random Effects (SRE) model, on which this package is based.--> `FRK` provides helper functions to model, fit, predict, and plot <!--using an SRE model -->with relative ease. 

* Cressie, N., & Johannesson, G. (2008). Fixed rank kriging for very large spatial data sets. Journal of the Royal Statistical Society: Series B, 70, 209–226.
* Diggle, P. J. & Tawn, J. A. & Moyeed, R. A. (1998). Model-based geostatistics. Journal of the Royal Statistical Society, 47:299–350
* Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M. (2016). TMB: Automatic differentiation and Laplace approximation. Journal of Statistical Software, 70:1–21.

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
plot_list <- plot(S, pred, zdf)
ggarrange(plotlist = plot_list, nrow = 1, legend = "top")

```

<!---
ggsave( 
  filename = "Gaussian_data.png", device = "png", 
  width = 10, height = 4,
  path = "~/Desktop/"
)
--->

![(Left) Gaussian data. (Centre) Predictions. (Right) Standard errors.](/man/figures/Gaussian_data.png?raw=true)

### Non-Gaussian data

Here we use analyse simulated Poisson data. We signify a Poisson data model with a mean response that is modelled using the square-root link function by setting `response = "poisson"` and `link = "square-root"` in `FRK()`. Other non-Gaussian data models available in `FRK` are the binomial, negative-binomial, gamma, and inverse-Gaussian. 

```r
## Setup
## Simulate Poisson data, using the data of the previous example to construct a mean 
RNGversion("3.6.0"); set.seed(1)                          
zdf$z <- rpois(m, lambda = zdf$z^2)

## Run FRK
S <- FRK(f = z ~ 1,                           # Formula to FRK
         list(zdf),                           # All datasets are supplied in list
         response = "poisson",                # Poisson data model
         link = "square-root")                # square-root link function
pred <- predict(S)                            # prediction stage


## Plotting
plot_list <- plot(S, pred, zdf)
ggarrange(plot_list$data, plot_list$p_mu, plot_list$interval_90_mu, 
          nrow = 1, legend = "top")

             
```    
<!---
ggsave( 
  filename = "Poisson_data.png", device = "png", 
  width = 10, height = 4,
  path = "~/Desktop/"
)
--->

![(Left) Poisson data. (Centre) Prediction of the mean response. (Right) Standard error of the mean response.](/man/figures/Poisson_data.png?raw=true)


```r
## Setup
data("NOAA_df_1990")
Tmax <- subset(NOAA_df_1990, month %in% 7 & year == 1993)
Tmax <- within(Tmax, {time = as.Date(paste(year,month,day,sep="-"))})
STObj <- stConstruct(x = Tmax, space = c("lon","lat"), time = "time", interval = TRUE)
STObj$std <- 2

## BAUs
## Choose spatial BAUs as 1x1 pixels, temporal BAUs as 4 day intervals
BAUs <- auto_BAUs(manifold = STplane(), 
                       cellsize = c(1, 1, 4),    
                       data=STObj, tunit = "days")
BAUs$fs <- 1 # scale fine-scale variance matrix, implicit in previous examples

## Basis functions
G <- auto_basis(manifold = STplane(), data = STObj, nres = 1, tunit = "days")

## Run FRK
S <- FRK(f = z ~ 1 + lat, data = list(STObj), 
         basis = G, BAUs = BAUs, est_error = FALSE, method = "TMB")
pred <- predict(S)

## Plotting
plot_list <- plot(S, pred) 
ggarrange(plot_list$p_mu, plot_list$interval_90_mu, align = "hv", legend = "top")


```

<!---
ggsave( 
  filename = "ST_data.png", device = "png", 
  width = 10, height = 4,
  path = "~/Desktop/"
)
--->



[//]: # (Currently `FRK` is not installing on OSX with `build_vignettes=TRUE` as it fails to find `texi2dvi`. Set `build_vignettes=FALSE` to ensure installation. Then download the `.Rnw` file in the `vignettes` folder and compile the pdf file separately in `RStudio` with `knitr`. )


Demonstrations
--------------

The package `FRK` is currently being used to generate spatio-temporal animations of fields observed by satellite data. [Here](https://www.youtube.com/watch?v=_kPa8VoeSdM) we show a daily prediction of CO2 using data from the NASA OCO-2 between September 2014 and June 2016.

[![alt tag](https://img.youtube.com/vi/ENx4CIZdoQk/0.jpg)](https://www.youtube.com/watch?v=ENx4CIZdoQk)

Acknowledgements
--------------

Thanks to [Michael Bertolacci](https://mbertolacci.github.io/) for designing the FRK hex logo!