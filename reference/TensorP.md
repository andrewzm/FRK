# Tensor product of basis functions

Constructs a new set of basis functions by finding the tensor product of
two sets of basis functions.

## Usage

``` r
TensorP(Basis1, Basis2)

# S4 method for class 'Basis,Basis'
TensorP(Basis1, Basis2)
```

## Arguments

- Basis1:

  first set of basis functions

- Basis2:

  second set of basis functions

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions and [`show_basis`](show_basis.md) for visualising basis
functions.

## Examples

``` r
library(spacetime)
library(sp)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
sim_data <- data.frame(lon = runif(20,-180,180),
                       lat = runif(20,-90,90),
                       t = 1:20,
                       z = rnorm(20),
                       std = 0.1)
time <- as.POSIXct("2003-05-01",tz="") + 3600*24*(sim_data$t-1)
space <- sim_data[,c("lon","lat")]
coordinates(space) = ~lon+lat # change into an sp object
slot(space, "proj4string") = CRS("+proj=longlat +ellps=sphere")
STobj <- STIDF(space,time,data=sim_data)
G_spatial <- auto_basis(manifold = sphere(),
                        data=as(STobj,"Spatial"),
                        nres = 1,
                        type = "bisquare",
                        subsamp = 20000)
G_temporal <- local_basis(manifold=real_line(),loc = matrix(c(1,3)),scale = rep(1,2))
G <- TensorP(G_spatial,G_temporal)
# show_basis(G_spatial)
# show_basis(G_temporal)
```
