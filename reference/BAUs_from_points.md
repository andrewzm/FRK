# Creates pixels around points

Takes a SpatialPointsDataFrame and converts it into
SpatialPolygonsDataFrame by constructing a tiny (within machine
tolerance) BAU around each SpatialPoint.

## Usage

``` r
BAUs_from_points(obj, offset = 1e-10)

# S4 method for class 'SpatialPoints'
BAUs_from_points(obj, offset = 1e-10)

# S4 method for class 'ST'
BAUs_from_points(obj, offset = 1e-10)
```

## Arguments

- obj:

  object of class `SpatialPointsDataFrame`

- offset:

  edge size of the mini-BAU (default 1e-10)

## Details

This function allows users to mimic standard geospatial analysis where
BAUs are not used. Since `FRK` is built on the concept of a BAU, this
function constructs tiny BAUs around the observation and prediction
locations that can be subsequently passed on to the functions `SRE` and
`FRK`. With `BAUs_from_points`, the user supplies both the data and
prediction locations accompanied with covariates.

## See also

[`auto_BAUs`](auto_BAUs.md) for automatically constructing generic BAUs.

## Examples

``` r
library(sp)
opts_FRK$set("parallel",0L)
df <- data.frame(x = rnorm(10),
                 y = rnorm(10))
coordinates(df) <- ~x+y
BAUs <- BAUs_from_points(df)
```
