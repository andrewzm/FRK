# Plot predictions from FRK analysis

This function acts as a wrapper around
[`plot_spatial_or_ST`](plot_spatial_or_ST.md). It plots the fields of
the `Spatial*DataFrame` or `STFDF` object corresponding to prediction
and prediction uncertainty quantification. It also uses the `@data` slot
of `SRE` object to plot the training data set(s), and generates
informative, latex-style legend labels for each of the plots.

## Usage

``` r
plot(x, y, ...)

# S4 method for class 'SRE,list'
plot(x, y, ...)

# S4 method for class 'SRE,STFDF'
plot(x, y, ...)

# S4 method for class 'SRE,SpatialPointsDataFrame'
plot(x, y, ...)

# S4 method for class 'SRE,SpatialPixelsDataFrame'
plot(x, y, ...)

# S4 method for class 'SRE,SpatialPolygonsDataFrame'
plot(x, y, ...)
```

## Arguments

- x:

  object of class `SRE`

- y:

  the `Spatial*DataFrame` or `STFDF` object resulting from the call
  `predict(x)`. Keep in mind that
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns a `list`
  when `method` = "TMB"; the element `$newdata` contains the required
  `Spatial`/`ST` object. If the list itself is passed, you will receive
  the error: "x" and "y" lengths differ.

- ...:

  optional arguments passed on to
  [`plot_spatial_or_ST`](plot_spatial_or_ST.md)

## Value

A list of `ggplot` objects consisting of the observed data, predictions,
and standard errors. This list can then be supplied to, for example,
[`ggpubr::ggarrange()`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html).

## Examples

``` r
## See example in the help file for SRE
```
