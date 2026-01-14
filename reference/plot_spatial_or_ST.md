# Plot a Spatial\*DataFrame or STFDF object

Takes an object of class `Spatial*DataFrame` or `STFDF`, and plots
requested data columns using `ggplot2`

## Usage

``` r
plot_spatial_or_ST(
  newdata,
  column_names,
  map_layer = NULL,
  subset_time = NULL,
  palette = "Spectral",
  plot_over_world = FALSE,
  labels_from_coordnames = TRUE,
  ...
)

# S4 method for class 'STFDF'
plot_spatial_or_ST(
  newdata,
  column_names,
  map_layer = NULL,
  subset_time = NULL,
  palette = "Spectral",
  plot_over_world = FALSE,
  labels_from_coordnames = TRUE,
  ...
)

# S4 method for class 'SpatialPointsDataFrame'
plot_spatial_or_ST(
  newdata,
  column_names,
  map_layer = NULL,
  subset_time = NULL,
  palette = "Spectral",
  plot_over_world = FALSE,
  labels_from_coordnames = TRUE,
  ...
)

# S4 method for class 'SpatialPixelsDataFrame'
plot_spatial_or_ST(
  newdata,
  column_names,
  map_layer = NULL,
  subset_time = NULL,
  palette = "Spectral",
  plot_over_world = FALSE,
  labels_from_coordnames = TRUE,
  ...
)

# S4 method for class 'SpatialPolygonsDataFrame'
plot_spatial_or_ST(
  newdata,
  column_names,
  map_layer = NULL,
  subset_time = NULL,
  palette = "Spectral",
  plot_over_world = FALSE,
  labels_from_coordnames = TRUE,
  ...
)
```

## Arguments

- newdata:

  an object of class `Spatial*DataFrame` or `STFDF`

- column_names:

  a vector of strings indicating the columns of the data to plot

- map_layer:

  (optional) a `ggplot` layer or object to add below the plotted layer,
  often a map

- subset_time:

  (optional) a vector of times to be included; applicable only for
  `STFDF` objects

- palette:

  the palette supplied to the argument `palette` of
  `scale_*_distiller()`. Alternatively, if `palette` = "nasa", a vibrant
  colour palette is created using `scale_*_gradientn()`

- plot_over_world:

  logical; if `TRUE`, `coord_map("mollweide")` and
  [`draw_world`](draw_world.md) are used to plot over the world

- labels_from_coordnames:

  logical; if `TRUE`, the coordinate names of `newdata` (i.e.,
  `coordnames(newdata)`) are used as the horizontal- and vertical-axis
  labels. Otherwise, generic names, s_1 and s_2, are used

- ...:

  optional arguments passed on to whatever geom is appropriate for the
  `Spatial*DataFrame` or `STFDF` object (`geom_point`, `geom_tile`,
  `geom_raster`, or `geom_polygon`)

## Value

A list of `ggplot` objects corresponding to the provided `column_names`.
This list can then be supplied to, for example,
[`ggpubr::ggarrange()`](https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html).

## See also

[`plot`](plot.md)

## Examples

``` r
## See example in the help file for FRK
```
