# Automatic BAU generation

This function calls the generic function `auto_BAU` (not exported) after
a series of checks and is the easiest way to generate a set of Basic
Areal Units (BAUs) on the manifold being used; see details.

## Usage

``` r
auto_BAUs(
  manifold,
  type = NULL,
  cellsize = NULL,
  isea3h_res = NULL,
  data = NULL,
  nonconvex_hull = TRUE,
  convex = -0.05,
  tunit = NULL,
  xlims = NULL,
  ylims = NULL,
  spatial_BAUs = NULL,
  ...
)
```

## Arguments

- manifold:

  object of class `manifold`

- type:

  either “grid” or “hex”, indicating whether gridded or hexagonal BAUs
  should be used. If `type` is unspecified, “hex” will be used if we are
  on the sphere, and “grid” will used otherwise

- cellsize:

  denotes size of gridcell when `type` = “grid”. Needs to be of length 1
  (square-grid case) or a vector of length `dimensions(manifold)`
  (rectangular-grid case)

- isea3h_res:

  resolution number of the isea3h DGGRID cells for when type is “hex”
  and manifold is the surface of a `sphere`

- data:

  object of class `SpatialPointsDataFrame`, `SpatialPolygonsDataFrame`,
  `STIDF`, or `STFDF`. Provision of `data` implies that the domain is
  bounded, and is thus necessary when the manifold is a
  `real_line, plane`, or `STplane`, but is not necessary when the
  manifold is the surface of a `sphere`

- nonconvex_hull:

  flag indicating whether to use `fmesher` to generate a non-convex
  hull. Otherwise a convex hull is used

- convex:

  convex parameter used for smoothing an extended boundary when working
  on a bounded domain (that is, when the object `data` is supplied); see
  details

- tunit:

  temporal unit when requiring space-time BAUs. Can be "secs", "mins",
  "hours", etc.

- xlims:

  limits of the horizontal axis (overrides automatic selection)

- ylims:

  limits of the vertical axis (overrides automatic selection)

- spatial_BAUs:

  object of class `SpatialPolygonsDataFrame` or `SpatialPixelsDataFrame`
  representing the spatial BAUs to be used in a spatio-temporal setting
  (if left `NULL`, the spatial BAUs are constructed automatically using
  the data)

- ...:

  currently unused

## Details

`auto_BAUs` constructs a set of Basic Areal Units (BAUs) used both for
data pre-processing and for prediction. As such, the BAUs need to be of
sufficienly fine resolution so that inferences are not affected due to
binning.

Two types of BAUs are supported by `FRK`: “hex” (hexagonal) and “grid”
(rectangular). In order to have a “grid” set of BAUs, the user should
specify a cellsize of length one, or of length equal to the dimensions
of the manifold, that is, of length 1 for `real_line` and of length 2
for the surface of a `sphere` and `plane`. When a “hex” set of BAUs is
desired, the first element of `cellsize` is used to determine the side
length by dividing this value by approximately 2. The argument `type` is
ignored with `real_line` and “hex” is not available for this manifold.

If the object `data` is provided, then automatic domain selection may be
carried out by employing the `fmesher` function
`fm_nonconvex_hull_inla`, which finds a (non-convex) hull surrounding
the data points (or centroids of the data polygons). This domain is
extended and smoothed using the parameter `convex`. The parameter
`convex` should be negative, and a larger absolute value for `convex`
results in a larger domain with smoother boundaries.

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions.

## Examples

``` r
## First a 1D example
library(sp)
set.seed(1)
data <- data.frame(x = runif(10)*10, y = 0, z= runif(10)*10)
coordinates(data) <- ~x+y
Grid1D_df <- auto_BAUs(manifold = real_line(),
                       cellsize = 1,
                       data=data)
if (FALSE) spplot(Grid1D_df) # \dontrun{}

## Now a 2D example
data(meuse)
coordinates(meuse) = ~x+y # change into an sp object

## Grid BAUs
GridPols_df <- auto_BAUs(manifold = plane(),
                         cellsize = 200,
                         type = "grid",
                         data = meuse,
                         nonconvex_hull = 0)
if (FALSE) plot(GridPols_df) # \dontrun{}

## Hex BAUs
HexPols_df <- auto_BAUs(manifold = plane(),
                        cellsize = 200,
                        type = "hex",
                        data = meuse,
                        nonconvex_hull = 0)
if (FALSE) plot(HexPols_df) # \dontrun{}
```
