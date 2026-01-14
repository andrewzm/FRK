# Evaluate basis functions

Evaluate basis functions at points or average functions over polygons.

## Usage

``` r
eval_basis(basis, s)

# S4 method for class 'Basis,matrix'
eval_basis(basis, s)

# S4 method for class 'Basis,SpatialPointsDataFrame'
eval_basis(basis, s)

# S4 method for class 'Basis,SpatialPolygonsDataFrame'
eval_basis(basis, s)

# S4 method for class 'Basis,STIDF'
eval_basis(basis, s)

# S4 method for class 'TensorP_Basis,matrix'
eval_basis(basis, s)

# S4 method for class 'TensorP_Basis,STIDF'
eval_basis(basis, s)

# S4 method for class 'TensorP_Basis,STFDF'
eval_basis(basis, s)
```

## Arguments

- basis:

  object of class `Basis`

- s:

  object of class `matrix`, `SpatialPointsDataFrame` or
  `SpatialPolygonsDataFrame` containing the spatial locations/footprints

## Details

This function evaluates the basis functions at isolated points, or
averages the basis functions over polygons, for computing the matrix
\\S\\. The latter operation is carried out using Monte Carlo integration
with 1000 samples per polygon. When using space-time basis functions,
the object must contain a field `t` containing a numeric representation
of the time, for example, containing the number of seconds, hours, or
days since the first data point.

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions.

## Examples

``` r
library(sp)

### Create a synthetic dataset
set.seed(1)
d <- data.frame(lon = runif(n=500,min = -179, max = 179),
                lat = runif(n=500,min = -90, max = 90),
                z = rnorm(500))
coordinates(d) <- ~lon + lat
slot(d, "proj4string") = CRS("+proj=longlat")

### Now create basis functions on sphere
G <- auto_basis(manifold = sphere(),data=d,
                nres = 2,prune=15,
                type = "bisquare",
                subsamp = 20000)
#> NOTE: Zero process variability is implicitly enforced in regions where basis functions are pruned. Please use the option prune carefully: regions of data paucity are generally not reflective of regions of low process variability. Please set prune = 0 if unsure what to do.

### Now evaluate basis functions at origin
S <- eval_basis(G,matrix(c(0,0),1,2))
```
