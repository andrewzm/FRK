# Return the number of resolutions

Return the number of resolutions from a basis function object.

## Usage

``` r
nres(b)

# S4 method for class 'Basis'
nres(b)

# S4 method for class 'TensorP_Basis'
nres(b)

# S4 method for class 'SRE'
nres(b)
```

## Arguments

- b:

  object of class `Basis` or `SRE`

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions and [`show_basis`](show_basis.md) for visualising basis
functions.

## Examples

``` r
library(sp)
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
nres(G)
#> [1] 2
```
