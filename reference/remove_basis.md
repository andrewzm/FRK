# Removes basis functions

Takes an object of class `Basis` and returns an object of class `Basis`
with selected basis functions removed

## Usage

``` r
remove_basis(Basis, rmidx)

# S4 method for class 'Basis,ANY'
remove_basis(Basis, rmidx)

# S4 method for class 'Basis,SpatialPolygons'
remove_basis(Basis, rmidx)
```

## Arguments

- Basis:

  object of class `Basis`

- rmidx:

  indices of basis functions to remove. Or a `SpatialPolygons` object;
  basis functions overlapping this `SpatialPolygons` object will be
  *retained*

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions and [`show_basis`](show_basis.md) for visualising basis
functions

## Examples

``` r
library(sp)
df <- data.frame(x = rnorm(10),
                 y = rnorm(10))
coordinates(df) <- ~x+y
G <- auto_basis(plane(),df,nres=1)
data.frame(G) # Print info on basis
#>          loc1       loc2    scale res
#> 1  -1.6219373 -1.2523559 2.055287   1
#> 2  -0.3118914 -1.2523559 2.055287   1
#> 3   0.9981544 -1.2523559 2.055287   1
#> 4  -1.6219373 -0.1562027 2.055287   1
#> 5  -0.3118914 -0.1562027 2.055287   1
#> 6   0.9981544 -0.1562027 2.055287   1
#> 7  -1.6219373  0.9399504 2.055287   1
#> 8  -0.3118914  0.9399504 2.055287   1
#> 9   0.9981544  0.9399504 2.055287   1
#> 10 -1.6219373  2.0361036 2.055287   1
#> 11 -0.3118914  2.0361036 2.055287   1
#> 12  0.9981544  2.0361036 2.055287   1

## Removing basis functions by index
G_subset <- remove_basis(G, 1:(nbasis(G)-1))
data.frame(G_subset)
#>         loc1     loc2    scale res
#> 12 0.9981544 2.036104 2.055287   1

## Removing basis functions using SpatialPolygons
x <- 1
poly <- Polygon(rbind(c(-x, -x), c(-x, x), c(x, x), c(x, -x), c(-x, -x)))
polys <- Polygons(list(poly), "1")
spatpolys <- SpatialPolygons(list(polys))
G_subset <- remove_basis(G, spatpolys)
data.frame(G_subset)
#>         loc1       loc2    scale res
#> 5 -0.3118914 -0.1562027 2.055287   1
#> 6  0.9981544 -0.1562027 2.055287   1
#> 8 -0.3118914  0.9399504 2.055287   1
#> 9  0.9981544  0.9399504 2.055287   1
```
