# Pre-configured distances

Useful objects of class `distance` included in package.

## Usage

``` r
measure(dist, dim)

Euclid_dist(dim = 2L)

gc_dist(R = NULL)

gc_dist_time(R = NULL)
```

## Arguments

- dist:

  a function taking two arguments `x1,x2`

- dim:

  the dimension of the manifold (e.g., 2 for a plane)

- R:

  great-circle radius

## Details

Initialises an object of class `measure` which contains a function
`dist` used for computing the distance between two points. Currently the
Euclidean distance and the great-circle distance are included with
`FRK`.

## Examples

``` r
M1 <- measure(distR,2)
D <- distance(M1,matrix(rnorm(10),5,2))
```
