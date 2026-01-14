# plane

Initialisation of a 2D plane.

## Usage

``` r
plane(measure = Euclid_dist(dim = 2L))
```

## Arguments

- measure:

  an object of class `measure`

## Details

A 2D plane is initialised using a `measure` object. By default, the
measure object (`measure`) is the Euclidean distance in 2 dimensions,
[Euclid_dist](distances.md).

## Examples

``` r
P <- plane()
print(type(P))
#> [1] "plane"
print(sp::dimensions(P))
#> [1] 2
```
