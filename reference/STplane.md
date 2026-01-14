# plane in space-time

Initialisation of a 2D plane with a temporal dimension.

## Usage

``` r
STplane(measure = Euclid_dist(dim = 3L))
```

## Arguments

- measure:

  an object of class `measure`

## Details

A 2D plane with a time component added is initialised using a `measure`
object. By default, the measure object (`measure`) is the Euclidean
distance in 3 dimensions, [Euclid_dist](distances.md).

## Examples

``` r
P <- STplane()
print(type(P))
#> [1] "STplane"
print(sp::dimensions(P))
#> [1] 3
```
