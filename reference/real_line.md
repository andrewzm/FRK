# real line

Initialisation of the real-line (1D) manifold.

## Usage

``` r
real_line(measure = Euclid_dist(dim = 1L))
```

## Arguments

- measure:

  an object of class `measure`

## Details

A real line is initialised using a `measure` object. By default, the
measure object (`measure`) describes the distance between two points as
the absolute difference between the two coordinates.

## Examples

``` r
R <- real_line()
print(type(R))
#> [1] "real_line"
print(sp::dimensions(R))
#> [1] 1
```
