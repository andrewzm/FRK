# sphere

Initialisation of the 2-sphere, S2.

## Usage

``` r
sphere(radius = 6371)
```

## Arguments

- radius:

  radius of sphere

## Details

The 2D surface of a sphere is initialised using a `radius` parameter.
The default value of the radius `R` is `R`=6371 km, Earth's radius,
while the measure used to compute distances on the sphere is the
great-circle distance on a sphere of radius `R`.

## Examples

``` r
S <- sphere()
print(sp::dimensions(S))
#> [1] 2
```
