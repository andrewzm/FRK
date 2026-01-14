# Space-time sphere

Initialisation of a 2-sphere (S2) with a temporal dimension

## Usage

``` r
STsphere(radius = 6371)
```

## Arguments

- radius:

  radius of sphere

## Details

As with the spatial-only sphere, the sphere surface is initialised using
a `radius` parameter. The default value of the radius `R` is `R`=6371,
which is the Earth's radius in km, while the measure used to compute
distances on the sphere is the great-circle distance on a sphere of
radius `R`. By default Euclidean geometry is used to factor in the time
component, so that dist((s1,t1),(s2,t2)) = sqrt(gc_dist(s1,s2)^2 + (t1 -
t2)^2). Frequently this distance can be used since separate correlation
length scales for space and time are estimated in the EM algorithm (that
effectively scale space and time separately).

## Examples

``` r
S <- STsphere()
print(sp::dimensions(S))
#> [1] 3
```
