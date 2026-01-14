# Type of manifold

Retrieve slot `type` from object

## Usage

``` r
type(.Object)

# S4 method for class 'manifold'
type(.Object)
```

## Arguments

- .Object:

  object of class `Basis` or `manifold`

## See also

[`real_line`](real_line.md), [`plane`](plane.md), [`sphere`](sphere.md),
[`STplane`](STplane.md) and [`STsphere`](STsphere.md) for constructing
manifolds.

## Examples

``` r
S <- sphere()
print(type(S))
#> [1] "surface of sphere"
```
