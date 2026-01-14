# Compute distance

Compute distance using object of class `measure` or `manifold`.

## Usage

``` r
distance(d, x1, x2 = NULL)

# S4 method for class 'measure'
distance(d, x1, x2 = NULL)

# S4 method for class 'manifold'
distance(d, x1, x2 = NULL)
```

## Arguments

- d:

  object of class `measure` or `manifold`

- x1:

  first coordinate

- x2:

  second coordinate

## See also

[`real_line`](real_line.md), [`plane`](plane.md), [`sphere`](sphere.md),
[`STplane`](STplane.md) and [`STsphere`](STsphere.md) for constructing
manifolds, and [`distances`](distances.md) for the type of distances
available.

## Examples

``` r
distance(sphere(),matrix(0,1,2),matrix(10,1,2))
#>          [,1]
#> [1,] 1568.521
distance(plane(),matrix(0,1,2),matrix(10,1,2))
#>          [,1]
#> [1,] 14.14214
```
