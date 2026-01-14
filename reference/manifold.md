# Retrieve manifold

Retrieve manifold from `FRK` object.

## Usage

``` r
manifold(.Object)

# S4 method for class 'Basis'
manifold(.Object)

# S4 method for class 'TensorP_Basis'
manifold(.Object)
```

## Arguments

- .Object:

  `FRK` object

## See also

[`real_line`](real_line.md), [`plane`](plane.md), [`sphere`](sphere.md),
[`STplane`](STplane.md) and [`STsphere`](STsphere.md) for constructing
manifolds.

## Examples

``` r
G <-  local_basis(manifold = plane(),
                   loc=matrix(0,1,2),
                   scale=0.2,
                   type="bisquare")
manifold(G)
#> Type of manifold: plane 
#> Dimension of manifold: 2 
#> Distance function:
#>  function (x1, x2)  distR(x1, x2) 
```
