# Combine basis functions

Takes a list of objects of class `Basis` and returns a single object of
class `Basis`.

## Usage

``` r
combine_basis(Basis_list)

# S4 method for class 'list'
combine_basis(Basis_list)
```

## Arguments

- Basis_list:

  a list of objects of class `Basis`. Each element of the list is
  assumed to represent a single resolution of basis functions

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions and [`show_basis`](show_basis.md) for visualising basis
functions

## Examples

``` r
## Construct two resolutions of basis functions using local_basis() 
Basis1 <- local_basis(manifold = real_line(), 
                     loc = matrix(seq(0, 1, length.out = 3), ncol = 1), 
                     scale = rep(0.4, 3))

Basis2 <- local_basis(manifold = real_line(), 
                     loc = matrix(seq(0, 1, length.out = 6), ncol = 1), 
                     scale = rep(0.2, 6))

## Combine basis-function resolutions into a single Basis object
combine_basis(list(Basis1, Basis2)) 
#> Number of basis functions: 9 
#> Number of resolutions: 2 
#> Regular: 0 
#> Type of manifold: real_line 
#> Dimension of manifold: 1 
#> First basis function:
#>  function (s)  {     stopifnot(ncol(s) == dimensions(manifold))     y <- distance(manifold, s, c)     (1 - (y/R)^2)^2 * (y < R) } 
```
