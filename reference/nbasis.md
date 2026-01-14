# Number of basis functions

Retrieve the number of basis functions from `Basis` or `SRE` object.

## Usage

``` r
nbasis(.Object)

# S4 method for class 'Basis_obj'
nbasis(.Object)

# S4 method for class 'SRE'
nbasis(.Object)
```

## Arguments

- .Object:

  object of class `Basis` or `SRE`

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions.

## Examples

``` r
library(sp)
data(meuse)
coordinates(meuse) = ~x+y # change into an sp object
G <- auto_basis(manifold = plane(),
                data=meuse,
                nres = 2,
                regular=1,
                type = "Gaussian")
print(nbasis(G))
#> [1] 129
```
