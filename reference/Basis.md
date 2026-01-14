# Generic basis-function constructor

This function is meant to be used for manual construction of arbitrary
basis functions. For \`local' basis functions, please use the function
[`local_basis`](local_basis.md) instead.

## Usage

``` r
Basis(manifold, n, fn, pars, df, regular = FALSE)
```

## Arguments

- manifold:

  object of class `manifold`, for example, `sphere`

- n:

  number of basis functions (should be an integer)

- fn:

  a list of functions, one for each basis function. Each function should
  be encapsulated within an environment in which the manifold and any
  other parameters required to evaluate the function are defined. The
  function itself takes a single input `s` which can be of class
  `numeric`, `matrix`, or `Matrix`, and returns a vector which contains
  the basis function evaluations at `s`.

- pars:

  A list containing a list of parameters for each function. For local
  basis functions these would correspond to location and scale
  parameters.

- df:

  A data frame containing one row per basis function, typically for
  providing informative summaries.

- regular:

  logical indicating if the basis functions (of each resolution) are in
  a regular grid

## Details

This constructor checks that all parameters are valid before
constructing the basis functions. The requirement that every function is
encapsulated is tedious, but necessary for FRK to work with a large
range of basis functions in the future. Please see the example below
which exemplifies the process of constructing linear basis functions
from scratch using this function.

## See also

[`auto_basis`](auto_basis.md) for constructing basis functions
automatically, [`local_basis`](local_basis.md) for constructing \`local'
basis functions, and [`show_basis`](show_basis.md) for visualising basis
functions.

## Examples

``` r
## Construct two linear basis functions on [0, 1]
manifold <- real_line()
n <- 2
lin_basis_fn <- function(manifold, grad, intercept) {
   function(s) grad*s + intercept
}
pars <- list(list(grad = 1, intercept = 0),
             list(grad = -1, intercept = 1))
fn <- list(lin_basis_fn(manifold, 1, 0),
           lin_basis_fn(manifold, -1, 1))
df <- data.frame(n = 1:2, grad = c(1, -1), m = c(1, -1))
G <- Basis(manifold = manifold, n = n, fn = fn, pars = pars, df = df)
if (FALSE) { # \dontrun{
eval_basis(G, s = matrix(seq(0,1, by = 0.1), 11, 1))} # }
```
