# Construct a set of local basis functions

Construct a set of local basis functions based on pre-specified location
and scale parameters.

## Usage

``` r
local_basis(
  manifold = sphere(),
  loc = matrix(c(1, 0), nrow = 1),
  scale = 1,
  type = c("bisquare", "Gaussian", "exp", "Matern32"),
  res = 1,
  regular = FALSE
)

radial_basis(
  manifold = sphere(),
  loc = matrix(c(1, 0), nrow = 1),
  scale = 1,
  type = c("bisquare", "Gaussian", "exp", "Matern32")
)
```

## Arguments

- manifold:

  object of class `manifold`, for example, `sphere`

- loc:

  a matrix of size `n` by `dimensions(manifold)` indicating centres of
  basis functions

- scale:

  vector of length `n` containing the scale parameters of the basis
  functions; see details

- type:

  either `"bisquare"`, `"Gaussian"`, `"exp"`, or `"Matern32"`

- res:

  vector of length `n` containing the resolutions of the basis functions

- regular:

  logical indicating if the basis functions (of each resolution) are in
  a regular grid

## Details

This functions lays out local basis functions in a domain of interest
based on pre-specified location and scale parameters. If `type` is
“bisquare”, then \$\$\phi(u) = \left(1- \left(\frac{\\ u
\\}{R}\right)^2\right)^2 I(\\u\\ \< R),\$\$ and `scale` is given by
\\R\\, the range of support of the bisquare function. If `type` is
“Gaussian”, then \$\$\phi(u) = \exp\left(-\frac{\\u
\\^2}{2\sigma^2}\right),\$\$ and `scale` is given by \\\sigma\\, the
standard deviation. If `type` is “exp”, then \$\$\phi(u) =
\exp\left(-\frac{\\u\\}{ \tau}\right),\$\$ and `scale` is given by
\\\tau\\, the e-folding length. If `type` is “Matern32”, then
\$\$\phi(u) = \left(1 +
\frac{\sqrt{3}\\u\\}{\kappa}\right)\exp\left(-\frac{\sqrt{3}\\ u
\\}{\kappa}\right),\$\$ and `scale` is given by \\\kappa\\, the
function's scale.

## See also

[`auto_basis`](auto_basis.md) for constructing basis functions
automatically, and [`show_basis`](show_basis.md) for visualising basis
functions.

## Examples

``` r
library(ggplot2)
G <-  local_basis(manifold = real_line(),
                   loc=matrix(1:10,10,1),
                   scale=rep(2,10),
                   type="bisquare")
if (FALSE) show_basis(G) # \dontrun{}
```
