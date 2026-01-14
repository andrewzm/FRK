# Show basis functions

Generic plotting function for visualising the basis functions.

## Usage

``` r
show_basis(basis, ...)

# S4 method for class 'Basis'
show_basis(basis, g = ggplot() + theme_bw() + xlab("") + ylab(""))

# S4 method for class 'TensorP_Basis'
show_basis(basis, g = ggplot())
```

## Arguments

- basis:

  object of class `Basis`

- ...:

  not in use

- g:

  object of class `gg` (a `ggplot` object) over which to overlay the
  basis functions (optional)

## Details

The function `show_basis` adapts its behaviour to the manifold being
used. With `real_line`, the 1D basis functions are plotted with colour
distinguishing between the different resolutions. With `plane`, only
local basis functions are supported (at present). Each basis function is
shown as a circle with diameter equal to the `scale` parameter of the
function. Linetype distinguishes the resolution. With `sphere`, the
centres of the basis functions are shown as circles, with larger sizes
corresponding to coarser resolutions. Space-time basis functions of
subclass `TensorP_Basis` are visualised by showing the spatial basis
functions and the temporal basis functions in two separate plots.

## See also

[`auto_basis`](auto_basis.md) for automatically constructing basis
functions.

## Examples

``` r
library(ggplot2)
library(sp)
data(meuse)
coordinates(meuse) = ~x+y # change into an sp object
G <- auto_basis(manifold = plane(),data=meuse,nres = 2,regular=2,prune=0.1,type = "bisquare")
#> NOTE: Zero process variability is implicitly enforced in regions where basis functions are pruned. Please use the option prune carefully: regions of data paucity are generally not reflective of regions of low process variability. Please set prune = 0 if unsure what to do.
if (FALSE) show_basis(G,ggplot()) + geom_point(data=data.frame(meuse),aes(x,y)) # \dontrun{}
```
