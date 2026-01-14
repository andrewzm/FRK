# Automatic basis-function placement

Automatically generate a set of local basis functions in the domain, and
automatically prune in regions of sparse data.

## Usage

``` r
auto_basis(
  manifold = plane(),
  data,
  regular = 1,
  nres = 3,
  prune = 0,
  max_basis = NULL,
  subsamp = 10000,
  type = c("bisquare", "Gaussian", "exp", "Matern32"),
  isea3h_lo = 2,
  bndary = NULL,
  scale_aperture = ifelse(is(manifold, "sphere"), 1, 1.25),
  verbose = 0L,
  buffer = 0,
  tunit = NULL,
  ...
)
```

## Arguments

- manifold:

  object of class `manifold`, for example, `sphere` or `plane`

- data:

  object of class `SpatialPointsDataFrame` or `SpatialPolygonsDataFrame`
  containing the data on which basis-function placement is based, or a
  list of these; see details

- regular:

  an integer indicating the number of regularly-placed basis functions
  at the first resolution. In two dimensions, this dictates the smallest
  number of basis functions in a row or column at the coarsest
  resolution. If `regular=0`, an irregular grid is used, one that is
  based on the triangulation of the domain with increased mesh density
  in areas of high data density; see details

- nres:

  the number of basis-function resolutions to use

- prune:

  a threshold parameter that dictates when a basis function is
  considered irrelevent or unidentifiable, and thus removed; see details
  \[deprecated\]

- max_basis:

  maximum number of basis functions. This overrides the parameter `nres`

- subsamp:

  the maximum amount of data points to consider when carrying out
  basis-function placement: these data objects are randomly sampled from
  the full dataset. Keep this number fairly high (on the order of 10^5),
  otherwise fine-resolution basis functions may be spuriously removed

- type:

  the type of basis functions to use; see details

- isea3h_lo:

  if `manifold = sphere()`, this argument dictates which ISEA3H
  resolution is the coarsest one that should be used for the first
  resolution

- bndary:

  a `matrix` containing points containing the boundary. If
  `regular == 0` this can be used to define a boundary in which
  irregularly-spaced basis functions are placed

- scale_aperture:

  the aperture (in the case of the bisquare, but similar interpretation
  for other basis) width of the basis function is the minimum distance
  between all the basis function centroids multiplied by
  `scale_aperture`. Typically this ranges between 1 and 1.5 and is
  defaulted to 1 on the sphere and 1.25 on the other manifolds.

- verbose:

  a logical variable indicating whether to output a summary of the basis
  functions created or not

- buffer:

  a numeric between 0 and 0.5 indicating the size of the buffer of basis
  functions along the boundary. The buffer is added by computing the
  number of basis functions in each dimension, and increasing this
  number by a factor of `buffer`. A buffer may be needed when the prior
  distribution of the basis-function coefficients is formulated in terms
  of a precision matrix

- tunit:

  temporal unit, required when constructing a spatio-temporal basis.
  Should be the same as used for the BAUs. Can be "secs", "mins",
  "hours", "days", "years", etc.

- ...:

  unused

## Details

This function automatically places basis functions within the domain of
interest. If the domain is a plane or the real line, then the object
`data` is used to establish the domain boundary.

Let \\\phi(u)\\ denote the value of a basis function evaluated at \\u =
s - c\\, where \\s\\ is a spatial coordinate and \\c\\ is the
basis-function centroid. The argument `type` can be either “Gaussian”,
in which case

*φ(u) = exp(-\|u\|²/2σ²),*

“bisquare”, in which case

*φ(u) = (1 -(\|u\|/R)²)²,*

“exp”, in which case

*φ(u) = exp(-\|u\|/τ),*

or “Matern32”, in which case

*φ(u) = (1 + √3\|u\|/κ)exp(-√3\|u\|/κ),*

where the parameters \\\sigma, R, \tau\\ and \\\kappa\\ are `scale`
arguments.

If the manifold is the real line, the basis functions are placed
regularly inside the domain, and the number of basis functions at the
coarsest resolution is dictated by the integer parameter `regular` which
has to be greater than zero. On the real line, each subsequent
resolution has thrice as many basis functions. The scale of the basis
function is set based on the minimum distance between the centre
locations following placement. The scale is equal to the minimum
distance if the type of basis function is Gaussian, exponential, or
Matern32, and is equal to 1.5 times this value if the function is
bisquare.

If the manifold is a plane, and `regular > 0`, then basis functions are
placed regularly within the bounding box of `data`, with the smallest
number of basis functions in each row or column equal to three times the
value of `regular` in the coarsest resolution. Subsequent resolutions
have thrice the number of basis functions in each row or column. If
`regular = 0`, then the function
[`fmesher::fm_nonconvex_hull_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull_inla.html)
is used to construct a (non-convex) hull around the data. The buffer and
smoothness of the hull is determined by the parameter `convex`. Once the
domain boundary is found,
[`fmesher::fm_mesh_2d_inla()`](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
is used to construct a triangular mesh such that the node vertices
coincide with data locations, subject to some minimum and maximum
triangular-side-length constraints. The result is a mesh that is dense
in regions of high data density and not dense in regions of sparse data.
Even basis functions are irregularly placed, the scale is taken to be a
function of the minimum distance between basis function centres, as
detailed above. This may be changed in a future revision of the package.

If the manifold is the surface of a sphere, then basis functions are
placed on the centroids of the discrete global grid (DGG), with the
first basis resolution corresponding to the third resolution of the DGG
(ISEA3H resolution 2, which yields 92 basis functions globally). It is
not recommended to go above `nres == 3` (ISEA3H resolutions 2–4) for the
whole sphere; `nres=3` yields a total of 1176 basis functions. Up to
ISEA3H resolution 6 is available with `FRK`; for finer resolutions;
please install `dggrids` from `https://github.com/andrewzm/dggrids`
using `devtools`.

Basis functions that are not influenced by data points may hinder
convergence of the EM algorithm when `K_type = "unstructured"`, since
the associated hidden states are, by and large, unidentifiable. We hence
provide a means to automatically remove such basis functions through the
parameter `prune`. The final set only contains basis functions for which
the column sums in the associated matrix \\S\\ (which, recall, is the
value/average of the basis functions at/over the data points/polygons)
is greater than `prune`. If `prune == 0`, no basis functions are removed
from the original design.

## See also

[`remove_basis`](remove_basis.md) for removing basis functions and
[`show_basis`](show_basis.md) for visualising basis functions

## Examples

``` r
if (FALSE) { # \dontrun{
library(sp)
library(ggplot2)

## Create a synthetic dataset
set.seed(1)
d <- data.frame(lon = runif(n=1000,min = -179, max = 179),
                lat = runif(n=1000,min = -90, max = 90),
                z = rnorm(5000))
coordinates(d) <- ~lon + lat
slot(d, "proj4string") = CRS("+proj=longlat +ellps=sphere")

## Now create basis functions over sphere
G <- auto_basis(manifold = sphere(),data=d,
                nres = 2,prune=15,
                type = "bisquare",
                subsamp = 20000)

## Plot
show_basis(G,draw_world())
} # }
```
