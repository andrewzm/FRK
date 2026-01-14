# Spatial Random Effects class

This is the central class definition of the `FRK` package, containing
the model and all other information required for estimation and
prediction.

## Details

The spatial random effects (SRE) model is the model employed in Fixed
Rank Kriging, and the `SRE` object contains all information required for
estimation and prediction from spatial data. Object slots contain both
other objects (for example, an object of class `Basis`) and matrices
derived from these objects (for example, the matrix \\S\\) in order to
facilitate computations.

## Slots

- `f`:

  formula used to define the SRE object. All covariates employed need to
  be specified in the object `BAUs`

- `data`:

  the original data from which the model's parameters are estimated

- `basis`:

  object of class `Basis` used to construct the matrix \\S\\

- `BAUs`:

  object of class `SpatialPolygonsDataFrame`, `SpatialPixelsDataFrame`
  of `STFDF` that contains the Basic Areal Units (BAUs) that are used to
  both (i) project the data onto a common discretisation if they are
  point-referenced and (ii) provide a BAU-to-data relationship if the
  data has a spatial footprint

- `S`:

  matrix constructed by evaluating the basis functions at all the data
  locations (of class `Matrix`)

- `S0`:

  matrix constructed by evaluating the basis functions at all BAUs (of
  class `Matrix`)

- `D_basis`:

  list of distance-matrices of class `Matrix`, one for each
  basis-function resolution

- `Ve`:

  measurement-error variance-covariance matrix (typically diagonal and
  of class `Matrix`)

- `Vfs`:

  fine-scale variance-covariance matrix at the data locations (typically
  diagonal and of class `Matrix`) up to a constant of proportionality
  estimated using the EM algorithm

- `Vfs_BAUs`:

  fine-scale variance-covariance matrix at the BAU centroids (typically
  diagonal and of class `Matrix`) up to a constant of proportionality
  estimated using the EM algorithm

- `Qfs_BAUs`:

  fine-scale precision matrix at the BAU centroids (typically diagonal
  and of class `Matrix`) up to a constant of proportionality estimated
  using the EM algorithm

- `Z`:

  vector of observations (of class `Matrix`)

- `Cmat`:

  incidence matrix mapping the observations to the BAUs

- `X`:

  design matrix of covariates at all the data locations

- `G`:

  list of objects of class Matrix containing the design matrices for
  random effects at all the data locations

- `G0`:

  list of objects of class Matrix containing the design matrices for
  random effects at all BAUs

- `K_type`:

  type of prior covariance matrix of random effects. Can be
  "block-exponential" (correlation between effects decays as a function
  of distance between the basis-function centroids), "unstructured" (all
  elements in `K` are unknown and need to be estimated), or "neighbour"
  (a sparse precision matrix is used, whereby only neighbouring basis
  functions have non-zero precision matrix elements).

- `mu_eta`:

  updated expectation of the basis-function random effects (estimated)

- `mu_gamma`:

  updated expectation of the random effects (estimated)

- `S_eta`:

  updated covariance matrix of random effects (estimated)

- `Q_eta`:

  updated precision matrix of random effects (estimated)

- `Khat`:

  prior covariance matrix of random effects (estimated)

- `Khat_inv`:

  prior precision matrix of random effects (estimated)

- `alphahat`:

  fixed-effect regression coefficients (estimated)

- `sigma2fshat`:

  fine-scale variation scaling (estimated)

- `sigma2gamma`:

  random-effect variance parameters (estimated)

- `fs_model`:

  type of fine-scale variation (independent or CAR-based). Currently
  only "ind" is permitted

- `info_fit`:

  information on fitting (convergence etc.)

- `response`:

  A character string indicating the assumed distribution of the response
  variable

- `link`:

  A character string indicating the desired link function. Can be "log",
  "identity", "logit", "probit", "cloglog", "reciprocal", or
  "reciprocal-squared". Note that only sensible link-function and
  response-distribution combinations are permitted.

- `mu_xi`:

  updated expectation of the fine-scale random effects at all BAUs
  (estimated)

- `Q_posterior`:

  updated joint precision matrix of the basis function random effects
  and observed fine-scale random effects (estimated)

- `log_likelihood`:

  the log likelihood of the fitted model

- `method`:

  the fitting procedure used to fit the SRE model

- `phi`:

  the estimated dispersion parameter (assumed constant throughout the
  spatial domain)

- `k_Z`:

  vector of known size parameters at the observation support level (only
  applicable to binomial and negative-binomial response distributions)

- `k_BAU`:

  vector of known size parameters at the observed BAUs (only applicable
  to binomial and negative-binomial response distributions)

- `include_fs`:

  flag indicating whether the fine-scale variation should be included in
  the model

- `include_gamma`:

  flag indicating whether there are gamma random effects in the model

- `normalise_wts`:

  if `TRUE`, the rows of the incidence matrices \\C_Z\\ and \\C_P\\ are
  normalised to sum to 1, so that the mapping represents a weighted
  average; if false, no normalisation of the weights occurs (i.e., the
  mapping corresponds to a weighted sum)

- `fs_by_spatial_BAU`:

  if `TRUE`, then each BAU is associated with its own fine-scale
  variance parameter

- `obsidx`:

  indices of observed BAUs

- `simple_kriging_fixed`:

  logical indicating whether one wishes to commit to simple kriging at
  the fitting stage: If `TRUE`, model fitting is faster, but the option
  to conduct universal kriging at the prediction stage is removed

## References

Zammit-Mangion, A. and Cressie, N. (2017). FRK: An R package for spatial
and spatio-temporal prediction with large datasets. Journal of
Statistical Software, 98(4), 1-48. doi:10.18637/jss.v098.i04.

## See also

[`SRE`](SRE.md) for details on how to construct and fit SRE models.
