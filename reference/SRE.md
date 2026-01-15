# Construct SRE object, fit and predict

The Spatial Random Effects (SRE) model is the central object in FRK. The
function `FRK()` provides a wrapper for the construction and estimation
of the SRE object from data, using the functions `SRE()` (the object
constructor) and `SRE.fit()` (for fitting it to the data). Please see
[`SRE-class`](SRE-class.md) for more details on the SRE object's
properties and methods.

## Usage

``` r
FRK(
  f,
  data,
  basis = NULL,
  BAUs = NULL,
  est_error = TRUE,
  average_in_BAU = TRUE,
  sum_variables = NULL,
  normalise_wts = TRUE,
  fs_model = "ind",
  vgm_model = NULL,
  K_type = c("block-exponential", "precision", "unstructured"),
  n_EM = 100,
  tol = 0.01,
  method = c("EM", "TMB"),
  lambda = 0,
  print_lik = FALSE,
  response = c("gaussian", "poisson", "gamma", "inverse-gaussian", "negative-binomial",
    "binomial"),
  link = c("identity", "log", "sqrt", "logit", "probit", "cloglog", "inverse",
    "inverse-squared"),
  optimiser = nlminb,
  fs_by_spatial_BAU = FALSE,
  known_sigma2fs = NULL,
  taper = NULL,
  simple_kriging_fixed = FALSE,
  ...
)

SRE(
  f,
  data,
  basis,
  BAUs,
  est_error = TRUE,
  average_in_BAU = TRUE,
  sum_variables = NULL,
  normalise_wts = TRUE,
  fs_model = "ind",
  vgm_model = NULL,
  K_type = c("block-exponential", "precision", "unstructured"),
  normalise_basis = TRUE,
  response = c("gaussian", "poisson", "gamma", "inverse-gaussian", "negative-binomial",
    "binomial"),
  link = c("identity", "log", "sqrt", "logit", "probit", "cloglog", "inverse",
    "inverse-squared"),
  include_fs = TRUE,
  fs_by_spatial_BAU = FALSE,
  ...
)

SRE.fit(
  object,
  n_EM = 100L,
  tol = 0.01,
  method = c("EM", "TMB"),
  lambda = 0,
  print_lik = FALSE,
  optimiser = nlminb,
  known_sigma2fs = NULL,
  taper = NULL,
  simple_kriging_fixed = FALSE,
  ...
)

# S4 method for class 'SRE'
predict(
  object,
  newdata = NULL,
  obs_fs = FALSE,
  pred_time = NULL,
  covariances = FALSE,
  nsim = 400,
  type = "mean",
  k = NULL,
  percentiles = c(5, 95),
  kriging = "simple",
  new_BAU_data = NULL
)

# S4 method for class 'SRE'
logLik(object)

# S4 method for class 'SRE'
nobs(object, ...)

# S4 method for class 'SRE'
coef(object, ...)

# S4 method for class 'SRE'
coef_uncertainty(
  object,
  percentiles = c(5, 95),
  nsim = 400,
  random_effects = FALSE
)

simulate(object, newdata = NULL, nsim = 400, conditional_fs = FALSE, ...)

# S4 method for class 'SRE'
fitted(object, ...)

# S4 method for class 'SRE'
residuals(object, type = "pearson")

# S4 method for class 'SRE'
AIC(object, k = 2)

# S4 method for class 'SRE'
BIC(object)
```

## Arguments

- f:

  `R` formula relating the dependent variable (or transformations
  thereof) to covariates

- data:

  list of objects of class `SpatialPointsDataFrame`,
  `SpatialPolygonsDataFrame`, `STIDF`, or `STFDF`. If using space-time
  objects, the data frame must have another field, `t`, containing the
  time index of the data point

- basis:

  object of class `Basis` (or `TensorP_Basis`)

- BAUs:

  object of class `SpatialPolygonsDataFrame`, `SpatialPixelsDataFrame`,
  `STIDF`, or `STFDF`. The object's data frame must contain covariate
  information as well as a field `fs` describing the fine-scale
  variation up to a constant of proportionality. If the function `FRK()`
  is used directly, then BAUs are created automatically, but only
  coordinates can then be used as covariates

- est_error:

  (applicable only if `response` = "gaussian") flag indicating whether
  the measurement-error variance should be estimated from variogram
  techniques. If this is set to 0, then `data` must contain a field
  `std`. Measurement-error estimation is currently not implemented for
  spatio-temporal datasets

- average_in_BAU:

  if `TRUE`, then multiple data points falling in the same BAU are
  averaged; the measurement error of the averaged data point is taken as
  the average of the individual measurement errors

- sum_variables:

  if `average_in_BAU == TRUE`, the string `sum_variables` indicates
  which data variables (can be observations or covariates) are to be
  summed rather than averaged

- normalise_wts:

  if `TRUE`, the rows of the incidence matrices ***C**_(Z)* and
  ***C**_(P)* are normalised to sum to 1, so that the mapping represents
  a weighted average; if false, no normalisation of the weights occurs
  (i.e., the mapping corresponds to a weighted sum)

- fs_model:

  if "ind" then the fine-scale variation is independent at the BAU
  level. Only the independent model is allowed for now, future
  implementation will include CAR/ICAR (in development)

- vgm_model:

  (applicable only if `response` = "gaussian") an object of class
  `variogramModel` from the package `gstat` constructed using the
  function `vgm`. This object contains the variogram model that will be
  fit to the data. The nugget is taken as the measurement error when
  `est_error = TRUE`. If unspecified, the variogram used is
  `gstat::vgm(1, "Lin", d, 1)`, where `d` is approximately one third of
  the maximum distance between any two data points

- K_type:

  the parameterisation used for the basis-function covariance matrix,
  `K`. If `method` = "EM", `K_type` can be "unstructured" or
  "block-exponential". If `method` = "TMB", `K_type` can be "precision"
  or "block-exponential". The default is "block-exponential", however if
  `FRK()` is used and `method` = "TMB", for computational reasons
  `K_type` is set to "precision"

- n_EM:

  (applicable only if `method` = "EM") maximum number of iterations for
  the EM algorithm

- tol:

  (applicable only if `method` = "EM") convergence tolerance for the EM
  algorithm

- method:

  parameter estimation method to employ. Currently "EM" and "TMB" are
  supported

- lambda:

  (applicable only if `K_type` = "unstructured") ridge-regression
  regularisation parameter (0 by default). Can be a single number, or a
  vector (one parameter for each resolution)

- print_lik:

  (applicable only if `method` = "EM") flag indicating whether to plot
  log-likelihood vs. iteration after convergence of the EM estimation
  algorithm

- response:

  string indicating the assumed distribution of the response variable.
  It can be "gaussian", "poisson", "negative-binomial", "binomial",
  "gamma", or "inverse-gaussian". If `method` = "EM", only "gaussian"
  can be used. Two distributions considered in this framework, namely
  the binomial distribution and the negative-binomial distribution, have
  an assumed-known ‘size’ parameter and a ‘probability of success’
  parameter; see the details below for the exact parameterisations used,
  and how to provide these ‘size’ parameters

- link:

  string indicating the desired link function. Can be "log", "identity",
  "logit", "probit", "cloglog", "reciprocal", or "reciprocal-squared".
  Note that only sensible link-function and response-distribution
  combinations are permitted. If `method` = "EM", only "identity" can be
  used

- optimiser:

  (applicable only if `method` = "TMB") the optimising function used for
  model fitting when `method` = "TMB" (default is `nlminb`). Users may
  pass in a function object or a string corresponding to a named
  function. Optional parameters may be passed to `optimiser` via `...`.
  The only requirement of `optimiser` is that the first three arguments
  correspond to the initial parameters, the objective function, and the
  gradient, respectively (this may be achieved by simply constructing a
  wrapper function)

- fs_by_spatial_BAU:

  (applicable only in a spatio-temporal setting and if `method` = "TMB")
  if `TRUE`, then each spatial BAU is associated with its own fine-scale
  variance parameter; otherwise, a single fine-scale variance parameter
  is used

- known_sigma2fs:

  known value of the fine-scale variance parameter. If `NULL` (the
  default), the fine-scale variance parameter is estimated as usual. If
  `known_sigma2fs` is not `NULL`, the fine-scale variance is fixed to
  the supplied value; this may be a scalar, or vector of length equal to
  the number of spatial BAUs (if `fs_by_spatial_BAU = TRUE`)

- taper:

  positive numeric indicating the strength of the
  covariance/partial-correlation tapering. Only applicable if `K_type` =
  "block-exponential", or if `K_type` = "precision" and the the
  basis-functions are irregular or the manifold is not the plane. If
  `taper` is `NULL` (default) and `method` = "EM", no tapering is
  applied; if `method` = "TMB", tapering must be applied (for
  computational reasons), and we set it to 3 if it is unspecified

- simple_kriging_fixed:

  commit to simple kriging at the fitting stage? If `TRUE`, model
  fitting is faster, but the option to conduct universal kriging at the
  prediction stage is removed

- ...:

  other parameters passed on to [`auto_basis()`](auto_basis.md) and
  [`auto_BAUs()`](auto_BAUs.md) when calling `FRK()`, or the user
  specified function `optimiser()` when calling `FRK()` or `SRE.fit()`

- normalise_basis:

  flag indicating whether to normalise the basis functions so that they
  reproduce a stochastic process with approximately constant variance
  spatially

- include_fs:

  (applicable only if `method` = "TMB") flag indicating whether the
  fine-scale variation should be included in the model

- object:

  object of class `SRE` returned from the constructor `SRE()` containing
  all the parameters and information on the SRE model

- newdata:

  object of class `SpatialPoylgons`, `SpatialPoints`, or `STI`,
  indicating the regions or points over which prediction will be carried
  out. The BAUs are used if this option is not specified.

- obs_fs:

  flag indicating whether the fine-scale variation sits in the
  observation model (systematic error; indicated by `obs_fs = TRUE`) or
  in the process model (process fine-scale variation; indicated by
  `obs_fs = FALSE`, default). For non-Gaussian data models, and/or
  non-identity link functions, if `obs_fs = TRUE`, then the fine-scale
  variation is removed from the latent process \\Y\\; however, they are
  re-introduced for prediction of the conditonal mean ***μ*** and
  simulated data ***Z**^(\*)*

- pred_time:

  vector of time indices at which prediction will be carried out. All
  time points are used if this option is not specified

- covariances:

  (applicable only for `method` = "EM") logical variable indicating
  whether prediction covariances should be returned or not. If set to
  `TRUE`, a maximum of 4000 prediction locations or polygons are allowed

- nsim:

  number of i) MC samples at each location when using `predict` or ii)
  response vectors when using `simulate`

- type:

  (applicable only if `method` = "TMB") vector of strings indicating the
  quantities for which inference is desired. If "link" is in `type`,
  inference on the latent Gaussian process \\Y(\cdot)\\ is included; if
  "mean" is in `type`, inference on the mean process \\\mu(\cdot)\\ is
  included (and the probability process, \\\pi(\cdot)\\, if applicable);
  if "response" is in `type`, inference on the noisy data ***Z**^(\*)*
  is included

- k:

  (applicable only if `response` is "binomial" or "negative-binomial")
  vector of size parameters at each BAU

- percentiles:

  (applicable only if `method` = "TMB") a vector of scalars in (0, 100)
  specifying the desired percentiles of the posterior predictive
  distribution; if `NULL`, no percentiles are computed

- kriging:

  (applicable only if `method` = "TMB") string indicating the kind of
  kriging: "simple" ignores uncertainty due to estimation of the fixed
  effects, while "universal" accounts for this source of uncertainty

- new_BAU_data:

  A `data.frame` containing updated covariate values at the BAU level.
  Must have the same number of rows, column names (in the same order),
  and column types as the original BAUs data. If `NULL` (default), the
  original BAUs data will be used.

- random_effects:

  logical; if set to true, confidence intervals will also be provided
  for the random effects random effects ***γ*** (see \`?SRE\` for
  details on these random effects)

- conditional_fs:

  condition on the fitted fine-scale random effects?

## Details

The following details provide a summary of the model and basic workflow
used in FRK. See Zammit-Mangion and Cressie (2021) and Sainsbury-Dale,
Zammit-Mangion and Cressie (2023) for further details.

**Model description**

The hierarchical model implemented in FRK is a spatial generalised
linear mixed model (GLMM), which may be summarised as

*Z_(j) \| **μ**_(Z), ψ ~ EF(μ_(Z_(j)) , ψ);       j = 1, ..., m,*

***μ**_(Z) = **C**_(Z)**μ**,*

*g(**μ**) = **Y**,*

***Y** = **Tα** + **Gγ** + **Sη** + **ξ**,*

***η** ~ N(**0**, **K**),*

***ξ** ~ N(**0**, **Σ**_(ξ) ),*

***γ** ~ N(**0**, **Σ**_(γ) ),*

  

where *Z_(j)* denotes a datum, \\EF\\ corresponds to a probability
distribution in the exponential family with dispersion parameter
\\\psi\\, ***μ**_(Z)* is the vector containing the conditional
expectations of each datum, ***C**_(Z)* is a matrix which aggregates the
BAU-level mean process over the observation supports, ***μ*** is the
mean process evaluated over the BAUs, \\g\\ is a link function, ***Y***
is a latent Gaussian process evaluated over the BAUs, the matrix ***T***
contains regression covariates at the BAU level associated with the
fixed effects ***α***, the matrix ***G*** is a design matrix at the BAU
level associated with random effects ***γ***, the matrix ***S***
contains basis-function evaluations over the BAUs associated with
basis-function random effects ***η***, and ***ξ*** is a vector
containing fine-scale variation at the BAU level.

The prior distribution of the random effects, ***γ***, is a mean-zero
multivariate Gaussian with diagonal covariance matrix, with each group
of random effects associated with its own variance parameter. These
variance parameters are estimated during model fitting.

The prior distribution of the basis-function coefficients, ***η***, is
formulated using either a covariance matrix ***K*** or precision matrix
***Q***, depending on the argument `K_type`. The parameters of these
matrices are estimated during model fitting.

The prior distribution of the fine-scale random effects, ***ξ***, is a
mean-zero multivariate Gaussian with diagonal covariance matrix,
***Σ**_(ξ)*. By default, ***Σ**_(ξ) = σ²_(ξ)**V***, where ***V*** is a
known, positive-definite diagonal matrix whose elements are provided in
the field `fs` in the BAUs. In the absence of problem specific
fine-scale information, `fs` can simply be set to 1, so that ***V** =
**I***. In a spatio-temporal setting, another model for ***Σ**_(ξ)* can
be used by setting `fs_by_spatial_BAU = TRUE`, in which case each
spatial BAU is associated with its own fine-scale variance parameter
(see Sainsbury-Dale et al., 2023, Sec. 2.6). In either case, the
fine-scale variance parameter(s) are either estimated during model
fitting, or provided by the user via the argument `known_sigma2fs`.

*Gaussian data model with an identity link function*

When the data is Gaussian, and an identity link function is used, the
preceding model simplifies considerably: Specifically,

***Z** = **C**_(Z)**Y** + **C**_(Z)**δ** + **e**,  
*

where ***Z*** is the data vector, ***δ*** is systematic error at the BAU
level, and ***e*** represents independent measurement error.

*Distributions with size parameters*

Two distributions considered in this framework, namely the binomial
distribution and the negative-binomial distribution, have an
assumed-known ‘size’ parameter and a ‘probability of success’ parameter.
Given the vector of size parameters associated with the data,
***k**_(Z)*, the parameterisation used in FRK assumes that *Z_(j)*
represents either the number of \`successes' from *k_(Z_(j))* trials
(binomial data model) or that it represents the number of failures
before *k_(Z_(j))* successes (negative-binomial data model).

When model fitting, the BAU-level size parameters ***k*** are needed.
The user must supply these size parameters either through the data or
though the BAUs. How this is done depends on whether the data are areal
or point-referenced, and whether they overlap common BAUs or not. The
simplest case is when each observation is associated with a single BAU
only and each BAU is associated with at most one observation support;
then, it is straightforward to assign elements from ***k**_(Z)* to
elements of ***k*** and vice-versa, and so the user may provide either
***k*** or ***k**_(Z)*. If each observation is associated with exactly
one BAU, but some BAUs are associated with multiple observations, the
user must provide ***k**_(Z)*, which is used to infer ***k*** ; in
particular, *k_(i) = Σ_(j∈a_(i)) k_(Z_(j))* , \\i = 1, \dots, N\\, where
*a_(i)* denotes the indices of the observations associated with BAU
*A_(i)*. If one or more observations encompass multiple BAUs, ***k***
must be provided with the BAUs, as we cannot meaningfully distribute
*k_(Z_(j))* over multiple BAUs associated with datum *Z_(j)*. In this
case, we infer ***k**_(Z)* using *k_(Z_(j)) = Σ_(i∈c_(j)) k_(i)* , \\j =
1, \dots, m\\, where *c_(j)* denotes the indices of the BAUs associated
with observation *Z_(j)*.

**Set-up**

`SRE()` constructs a spatial random effects model from the user-defined
formula, data object (a list of spatially-referenced data), basis
functions and a set of Basic Areal Units (BAUs). It first takes each
object in the list `data` and maps it to the BAUs – this entails binning
point-referenced data into the BAUs (and averaging within the BAU if
`average_in_BAU = TRUE`), and finding which BAUs are associated with
observations. Following this, the incidence matrix, ***C**_(Z)*, is
constructed. All required matrices (***S***, ***T***, ***C**_(Z)*, etc.)
are constructed within `SRE()` and returned as part of the `SRE` object.
`SRE()` also intitialises the parameters and random effects using
sensible defaults. Please see [`SRE-class`](SRE-class.md) for more
details. The functions [`observed_BAUs()`](observed_BAUs.md) and
[`unobserved_BAUs()`](observed_BAUs.md) return the indices of the
observed and unobserved BAUs, respectively.

To include random effects in FRK please follow the notation as used in
lme4. For example, to add a random effect according to a variable `fct`,
simply add \``(1 | fct)`' to the formula used when calling `FRK()` or
`SRE()`. Note that FRK only supports simple, uncorrelated random effects
and that a formula term such as '`(1 + x | fct)`' will throw an error
(since in lme4 parlance this implies that the random effect
corresponding to the intercept and the slope are correlated). If one
wishes to model a an intercept and linear trend for each level in `fct`,
then one can force the intercept and slope terms to be uncorrelated by
using the notation "`(x || fct)`", which is shorthand for
"`(1 | fct) + (x - 1 | x2)`".

**Model fitting**

`SRE.fit()` takes an object of class `SRE` and estimates all unknown
parameters, namely the covariance matrix ***K***, the fine scale
variance (*σ²_(ξ)* or *σ²_(δ)*, depending on whether Case 1 or Case 2 is
chosen; see the vignette "FRK_intro") and the regression parameters
***α***. There are two methods of model fitting currently implemented,
both of which implement maximum likelihood estimation (MLE).

- MLE via the expectation maximisation (EM) algorithm. :

  This method is implemented only for Gaussian data and an identity link
  function. The log-likelihood (given in Section 2.2 of the vignette) is
  evaluated at each iteration at the current parameter estimate.
  Optimation continues until convergence is reached (when the
  log-likelihood stops changing by more than `tol`), or when the number
  of EM iterations reaches `n_EM`. The actual computations for the
  E-step and M-step are relatively straightforward. The E-step contains
  an inverse of an \\r \times r\\ matrix, where \\r\\ is the number of
  basis functions which should not exceed 2000. The M-step first updates
  the matrix ***K***, which only depends on the sufficient statistics of
  the basis-function coefficients ***η***. Then, the regression
  parameters ***α*** are updated and a simple optimisation routine (a
  line search) is used to update the fine-scale variance *σ²_(δ)* or
  *σ²_(ξ)*. If the fine-scale errors and measurement random errors are
  homoscedastic, then a closed-form solution is available for the update
  of *σ²_(ξ)* or *σ²_(δ)*. Irrespectively, since the updates of ***α***,
  and *σ²_(δ)* or *σ²_(ξ)*, are dependent, these two updates are
  iterated until the change in *σ²_(.)* is no more than 0.1%.

- MLE via `TMB`. :

  This method is implemented for all available data models and link
  functions offered by FRK. Furthermore, this method facilitates the
  inclusion of many more basis function than possible with the EM
  algorithm (in excess of 10,000). `TMB` applies the Laplace
  approximation to integrate out the latent random effects from the
  complete-data likelihood. The resulting approximation of the marginal
  log-likelihood, and its derivatives with respect to the parameters,
  are then called from within `R` using the optimising function
  `optimiser` (default
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)).

*Wrapper for set-up and model fitting*

The function `FRK()` acts as a wrapper for the functions `SRE()` and
`SRE.fit()`. An added advantage of using `FRK()` directly is that it
automatically generates BAUs and basis functions based on the data.
Hence `FRK()` can be called using only a list of data objects and an `R`
formula, although the `R` formula can only contain space or time as
covariates when BAUs are not explicitly supplied with the covariate
data.

**Prediction**

Once the parameters are estimated, the `SRE` object is passed onto the
function [`predict()`](https://rdrr.io/r/stats/predict.html) in order to
carry out optimal predictions over the same BAUs used to construct the
SRE model with `SRE()`. The first part of the prediction process is to
construct the matrix ***S*** over the prediction polygons. This is made
computationally efficient by treating the prediction over polygons as
that of the prediction over a combination of BAUs. This will yield valid
results only if the BAUs are relatively small. Once the matrix ***S***
is found, a standard Gaussian inversion (through conditioning) using the
estimated parameters is used for prediction.

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the BAUs (or
an object specified in `newdata`), which are of class
`SpatialPixelsDataFrame`, `SpatialPolygonsDataFrame`, or `STFDF`, with
predictions and uncertainty quantification added. If `method` = "TMB",
the returned object is a list, containing the previously described
predictions, and a list of Monte Carlo samples. The predictions and
uncertainties can be easily plotted using [`plot`](plot.md) or `spplot`
from the package `sp`.

## References

Zammit-Mangion, A. and Cressie, N. (2021). FRK: An R package for spatial
and spatio-temporal prediction with large datasets. Journal of
Statistical Software, 98(4), 1-48. doi:10.18637/jss.v098.i04.

Sainsbury-Dale, M. and Zammit-Mangion, A. and Cressie, N. (2024)
Modelling Big, Heterogeneous, Non-Gaussian Spatial and Spatio-Temporal
Data using FRK. Journal of Statistical Software, 108(10), 1–39.
doi:10.18637/jss.v108.i10.

## See also

[`SRE-class`](SRE-class.md) for details on the SRE object internals,
[`auto_basis`](auto_basis.md) for automatically constructing basis
functions, and [`auto_BAUs`](auto_BAUs.md) for automatically
constructing BAUs.

## Examples

``` r
library("FRK")
library("sp")
## Generate process and data
m <- 250                                                   # Sample size
zdf <- data.frame(x = runif(m), y= runif(m))               # Generate random locs
zdf$Y <- 3 + sin(7 * zdf$x) + cos(9 * zdf$y)               # Latent process
zdf$z <- rnorm(m, mean = zdf$Y)                            # Simulate data
coordinates(zdf) = ~x+y                                    # Turn into sp object

## Construct BAUs and basis functions
BAUs <- auto_BAUs(manifold = plane(), data = zdf,
                  nonconvex_hull = FALSE, cellsize = c(0.03, 0.03), type="grid")
BAUs$fs <- 1 # scalar fine-scale covariance matrix
basis <- auto_basis(manifold =  plane(), data = zdf, nres = 2)

## Construct the SRE model
S <- SRE(f = z ~ 1, list(zdf), basis = basis, BAUs = BAUs)
#> Loading required namespace: gstat

## Fit with 2 EM iterations so to take as little time as possible
S <- SRE.fit(S, n_EM = 2, tol = 0.01, print_lik = TRUE)
#> NOTE: In FRK >2.0 simple_kriging_fixed = FALSE by default, and hence
#>       universal kriging is done by default. However this is only the case when
#>       method = 'TMB'. When method = 'EM', simple kriging is done, irrespective of
#>       what the argument simple_kriging_fixed is set to.
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
#> Maximum EM iterations reached


## Check fit info, final log-likelihood, and estimated regression coefficients
info_fit(S)
#> $method
#> [1] "EM"
#> 
#> $num_iterations
#> [1] 2
#> 
#> $sigma2fshat_equal_0
#> [1] 1
#> 
#> $converged
#> [1] 0
#> 
#> $plot_lik
#> $plot_lik$x
#> [1] 1 2
#> 
#> $plot_lik$llk
#> [1] -356.1047 -350.2567
#> 
#> $plot_lik$ylab
#> [1] "log likelihood"
#> 
#> $plot_lik$xlab
#> [1] "EM iteration"
#> 
#> 
#> $time
#>    user  system elapsed 
#>   0.240   0.372   0.205 
#> 
logLik(S)
#> [1] -349.8645
coef(S)
#> Intercept 
#>  3.097504 

## Predict over BAUs
pred <- predict(S)

## Plot
if (FALSE) { # \dontrun{
plotlist <- plot(S, pred)
ggpubr::ggarrange(plotlist = plotlist, nrow = 1, align = "hv", legend = "top")} # }
```
