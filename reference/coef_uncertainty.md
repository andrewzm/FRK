# Uncertainty quantification of the fixed effects

Compute confidence intervals for the fixed effects (upper and lower
bound specifed by percentiles; default 90% confidence central interval)

## Usage

``` r
coef_uncertainty(
  object,
  percentiles = c(5, 95),
  nsim = 400,
  random_effects = FALSE
)
```

## Arguments

- object:

  object of class `SRE` returned from the constructor [`SRE()`](SRE.md)
  containing all the parameters and information on the SRE model

- percentiles:

  (applicable only if `method` = "TMB") a vector of scalars in (0, 100)
  specifying the desired percentiles of the posterior predictive
  distribution; if `NULL`, no percentiles are computed

- nsim:

  number of Monte Carlo samples used to compute the confidence intervals

- random_effects:

  logical; if set to true, confidence intervals will also be provided
  for the random effects random effects ***γ*** (see \`?SRE\` for
  details on these random effects)
