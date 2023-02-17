# FRK: An R Software package for spatial and spatio-temporal prediction
# with large datasets.
# Copyright (c) 2017 University of Wollongong
# Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' @name SRE
#'
#' @title Construct SRE object, fit and predict
#' @description The Spatial Random Effects (SRE) model is the central object in \pkg{FRK}. The function \code{FRK()} provides a wrapper for the construction and estimation of the SRE object from data, using the functions \code{SRE()} (the object constructor) and \code{SRE.fit()} (for fitting it to the data). Please see \code{\link{SRE-class}} for more details on the SRE object's properties and methods.
#' @param f \code{R} formula relating the dependent variable (or transformations thereof) to covariates
#' @param data list of objects of class \code{SpatialPointsDataFrame}, \code{SpatialPolygonsDataFrame}, \code{STIDF}, or  \code{STFDF}. If using space-time objects, the data frame must have another field, \code{t}, containing the time index of the data point
#' @param basis object of class \code{Basis} (or \code{TensorP_Basis})
#' @param BAUs object of class \code{SpatialPolygonsDataFrame}, \code{SpatialPixelsDataFrame}, \code{STIDF}, or \code{STFDF}. The object's data frame must contain covariate information as well as a field \code{fs} describing the fine-scale variation up to a constant of proportionality. If the function \code{FRK()} is used directly, then BAUs are created automatically, but only coordinates can then be used as covariates
#' @param est_error (applicable only if \code{response} = "gaussian") flag indicating whether the measurement-error variance should be estimated from variogram techniques. If this is set to 0, then \code{data} must contain a field \code{std}. Measurement-error estimation is currently not implemented for spatio-temporal datasets
#' @param include_fs (applicable only if \code{method} = "TMB") flag indicating whether the fine-scale variation should be included in the model
#' @param average_in_BAU if \code{TRUE}, then multiple data points falling in the same BAU are averaged; the measurement error of the averaged data point is taken as the average of the individual measurement errors
#' @param sum_variables if \code{average_in_BAU == TRUE}, the string \code{sum_variables} indicates which data variables (can be observations or covariates) are to be summed rather than averaged
#' @param normalise_wts if \code{TRUE}, the rows of the incidence matrices \ifelse{html}{\out{<i><b>C</b><sub>Z</sub></i>}}{\eqn{\boldsymbol{C}_Z}{C_Z}} and \ifelse{html}{\out{<i><b>C</b><sub>P</sub></i>}}{\eqn{\boldsymbol{C}_P}{C_P}} are normalised to sum to 1, so that the mapping represents a weighted average; if false, no normalisation of the weights occurs (i.e., the mapping corresponds to a weighted sum)
#' @param fs_by_spatial_BAU (applicable only in a spatio-temporal setting and if \code{method} = "TMB") if \code{TRUE}, then each spatial BAU is associated with its own fine-scale variance parameter; otherwise, a single fine-scale variance parameter is used
#' @param fs_model if "ind" then the fine-scale variation is independent at the BAU level. Only the independent model is allowed for now, future implementation will include CAR/ICAR (in development)
# #' If "ICAR", then an ICAR model for the fine-scale variation is placed on the BAUs
#' @param vgm_model (applicable only if \code{response} = "gaussian") an object of class \code{variogramModel} from the package \code{gstat} constructed using the function \code{vgm}. This object contains the variogram model that will be fit to the data. The nugget is taken as the measurement error when \code{est_error = TRUE}. If unspecified, the variogram used is \code{gstat::vgm(1, "Lin", d, 1)}, where \code{d} is approximately one third of the maximum distance between any two data points
#' @param K_type the parameterisation used for the basis-function covariance matrix, \code{K}. If \code{method} = "EM", \code{K_type} can be "unstructured" or "block-exponential". If \code{method} = "TMB", \code{K_type} can be "precision" or "block-exponential". The default is "block-exponential", however if \code{FRK()} is used and \code{method} = "TMB", for computational reasons \code{K_type} is set to "precision"
#' @param normalise_basis flag indicating whether to normalise the basis functions so that they reproduce a stochastic process with approximately constant variance spatially
#' @param object object of class \code{SRE} returned from the constructor \code{SRE()} containing all the parameters and information on the SRE model
#' @param n_EM (applicable only if \code{method} = "EM") maximum number of iterations for the EM algorithm
#' @param tol (applicable only if \code{method} = "EM") convergence tolerance for the EM algorithm
#' @param method parameter estimation method to employ. Currently "EM" and "TMB" are supported
#' @param lambda (applicable only if \code{K_type} = "unstructured") ridge-regression regularisation parameter (0 by default). Can be a single number, or a vector (one parameter for each resolution)
#' @param print_lik (applicable only if \code{method} = "EM") flag indicating whether to plot log-likelihood vs. iteration after convergence of the EM estimation algorithm
# #' @param use_centroid flag indicating whether the basis functions are averaged over the BAU, or whether the basis functions are evaluated at the BAUs centroid in order to construct the matrix \ifelse{html}{\out{<i> <b>S</b> </i>}}{\eqn{\boldsymbol{S}}{S}{S}}. The flag can safely be set when the basis functions are approximately constant over the BAUs in order to reduce computational time
#' @param newdata object of class \code{SpatialPoylgons}, \code{SpatialPoints}, or \code{STI}, indicating the regions or points over which prediction will be carried out. The BAUs are used if this option is not specified.
#' @param obs_fs flag indicating whether the fine-scale variation sits in the observation model (systematic error; indicated by \code{obs_fs = TRUE}) or in the process model (process fine-scale variation; indicated by \code{obs_fs = FALSE}, default). For non-Gaussian data models, and/or non-identity link functions, if \code{obs_fs = TRUE}, then the fine-scale variation is removed from the latent process \eqn{Y}; however, they are re-introduced for prediction of the conditonal mean \ifelse{html}{\out{<i> <b> &mu; </b> </i>}}{\eqn{\boldsymbol{\mu}}{mu}} and simulated data \ifelse{html}{\out{<i> <b>Z</b><sup>*</sup> </i>}}{\eqn{\boldsymbol{Z}^*}{Z*}}
#' @param pred_time vector of time indices at which prediction will be carried out. All time points are used if this option is not specified
#' @param covariances (applicable only for \code{method} = "EM") logical variable indicating whether prediction covariances should be returned or not. If set to \code{TRUE}, a maximum of 4000 prediction locations or polygons are allowed
#' @param response string indicating the assumed distribution of the response variable. It can be "gaussian", "poisson", "negative-binomial", "binomial", "gamma", or "inverse-gaussian". If \code{method} = "EM", only "gaussian" can be used. Two distributions considered in this framework, namely the binomial distribution and the negative-binomial distribution, have an assumed-known ‘size’ parameter and a ‘probability of success’ parameter; see the details below for the exact parameterisations used, and how to provide these ‘size’ parameters
#' @param link  string indicating the desired link function. Can be "log", "identity", "logit", "probit", "cloglog", "reciprocal", or "reciprocal-squared". Note that only sensible link-function and response-distribution combinations are permitted. If \code{method} = "EM", only "identity" can be used
#' @param taper positive numeric indicating the strength of the covariance/partial-correlation tapering. Only applicable if \code{K_type} = "block-exponential", or if \code{K_type} = "precision" and the the basis-functions are irregular or the manifold is not the plane. If \code{taper} is \code{NULL} (default) and \code{method} = "EM", no tapering is applied; if \code{method} = "TMB", tapering must be applied (for computational reasons), and we set it to 3 if it is unspecified
#' @param optimiser (applicable only if \code{method} = "TMB") the optimising function used for model fitting when \code{method} = "TMB" (default is \code{nlminb}). Users may pass in a function object or a string corresponding to a named function. Optional parameters may be passed to \code{optimiser} via \code{...}. The only requirement of \code{optimiser} is that the first three arguments correspond to the initial parameters, the objective function, and the gradient, respectively (this may be achieved by simply constructing a wrapper function)
#' @param known_sigma2fs known value of the fine-scale variance parameter. If \code{NULL} (the default), the fine-scale variance parameter is estimated as usual. If \code{known_sigma2fs} is not \code{NULL}, the fine-scale variance is fixed to the supplied value; this may be a scalar, or vector of length equal to the number of spatial BAUs (if \code{fs_by_spatial_BAU = TRUE})
#' @param kriging (applicable only if \code{method} = "TMB") string indicating the kind of kriging: "simple" ignores uncertainty due to estimation of the fixed effects, while "universal" accounts for this source of uncertainty
#' @param type (applicable only if \code{method} = "TMB") vector of strings indicating the quantities for which inference is desired. If "link" is in \code{type}, inference on the latent Gaussian process \eqn{Y(\cdot)}{Y(.)} is included; if "mean" is in \code{type}, inference on the mean process \eqn{\mu(\cdot)}{\mu(.)} is included (and the probability process, \eqn{\pi(\cdot)}{\pi(.)},  if applicable); if "response" is in \code{type}, inference on the noisy data \ifelse{html}{\out{<i> <b>Z</b><sup>*</sup> </i>}}{\eqn{\boldsymbol{Z}^*}{Z*}} is included
#' @param nsim number of i) MC samples at each location when using \code{predict} or ii) response vectors when using \code{simulate}
#' @param k (applicable only if \code{response} is "binomial" or "negative-binomial") vector of size parameters at each BAU
#' @param percentiles (applicable only if \code{method} = "TMB") a vector of scalars in (0, 100) specifying the desired percentiles of the posterior predictive distribution; if \code{NULL}, no percentiles are computed
#' @param simple_kriging_fixed commit to simple kriging at the fitting stage? If \code{TRUE}, model fitting is faster, but the option to conduct universal kriging at the prediction stage is removed
#' @param conditional_fs condition on the fitted fine-scale random effects?
#' @param ... other parameters passed on to \code{auto_basis()} and \code{auto_BAUs()} when calling \code{FRK()}, or the user specified function \code{optimiser()} when calling \code{FRK()} or \code{SRE.fit()}
#' @details
#'
#' The following details provide a summary of the model and basic workflow
#' used in \pkg{FRK}. See Zammit-Mangion and Cressie
#' (2021) and Sainsbury-Dale, Zammit-Mangion and Cressie (2021) for further details.
#'
#' \strong{Model description}
#'
#' The hierarchical model implemented in \pkg{FRK} is a spatial generalised
#' linear mixed model (GLMM), which may be summarised as
#'
#' \ifelse{html}{\out{<div style="text-align:center"> <i> Z<sub>j</sub> | <b>&mu;</b><sub>Z</sub>, &psi; ~ EF(&mu;<sub>Z<sub>j</sub></sub> , &psi;); &nbsp; &nbsp; &nbsp; j = 1, ..., m, </i></div>}}{\deqn{Z_j \mid \boldsymbol{\mu}_{Z}, \psi \sim EF(\mu_{Z_j}, \psi); \quad j = 1, \dots, m,}{Z_j | \mu_{Z}, \psi ~ EF(\mu_{Z_j}, \psi);  j = 1, \dots, m,}}
#' \ifelse{html}{\out{<div style="text-align:center"> <i> <b>&mu;</b><sub>Z</sub> = <b>C</b><sub>Z</sub> <b>&mu;</b>, </i></div>}}{\deqn{\boldsymbol{\mu}_Z = \boldsymbol{C}_Z\boldsymbol{\mu}}{\mu_Z = C_Z \mu}}
#' \ifelse{html}{\out{<div style="text-align:center"> <i> g(<b>&mu;</b>) = <b>Y</b>, </i></div>}}{\deqn{g(\boldsymbol{\mu}) = \boldsymbol{Y}}{g(\mu) = Y}}
#' \ifelse{html}{\out{<div style="text-align:center"> <i> <b>Y</b> = <b>T&alpha;</b> + <b>S&eta;</b> + <b>&xi;</b>, </i></div>}}{\deqn{\boldsymbol{Y} = \boldsymbol{T}\boldsymbol{\alpha} + \boldsymbol{S}\boldsymbol{\eta} + \boldsymbol{\xi}}{Y = T\alpha + S\eta + \xi}}
#' \ifelse{html}{\out{<div style="text-align:center"> <i> <b>&eta;</b> | <b>&vartheta;</b> ~ N(<b>0</b>, <b>K</b>),</i></div>}}{\deqn{\boldsymbol{\eta} \mid \boldsymbol{\vartheta} \sim N(\boldsymbol{0}, \boldsymbol{K})}{\eta | \vartheta ~ N(0, K)}}
#' \ifelse{html}{\out{<div style="text-align:center"> <i> <b>&xi;</b> ~ N(<b>0</b>, <b>&Sigma;</b><sub>&xi;</sub>),</i></div><br>}}{\deqn{\boldsymbol{\xi} \mid \sigma^2_\xi \sim N(\boldsymbol{0}, \boldsymbol{\Sigma}_\xi),}{\xi | \sigma^2_\xi ~ N(0, \Sigma_\xi),}}
#'
#' where \ifelse{html}{\out{<i>Z<sub>j</sub></i>}}{\eqn{Z_j}} denotes a datum, \eqn{EF}  corresponds to a probability
#' distribution in the exponential family with dispersion parameter \eqn{\psi},
#' \ifelse{html}{\out{<i> <b>&mu;</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{\mu}_Z}{\mu_Z}} is the vector containing the conditional expectations of each datum,
#' \ifelse{html}{\out{<i> <b>C</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{C}_Z}{C_Z}} is a matrix which aggregates the BAU-level mean process over the observation supports,
#' \ifelse{html}{\out{<i> <b>&mu;</b> </i>}}{\eqn{\boldsymbol{\mu}}{\mu}} is the mean process evaluated over the BAUs, \eqn{g} is a link function,
#' \ifelse{html}{\out{<i> <b>Y</b> </i>}}{\eqn{\boldsymbol{Y}}{Y}} is a latent Gaussian process evaluated over the BAUs,
#' the matrix \ifelse{html}{\out{<i> <b>T</b> </i>}}{\eqn{\boldsymbol{T}}{T}} contains regression covariates at the BAU level associated with the fixed effects
#' \ifelse{html}{\out{<i> <b>&alpha;</b> </i>}}{\eqn{\boldsymbol{\alpha}}{\alpha}},
#' the matrix \ifelse{html}{\out{<i> <b>S</b> </i>}}{\eqn{\boldsymbol{S}}{S}} contains basis function evaluations over the BAUs,
#' \ifelse{html}{\out{<i> <b>&eta;</b> </i>}}{\eqn{\boldsymbol{\eta}}{\eta}} are the random coefficients associated with the basis functions, and \ifelse{html}{\out{<i> <b>&xi;</b> </i>}}{\eqn{\boldsymbol{\xi}}{\xi}} is a vector containing fine-scale variation at the BAU level.
#'
#' The prior distribution of the basis-function coefficients, \ifelse{html}{\out{<i> <b>&eta;</b> </i>}}{\eqn{\boldsymbol{\eta}}{\eta}}, are formulated
#' using either a covariance matrix \ifelse{html}{\out{<i> <b>K</b> </i>}}{\eqn{\boldsymbol{K}}{K}} or precision matrix \ifelse{html}{\out{<i> <b>Q</b> </i>}}{\eqn{\boldsymbol{Q}}{Q}}, depending on the argument
#' \code{K_type}; the parameters of these matrices, \ifelse{html}{\out{<i> <b>&vartheta;</b> </i>}}{\eqn{\boldsymbol{\vartheta}}{\vartheta}}, are estimated during model
#' fitting.
#' The covariance matrix of \ifelse{html}{\out{<i> <b>&xi;</b> </i>}}{\eqn{\boldsymbol{\xi}}{\xi}},
#' \ifelse{html}{\out{<i> <b>&Sigma;</b><sub>&xi;</sub> </i>}}{\eqn{\boldsymbol{\Sigma}_\xi}{\Sigma_\xi}},
#' is diagonal.
#' By default, \ifelse{html}{\out{<i> <b>&Sigma;</b><sub>&xi;</sub> = &sigma;<sup>2</sup><sub>&xi;</sub><b>V</b> </i>}}{\eqn{\boldsymbol{\Sigma}_\xi = \sigma^2_\xi \boldsymbol{V}}{\Sigma_\xi = \sigma^2_\xi V}}, where \ifelse{html}{\out{<i> <b>V</b> </i>}}{\eqn{\boldsymbol{V}}{V}} is a
#' known, positive-definite diagonal matrix whose elements are provided in the
#' field `fs' in the BAUs; in the absence of problem
#' specific fine-scale information, `fs' can simply be set to 1, so that
#' \ifelse{html}{\out{<i> <b>V</b> = <b>I</b> </i>}}{\eqn{\boldsymbol{V} = \boldsymbol{I}}{V = I}}.
#' In a spatio-temporal setting, another model for \ifelse{html}{\out{<i> <b>&Sigma;</b><sub>&xi;</sub> </i>}}{\eqn{\boldsymbol{\Sigma}_\xi}{\Sigma_\xi}}
#' can be used by setting \code{fs_by_spatial_BAU = TRUE}, in which case each
#' spatial BAU is associated with its own fine-scale variance parameter (see
#' Section 2.6 of Sainsbury-Dale, Zammit-Mangion and Cressie (2021) for details).
#' In either case, the fine-scale variance parameter(s) are either estimated during model fitting, or provided by
#' the user via the argument \code{known_sigma2fs}.
#'
#' \emph{Gaussian data model with an identity link function}
#'
#'
#' When the data is Gaussian, and an identity link function is used, the preceding
#' model simplifies considerably: Specifically,
#'
#' \ifelse{html}{\out{
#' <div style="text-align:center"><i> <b>Z</b> = <b>C</b><sub>Z</sub><b>Y</b> + <b>C</b><sub>Z</sub><b>&delta;</b> + <b>e</b>,<br> </i></div>
#' }}{
#' \deqn{\boldsymbol{Z} = \boldsymbol{C}_Z\boldsymbol{Y} + \boldsymbol{C}_Z\boldsymbol{\delta} + \boldsymbol{e},}{Z = C_ZY + C_Z \delta + e,}
#' }
#'
#' where
#' \ifelse{html}{\out{<i><b>Z</b></i>}}{\eqn{\boldsymbol{Z}}{Z}} is the data vector,
#' \ifelse{html}{\out{<i><b>&delta;</b></i>}}{\eqn{\boldsymbol{\delta}}{\delta}} is systematic error at the BAU level, and
#' \ifelse{html}{\out{<i><b>e</b></i>}}{\eqn{\boldsymbol{e}}{e}} represents independent measurement error.
#'
#' \emph{Distributions with size parameters}
#'
#' Two distributions considered in this framework, namely the binomial
#' distribution and the negative-binomial distribution, have an assumed-known
#' ‘size’ parameter and a ‘probability of success’ parameter.
#' Given the vector of size parameters associated with the data,
#' \ifelse{html}{\out{<i> <b>k</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{k}_Z}{k_Z}}, the parameterisation used in \pkg{FRK} assumes that
#' \ifelse{html}{\out{<i>Z<sub>j</sub></i>}}{\eqn{Z_j}} represents either the number of `successes' from
#' \ifelse{html}{\out{<i>k<sub>Z<sub>j</sub></sub></i>}}{\eqn{k_{Z_j}}} trials (binomial data model) or that it represents the number of failures before
#' \ifelse{html}{\out{<i>k<sub>Z<sub>j</sub></sub></i>}}{\eqn{k_{Z_j}}} successes (negative-binomial data model).
#'
#' When model fitting, the BAU-level size parameters
#' \ifelse{html}{\out{<i> <b> k </b> </i>}}{\eqn{\boldsymbol{k}}{k}} are needed.
#' The user must supply these size parameters either through the data or though
#' the BAUs. How this is done depends on whether the data are areal or
#' point-referenced, and whether they overlap common BAUs or not.
#' The simplest case is when each observation is associated with a single BAU
#' only and each BAU is associated with at most one observation support; then,
#' it is straightforward to assign elements from
#' \ifelse{html}{\out{<i> <b>k</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{k}_Z}{k_Z}} to elements of
#' \ifelse{html}{\out{<i> <b> k </b> </i>}}{\eqn{\boldsymbol{k}}{k}} and vice-versa, and so the user may provide either
#' \ifelse{html}{\out{<i> <b> k </b> </i>}}{\eqn{\boldsymbol{k}}{k}} or
#' \ifelse{html}{\out{<i> <b>k</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{k}_Z}{k_Z}}.
#' If each observation is associated with
#' exactly one BAU, but some BAUs are associated with multiple observations,
#' the user must provide \ifelse{html}{\out{<i> <b>k</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{k}_Z}{k_Z}}, which is used to infer
#' \ifelse{html}{\out{<i> <b> k </b> </i>}}{\eqn{\boldsymbol{k}}{k}}; in
#' particular,
#' \ifelse{html}{\out{<i>k<sub>i</sub> = &Sigma;<sub>j&isin;a<sub>i</sub></sub> k<sub>Z<sub>j</sub></sub> </i> }}{\eqn{k_i = \sum_{j \in a_i} k_{Z_j}}},
#' \eqn{i = 1, \dots, N}, where
#' \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}}
#' denotes the indices of the observations associated with BAU
#' \ifelse{html}{\out{<i>A<sub>i</sub></i>}}{\eqn{A_i}}.
#' If one or more observations encompass multiple BAUs,
#' \ifelse{html}{\out{<i> <b> k </b> </i>}}{\eqn{\boldsymbol{k}}{k}}
#' must be provided with the BAUs, as we cannot meaningfully
#' distribute
#' \ifelse{html}{\out{<i>k<sub>Z<sub>j</sub></sub></i>}}{\eqn{k_{Z_j}}}
#' over multiple BAUs associated with datum
#' \ifelse{html}{\out{<i>Z<sub>j</sub></i>}}{\eqn{Z_j}}.
#' In this case, we infer
#' \ifelse{html}{\out{<i> <b>k</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{k}_Z}{k_Z}} using
#' \ifelse{html}{\out{<i>k<sub>Z<sub>j</sub></sub> = &Sigma;<sub>i&isin;c<sub>j</sub></sub> k<sub>i</sub> </i> }}{\eqn{k_{Z_j} = \sum_{i \in c_j} k_i}},
#' \eqn{j = 1, \dots, m}, where
#' \ifelse{html}{\out{<i>c<sub>j</sub></i>}}{\eqn{c_j}}
#' denotes the indices of the BAUs associated with observation
#' \ifelse{html}{\out{<i>Z<sub>j</sub></i>}}{\eqn{Z_j}}.
#'
#'
#' \strong{Set-up}
#'
#' \code{SRE()}
#' constructs a spatial random effects model from the user-defined formula, data object (a list
#' of spatially-referenced data), basis functions and a set of Basic Areal Units (BAUs).
#' It first takes each object in the list \code{data} and maps it to the BAUs -- this
#' entails binning point-referenced data into the BAUs (and averaging within the
#' BAU if \code{average_in_BAU = TRUE}), and finding which BAUs are associated
#' with observations. Following this, the incidence matrix, \ifelse{html}{\out{<i> <b>C</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{C}_Z}{C_Z}}, is
#' constructed.
#' All required matrices (\ifelse{html}{\out{<i><b>S</b></i>}}{\eqn{\boldsymbol{S}}{S}}, \ifelse{html}{\out{<i> <b>T</b> </i>}}{\eqn{\boldsymbol{T}}{T}}, \ifelse{html}{\out{<i> <b>C</b><sub>Z</sub> </i>}}{\eqn{\boldsymbol{C}_Z}{C_Z}}, etc.)
#' are constructed within \code{SRE()} and returned as part of the \code{SRE} object.
#' \code{SRE()} also intitialises the parameters and random effects using
#' sensible defaults. Please see
#' \code{\link{SRE-class}} for more details.
#' The functions \code{observed_BAUs()} and \code{unobserved_BAUs()} return the
#' indices of the observed and unobserved BAUs, respectively.
#'
#'
#' \strong{Model fitting}
#'
#' \code{SRE.fit()} takes an object of class \code{SRE} and estimates all unknown
#' parameters, namely the covariance matrix \ifelse{html}{\out{<i><b>K</b></i>}}{\eqn{\boldsymbol{K}}{K}}, the fine scale variance
#' (\ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&xi;</sub></i>}}{\eqn{\sigma^2_{\xi}}} or \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&delta;</sub></i>}}{\eqn{\sigma^2_{\delta}}}, depending on whether Case 1
#' or Case 2 is chosen; see the vignette "FRK_intro") and the regression parameters \ifelse{html}{\out{<i> <b>&alpha;</b> </i>}}{\eqn{\boldsymbol{\alpha}}{\alpha}}.
#' There are two methods of model fitting currently implemented, both of which
#' implement maximum likelihood estimation (MLE).
#' \itemize{
#'  \item{MLE via the expectation maximisation
#'  (EM) algorithm. }{This method is implemented only
#'  for Gaussian data and an identity link function.
#'  The log-likelihood (given in Section 2.2 of the vignette) is evaluated at each
#' iteration at the current parameter estimate. Optimation continues until
#' convergence is reached (when the log-likelihood stops changing by more than
#' \code{tol}), or when the number of EM iterations reaches \code{n_EM}.
#' The actual computations for the E-step and M-step are relatively straightforward.
#' The E-step contains an inverse of an \eqn{r \times r}{r x r} matrix, where \eqn{r}
#' is the number of basis functions which should not exceed 2000. The M-step
#' first updates the matrix \ifelse{html}{\out{<i><b>K</b></i>}}{\eqn{\boldsymbol{K}}{K}}, which only depends on the sufficient
#' statistics of the basis-function coefficients \ifelse{html}{\out{<i> <b>&eta;</b> </i>}}{\eqn{\boldsymbol{\eta}}{\eta}}. Then, the regression
#' parameters \ifelse{html}{\out{<i> <b>&alpha;</b> </i>}}{\eqn{\boldsymbol{\alpha}}{\alpha}} are updated and a simple optimisation routine
#' (a line search) is used to update the fine-scale variance
#' \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&delta;</sub></i>}}{\eqn{\sigma^2_{\delta}}} or \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&xi;</sub></i>}}{\eqn{\sigma^2_{\xi}}}. If the fine-scale errors and
#' measurement random errors are homoscedastic, then a closed-form solution is
#' available for the update of \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&xi;</sub></i>}}{\eqn{\sigma^2_{\xi}}} or \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&delta;</sub></i>}}{\eqn{\sigma^2_{\delta}}}.
#' Irrespectively, since the updates of \ifelse{html}{\out{<i> <b>&alpha;</b> </i>}}{\eqn{\boldsymbol{\alpha}}{\alpha}}, and \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&delta;</sub></i>}}{\eqn{\sigma^2_{\delta}}}
#' or \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>&xi;</sub></i>}}{\eqn{\sigma^2_{\xi}}}, are dependent, these two updates are iterated until
#' the change in \ifelse{html}{\out{<i>&sigma;<sup>2</sup><sub>.</sub></i>}}{\eqn{\sigma^2_{\cdot}}} is no more than 0.1\%.}
#'  \item{MLE via \code{TMB}. }{This method is implemented for
#'  all available data models and link functions offered by \pkg{FRK}. Furthermore,
#'  this method faciliates the inclusion of many more basis function than possible
#'  with the EM algorithm (in excess of 10,000). \code{TMB} applies
#'  the Laplace approximation to integrate out the latent random effects from the
#'  complete-data likelihood. The resulting approximation of the marginal
#'  log-likelihood, and its derivatives with respect to the parameters, are then
#'  called from within \code{R} using the optimising function \code{optimiser}
#'  (default \code{nlminb()}).}
#' }
#'
#'
#' \emph{Wrapper for set-up and model fitting}
#'
#' The function \code{FRK()} acts as a wrapper for the functions \code{SRE()} and
#' \code{SRE.fit()}. An added advantage of using \code{FRK()} directly is that it
#' automatically generates BAUs and basis functions based on the data. Hence
#' \code{FRK()} can be called using only a list of data objects and an \code{R}
#' formula, although the \code{R} formula can only contain space or time as
#' covariates when BAUs are not explicitly supplied with the covariate data.
#'
#'
#' \strong{Prediction}
#'
#' Once the parameters are estimated, the \code{SRE} object is passed onto the
#' function \code{predict()} in order to carry out optimal predictions over the
#' same BAUs used to construct the SRE model with \code{SRE()}. The first part
#' of the prediction process is to construct the matrix \ifelse{html}{\out{<i> <b>S</b> </i>}}{\eqn{\boldsymbol{S}}{S}} over the
#' prediction polygons. This is made computationally efficient by treating the
#' prediction over polygons as that of the prediction over a combination of BAUs.
#' This will yield valid results only if the BAUs are relatively small. Once the
#' matrix \ifelse{html}{\out{<i> <b>S</b> </i>}}{\eqn{\boldsymbol{S}}{S}} is found, a standard Gaussian inversion (through conditioning)
#' using the estimated parameters is used for prediction.
#'
#' \code{predict()} returns the BAUs (or an object specified in \code{newdata}),
#' which are of class \code{SpatialPixelsDataFrame}, \code{SpatialPolygonsDataFrame},
#' or \code{STFDF}, with predictions and
#' uncertainty quantification added.
#' If \code{method} = "TMB", the returned object is a list, containing the
#' previously described predictions, and a list of Monte Carlo samples.
#' The predictions and uncertainties can be easily plotted using \code{\link{plot}}
#' or \code{spplot} from the package \code{sp}.
#' @seealso \code{\link{SRE-class}} for details on the SRE object internals,
#' \code{\link{auto_basis}} for automatically constructing basis functions, and
#' \code{\link{auto_BAUs}} for automatically constructing BAUs.
#' @references
#' Zammit-Mangion, A. and Cressie, N. (2021). FRK: An R package for spatial and spatio-temporal prediction with large datasets. Journal of Statistical Software, 98(4), 1-48. doi:10.18637/jss.v098.i04.
#'
#' Sainsbury-Dale, M. and Zammit-Mangion, A. and Cressie, N. (2021) Modelling with Non-Gaussian Spatial and Spatio-Temporal Data using FRK, arXiv:2110.02507
#' @export
#' @examples
#' library("FRK")
#' library("sp")
#' ## Generate process and data
#' m <- 250                                                   # Sample size
#' zdf <- data.frame(x = runif(m), y= runif(m))               # Generate random locs
#' zdf$Y <- 3 + sin(7 * zdf$x) + cos(9 * zdf$y)               # Latent process
#' zdf$z <- rnorm(m, mean = zdf$Y)                            # Simulate data
#' coordinates(zdf) = ~x+y                                    # Turn into sp object
#'
#' ## Construct BAUs and basis functions
#' BAUs <- auto_BAUs(manifold = plane(), data = zdf,
#'                   nonconvex_hull = FALSE, cellsize = c(0.03, 0.03), type="grid")
#' BAUs$fs <- 1 # scalar fine-scale covariance matrix
#' basis <- auto_basis(manifold =  plane(), data = zdf, nres = 2)
#'
#' ## Construct the SRE model
#' S <- SRE(f = z ~ 1, list(zdf), basis = basis, BAUs = BAUs)
#'
#' ## Fit with 2 EM iterations so to take as little time as possible
#' S <- SRE.fit(S, n_EM = 2, tol = 0.01, print_lik = TRUE)
#'
#' ## Check fit info, final log-likelihood, and estimated regression coefficients
#' info_fit(S)
#' logLik(S)
#' coef(S)
#'
#' ## Predict over BAUs
#' pred <- predict(S)
#'
#' ## Plot
#' \dontrun{
#' plotlist <- plot(S, pred)
#' ggpubr::ggarrange(plotlist = plotlist, nrow = 1, align = "hv", legend = "top")}
SRE <- function(f, data,basis,BAUs, est_error = TRUE, average_in_BAU = TRUE,
                sum_variables = NULL,
                normalise_wts = TRUE,
                fs_model = "ind", vgm_model = NULL,
                K_type = c("block-exponential", "precision", "unstructured"),
                normalise_basis = TRUE,
                response = c("gaussian", "poisson", "gamma",
                             "inverse-gaussian", "negative-binomial", "binomial"),
                link = c("identity", "log", "sqrt", "logit", "probit", "cloglog", "inverse", "inverse-squared"),
                include_fs = TRUE, fs_by_spatial_BAU = FALSE,
                ...) {



  ## Strings that must be lower-case (this allows users to enter
  ## response = "Gaussian", for example, without causing issues)
  response  <- tolower(response)
  link      <- tolower(link)
  K_type    <- tolower(K_type)

  ## Allow partial matching of string arguments
  K_type   <- match.arg(K_type)
  response <- match.arg(response)
  link     <- match.arg(link)


  ## data should be a list: Instead of checking this condition in .check_args1()
  ## and throwing an error if it doesn't hold, here we just force it to be a list
  if (!is.list(data)) data <- list(data)

  ## When the response has a size parameter, restrict the incidence matrices
  ## (Cz and Cp) to represent simple sums only. (See Appendix B, last paragraph,
  ## of the FRK v2 paper).
  ## Note that this code has to go here, and not in .check_args1(), because we
  ## are altering some of the variables.
  ## Note also that we inform the user of these decisions only if some of the
  ## observations are associated with multiple BAUs (these terms are irrelevant
  ## otherwise): To check this requires the incidence matrix.
  if (response %in% c("binomial", "negative-binomial")) {

    normalise_wts <- FALSE
    BAUs$wts      <- 1

    if(any(sapply(data, function(x) is(x, "SpatialPoints")))) {

      ## NOTE: The users are informed of these decisions later, only if they are
      ## found to be relevant.

      ## Enforce average_in_BAU = TRUE for simplicity
      average_in_BAU <- TRUE
      ## sum_variables is a vector of variable names that are to be summed
      ## rather than averaged when average_in_BAU = TRUE
      sum_variables  <- c(all.vars(f)[1], "k_Z", sum_variables)
      ## Remove possible duplicates of all.vars(f)[1] and "k_Z"
      sum_variables  <- unique(sum_variables)
    }
  }

  ## Check that the arguments are OK
  .check_args1(f = f, data = data, basis = basis, BAUs = BAUs, est_error = est_error,
               response = response, link = link, K_type = K_type,
               fs_by_spatial_BAU = fs_by_spatial_BAU, normalise_wts = normalise_wts,
               sum_variables = sum_variables, average_in_BAU = average_in_BAU)

  ## The weights of the BAUs only really matter if the data are SpatialPolygons.
  ## However, we still need a 1 for all other kinds; only produce a warning if
  ## we have areal data.
  if (is.null(BAUs$wts)) {
    BAUs$wts <- 1
    if (any(sapply(data, function(x) is(x, "SpatialPolygons"))) &&
        !response %in% c("binomial", "negative-binomial")) # wts doesn't come into play for binomial or neg. binomial data, as it is forced to 1
      cat("SpatialPolygons were provided for the data support. No 'wts' field was found in the BAUs, so all BAUs are assumed to be of equal weight; if this is not the case, set the 'wts' field in the BAUs accordingly.\n")
  }

  ## Extract the dependent variable from the formula
  av_var <- all.vars(f)[1]

  ## Number of data objects
  ndata <- length(data)

  ## Initialise list of matrices (We construct one for every data object then concatenate)
  S <- Ve <- Vfs <- X <- Z <- Cmat <- k_Z <- list()

  ## Number of spatial BAUs and basis functions
  ns <- dim(BAUs)[1]
  n_basis_spatial <- if(is(basis,"TensorP_Basis")) nbasis(basis@Basis1) else nbasis(basis)

  ## Evaluate the basis functions over the BAUs. If we have fewer spatial BAUs
  ## than basis functions, then we average the basis functions over the BAUs
  ## using Monte Carlo integration with 1000 samples per BAU.
  ## Otherwise, evaluate the basis functions over the BAU centroids.
  S0 <- eval_basis(basis, if(ns < n_basis_spatial) BAUs else .polygons_to_points(BAUs))

  ## Normalise basis functions for the prior process to have constant variance. This was seen to pay dividends in
  ## latticekrig, however we only do it once initially
  if(normalise_basis) {
    if(opts_FRK$get("verbose")) cat("Normalising basis function evaluations at BAU level...\n")
    xx <- sqrt(rowSums((S0) * S0))                        # Find the standard deviation (assuming unit basis function weight)
    xx <- xx + 1*(xx == 0)                                # In the rare case all basis functions evaluate to zero don't do anything
    S0 <- S0 / (as.numeric(xx))                           # Normalise the S matrix
  }

  ## Find the distance matrix associated with the basis-function centroids
  D_basis <- BuildD(basis)

  ## For each data object
  for(i in 1:ndata) {

    ## If we are estimating measurement error
    if(est_error && response == "gaussian") {
      ## Algorithm for estimating measurement error in space-time objects still not implemented
      if(is(data[[i]],"ST"))
        stop("Estimation of error not yet implemented for spatio-temporal data")
      data[[i]]$std <- 0            # Set it to zero initially
      this_data <- data[[i]]        # Allocate current data object

      ## Now estimate the measurement error using variogram methods
      this_data <- .est_obs_error(this_data,variogram.formula=f,
                                  vgm_model = vgm_model)

      ## Allocate the measurement standard deviation from variogram analysis
      data[[i]]$std <- this_data$std
    } else if (est_error && response != "gaussian") {
        cat("Not estimating the measurement error variance since the response model is not Gaussian.\n")
    }


    ## The next step is to allocate all data (both point and polygon referenced) to BAUs.
    ## We can either average data points falling in the same
    ## BAU (average_in_BAU == TRUE) or not (average_in_BAU == FALSE).
    if(opts_FRK$get("verbose")) cat("Binning data...\n")

    data_proc <- map_data_to_BAUs(data[[i]],
                                  BAUs,
                                  average_in_BAU = average_in_BAU,
                                  sum_variables = sum_variables)

    ## The mapping can fail if not all data are covered by BAUs. Throw an error message if this is the case
    if(any(is.na(data_proc@data[av_var])))
      stop("NAs found when mapping data to BAUs. Do you have NAs in your data?
                 If not, are you sure all your data are covered by BAUs?")

    ## Extract information from the data using the .extract.from.formula internal function
    L <- .extract.from.formula(f,data=data_proc)
    X[[i]] <- as(L$X,"Matrix")                # covariate information
    Z[[i]] <- Matrix(L$y)                     # data values
    Ve[[i]] <- Diagonal(x=data_proc$std^2)    # measurement-error variance

    ## Construct the incidence matrix mapping data to BAUs. This just returns
    ## indices and values which then need to be assembled into a sparse matrix.
    C_idx <- BuildC(data_proc,BAUs)

    ## Construct the sparse incidence Matrix from above indices. This is the
    ## matrix C_Z in the vignette
    Cmat[[i]] <- sparseMatrix(i=C_idx$i_idx,
                              j=C_idx$j_idx,
                              x=C_idx$x_idx,
                              dims=c(length(data_proc),  # ensure dimensions of C are good
                                     length(BAUs)))

    ## Every data should be affected by at least one BAU. If this is not the case throw an error message.
    if(any(rowSums(Cmat[[i]])==0))
      stop("I have found difficulty in associating the data with the BAUs.
                 If you have point-referenced data
                 then this could be because you have data outside BAUs. If you have
                 polygon data, then this could be because no BAUs centroids are
                 within the polygons. For polygon data, influence on a BAU is determined from
                 whether the BAU centroid falls within the polygon or not.")

    ## We ensure the polygon observations are a weighted average over the
    ## BAUs. This just means dividing each row by its row sum, so that each
    ## entry is between 0 and 1, and the row sums are all equal to 1.
    if (normalise_wts) Cmat[[i]] <- Cmat[[i]] / rowSums(Cmat[[i]])

    ## Only the independent model is allowed for now, future implementation will include CAR/ICAR (in development)
    if(fs_model == "ind") {
      Vfs[[i]] <- tcrossprod(Cmat[[i]] %*% Diagonal(x=sqrt(BAUs$fs)))
    } else stop("No other fs-model implemented yet")

    ## S0 is the matrix S in the vignette. Here S is the matrix SZ in the vignette.
    S[[i]] <- Cmat[[i]] %*% S0

    ## Construct k_Z in the same way that Z is constructed
    if("k_Z" %in% names(data_proc@data)) k_Z[[i]] <- Matrix(data_proc$k_Z)
  }

  if(fs_model == "ind") {
    Qfs_BAUs <- Diagonal(x=1/BAUs$fs)
    Vfs_BAUs <- Diagonal(x=BAUs$fs)
  } else stop("No other fs-model implemented yet")

  ## Now concatenate the matrices obtained from all the observations together
  S    <- do.call("rbind",S)
  X    <- do.call("rbind",X)
  Cmat <- do.call("rbind",Cmat)
  Z    <- do.call("rbind",Z)
  Ve   <- do.call("bdiag",Ve)
  Vfs  <- do.call("bdiag",Vfs)

  ## Indices of observed BAUs
  obsidx <- .observed_BAUs_from_Cmat(Cmat)
  
  # Size parameter
  if(response %in% c("binomial", "negative-binomial")) {

    Cmat_dgT <- .as(Cmat, "dgTMatrix")
    num_obs_each_BAU <- table(Cmat_dgT@j)
    if (!all(num_obs_each_BAU == 1)) {
      cat("For binomial and negative-binomial data, FRK enforces average_in_BAU = TRUE, and the size parameter and response variable of observations falling into the same BAU are summed rather than averaged (i.e., the response and the size parameter are included in the argument 'sum_variables').\n")
    }

    num_BAUs_each_obs <- table(Cmat_dgT@i)

    if (!all(num_BAUs_each_obs == 1)) {

      ## At least one observation is associated with multiple BAUs.
      ## Note: k_Z is not used in this scenario; i.e., it does not need
      ## to be provided by the user.

      ## Inform the user of our restrictions to BAUs$wts and normalise_wts
      ## (These terms are only applicable in this scenario)
      cat("Since the response distribution is binomial or negative-binomial, FRK enforces the non-zero elements of the incidence matrices (Cz and Cp in the papers) to be 1 and normalise_wts = FALSE: This means that aggregation over the BAUs is a simple, unweighted sum.\n")

      ## We require the size parameter in a field of the BAUs:
      if (!("k_BAU" %in% names(BAUs@data))) {
        stop("If the response distribution is binomial or negative-binomial and some observations supports are associated with multiple BAUs (e.g., areal data), the size parameter must be provided in the BAUs objects in a field named 'k_BAU'.")
      } else {
        k_BAU_O <- BAUs$k_BAU[obsidx]
      }

      k_Z <- Cmat %*% BAUs$k_BAU # construct k_Z by aggregating over the BAUs

    } else {

      ## Each observation is associated with a single BAU only (e.g., we have
      ## point-referenced data, or areal data where the observation supports and
      ## BAUs coincide). In this case, the user need only provide the size
      ## parameters associated with the observation supports, k_Z: See Section
      ## 2.5 of the FRK v2 paper for details.


      ## Check that k_Z is provided
      if (!all(sapply(data, function(l) "k_Z" %in% names(l@data)))) {

        ## If k_Z was not provided,
        ## as a back up convenience for the user, in the special case that each
        ## observation is associated with a single BAU only AND
        ## each BAU is associated with at most one observation support,
        ## k_Z and k are essentially the same, and we can use k_BAU if it
        ## was provided:
        if (all(num_obs_each_BAU == 1) && "k_BAU" %in% names(BAUs@data)) {
          k_Z <- Cmat %*% BAUs$k_BAU # construct k_Z by aggregating over the BAUs
        } else {
          stop("For binomial or negative-binomial data where each observation is associated with a single BAU only, the known constant size parameter must be provided for each observation. Please provide this in the data object, in a field called 'k_Z'.")
        }

      } else {
        ## k_Z was provided, and here we construct it from the binning process:
        k_Z <- as.numeric(do.call("rbind", k_Z))
      }

      ## Note: it is OK if non-observed BAUs do not have a size parameter,
      ## because these non-observed BAUs are "zeroed" out by the incidence
      ## matrix, Cz.
      ## Note: The aggregation described above was already done during the data
      ## binning stage (i.e., k_Z is treated as a covariate and summed).
      ## Note: Below, we must first re-order the observation size parameters,
      ## so that element i of k_Z is associated with the same BAU as element i
      ## of k_BAU_O; Cmat_dgT@i contains the index of the observation associated
      ## with BAU Cmat_dgT@j
      k_BAU_O <- k_Z[Cmat_dgT@i + 1]
      BAUs@data[Cmat_dgT@j + 1, "k_BAU"] <- k_BAU_O
    }

    if (any(is.na(k_BAU_O)))
      stop("The size parameter is required at all observed BAUs")

  }

  ## Irrespective of whether it was provided by the user or infered from the
  ## BAU-level size parameters, from this point forward we have a value for
  ## k_Z, the size-parameters associated with the observations.

  ## If we are estimating a unique fine-scale variance at each spatial BAU,
  ## simply replicate the initialisation of sigma2fs ns times.
  ## Also check a few things that we couldn't check before computing the matrices.
  if (fs_by_spatial_BAU) {
    ## Check each spatial BAU is observed enough times for a unique fine-scale
    ## variance parameter to be associated with each spatial BAU
    fewest_obs <-  min(table(obsidx %% ns)) # count of observations from spatial BAU with fewest observations
    if(fewest_obs == 0) {
      stop("A unique fine-scale variance at each spatial BAU can only be fit (i.e., fs_by_spatial_BAU = TRUE) if all spatial BAUs are observed, which is not the case for the provided data and BAUs. Please set fs_by_spatial_BAU = FALSE.")
    } else if(fewest_obs < 10) {
      warning(paste0("The smallest number of observations associated with a spatial BAUs is: ", fewest_obs,
                     ". As you have selected to fit a unique fine-scale variance at each spatial BAU (i.e., fs_by_spatial_BAU = TRUE), please consider if this is a sufficient number of observations."))
    }

    ## Throw a warning if the number of spatial BAUs (and hence number
    ## of fine-scale variance parameters) is very large
    if(ns > 500) warning(paste0("The number of spatial BAUs is relatively large (",ns,"). As you have chosen to fit a separate fine-scale variance parameter at each spatial BAU, there will be ",ns," fine-scale variance parameters to estimate, which may result in difficulties in model fitting. (However, it may not be an issue.)"))
  }

  ## If the response should be an integer, round to be safe
  ## (otherwise factorials will not make sense).
  if (response %in% c("poisson", "binomial", "negative-binomial"))
    Z <- round(Z)

  ## Information on fitting
  info_fit <- list("SRE not fitted to data yet. Please use SRE.fit()")

  ## Initialise the fixed effects and parameters for the EM algorithm.
  ## Note that this is done here for backwards compatability, and that
  ## the initalisations are not computationally demanding, so it is not a major
  ## issue that we always do them irrespective of which method is used.
  ## The initialisations for method = 'TMB' is more complicated, and is
  ## performed in SRE.fit(); doing it there means we know the method, and
  ## we do not have to create extra slots to pass the initialised parameters
  ## from SRE() to SRE.fit().
  l <- .EM_initialise(basis, Z, X, Ve, mstar = length(obsidx))

  ## Dummy values:
  if(!(response %in% c("binomial", "negative-binomial"))) k_BAU_O <- k_Z <- -1

  ## Construct the SRE object
  new("SRE",
      data=data,
      basis=basis,
      BAUs=BAUs,
      f = f,
      S = S,
      S0 = S0,
      D_basis = D_basis,
      Ve = Ve,
      Vfs = Vfs,
      Vfs_BAUs = Vfs_BAUs,
      Qfs_BAUs = Qfs_BAUs,
      Z = Z,
      Cmat = Cmat,
      X = X,
      mu_eta = l$mu_eta_init,
      S_eta = l$S_eta_init,
      Q_eta = l$Q_eta_init,
      K_type = K_type,
      Khat = l$K_init,
      Khat_inv = l$K_inv_init,
      alphahat = l$alphahat_init,
      sigma2fshat = l$sigma2fshat_init,
      fs_model = fs_model,
      info_fit = info_fit,
      response = response,
      link = link,
      mu_xi = l$mu_xi_init,
      k_Z = as.numeric(k_Z),
      k_BAU_O = as.numeric(k_BAU_O),
      include_fs = include_fs,
      normalise_wts = normalise_wts,
      fs_by_spatial_BAU = fs_by_spatial_BAU,
      obsidx = obsidx)
}

## Initialise the fixed effects and parameters for method = 'EM'
.EM_initialise <- function(basis, Z, X, Ve, mstar) {

  l    <- list()              # list of initial values
  nres <- max(basis@df$res)   # number of basis function resolutions

  ## Initialise the expectations and covariances from E-step to reasonable values
  l$mu_eta_init <- Matrix(0,nbasis(basis),1)
  l$mu_xi_init <- Matrix(0,mstar,1)
  l$S_eta_init <- Diagonal(x = rep(1,nbasis(basis)))
  l$Q_eta_init <- Diagonal(x = rep(1,nbasis(basis)))

  ## Start with reasonable parameter estimates (that will be updated in M-step)
  l$K_init = Diagonal(n=nbasis(basis),x = 1/(1/var(Z[,1])))
  l$K_inv_init = solve(l$K_init)

  if(!is.finite(determinant(t(X) %*% X)$modulus))
    stop("Matrix of covariates has columns that are linearly dependent. Please change formula or covariates.")

  if (ncol(X) == 0) {
    stop("We need at least one covariate in the model")
  } else {
    l$alphahat_init <- solve(t(X) %*% X) %*% t(X) %*% Z
  }
  l$sigma2fshat_init <- mean(diag(Ve)) / 4

  return(l)
}
