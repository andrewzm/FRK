#' Link and inverse-link function generator.
#'
#' Create the link function \eqn{g}  and inverse-link function \eqn{\psi},
#' which link the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y} such that \eqn{g(\mu) = Y} and
#' \eqn{\mu = \psi(Y)}, respectively.
#'
#' For families lacking a "known constant" parameter,
#' \code{.link_fn} generates the functions linking
#' the conditional mean of the data \eqn{\mu} to the latent
#' geostatistical process \eqn{Y}. However, for families with a "known constant"
#' parameter (such as the "binomial" and "negative-binomial"), \code{.link_fn} generate the functions linking
#' the probability parameter \eqn{p} to the latent
#' geostatistical process \eqn{Y}, and then the parameter \eqn{p} to the conditional mean \eqn{\mu}.
#'
#' @param kind A character string indicating which kind of link function is desired. 
#' Valid values are 
#' \describe{
#'   \item{\code{"Y_to_mu"}}{Provides the function \eqn{\psi} such that \eqn{\mu = \psi(Y)}.}
#'   \item{\code{"mu_to_Y"}}{Provides the function \eqn{g} such that \eqn{g(\mu) = Y}.}
#'   \item{\code{"Y_to_prob"}}{Provides the function \eqn{\zeta} such that \eqn{p = \zeta(Y)}.}
#'   \item{\code{"prob_to_Y"}}{Provides the function \eqn{h} such that \eqn{h(p) = Y}.}
#'   \item{\code{"prob_to_mu"}}{Provides the function \eqn{\chi} such that \eqn{\mu = \chi(p)}.}
#'   \item{\code{"mu_to_prob"}}{Provides the function \eqn{f} such that \eqn{f(\mu) = p}.}
#' }
#' Note that the latter four values are relevant only to the binomial and negative-binomial distributions with logit, probit, or cloglog link functions. 
#' @param link A character string indicating the assumed link function. \emph{Not} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}. 
#' @param response A character string indicating the assumed response distribution. \emph{Only} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}.
#' @return A function.
.link_fn <- function (kind, link, response) {
  
  if (kind == "Y_to_mu") {
    if (link == "log")             psi <- function(Y) exp(Y)
    if (link == "identity")        psi <- function(Y) Y
    if (link == "logit")           psi <- function(Y) 1/(1 + exp(-Y))
    if (link == "probit")          psi <- function(Y) pnorm(Y)
    if (link == "cloglog")         psi <- function(Y) 1 - exp(-exp(Y))
    if (link == "inverse")         psi <- function(Y) 1/Y
    if (link == "inverse-squared") psi <- function(Y) 1/(sqrt(Y))
    if (link == "square-root")     psi <- function(Y) Y^2
    return(psi)
    
  } else if (kind == "mu_to_Y") {
    if (link == "log")             g <- function(mu) log(mu)
    if (link == "identity")        g <- function(mu) mu
    if (link == "logit")           g <- function(mu) log( mu /(1 - mu))
    if (link == "probit")          g <- function(mu) qnorm(mu)
    if (link == "cloglog")         g <- function(mu) log(-log(1 - mu))
    if (link == "inverse")         g <- function(mu) 1/mu
    if (link == "inverse-squared") g <- function(mu) 1/(mu^2)
    if (link == "square-root")     g <- function (mu) sqrt(mu)
    return(g)
    
  } else if (kind == "Y_to_prob") {
    if (link == "logit")      zeta <- function(Y) 1/(1 + exp(-Y))
    if (link == "probit")     zeta <- function(Y) pnorm(Y)
    if (link == "cloglog")    zeta <- function(Y) 1 - exp(-exp(Y))
    return(zeta)
    
  } else if (kind == "prob_to_Y") {
    if (link == "logit")      h <- function(p) log(p /(1 - p))
    if (link == "probit")     h <- function(p) qnorm(p)
    if (link == "cloglog")    h <- function(p) log(-log(1 - p))
    return(h)
    
  } else if (kind == "prob_to_mu") {
    if (response == "binomial")          chi <- function(p, k) k * p
    if (response == "negative-binomial") chi <- function(p, k) k * (1 / p - 1)
    return(chi)
    
  } else if (kind == "mu_to_prob") {
    if (response == "binomial")          f <- function(mu, k) mu / k
    if (response == "negative-binomial") f <- function(mu, k) k / (k + mu)
    return(f)
  }
  
  stop("Invalid arguments.")
}


## Exctracts the covariate design matrix of the fixed effects at the BAU level.
.extract_BAU_X_matrix <- function (formula, BAUs) {
  ## Retrieve the dependent variable name
  depname <- all.vars(formula)[1]
  
  ## Set the dependent variable in BAUs to something just so that 
  ## .extract.from.formula doesn't throw an error. We are not returning M, 
  ## so we don't need to worry about NULL-ing afterwards.
  BAUs[[depname]] <- 0.1
  
  ## Extract covariates from BAUs
  L <- .extract.from.formula(formula, data = BAUs)
  X <- as(L$X,"Matrix")
  
  return(X)
}



#' @rdname observed_BAUs
setMethod("observed_BAUs", signature(SRE_model = "SRE"), function (SRE_model) {
  
  ## Note that Cmat maps BAUs to the observations. The dimension of SRE_model@Cmat is
  ## (number of observations) * (number of BAUs).
  Cmat <- as(SRE_model@Cmat, "dgTMatrix")
  obsidx <- unique(Cmat@j) + 1 # unique() probably not necessary but want to be safe
  
  return(obsidx)
})


#' @rdname unobserved_BAUs
setMethod("unobserved_BAUs",signature(SRE_model = "SRE"), function (SRE_model) {
  
  ## Id of observed BAUs:
  obsidx <- observed_BAUs(SRE_model)
  
  ## Id of unobserved BAUs (ncol(SRE_model@Cmat) is the total number of BAUs):
  unobsidx <- (1:ncol(SRE_model@Cmat))[-obsidx]
  
  return(unobsidx)
})


## Determine if the basis functions are in a regular rectangular grid
## Inspired by: https://stackoverflow.com/a/32916893
.test_regular_grid <- function(x, y, rectangular) {
  
  a <- sort(unique(x))
  b <- sort(unique(y))
  
  condition <- .zero_range(diff(a)) && .zero_range(diff(b))
  
  if (rectangular)
    condition <- condition & ((length(a) * length(b)) == length(x))
  
  return(condition)
}



## Determine if range of vector is FP 0.
## (tests whether all elements of a vector are equal, with a tolerance) 
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}


#' Test equality within a given tolerance.
#'
#' Tests equality between objects \code{x} and \code{y} (within a given tolerance, \code{tol}).
#' The primary purpose of this function is to avoid deeming objects to be unequal if they differ
#' by some very tiny amount due to floating point inaccuracies.
#' Of particular note is that the function can accept a matrix argument for \code{x} and a single numeric
#' for \code{y}, and the output will be a matrix with elements 1 and 0 if elements of \code{x} are equal to
#' \code{y} or not, respectively; i.e., it does elementwise comparisons.
#'
#' @param x \code{R} object.
#' @param y \code{R} object we wish to compare to \code{x}.
#' @param tol Tolerance.
#' @return If \code{x} and \code{y} are single numbers, then the function
#' returns 1 if \code{x} and \code{y} are equal (within \code{tol}), and 0
#' otherwise. However matrices may also be passed, in which case the function
#' returns a matrix of equal size with elementwise comparisons.
## Examples:
## Note that .equal_within_tol is not exported, so it cannot be used in the 
## usual roxygen @examples format. We could use FRK:::.equal_within_tol, 
## but this produces a warning. For unexported functions, it is best to leave
## examples as regular comments.
# .equal_within_tol(2, 2 + 0.00000000000001)
# .equal_within_tol(2, 2 + 0.1)
# 
# A <- matrix(1:4, ncol = 2, byrow = TRUE)
# B <- matrix(c(1:3, 5), ncol = 2, byrow = TRUE)
# .equal_within_tol(A, B)
# .equal_within_tol(A, 3)
.equal_within_tol <- function(x, y, tol = 1e-8) {
  return(1 * (abs(x - y) < tol))
}



#' Computation and concatenation of percentiles to a dataframe.
#'
#' Computes the percentiles or HPD interval bounds at each prediction location and appends 
#' the result to \code{data}. Note that we use percentiles rather than quantiles
#' because we including a "dot" (corresponding to the decimal place) in the 
#' dataframe column name may cause issues. 
#'
#' @param data The dataframe we will append percentiles to; the number of rows of the matrices in \code{MC} and  in \code{data} must be equal
#' @param MC List of matrices containing Monte Carlo samples
#' @param percentiles a vector of scalars in [0, 100] specifying the desired percentiles; if \code{percentiles = NULL}, no percentiles are computed 
#' @return The dataframe \code{data} with appended percentiles
.concat_percentiles_to_df <- function (data, MC, percentiles) {
  
  if (!is.null(percentiles)) {
    QOI <- gsub("_.*", "", names(MC))
    
    for(i in seq_along(MC)) {
      Q           <- t(apply(as.matrix(MC[[i]]), 1, quantile, percentiles / 100))
      colnames(Q) <- paste(QOI[i], "percentile", as.character(percentiles), sep = "_")
      data        <- cbind(data, Q)
    }
  }
  
  return(data)
}


## Since we will use ggplot2 we will first convert our objects to data frames.
## The fortify command in ggplot2 was found to not work well with the sp polys.
## The following function converts a SpatialPolygonsDataFrame to a data.frame.
## sp_polys: Object of class SpatialPolygonsDataFrame
## vars: vector of characters identifying the field names to retain
.SpatialPolygonsDataFrame_to_df <- function (sp_polys, vars = names(sp_polys)) {
  ## Extract the names of the polygons
  polynames <- as.character(row.names(sp_polys))
  
  ## For each polygon do the following
  list_polys <- lapply(1:length(sp_polys), function(i) {
    ## Extract coordinates
    coords <- sp_polys@polygons[[i]]@Polygons[[1]]@coords
    row.names(coords) <- NULL
    coords <- data.frame(coords)
    
    ## Create data frame with coordinates and with "id" as the polygon name
    poldf <- cbind(coords, id = polynames[i], stringsAsFactors = FALSE)
    rownames(poldf) <- NULL
    poldf
  })
  
  ## Bind all the returned data frames together
  df_polys <- bind_rows(list_polys)
  
  ## Make sure the data frame and SpatialPolygonsDataFrame have the same "key" id var
  df_polys$id <- as.character(df_polys$id)
  sp_polys$id <- row.names(sp_polys)
  
  cnames <- coordnames(sp_polys)                     # coordinate names
  vars_no_coords <- vars[which(!vars %in% cnames)]   # all selected var names
  
  ## Now bind the variables to the data frame
  if (length(vars_no_coords) > 0)
    df_polys <- left_join(df_polys, sp_polys@data[c("id",
                                                    vars_no_coords)], by = "id")
  df_polys
}

# ---- Various matrix construction functions ----

#' Covariance tapering based on distances.
#'
#' Computes the covariance tapering parameters \eqn{\beta} (which are dependent
#' on resolution), number of non-zeros in each block of the tapered
#' covariance matrix K_tap, and the tapered distance matrix whereby some distances
#' have been set to zero post tapering (although the remaining non-zero distances 
#' are unchanged).
#'
#' \code{taper} determines how strong the covariance tapering is; the ith taper parameter
#' \eqn{\beta_i} is equal to \code{taper[i]} * \code{minDist[i]}, where
#' \code{minDist[i]} is the minimum distance between basis functions at the
#' \code{i}th resolution.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param taper The strength of the taper (either a vector or a single number).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{A vector of taper parameters.}
#'   \item{nnz}{A vector containing the number of non-zeros in each block of the 
#'   tapered prior covariance matrix, K_tap.}
#'   \item{D_tap}{A sparse block-diagonal matrix containing the distances with 
#'   some distances set to zero post tapering. }
#' }
#' @seealso \code{\link{.K_matrix}}
.cov_tap <- function(D_matrices, taper = 8){
  
  ri <- sapply(D_matrices, nrow) # No. of basis functions at each res
  nres <- length(ri)
  
  ## Minimum distance between neighbouring basis functions.
  ## (add a large number to the diagonal, which would otherwise be 0)
  ## FIXME: This approach is flawed. What if we have a very skinny rectangle for the domain of interest?
  minDist <- vector()
  for(i in 1:nres) minDist[i] <- min(D_matrices[[i]] + 10^8 * diag(ri[i]))
  
  ## Taper parameters
  beta   <- taper * minDist
  
  ## Construct D matrix with elements set to zero after tapering
  D_tap <- list()
  nnz <- c()
  for (i in 1:nres) {
    indices    <- D_matrices[[i]] < beta[i]
    D_tap[[i]] <- as(D_matrices[[i]] * indices, "sparseMatrix")
    
    # Add explicit zeros
    D_tap[[i]] <- D_tap[[i]] + sparseMatrix(i = 1:ri[i], j = 1:ri[i], x = 0) 
    
    # Number of non-zeros in tapered covariance matrix at each resolution
    nnz[i] <- length(D_tap[[i]]@x) # Note that this approach DOES count explicit zeros
  }
  
  D_tap <- Matrix::bdiag(D_tap)
  
  return(list("beta" = beta, "nnz" = nnz, "D_tap" = D_tap))
}



#' K matrix, the un-tapered prior covariance matrix.
#'
#' Construct the prior covariance matrix of the random effects (NOT-TAPERED).
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param sigma2 The variance components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @param tau The correlation components: a vector with length equal to the number
#' of resolutions i.e. \code{length(D_matrices)}.
#' @return The block-diagonal (class \code{dgCmatrix}) prior covariance matrix, K.
.K_matrix <- function(D_matrices, sigma2, tau){
  
  exp_cov_function <- function(dist_matrix, sigma2, tau) {
    sigma2 * exp(- dist_matrix / tau)
  }
  
  if (class(D_matrices) == "list") {
    
    ## Construct the individual blocks of the K-matrix
    K_matrices <- mapply(exp_cov_function,
                         dist_matrix = D_matrices,
                         sigma2 = sigma2, tau = tau,
                         SIMPLIFY = FALSE)
    
    ## Construct the full K-matrix
    K <- Matrix::bdiag(K_matrices)
    
  } else {
    K <- exp_cov_function(D_matrices, sigma2 = sigma2, tau  = tau)
  }
  
  return(K)
}

#' \eqn{K_\beta} matrix, the taper matrix.
#'
#' Constructs \eqn{K_\beta}, the taper matrix, using the Spherical taper.
#'
#' Denoting by \eqn{d(s,s^*)} the distance between two spatial locations,
#' the Spherical taper is given by:
#' \deqn{C_\beta(s, s^*) = (1-\frac{d(s,s^*)}{\beta})^2_{+}  (1+\frac{d(s,s^*)}{2\beta}), }
#' where \eqn{x_+ = max(x, 0)}.
#'
#' @param D_matrices A list of distance matrices corresponding to each resolution.
#' @param beta The taper parameters (which vary based on resolution).
#' @return A taper matrix, which is block-diagonal and of class \code{dgCmatrix}.
.K_beta_matrix <- function(D_matrices, beta){
  
  spherical_taper <- function(matrix, beta) {
    apply(matrix, c(1,2), function(h){(1 + h / (2 * beta)) * (max(1 - h / beta, 0))^2})
  }
  
  if (class(D_matrices) == "list") {
    K_beta <- mapply(spherical_taper, matrix = D_matrices, beta = beta, SIMPLIFY = FALSE)
  } else {
    K_beta <- spherical_taper(matrix = D_matrices, beta = beta)
  }
  
  K_beta <- Matrix::bdiag(K_beta)
  
  return(K_beta)
}

#' K_tap, the tapered prior covariance matrix.
#'
#' Construct K_tap, the \emph{tapered} prior covariance matrix, using the Spherical taper.
#'
#' @inheritParams .K_matrix
#' @inheritParams .K_beta_matrix
#' @return The \emph{tapered} prior covariance matrix, which is
#' block-diagonal and of class \code{dgCmatrix}.
.K_tap_matrix <- function(D_matrices, beta, sigma2, tau){
  
  K <- .K_matrix(D_matrices, sigma2, tau)
  K_beta <- .K_beta_matrix(D_matrices, beta)
  
  K_tap <- K * K_beta
  
  return(K_tap)
}


#' Neighbour matrix.
#'
#' Creates a matrix \eqn{A} with elements \eqn{A_{i, j}} equal to 1 if basis
#' functions i and j (i not equal to j) are first order neighbours, 1/2 if they are second order neighbours, and so on.
#' Neighbours with a larger order than specified by \code{order} have 0 in this matrix.
#' The diagonal elements of \eqn{A} (i.e. \eqn{A_{i, i}}) indicate the row sums (if 
#' the \code{order == 1}, then it is the totalnumber of first order neighbours associated with basis function i).
#'
#' This function is only designed for basis functions
#' at \emph{one resolution}. It also assumes the basis functions are in a
#' regularly-spaced lattice; the shape of the lattice is not important, however
#' there \emph{must} be constant spacing between basis functions in a given
#' direction (horizontal and vertical spacing can be different).
#'
#' @seealso \code{\link{.sparse_Q}}, \code{\link{.sparse_Q_block_diag}}
#'
#' @param df a dataframe containing the spatial coordinates
#' @param order If order == 1, only first order neighbours are considered. If order == 2, second order neighbours are also considered, and so on
#' @param diag_neighbours Indicates whether to consider the diagonal neighbours. If FALSE (default), only the horizontal and vertical neighbours are considered
#' @return A "neighbour" matrix with element (i, j), for i not equal to j, equal to 1/l if basis functions i and j are lth order neighbours (provided \code{l <= order}), and 0 otherwise. Diagonal elements indicate the row sums
.neighbour_matrix <- function(df, order = 1, diag_neighbours = FALSE) {
   
  A <- matrix(0, nrow = nrow(df), ncol = nrow(df))
  
  ## absolute difference in each dimension for each knot
  abs_diff <- lapply(df, function(x) abs(outer(x, x, "-")))
  
  ## Vectors containing all x and y distances. 
  ## Note that the first elements in each vector is 0. 
  ## Note also that we only use 1 row from the abs_diff matrices (this helps 
  ## to prevents problems with unique() and floating point accuracy and is a 
  ## bit faster)
  distances <- lapply(abs_diff, function(X) sort(unique(X[1, ])))

  ## TODO: I'd like to think about how to use irregular basis functions. 
  for (current_order in 1:order) { ## Order is experimental, and is hard-coded to be 1 for now
    
    ## Extract the smallest distance which is not zero, provided the basis functions 
    ## are not in a straight line. 
    ## If the basis functions have the same coordinate in a given dimension, the 
    ## corresponding entry in distances will contain only a single element (0), 
    ## so x[i + 1] doesn't work (results in NA). 
    min_d <- lapply(distances, function(x) if(length(x) == 1) 0 else x[current_order + 1])

    ## Find the neighbours. This is based on the idea that, provided the basis
    ## functions are regularly spaced, neighbours in a given dimension should be 
    ## separated by the minimum distance between basis functions of that dimension (condition1), 
    ## and have the same coordinates for the other dimensions (condition2). 
    d <- length(abs_diff) # number of spatial coordinates
    neighbours <- lapply(seq_along(abs_diff), function(i) {
      condition1 <- .equal_within_tol(abs_diff[[i]], min_d[[i]])
      if (d > 1) {
        condition2 <- lapply((1:d)[-i], function(j) .equal_within_tol(abs_diff[[j]], 0))
        condition2 <- Reduce("&", condition2) # Convert list of matrices to a single matrix
      } else {
        condition2 <- TRUE
      }
      return(condition1 & condition2)
    })
    
    ## Consider the diagonal neighbours, if specified.
    if (diag_neighbours == TRUE) {
      ## For basis-functions to be diagonal neighbours, they need to have the minimum
      ## (non-zero distance) between them in each spatial coordinate.
      diagonal_neighbours <- lapply(seq_along(abs_diff), function(i) {
        .equal_within_tol(abs_diff[[i]], min_d[[i]])
      })
     
      diagonal_neighbours <-  Reduce("&", diagonal_neighbours) # Convert list of matrices to a single matrix
      ## Update A
      A <- A + 1/current_order * diagonal_neighbours
    }
    ## Update neighbour matrix (with zeros along the diagonal)
    ## We weight the neighbours by their order.
    A <- A + 1/current_order * Reduce("+", neighbours)
  }
  
  ## Add the sums of each row to the diagonal (required for use in the 
  ## precision matrix computation later)
  diag(A) <- rowSums(A)
  
  return(A)
}



#' Sparse precision matrix.
#'
#' Creates a sparse precision matrix \eqn{Q} with off diagonal elements equal
#' to -1 if the basis functions are neighbours, and zero otherwise.
#' The diagonal elements are equal to the number of neighbours for that basis
#' function, plus some amount given by \code{kappa}.
#'
#' @param A "neighbour" matrix with element (i, j), for i not equal to j,
#' equal to 1 if basis functions i and j are neighbours, and 0 otherwise,
#' diagonal elements indicating the number of neighbours for that basis function.
#' @param kappa Quantity to add to the diagonal elements. This must be positive if Q is to be positive definite.
#' @param rho Quantity to multiply matrix by. This must be positive if Q is to be positive definite.
#' @return A sparse precision matrix of class \code{dgCMatrix}.
#' @seealso \code{\link{.neighbour_matrix}}, \code{\link{.sparse_Q_block_diag}}
.sparse_Q <- function(A, kappa, rho) {
  
  Q <- -A
  diag(Q) <- diag(A) + kappa
  Q <- rho * Q
  Q <- as(Q, "sparseMatrix")
  
  return(Q)
}



#' Block-diagonal sparse precision matrix.
#'
#' Creates a block-diagonal sparse precision matrix, where the blocks are created
#' using \code{sparse_Q}.
#'
#' @inheritParams .sparse_Q
#' @inheritParams .neighbour_matrix
#' @param df dataframe containing the spatial coordinates (named "loc1" and "loc2", etc.) and a column indicating the resolution of each basis function (named "res").
#' @return list containing the sparse block-diagonal precision matrix (Q) of class "dgCMatrix", and the number of non-zero elements (nnz) at each resolution.
#' @seealso \code{\link{.sparse_Q}}, \code{\link{.neighbour_matrix}}
.sparse_Q_block_diag <- function(df, kappa, rho, order = 1, diag_neighbours = FALSE) {
  
  if (!("res" %in% names(df)))
    stop("To construct the sparse precision matrix, the basis-function dataframe must contain a column named res, indicating the resolution of the corresponding basis function.")
  nres <- length(unique(df$res))
  if (length(kappa) == 1) kappa <- rep(kappa, nres)
  if (length(rho) == 1) rho <- rep(rho, nres)
  
  ## Find the location columns (should be called loc1, loc2, etc.)
  loc_idx <- grep("loc", names(df))
  
  ## Construct the blocks
  Q_matrices  <- list()
  nnz <- c()
  for (i in unique(df$res)) { 
    A_i <- .neighbour_matrix(df[df$res == i, loc_idx], order = order, diag_neighbours = diag_neighbours)
    Q_matrices[[i]] <- .sparse_Q(A = A_i,
                                 kappa = kappa[i],
                                 rho = rho[i])
    nnz[i] <- Matrix::nnzero(Q_matrices[[i]]) # note that nnzero does not count explicit zeros
  }
  
  ## Block diagonal
  Q <- Matrix::bdiag(Q_matrices)
  
  return(list(Q = Q, nnz = nnz))
}


# ---- Basis manipulation ----

remove_basis_outside_polygon <- function(basis_object, polygon) {
  sp_df <- basis_object@df
  coordinates(sp_df) <- ~ loc1 + loc2
  keep <- which(!is.na(over(sp_df, polygon)))
  
  basis_object@fn <- basis_object@fn[keep]
  basis_object@pars <- basis_object@pars[keep]
  basis_object@df <- basis_object@df[keep, ]
  basis_object@n <- length(keep)
  
  return(basis_object)
}


# ---- BAU manipulation ----


setMethod("reverse_spatial_coords",signature(BAUs="SpatialPixelsDataFrame"),function(BAUs) {
  
  ## Number of dimensions of the BAUs
  n_coord <- dimensions(BAUs)
  
  ## Reverse the slots (must do this for each slot):
  ## First, deal with the @data slot.
  ## It is possible that there is additional data columns present, 
  ## so we need to cater for this.
  ## Desired coordinate order and their indices in the data:
  new_coord_order <- rev(coordnames(BAUs))
  coord_idx <- match(new_coord_order, names(BAUs@data), nomatch = 0)
  
  ## Indices of the remaining columns
  remaining_idx <- which(!names(BAUs@data) %in% new_coord_order)
  
  ## Subset the data
  BAUs@data <- BAUs@data[, c(coord_idx, remaining_idx)]
  
  ## Reverse the remaining slots
  BAUs@coords <- BAUs@coords[, new_coord_order]
  BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
  
  ## Note that we cannot subset BAUs@grid (S4 method). 
  ## Instead, we will manipulate the slots directly.
  tmp <- BAUs@grid
  tmp@cellcentre.offset <- tmp@cellcentre.offset[n_coord:1]
  tmp@cellsize <- tmp@cellsize[n_coord:1]
  tmp@cells.dim <- tmp@cells.dim[n_coord:1]
  ## NB: note sure how to access the coordinate names of tmp directly as of now (GridTopology
  ## objects have only an assignment method for coordnames, not a retrieval method). 
  ## For now I will just use the coordinate names of the BAUs.
  ## It doesn't seem to matter anyway, though, as it adjusts automatically.
  coordnames(tmp) <- coordnames(BAUs)[n_coord:1]
  
  BAUs@grid <- tmp
  
  return(BAUs)
})


setMethod("reverse_spatial_coords",signature(BAUs="SpatialPointsDataFrame"),function(BAUs) {
  
  ## The documentation says that the @data slot may in fact contain the coordinates.
  ## So, we need to allow for this possibility.
  ## Desired coordinate order and their indices in the data:
  new_coord_order <- rev(coordnames(BAUs))
  coord_idx <- match(new_coord_order, names(BAUs@data), nomatch = 0)
  
  ## Indices of the remaining columns
  remaining_idx <- which(!names(BAUs@data) %in% new_coord_order)
  
  ## Subset the data
  BAUs@data <- BAUs@data[, c(coord_idx, remaining_idx), drop = FALSE]
  
  ## Reverse the order of the remaining slots:
  BAUs@coords <- BAUs@coords[, new_coord_order]
  BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
  
  return(BAUs)
})


setMethod("reverse_spatial_coords",signature(BAUs="SpatialPoints"),function(BAUs) {
  
  ## The documentation does not explicitly say that coords can or cannot contain data, 
  ## but for safety I will assume it can. 
  ## For the spatial points object, I will simply reverse the coordinates matrix. 
  ## Hence, if we have (lon, lat, z), the reversed order will be (z, lat, lon).
  ## Desired coordinate order and their indices in the data:
  new_coord_order <- rev(coordnames(BAUs))
  
  ## Subset the data
  BAUs@coords <- BAUs@coords[, new_coord_order, drop = FALSE]
  
  ## Reverse the order of the remaining slot:
  BAUs@bbox   <- BAUs@bbox[new_coord_order, ]
  
  return(BAUs)
})


setMethod("reverse_spatial_coords",signature(BAUs="SpatialPolygonsDataFrame"),function(BAUs) {
  
  
  ## First, the @data slot:
  ## (Note that both coordinate and data columns may be present - we will just reverse all)
  new_coord_order <- 
    BAUs@data <- BAUs@data[, rev(colnames(BAUs@data))]
  
  ## Second, the @bbox:
  BAUs@bbox   <- BAUs@bbox[rev(row.names(BAUs@bbox)), ]
  
  ## Third, the @polygons
  ## FIXME: Have to do this for each polygons list (use lapply or something)
  tmp <- BAUs@polygons$`1` 
  class(BAUs@polygons$`1`)
  ## Swap the order of the label point ($labpt) slot:
  tmp@labpt <- rev(tmp@labpt)
  ## Swap the order of the coordinates ($coords) slot:
  coordinates(tmp@Polygons)
  coordinates(tmp)
  ## FIXME: not sure how to access the coordinates of the polygons.
  ## Don't want to spend any more time on this until Andrew says its worthwhile.
  
  return(BAUs)
})


setMethod("reverse_spatial_coords",signature(BAUs="STFDF"),function(BAUs) {
  
  ## FIXME: should allow for SpatialPolygonsDF and SpatialPolygons too
  if(is(BAUs@sp, "SpatialPointsDataFrame") | 
     is(BAUs@sp, "SpatialPixelsDataFrame") | 
     is(BAUs@sp, "SpatialPoints"))
    BAUs@sp <- reverse_spatial_coords(BAUs@sp)
  else
    stop("The underlying spatial object should be of class 'SpatialPointsDataFrame' or 'SpatialPixelsDataFrame'.")
  
  return(BAUs)
})


setMethod("reverse_spatial_coords",signature(BAUs="STIDF"),function(BAUs) {
  
  ## FIXME: should allow for SpatialPolygonsDF and SpatialPolygons too
  if(is(BAUs@sp, "SpatialPointsDataFrame") | 
     is(BAUs@sp, "SpatialPixelsDataFrame") | 
     is(BAUs@sp, "SpatialPoints"))
    BAUs@sp <- reverse_spatial_coords(BAUs@sp)
  else
    stop("The underlying spatial object should be of class 'SpatialPointsDataFrame', 'SpatialPixelsDataFrame', or 'SpatialPoints'.")
  
  return(BAUs)
})


## Code common to the removal of spatial BAUs across methods
.remove_spatial_BAUs <- function(BAUs, rmidx, redefine_index = FALSE) {
  ntot <- nrow(BAUs@coords)
  if(!all(rmidx %in% 1:ntot))
    stop("Please ensure indices are numeric and within
             1 and the number of spatial BAUs.")
  
  BAUs_orig  <- BAUs
  BAUs <- BAUs[-rmidx, ]
  
  ## Check for indexing columns in @data slot of BAU object
  ## Note that the first check is for SpatialPoints and SpatialPolygons, each of which 
  ## do not have a @data slot.
  
  if ("data" %in% slotNames(class(BAUs_orig)) & redefine_index) {
    BAUs@data <- .redefine_indexing_variables(BAUs_orig@data, BAUs@data)
  }
  
  
  return(BAUs)
}


## Function to check if columns are "indexing variables".
## Returns the indices of these variables in the column names of df.
.redefine_indexing_variables <- function(data_orig, data_new) {
  
  ## Find the column indices of the indexing variables (if present)
  indexing_variables_idx <- apply(data_orig, 2, function(x) all(x == 1:nrow(data_orig))) %>%
    which()
  
  
  ## if indexing variables are present, redefine the corresponding column.
  if (length(indexing_variables_idx)) { 
    print(paste0("Possible indexing variables present, which will be redefined to maintain indexing following the removal of BAUs: ", names(indexing_variables_idx)))
    data_new[, indexing_variables_idx] <- 1:nrow(data_new)
  } 
  
  return(data_new)
}

setMethod("remove_BAUs",signature(BAUs="SpatialPoints"),function(BAUs, rmidx, redefine_index = FALSE) {
  return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})


setMethod("remove_BAUs",signature(BAUs="SpatialPointsDataFrame"),function(BAUs, rmidx, redefine_index = FALSE) {
  return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})

setMethod("remove_BAUs",signature(BAUs="SpatialPixelsDataFrame"),function(BAUs, rmidx, redefine_index = FALSE) {
  return(.remove_spatial_BAUs(BAUs, rmidx, redefine_index))
})

setMethod("remove_BAUs",signature(BAUs="STFDF"),function(BAUs, rmidx, redefine_index = FALSE) {
  
  BAUs_orig <- BAUs # Backup of BAUs for checks later
  
  ## Remove the spatial BAUs:
  BAUs@sp <- .remove_spatial_BAUs(BAUs@sp, rmidx)
  
  ## SUBSET @data to adjust for the new spatial BAUs:
  ## From the documentation, it says that @data is a data.frame containing
  ## measured values; space index cycling first, and time order preserved. 
  ## Hence, we should be able to remove the rows which correspond to rmidx
  ## by adding time to rmidx. 
  ## Also need to maintain space index cycling first. Achieve this by sorting the indcies.
  n_time     <- length(BAUs@time) # number of temporal frames
  n_spat     <- nrow(BAUs_orig@sp)      # ORIGINAL number of spatial BAUs
  increment   <- n_spat * 0:(n_time - 1) # Terms to add to each rmidx
  data_rmidx<- sort(c(outer(rmidx, increment, "+")))
  
  BAUs@data <- BAUs@data[-data_rmidx, ]
  
  ## Check for indexing columns in @data slot of BAU object
  ## Note that the first check is for SpatialPoints and SpatialPolygons, each of which 
  ## do not have a @data slot.
  if ("data" %in% slotNames(class(BAUs_orig)) & redefine_index)
    BAUs@data <- .redefine_indexing_variables(BAUs_orig@data, BAUs@data)
  
  return(BAUs)
})