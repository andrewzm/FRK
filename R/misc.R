## #' Link and inverse-link function generator.
## #'
## #' Create the link function \eqn{g}  and inverse-link function \eqn{\psi},
## #' which link the conditional mean of the data \eqn{\mu} to the latent
## #' geostatistical process \eqn{Y} such that \eqn{g(\mu) = Y} and
## #' \eqn{\mu = \psi(Y)}, respectively.
## #'
## #' For families lacking a "known constant" parameter,
## #' \code{.link_fn} generates the functions linking
## #' the conditional mean of the data \eqn{\mu} to the latent
## #' geostatistical process \eqn{Y}. However, for families with a "known constant"
## #' parameter (such as the "binomial" and "negative-binomial"), \code{.link_fn} generate the functions linking
## #' the probability parameter \eqn{p} to the latent
## #' geostatistical process \eqn{Y}, and then the parameter \eqn{p} to the conditional mean \eqn{\mu}.
## #'
## #' @param kind A character string indicating which kind of link function is desired. 
## #' Valid values are 
## #' \describe{
## #'   \item{\code{"Y_to_mu"}}{Provides the function \eqn{\psi} such that \eqn{\mu = \psi(Y)}.}
## #'   \item{\code{"mu_to_Y"}}{Provides the function \eqn{g} such that \eqn{g(\mu) = Y}.}
## #'   \item{\code{"Y_to_prob"}}{Provides the function \eqn{\zeta} such that \eqn{p = \zeta(Y)}.}
## #'   \item{\code{"prob_to_Y"}}{Provides the function \eqn{h} such that \eqn{h(p) = Y}.}
## #'   \item{\code{"prob_to_mu"}}{Provides the function \eqn{\chi} such that \eqn{\mu = \chi(p)}.}
## #'   \item{\code{"mu_to_prob"}}{Provides the function \eqn{f} such that \eqn{f(\mu) = p}.}
## #' }
## #' Note that the latter four values are relevant only to the binomial and negative-binomial distributions with logit, probit, or cloglog link functions. 
## #' @param link A character string indicating the assumed link function. \emph{Not} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}. 
## #' @param response A character string indicating the assumed response distribution. \emph{Only} required if \code{kind} is \code{"prob_to_mu"} or \code{"mu_to_prob"}.
## #' @return A function.
.link_fn <- function (kind, link, response) {
  
  if (kind == "Y_to_mu") {
    if (link == "log")             psi <- function(Y) exp(Y)
    if (link == "identity")        psi <- function(Y) Y
    if (link == "logit")           psi <- function(Y) 1/(1 + exp(-Y))
    if (link == "probit")          psi <- function(Y) pnorm(Y)
    if (link == "cloglog")         psi <- function(Y) 1 - exp(-exp(Y))
    if (link == "inverse")         psi <- function(Y) 1/Y
    if (link == "inverse-squared") psi <- function(Y) 1/(sqrt(Y))
    if (link == "sqrt")     psi <- function(Y) Y^2
    return(psi)
    
  } else if (kind == "mu_to_Y") {
    if (link == "log")             g <- function(mu) log(mu)
    if (link == "identity")        g <- function(mu) mu
    if (link == "logit")           g <- function(mu) log( mu /(1 - mu))
    if (link == "probit")          g <- function(mu) qnorm(mu)
    if (link == "cloglog")         g <- function(mu) log(-log(1 - mu))
    if (link == "inverse")         g <- function(mu) 1/mu
    if (link == "inverse-squared") g <- function(mu) 1/(mu^2)
    if (link == "sqrt")     g <- function (mu) sqrt(mu)
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

# #' NASA colour palette
# #' @export
nasa_palette <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                  "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d","#f6da0c","#f5a009",
                  "#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

.constructS_O <- function(object) {
  obsidx <- observed_BAUs(object)
  return(drop0(object@S0[obsidx, , drop = FALSE]))
}

.constructC_O <- function(object) {
  obsidx <- observed_BAUs(object)
  return(object@Cmat[, obsidx, drop = FALSE])
}

.constructX_O <- function(object) {
  obsidx <- observed_BAUs(object)
  X_BAU <- as(.extract_BAU_X_matrix(object@f, object@BAUs), "matrix") # fixed-effect design matrix at BAU level
  return(X_BAU[obsidx, , drop = FALSE])
}




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


## #' Test equality within a given tolerance.
## #'
## #' Tests equality between objects \code{x} and \code{y} (within a given tolerance, \code{tol}).
## #' The primary purpose of this function is to avoid deeming objects to be unequal if they differ
## #' by some very tiny amount due to floating point inaccuracies.
## #' Of particular note is that the function can accept a matrix argument for \code{x} and a single numeric
## #' for \code{y}, and the output will be a matrix with elements 1 and 0 if elements of \code{x} are equal to
## #' \code{y} or not, respectively; i.e., it does elementwise comparisons.
## #'
## #' @param x \code{R} object.
## #' @param y \code{R} object we wish to compare to \code{x}.
## #' @param tol Tolerance.
## #' @return If \code{x} and \code{y} are single numbers, then the function
## #' returns 1 if \code{x} and \code{y} are equal (within \code{tol}), and 0
## #' otherwise. However matrices may also be passed, in which case the function
## #' returns a matrix of equal size with elementwise comparisons.
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



## #' Computation and concatenation of percentiles to a dataframe.
## #'
## #' Computes the percentiles or HPD interval bounds at each prediction location and appends 
## #' the result to \code{data}. Note that we use percentiles rather than quantiles
## #' because we including a "dot" (corresponding to the decimal place) in the 
## #' dataframe column name may cause issues. 
## #'
## #' @param data The dataframe we will append percentiles to; the number of rows of the matrices in \code{MC} and  in \code{data} must be equal
## #' @param MC List of matrices containing Monte Carlo samples
## #' @param percentiles a vector of scalars in [0, 100] specifying the desired percentiles; if \code{percentiles = NULL}, no percentiles are computed 
## #' @return The dataframe \code{data} with appended percentiles
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

.tapering_params <- function(D_matrices, taper) {
  
  ## Minimum distance between neighbouring basis functions.
  ## (add a large number to the diagonal, which would otherwise be 0)
  minDist <- sapply(D_matrices, function(D) min(D + 10^8 * diag(nrow(D))))

  return(taper * minDist)
}

## Covariance tapering based on distances.
##
## Computes the covariance tapering parameters \eqn{\beta} (which are dependent
## on resolution), number of non-zeros in each block of the tapered
## covariance matrix K_tap, and the tapered distance matrix whereby some distances
## have been set to zero post tapering (although the remaining non-zero distances 
## are unchanged).
##
## \code{taper} determines how strong the covariance tapering is; the ith taper parameter
## \eqn{\beta_i} is equal to \code{taper[i]} * \code{minDist[i]}, where
## \code{minDist[i]} is the minimum distance between basis functions at the
## \code{i}th resolution.
##
## @param D_matrices a list of distance matrices corresponding to each resolution of basis function.
## @param taper the strength of the taper (either a vector or a single number).
##
## @return A list containing:
## \describe{
##   \item{beta}{A vector of taper parameters.}
##   \item{nnz}{A vector containing the number of non-zeros in each block of the 
##   tapered prior covariance matrix, K_tap.}
##   \item{D_tap}{A sparse block-diagonal matrix containing the distances with 
##   some distances set to zero post tapering.}
## }
.cov_tap <- function(D_matrices, taper){

  if (class(D_matrices) != "list") 
    stop("D_matrices should be a list of matrices giving the distance between basis functions at each resolution")
  
  ## Taper parameters
  beta   <- .tapering_params(D_matrices = D_matrices, taper = taper)
  
  ## Construct D matrix with elements set to zero after tapering (non-zero entries
  ## are the untapered distances)
  D_tap <- .tapered_dist_matrices_nonzeroes_untouched(D_matrices, beta)
  
  ## Number of non-zeros in tapered covariance matrix at each resolution
  nnz <- sapply(D_tap, function(D) length(D@x))
  
  D_tap <- Matrix::bdiag(D_tap)

  return(list("beta" = beta, "nnz" = nnz, "D_tap" = D_tap))
}

.tapered_dist_matrices_nonzeroes_untouched <- function(D_matrices, beta) {
  
  D_tap <- lapply(1:length(D_matrices), function(i, D, beta) {
    D[[i]] <- as(D[[i]] * (D[[i]] < beta[i]), "sparseMatrix")
    
    ## Add explicit zeros to diagonal
    D <- D[[i]] + sparseMatrix(i = 1:nrow(D[[i]]), j = 1:nrow(D[[i]]), x = 0) 
  }, D = D_matrices, beta = beta)
  
  return(D_tap)
}

## Constructs \eqn{T_\beta}, the taper matrix, using the Spherical taper.
##
## Takes a list of distance matrices, and compute the spherical taper function
## for each matrix. Returns a list object. This is only used in the EM side; 
## for tapering within TMB, we compute the taper function within the C++ template.
## Denoting by \eqn{d(s,s^*)} the distance between two spatial locations,
## the Spherical taper is given by:
## \deqn{C_\beta(s, s^*) = (1-\frac{d(s,s^*)}{\beta})^2_{+}  (1+\frac{d(s,s^*)}{2\beta}), }
##  where \eqn{x_+ = max(x, 0)}.
## 
## @param beta The taper parameters (which vary based on resolution).
## @return A taper matrix, which is block-diagonal and of class \code{dgCmatrix}.
.T_beta_taper_matrix <- function(D_matrices, beta) {
  
  if (class(D_matrices) != "list") 
    stop("D_matrices should be a list of matrices giving the distance between basis functions at each resolution")
  
  spherical_taper <- function(D, beta) {
    apply(D, c(1,2), function(d){(1 + d / (2 * beta)) * (max(1 - d / beta, 0))^2})
  }
  
  T_beta <- mapply(spherical_taper, D = D_matrices, beta = beta, SIMPLIFY = FALSE)
  
  T_beta <- Matrix::bdiag(T_beta)
  
  ## Sanity check:
  # nnzero(T_beta) == nnzero(zapsmall(T_beta))
  
  return(T_beta)
}


## #' Neighbour matrix.
## #'
## #' Creates a matrix \eqn{A} with elements \eqn{A_{i, j}} equal to 1 if basis
## #' functions i and j (i not equal to j) are first order neighbours, 1/2 if they are second order neighbours, and so on.
## #' Neighbours with a larger order than specified by \code{order} have 0 in this matrix.
## #' The diagonal elements of \eqn{A} (i.e. \eqn{A_{i, i}}) indicate the row sums (if 
## #' the \code{order == 1}, then it is the totalnumber of first order neighbours associated with basis function i).
## #'
## #' This function is only designed for basis functions
## #' at \emph{one resolution}. It also assumes the basis functions are in a
## #' regularly-spaced lattice; the shape of the lattice is not important, however
## #' there \emph{must} be constant spacing between basis functions in a given
## #' direction (horizontal and vertical spacing can be different).
## #'
## #' @seealso \code{\link{.sparse_Q}}, \code{\link{.sparse_Q_block_diag}}
## #'
## #' @param df a dataframe containing the spatial coordinates
## #' @param order If order == 1, only first order neighbours are considered. If order == 2, second order neighbours are also considered, and so on
## #' @param diag_neighbours Indicates whether to consider the diagonal neighbours. If FALSE (default), only the horizontal and vertical neighbours are considered
## #' @return A "neighbour" matrix with element (i, j), for i not equal to j, equal to 1/l if basis functions i and j are lth order neighbours (provided \code{l <= order}), and 0 otherwise. Diagonal elements indicate the row sums
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



## #' Sparse precision matrix.
## #'
## #' Creates a sparse precision matrix \eqn{Q} with off diagonal elements equal
## #' to -1 if the basis functions are neighbours, and zero otherwise.
## #' The diagonal elements are equal to the number of neighbours for that basis
## #' function, plus some amount given by \code{kappa}.
## #'
## #' @param A "neighbour" matrix with element (i, j), for i not equal to j,
## #' equal to 1 if basis functions i and j are neighbours, and 0 otherwise,
## #' diagonal elements indicating the number of neighbours for that basis function.
## #' @param kappa Quantity to add to the diagonal elements. This must be positive if Q is to be positive definite.
## #' @param rho Quantity to multiply matrix by. This must be positive if Q is to be positive definite.
## #' @return A sparse precision matrix of class \code{dgCMatrix}.
## #' @seealso \code{\link{.neighbour_matrix}}, \code{\link{.sparse_Q_block_diag}}
.sparse_Q <- function(A, kappa, rho) {
  
  Q <- -A
  diag(Q) <- diag(A) + kappa
  Q <- rho * Q
  Q <- as(Q, "sparseMatrix")
  
  return(Q)
}



## #' Block-diagonal sparse precision matrix.
## #'
## #' Creates a block-diagonal sparse precision matrix, where the blocks are created
## #' using \code{sparse_Q}.
## #'
## #' @inheritParams .sparse_Q
## #' @inheritParams .neighbour_matrix
## #' @param df dataframe containing the spatial coordinates (named "loc1" and "loc2", etc.) and a column indicating the resolution of each basis function (named "res").
## #' @return list containing the sparse block-diagonal precision matrix (Q) of class "dgCMatrix", and the number of non-zero elements (nnz) at each resolution.
## #' @seealso \code{\link{.sparse_Q}}, \code{\link{.neighbour_matrix}}
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

