#' @rdname eval_basis
#' @aliases eval_basis,Basis-matrix-method
setMethod("eval_basis",signature(basis="Basis",s="matrix"),function(basis,s,output = "matrix"){
        stopifnot(output %in% c("list","matrix"))
        x <- lapply(basis@fn,function(f) f(s))
        if (output=="matrix") {
            x <- Reduce("cbind",x)
        }
        x
    })


#' @rdname mass_matrix
#' @aliases mass_matrix,FEBasis-method
setMethod("mass_matrix",signature(B = "FEBasis"),
          function(B) {
            return(B@pars$M)
          })

#' @rdname stiffness_matrix
#' @aliases stiffness_matrix,FEBasis-method
setMethod("stiffness_matrix",signature(B = "FEBasis"),
          function(B) {
            return(B@pars$K)
          })

#' Load meshes
#'
#' Imports meshes from folders containing the files \code{p.csv, t.csv, M.csv} and \code{K.csv} which are csv files containing the vertices, triangulations, mass matrix and stiffness matrix details respectively. The matrices should be stored in the three-column format \code{i,j,k}.
#' @param paths a list of path names to the relevant folders
#' @export
#' @return a list of objects of class \code{FEBasis}. Please refer to the vignette for and example.
Load_meshes <- function(paths) {

  Meshes <- lapply(paths,function(x) {
    p <-  round(as.matrix(read.table(file.path(x,"p.csv"),sep=",")))
    tri <- as.matrix(read.table(file.path(x,"t.csv"),sep=","))
    M <-  as.matrix(read.table(file.path(x,"M.csv"),sep=","))
    K <-  as.matrix(read.table(file.path(x,"K.csv"),sep=","))
    n <- nrow(p)
    M <- sparseMatrix(i=M[,1],j=M[,2],x=M[,3],dims=c(n,n))
    K <- sparseMatrix(i=K[,1],j=K[,2],x=K[,3],dims=c(n,n))

    return(initFEbasis(p=p, t = tri, M = M, K = K))

  })
  names(Meshes) <- names(paths)

  return(Meshes)

}

# GRBF <- function(c,std,s) {
#   c_ext <- matrix(c,nrow(s),2,byrow=T)
#   dist_sq <- (rowSums((s - c_ext)^2))
#   return(exp(-0.5* dist_sq/(std^2) ))
# }



GRBF_wrapper <- function(manifold,mu,std) {
    stopifnot(is.matrix(mu))
    stopifnot(dimensions(manifold) == ncol(mu))
    stopifnot(is.numeric(std))
    stopifnot(std > 0)
    function(s) {
        stopifnot(ncol(s) == dimensions(manifold))
        dist_sq <- distance(manifold,s,mu)^2
        exp(-0.5* dist_sq/(std^2) )
    }
}

bisquare_wrapper <- function(manifold,c,R) {
    stopifnot(dimensions(manifold) == ncol(c))
    stopifnot(is.numeric(R))
    stopifnot(R > 0)
    function(s) {
        y <- distance(manifold,s,c)
        (1-(y/R)^2)^2 * (y < R)
    }
}


my_RBF <- function(r,mu=matrix(0,1,2),sigma2=1,A=1) {
  y <- A*exp(-r^2/(2*sigma2))
  return(y)
}


setMethod("basisinterp",signature(G = "Basis"),  # GRBF basis with mean offset as last weight
          function(G,s,weights=NULL) {
            y <- matrix(0,nrow(s),1)
            n <- length(weights)
            if (is.null(weights)) {
              weights <- rep(1,n)
            }
            for (i in 1:n) {
              y <- y + weights[i]*(G@fn[[i]](G@pars[[i]],s))
            }
            return(as.vector(y))
          })

tsearch2 <- function(x,y,t,xi,yi,bary=FALSE) {
  if (!is.vector(x)) {
    stop(paste(deparse(substitute(x)), "is not a vector"))
  }
  if (!is.vector(y)) {
    stop(paste(deparse(substitute(y)), "is not a vector"))
  }
  if (!is.matrix(t)) {
    stop(paste(deparse(substitute(t)), "is not a matrix"))
  }
  if (!is.vector(xi)) {
    stop(paste(deparse(substitute(xi)), "is not a vector"))
  }
  if (!is.vector(yi)) {
    stop(paste(deparse(substitute(yi)), "is not a vector"))
  }
  if (length(x) != length(y)) {
    stop(paste(deparse(substitute(x)), "is not same length as ",
               deparse(substitute(y))))
  }
  if (length(xi) != length(yi)) {
    stop(paste(deparse(substitute(xi)), "is not same length as ",
               deparse(substitute(yi))))
  }
  if (ncol(t) != 3) {
    stop(paste(deparse(substitute(t)), "does not have three columns"))
  }
  storage.mode(t) <- "integer"
  out <- .Call("tsearch", as.double(x), as.double(y), t, as.double(xi),
               as.double(yi), as.logical(bary),PACKAGE="geometry")
  if (bary) {
    names(out) <- c("idx", "p")
  }
  return(out)


}
