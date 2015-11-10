#' @title Sparse Cholesky Factorisation with fill-in reducing permutations
#'
#' @noRd
#' @description This function is similar to chol(A,pivot=T) when A is a sparse matrix. The fill-in reduction permutation is the approximate minimum degree permutation of
#' Davis' SuiteSparse package configured to be slightly more aggressive than that in the Matrix package. If the Cholesky factor fails, the matrix is coerced to be symmetric.
#'
#' @param Q matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @param method If "amd", Timothy Davis SuiteSparse algorithm is used, if not that in the R Matrix package is employed
#' @return A list with two elements, Qpermchol (the permuted Cholesky factor) and P (the pivoting order matrix)
#' @keywords Cholesky factor
#' @examples
#' require(Matrix)
#' cholPermute(sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1)))
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholPermute <- function(Q,method="amd")  {
  n <- nrow(Q)

  if(method == "amd") {
    P <- amd_Davis(Q)
    Qp <- Q[P,P]
    Qpermchol  <- t(chol(Qp))
    P <- sparseMatrix(i=P,j=1:n,x=1)
    return(list(Qpermchol=Qpermchol,P=P))

  } else {
    e <-tryCatch({ symchol <- Cholesky(Q)},error= function(temp) {print("Cholesky failed, coercing to symmetric")},finally="Cholesky successful")
    if (class(e) == "character")  {
      symchol <- Cholesky(forceSymmetric(Q))
    }


    j <- 1:n
    i <- symchol@perm + 1
    P <- sparseMatrix(i,j,x=rep(1,n))
    if (class(e) == "character")  {
      Qpermchol <- t(chol(forceSymmetric(t(P)%*%Q%*%P)))
    } else { Qpermchol <- t(chol(t(P)%*%Q%*%P)) }
    return(list(Qpermchol=Qpermchol,P=P))
  }

}

#' @title Solve the equation Qx = y
#'
#' @noRd
#' @description This function is similar to \code{solve(Q,y)} but with the added benefit that it allows for permuted matrices. This function does the job in order to minimise
#' user error when attempting to re-permute the matrices prior or after solving. The user also has an option for the permuted Cholesky factorisation of Q to be carried out
#' internally.
#'
#' @param Q matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @param y matrix with the same number of rows as Q
#' @param perm if F no permutation is carried out, if T permuted Cholesky factors are used
#' @param cholQ the Cholesky factor of Q (if known already)
#' @param cholQp the permuted Cholesky factor of Q (if known already)
#' @param P the pivot matrix (if known already)
#' @return x solution to Qx = y
#' @keywords Cholesky factor, linear solve
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' y = matrix(c(1,2),2,1)
#' cholsolve(Q,y)
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholsolve <- function(Q,y,perm=F,cholQ = matrix(1,0,0),cholQp = matrix(1,0,0),P=NA)  {
  ## Solve Qx = y
  if (perm == F) {
    if (dim(cholQ)[1] == 0) {
      e <-tryCatch({L <- t(chol(Q))},error= function(temp) {print("Cholesky failed, coercing to symmetric")},finally="Cholesky successful")
      if (class(e) == "character") {
        L <- t(chol(forceSymmetric(Q))) }
    }  else {
      L <- cholQ
    }

    v <- solve(L,y)
    x <- solve(t(L),v)
  }
  if (perm == T) {
    if (dim(cholQp)[1] == 0) {
      QP <- cholPermute(Q)
      Lp <- QP$Qpermchol
      P <- QP$P
    } else {
      Lp <- cholQp
    }

    v <- solve(Lp,t(P)%*%y)
    w <- solve(t(Lp),v)
    x <- P%*%w
  }
  return(x)
}

#' @title Solve the equation X = AQ^{-1}t(A) under permutations
#' @noRd
#' @description This function is a wrapper of solve() for finding \code{X = AQ^{-1}t(A)} when the permuted Cholesky factor of Q is known.
#' #'
#' @param Q ignored (deprecated)
#' @param A matrix
#' @param Lp Permuted Cholesky factor of Q
#' @param P the pivot matrix
#' @return x solution to \code{X = AQ^{-1}t(A)}
#' @keywords Cholesky factor, linear solve
#' @examples
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' X <- cholPermute(Q)
#' y <- matrix(c(1,2),2,1)
#' A <- y %*% t(y)
#' cholsolveAQinvAT(Q,A,X$Qpermchol,X$P)
cholsolveAQinvAT <- function(Q,A,Lp,P) {
  #Solve X = AQ^{-1}t(A)
  W <- t(solve(Lp,t(P)%*%t(A)))
  return(W %*% t(W))

}


#' @title Compute the Takahashi equations
#' @noRd
#' @description This function is wrapper for the Takahashi equations required to compute the marginal variances from the Cholesky factor of a precision matrix.
#' The equations themselves are implemented in C using the SparseSuite package of Timothy Davis.
#'
#' @param Q precision matrix (sparse or dense)
#' @param return_perm_chol if 1 returns the permuted Cholesky factor (not advisable for large systems)
#' @param cholQp the permuted Cholesky factor of Q (if known already)
#' @param P the pivot matrix (if known already)
#' @return if return_perm_chol == 0, returns the partial matrix inverse of Q, where the non-zero elements correspond to those in the Cholesky factor.
#' If !(return_perm_chol  == 0), returns a list with three elements, S (the partial matrix inverse), Lp (the Cholesky factor of the permuted matrix) and P (the
#' permutation matrix)
#' @keywords Cholesky factor, linear solve
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' X <- cholPermute(Q)
#' S_partial = Takahashi_Davis(Q,cholQp = X$Qpermchol,P=X$P)
#' @references Yogin E. Campbell and Timothy A Davis (1995). Computing the sparse inverse subset: an inverse multifrontal approach. \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.37.9276&rep=rep1&type=pdf}
Takahashi_Davis <- function(Q,return_perm_chol = 0,cholQp = matrix(0,0,0),P=0) {

  n <- nrow(Q)


  if (dim(cholQp)[1] == 0) {
    symchol <- Cholesky(forceSymmetric(Q))
    j <- 1:n
    i <- symchol@perm + 1
    P <- sparseMatrix(i,j,x=rep(1,n))
    Lperm <- L <- t(chol(t(P)%*%Q%*%P))
  } else {
    L <- cholQp
    P <- P
  }
  rm(Q)
  if (return_perm_chol == 0) rm(cholQp)

  d <- diag (L)
  L <- tril(L%*%sparseMatrix(i=1:n,j=1:n,x=1/d),-1)
  d <- d^2
  D <- sparseMatrix(i=1:n,j=1:n,x=d)

  #ii <- L@i + 1 # in {1,...,n}
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing

  gc()
  Zpattern <- sparseMatrix(c(L@i + 1,jj,1:n),c(jj,L@i + 1,1:n))
  rm(dp,jj)

  gc()
  Z <- sparseinv_wrapper(L,d,L,Zpattern)
  if (return_perm_chol == 0) {
    return(P%*%Z%*%t(P))
  } else {
    return(list(S=P%*%Z%*%t(P),Lp = cholQp,P=P)) # Only possible for small problems
  }

}


amd_Davis <- function(Q) {
  n <- nrow(Q)
  Ap <- Q@p
  Ai <- Q@i

  X <- .C("AMD_order_wrapper",as.integer(n),as.integer(Ap),as.integer(Ai),
          P = integer(n), Control=double(5),Info=double(20))
  return(X$P + 1)
}
amd_test <- function() {
  n=24
  Ap = c( 0, 9, 15, 21, 27, 33, 39, 48, 57, 61, 70, 76, 82, 88, 94, 100,
          106, 110, 119, 128, 137, 143, 152, 156, 160 )

  Ai = c(0, 5, 6, 12, 13, 17, 18, 19, 21,
         1, 8, 9, 13, 14, 17,
         2, 6, 11, 20, 21, 22,
         3, 7, 10, 15, 18, 19,
         4, 7, 9, 14, 15, 16,
         0, 5, 6, 12, 13, 17,
         0, 2, 5, 6, 11, 12, 19, 21, 23,
         3, 4, 7, 9, 14, 15, 16, 17, 18,
         1, 8, 9, 14,
         1, 4, 7, 8, 9, 13, 14, 17, 18,
         3, 10, 18, 19, 20, 21,
         2, 6, 11, 12, 21, 23,
         0, 5, 6, 11, 12, 23,
         0, 1, 5, 9, 13, 17,
         1, 4, 7, 8, 9, 14,
         3, 4, 7, 15, 16, 18,
         4, 7, 15, 16,
         0, 1, 5, 7, 9, 13, 17, 18, 19,
         0, 3, 7, 9, 10, 15, 17, 18, 19,
         0, 3, 6, 10, 17, 18, 19, 20, 21,
         2, 10, 19, 20, 21, 22,
         0, 2, 6, 10, 11, 19, 20, 21, 22,
         2, 20, 21, 22,
         6, 11, 12, 23 )
  Q <- as(sparseMatrix(i=Ai,p=Ap,index1=F,x=1),"dgTMatrix")
  write.table(data.frame(i=Q@i,j=Q@j,x=1),file="Chol_test.csv")
  X <- .C("AMD_order_wrapper",as.integer(n),as.integer(Ap),as.integer(Ai),
          P = integer(n), Control=double(5),Info=double(20))
}



sparseinv_wrapper <- function(L,d,U,Zpattern) {

  n <- nrow(L)
  Lp <- L@p
  Li <- L@i
  Lx <- L@x

  Up <- U@p
  Uj <- U@i
  Ux <- U@x

  Zpatp <- Zpattern@p
  Zpati <- Zpattern@i
  znz = Zpatp [n+1]


  X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))
  X <- X$result

  rm(U,L,Zpattern,Ux,Uj,Up,Lp,Li,Lx)
  Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X,index1=F)

  return(Z)
}


tr <- function(X) {
    sum(diag(X))
}

diag2 <- function(X,Y) {
    rowSums(X * t(Y))
}
