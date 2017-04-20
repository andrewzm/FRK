
##########################################
######### NOT EXPORTED ###################
##########################################

#' @title Sparse Cholesky Factorisation with fill-in reducing permutations
#'
#' @noRd
#' @description This function is similar to chol(A,pivot=T) when A is a sparse matrix. By default, the fill-in reduction permutation is the approximate minimum degree permutation of Davis' SuiteSparse package configured to be slightly more aggressive than that in the Matrix package. When using the \code{R} Cholesyk decomposition, if the Cholesky factor fails because of lack of symmetry, the matrix is coerced to be symmetric using \code{forceSymmetric()}.
#' @param Q matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @param method If "amd", the SuiteSparse amd algorithm is used, if not that in the R Matrix package is employed
#' @return A list with two elements, Qpermchol (the permuted Cholesky factor) and P (the pivoting order matrix)
#' @keywords Cholesky factor
#' @examples
#' require(Matrix)
#' cholPermute(sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1)))
#' @references Davis T (2011). ``SPARSEINV: a MATLAB toolbox for computing the sparse inverse subset using the Takahashi equations.'' http://faculty.cse.tamu.edu/davis/suitesparse.html, Online: Last accessed 01 February 2016.
#' Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press, Boca Raton, FL.
cholPermute <- function(Q,method="amd")  {
  n <- nrow(Q)   # matrix dimension

  if(method == "amd") {  # is we will use the SuiteSparse amd
    P <- amd_Davis(Q)    # call the permutation algorithm
    Qp <- Q[P,P]         # permute the matrix
    Qpermchol  <- t(chol(Qp))                # do the Cholesky decomposition
    P <- sparseMatrix(i=P,j=1:n,x=1)         # construct the permutation matrix
    return(list(Qpermchol=Qpermchol,P=P))    # return both in a list

  } else {
    ## Try to do the Cholesky decomposition with R (usually this is done just for testing purposes)
    e <-tryCatch({ symchol <- Cholesky(Q)},
                 error= function(temp) {print("Cholesky failed, coercing to symmetric")},
                 finally="Cholesky successful")

    ## If an error was returned try to symmetrise it first using forceSymmetric().
    if (class(e) == "character")  {
      symchol <- Cholesky(forceSymmetric(Q))
    }

    j <- 1:n                               # column indicies
    i <- symchol@perm + 1                  # row indices
    P <- sparseMatrix(i,j,x=rep(1,n))      # construct permutation matrix
    if (class(e) == "character")  {
      Qpermchol <- t(chol(forceSymmetric(t(P)%*%Q%*%P)))   # do the Cholesky decomposition on the
    } else { Qpermchol <- t(chol(t(P)%*%Q%*%P)) }          # permuted matrix, possibly after a
    return(list(Qpermchol=Qpermchol,P=P))                  # forceSymmetric, then return list
  }

}

#' @title Solve the equation Qx = y
#'
#' @noRd
#' @description This function is similar to \code{solve(Q,y)} but with the added benefit that it allows for permuted matrices. This function does the job in order to minimise user error when attempting to re-permute the matrices prior to after solving. The user also has an option for the permuted Cholesky factorisation of Q to be carried out internally. This function is not exported.
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
  if (perm == F) {                                    # of there is no permutation
    if (dim(cholQ)[1] == 0) {                         # and the Cholesky is not already supplied
      e <-tryCatch({L <- t(chol(Q))},                 # try to do the Cholesky decomposition without permuting
                   error= function(temp) {            # possibly after attempting a forceSymmetric
                       print("Cholesky failed, coercing to symmetric")},
                   finally="Cholesky successful")
      if (class(e) == "character") {
        L <- t(chol(forceSymmetric(Q))) }
    }  else {
      L <- cholQ
    }

    v <- solve(L,y)                # standard solving of Ax=b using Cholesky
    x <- solve(t(L),v)
  }
  if (perm == T) {                 # If we wish to use permutations
    if (dim(cholQp)[1] == 0) {     # and the Cholesky factor was not supplied
      QP <- cholPermute(Q)         # permute and find the Cholesky
      Lp <- QP$Qpermchol           # Cholesky of permuted Q
      P <- QP$P                    # Permutation matrix
    } else {
      Lp <- cholQp                 # If supplied, just assign
    }

    v <- solve(Lp,t(P)%*%y)        # Standard solving for Ax=b under a permutation
    w <- solve(t(Lp),v)            # of the matrix A
    x <- P%*%w
  }
  return(x)
}

#' @title Solve the equation X = AQ^{-1}t(A) under permutations
#' @noRd
#' @description This function is a wrapper of solve() for finding \code{X = AQ^{-1}t(A)} when the permuted Cholesky factor of Q is known.
#' #'
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
#' cholsolveAQinvAT(A,X$Qpermchol,X$P)
cholsolveAQinvAT <- function(A,Lp,P) {
  ## Solve X = AQ^{-1}t(A) using the permuted Cholesky factor
  W <- t(solve(Lp,t(P)%*%t(A)))
  return(W %*% t(W))
}

#' @title Compute the Takahashi equations
#' @noRd
#' @description This function is wrapper for the Takahashi equations required to compute the marginal variances from the Cholesky factor of a precision matrix. The equations themselves are implemented in C using the SparseSuite package.
#' @param Q precision matrix (sparse or dense)
#' @param return_perm_chol if 1 returns the permuted Cholesky factor (not advisable for large systems)
#' @param cholQp the permuted Cholesky factor of Q (if known already)
#' @param P the pivot matrix (if known already)
#' @return If return_perm_chol == 0, returns the partial matrix inverse of Q, where the non-zero elements correspond to those in the Cholesky factor.
#' If !(return_perm_chol  == 0), returns a list with three elements, S (the partial matrix inverse), Lp (the Cholesky factor of the permuted matrix) and P (the permutation matrix).
#' @keywords Cholesky factor, linear solve
#' @examples
#' require(Matrix)
#' Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' X <- cholPermute(Q)
#' S_partial = Takahashi_Davis(Q,cholQp = X$Qpermchol,P=X$P)
#' @references Takahashi, K., Fagan, J., Chin, M.-S., 1973. Formation of a sparse bus impedance matrix and
#' its application to short circuit study. 8th PICA Conf. Proc.June 4--6, Minneapolis, Minn.
#' Davis T (2011). ``SPARSEINV: a MATLAB toolbox for computing the sparse inverse subset using the Takahashi equations.'' http://faculty.cse.tamu.edu/davis/suitesparse.html, Online: Last accessed 01 February 2016.
Takahashi_Davis <- function(Q,return_perm_chol = 0,cholQp = matrix(0,0,0),P=0) {

  n <- nrow(Q)  # matrix dimension

  if (dim(cholQp)[1] == 0) {                    # if permuted Cholesky factor not supplied
    symchol <- Cholesky(forceSymmetric(Q))      # find symbolic Cholesky decomposition
    j <- 1:n                                    # column indices
    i <- symchol@perm + 1                       # row indices
    P <- sparseMatrix(i,j,x=rep(1,n))           # Permutation matrix
    Lperm <- L <- t(chol(t(P)%*%Q%*%P))         # Cholesky factor
  } else {
    L <- cholQp                                 # else just assign
    P <- P
  }
  rm(Q)                                         # we don't need Q anymore, remove it
  if (return_perm_chol == 0) rm(cholQp)         # we also don't need the factor if it was supplied


  ## The following commands are adapted from sparseinv.m
  ## See https://au.mathworks.com/matlabcentral/fileexchange/33966-sparseinv-sparse-inverse-subset/content/sparseinv/sparseinv.m
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

## The AMD algrithm in the SuiteSparse package
amd_Davis <- function(Q) {
  n <- nrow(Q)  # matrix dimension
  Ap <- Q@p     # indices of compressed format matrix
  Ai <- Q@i

  ## Call the C functions in the SuiteSparse library
  X <- .C("AMD_order_wrapper",as.integer(n),as.integer(Ap),as.integer(Ai),
          P = integer(n), Control=double(5),Info=double(20))
  return(X$P + 1)
}

## Simple test function to ensure AMD works as it should
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
  #write.table(data.frame(i=Q@i,j=Q@j,x=1),file="Chol_test.csv")
  X <- .C("AMD_order_wrapper",as.integer(n),as.integer(Ap),as.integer(Ai),
          P = integer(n), Control=double(5),Info=double(20))
}


## Wrapper for the sparse inverse function in the SuiteSparse package
sparseinv_wrapper <- function(L,d,U,Zpattern) {

  n <- nrow(L)  # number of columns
  Lp <- L@p     # matrix description of Cholesky factor in column-compressed format
  Li <- L@i
  Lx <- L@x

  Up <- U@p     # same as above -- in our case U = L
  Uj <- U@i
  Ux <- U@x

  ## Set up Zpattern matrix (see sparseinv.m for details)
  Zpatp <- Zpattern@p
  Zpati <- Zpattern@i
  znz = Zpatp [n+1]

  ## Call SuiteSparse package
  X <- .C("sparseinv",as.integer(n),as.integer(Lp),as.integer(Li),as.double(Lx),as.double(d),as.integer(Up),as.integer(Uj),as.double(Ux),as.integer(Zpatp),as.integer(Zpati),result = double(znz))

  ## Retrive result
  X <- X$result

  ## Remove other variables (this seemed to help for enormous systems)
  rm(U,L,Zpattern,Ux,Uj,Up,Lp,Li,Lx)

  ## Construct matrix from returned results
  Z <- sparseMatrix(p = Zpatp, i =Zpati, x = X,index1=F)

  return(Z)
}

## Return trace of matrix
tr <- function(X) {
    sum(diag(X))
}

## Efficient method of finding the diagonal of the product of two matrices
## If the matrix is symmetric we don't need to transpose and save some time
diag2 <- function(X,Y,symm=FALSE) {
    if(!symm) rowSums(X * t(Y)) else rowSums(X * Y)
}

## Compute the log determinant from a Cholesky factor L
logdet <- function (L)
{
    diagL <- diag(L)
    return(2 * sum(log(diagL)))
}

## quickBinds on columns
quickcBind <- function(L) {
  quickBind(L,"c")
}

## quickBinds on rows
quickrBind <- function(L) {
  quickBind(L,"r")
}

## Performs a quick binding of sparse matrices by extract the indices and rearranging. This code was adapted fro
## http://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors
## This function should probably be implemented in C at some point
quickBind <- function(L,rc = "c") {

  ## L list a list of sparseMatrices
  nzCount<-lapply(L, function(x) length(as(x,"dgTMatrix")@x));    # number off non-zeros in each matrix
  nz<-sum(do.call(rbind,nzCount));                                # total number of non-zeros
  r<-vector(mode="integer",length=nz);                            # row indices
  c<-vector(mode="integer",length=nz);                            # column indices
  v<-vector(mode="double",length=nz);                             # values to go in matrix
  ind <- 1                                                        # starting
  nc  <- 0                                                        # column number
  nr  <- 0                                                        # row number
  for(i in 1:length(L)){                                          # for each matrix
    tempMat <- as(L[[i]],"dgTMatrix")                             # convert to row-column storage format
    ln<-length(tempMat@x);                                        # number of nonzeros for this matrix
    if(ln>0){                                                     # if there is at least one non-zero
      if(rc == "c") {                                             # if column bind
        r[ind:(ind+ln-1)] <- tempMat@i + 1;                       # add to row indices
        c[ind:(ind+ln-1)] <- tempMat@j+ nc + 1                    # add to column indices
      } else if (rc == "r") {                                     # if row bind
        r[ind:(ind+ln-1)] <- tempMat@i + nr + 1;                  # add to row indices
        c[ind:(ind+ln-1)] <- tempMat@j + 1                        # add to column indices
      }
      v[ind:(ind+ln-1)] <- tempMat@x                              # add to final matrix values
      ind<-ind+ln;                                                # update "starting index"
    }

    ## Adjust number of rows and columns so far in matrix
    if(rc == "c") {
      nc <- nc + ncol(tempMat)
      nr <- nrow(tempMat)
    } else if (rc == "r") {
      nr <- nr + nrow(tempMat)
      nc <- ncol(tempMat)
    }
  }

  ## Return final sparse matrix
  return (sparseMatrix(i=r,j=c,x=v,dims = c(nr,nc)));
}

## Given a matrix X returns Y such that Y[idx,idx] = X
reverse_permute <- function(X,idx) {
  X <- as(X,"dgTMatrix")
  dict <- data.frame(from = 1:length(idx),to = idx)
  i_idx <- data.frame(from = X@i+1) %>% left_join(dict,by="from")
  j_idx <- data.frame(from = X@j+1) %>% left_join(dict,by="from")
  sparseMatrix(i=i_idx$to, j=j_idx$to, x = X@x,dims = dim(X))
}
