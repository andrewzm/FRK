#' @title Pivoted Cholesky factorisation
#'
#' @description Pivoted Cholesky (same as pivot=T when carrying out Cholesky with dense matrices)
#' 
#' @param A matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @details This function should just be used to verify the pivot=T with dense matrices. Since it is an R implementation it is much slower than that in the base package.
#' @return A list with two elements, R (the Cholesky factor) and piv (the pivoting order)
#' @keywords Cholesky factor
#' @export
#' @examples
#' find_pivot_Bastos(matrix(c(0.1,0.2,0.2,1),2,2))
#' @references Leonardo S. Bastos and A. O'Hagan (2007). Diagnostics for Gaussian Process Emulators. \url{www.tonyohagan.co.uk/academic/pdf/diagtech.pdf}
find_pivot_Bastos <- function(A) {
 n <- nrow(A)
 R <- matrix(0,n,n)
 piv=1:n
 for (k in 1:(n-1)) {
  q <- which.max(diag(A[k:n,k:n])) + k - 1
  A <- swapcol(A,k,q)
  R <- swapcol(R,k,q)
  A <- swaprow(A,k,q)
  piv <- swap(piv,k,q)
  R[k,k] <- sqrt(A[k,k])
  R[k,((k+1):n)] <- (1/R[k,k]) * A[k,((k+1):n)]
  A[((k+1):n),((k+1):n)] <- A[((k+1):n),((k+1):n)] - outer(R[k,((k+1):n)],R[k,((k+1):n)])
 }
 R[n,n] <- sqrt(A[n,n])
 return(list(R=R,piv=piv))
}

#' @title Sparse Cholesky Factorisation with fill-in reducing permutations
#'
#' @description This function is similar to chol(A,pivot=T) when A is a sparse matrix. The fill-in reduction permutation is the approximate minimum degree permutation of 
#' Davis' SuiteSparse package configured to be slightly more aggressive than that in the Matrix package. If the Cholesky factor fails, the matrix is coerced to be symmetric.
#' 
#' @param Q matrix (sparse or dense), the Cholesky factor of which needs to be found
#' @param method If "amd", Timothy Davis SuiteSparse algorithm is used, if not that in the R Matrix package is employed
#' @param matlab_server A matlab server initiated using the R.matlab package by Henrik Bengtsson. Sparse Cholesky factorisation in MATLAB is generally much faster than in R.
#' @return A list with two elements, Qpermchol (the permuted Cholesky factor) and P (the pivoting order matrix)
#' @keywords Cholesky factor
#' @export
#' @examples 
#' require(Matrix)
#' cholPermute(sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1)))
#' @references Havard Rue and Leonhard Held (2005). Gaussian Markov Random Fields: Theory and Applications. Chapman & Hall/CRC Press
cholPermute <- function(Q,method="spam",matlab_server=NULL)  {
  n <- nrow(Q)
  
  if(method == "amd") {
    P <- amd_Davis(Q)
    Qp <- Q[P,P]
    if(!(is.null(matlab_server))) {
      Qpermchol  <- t(cholMATLAB(Qp,matlab_server))
    } else {
      Qpermchol  <- t(chol(Qp))
    }
    P <- sparseMatrix(i=P,j=1:n,x=1)
    return(list(Qpermchol=Qpermchol,P=P))
    
  } else if (method == "spam")   {
    return(.spam_chol(Q,amd=0))
    
  } else if (method == "spam_amd")   {
    return(.spam_chol(Q,amd=1))
    
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
#' @export
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
#'
#' @description This function is a wrapper of solve() for finding \code{X = AQ^{-1}t(A)} when the permuted Cholesky factor of Q is known.
#' #' 
#' @param Q ignored (deprecated)
#' @param A matrix
#' @param Lp Permuted Cholesky factor of Q
#' @param P the pivot matrix
#' @return x solution to \code{X = AQ^{-1}t(A)}
#' @keywords Cholesky factor, linear solve
#' @export
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
#'
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
#' @export
#' @examples 
#' require(Matrix)
#' Q = sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' X <- cholPermute(Q)
#' S_partial = Takahashi_Davis(Q,cholQp = X$Qpermchol,P=X$P)
#' @references Yogin E. Campbell and Timothy A Davis (1995). Computing the sparse inverse subset: an inverse multifrontal approach. \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.37.9276&rep=rep1&type=pdf}
Takahashi_Davis <- function(Q,return_perm_chol = 0,cholQp = matrix(0,0,0),P=0) {
  
  n <- nrow(Q)
  rm(Q)
  
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


#' @title Create an empty matrix
#'
#' @description Creates an empty sparse matrix of size 0 x 0
#' @export
#' @examples 
#' require(Matrix)
#' Q <- emptySp()
emptySp <- function() {
  as(matrix(0,0,0),"dgCMatrix")
}

#' @title Create a sparse identity matrix
#'
#' @description Creates a sparse identity matrix of size n x n
#' @param n size of matrix
#' @export
#' @examples 
#' require(Matrix)
#' Q <- Imat(4)
Imat <- function(n) {
  sparseMatrix(i=1:n, j=1:n,x=1)
}

#' @title Create an empty sparse matrix
#'
#' @description Creates an empty sparse matrix of size ni x nj
#' @param ni number of rows
#' @param nj number of columns. If NULL a square matrix is produced
#' @export
#' @examples 
#' require(Matrix)
#' Q <- Zeromat(2,5)
Zeromat <- function(ni,nj=NULL) {
  if(is.null(nj)) nj <- ni
   return(as(sparseMatrix(i={},j={},dims=c(ni,nj)),"dgCMatrix"))
}


#' @title Find the log determinant
#'
#' @description Find the log determinant of a matrix Q from its Cholesky factor L (which could be permutated or not)
#' @param L the Cholesky factor of Q
#' @examples 
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' logdet(chol(Q))
logdet <- function(L)  {
  ## Find the log-determinant of Q from its Cholesky L
  diagL <- diag(L)
  return(2*sum(log(diagL)))
}


#' @title Create a sparse diagonal matrix
#'
#' @description Creates a sparse diagonal matrix of length(xx) where xx is the vector containing the elements on the diagonal
#' @param xx diagonal vector
#' @export
#' @examples 
#' require(Matrix)
#' Q <- sparsediag(c(1,2,3,4))
sparsediag <- function(xx) {
  n <- length(xx)
  return(sparseMatrix(i=1:n,j=1:n,x=xx))
}


# Internal functions
swapcol <- function(A,p,q) {
 temp <- A[,p]
 A[,p] <- A[,q]
 A[,q] <- temp
 return(A)
}

swap <- function(v,p,q) {
  temp <- v[p]
  v[p] <- v[q]
  v[q] <- temp
  return(v)
}

swaprow <- function(A,p,q) {
  temp <- A[p,]
  A[p,] <- A[q,]
  A[q,] <- temp
  return(A)
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


# Find cholPermute using the spam package
.spam_chol <- function(Q,amd=T) {
  Qspam <- as.spam.dgCMatrix(Q)
  if(amd) {
    P <- linalg:::amd_Davis(Q)
    X  <- spam::chol(Qspam,pivot=P)
  } else {
    X  <- spam::chol(Qspam)
  }
  P <- sparseMatrix(i=X@pivot,j=1:nrow(X),x=1)
  Qpermchol <- as(as.dgCMatrix.spam(t(X)),"dtCMatrix")
  return(list(Qpermchol = Qpermchol,P=P))
}


cholMATLAB <- function(Q,matlab_server) {
  Q <- as(Q,"dgTMatrix")
  i = Q@i+1
  j = Q@j+1
  x = Q@x
  setVariable(matlab_server, i=i,j=j,x=x)
  cat("Doing Cholesky in MATLAB",sep="\n")
  evaluate(matlab_server, paste("Q = sparse(i,j,x);",
                                "L = (chol(Q));"))
  L <- getVariable(matlab_server,"L")$L
  L <- as(L,"dtCMatrix")
  cat("Finished Cholesky in MATLAB",sep="\n")
  evaluate(matlab_server, "clearvars Q L i j x")
  return(L)
  
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

# deprecated Takahashi trials
#----------------------------
Takahashidiag <- function(L) {
  ## Takahashi diag: Find marginal variance from Cholesky factor
  # For now convert to full matrices for indexing. We need to do this intelligently in the future for matrices > 10000x10000 which would fill up memory
  n <- dim(L)[1]
  X <- which(L!=0,arr.ind=T)
  i_ind <- X[,1]
  j_ind <- X[,2]
  Sigma <- Sigma + t(Sigma) - sparseMatrix(i=1:n,j=1:n,x=diag(Sigma))
  Sigma <- as(Sigma,"dgTMatrix") # coerce to i,j format
  Sigma[n,n] <- 1/L[n,n]^2
  numcomp = 0
  
  # Sigma <- as.matrix(Sigma)
  # L <- as.matrix(L)
  
  tic()
  for (i in seq(n-1,1,-1)) {
    Lii <- L[i,i]
    nz_indices <- intersect(i_ind[j_ind == i],(i+1):n)
    Lip1n_i<-L[nz_indices,i]
    for (j in intersect(seq(n,i,-1),which(abs(L[,i])>0))) {
      Sigma[i,j] <-  Sigma[j,i] <- (i==j)/(Lii^2) - 1/Lii*sum(Lip1n_i*Sigma[j,nz_indices ])
      numcomp = numcomp + 1
    }
  }
  toc()
  
  
  
  return(diag(Sigma))
}
Takahashidiag_Cseke <- function(L) {
  
  n <- dim(L)[1]
  invdiagL2 <- 1/diag(L)^2
  S <- L + t(L)
  S[S >0] = 1
  S[n,n] =invdiagL2[n]
  
  for (i in seq(n-1,1,-1)) {
    I   = i+which(abs(L[seq(i+1,n,1),i]) > 0)
    S[I,i] = -(S[I,I]%*%L[I,i])/L[i,i]
    S[i,I] = t(S[I,i])
    S[i,i] = invdiagL2[i] - (S[i,I]%*%L[I,i])/L[i,i];
    
  }
  return(diag(S))
}

Takahashi <- function(L, diag = TRUE, method = c("R", "C")) {
  #### R and C implementation of Takahashi diagonal equations by Jonty
  
  method <- match.arg(method)
  stopifnot(inherits(L, "CsparseMatrix")) # has @i and @p
  
  
  n <- nrow(L)
  ii <- L@i + 1 # in {1,...,n}
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp) # in {1,...,n}, non-decreasing
  N = length(ii)
  
  stopifnot(ii >= jj,              # lower triangular
            1:n %in% ii[ii == jj]) # full diagonal
  
  
  if (method == "C") {
    tic();
    dyn.load("Test.so")
    X <- .C("TakahashiC",as.integer(n),as.integer(N),as.integer(ii),as.integer(jj),as.double(L@x),results = double(N))
    toc();
    return(X$results[ii==jj])
  } else if (method == "R") {
    
    if (diag) { # this to speed up the calculation
      tic()
      S <- L; S@x[] <- -999
      
      S@x[ii == n & jj == n] <- 1 / L[n, n]^2
      
      if (n > 1)
        for (i in (n-1):1) {
          
          k <- ii[ii > i & jj == i]      # Find row numbers with non zero indices at this column
          if (length(k) == 0) {
            S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]
            
            js <- rev(jj[ii %in% k & jj >= i]) # going backwards
            #for (j in js) {
            #  skj <- S@x[ii == pmax(k, j) & jj == pmin(k, j)] # select from lower triangle
            #  S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            #}
            js = unique(js)
            js <- c(k,i)
            for(j in js) {
              skj <- apply(matrix(k),1,function(ind){ S@x[ii == max(ind,j) & jj == min(j,ind)] } )
              S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }
            
          }
        }
      toc()
      return(diag(S))
      
    } else { # diag = FALSE : use full-size S but only fill lower triangle
      
      S <- matrix(NA, n, n)
      
      S[n, n] <- 1 / L[n, n]^2
      
      if (n > 1)
        for (i in (n-1):1) {
          
          k <- ii[ii > i & jj == i]
          if (length(k) == 0) {
            S@x[ii == i & jj == i] <- 1 / L[i, i]^2
          } else {
            Lii <- L[i, i]
            Lki <- L[k, i]
            
            js <- n:i # going backwards
            
            for (j in js) {
              skj <- S[pmax(k, j), pmin(k, j)] # select from lower triangle
              S[j, i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
            }
          }
        }
      
      return(ifelse(is.na(S), t(S), S))
    }
    
  } else stop("Never get here!")
}
