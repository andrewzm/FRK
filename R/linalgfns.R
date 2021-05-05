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

##########################################
######### NOT EXPORTED ###################
##########################################

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

## quickbinds on columns
quickcbind <- function(L) {
  quickbind(L,"c")
}

## quickbinds on rows
quickrbind <- function(L) {
  quickbind(L,"r")
}

## Performs a quick binding of sparse matrices by extract the indices and rearranging. This code was adapted fro
## https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors
## This function should probably be implemented in C at some point
quickbind <- function(L, rc = "c") {

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
