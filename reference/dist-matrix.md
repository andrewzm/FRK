# Distance Matrix Computation from Two Matrices

This function extends `dist` to accept two arguments.

## Usage

``` r
distR(x1, x2 = NULL)
```

## Arguments

- x1:

  matrix of size N1 x n

- x2:

  matrix of size N2 x n

## Value

Matrix of size N1 x N2

## Details

Computes the distances between the coordinates in `x1` and the
coordinates in `x2`. The matrices `x1` and `x2` do not need to have the
same number of rows, but need to have the same number of columns (e.g.,
manifold dimensions).

## Examples

``` r
A <- matrix(rnorm(50),5,10)
D <- distR(A,A[-3,])
```
