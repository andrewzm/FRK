#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix distR_C(NumericMatrix x,NumericMatrix y) {
    int nrow_x = x.nrow(), ncol_x = x.ncol();
    int nrow_y = y.nrow();
    NumericMatrix out(nrow_x,nrow_y);

    for (int i = 0; i < nrow_x; i++) {
        for (int j = 0; j < nrow_y; j++) {
            for (int k = 0; k < ncol_x; k++) {
                out(i,j) = out(i,j) + pow(x(i,k) - y(j,k),2);
            }
            out(i,j) = sqrt(out(i,j));
        }
    }
    return out;
}
