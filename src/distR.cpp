/*FRK: An R Software package for spatial and spatio-temporal prediction
  with large datasets.
  Copyright (c) 2017 University of Wollongong
  Author: Andrew Zammit-Mangion, azm (at) uow.edu.au
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details. */

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
