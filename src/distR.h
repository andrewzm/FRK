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

#ifndef _DISTR_H_
#define _DISTR_H_
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>

SEXP _FRK_distR_C(SEXP xSEXP, SEXP ySEXP);
SEXP overhead_cpp(SEXP a, SEXP b);

#endif
