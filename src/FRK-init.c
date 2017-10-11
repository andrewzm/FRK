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

#include "amd_order_wrapper.h"
#include "sparseinvR.h"
#include "distR.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


void attribute_visible R_init_FRK(DllInfo *info) {

  static R_NativePrimitiveArgType sparseinv_t[] = {
      INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,INTSXP,REALSXP
  };

  static R_NativePrimitiveArgType AMD_order_wrapper_t[] = {
        INTSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP
  };

  static R_CMethodDef cMethods[] = {
        {"sparseinv", (DL_FUNC) &sparseinv, 11, sparseinv_t},
        {"AMD_order_wrapper", (DL_FUNC) &AMD_order_wrapper, 6, AMD_order_wrapper_t},
        {NULL, NULL, 0}
  };


  static R_CallMethodDef callMethods[]  = {
              {"_FRK_distR_C", (DL_FUNC) &_FRK_distR_C, 2},
              {NULL, NULL, 0}
            };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}

