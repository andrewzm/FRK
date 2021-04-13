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

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifndef TMB_CALLDEFS
#define TMB_CALLDEFS                                            \
{"MakeADFunObject",     (DL_FUNC) &MakeADFunObject,     4},     \
{"InfoADFunObject",     (DL_FUNC) &InfoADFunObject,     1},     \
{"EvalADFunObject",     (DL_FUNC) &EvalADFunObject,     3},     \
{"MakeDoubleFunObject", (DL_FUNC) &MakeDoubleFunObject, 3},     \
{"EvalDoubleFunObject", (DL_FUNC) &EvalDoubleFunObject, 3},     \
{"getParameterOrder",   (DL_FUNC) &getParameterOrder,   3},     \
{"MakeADGradObject",    (DL_FUNC) &MakeADGradObject,    3},     \
{"MakeADHessObject2",   (DL_FUNC) &MakeADHessObject2,   4},     \
{"usingAtomics",        (DL_FUNC) &usingAtomics,        0},     \
{"TMBconfig",           (DL_FUNC) &TMBconfig,           2}
#endif


extern "C" {

   const static R_CallMethodDef callMethods[]  = {
     /*Currently not using C code for distances*/
     #if 0
        {"_FRK_distR_C", (DL_FUNC) &_FRK_distR_C, 2},
     #endif
     TMB_CALLDEFS,
     {NULL, NULL, 0}
   };

    void R_init_FRK(DllInfo *info) {

        R_registerRoutines(info, NULL, callMethods, NULL, NULL);
        R_useDynamicSymbols(info, FALSE);
        #ifdef TMB_CCALLABLES
        TMB_CCALLABLES("FRK");
        #endif
  }
}
