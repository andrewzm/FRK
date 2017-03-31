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
              {"FRK_distR_C", (DL_FUNC) &FRK_distR_C, 2},
              {NULL, NULL, 0}
            };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}

