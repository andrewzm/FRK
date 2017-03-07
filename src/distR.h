#ifndef _DISTR_H_
#define _DISTR_H_
#include <stddef.h>
#include <R.h>
#include <Rinternals.h>

SEXP FRK_distR_C(SEXP xSEXP, SEXP ySEXP);
SEXP overhead_cpp(SEXP a, SEXP b);

#endif
