#ifndef _AMD_ORDER_WRAPPER_H_
#define _AMD_ORDER_WRAPPER_H_
#include <stddef.h>
#define Int int

void AMD_order_wrapper
    (
            Int *n,
            const Int *Ap,
            const Int *Ai,
            Int *P,
            double *Control,
            double *Info
    );

#endif
