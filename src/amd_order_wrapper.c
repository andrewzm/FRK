#include "amd_internal.h"
#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  
GLOBAL void AMD_order_wrapper
(
    Int *n,
    const Int *Ap,
    const Int *Ai,
    Int *P,
    double *Control,
    double *Info
) {
amd_defaults (Control) ;		
amd_order(*n,Ap,Ai,P,Control,Info);
return(NULL);
}
