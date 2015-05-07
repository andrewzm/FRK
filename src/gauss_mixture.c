/* Call this from R using 
   R CMD SHLIB gauss_mixture.c
   Once compiled load the library in R using dyn.load("gauss_mixture.so")
   Run as .C("gauss_mixture",...)
   Like this all the functions here are also in R */

#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define verbose 0



void gauss_mixture(double *x, double *w, double *mu, double *stdev, int *n, int *nxpts, double *result)
{
	int i,j;
	double ksum;
	for(i=0; i<*nxpts;i++) {
	   ksum =0;
	   for(j=0;j<*n;j++) {
		ksum += w[j]*dnorm(x[i],mu[j], stdev[j], 0);
	   }
	   result[i] = ksum;
	}
}


