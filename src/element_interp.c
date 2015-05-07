#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  

/* The 3x3 inverse was obtained from http://www.cquestions.com/2011/09/c-program-to-find-inverse-of-matrix.html */


void element_interp(int *n, int *t_num, int *tri1, int *tri2, int *tri3, double *p1, double *p2, double *locs1, double *locs2, double *i_ind1, double *i_ind2, double *i_ind3, double *j_ind1, double *j_ind2, double *j_ind3, double *z1, double *z2, double *z3)

{
	int i,k,l,t_num_i;
	int this_tri1,this_tri2,this_tri3;
     	double A[3][3];
	double X[3][3];
	double r[3];
	double determinant;
	int b[3];

	A[0][0] = 1;
	A[1][0] = 1;
	A[2][0] = 1;

	for(i=1; i<=*n; i++) {
	        t_num_i = t_num[i-1];
	        this_tri1 = tri1[t_num_i-1];
	        this_tri2 = tri2[t_num_i-1];
	        this_tri3 = tri3[t_num_i-1];
	  
      		A[0][1] = p1[this_tri1-1];
      		A[0][2]= p2[this_tri1-1];
      	  A[1][1] = p1[this_tri2-1];
      	  A[1][2]= p2[this_tri2-1];
      	  A[2][1]= p1[this_tri3-1];
        	A[2][2]= p2[this_tri3-1];

	        determinant = 0;
	        for(l=0;l<3;l++)
      			determinant = determinant + (A[0][l]*(A[1][(l+1)%3]*A[2][(l+2)%3] - A[1][(l+2)%3]*A[2][(l+1)%3]));

      		for(k=0;k<3;k++){
      		      for(l=0;l<3;l++)
      			X[l][k] = ((A[(k+1)%3][(l+1)%3] * A[(k+2)%3][(l+2)%3]) - (A[(k+1)%3][(l+2)%3]*A[(k+2)%3][(l+1)%3]))/determinant;
      	        }          
		
		  
    		   r[0] = X[0][0] + locs1[i-1]*X[1][0] + locs2[i-1]*X[2][0];
    		   r[1] = X[0][1] + locs1[i-1]*X[1][1] + locs2[i-1]*X[2][1];
    		   r[2] = X[0][2] + locs1[i-1]*X[1][2] + locs2[i-1]*X[2][2];
    
    		   b[0] = 1;
    		   b[1] = 0;
    		   b[2] = 0;
    		   i_ind1[i-1] = i;
    		   j_ind1[i-1] = this_tri1;
    		   z1[i-1] = r[0]*b[0] + r[1]*b[1] + r[2]*b[2];
    	          			
    
    
    		   b[0] = 0;
    		   b[1] = 1;
    		   i_ind2[i-1] = i;
    		   j_ind2[i-1] = this_tri2;
    		   z2[i-1] = r[0]*b[0] + r[1]*b[1] + r[2]*b[2];
    
    		   b[1] = 0;
    		   b[2] = 1;
    		   i_ind3[i-1] = i;
    		   j_ind3[i-1] = this_tri3;
    		   z3[i-1] = r[0]*b[0] + r[1]*b[1] + r[2]*b[2];
    	
     }

}


