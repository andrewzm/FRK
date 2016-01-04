/*-----------------------------------------------------------------------------
  File:   thin.c
   
  Implementation of polygon thinning routine(s).

  History:
    18 May 1999 [Paul Wessel]
      - original code copied from gshhs_dp.c to thin.c

    04 May 2004 [Nicholas Boers]
      - removed main() function
      - moved Douglas_Peucker_i() declaration to thin.h

    02 Jun 2004 [Nicholas Boers]
      - additional comments
      - minor changes to code to support UTM as well as LL

    17 Jun 2004 [Nicholas Boers]
      - updated comments

    14 Nov 2004 [Nicholas Boers]
      - updated for PBSINT type
  ---------------------------------------------------------------------------*/

/* original header comments */
/*      @(#)gshhs_dp.c  1.1  05/18/99
 *
 * gshhs_dp applies the Douglas-Peucker algorithm to simplify a line
 * segment given a tolerance.  The algorithm is based on the paper
 * Douglas, D. H., and T. K. Peucker, Algorithms for the reduction
 *   of the number of points required to represent a digitized line
 *   of its caricature, Can. Cartogr., 10, 112-122, 1973.
 * The impmementation of this algorightm has been kindly provided by
 * Dr. Gary J. Robinson, Environmental Systems Science Centre,
 * University of Reading, Reading, UK
 *
 * Paul Wessel, 18-MAY-1999
 * Version: 1.1 Added byte flipping
 *          1.2 Explicit binary read for DOS.  POSIX compliance
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "globals.h"
#include "thin.h"

#define sqr(x) ((x)*(x))
#define D2R (M_PI/180.0)
#define F (D2R * 0.5)
#define FALSE 0
#define VNULL (void *)NULL

void *get_memory (void *prev_addr, PBSINT n, size_t size);

/* Stack-based Douglas Peucker line simplification routine */
/* returned value is the number of output points, or -1 on error */
/* Kindly provided by  Dr. Gary J. Robinson,
 *      Environmental Systems Science Centre,
 *      University of Reading, Reading, UK
 */

/* returns -1 if insufficient memory */
int Douglas_Peucker_i (PBSINT *x_source, PBSINT *y_source, PBSINT n_source,
                       double band, PBSINT *index, short units)
{
  const short LL = 0;   /* "LL" units */
  PBSINT        n_stack, n_dest, start, end, i, sig;
  PBSINT        *sig_start, *sig_end;   /* indices of start&end of working section */

  double dev_sqr, max_dev_sqr, band_sqr;
  double  x12, y12, d12, x13, y13, d13, x23, y23, d23;

  double factor; /* conversion factor */

  /* check for simple cases */

  if ( n_source < 3 ) return(0);    /* one or two points */

  /* more complex case. initialise stack */

  sig_start = (PBSINT *) get_memory (VNULL, n_source, sizeof (PBSINT));
  sig_end   = (PBSINT *) get_memory (VNULL, n_source, sizeof (PBSINT));

  if (sig_start == NULL || sig_end == NULL) {
    if (sig_start != NULL) free(sig_start);
    if (sig_end != NULL) free(sig_end);

    return (-1);
  }

  /* unit specific initialization */
  if (units == LL)
    {
      /* micro-degrees >> degrees */
      factor = 1.0e-6;
      /* kilometers >> degrees */
      band *= 360.0 / (2.0 * M_PI * MEAN_RADIUS);
    }
  else /* UTM */
    {
      /* meters >> meters */
      factor = 1.0;
      /* kilometers >> meters */
      band *= 1.0e3;
    }

  band_sqr = sqr(band);
                
  n_dest = 0;

  sig_start[0] = 0;
  sig_end[0] = n_source-1;

  n_stack = 1;

  /* while the stack is not empty  ... */

  while ( n_stack > 0 )
    {
      /* ... pop the top-most entries off the stacks */

      start = sig_start[n_stack-1];
      end = sig_end[n_stack-1];

      n_stack--;

      if ( end - start > 1 )  /* any intermediate points ? */
        {
          /* ... yes, so find most deviant intermediate point to
             either side of line joining start & end points */

          /* micro-degrees >> decimal degrees */
          x12 = factor * (x_source[end] - x_source[start]);
          if ((units == LL) && fabs (x12) > 180.0) x12 = 360.0 - fabs (x12);
          y12 = factor * (y_source[end] - y_source[start]);
          /* scale x12 for latitude (same as in .initPlotRegion() )*/
          /* (mean y in micro-degrees >> decimal degrees >> radians) */
          if (units == LL)
            x12 *= cos (F * factor * (y_source[end] + y_source[start]));
          /* (degrees)^2 between 1 and 2 */
          d12 = sqr(x12) + sqr(y12);

          for ( i = start + 1, sig = start, max_dev_sqr = -1.0; i < end; i++ )
            {
              /* same as above, different points: intermediate and start */
              x13 = factor * (x_source[i] - x_source[start]);
              if ((units == LL) && fabs (x13) > 180.0) x13 = 360.0 - fabs (x13);
              y13 = factor * (y_source[i] - y_source[start]);
              if (units == LL)
                x13 *= cos (F * factor * (y_source[i] + y_source[start]));
              d13 = sqr(x13) + sqr(y13);

              /* same as above, different points: intermediate and end */
              x23 = factor * (x_source[i] - x_source[end]);
              if ((units == LL) && fabs (x23) > 180.0) x23 = 360.0 - fabs (x23);
              y23 = factor * (y_source[i] - y_source[end]);
              if (units == LL)
                x23 *= cos (F * factor * (y_source[i] + y_source[end]));
              d23 = sqr(x23) + sqr(y23);

              /* compare the distances -- since relative, all done based
                 on (degrees)^2 */
              if ( d13 >= ( d12 + d23 ) )
                dev_sqr = d23;
              else if ( d23 >= ( d12 + d13 ) )
                dev_sqr = d13;
              else
                dev_sqr =  sqr( x13 * y12 - y13 * x12 ) / d12;

              if ( dev_sqr > max_dev_sqr  )
                {
                  sig = i;
                  max_dev_sqr = dev_sqr;
                }
            }

          /* since just a comparison, compare the two squares rather than
             taking a square root */
          if ( max_dev_sqr < band_sqr )   /* is there a sig. intermediate point ? */
            {
              /* ... no, so transfer current start point */

              index[n_dest] = start;
              n_dest++;
            }
          else
            {
              /* ... yes, so push two sub-sections on stack for further processing */

              n_stack++;

              sig_start[n_stack-1] = sig;
              sig_end[n_stack-1] = end;

              n_stack++;

              sig_start[n_stack-1] = start;
              sig_end[n_stack-1] = sig;
            }
        }
      else
        {
          /* ... no intermediate points, so transfer current start point */

          index[n_dest] = start;
          n_dest++;
        }
    }


  /* transfer last point */

  index[n_dest] = n_source-1;
  n_dest++;

  free ((void *)sig_start);
  free ((void *)sig_end);
        
  return (n_dest);
}

void *get_memory (void *prev_addr, PBSINT n, size_t size)
{
  if (n == 0)
    return(VNULL); /* Take care of n = 0 */

  if (prev_addr)
    return (realloc ((void *) prev_addr, (size_t) (n * size)));
  else
    return (calloc ((size_t) n, (size_t) size));
}
