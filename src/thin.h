/*-----------------------------------------------------------------------------
  File:   thin.h
   
  Interface for polygon thinning routine(s).

  History:
    ?? ??? ????? [Paul Wessel]
      - original code copied from gshhs_dp.c to thin.h

    17 Jun 2004 [Nicholas Boers]
      - cleaned up comments
  ---------------------------------------------------------------------------*/

#include "globals.h"

/* Stack-based Douglas Peucker line simplification routine */
/* returned value is the number of output points */
/* Kindly provided by  Dr. Gary J. Robinson,
 *      Environmental Systems Science Centre,
 *      University of Reading, Reading, UK
 *      gazza@mail.nerc-essc.ac.uk
 */
/* Returns -1 if insufficient memory. */
int
Douglas_Peucker_i(PBSINT *x_source, PBSINT *y_source, PBSINT n_source, double band,
                  PBSINT *index, short units);
     /* *x_source:      Input coordinates in micro-degrees (LL) */
     /* *y_source:                     OR in meters (UTM)       */
     /* n_source:       Number of points                        */
     /* band:           tolerance in kilometers                 */
     /* *index:         output co-ordinates indeces             */
     /* units:          0 = "LL", 1 = "UTM"                     */

