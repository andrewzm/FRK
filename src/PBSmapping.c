/*=============================================================================
  Copyright (C) 2003-2013 Fisheries and Oceans Canada

  This file is part of PBS Mapping.

  PBS Mapping is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PBS Mapping is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PBS Mapping; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
=============================================================================*/
/*-----------------------------------------------------------------------------
  File: PBSmapping.c

  Interface between R/S-PLUS and C code.

  Version:
    2.53

  Authors:
    Specifications: Jon Schnute, Rowan Haigh, Nicholas Boers
    Code:           Nicholas Boers

  To do:

  History:
    14 Jul 2003 [Nicholas Boers]
      - started the change log
    16 Jul 2003 [Nicholas Boers]
      - changed the interface to fixPOS() so that it will take
        floating-point POSs as input
      - moved fixBound() to R
    21 Jul 2003 [Nicholas Boers]
      - removed lines_intersect() from the file (not used)
      - moved calcPolygonArea to it's own files, polygons.c/.h
      - in time, polygon clipping should be moved to the polygons.c/.h
        files
    23 Jul 2003 [Nicholas Boers]
      - moved `pointInPolygon' to polygons.c/.h
      - renamed ulconv.c/.h to conversions.c/.h (in case we add more
        conversions)
      - creating a "floating.h" for floating point equality test, etc.
      - moved `clipPolygon' to polygons.c/.h
    27 Jul 2003 [Nicholas Boers]
      - rewrote closePolys()
    05 Aug 2003 [Nicholas Boers]
      - fixed epsilon
    14 Aug 2003 [Nicholas Boers]
      - fixed bug in calcArea() where polygons with multiple
        components were not being summed
    15 Aug 2003 [Nicholas Boers]
      - renamed functions:
        - ulConv() --> convUL()
        - locateEvents() --> findPolys()
        - clipPolys() --> clip()
    22 Dec 2003 [Nicholas Boers]
      - bug fix: previously, clip() assumed that if the first vertex
        in a polygon was not equal to 1, then the polygon must be a
        hole; now, it looks at the first two vertices to determine
        increasing or decreasing order
    04 Jan 2004 [Nicholas Boers]
      - cosmetic changes to the code
    04 Jan 2004 [Nicholas Boers]
       ## VERSION 1.00 ########################################################
    04 May 2004 [Nicholas Boers]
      - added thinPolys() function
    08 May 2004 [Nicholas Boers]
      - added clipPoly() function to clip a polygon against another polygon
    10 May 2004 [Nicholas Boers]
      - continued use of the Greiner/Hormann polygon clip routine would
        requires accounting for numerous special cases; switching to the
        gpc (General Polygon Clipper) library
      - fully implemented an interface to the gpc library (holes and all)
        - we will need to rename the R function... from clipPoly to something
          more intuitive
    11 May 2004 [Nicholas Boers]
      - renamed clipPoly() to clipToPoly()
    12 May 2004 [Nicholas Boers]
      - renamed clipToPoly() to joinPolys()
      - extensively cleaned up joinPolys() code
      - added `extern gpc_success` to catch errors that occur within GPC
      - discovered setjmp.h (setjmp() and longjmp())!!!  Now I can reliably
        detect memory allocation errors that occur in GPC!!!
    16 May 2004 [Nicholas Boers]
       ## VERSION 1.10 ########################################################
    20 May 2004 [Nicholas Boers]
      - added convexHull()
    24 May 2004 [Nicholas Boers]
      - fixed huge bug in joinPolys(), where an error would be signalled where
        one didn't really exist (check if all holes were processed was placed
        in the wrong spot...)
    26 May 2004 [Nicholas Boers]
      - fixed huge memory leaks in joinPolys()
      - added support for all R-accessible functions to interface with
        MemCheckDeluxe
    30 May 2004 [Nicholas Boers]
       ## VERSION 1.20 ########################################################
    02 Jun 2004 [Nicholas Boers]
      - updated thinPolys() to support UTM
    08 Jun 2004 [Nicholas Boers]
      - added thickenPolys() (only supports UTM)
    09 Jun 2004 [Nicholas Boers]
      - performed preliminary memory tests on all functions, and none appeared
        to have memory leaks
      - added calcOrientation() to determine whether vertices are clockwise
        or counter-clockwise
    10 Jun 2004 [Nicholas Boers]
      - gave fixPOS() the ability to update X/Y values so that it can make
        polygons that follow the GIS standards for clockwise/counter-clockwise
      - renamed calcPolygonArea() to calcPolyArea()    
    14 Jun 2004 [Nicholas Boers]
       ## VERSION 1.30 ########################################################
    15 Jun 2004 [Nicholas Boers]
      - renamed integrateHoles() to mergePolys(); now merges in all SIDs when
        specified as well as those representing holes
      - fixed huge memory leak and bug in joinPolys() -- the clip polygon
        wasn't handled correctly
    16 Jun 2004 [Nicholas Boers]
      - adjusted paths to #include's
    17 Jun 2004 [Nicholas Boers]
      - cleaning
      - added rollupPolys() to replace fixPOS() and integrateHoles() -- but it
        does much more than those two
    18 Jun 2004 [Nicholas Boers]
      - removed fixPOS() and integrateHoles()
    21 Jun 2004 [Nicholas Boers]
      - changed calcArea() to match new interface for calcPolyArea(); can now
        return negative areas for holes; substantially simplified function
    22 Jun 2004 [Nicholas Boers]
       ## VERSION 1.90 ########################################################
    28 Jun 2004 [Nicholas Boers]
      - added isIntersecting() to determine if a polygon is self-intersecting
    29 Jun 2004 [Nicholas Boers]
      - added isRetrace() to support isIntersecting()
      - added lineIntersect() to support isIntersecting()
    30 Jun 2004 [Nicholas Boers]
      - added isConvex()
      - added nPolyIntersects()
    05 Jul 2004 [Nicholas Boers]
      - renamed `convexHull()' to `calcConvexHull()'; it called
        calcConvexHull(), and renamed this call to convexHull()
    12 Jul 2004 [Nicholas Boers]
       ## VERSION 1.91 ########################################################
    18 Jul 2004 [Nicholas Boers]
       ## VERSION 1.92 ########################################################
    21 Jul 2004 [Nicholas Boers]
       ## VERSION 1.93 ########################################################
    28 Jul 2004 [Nicholas Boers]
      - rewrote joinPolys() and added several functions to support it
        (gpcOutputPoly(), gpcCreatePoly()); the code is much cleaner (and
        more powerful now)
    29 Jul 2004 [Nicholas Boers]
      - updated `joinPolys()' to use existing PIDs when either `polysA' or
        `polysB' contains only one "generic" polygon
       ## VERSION 1.99 ########################################################
    03 Aug 2004 [Nicholas Boers]
      - added note to isPolyConvex()
    16 Aug 2004 [Nicholas Boers]
      - moved isPolyConvex() to polygons.c
      - moved nPolyIntersects() to polygons.c
      - moved isRetrace() to polygons.c
      - moved linesIntersect() to polygons.c
    23 Aug 2004 [Nicholas Boers]
       ## VERSION 2.00 ########################################################
    14 Nov 2004 [Nicholas Boers]
      - replaced all instances of long with PBSINT to quickly support R's
        use of `int' rather than `long'
      - will make similar replacements in the dependent files
      - checked GPC functions to ensure appropriate types were being passed
        into the GPC, because I won't make similar replacements within that
        module
    For future changes, see the ChangeLog file.
  ---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <setjmp.h>

#include <R.h>
#include <Rdefines.h>

#include "globals.h"
#include "conversions.h"
#include "polygons.h"
#include "thin.h"
#include "floating.h"

#ifdef _MCD_CHECK
#include "mcd.h"
#endif /* _MCD_CHECK */

/* status options to pass back to R/S-PLUS */
#define PBS_SUCCESS     0       /* success */
#define PBS_ERR_MEM     1       /* insufficient memory */
#define PBS_ERR_OUT     2       /* output full */
#define PBS_ERR_OT1     3       /* other error */
#define PBS_ERR_OT2     4       /* other error */
#define PBS_ERR_OT3     5       /* other error */
#define PBS_ERR_OT4     6       /* other error */
#define PBS_ERR_OT5     7       /* other error */
#define PBS_ERR_OT6     8       /* other error */

/* add a point to the output polygon set
   NOTE:
   This macro relies on the existance of: outPID, outSID, outPOS, outX, outY,
   and outVerts. */
#define ADD_PT(pid, sid, pos, x, y)                                     \
        {                                                               \
          outPID[*outVerts] = pid;                                      \
          outSID[*outVerts] = sid;                                      \
          outPOS[*outVerts] = pos;                                      \
          outX[*outVerts] = x;                                          \
          outY[*outVerts] = y;                                          \
          (*outVerts)++;                                                \
        }

#define FREE(ptr)       { if (ptr) { free(ptr); (ptr) = NULL; }}

/* order for N, E, S, W must match `limits' arg. passed into these functions
   from S-Plus/R */
typedef enum {
  EDGE_W = 0,
  EDGE_E,
  EDGE_S,
  EDGE_N,
  NUM_EDGES,
  NO_EDGE
} edge;

typedef enum {
  C_NW = 0,
  C_NE,
  C_SE,
  C_SW,
  NUM_CORNERS,
  C_NONE
} corner;

int cornerXY[NUM_CORNERS][2] = {
  {EDGE_W, EDGE_N},
  {EDGE_E, EDGE_N},
  {EDGE_E, EDGE_S},
  {EDGE_W, EDGE_S}
};

jmp_buf returnLocation;


/*-----------------------------------------------------------------------------
  polyStartsEnds:
    This function determines the indices (0 .. n-1) where polygons
    start and end.
  
  Author:  Nicholas Boers (June 11, 2003)
  
  Returns: 
    the number of polygons (start/end pairs).
  ---------------------------------------------------------------------------*/
static PBSINT
polyStartsEnds(PBSINT *polyStarts, PBSINT *polyEnds,
               PBSINT *inPID, PBSINT *inSID, PBSINT *inVerts)
{
  PBSINT curPID, curSID;
  PBSINT count = 0;
  PBSINT i;

  /* ensure at least one polygon must exist */
  if (*inVerts == 0)
    return 0;

  /* set up for the first one */
  curPID = inPID[0];
  curSID = inSID[0];

  /* add the start of the first one */
  *polyStarts++ = 0;
  count++;
  
  /* walk through all the vertices (less the first one)... */
  for (i = 1; i < *inVerts; i++) {
    if (inPID[i] != curPID ||
        inSID[i] != curSID) {
      curPID = inPID[i];
      curSID = inSID[i];

      *polyEnds++ = i - 1;
      *polyStarts++ = i;
      count++;
    }
  }

  /* add the end of the last one */
  *polyEnds = i - 1;

  return count;
}

/*-----------------------------------------------------------------------------
  clip:
    This function clips polygons to a rectangular viewing window.
  
  Author:  Nicholas Boers (June 11, 2003)
  
  Implementation Notes:
    For each pair of points that are tested, the _first_ point of the
    pair is added to the final output.  If necessary, an intersection
    with the border is added as well.
  
  Notes:
    Recommended allocated space for out*: 2 x *inVerts
  ---------------------------------------------------------------------------*/
void
clip(PBSINT *inID, double *inXY, PBSINT *inVerts, PBSINT *polygons,
     double *limits, PBSINT *outID, double *outXY, PBSINT *outVerts,
     PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  PBSINT *inPOS = inID + (2 * (*inVerts));
  double *inX = inXY;
  double *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  PBSINT *outPOS = outID + (2 * (*outVerts));
  PBSINT *outOLD = outID + (3 * (*outVerts));
  double *outX = outXY;
  double *outY = outXY + (*outVerts);

  /* miscellaneous variables */
  PBSINT nVerts;
  PBSINT nPolys;
  PBSINT i;
  PBSINT j;
  PBSINT pos;
  PBSINT firstVertex;
  PBSINT tempOutVerts;
  PBSINT allocatedSpace = *outVerts;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts;
  PBSINT *polyEnds;

  short isHole;

  polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (!polyStarts || !polyEnds) {
    (*status) = PBS_ERR_MEM;
    goto CLIP_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    firstVertex = *outVerts;
    nVerts = (polyEnds[i] - polyStarts[i] + 1);
    if (nVerts > 1)
      isHole = (inPOS[polyStarts[i]] > inPOS[polyStarts[i] + 1]);
    else
      isHole = FALSE;

    tempOutVerts = allocatedSpace - *outVerts;
    /* walk through the four edges */
    clipPolygon(&inX[polyStarts[i]], &inY[polyStarts[i]], 
                &inPOS[polyStarts[i]], nVerts, 
                &outX[*outVerts], &outY[*outVerts], 
                &outOLD[*outVerts], &tempOutVerts, limits, (short)(*polygons));

    /* clipPolygon() returns whether or not it was successful using *outVerts:
       -1 on insufficient memory
       -2 on insufficient space in out* */
    if (*outVerts < 0) {
      (*status) = (*outVerts == -1) ? PBS_ERR_MEM : PBS_ERR_OUT;
      goto CLIP_FREE_MEM;
    }
    
    *outVerts += tempOutVerts;

    pos = (isHole) ? tempOutVerts : 1;
    /* update the PIDs/SIDs/POSs for the clipped polygon*/
    for (j = firstVertex; j < *outVerts; j++) {
      outPID[j] = inPID[polyStarts[i]];
      outSID[j] = inSID[polyStarts[i]];
      outPOS[j] = (isHole) ? pos-- : pos++;
    }
  }

  (*status) = PBS_SUCCESS;

 CLIP_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

#define ROLL_ADD_PT(pid, sid, pos, x, y)                                      \
          {                                                                   \
            outPID[*outVerts] = pid;                                          \
            outSID[*outVerts] = sid;                                          \
            outPOS[*outVerts] = pos;                                          \
            outX[*outVerts] = x;                                              \
            outY[*outVerts] = y;                                              \
            (*outVerts)++;                                                    \
          }
 
