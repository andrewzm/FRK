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
  pnpoly:
    This function computes whether or not an array of points lie
    outside, on the boundary of, or inside a polygon.  The results are
    returned in `results' (an array).
  
    The polygon can optionally be closed: (x_0, y_0) == (x_n, y_n).
  
  Author:  Nicholas Boers (June 13, 2003)
  ---------------------------------------------------------------------------*/
void
pnpoly(PBSINT *polyPts, double *polyX, double *polyY,
       PBSINT *pts, double *x, double *y,
       PBSINT *results)
{
  PBSINT i;
  double limits[NUM_EDGES];

  if (*polyPts > 0) {
    /* calculate the limits of the polygon */
    limits[EDGE_W] = limits[EDGE_E] = polyX[0];
    limits[EDGE_S] = limits[EDGE_N] = polyY[0];
    for (i = 1; i < *polyPts; i++) {
      if (DBL_LT(polyX[i], limits[EDGE_W])) limits[EDGE_W] = polyX[i];
      if (DBL_GT(polyX[i], limits[EDGE_E])) limits[EDGE_E] = polyX[i];
      if (DBL_LT(polyY[i], limits[EDGE_S])) limits[EDGE_S] = polyY[i];
      if (DBL_GT(polyY[i], limits[EDGE_N])) limits[EDGE_N] = polyY[i];
    }

    /* walk through the points */
    for (i = 0; i < *pts; i++) {
      /* do a quick test to see if the event is outside the limits */
      if (DBL_LT(x[i], limits[EDGE_W]) || DBL_GT(x[i], limits[EDGE_E]) ||
	  DBL_LT(y[i], limits[EDGE_S]) || DBL_GT(y[i], limits[EDGE_N]))
	results[i] = OUTSIDE;
      else
	results[i] = pointInPolygon(polyX, polyY, *polyPts, x[i], y[i]);
    }
  } else {
    /* when the polygon consists of 0 points, set all points
       as being outside (can this really happen?) */
    for (i = 0; i < *pts; i++) {
      results[i] = OUTSIDE;
    }
  }
}

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
 
/*-----------------------------------------------------------------------------
  rollupPolys:
    Performs several operations on a PolySet.  See the arguments below for
    further details.

  Non-standard parameters:
    rollupMode: method for rolling up the PolySet
       1 = roll-up to the PID level (only PIDs in the result)
       2 = roll-up to the outer contour level (only outer contours in the
           result)
       3 = do not roll-up

    exteriorCCW: modify vertices orientation (CW/CCW)?
      -1 = don't modify
       0 = exterior should be CW
      +1 = exterior should be CCW
    
    closedPolys: whether the last and first vertices should be the same
      -1 = don't modify
       0 = ensure polygons do not close
      +1 = close the polygons

    addRetrace: determines whether it adds retrace lines to the first vertex
    of the parent after outputting a child
       0 = don't add
       1 = add

  Author:  Nicholas Boers (June 17, 2004)
  
  Status values:
    PBS_SUCCESS: everything OK
    PBS_ERR_MEM: insufficient memory
    PBS_ERR_OUT: output array full
    PBS_ERR_OT1: encountered child where not allowed

  Notes:
    - maximum allocated space required for out*:
      *inVerts + 1 (for each it needs to close) + 1 (for each retrace line)
    - recalculates the "POS" column
  ---------------------------------------------------------------------------*/
void rollupPolys(PBSINT *inID, double *inPOS, double *inXY, PBSINT *inVerts,
                 PBSINT *outID, double *outXY, PBSINT *outVerts,
                 PBSINT *rollupMode, PBSINT *exteriorCCW, PBSINT *closedPolys, 
                 PBSINT *addRetrace, PBSINT *status)
{
  /* 3 rollup modes */
  const PBSINT MODE_SID = 1;
  const PBSINT MODE_HOLE = 2;
  const PBSINT MODE_NONE = 3;

  /* make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  /* double *inPOS (argument) */
  double *inX = inXY;
  double *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  PBSINT *outPOS = outID + (2 * (*outVerts));
  double *outX = outXY;
  double *outY = outXY + (*outVerts);

  PBSINT allocatedMemory = (*outVerts);

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;

  /* for rollup == 1 and 2; ensure first PID != parentPID */
  PBSINT parentPID = (*inPID) - 1, parentSID = 0;
  short allowedChild = FALSE;

  /* for addRetrace */
  double parentX = 0, parentY = 0;

  /* miscellaneous variables */
  PBSINT i;
  PBSINT lastPOS = 0;

  /* set outVerts now in case we return early */
  (*outVerts) = 0;
  if (!polyStarts || !polyEnds) {
    (*status) = PBS_ERR_MEM;
    goto ROLLUPPOLYS_FREE_MEM;
  }

  /* VALIDATE THE INPUT ARGUMENTS */

  /* calculate polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through polygons */
  for (i = 0; i < nPolys; i++) 
    {
      PBSINT nVerts  = polyEnds[i] - polyStarts[i] + 1;
      short isHole = ((*inVerts) < 2) ? FALSE :
        (DBL_GT(inPOS[polyStarts[i]], inPOS[polyStarts[i] + 1]));
      short revVerts = FALSE;
      short isParent = FALSE;
      PBSINT j;

      double *X, *Y, *firstX, *firstY;
      short addFirst = FALSE;

      /* deal with `exteriorCCW': set `revVerts' to TRUE if we need
         to reverse the vertices */
      if ((*exteriorCCW) != -1)
        {
          short orientation = calcPolyOrientation(&inX[polyStarts[i]], 
                                                  &inY[polyStarts[i]], nVerts);
          /* determine whether it is oriented like a hole */
          /* exterior goes CW and this one goes CCW
             OR exterior goes CCW and this one goes CW */
          short likeHole = ((!(*exteriorCCW) && (orientation == -1)) ||
                            ((*exteriorCCW) && (orientation == 1)));
          /* if oriented like a hole but not a hole, then reverse vertices */
          if (likeHole != isHole)
            revVerts = TRUE;
        }
      /* revVerts == TRUE if we need to reverse vertices */

      isParent = (((*rollupMode) == MODE_NONE)
                  || (((*rollupMode) == MODE_SID) 
                      && inPID[polyStarts[i]] != parentPID) 
                  || (((*rollupMode) == MODE_HOLE) && !isHole));
      
      /* deal with `closedPolys': shrink `nVerts' to not include points,
         or set `addFirst' to TRUE to add the first one */
      if ((*closedPolys) != -1)
        {
          /* if ensure they don't close... */
          if ((*closedPolys) == 0)
            {
              PBSINT start = polyStarts[i];
              
              while ((inX[start] == inX[start + nVerts - 1])
                     && (inY[start] == inY[start + nVerts - 1])
                     && (nVerts > 1))
                nVerts--;
            }
          /* if ensure they close... */
          else if ((*closedPolys) == 1)
            {
              if ((inX[polyStarts[i]] != inX[polyEnds[i]])
                  || (inY[polyStarts[i]] != inY[polyEnds[i]]))
                addFirst = TRUE;
            }
        }
      /* nVerts only includes necessary points,
         addFirst == TRUE if we need to add the first point */
      
      firstX = X = (revVerts) ? &inX[polyEnds[i]] : &inX[polyStarts[i]];
      firstY = Y = (revVerts) ? &inY[polyEnds[i]] : &inY[polyStarts[i]];
      
      /* set parent variables */
      if (isParent)
        {
          parentPID = inPID[polyStarts[i]];
          parentSID = inSID[polyStarts[i]];
          parentX = *X;
          parentY = *Y;

          lastPOS = (isHole && ((*rollupMode) == MODE_NONE)) ? nVerts : 1;
          if (addFirst && isHole && ((*rollupMode) == MODE_NONE)) lastPOS++;

          allowedChild = TRUE;
        }
      else if (!allowedChild)
        {
          (*status) = PBS_ERR_OT1;
          goto ROLLUPPOLYS_FREE_MEM;
        }
      
      /* copy over polygon */
      for (j = 0; j < nVerts; j++)
        {
          if ((*outVerts) == allocatedMemory)
            {
              (*status) = PBS_ERR_OUT;
              goto ROLLUPPOLYS_FREE_MEM;
            }
          ROLL_ADD_PT(parentPID, parentSID, 
                      (isHole && ((*rollupMode) == MODE_NONE))
                      ? lastPOS-- : lastPOS++,
                      *X, *Y);
          if (revVerts) { X--; Y--; } else { X++; Y++; }
        }
      
      /* deal with special cases */
      if (addFirst)
        {
          if ((*outVerts) == allocatedMemory)
            {
              (*status) = PBS_ERR_OUT;
              goto ROLLUPPOLYS_FREE_MEM;
            }
          ROLL_ADD_PT(parentPID, parentSID, 
                      (isHole && ((*rollupMode) == MODE_NONE)) 
                      ? lastPOS-- : lastPOS++,
                      *firstX, *firstY);
        }
      if (!isParent && (*addRetrace))
        {
          if ((*outVerts) == allocatedMemory)
            {
              (*status) = PBS_ERR_OUT;
              goto ROLLUPPOLYS_FREE_MEM;
            }
          ROLL_ADD_PT(parentPID, parentSID, 
                      (isHole && ((*rollupMode) == MODE_NONE))
                      ? lastPOS-- : lastPOS++,
                      parentX, parentY);
        }
    }
  
  (*status) = PBS_SUCCESS;
  
 ROLLUPPOLYS_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  closestCorner:
    Given two points (stored in X and Y), returns the corner that is
    closest to one of those points.
  
  Author:  Nicholas Boers
  ---------------------------------------------------------------------------*/
static corner
findClosestCorner(double x[2], double y[2], double *limits)
{
  int i, j;
  double temp;

  double dist = DIST(x[0], y[0], limits[EDGE_W], limits[EDGE_N]);
  corner closest = C_NW;

  for (i = 0; i < NUM_CORNERS; i++) {
    for (j = 0; j < 2; j++) {
      temp = DIST(x[j], y[j], 
                  limits[cornerXY[i][0]], limits[cornerXY[i][1]]);
      if (temp < dist) {
        dist = temp;
        closest = i;
      }
    }
  }

  return(closest);
}

/*-----------------------------------------------------------------------------
  closePolys: 
    "Fix" the closure of open polygons.
 
  Author:  Nicholas Boers
  
  Notes:
    - recommended allocated space for out*: 

    - cornerToAdd was developed as follows:
                      EDGE_N (3)
           C_NW (0)                 C_NE (1)
                    ---------------
                    |             |
        EDGE_W (0)  |             |  EDGE_E (1)
                    |             |     
                    ---------------
           C_SW (3)                 C_SE (2)
                       EDGE_S (2)
   
     - from the above picture, create a table:
       START EDGE | END EDGE   | CLOSEST CORNER || 1ST TO ADD | 2ND TO ADD
       -----------|------------|----------------||------------|-----------
       EDGE_W     | EDGE_W     | C_NW           || C_NONE     | C_NONE     
       EDGE_W     | EDGE_W     | C_NE           || C_NONE     | C_NONE     
       EDGE_W     | EDGE_W     | C_SE           || C_NONE     | C_NONE     
       EDGE_W     | EDGE_W     | C_SW           || C_NONE     | C_NONE     
       EDGE_W     | EDGE_E     | C_NW           || C_NE       | C_NW       
       EDGE_W     | EDGE_E     | C_NE           || C_NE       | C_NW       
       EDGE_W     | EDGE_E     | C_SE           || C_SE       | C_SW       
       EDGE_W     | EDGE_E     | C_SW           || C_SE       | C_SW       
       ...        | ...        | ...            || ...        | ...        
  
     - the table should be completed for all combinations of START EDGE,
       END EDGE, and CLOSEST CORNER
  
     - the cornerToAdd array must follow the order in the table above, where
       CLOSEST corner is cycled through first for each (START EDGE, END EDGE),
       and then EDGE EDGE is cycled through for each START EDGE
   ---------------------------------------------------------------------------*/
void
closePolys(PBSINT *inID, double *inXY, PBSINT *inVerts, double *limits, 
           PBSINT *outID, double *outXY, PBSINT *outVerts,
           PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (1 * (*inVerts));
  PBSINT *inPOS = inID + (2 * (*inVerts));
  double *inX = inXY;
  double *inY = inXY + (1 * (*inVerts));
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (1 * (*outVerts));
  PBSINT *outPOS = outID + (2 * (*outVerts));
  double *outX = outXY;
  double *outY = outXY + (1 * (*outVerts));

  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;
  PBSINT nVerts;

  short isHole;

  edge edges[2];
  corner closestCorner;
  corner cornerToAdd1, cornerToAdd2;

  double x[2], y[2];

  PBSINT i, j;

  const int START = 0;
  const int END   = 1;
 
  const int X = 0;
  const int Y = 1;

  /* see the note above for how this array was developed... */
  corner cornerToAdd[2][64] = {
    {C_NONE, C_NONE, C_NONE, C_NONE, C_NE, C_NE, C_SE, C_SE, C_SW,
     C_SW, C_SW, C_SW, C_NW, C_NW, C_NW, C_NW, C_NW, C_NW, C_SW, C_SW,
     C_NONE, C_NONE, C_NONE, C_NONE, C_SE, C_SE, C_SE, C_SE, C_NE,
     C_NE, C_NE, C_NE, C_SW, C_SW, C_SW, C_SW, C_SE, C_SE, C_SE, C_SE,
     C_NONE, C_NONE, C_NONE, C_NONE, C_NW, C_NE, C_NE, C_NW, C_NW,
     C_NW, C_NW, C_NW, C_NE, C_NE, C_NE, C_NE, C_SW, C_SE, C_SE, C_SW,
     C_NONE, C_NONE, C_NONE, C_NONE},
    {C_NONE, C_NONE, C_NONE, C_NONE, C_NW, C_NW, C_SW, C_SW, C_NONE,
     C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NE,
     C_NE, C_SE, C_SE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE,
     C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE,
     C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NONE,
     C_NONE, C_NONE, C_SW, C_SE, C_SE, C_SW, C_NONE, C_NONE, C_NONE,
     C_NONE, C_NONE, C_NONE, C_NONE, C_NONE, C_NW, C_NE, C_NE, C_NW,
     C_NONE, C_NONE, C_NONE, C_NONE}
  };

  PBSINT allocatedMemory = (*outVerts);

  *outVerts = 0;

  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto CLOSEPOLYS_FREE_MEM;
  }

  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  for (i = 0; i < nPolys; i++) {
    nVerts = (polyEnds[i] - polyStarts[i]) + 1;
    if (nVerts > 1)
      isHole = (inPOS[polyStarts[i]] > inPOS[polyStarts[i] + 1]);
    else
      isHole = FALSE;
    
    /* calculate which edge each one falls on */
    edges[START] = NO_EDGE;
    if (DBL_EQ(inX[polyStarts[i]], limits[EDGE_W])) edges[START] = EDGE_W;
    if (DBL_EQ(inX[polyStarts[i]], limits[EDGE_E])) edges[START] = EDGE_E;
    if (DBL_EQ(inY[polyStarts[i]], limits[EDGE_S])) edges[START] = EDGE_S;
    if (DBL_EQ(inY[polyStarts[i]], limits[EDGE_N])) edges[START] = EDGE_N;
    
    edges[END] = NO_EDGE;
    if (DBL_EQ(inX[polyEnds[i]], limits[EDGE_W])) edges[END] = EDGE_W;
    if (DBL_EQ(inX[polyEnds[i]], limits[EDGE_E])) edges[END] = EDGE_E;
    if (DBL_EQ(inY[polyEnds[i]], limits[EDGE_S])) edges[END] = EDGE_S;
    if (DBL_EQ(inY[polyEnds[i]], limits[EDGE_N])) edges[END] = EDGE_N;

    /* calculate which corner is closest */
    x[START] = inX[polyStarts[i]];
    y[START] = inY[polyStarts[i]];
    x[END] = inX[polyEnds[i]];
    y[END] = inY[polyEnds[i]];

    closestCorner = C_NONE;
    /* findClosestCorner() doesn't dynamically allocate memory */
    closestCorner = findClosestCorner(x, y, limits);

    /* create indices into ugly 2 * 64 array */
    cornerToAdd1 = cornerToAdd[START][16*edges[START] + 
                                      4*edges[END] + closestCorner];
    cornerToAdd2 = cornerToAdd[END][16*edges[START] + 
                                    4*edges[END] + closestCorner];
    
    /* if it's a hole */
    if (isHole && edges[START] != NO_EDGE && edges[END] != NO_EDGE) {
      if (cornerToAdd1 != C_NONE && cornerToAdd2 != C_NONE) {
        if (((*outVerts) + 1) >= allocatedMemory) {
          (*status) = PBS_ERR_OUT;
          goto CLOSEPOLYS_FREE_MEM;
        }
        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyStarts[i]] + 2,
               limits[cornerXY[cornerToAdd1][X]], 
               limits[cornerXY[cornerToAdd1][Y]]);

        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyStarts[i]] + 1,
               limits[cornerXY[cornerToAdd2][X]], 
               limits[cornerXY[cornerToAdd2][Y]]);
      } else if (cornerToAdd1 != C_NONE) {
        if ((*outVerts) >= allocatedMemory) {
          (*status) = PBS_ERR_OUT;
          goto CLOSEPOLYS_FREE_MEM;
        }
        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyStarts[i]] + 1,
               limits[cornerXY[cornerToAdd1][X]], 
               limits[cornerXY[cornerToAdd1][Y]]);
      }
    }

    /* copy the current points over */
    for (j = 0; j < nVerts; j++) {
      if ((*outVerts) >= allocatedMemory) {
        (*status) = PBS_ERR_OUT;
        goto CLOSEPOLYS_FREE_MEM;
      }
      ADD_PT(inPID[polyStarts[i] + j], inSID[polyStarts[i] + j],
                inPOS[polyStarts[i] + j], inX[polyStarts[i] + j],
                inY[polyStarts[i] + j]);
    }

    /* if it isn't a hole */
    if (!isHole && edges[START] != NO_EDGE && edges[END] != NO_EDGE) {
      if (cornerToAdd1 != C_NONE && cornerToAdd2 != C_NONE) {
        if (((*outVerts) + 1) >= allocatedMemory) {
          (*status) = PBS_ERR_OUT;
          goto CLOSEPOLYS_FREE_MEM;
        }
        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyEnds[i]] + 1,
               limits[cornerXY[cornerToAdd1][X]], 
               limits[cornerXY[cornerToAdd1][Y]]);

        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyEnds[i]] + 2,
               limits[cornerXY[cornerToAdd2][X]],
               limits[cornerXY[cornerToAdd2][Y]]);
      } else if (cornerToAdd1 != C_NONE) {
        if ((*outVerts) >= allocatedMemory) {
          (*status) = PBS_ERR_OUT;
          goto CLOSEPOLYS_FREE_MEM;
        }
        ADD_PT(inPID[polyStarts[i]], inSID[polyStarts[i]],
               inPOS[polyEnds[i]] + 1,
               limits[cornerXY[cornerToAdd1][X]],
               limits[cornerXY[cornerToAdd1][Y]]);
      }
    }
  }

  (*status) = PBS_SUCCESS;

 CLOSEPOLYS_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  findCells:
    Identify integers within a vector of break points.

  Author:  Nicholas Boers (Mar. 29, 2006)

  Notes:
    Some ideas from "findInterval()" in R source code.
  ---------------------------------------------------------------------------*/
void
findCells(double *inPt, PBSINT *inPts,
          double *inBrk, PBSINT *inBrks,
          PBSINT *outCell, PBSINT *outBdry,
          PBSINT *status)
{
  double *pts = inPt;
  int nPts = (*inPts);

  double *brks = inBrk;
  int n = (*inBrks);

  int i;
  int ilo, imid, ihi;

#define outside { outBdry[i] = 0; outCell[i] = -1; continue; }
#define inside { outBdry[i] = DBL_EQ(pt, brks[ilo]); outCell[i] = ilo;  continue; }

  --brks;       /* 1-based indexing */

  for (i = 0; i < nPts; i++) {
    double pt = pts[i];

    /* out of range? */
    if (DBL_LT(pt, brks[1]))
      outside;
    if (DBL_GT(pt, brks[n]))
      outside;

    ilo = 1;
    ihi = n;
    
    /* narrow range with essentially a binary search */
    for (; (ihi - ilo) > 1;) {
      imid = (ilo + ihi) / 2;
      /* "pt" less than "brks[imid]" */
      if (DBL_LT(pt, brks[imid])) {
        ihi = imid - 1;
      } 
      /* "pt" greater than or equal to "brks[imid]" */
      else {
        ilo = imid;
      }
    }

    /* if two cells, determine correct one */
    if ((ilo != ihi) && !(DBL_LT(pt, brks[ihi])))
      ilo++;

    inside
  }

  (*status) = PBS_SUCCESS;
}

/*-----------------------------------------------------------------------------
  findPolys:
    Locate events within a polyset.
  
  Author:  Nicholas Boers
  
  Notes:
    - recommended allocated space for out*: 
  ---------------------------------------------------------------------------*/
void
findPolys(PBSINT *inEventsID, double *inEventsXY, PBSINT *inEvents,
          PBSINT *inPolysID, double *inPolysXY, PBSINT *inPolys,
          PBSINT *outID, PBSINT *outIDs,
          PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  double *inEventsX = inEventsXY;
  double *inEventsY = inEventsXY + (1 * (*inEvents));

  PBSINT *inPolysPID = inPolysID;
  PBSINT *inPolysSID = inPolysID + (1 * (*inPolys));
  PBSINT *inPolysPOS = inPolysID + (2 * (*inPolys));
  double *inPolysX = inPolysXY;
  double *inPolysY = inPolysXY + (1 * (*inPolys));

  PBSINT *outEID = outID;
  PBSINT *outPID = outID + (1 * (*outIDs));
  PBSINT *outSID = outID + (2 * (*outIDs));
  PBSINT *outBdry = outID + (3 * (*outIDs));

  PBSINT i, j;

  PBSINT nVerts;
  PBSINT nPolys;
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inPolys));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inPolys));

  PBSINT *results = (PBSINT *) malloc(sizeof(PBSINT) * (*inEvents));
  PBSINT *resultsTemp = (PBSINT *) malloc(sizeof(PBSINT) * (*inEvents));
  short *boundary = (short *) malloc(sizeof(short) * (*inEvents));

  PBSINT parentPID = -1, parentSID = -1;
  short isHole, isNext;

  PBSINT allocatedMemory = (*outIDs);
  *outIDs = 0;

  if (results == NULL || polyStarts == NULL || polyEnds == NULL ||
      resultsTemp == NULL || boundary == NULL) {
    (*status) = PBS_ERR_MEM;
    goto FINDPOLYS_FREE_MEM;
  }

  nPolys = polyStartsEnds(polyStarts, polyEnds, inPolysPID, 
                          inPolysSID, inPolys);

  /* walk through all the polygons */
  for (i = 0; i < nPolys; i++) {
    /* run points-in-polygons */
    nVerts = polyEnds[i] - polyStarts[i] + 1;
    /* does not dynamically allocate memory */
    pnpoly(&nVerts, &inPolysX[polyStarts[i]], &inPolysY[polyStarts[i]],
           inEvents, inEventsX, inEventsY,
           results);

    /* if !hole, then set "parent" variables */
    if (!(isHole = (inPolysPOS[polyStarts[i]] > 
                    inPolysPOS[polyStarts[i] + 1]))) {
      parentPID = inPolysPID[polyStarts[i]];
      parentSID = inPolysSID[polyStarts[i]];
    }

    /* process the results from points-in-polygons */
    if (!isHole) {
      for (j = 0; j < *inEvents; j++) {
        resultsTemp[j] = (results[j] == INSIDE || results[j] == BOUNDARY);
        boundary[j] = (results[j] == BOUNDARY);
      }
    }
    else {
      for (j = 0; j < *inEvents; j++) {
        resultsTemp[j] = (resultsTemp[j] && !(results[j] == INSIDE));
        boundary[j] = boundary[j] || (results[j] == BOUNDARY);
      }
    }

    /* determine if the next one is a hole */
    if ((isNext = ((i + 1) < nPolys)))
      /* if the polySet is valid, the PID check is useless... */
      isHole = ((inPolysPID[polyStarts[i]] == 
                 inPolysPID[polyStarts[i + 1]]) &&
                (inPolysPOS[polyStarts[i + 1]] > 
                 inPolysPOS[polyStarts[i + 1] + 1]));

    /* output when necessary */
    if ((isNext && !isHole) || (!isNext)) {
      for (j = 0; j < *inEvents; j++) {
        if (resultsTemp[j]) {
          if ((*outIDs) >= allocatedMemory) {
            (*status) = PBS_ERR_OUT;
            goto FINDPOLYS_FREE_MEM;
          }
          outEID[*outIDs] = inEventsID[j];
          outPID[*outIDs] = parentPID;
          outSID[*outIDs] = parentSID;
          outBdry[*outIDs] = boundary[j];
          (*outIDs)++;
        }
      }
    }
  }

  (*status) = PBS_SUCCESS;

 FINDPOLYS_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);
  FREE(results);
  FREE(resultsTemp);
  FREE(boundary);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  convUL:
    Convert Lon/Lat <--> UTME/UTMN.
  
  Author:  Nicholas Boers
  
  Notes:
    Maximum space required for out*:
  ---------------------------------------------------------------------------*/
void
convUL(double *inXY, PBSINT *inVerts, PBSINT *toUTM, PBSINT *zone,
       PBSINT *southern, double *outXY, PBSINT *outVerts,
       PBSINT *status)
{
  double *inX = inXY;
  double *inY = inXY + (*inVerts);

  double *outX = outXY;
  double *outY = outXY + (*outVerts);

  double num1, num2;
  struct pair results;
  PBSINT lineno = 0;
  PBSINT i;

  if ((*inVerts) > (*outVerts)) {
    (*status) = PBS_ERR_OUT;
    return;
  }

  for (i = 0; i < *inVerts; i++) {
    lineno++;
    num1 = inX[i];
    num2 = inY[i];

    if(*toUTM) {
      num1 *= DEG_TO_RAD;
      num2 *= DEG_TO_RAD;
      /* no dynamic memory allocation */
      lonlat_to_utm(num1, num2, (int)(*zone), &results);
    }
    else {
      /* no dynamic memory allocation */
      utm_to_lonlat(num1, num2, (*southern ? 'S' : 'N'), (int)(*zone), &results);
      results.x *= RAD_TO_DEG;
      results.y *= RAD_TO_DEG;
    }
    
    outX[i] = results.x;
    outY[i] = results.y;
  }

  *outVerts = (lineno == *inVerts) ? *inVerts : 0;
  (*status) = PBS_SUCCESS;

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  calcArea:
    This function calculates the areas of a set of polygons.
    It handles holes, but has no concept of projection.
  
  Author:  Nicholas Boers (July 11, 2003)

  Capabilities:
    [Y] Holes
    [>] Projection
        [N] LL
        [Y] UTM
        [Y] 1:1
  Robustness:
    [?] Floating-point
  Legend:  
    [Y] yes   [N] no   [-] N/A   [?] good question   [>] details

  Notes:
    - maximum space required for out*: length(unique(paste(polys$PID, 
                                                           polys$SID)))
  ---------------------------------------------------------------------------*/
void
calcArea(PBSINT *inID, double *inXY, PBSINT *inVerts,
         PBSINT *outID, double *outArea, PBSINT *outVerts,
         PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  /* PBSINT *inPOS = inID + (2 * (*inVerts)); */
  double *inX = inXY;
  double *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);

  /* miscellaneous variables */
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;
  PBSINT nVerts;

  double area;

  PBSINT allocatedMemory = (*outVerts);

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto CALCAREA_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    nVerts = polyEnds[i] - polyStarts[i] + 1;

    /* calculate the polygon's area */
    if (calcPolyArea(&inX[polyStarts[i]], &inY[polyStarts[i]], &area, nVerts)
        < 0)
      {
        (*status) = PBS_ERR_MEM;
        goto CALCAREA_FREE_MEM;
      }

    /* output area */
    if ((*outVerts) >= allocatedMemory) {
      (*status) = PBS_ERR_OUT;
      goto CALCAREA_FREE_MEM;
    }
    outPID[*outVerts] = inPID[polyStarts[i]];
    outSID[*outVerts] = inSID[polyStarts[i]];
    outArea[*outVerts] = area;
    (*outVerts)++;
  }

  (*status) = PBS_SUCCESS;

 CALCAREA_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  calcCentroid:
    This function calculates the centroids of a set of polygons.
    Since it uses signed area calculations, it can handle holes if they are
    properly described (in a CW/CCW GIS sense) in the PolySet, and they are
    "integrated" (holes have the same PID/SID as their parent). It has no
    concept of projection.
  
  Author:  Nicholas Boers (June 10, 2004)

  Capabilities:
    [>] Holes
        Yes in certain circumstances (see note above).
    [-] Projection
  Robustness:
    [?] Floating-point
  Legend:  
    [Y] yes   [N] no   [-] N/A   [?] good question   [>] details

  Notes:
    - maximum space required for out*: length(unique(paste(polys$PID, 
                                                           polys$SID)))
  ---------------------------------------------------------------------------*/
void
calcCentroid(PBSINT *inID, double *inXY, PBSINT *inVerts,
             PBSINT *outID, double *outXY, PBSINT *outVerts,
             PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  /* PBSINT *inPOS = inID + (2 * (*inVerts)); */
  double *inX = inXY;
  double *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  double *outX = outXY;
  double *outY = outXY + (*outVerts);

  /* miscellaneous variables */
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;
  PBSINT nVerts;

  PBSINT allocatedMemory = (*outVerts);

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (!polyStarts || !polyEnds) {
    (*status) = PBS_ERR_MEM;
    goto CALCCENTROID_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    short rValue;

    nVerts = polyEnds[i] - polyStarts[i] + 1;

    if ((*outVerts) >= allocatedMemory) {
      (*status) = PBS_ERR_OUT;
      goto CALCCENTROID_FREE_MEM;
    }
    rValue = calcPolyCentroid(&inX[polyStarts[i]], &inY[polyStarts[i]], nVerts,
                              &outX[*outVerts], &outY[*outVerts]);
    if (rValue == -2) {
      (*status) = PBS_ERR_MEM;
      goto CALCCENTROID_FREE_MEM;
    }
    if (rValue == -1) {
      continue;
    }
    outPID[*outVerts] = inPID[polyStarts[i]];
    outSID[*outVerts] = inSID[polyStarts[i]];
    (*outVerts)++;
  }

  (*status) = PBS_SUCCESS;

 CALCCENTROID_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  calcOrientation:
    This function calculates the orientation of polygons (i.e., clockwise or
    counter-clockwise).  Holes are irrelevant and it needs no concept of
    projection.
  
  Author:  Nicholas Boers (June 9, 2004)

  Capabilities:
    [-] Holes
    [-] Projection
  Robustness:
    [?] Floating-point
  Legend:  
    [Y] yes   [N] no   [-] N/A   [?] good question   [>] details

  Notes:
    - space required for out*: length(unique(paste(polys$PID, polys$SID)))

  Returns:
    -1 when counter-clockwise
     0 when N/A
    +1 when clockwise
  ---------------------------------------------------------------------------*/
void
calcOrientation(PBSINT *inID, double *inXY, PBSINT *inVerts,
                PBSINT *outID, double *outOrientation, PBSINT *outVerts,
                PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  /* PBSINT *inPOS = inID + (2 * (*inVerts)); */
  double *inX = inXY;
  double *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);

  /* miscellaneous variables */
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;
  PBSINT nVerts;

  PBSINT allocatedMemory = (*outVerts);

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto CALCORIENT_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    nVerts = polyEnds[i] - polyStarts[i] + 1;

    /* calculate the polygon's area */
    if ((*outVerts) == allocatedMemory) {
      (*status) = PBS_ERR_OUT;
      goto CALCORIENT_FREE_MEM;
    }
    outOrientation[*outVerts] = 
      (PBSINT)calcPolyOrientation(&inX[polyStarts[i]], &inY[polyStarts[i]], nVerts);
    outPID[*outVerts] = inPID[polyStarts[i]];
    outSID[*outVerts] = inSID[polyStarts[i]];
    (*outVerts)++;
  }

  (*status) = PBS_SUCCESS;

 CALCORIENT_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  isConvex:
    Determines whether a PolySet contains convex polygons.

  Author:  Nicholas Boers (June 30, 2004)
  
  Status values:
    PBS_SUCCESS:    everything OK
    PBS_ERR_MEM:    insufficient memory
    PBS_ERR_OUT:    output array full

  Notes:
    - maximum allocated space required for out*: number of polygons
  ---------------------------------------------------------------------------*/
void
isConvex(PBSINT *inID, double *inXY, PBSINT *inVerts,
         PBSINT *outID, PBSINT *outResult, PBSINT *outVerts,
         PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  double *inX = inXY;
  double *inY = inXY + (*inVerts);

  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  
  /* miscellaneous variables */
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;

  PBSINT allocatedMemory = (*outVerts);

  /* set `outVerts' now rather than later in case we return early (below) */
  *outVerts = 0;
  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto ISCONVEX_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);
  /* check for sufficient memory for output*/
  if (allocatedMemory < nPolys) {
    (*status) = PBS_ERR_OUT;
    goto ISCONVEX_FREE_MEM;
  }

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    PBSINT nVerts = polyEnds[i] - polyStarts[i] + 1;
    short result;

    result = isPolyConvex(&inX[polyStarts[i]], &inY[polyStarts[i]], nVerts);

    /* add result to output */
    outPID[*outVerts] = inPID[polyStarts[i]];
    outSID[*outVerts] = inSID[polyStarts[i]];
    outResult[*outVerts] = result;
    (*outVerts)++;
  }

  (*status) = PBS_SUCCESS;

 ISCONVEX_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  isIntersecting:
    Determines whether the polygons in a PolySet self-intersect.

  Non-standard parameters:
    numericResult: if TRUE, returns a numeric result (a count); otherwise,
      returns a Boolean for whether or not the polygon self-intersects

  Author:  Nicholas Boers (June 28, 2004)
  
  Status values:
    PBS_SUCCESS:    everything OK
    PBS_ERR_MEM:    insufficient memory
    PBS_ERR_OUT:    output array full

  Notes:
    - maximum allocated space required for out*: number of polygons
    - counts certain types of intersections (i.e., those involving vertices
      and those where an edge retraces over an edge) more than once
  ---------------------------------------------------------------------------*/
void
isIntersecting(PBSINT *inID, double *inXY, PBSINT *inVerts,
               PBSINT *numericResult,
               PBSINT *outID, PBSINT *outResult, PBSINT *outVerts,
               PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  /* PBSINT *inPOS = inID + (2 * (*inVerts)); */
  double *inX = inXY;
  double *inY = inXY + (*inVerts);

  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  
  /* miscellaneous variables */
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT nPolys;

  PBSINT allocatedMemory = (*outVerts);

  /* set `outVerts' now rather than later in case we return early (below) */
  *outVerts = 0;
  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto ISINTERSECTING_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);
  /* check for sufficient memory for output*/
  if (allocatedMemory < nPolys) {
    (*status) = PBS_ERR_OUT;
    goto ISINTERSECTING_FREE_MEM;
  }

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    PBSINT nVerts = polyEnds[i] - polyStarts[i] + 1;
    PBSINT count = 0;

    count = nPolyIntersects(&inX[polyStarts[i]], &inY[polyStarts[i]],
                            nVerts, (short)(*numericResult));

    /* add result to output */
    outPID[*outVerts] = inPID[polyStarts[i]];
    outSID[*outVerts] = inSID[polyStarts[i]];
    outResult[*outVerts] = (*numericResult) ? count : (count > 0);
    (*outVerts)++;
  }

  (*status) = PBS_SUCCESS;

 ISINTERSECTING_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  thickenPolys:
    This function thickens polygons.
  
  Author:  Nicholas Boers (June 8, 2004)
  
  Notes:
    - if units == 0, "LL": tolerance in kilometers and inXY in decimal-degrees
    - if units == 1, other: tolerance and inXY in same units
  ---------------------------------------------------------------------------*/
void
thickenPolys(PBSINT *inID, double *inXY, PBSINT *inVerts,
             double *tolerance, PBSINT *filter, PBSINT *units, PBSINT *keepOrig,
             PBSINT *close,
             PBSINT *outID, double *outXY, PBSINT *outVerts,
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
  double *outX = outXY;
  double *outY = outXY + (*outVerts);

  /* miscellaneous variables */
  PBSINT nPolys;
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));

  PBSINT allocatedSpace = (*outVerts);
  PBSINT tempOutVerts;
  PBSINT pos;

  short isHole;

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto THICKENPOLYS_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    PBSINT nVerts = polyEnds[i] - polyStarts[i] + 1;

    if (nVerts > 1)
      isHole = (inPOS[polyStarts[i]] > inPOS[polyStarts[i] + 1]);
    else
      isHole = FALSE;
    
    tempOutVerts = thickenPoly(&inX[polyStarts[i]], &inY[polyStarts[i]], 
                               nVerts, &outX[*outVerts], &outY[*outVerts],
                               (allocatedSpace - *outVerts),
                               (*tolerance), (short)(*units),
                               (short)(*keepOrig), (short)(*close));
    if (tempOutVerts < 0) {
      (*status) = PBS_ERR_OUT;
      goto THICKENPOLYS_FREE_MEM;
    }
    else if (tempOutVerts >= (*filter)) {
      PBSINT j;

      pos = (isHole) ? tempOutVerts : 1;
      for (j = 0; j < tempOutVerts; j++) {
        outPID[*outVerts] = inPID[polyStarts[i]];
        outSID[*outVerts] = inSID[polyStarts[i]];
        outPOS[*outVerts] = (isHole) ? pos-- : pos++;
        (*outVerts)++;
      }
    }
  }

  (*status) = PBS_SUCCESS;

 THICKENPOLYS_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}

/*-----------------------------------------------------------------------------
  thinPolys:
    This function thins polygons.
  
  Author:  Nicholas Boers (May 4, 2004)
  
  Notes:
    - X and Y are `PBSINT,' rather than the usual double
    - recommended allocated space for out*:
      *inVerts
    - does not renumber POS
    - if units == 0, "LL", and inXY should be in micro-degrees
    - if units == 1, "UTM", and inXY should be in meters
  ---------------------------------------------------------------------------*/
void
thinPolys(PBSINT *inID, PBSINT *inXY, PBSINT *inVerts, double *tolerance,
          PBSINT *filter, PBSINT *units,
          PBSINT *outID, PBSINT *outXY, PBSINT *outVerts,
          PBSINT *status)
{
  /* declarations to make accessing values user-friendly */
  PBSINT *inPID = inID;
  PBSINT *inSID = inID + (*inVerts);
  PBSINT *inPOS = inID + (2 * (*inVerts));
  PBSINT *inX = inXY;
  PBSINT *inY = inXY + (*inVerts);
  
  PBSINT *outPID = outID;
  PBSINT *outSID = outID + (*outVerts);
  PBSINT *outPOS = outID + (2 * (*outVerts));
  PBSINT *outX = outXY;
  PBSINT *outY = outXY + (*outVerts);

  /* miscellaneous variables */
  PBSINT nPolys;
  PBSINT i;

  /* keep track of where polygons start and end */
  PBSINT *polyStarts = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *polyEnds = (PBSINT *) malloc(sizeof(PBSINT) * (*inVerts));
  PBSINT *index = NULL;

  PBSINT allocatedMemory = (*outVerts);

  /* set it now rather than later in case we return early (below) */
  *outVerts = 0;

  if (polyStarts == NULL || polyEnds == NULL) {
    (*status) = PBS_ERR_MEM;
    goto THINPOLYS_FREE_MEM;
  }

  /* calculate the polygon starts and ends */
  nPolys = polyStartsEnds(polyStarts, polyEnds, inPID, inSID, inVerts);

  /* walk through the various polygons */
  for (i = 0; i < nPolys; i++) {
    int n, j;
    PBSINT nVerts;

    nVerts = polyEnds[i] - polyStarts[i] + 1;
    
    index = (PBSINT *)malloc(sizeof(PBSINT) * nVerts);
    if (index == NULL) {
      (*status) = PBS_ERR_MEM;
      goto THINPOLYS_FREE_MEM;
    }

    n = Douglas_Peucker_i(&inX[polyStarts[i]], &inY[polyStarts[i]],
                          nVerts, *tolerance, index, (short)(*units));
    if (n < 0) {
      (*status) = PBS_ERR_MEM;
      goto THINPOLYS_FREE_MEM;
    }

    if (n >= *filter) {
      for (j = 0; j < n; j++) {
        if ((*outVerts) >= allocatedMemory) {
          (*status) = PBS_ERR_OUT;
          goto THINPOLYS_FREE_MEM;
        }
        outPID[*outVerts] = inPID[polyStarts[i] + index[j]];
        outSID[*outVerts] = inSID[polyStarts[i] + index[j]];
        outPOS[*outVerts] = inPOS[polyStarts[i] + index[j]];
        outX[*outVerts] = inX[polyStarts[i] + index[j]];
        outY[*outVerts] = inY[polyStarts[i] + index[j]];
        
        (*outVerts)++;
      }
    }

    FREE(index);
  }

  (*status) = PBS_SUCCESS;

 THINPOLYS_FREE_MEM:
  FREE(polyStarts);
  FREE(polyEnds);
  FREE(index);

#ifdef _MCD_CHECK
  showMemStats();
#endif /* _MCD_CHECK */
}
