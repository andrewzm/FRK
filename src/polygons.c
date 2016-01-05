/*=============================================================================
  Copyright (C) 2003-2013  Fisheries and Oceans Canada

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
  File: polygons.c
   
  Implementation of polygon manipulation routines.

  History:
    08 Jun 2004 [Nicholas Boers]
      - added thickenPoly()
    09 Jun 2004 [Nicholas Boers]
      - modifications to calcPolygonArea(); make code easier to understand
      - added calcPolyOrientation() to determine whether vertices are
        clockwise or counter-clockwise
    10 Jun 2004 [Nicholas Boers]
      - finished calcPolyOrientation()
      - renamed calcPolygonArea() to calcPolyArea()
      - with the addition of calcPolyOrientation(), calcPolyArea() initially
        appeared redundant; after comparing the speed of the computations done
        natively done in R with the C routine, the C routine was substantially
        faster
    11 Jun 2004 [Nicholas Boers]
      - added LL support to thickenPolys()
    16 Jun 2004 [Nicholas Boers]
      - merged convexhull.c into this file
      - cleaned up comments (most are now in polygons.h to reduce duplication)
    17 Jun 2004 [Nicholas Boers]
      - additional cleaning
    21 Jun 2004 [Nicholas Boers]
      - modified interface to calcPolyArea() -- added `area' as an argument,
        and changed return type to `short'; now it can return negative areas
        if polygon goes CW
    05 Jul 2004 [Nicholas Boers]
      - renamed calcConvexHull() to convexHull()
    13 Jul 2004 [Nicholas Boers]
      - fixed huge bug in calcPolyCentroid()
    16 Aug 2004 [Nicholas Boers]
      - added isPolyConvex(), from PBSmapping.c
      - added nPolyIntersect(), from PBSmapping.c
      - added isRetrace(), from PBSmapping.c
      - added linesIntersect(), from PBSmapping.c
    14 Nov 2004 [Nicholas Boers]
      - updated for PBSINT type
    22 Jan 2008 [Rowan Haigh]
      - removed convexHull(), sortPointList(), and rightTurn()
    22 Jun 2010 [Nicholas Boers]
      - require *inVerts > 0 in sutherlandHodgmanPolygonClip
  ---------------------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "floating.h"
#include "polygons.h"

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

/* used in convexHull() */
typedef struct _point {
  double x;
  double y;
} point;

/* used in convexHull() */
#define APPEND(xList, yList, n, x, y)  { xList[n] = x; yList[n] = y; n++; }

/* calculate the X-value for the intersection of the line
   (x1, y1) -> (x2, y2) with the horizontal line at boundary */
#define X_INTERCEPT(x1, x2, y1, y2, boundary)                            \
                (((boundary - y1) * (x1 - x2) / (y1 - y2)) + x1)

/* calculate the Y-value for the intersection of the line
   (x1, y1) -> (x2, y2) with the vertical line at boundary */
#define Y_INTERCEPT(x1, x2, y1, y2, boundary)                            \
                (((boundary - x1) * (y1 - y2) / (x1 - x2)) + y1)

/* calculate intersection with an edge */
#define EDGE_INTERSECTION(x1, y1, x2, y2, xOut, yOut, limits, e)         \
        {                                                                \
          /* calculate the intersection */                               \
          switch(e) {                                                    \
          case EDGE_N:                                                   \
          case EDGE_S:                                                   \
            xOut = X_INTERCEPT(x1, x2, y1, y2, limits[e]);               \
            yOut = limits[e];                                            \
            break;                                                       \
          case EDGE_E:                                                   \
          case EDGE_W:                                                   \
            xOut = limits[e];                                            \
            yOut = Y_INTERCEPT(x1, x2, y1, y2, limits[e]);               \
            break;                                                       \
          case NO_EDGE:                                                  \
          default:                                                       \
            break;                                                       \
          }                                                              \
        }

/*-----------------------------------------------------------------------------
  inside:
    Returns whether or not a point is 'within' a given edge.
    Used in sutherlandHodgmanPolygonClip().
  ---------------------------------------------------------------------------*/
static short 
inside(double x, double y, double *limits, edge e)
{
  switch(e) 
    {
    case EDGE_N:
      return (DBL_LTEQ(y, limits[EDGE_N]));
    case EDGE_E:
      return (DBL_LTEQ(x, limits[EDGE_E]));
    case EDGE_S:
      return (DBL_GTEQ(y, limits[EDGE_S]));
    case EDGE_W:
      return (DBL_GTEQ(x, limits[EDGE_W]));
    case NO_EDGE:
    default:
      return FALSE;
    }
}

/*-----------------------------------------------------------------------------
  sutherlandHodgmanPolygonClip: 
    Clip a single polygon against one edge.
    Used in clipPolygon().
   
  Author:  Nicholas Boers (June 18, 2003)

  Algorithm source:
    Computer Graphics: Principles and Practice
  
  Notes:
    - recommended allocated space for out*: 2 x *inVerts

  Returns:
    *outVerts >= 0 on success
    *outVerts = -1 on insufficient memory allocated in out*
  ---------------------------------------------------------------------------*/
static void
sutherlandHodgmanPolygonClip(double *inX, double *inY, PBSINT *inPOS,
                             PBSINT *inVerts,
                             double *outX, double *outY, PBSINT *outOLD,
                             PBSINT *outVerts,
                             double *limits, edge e, short polygons)
{
  double sX, sY, pX, pY;   /* for walking through polygon */
  double iX = 0, iY = 0;   /* for intersections */
  PBSINT allocatedMemory = *outVerts;
  PBSINT i;
 
  *outVerts = 0;

  if (*inVerts > 0) {
    sX = inX[(*inVerts) - 1];
    sY = inY[(*inVerts) - 1];
  }
  
  /* Idea:
     - add p and any intersections on each iteration of the loop

     The points:
       - p is the 'current' point
       - s in the 'previous' point */

  for (i = 0; i < *inVerts; i++) {
    pX = inX[i];
    pY = inY[i];

    /* if clipping lines, don't test the edge connecting the end -> start */
    if (i == 0 && !polygons) {
      /* because we always add the segment's end point (rather than start
         point, simply add the end point here if necessary */
      if (inside(pX, pY, limits, e)) {
        if (*outVerts == allocatedMemory) {
          *outVerts = -1;
          return;
        }
        outX[*outVerts] = pX;
        outY[*outVerts] = pY;
        outOLD[*outVerts] = inPOS[i];
        (*outVerts)++;
      }
      sX = pX;
      sY = pY;
      continue;
    }

    /* 4 cases:
       #1: s inside     p inside
       #2: s inside     p outside
       #3: s outside    p outside
       #4: s outside    p inside */

    if (inside(pX, pY, limits, e)) {            /* case #1 & #4 */
      if (inside(sX, sY, limits, e)) {          /* case #1 */
        if (*outVerts == allocatedMemory) {
          *outVerts = -1;
          return;
        }
        outX[*outVerts] = pX;
        outY[*outVerts] = pY;
        outOLD[*outVerts] = inPOS[i];
        (*outVerts)++;
      }
      else {                                    /* case #4 */
        if (*outVerts == allocatedMemory) {
          *outVerts = -1;
          return;
        }
        EDGE_INTERSECTION(sX, sY, pX, pY, iX, iY, limits, e);
        outX[*outVerts] = iX;
        outY[*outVerts] = iY;
        outOLD[*outVerts] = NO_OLD_POS;
        (*outVerts)++;

        if (*outVerts == allocatedMemory) {
          *outVerts = -1;
          return;
        }
        outX[*outVerts] = pX;
        outY[*outVerts] = pY;
        outOLD[*outVerts] = inPOS[i];
        (*outVerts)++;
      }
    }
    else {                                      /* case #2 & #3 */
      if (inside(sX, sY, limits, e)) {          /* case #2 */
        if (*outVerts == allocatedMemory) {
          *outVerts = -1;
          return;
        }
        EDGE_INTERSECTION(sX, sY, pX, pY, iX, iY, limits, e);
        outX[*outVerts] = iX;
        outY[*outVerts] = iY;
        outOLD[*outVerts] = NO_OLD_POS;
        (*outVerts)++;
      }                                         /* no action for case #3 */
    }
    
    sX = pX;
    sY = pY;
  }
}

/*-----------------------------------------------------------------------------
  clipPolygon:
  ---------------------------------------------------------------------------*/
void
clipPolygon(double *inX, double *inY, PBSINT *inPOS, PBSINT inVerts,
            double *outX, double *outY, PBSINT *outOLD, PBSINT *outVerts,
            double *limits, short polygons)
{
  /* create a second copy of the output variables, for temporary work */
  double *tempX = (double *) malloc(sizeof(double) * (*outVerts));
  double *tempY = (double *) malloc(sizeof(double) * (*outVerts));
  PBSINT *tempOLD = (PBSINT *) malloc(sizeof(PBSINT) * (*outVerts));
  PBSINT tempVerts;

  PBSINT allocatedMemory = (*outVerts);
  edge e;

  const short ERROR_NO_MEM = -1;
  const short ERROR_NO_OUT = -2;

  if (!tempX || !tempY || !tempOLD) {
    if (!tempX) free(tempX);
    if (!tempY) free(tempY);
    if (!tempOLD) free(tempOLD);

    (*outVerts) = ERROR_NO_MEM;
    return;
  }

  if (inVerts > (*outVerts)) {
    free(tempX);
    free(tempY);
    free(tempOLD);

    (*outVerts) = ERROR_NO_OUT;
    return;
  }
  memcpy(tempX, inX, sizeof(double) * inVerts);
  memcpy(tempY, inY, sizeof(double) * inVerts);
  memcpy(tempOLD, inPOS, sizeof(PBSINT) * inVerts);
  tempVerts = inVerts;

  /* clip against each edge */
  for (e = 0; e < NUM_EDGES; e++) {
    (*outVerts) = allocatedMemory;

    sutherlandHodgmanPolygonClip(tempX, tempY, tempOLD, &tempVerts,
                                 outX, outY, outOLD, outVerts,
                                 limits, e, polygons);

    /* sutherlandHodgmanPolygonClip() returns whether or not it was
       successful using tempOutVerts:
       -1 on insufficient memory allocated to tempOut* */
    if ((*outVerts) == -1) {
      free(tempX);
      free(tempY);
      free(tempOLD);

      *outVerts = ERROR_NO_OUT; /* out* were too small */
      return;
    }

    /* only copy it over if we'll need it again */
    if ((e + 1) < NUM_EDGES) {
      memcpy(tempX, outX, sizeof(double) * (*outVerts));
      memcpy(tempY, outY, sizeof(double) * (*outVerts));
      memcpy(tempOLD, outOLD, sizeof(PBSINT) * (*outVerts));
      tempVerts = (*outVerts);
    }
  }
  
  free(tempX);
  free(tempY);
  free(tempOLD);
}
