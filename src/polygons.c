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

/*-----------------------------------------------------------------------------
  pointInPolygon:
  ---------------------------------------------------------------------------*/
short
pointInPolygon(double *inX, double *inY, PBSINT inVerts, double x, 
               double y)
{
  short result, inside;
  PBSINT i, j;
  double yCalc;

  /* each time we start a new point, reset */
  result = UNKNOWN;
  inside = FALSE;
    
  /* walk through the polygon */
  for (i = 0, j = inVerts-1; i < inVerts; j = i++) {
    /* in X-range */
    if ((DBL_LTEQ(inX[i], x) && DBL_LT(x, inX[j])) ||
        (DBL_LTEQ(inX[j], x) && DBL_LT(x, inX[i]))) {

      yCalc = (inY[i] - inY[j]) / (inX[i] - inX[j]) * (x - inX[j]) + inY[j];

      if (DBL_LT(y, yCalc))
        inside = !inside;
      else if (DBL_EQ(y, yCalc))
        result = BOUNDARY;

    }      
    else if ((DBL_LTEQ(inX[i], x) &&  DBL_LTEQ(x, inX[j])) ||
             (DBL_LTEQ(inX[j], x) &&  DBL_LTEQ(x, inX[i]))) {
      
      if ((DBL_EQ(x, inX[i]) && DBL_EQ(y, inY[i])) ||
          (DBL_EQ(x, inX[j]) && DBL_EQ(y, inY[j])))
        result = BOUNDARY;
      else if (DBL_EQ(inX[i], inX[j]) &&
               ((DBL_LTEQ(inY[i], y) && DBL_LTEQ(y, inY[j])) ||
                (DBL_LTEQ(inY[j], y) && DBL_LTEQ(y, inY[i]))))
        result = BOUNDARY;

    }
  }

  if (result == UNKNOWN)
    result = (inside) ? INSIDE : OUTSIDE;

  return(result);
}

/*-----------------------------------------------------------------------------
  calcPolyArea:
  ---------------------------------------------------------------------------*/
short
calcPolyArea(double *inX, double *inY, double *area, PBSINT inVerts)
{
  double *x = (double *)malloc(sizeof(double) * inVerts);
  double *y = (double *)malloc(sizeof(double) * inVerts);
  double *x0, *x1, *y0, *y1;
  PBSINT i;

  (*area) = 0;

  if ((x == NULL) || (y == NULL)) {
    if (x != NULL) free(x);
    if (y != NULL) free(y);

    return -1;
  }
  if (inVerts < 3) {
    (*area) = 0;
    free(x);
    free(y);

    return 0;
  }

  /* move the origin temporarily to a location near the vertices to
     improve precision (especially applicable with large numbers) */
  for (i = 0; i < inVerts; i++) {
    x[i] = inX[i] - inX[0];
    y[i] = inY[i] - inY[0];
  }

  x0 = &x[0];  x1 = &x[1];
  y0 = &y[0];  y1 = &y[1];

  /* SIMILAR statements to calcPolyCentroid() */
  for (i = 0; i < (inVerts - 1); i++)
    (*area) += (x0[i] * y1[i] - x1[i] * y0[i]);
  /* if first and last points are different, process closing edge */
  if (!DBL_EQ(x[0], x[inVerts-1]) || !DBL_EQ(y[0], y[inVerts-1]))
    (*area) += (x[inVerts - 1] * y[0] - x[0] * y[inVerts - 1]);
  (*area) = (*area) / 2.0;

  free(x);
  free(y);

  return 0;
}

/*-----------------------------------------------------------------------------
  calcPolyCentroid:
  ---------------------------------------------------------------------------*/
short
calcPolyCentroid(double *inX, double *inY, PBSINT inVerts,
                 double *outX, double *outY)
{
  double *term, *x, *y;
  double *x0, *x1, *y0, *y1;
  double area = 0.0;
  double sumX = 0.0, sumY = 0.0;
  PBSINT i;

  /* special cases */
  if (inVerts < 1)
    return -1;
  if (inVerts == 1) {
    (*outX) = inX[0];
    (*outY) = inY[0];
    return 0;
  }
  if (inVerts == 2) {
    (*outX) = (inX[0] + inX[1]) / 2.0;
    (*outY) = (inY[0] + inY[1]) / 2.0;
    return 0;
  }

  /* get memory */
  term = (double *)malloc(sizeof(double) * inVerts);
  /* +1 for case where first and last points are different */
  x = (double *)malloc(sizeof(double) * inVerts);
  y = (double *)malloc(sizeof(double) * inVerts);
  if (!x || !y || !term) {
    if (x) free(x);
    if (y) free(y);
    if (term) free(term);
    return -2;
  }

  /* translate the origin temporarily to a location near the vertices to
     improve precision (especially applicable with large numbers) */
  for (i = 0; i < inVerts; i++) {
    x[i] = inX[i] - inX[0]; 
    y[i] = inY[i] - inY[0]; 
  }
  /* set up pointers for easier access */
  x0 = &x[0];  x1 = &x[1];
  y0 = &y[0];  y1 = &y[1];

  /* SIMILAR statements to calcPolyArea() */
  for (i = 0; i < (inVerts - 1); i++)
    {
      area += term[i] = (x0[i] * y1[i] - x1[i] * y0[i]);
      sumX += (x0[i] + x1[i]) * term[i];
      sumY += (y0[i] + y1[i]) * term[i];
    }
  /* if first and last points are different, process closing edge */
  if (!DBL_EQ(x[0], x[inVerts-1]) || !DBL_EQ(y[0], y[inVerts-1])) {
    area += term[inVerts-1] = (x[inVerts-1] * y[0] - x[0] * y[inVerts-1]);
    sumX += (x[inVerts-1] + x[0]) * term[inVerts-1];
    sumY += (y[inVerts-1] + y[0]) * term[inVerts-1];
  }
  area /= 2.0;

  (*outX) = sumX / (6.0 * area);
  (*outY) = sumY / (6.0 * area);

  /* translate centroid back */
  (*outX) += inX[0];
  (*outY) += inY[0];

  free(x);
  free(y);
  free(term);

  return 0;
}

/*-----------------------------------------------------------------------------
  calcPolyOrientation:
  ---------------------------------------------------------------------------*/
short
calcPolyOrientation(double *inX, double *inY, PBSINT inVerts)
{
  double minY;
  PBSINT minIdx;
  PBSINT next, prev;
  PBSINT update;
  PBSINT i;

  minY = inY[0];
  minIdx = 0;

  /* find a minimum (may occur more than once) */
  for (i = 1; i < inVerts; i++)
    {
      if (inY[i] < minY)
        {
          minY = inY[i];
          minIdx = i;
        }
    }

  /* search left then right in array for right-most bottom point */
  next = minIdx;
  for (update = -1; update < 2; update += 2)
    {
      /* only loop up to a maximum to prevent infinite loop if line */
      for (i = 0; i < inVerts; i++)
        {
          next += update;
          if (next == inVerts) next = 0;
          if (next < 0) next = inVerts - 1;

          if (!DBL_EQ(minY, inY[next]) ||
              DBL_LT(inX[next], inX[minIdx]))
            break;

          minIdx = next;
        }
      if (i == inVerts)
        return 0;
    }

  /* now have a minimum, right-most point; search for a point "previous"
     to this point */
  prev = minIdx - 1;
  if (prev < 0) prev = inVerts - 1;
  while (DBL_EQ(inX[minIdx], inX[prev]) &&
         DBL_EQ(inY[minIdx], inY[prev]))
    {
      prev--;
      if (prev < 0) prev = inVerts - 1;
    }

  /* have a minimum and the point before it; search for a point "next" to
     this point that is not the previous point or the current point */
  next = minIdx + 1;
  if (next == inVerts) next = 0;
  /* only loop up to a maximum to prevent infinite loop if line */
  for (i = 0; i < inVerts; i++)
    {
      if (((!DBL_EQ(inX[minIdx], inX[next]))
           || (!DBL_EQ(inY[minIdx], inY[next])))
          && ((!DBL_EQ(inX[prev], inX[next]))
              || (!DBL_EQ(inY[prev], inY[next]))))
        break;

      next++;
      if (next == inVerts) next = 0;
    }
  if (i == inVerts)
    return 0;

  /* calculate whether the cross-product z-vector points up (positive) or
     down (negative). If it points up, then clockwise. */
  if (((inX[prev] - inX[minIdx]) * (inY[next] - inY[minIdx])
       - (inX[next] - inX[minIdx]) * (inY[prev] - inY[minIdx])) > 0)
    return 1;
  else
    return -1;
}

/*-----------------------------------------------------------------------------
  isPolyConvex:
  ---------------------------------------------------------------------------*/
short
isPolyConvex(double *inX, double *inY, PBSINT inVerts)
{
  /* if no self-intersections, check for convex */
  if ((inVerts > 2) && (nPolyIntersects(inX, inY, inVerts, FALSE) == 0)) {
    PBSINT lastIdx = inVerts - 1;
    PBSINT f;                   /* first non-straight */
    double temp;
    PBSINT turn;
    PBSINT i;

    /* if same first and last points, decrement the last point */
    if (DBL_EQ(inX[0], inX[lastIdx]) && DBL_EQ(inY[0], inY[lastIdx]))
      lastIdx--;

    /* Walk the polygon, calculating the Z component of the
       cross-product for each corner.  If the first corner is a right
       turn, and all subsequent corners are either straight or a right
       turn, then it is convex. The same pattern applies to left
       turns. */
    /* Test all the corners from 1 to nVerts-1; for the corner at 0, and
       the (possible) corner at nVerts, handle as special cases. */

    /* Calculate the first turn before starting the loop to determine
       what all subsequent turns must be. */
    temp = ((inX[lastIdx] - inX[0]) * (inY[1] - inY[0]) - 
            (inY[lastIdx] - inY[0]) * (inX[1] - inX[0]));
    /* extra tests in case the first three points are colinear */
    f = 1;
    while (DBL_EQ(temp, 0) && (f < lastIdx)) {
      temp = (((inX[f-1] - inX[f]) * (inY[f+1] - inY[f])) - 
              ((inY[f-1] - inY[f]) * (inX[f+1] - inX[f])));
      f++;
    }
    turn = (DBL_GT(temp, 0));

    /* check the turns between the first and the last */
    for (i = f; i < lastIdx; i++) {
      /* Calculate whether the cross-product z-vector points up (positive) or
         down (negative). If it points up, then clockwise. */
      temp = ((inX[i-1] - inX[i]) * (inY[i+1] - inY[i]) - 
              (inY[i-1] - inY[i]) * (inX[i+1] - inX[i]));
      if (!DBL_EQ(temp, 0) && (DBL_GT(temp, 0) != turn))
        return FALSE;
    }  

    /* check the last turn */
    temp = ((inX[lastIdx-1] - inX[lastIdx]) * (inY[0] - inY[lastIdx]) - 
            (inY[lastIdx-1] - inY[lastIdx]) * (inX[0] - inX[lastIdx]));
    if (!DBL_EQ(temp, 0) && (DBL_GT(temp, 0) != turn))
      return FALSE;
  }
  else {
    return FALSE;
  }
  
  return TRUE;
}

/*-----------------------------------------------------------------------------
  linesIntersect:
    Determines whether two line segments intersect.

  Source:
    Adapted from `insectc.c'.
    (http://http://www.acm.org/pubs/tog/GraphicsGems/category.html

  Author:  Franklin Antonio

  Modifications:
    29 Jun 2004 [Nicholas Boers]
      - converted to double precision (from PBSINT)

  Returns:
    0 if no intersection
    1 if intersection
    2 if parallel
  ---------------------------------------------------------------------------*/
static int 
linesIntersect(double x1, double y1, double x2, double y2,
               double x3, double y3, double x4, double y4)
{
  const int NO_INTERSECT = 0;
  const int INTERSECT = 1;
  const int PARALLEL = 2;

  double Ax, Bx, Cx;
  double Ay, By, Cy;

  double x1lo, x1hi, y1lo, y1hi;

  double d, e, f;

  Ax = x2 - x1;
  Bx = x3 - x4;

  /* X bound box test*/
  if (DBL_LT(Ax, 0))
    {
      x1lo = x2; 
      x1hi = x1;
    }
  else
    {
      x1hi = x2;
      x1lo = x1;
    }
  if (DBL_GT(Bx, 0)) 
    {
      if (DBL_LT(x1hi, x4) || DBL_LT(x3, x1lo))
        return NO_INTERSECT;
    }
  else
    {
      if (DBL_LT(x1hi, x3) || DBL_LT(x4, x1lo))
        return NO_INTERSECT;
    }
  
  Ay = y2 - y1;
  By = y3 - y4;

  /* Y bound box test*/
  if (DBL_LT(Ay, 0)) 
    {
      y1lo = y2;
      y1hi = y1;
    }
  else
    {
      y1hi = y2;
      y1lo = y1;
    }
  if (DBL_GT(By, 0))
    {
      if (DBL_LT(y1hi, y4) || DBL_LT(y3, y1lo))
        return NO_INTERSECT;
    }
  else
    {
      if (DBL_LT(y1hi, y3) || DBL_LT(y4, y1lo))
        return NO_INTERSECT;
    }

  Cx = x1 - x3;
  Cy = y1 - y3;

  d = By * Cx - Bx * Cy;                                /* alpha numerator*/
  f = Ay * Bx - Ax * By;                                /* both denominator*/

  /* alpha tests*/
  if (DBL_GT(f, 0))
    {
      if (DBL_LT(d, 0) || DBL_GT(d, f))
        return NO_INTERSECT;
    } 
  else 
    {
      if (DBL_GT(d, 0) || DBL_LT(d, f))
        return NO_INTERSECT;
    }

  e = Ax*Cy - Ay*Cx;                                    /* beta numerator*/

  /* beta tests */
  if (DBL_GT(f, 0))
    {
      if (DBL_LT(e, 0) || DBL_GT(e, f))
        return NO_INTERSECT;
    } 
  else
    {
      if (DBL_GT(e, 0) || DBL_LT(e, f))
        return NO_INTERSECT;
    }

  if (DBL_EQ(f, 0))
    return PARALLEL;

#ifdef UNDEF
  /*compute intersection coordinates*/
  num = d*Ax;                                   /* numerator */
  offset = SAME_SIGNS(num,f) ? f/2 : -f/2;      /* round direction*/
  *x = x1 + (num+offset) / f;                   /* intersection x */

  num = d*Ay;
  offset = SAME_SIGNS(num,f) ? f/2 : -f/2;
  *y = y1 + (num+offset) / f;                   /* intersection y */
#endif /* UNDEF */

  return INTERSECT;
}

/*-----------------------------------------------------------------------------
  isRetrace:
    A line segment goes from P0 to P1.  A second line segment goes from
    P1 to Pt.  Does the second line segment (P1 to Pt) retrace over the first?

  Non-standard parameters:
    P0x, P0y: coordinates of P0

    P1x, P1y: coordinates of P1

    Ptx, Pty: coorindates of Pt

  Author:  Nicholas Boers (June 29, 2004)

  Return values:
    TRUE if Pt retraces over P0 --> P1.  Otherwise, FALSE.

  Notes:
    - based on the parametric equation of a line
      P(t) = P0 + (P1 - P0)*t
      Since we know P(t), rearrange the equation like:
      P(t) - P0
      --------- = t
       P1 - P0
      Calculate the above equation separately for X and Y coordinates.
    - if P0x/y and P1x/y are the same point, returns TRUE if Ptx/y is also
      the same point
    - if Ptx/y is the same as P1x/y, or it otherwise 'retraces' over
      the segment P0->P1, then it returns TRUE
  ---------------------------------------------------------------------------*/
static short 
isRetrace(double P0x, double P0y, double P1x, double P1y,
          double Ptx, double Pty)
{
  double numX =   Ptx - P0x;
  double numY =   Pty - P0y;
  double denomX = P1x - P0x;
  double denomY = P1y - P0y;

  if (DBL_EQ(denomX, 0) || DBL_EQ(denomY, 0)) {
    /* if the line doesn't leave this point, then it's a retrace */
    if (DBL_EQ(P1x, Ptx) && DBL_EQ(P1y, Pty)) {
      return TRUE;
    }
    /* if no change in X, check for retrace in Y */
    else if (DBL_EQ(denomX, 0) && DBL_EQ(P1x, Ptx)) {
      if ((DBL_GT(P1y, P0y) && !DBL_GT(Pty, P1y)) ||
          (DBL_LT(P1y, P0y) && !DBL_LT(Pty, P1y)))
        return TRUE;
    }
    /* if no change in Y, check for retrace in X */
    else if (DBL_EQ(denomY, 0) && DBL_EQ(P1y, Pty)) {
      if ((DBL_GT(P1x, P0x) && !DBL_GT(Ptx, P1x)) ||
          (DBL_LT(P1x, P0x) && !DBL_LT(Ptx, P1x)))
        return TRUE;
    }
  }
  else {
    double tX = numX / denomX;
    double tY = numY / denomY;

    if (DBL_EQ(tX, tY) && !DBL_GT(tX, 1))
      return TRUE;
  }
  
  return FALSE;
}

/*-----------------------------------------------------------------------------
  nPolyIntersects:
  ---------------------------------------------------------------------------*/
PBSINT
nPolyIntersects(double *inX, double *inY, PBSINT inVerts, short numericResult)
{
  PBSINT lastIdx = inVerts - 1;
  PBSINT count = 0;
  PBSINT j, k;

  /* since the polygon may be 'optionally' closed... */
  if (DBL_EQ(inX[0], inX[lastIdx]) && DBL_EQ(inX[0], inX[lastIdx]))
    lastIdx--;
  
  /* lastIdx == 0 means one point; not self-intersecting */
  if (lastIdx == 0)
    return 0;
  /* lastIdx == 1 means two points; self-intersecting */
  if (lastIdx == 1)
    return 1;
  
  /*
    Algorithm:
    
    Test all edges against the edge that closes the polygon.  The
    edge immediately before and immediately after this closing edge
    require special handling (as they share vertices).  If either
    of these two immediately neighbouring edges retraces over the
    closing edge, count an intersection.
    
    After checking the closing edge, 'walk' the polygon's edges to
    compare the current edge with all the edges that follow the
    current edge.  When checking the next edge, ensure that it
    doesn't colinearly retrace over the current edge (treat the
    next edge as a special case). 
  */
  
  /* handle the closing edge */
  /* P0 = lastIdx, P1 = 0, P(t) = 1 */
  if (isRetrace(inX[lastIdx],   inY[lastIdx],
                inX[0],         inY[0],
                inX[1],         inY[1])) {
    count++;
  }
  if (isRetrace(inX[0],         inY[0],
                inX[lastIdx],   inY[lastIdx],
                inX[lastIdx-1], inY[lastIdx-1])) {
    count++;
  }
  for (j = 1; j < lastIdx-1; j++) {
    /* compare the closing edge to each point (except it's nearest 
       neighbour on each side) */
    if (linesIntersect(inX[lastIdx],    inY[lastIdx],
                       inX[0],          inY[0],
                       inX[j],          inY[j],
                       inX[j+1],        inY[j+1])
        != 0) {
      count++;
    }
  }

  if (numericResult || (count == 0)) {
    /* now that we processed the closing edge, process the remaining
       edges */
    for (j = 0; j < lastIdx; j++) {
      for (k = j + 1; k < lastIdx; k++) {
        /* next edge is a special case */
        if (k == (j + 1)) {
          if (isRetrace(inX[j], inY[j], inX[k], inY[k], inX[k+1], inY[k+1])) {
            count++;
            if (!numericResult)
              j = k = lastIdx;
          }
        }
        else {
          if (linesIntersect(inX[j], inY[j], inX[j+1], inY[j+1],
                             inX[k], inY[k], inX[k+1], inY[k+1]) != 0) {
            count++;
            if (!numericResult)
              j = k = lastIdx;
          }
        } /* else */
      } /* for (k) */
    } /* for (j) */
  } /* if */

  return count;
}

/*-----------------------------------------------------------------------------
  thickenPoly:
  - uses Great Circle distances for LL
  - uses Pythagoras' Theorem for UTM and 1
  ---------------------------------------------------------------------------*/
PBSINT
thickenPoly(double *inX, double *inY, PBSINT inN,
            double *outX, double *outY, PBSINT outN,
            double tol, short units, short keepOrig, short close)
{
  const short LL = 0;
  double accDist = 0.0;
  PBSINT allocatedSpace = outN;
  PBSINT i;
  PBSINT ptLast;

  outN = 0;

  /* we will loop back to vertex 0 to close if necessary */
  if ((close) && ((inX[0] != inX[inN-1]) || (inY[0] != inY[inN-1])))
    ptLast = inN;
  else
    ptLast = inN - 1;

  /* SPECIAL CASE: first point, and !keepOrig */
  if (!keepOrig && (inN > 0))
    {
      /* add the first point, which is otherwise skipped */
      if (outN == allocatedSpace)
        return -1;
      outX[outN] = inX[0];
      outY[outN] = inY[0];
      outN++;
    }

  /* process all points */
  for (i = 0; i < ptLast; i++)
    {
      /* assign indices for points */
      PBSINT ptCurr = i;
      PBSINT ptNext = ((i + 1) == inN) ? 0 : (i + 1);
      /* calculate distance from current point (i) to next point (i+1) */
      double dX = inX[ptNext] - inX[ptCurr];
      double dY = inY[ptNext] - inY[ptCurr];
      double dist;
      /* when LL, calculate distances in kilometers rather than degrees, using
         Great Circles */
     if (units == LL)
        {
          /* added consideration for LL here (rather than in the tolerance) to
             improve the accuracy (compared results with UTM, and doing the
             computations here is more accurate */
          /* Formula source: http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
             Uses the Haversine Formula. */
          double dlon = (inX[ptNext] - inX[ptCurr]) * M_PI/180;
          double dlat = (inY[ptNext] - inY[ptCurr]) * M_PI/180;
          double a = (sin(dlat/2)) * (sin(dlat/2))
            + cos(inY[ptCurr]*M_PI/180) * cos(inY[ptNext]*M_PI/180)
            * (sin(dlon/2)) * (sin(dlon/2));
          a = sqrt(a);
          /* don't exceed the domain of asin() due to rounding */
          if (a > 1) a = 1;
          dist = MEAN_RADIUS * (2.0 * asin(a));
        }
      else
        {
          dist = sqrt((dX * dX) + (dY * dY));
        }

      if (keepOrig)
        {
          /* add the current point */
          if (outN == allocatedSpace)
            return -1;
          outX[outN] = inX[ptCurr];
          outY[outN] = inY[ptCurr];
          outN++;

          /* add intermediate points if necessary */
          if (dist > tol)
            {
              /* space them evenly across the segment */
              PBSINT segments = ceil(dist / tol);
              PBSINT j;
              
              for (j = 1; j < segments; j++)
                {
                  if (outN == allocatedSpace)
                    return -1;
                  outX[outN] = inX[ptCurr] + (dX * ((double)j / segments));
                  outY[outN] = inY[ptCurr] + (dY * ((double)j / segments));
                  outN++;
                }
            }
        }
      else /* !keepOrig */
        {
          accDist += dist;

          /* add intermediate points; for !keepOrig case */
          while (accDist >= tol)
            {
              double frac;
              
              accDist -= tol;
              frac = ((dist - accDist) / dist);
              
              if (outN == allocatedSpace)
                return -1;
              outX[outN] = inX[ptCurr] + (dX * frac);
              outY[outN] = inY[ptCurr] + (dY * frac);
              outN++;
            }
        }
    } /* process all points */

  return outN;
}

