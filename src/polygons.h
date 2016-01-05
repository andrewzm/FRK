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
  File: polygons.h
   
  Interfaces for polygon manipulation routines.

  History:
    08 Jun 2004 [Nicholas Boers]
      - added thickenPoly()
    09 Jun 2004 [Nicholas Boers]
      - modifications to calcPolygonArea(); make code easier to understand
      - added calcPolyOrientation() to determine whether vertices are clockwise
        or counter-clockwise
    16 Jun 2004 [Nicholas Boers]
      - moved several comments from polygons.c to polygons.h to avoid 
        duplication
    17 Jun 2004 [Nicholas Boers]
      - additional cleaning
    21 Jun 2004 [Nicholas Boers]
      - modified interface to calcPolyArea() -- added `area' as an argument,
        and changed return type to `short'
    07 Jul 2004 [Nicholas Boers]
      - renamed calcConvexHull() to convexHull()
    14 Nov 2004 [Nicholas Boers]
      - updated for PBSINT type
  ---------------------------------------------------------------------------*/

#ifndef _POLYGONS_H
#define _POLYGONS_H

/* possible return values from pointInPolygon */
#define UNKNOWN          -2
#define OUTSIDE          -1
#define BOUNDARY          0
#define INSIDE            1

/* value of oldPOS when clipPolygon generates a point */
#define NO_OLD_POS       -1

/* the distance between two points (x1, y1) -> (x2, y2) */
#define DIST(x1, y1, x2, y2)                                             \
                (sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)))

/*-----------------------------------------------------------------------------
  calcPolyArea:
    This function has no concept of holes.  The routine that calls it must
    support holes.  Additionally, it has no concept of projection.  If the
    vertices are oriented in a CCW fashion, area is positive.  Otherwise, it
    is negative.

  Source for algorithm:
    http://www.geog.ubc.ca/courses/klink/gis.notes/ncgia/u33.html

  Returns:
     0 on success.
    -1 on insufficient memory.
  ---------------------------------------------------------------------------*/
short calcPolyArea(double *inX, double *inY, double *area, PBSINT inVerts);

/*-----------------------------------------------------------------------------
  calcPolygonCentroid:
    This function has no concept of holes.  The routine that calls it must
    support holes.  Additionally, it has no concept of projection.

  Algorithm source:
    http://astronomy.swin.edu.au/~pbourke/geometry/polyarea/

  Returns:
    -2 on insufficient memory
    -1 on no polygon!
     0 on success
  ---------------------------------------------------------------------------*/
short calcPolyCentroid(double *inX, double *inY, PBSINT inVerts,
                       double *outX, double *outY);

/*-----------------------------------------------------------------------------
  calcPolyOrientation:
    Calculates whether the vertices of a polygon go in a clockwise or
    counter-clockwise order.

  Returns:
    -1 when counter-clockwise
     0 when N/A
    +1 when clockwise
 ---------------------------------------------------------------------------*/
short calcPolyOrientation(double *inX, double *inY, PBSINT inVerts);

/*-----------------------------------------------------------------------------
  clipPolygon:
  
  Author:  Nicholas Boers (June 18, 2003)

  Algorithm source:
    - Computer Graphics: Principles and Practice 
  
  Description:
    This function clips a single polygon to the limits described in `limits'.
  
  Notes:
    - limits is a 4-element array of doubles (left, right, bottom, top)
    - polygons is a boolean value specifying whether to clip as polygons or lines
    - returns value in out*
    - recommended allocated space for out*: 2 x *inVerts
  
  Returns:
    *outVerts >= 0 on success
    *outVerts = -1 on insufficient memory
    *outVerts = -2 on insufficient space in out*
---------------------------------------------------------------------------*/
void clipPolygon(double *inX, double *inY, PBSINT *inPOS, PBSINT inVerts,
                 double *outX, double *outY, PBSINT *outOLD, PBSINT *outVerts,
                 double *limits, short polygons);

/*-----------------------------------------------------------------------------
  isPolyConvex:
    Determines whether a polygon is convex.

  Author:  Nicholas Boers (June 30, 2004)
  
  Notes:
    - checks for self-intersections
    - polygons consisting of one or two points cannot possibly be convex
  ---------------------------------------------------------------------------*/
short isPolyConvex(double *inX, double *inY, PBSINT inVerts);

/*-----------------------------------------------------------------------------
  nPolyIntersects:
    Counts the number of times a polygon self-intersects itself.

  Author:  Nicholas Boers (June 30, 2004)
  
  Notes:
    - counts certain types of intersections (i.e., those involving vertices
      and those where an edge retraces over an edge) more than once
    - returns the number of times a polygon intersects itself;
      if !numericResult, stops counting at 1 intersection
    - if a polygon only has two points, it's self-intersecting
    - if it only has one point, it isn't self-intersecting
  ---------------------------------------------------------------------------*/
PBSINT nPolyIntersects(double *inX, double *inY, PBSINT inVerts,
                     short numericResult);


/*-----------------------------------------------------------------------------
  thickenPoly:
    Thickens a polygon (adds additional points to it).

  Note:
    - if units == 0, "LL": tolerance in kilometers and inXY in decimal-degrees
    - if units == 1, other: tolerance and inXY in same units
    - there is one special case when !keepOrig (commented in code, before the
      loops)
    - could improve accuracy for LL by changing the calculation that determines
      intermediate LL points (should the point be on an arc?)

  Returns:
    -1 if out* is too small
     # of output vertices on success
  ---------------------------------------------------------------------------*/
PBSINT thickenPoly(double *inX, double *inY, PBSINT inN,
                 double *outX, double *outY, PBSINT outN,
                 double tol, short units, short keepOrig, short close);

#endif
