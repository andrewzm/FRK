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
  File: convGSHHS.c

  Interface between PBS Mapping and GSHHS.  This file defines a number of
  callbacks that are called by our slightly modified version of GSHHS. 
  Although this approach isn't as efficient as directly modifying the
  GSHHS source, it should allow us to easily adapt to new versions of GSHHS.

  Given a new version of the GSHHS source, search through the current version
  for the STANDALONE defines and apply those changes to the GSHHS source.
  ---------------------------------------------------------------------------*/
#include <R.h>
#include <Rdefines.h>
#include <errno.h>

#include "globals.h"
#include "floating.h"
#include "polygons.h"
#include "gshhs.h"

#ifndef FALSE
#define FALSE	0
#endif /* FALSE */

#ifndef TRUE
#define TRUE	1
#endif /* TRUE */

/* number of points/polygons/lines to extract (calculated during phase 1) */
static int numPoints, numPolys, numLines;

/* arguments supplied with the call from R */
static double xlim[2], ylim[2];
static int maximumLevel, minimumVerts;

/* current state of the extraction */
static int pid = -1, sid = -1, pos = -1;
static int extractCur = FALSE;	/* Boolean; extracting current polygon? */

/* vectors for extracted PolySet and PolyData */
static SEXP psetPID, psetSID, psetPOS, psetX, psetY;
static SEXP pdatPID, pdatSID, pdatLevel, pdatSource;

/* if a polygon is (a) outside of the limits, (b) has too few vertices,
   or (c) is at too high of a level, skip it */
#define SET_POLY_INCLUSION	do {					\
	if ((e < xlim[0] || w > xlim[1] || n < ylim[0] || s > ylim[1])	\
	    || h_n < minimumVerts || level > maximumLevel)		\
	    extractCur = FALSE;						\
	else								\
	    extractCur = TRUE;						\
    } while (0)

/* if a line is (a) outside of the limits or (b) has too few vertices,
   skip it */
#define SET_LINE_INCLUSION	do {					\
	if ((e < xlim[0] || w > xlim[1] || n < ylim[0] || s > ylim[1])	\
	    || h_n < minimumVerts)					\
	    extractCur = FALSE;						\
	else								\
	    extractCur = TRUE;						\
    } while (0)

/* polyCountPoints(...): Callback for gshhs; prior to extracting, this
   routine is used to determine the number of points that we'll need to
   extract.
  
   Arguments:
     c		kind {P,L} for polygon or line, respectively
     h_id	identifier
     n		number of points
     level	1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
     source	0 = CIA WDBII, 1 = WVS
     area	area of polygon in 1/10 km^2
     f_area	area of original full-resolution polygon in 1/10 km^2
     w,e,s,n	min/max extent
     container	identifier of parent (at level-1)
     ancestor	polygon's full-resolution ancestor
 */
void polyCountPoints (char c, int  h_id, int h_n, int level, char source,
		      double area, double f_area, double w, double e,
		      double s, double n, int container, int ancestor)
{
    SET_POLY_INCLUSION;

    if (extractCur) {
	numPolys++;
	numPoints += h_n;
    }
}

/* polyExtract(...): Callback for gshhs; when extracting polygons,
   this routine is used to set the PID/SID/POS field appropriately.
  
   The level indicates forward or reverse pos:
     if level 1, PID = h_id, SID = 0, POS = normal
     if level 2, PID = container, SID = h_id, POS = reverse
     if level 3, PID = h_id, SID = 0, POS = normal
     if level 4, PID = container, SID = h_id, POS = reverse
*/
void polyExtract (char c, int  h_id, int h_n, int level, char source,
		  double area, double f_area, double w, double e,
		  double s, double n, int container, int ancestor)
{
    SET_POLY_INCLUSION;

    if (extractCur) {
	if (level == 2 || level == 4) {
	    /* a hole */
	    pid = container;
	    sid = h_id;		/* non-zero => decreasing pos */
	    pos = h_n;
	} else {
	    /* an outer boundary */
	    pid = h_id;
	    sid = 0;		/* 0 => increasing pos */
	    pos = 0;
	}

	/* save the polygon in the PolyData */
	INTEGER_POINTER(pdatPID)[numPolys] = pid;
	INTEGER_POINTER(pdatSID)[numPolys] = sid;
	INTEGER_POINTER(pdatLevel)[numPolys] = level;
	INTEGER_POINTER(pdatSource)[numPolys] = (source == 'W' ? 1 : 0);

	numPolys++;
    }
}

/* lineCountPoints(...): Callback for gshhs; similar to one for polys. */
void lineCountPoints (char c, int h_id, int h_n, int level, char source, 
		      double w, double e, double s, double n)
{
    SET_LINE_INCLUSION;

    if (extractCur) {
	numLines++;
	numPoints += h_n;
    }
}

/* lineExtract(...): Callback for gshhs; similar to previous, but for
   lines rather than polygons. */
void lineExtract (char c, int h_id, int h_n, int level, char source, 
		  double w, double e, double s, double n)
{
    SET_LINE_INCLUSION;

    if (extractCur) {
	pid = h_id;
	sid = 0;		/* 0 => increasing pos */
	pos = 0;

	/* save the line in the PolyData */
	INTEGER_POINTER(pdatPID)[numLines] = pid;
	INTEGER_POINTER(pdatSID)[numLines] = sid;
	INTEGER_POINTER(pdatLevel)[numLines] = level;
	INTEGER_POINTER(pdatSource)[numLines] = (source == 'W' ? 1 : 0);

	numLines++;
    }
}

/* extractPoint(...): Callback for gshhs; when extracting lines/polygons,
   this routine handles the extraction of a single point.
*/
void extractPoint (double lon, double lat)
{
    if (extractCur) {
	/* save the point in the PolySet */
	INTEGER_POINTER(psetPID)[numPoints] = pid;
	INTEGER_POINTER(psetSID)[numPoints] = sid;
	INTEGER_POINTER(psetPOS)[numPoints] = pos;
	NUMERIC_POINTER(psetX)[numPoints] = lon;
	NUMERIC_POINTER(psetY)[numPoints] = lat;
	
	pos += (sid ? -1 : 1);	/* 0 => increasing; non-zero => decreasing */
	numPoints++;
    }
}

/* importGSHHS(...): Routine called by R to extract a region.
   Arguments:
     gshhsFileName	file name (string)
     clipLimits		longitude/latitude limits for extraction
     levels		maximum level to extract
     minVerts		minimum vertices in an extracted polygon
*/
SEXP importGSHHS(SEXP gshhsFileName, SEXP clipLimits, SEXP levels,
		 SEXP minVerts)
{
    char *cmdList[] = {"", "-L", "" };	/* for calling gshhs to list polygons */
    char *cmdExtract[] = {"", "" };	/* for calling gshhs to extract 
					   polygons */

    SEXP polySet, polyData, temp;

    /* obtain arguments */
    cmdList[2] = cmdExtract[1] = (char *)CHAR(STRING_ELT(gshhsFileName,0));
    xlim[0] = NUMERIC_POINTER(clipLimits)[0];
    xlim[1] = NUMERIC_POINTER(clipLimits)[1];
    ylim[0] = NUMERIC_POINTER(clipLimits)[2];
    ylim[1] = NUMERIC_POINTER(clipLimits)[3];
    maximumLevel = INTEGER_POINTER(levels)[0];
    minimumVerts = INTEGER_POINTER(minVerts)[0];

    /* calculate the number of points to extract */
    Rprintf ("importGSHHS status:\n");
    Rprintf ("--> Pass 1: ");
    numPoints = numPolys = numLines = 0;
    polygonHeader = &polyCountPoints;
    lineHeader = &lineCountPoints;
    point = NULL;
    if (gshhs (3, cmdList) != EXIT_SUCCESS)
	error("call to gshhs failed.\n");
    if (numPolys && numLines)
	error ("encountered both polygons and lines: unsupported.\n");
    Rprintf ("complete: %d bounding boxes within limits.\n",
	     numPolys + numLines);

    /* set up the polySet that we'll return */
    PROTECT(polySet = allocVector(VECSXP, 5));

    /* add attribute to indicate whether we later need to clip */
    PROTECT(temp = allocVector(INTSXP, 1));
    INTEGER_POINTER(temp)[0] = (numPolys ? 1 : 0);    
    setAttrib(polySet, install("clipAsPolys"), temp);
    UNPROTECT(1);

    /* add column names to polySet */
    PROTECT(temp = allocVector(STRSXP, 5));
    SET_STRING_ELT(temp, 0, mkChar("PID"));
    SET_STRING_ELT(temp, 1, mkChar("SID"));
    SET_STRING_ELT(temp, 2, mkChar("POS"));
    SET_STRING_ELT(temp, 3, mkChar("X"));
    SET_STRING_ELT(temp, 4, mkChar("Y"));
    setAttrib(polySet, R_NamesSymbol, temp);
    UNPROTECT(1); // temp now protected by polySet

    /* add columns to polySet */
    PROTECT(psetPID = allocVector(INTSXP, numPoints));
    PROTECT(psetSID = allocVector(INTSXP, numPoints));
    PROTECT(psetPOS = allocVector(INTSXP, numPoints));
    PROTECT(psetX = allocVector(REALSXP, numPoints));
    PROTECT(psetY = allocVector(REALSXP, numPoints));
    SET_VECTOR_ELT(polySet, 0, psetPID);
    SET_VECTOR_ELT(polySet, 1, psetSID);
    SET_VECTOR_ELT(polySet, 2, psetPOS);
    SET_VECTOR_ELT(polySet, 3, psetX);
    SET_VECTOR_ELT(polySet, 4, psetY);
    UNPROTECT(5); // vectors now protected by polySet

    /* set up the PolyData that we'll return */
    PROTECT(polyData = allocVector(VECSXP, 4));

    /* add column names to polyData */
    PROTECT(temp = allocVector(STRSXP, 4));
    SET_STRING_ELT(temp, 0, mkChar("PID"));
    SET_STRING_ELT(temp, 1, mkChar("SID"));
    SET_STRING_ELT(temp, 2, mkChar("Level"));
    SET_STRING_ELT(temp, 3, mkChar("Source"));
    setAttrib(polyData, R_NamesSymbol, temp);
    UNPROTECT(1); // temp now protected by polyData
    
    /* add columns to PolyData */
    PROTECT(pdatPID = allocVector(INTSXP, numPolys + numLines));
    PROTECT(pdatSID = allocVector(INTSXP, numPolys + numLines));
    PROTECT(pdatLevel = allocVector(INTSXP, numPolys + numLines));
    PROTECT(pdatSource = allocVector(INTSXP, numPolys + numLines));
    SET_VECTOR_ELT(polyData, 0, pdatPID);
    SET_VECTOR_ELT(polyData, 1, pdatSID);
    SET_VECTOR_ELT(polyData, 2, pdatLevel);
    SET_VECTOR_ELT(polyData, 3, pdatSource);
    UNPROTECT(4); // vectors now protected by polyData

    /* attach polyData to polySet as "PolyData" attribute */
    setAttrib(polySet, install("PolyData"), polyData);
    UNPROTECT(1); // polyData now protected by polySet

    /* extract the polygons */
    Rprintf ("--> Pass 2: ");
    numPoints = numPolys = numLines = 0;
    polygonHeader = &polyExtract;
    lineHeader = &lineExtract;
    point = &extractPoint;
    if (gshhs (2, cmdExtract) != EXIT_SUCCESS)
	error("call to gshhs failed.\n");
    Rprintf ("complete.\n");
    Rprintf ("--> Clipping...\n");

    UNPROTECT(1); // polySet is the last one...

    return (polySet);
}
