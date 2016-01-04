/*------------------------------------------------------------------------------
 * File:   gshhs.c
 *
 * This code is slightly modified from the original code from Paul Wessel
 * that was downloaded from 
 *   http://www.soest.hawaii.edu/wessel/gshhs/
 * and made available under GNU GPL.
 *
 * Modified by Nicholas Boers.
 *
 * Original header follows.
 * -----------------------------------------------------------------------------
 *	$Id: gshhs.c 9656 2012-01-13 18:29:02Z pwessel $
 *
 *	Copyright (c) 1996-2013 by P. Wessel and W. H. F. Smith
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 * PROGRAM:	gshhs.c
 * AUTHOR:	Paul Wessel (pwessel@hawaii.edu)
 * CREATED:	JAN. 28, 1996
 * PURPOSE:	To extract ASCII data from the binary GSHHS shoreline data
 *		as described in the 1996 Wessel & Smith JGR Data Analysis Note.
 * VERSION:	1.1 (Byte flipping added)
 *		1.2 18-MAY-1999:
 *		   Explicit binary open for DOS systems
 *		   POSIX.1 compliant
 *		1.3 08-NOV-1999: Released under GNU GPL
 *		1.4 05-SEPT-2000: Made a GMT supplement; FLIP no longer needed
 *		1.5 14-SEPT-2004: Updated to deal with latest GSHHS database (1.3)
 *		1.6 02-MAY-2006: Updated to deal with latest GSHHS database (1.4)
 *		1.7 11-NOV-2006: Fixed bug in computing level (&& vs &)
 *		1.8 31-MAR-2007: Updated to deal with latest GSHHS database (1.5)
 *		1.9 27-AUG-2007: Handle line data as well as polygon data
 *		1.10 15-FEB-2008: Updated to deal with latest GSHHS database (1.6)
 *		1.11 15-JUN-2009: Now contains information on container polygon,
 *				the polygons ancestor in the full resolution, and
 *				a flag to tell if a lake is a riverlake.
 *				Updated to deal with latest GSHHS database (2.0)
 *		1.12 24-MAY-2010: Deal with 2.1 format.
 *		1.13 15-JUL-2011: Now contains improved area information (2.2.0),
 *				 and revised greenwhich flags (now 2-bit; see gshhs.h).
 *				 Also added -A and -G as suggested by José Luis García Pallero,
 *				 as well as -Qe|i to control river-lake output, and -N to
 *				 get a particular level.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: www.soest.hawaii.edu/pwessel */

#include "gshhs.h"

#ifndef STANDALONE
/* when running through a call to 'gshhs' rather than as a stand-alone application,
   the code will call the following functions rather than printing to standard
   output */
void (*polygonHeader) (char c, int h_id, int h_n, int level, char source, double area, double f_area,
		       double w, double e, double s, double n, int container, int ancestor) = NULL;
void (*lineHeader) (char c, int h_id, int h_n, int level, char source, double w, double e, double s,
		    double n) = NULL;
void (*point) (double lon, double lat) = NULL;

/* return rather than exit */
#define exit			return
/* effectively disable printing to stderr */
#define fprintf(...)
#endif /* STANDALONE */

#ifdef STANDALONE
int main (int argc, char **argv)
#else /* !STANDALONE */
int gshhs (int argc, char **argv)
#endif /* !STANDALONE */
{
	double w, e, s, n, area, f_area, lon, lat, scale = 10.0, minarea = 0.0;
	char source, kind[2] = {'P', 'L'}, c = '>', *file = NULL;
#ifdef STANDALONE
	char *name[2] = {"polygon", "line"}, container[8], ancestor[8];
#endif /* STANDALONE */
	FILE *fp = NULL;
	int k, line, max_east = 270000000, info, single, error, ID, flip, m, octave = 0, rank = -1;
	int  OK, level, version, greenwich, river, src, msformat = 0, qmode = 0, first = 1;
	size_t n_read;
	struct GSHHS h;
 	struct POINT p;
       
	info = single = error = ID = 0;
	for (k = 1; k < argc; k++) {
		if (argv[k][0] == '-') {	/* Option switch */
			switch (argv[k][1]) {
				case 'A':
					minarea = atof (&argv[k][2]);
					break;
				case 'G':
					octave = 1;
					break;
				case 'L':
					info = 1;
					break;
				case 'M':
					msformat = 1;
					break;
				case 'N':
					rank = atoi (&argv[k][2]);
					break;
				case 'I':
					if (argv[k][2] == 'c')
						single = 2;
					else {
						single = 1;
						ID = atoi (&argv[k][2]);
					}
					break;
				case 'Q':
					if (argv[k][2] == 'e')
						qmode = 1;
					else if (argv[k][2] == 'i')
						qmode = 2;
					else
						error++;
					break;
				default:
					error++;
					break;
			}
		}
		else
			file = argv[k];
	}

	if (!file) error++;
	
	if (argc == 1 || error) {
		fprintf (stderr, "gshhs %s - ASCII export of GSHHS/WDBII %s data\n\n", GSHHS_PROG_VERSION, GSHHS_DATA_VERSION);
		fprintf (stderr, "usage:  gshhs gshhs|wdbii_rivers|wdb_borders_[f|h|i|l|c].b [-A<area>] [-G] [-I<id>] [-L] [-M] [-N<level>] [-Qe|i] > ascii.txt\n\n");
		fprintf (stderr, "\t-A will extract only polygons whose area is greater or equal than <area>.\n");
		fprintf (stderr, "\t-G will write '%%' at start of each segment header [P or L] (overwrites -M)\n"
				 "\t   and write 'NaN NaN' after each segment to enable import by GNU Octave or Matlab.\n");
		fprintf (stderr, "\t-I Output data for polygon number <id>.  Use -Ic to get all continent polygons\n");
		fprintf (stderr, "\t   [Default is all polygons].\n");
		fprintf (stderr, "\t-L will only list headers (no data output).\n");
		fprintf (stderr, "\t-M will write '>' at start of each segment header [P or L].\n");
		fprintf (stderr, "\t-N will only output features whose level matches <level> [Default is all levels].\n");
		fprintf (stderr, "\t-Q controls river-lakes.  By default all polygons are output.\n");
		fprintf (stderr, "\t   use -Qe to exclude river-lakes, and -Qi to ONLY get river-lakes.\n");
		exit (EXIT_FAILURE);
	}
	if (!file) {
		fprintf (stderr, "gshhs: No data file given!\n");
		error++;
	}

	if ((fp = fopen (file, "rb")) == NULL ) {
		fprintf (stderr, "gshhs:  Could not find file %s.\n", file);
		exit (EXIT_FAILURE);
	}
		
	n_read = fread ((void *)&h, (size_t)sizeof (struct GSHHS), (size_t)1, fp);
	version = (h.flag >> 8) & 255;
	flip = (version != GSHHS_DATA_RELEASE);	/* Take as sign that byte-swabbing is needed */
	
	while (n_read == 1) {
		if (flip) {
			h.id = swabi4 ((unsigned int)h.id);
			h.n  = swabi4 ((unsigned int)h.n);
			h.west  = swabi4 ((unsigned int)h.west);
			h.east  = swabi4 ((unsigned int)h.east);
			h.south = swabi4 ((unsigned int)h.south);
			h.north = swabi4 ((unsigned int)h.north);
			h.area  = swabi4 ((unsigned int)h.area);
			h.area_full  = swabi4 ((unsigned int)h.area_full);
			h.flag  = swabi4 ((unsigned int)h.flag);
			h.container = swabi4 ((unsigned int)h.container);
			h.ancestor  = swabi4 ((unsigned int)h.ancestor);
		}
		level = h.flag & 255;				/* Level is 1-4 */
		version = (h.flag >> 8) & 255;			/* Version is 1-255 */
		if (first) fprintf (stderr, "gshhs %s - Found GSHHS/WDBII version %d in file %s\n", GSHHS_PROG_VERSION, version, file);
		greenwich = (h.flag >> 16) & 3;			/* Greenwich is 0-3 */
		src = (h.flag >> 24) & 1;			/* Greenwich is 0 (WDBII) or 1 (WVS) */
		river = (h.flag >> 25) & 1;			/* River is 0 (not river) or 1 (is river) */
		w = h.west  * GSHHS_SCL;			/* Convert from microdegrees to degrees */
		e = h.east  * GSHHS_SCL;
		s = h.south * GSHHS_SCL;
		n = h.north * GSHHS_SCL;
		source = (src == 1) ? 'W' : 'C';		/* Either WVS or CIA (WDBII) pedigree */
		if (river) source = tolower ((int)source);	/* Lower case c means river-lake */
		line = (h.area) ? 0 : 1;			/* Either Polygon (0) or Line (1) (if no area) */
		if (version >= 9) {				/* Variable magnitude for area scale */
			m = h.flag >> 26;
			scale = pow (10.0, (double)m);		/* Area scale */
		}
		area = h.area / scale;				/* Now im km^2 */
		f_area = h.area_full / scale;			/* Now im km^2 */

		OK = ((!single || ((single == 1 && h.id == ID) || (single == 2 && h.id <= 5))) && area >= minarea);	/* Skip if not the one (-I) or too small (-A) */
		if (OK && qmode && ((river && qmode == 1) || (!river && qmode == 2))) OK = 0;	/* Skip if riverlake/not riverlake (-Q) */
		if (OK && rank >= 0 && level != rank) OK = 0;		/* Skip if not the right level (-N) */
		first = 0;
		
		if (!msformat) c = kind[line];
	        if (octave) c = '%';
		if (OK) {
			if (line) {
#ifdef STANDALONE
				printf ("%c %6d%8d%3d%2c%10.5f%10.5f%10.5f%10.5f\n", c, h.id, h.n, level, source, w, e, s, n);
#else /* !STANDALONE */
				if (lineHeader)
					(*lineHeader)(c, h.id, h.n, level, source, w, e, s, n);
#endif /* !STANDALONE */
			} else {
#ifdef STANDALONE
				(h.container == -1) ? sprintf (container, "-") : sprintf (container, "%6d", h.container);
				(h.ancestor == -1) ? sprintf (ancestor, "-") : sprintf (ancestor, "%6d", h.ancestor);
				printf ("%c %6d%8d%2d%2c %.12g %.12g%11.5f%11.5f%10.5f%10.5f %s %s\n", c, h.id, h.n, level, source, area, f_area, w, e, s, n, container, ancestor);
#else /* !STANDALONE */
				if (polygonHeader)
					(*polygonHeader)(c, h.id, h.n, level, source, area, f_area, w, e, s, n, h.container, h.ancestor);
#endif /* !STANDALONE */
			}
		}

		if (info || !OK) {	/* Skip data, only want headers */
			fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
		}
		else {
			for (k = 0; k < h.n; k++) {

				if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
					fprintf (stderr, "gshhs:  Error reading file %s for %s %d, point %d.\n", argv[1], name[line], h.id, k);
					exit (EXIT_FAILURE);
				}
				if (flip) {
					p.x = swabi4 ((unsigned int)p.x);
					p.y = swabi4 ((unsigned int)p.y);
				}
				lon = p.x * GSHHS_SCL;
				if (((greenwich & 1) && p.x > max_east) || (h.west > 180000000)) lon -= 360.0;
				lat = p.y * GSHHS_SCL;
#ifdef STANDALONE
				printf ("%11.6f%11.6f\n", lon, lat);
#else /* !STANDALONE */
				if (point)
					(*point) (lon, lat);
#endif /* !STANDALONE */
			}
#ifdef STANDALONE
			if (octave) printf ("NaN\tNaN\n");
#endif /* STANDALONE */
		}
		max_east = 180000000;	/* Only Eurasia needs 270 */
		n_read = fread((void *)&h, (size_t)sizeof (struct GSHHS), (size_t)1, fp);
	}
		
	fclose (fp);

	exit (EXIT_SUCCESS);
}
