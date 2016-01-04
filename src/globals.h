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
  File:   globals.h
   
  Constants global to many files.

  History:
    ?? ??? ???? [Nicholas Boers]
      - original file
    17 Jun 2004 [Nicholas Boers]
      - cleaned up
    19 Aug 2004 [Nicholas Boers]
      - changed Earth radius to 6371.3
    14 Nov 2004 [Nicholas Boers]
      - updated for PBSINT type
  ---------------------------------------------------------------------------*/

#ifndef _GLOBALS_H
#define _GLOBALS_H

#if defined(_SPLUS_) || defined(_STANDALONE_)
#define PBSINT          long
#else
#define PBSINT          int
#endif

#ifndef TRUE
#define TRUE            1
#endif /* TRUE */

#ifndef FALSE
#define FALSE           0
#endif /* FALSE */

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif /* M_PI */

/* Equatorial radius 6,378.14 km 
   Polar radius 6,356.78 km 
   Mean radius 6,371.3 km 

   Sources: 
     http://en.wikipedia.org/wiki/Earth
     http://en.wikipedia.org/wiki/Earth_radius */
#ifndef MEAN_RADIUS
#define MEAN_RADIUS     6371.3
#endif /* MEAN_RADIUS */

#endif /* _GLOBALS_H */
