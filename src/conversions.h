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
  File: conversions.h
   
  Interfaces for converting units/projections.

  History:
    ?? ??? ???? [Nicholas Boers]
      - created

    17 Jun 2004 [Nicholas Boers]
      - cleaned up comments
  ---------------------------------------------------------------------------*/

#ifndef _CONVERSIONS_H
#define _CONVERSIONS_H

#define DEG_TO_RAD      (M_PI / 180.0)      /* degrees to radians */
#define RAD_TO_DEG      (180.0 / M_PI)      /* radians to degrees */

/* the return types from UTM<-->LL conversion routines */
struct pair
{
  double x;
  double y;
};

/*-----------------------------------------------------------------------------
  lonlat_to_utm:
    Convert a point from Lat/Lon to UTM. All angles are in radians.
    Memory should be allocated for the structure before the function
    is called. 

  Author:  Chris Grandin

  Algorithm Source:
    National Mapping Agency of Great Britain Ordnance Survey 
    <http://www.ordsvy.gov.uk>
  ---------------------------------------------------------------------------*/
void lonlat_to_utm(double lambda, double phi, int utmZone, struct pair *result);

/*-----------------------------------------------------------------------------
  utm_to_lonlat:
    Convert a point from UTM to Lat/Lon. All angles are in radians.
    Memory should be allocated for the structure before the function
    is called.

  Author:  Chris Grandin

  Algorithm Source:
    National Mapping Agency of Great Britain Ordnance Survey 
    <http://www.ordsvy.gov.uk>
  ---------------------------------------------------------------------------*/
void utm_to_lonlat(double E, double N, char hem, int utmZone, 
		   struct pair *result);

#endif /* _CONVERSIONS_H */
