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
  File: floating.h
   
  #defines and macros for floating point calculations.

  History:
    ?? ??? ???? [Nicholas Boers]
      - original file

    17 Jun 2004 [Nicholas Boers]
      - cleaned up
  ---------------------------------------------------------------------------*/

#ifndef _FLOATING_H
#define _FLOATING_H

#include <math.h>
#include <float.h>

#define EPSILON         (DBL_EPSILON)

/* for checking the equality of floating point numbers */
#define FABS(n)                                                          \
                (((n) < 0) ? -(n) : (n))

#define DBL_EQ(n1, n2)                                                   \
                (((n1) == 0 && (n2) == 0)                                \
                 || (((n1) != 0)                                         \
                     && FABS((n1) - (n2))/FABS((n1)) <= DBL_EPSILON)     \
                 || FABS((n1) - (n2)) <= DBL_EPSILON)
#define DBL_LTEQ(n1, n2)                                                 \
                (((n1) < (n2)) || DBL_EQ((n1), (n2)))
#define DBL_GTEQ(n1, n2)                                                 \
                (((n1) > (n2)) || DBL_EQ((n1), (n2)))
#define DBL_LT(n1, n2)                                                   \
                (((n1) < (n2)) && !DBL_EQ((n1), (n2)))
#define DBL_GT(n1, n2)                                                   \
                (((n1) > (n2)) && !DBL_EQ((n1), (n2)))

#endif /* _FLOATING_H */
