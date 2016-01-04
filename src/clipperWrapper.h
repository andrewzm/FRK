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
  File: clipperWrapper.hpp

  Interface between R and the Clipper library.

  Author: Nicholas Boers
  ---------------------------------------------------------------------------*/
#ifndef _CLIPPERWRAPPER_H_
#define _CLIPPERWRAPPER_H_

#ifndef STANDALONE

#include <R.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

SEXP clipperWrapper (SEXP operation,
		     SEXP sPID, SEXP sSID, SEXP sPOS, SEXP sX, SEXP sY,
		     SEXP cPID, SEXP cSID, SEXP cPOS, SEXP cX, SEXP cY);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* not defined STANDALONE */

#endif /* _CLIPPERWRAPPER_H_ */
