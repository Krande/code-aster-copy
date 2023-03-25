/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#ifndef ASTER_FORT_BLAS_H_
#define ASTER_FORT_BLAS_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif
/*
// #define CALL_DAXPY( a, b, c, e, f, g ) CALLPSSSSSS( DAXPY, daxpy, a, b, c, e, f, g )
// #define CALLO_DAXPY( a, b, c, e, f, g ) CALLPOOOOOO( DAXPY, daxpy, a, b, c, e, f, g )
// extern void DEFPSSSSSS( DAXPY, daxpy, ASTERINTEGER *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERINTEGER
// *,
//                        ASTERDOUBLE *, ASTERINTEGER * );
*/

#define CALL_DSCAL( a, b, c, d ) CALLPPPP( DSCAL, dscal, a, b, c, d )
extern void DEFPPPP( DSCAL, dscal, const ASTERINTEGER *, const ASTERDOUBLE *, ASTERDOUBLE *,
                     const ASTERINTEGER * );

#define CALL_ZSCAL( a, b, c, d ) CALLPPPP( ZSCAL, zscal, a, b, c, d )
extern void DEFPPPP( ZSCAL, zscal, const ASTERINTEGER *, const ASTERCOMPLEX *, ASTERCOMPLEX *,
                     const ASTERINTEGER * );

#define CALL_ZDSCAL( a, b, c, d ) CALLPPPP( ZDSCAL, zdscal, a, b, c, d )
extern void DEFPPPP( ZDSCAL, zdscal, const ASTERINTEGER *, const ASTERDOUBLE *, ASTERCOMPLEX *,
                     const ASTERINTEGER * );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_BLAS_H_ */
#endif
