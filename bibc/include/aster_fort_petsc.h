/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
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

#ifndef ASTER_FORT_PETSC_H_
#define ASTER_FORT_PETSC_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef ASTER_HAVE_PETSC
#include "petscmat.h"
#include "petscvec.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ASTER_HAVE_PETSC
#define CALLO_MATASS2PETSC( a, b, c, d ) CALLOPPP( MATASS2PETSC, matass2petsc, a, b, c, d )
#define CALL_MATASS2PETSC( a, b, c, d ) CALLSPPP( MATASS2PETSC, matass2petsc, a, b, c, d )
void DEFSPPP( MATASS2PETSC, matass2petsc, const char *, STRING_SIZE, const ASTERINTEGER *, Mat *,
              PetscErrorCode * );

#define CALLO_AP_ON_OFF( a, b ) CALLOO( AP_ON_OFF, ap_on_off, a, b )
void DEFSS( AP_ON_OFF, ap_on_off, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VECT_ASSE_FROM_PETSC( a, b, c, d, e )                                                \
    CALLOOPPP( VECT_ASSE_FROM_PETSC, vect_asse_from_petsc, a, b, c, d, e )
extern void DEFSSPPP( VECT_ASSE_FROM_PETSC, vect_asse_from_petsc, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const Vec *, const ASTERDOUBLE *,
                      const ASTERINTEGER * );

#endif

#define CALLO_VECT_ASSE_UPDATE_GHOST_VALUES( a, b )                                                \
    CALLOO( VECT_ASSE_UPDATE_GHOST_VALUES, vect_asse_update_ghost_values, a, b )
extern void DEFSS( VECT_ASSE_UPDATE_GHOST_VALUES, vect_asse_update_ghost_values, const char *,
                   STRING_SIZE, const char *, STRING_SIZE );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_PETSC_H */
#endif
