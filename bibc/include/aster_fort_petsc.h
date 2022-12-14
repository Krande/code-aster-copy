/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org             */
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
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ASTER_HAVE_PETSC
#define CALLO_MATASS2PETSC( a, b, c ) CALLOPP( MATASS2PETSC, matass2petsc, a, b, c )
#define CALL_MATASS2PETSC( a, b, c ) CALLSPP( MATASS2PETSC, matass2petsc, a, b, c )
void DEFSPP( MATASS2PETSC, matass2petsc, const char *, STRING_SIZE, Mat *, PetscErrorCode * );

#define CALLO_AP_ON_OFF( a, b ) CALLOO( AP_ON_OFF, ap_on_off, a , b)
void DEFSS( AP_ON_OFF, ap_on_off, const char *, STRING_SIZE , const char *, STRING_SIZE );
#endif

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_PETSC_H */
#endif
