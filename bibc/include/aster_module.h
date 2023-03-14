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

/* person_in_charge: mathieu.courtois at edf.fr */

#ifndef ASTER_MODULE_H_
#define ASTER_MODULE_H_

#include "aster.h"
/*
 *   PUBLIC FUNCTIONS
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __FILE_NAME__
#define PATH __FILE_NAME__
#else
#define PATH __FILE__
#endif

extern PyObject *PyInit_aster( void );

extern void PRE_myabort( _IN const char *nomFichier, _IN const int numeroLigne,
                         _IN const char *message );
#define MYABORT( message ) PRE_myabort( PATH, __LINE__, message )

#define CALL_ISJVUP() CALL0( ISJVUP, isjvup )
extern ASTERINTEGER DEF0( ISJVUP, isjvup );

extern void CALLP( XFINI, xfini, _IN ASTERINTEGER * );

extern void getvpy( const char *, const char *, const int, PyObject ** );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_MODULE_H */
#endif
