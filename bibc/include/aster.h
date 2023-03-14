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

#ifndef ASTER_H_
#define ASTER_H_

/* NOTE:
 *  Since Python may define some pre-processor definitions which affect the
 *  standard headers on some systems, you must include "Python.h" before any
 *  standard headers are included.
 *  The warning on _POSIX_C_SOURCE redefinition must not occur.
 *
 *  source: http://docs.python.org/c-api/intro.html
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "aster_depend.h"
#include "asterc_debug.h"
#include "definition.h"

#include <signal.h>
#include <stdio.h>

/* pour indiquer le statut des arguments des fonctions. */

#define _IN
#define _OUT
#define _INOUT
#define _UNUSED

/* AS_ASSERT is similar to the ASSERT macro used in fortran */
#define AS_ASSERT( cond )                                                                          \
    if ( !( cond ) ) {                                                                             \
        DEBUG_LOC;                                                                                 \
        DBGV( "Assertion failed: %s", #cond );                                                     \
        INTERRUPT( 17 );                                                                           \
    }

/* deprecated functions on Windows */
#ifdef ASTER_PLATFORM_WINDOWS
#define strdup _strdup
#endif

#endif
