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

#ifndef ASTERC_DEBUG_H_
#define ASTERC_DEBUG_H_

#include "aster.h"

#include <stdio.h>
#include <stdlib.h>

/*! Here are defined some flags to add debugging informations.

If the flag is defined, a function prints informations on stdout.
If the flag is not defined, the function must be empty macro.

to enable DEBUG_ASSERT
#define ASTER_DEBUG_ASSERT

to add all traces
#define ASTER_DEBUG_ALL

to debug MPI (add to asterf_debug.h too)
#define ASTER_DEBUG_MPI
*/

#ifdef __FILE_NAME__
#define PATH __FILE_NAME__
#else
#define PATH __FILE__
#endif

/*! print the filename and line where the error occurred */
#define DEBUG_LOC                                                                                  \
    fprintf( stdout, "DEBUG: %s #%d: ", PATH, __LINE__ );                                          \
    fflush( stdout );

/*! interrupt the execution, return SIGABRT */
#define INTERRUPT( code )                                                                          \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        fprintf( stdout, "ABORT - exit code %d\n", code );                                         \
        fflush( stdout );                                                                          \
        abort();                                                                                   \
    }

/*! internal utility to print a PyObject */
#define PYDBG( label, pyobj )                                                                      \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        fprintf( stdout, label );                                                                  \
        PyObject_Print( pyobj, stdout, 0 );                                                        \
        printf( "\n" );                                                                            \
        fflush( stdout );                                                                          \
    }

#define DBG( label )                                                                               \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        printf( label );                                                                           \
        printf( "\n" );                                                                            \
        fflush( stdout );                                                                          \
    }
#define DBGV( fmt, a )                                                                             \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        printf( fmt, a );                                                                          \
        printf( "\n" );                                                                            \
        fflush( stdout );                                                                          \
    }
#define DBGVV( fmt, a, b )                                                                         \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        printf( fmt, a, b );                                                                       \
        printf( "\n" );                                                                            \
        fflush( stdout );                                                                          \
    }

/*! enable DEBUG_ASSERT */
#if defined( ASTER_DEBUG_ASSERT ) || defined( ASTER_DEBUG_ALL )
#define DEBUG_ASSERT( cond ) AS_ASSERT( cond )
#else
#define DEBUG_ASSERT( cond )
#endif

/*! enable DEBUG_DLL */
#if defined( ASTER_DEBUG_DLL ) || defined( ASTER_DEBUG_ALL )
#define DEBUG_DLL_PYOB( label, pyobj ) PYDBG( label, pyobj )
#define DEBUG_DLL_VV( fmt, a, b ) DBGVV( fmt, a, b )
#else
#define DEBUG_DLL_PYOB( label, pyobj )
#define DEBUG_DLL_VV( fmt, a, b )
#endif

/*! enable DEBUG_EXCEPT, not in ASTER_DEBUG_ALL */
#if defined( ASTER_DEBUG_EXCEPT )
#define DEBUG_EXCEPT( fmt, a ) DBGV( fmt, a )
#else
#define DEBUG_EXCEPT( fmt, a )
#endif

/*! debug MPI communicator as aster_comm_t */
#if defined( ASTER_DEBUG_MPICOM ) || defined( ASTER_DEBUG_ALL )
#define COMM_DEBUG( c )                                                                            \
    {                                                                                              \
        DEBUG_LOC;                                                                                 \
        printf( "%-8s #%d (%d/@", ( c ).name, (int)MPI_Comm_c2f( ( c ).id ), ( c ).level );        \
        if ( ( c ).parent ) {                                                                      \
            printf( "%-8s", ( c ).parent->name );                                                  \
        } else {                                                                                   \
            printf( "        " );                                                                  \
        }                                                                                          \
        printf( ")\n" );                                                                           \
        fflush( stdout );                                                                          \
    }
#else
#define COMM_DEBUG( c )
#endif

/*! debug MPI communications */
#if defined( ASTER_DEBUG_MPI ) || defined( ASTER_DEBUG_ALL )
#define DEBUG_MPI( fmt, a, b ) DBGVV( fmt, a, b )
#else
#define DEBUG_MPI( fmt, a, b )
#endif

/*! enable DEBUG_ASTER_FONCTIONS */
#if defined( ASTER_DEBUG_ALL )
#define ASTER_DEBUG_FONCTIONS
#endif

/*! enable DEBUG_IODR */
#if defined( ASTER_DEBUG_IODR ) || defined( ASTER_DEBUG_ALL )
#define DEBUG_IODR( fmt, a, b ) DBGVV( fmt, a, b )
#else
#define DEBUG_IODR( fmt, a, b )
#endif

#endif
