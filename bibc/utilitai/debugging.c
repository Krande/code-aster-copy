/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org             */
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

#include "aster.h"

#ifndef ASTER_PLATFORM_MINGW
#include <execinfo.h>
#else
#include <windows.h>
#
#include <dbghelp.h>
#endif

#include <stdio.h>
#include <stdlib.h>

/*
 * This module defines a wrapper to call GNU libc backtrace functions.
 */
#define LEVEL 25

/* Obtain a backtrace and print it to stdout. */
void DEF0( PRINT_TRACE, print_trace ) {
#ifdef ASTER_HAVE_BACKTRACE

#ifndef ASTER_PLATFORM_MINGW
    void *array[LEVEL];
    size_t size;
    char **strings;
    size_t i;

    size = backtrace( array, LEVEL );
    strings = backtrace_symbols( array, size );

    fprintf( stderr, "Traceback returned by GNU libc (last %zd stack frames):\n", size );

    for ( i = 0; i < size; i++ )
        fprintf( stderr, "%s\n", strings[i] );

    free( strings );
#else
    HANDLE process = GetCurrentProcess();
    HANDLE thread = GetCurrentThread();

    CONTEXT context;
    memset( &context, 0, sizeof( CONTEXT ) );
    context.ContextFlags = CONTEXT_FULL;
    RtlCaptureContext( &context );

    SymInitialize( process, NULL, TRUE );

    DWORD image;
    STACKFRAME64 stackframe;
    ZeroMemory( &stackframe, sizeof( STACKFRAME64 ) );

    image = IMAGE_FILE_MACHINE_AMD64;
    stackframe.AddrPC.Offset = context.Rip;
    stackframe.AddrPC.Mode = AddrModeFlat;
    stackframe.AddrFrame.Offset = context.Rsp;
    stackframe.AddrFrame.Mode = AddrModeFlat;
    stackframe.AddrStack.Offset = context.Rsp;
    stackframe.AddrStack.Mode = AddrModeFlat;

    size_t i;
    printf( "Traceback returned by pdbhelp:\n" );
    for ( i = 0; i < LEVEL; i++ ) {

        BOOL result = StackWalk64( image, process, thread, &stackframe, &context, NULL,
                                   SymFunctionTableAccess64, SymGetModuleBase64, NULL );

        if ( !result ) {
            break;
        }

        char buffer[sizeof( SYMBOL_INFO ) + MAX_SYM_NAME * sizeof( TCHAR )];
        PSYMBOL_INFO symbol = (PSYMBOL_INFO)buffer;
        symbol->SizeOfStruct = sizeof( SYMBOL_INFO );
        symbol->MaxNameLen = MAX_SYM_NAME;

        DWORD64 displacement = 0;
        if ( SymFromAddr( process, stackframe.AddrPC.Offset, &displacement, symbol ) ) {
            printf( "[%i] %s\n", i, symbol->Name );
        } else {
            printf( "[%i] ???\n", i );
        }
    }
    fflush( stdout );
    SymCleanup( process );
#endif

#endif
}
