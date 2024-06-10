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

/* ------------------------------------------------------------------ */
/*
** Initialisation de l'interception de certains signaux
** Actuellement sont traites les signaux :
**    CPULIM  : plus de temps CPU
**    FPE     : Floating point exception
*/

#include "aster.h"

#include "aster_fort_utils.h"

#include <signal.h>

void abort();

void hancpu( int sig );

#if defined ASTER_PLATFORM_SOLARIS
#include <siginfo.h>
#include <ucontext.h>
void hanfpe( int sig, siginfo_t *sip, ucontext_t *uap );

#elif defined ASTER_PLATFORM_MINGW
#include <float.h>
void hanfpe( int sig );
void sigsegv( int sig );

#elif defined ASTER_PLATFORM_POSIX
void hanfpe( int sig );
void stpusr1( int sig );

#elif defined ASTER_PLATFORM_MSVC64
#include <windows.h>
#include <float.h>
#include <excpt.h>
void hanfpe(int sig);
#endif

#if defined ASTER_PLATFORM_LINUX
#define _GNU_SOURCE 1
#include <fenv.h>
#endif

void DEF0( INISIG, inisig ) {
#if defined ASTER_PLATFORM_POSIX
    struct sigaction action_CPU_LIM;
#else
    unsigned int cw, cwOrig;
#endif

/*            */
/* CPU LIMITE */
/*            */
#if defined ASTER_PLATFORM_POSIX
    action_CPU_LIM.sa_handler = hancpu;
    sigemptyset( &action_CPU_LIM.sa_mask );
    action_CPU_LIM.sa_flags = 0;
    sigaction( SIGXCPU, &action_CPU_LIM, NULL );
#endif

/*                          */
/* Floating point exception */
/*                          */
#if defined ASTER_PLATFORM_SOLARIS
    ieee_handler( "set", "common", hanfpe );
    ieee_handler( "clear", "invalid", hanfpe );

#elif defined ASTER_PLATFORM_LINUX

    /* Enable some exceptions. At startup all exceptions are masked. */
    feenableexcept( FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID );

    signal( SIGFPE, hanfpe );

#elif defined ASTER_PLATFORM_MINGW
    _clearfp();
    cw = _controlfp( 0, 0 );
    cw &= ~( _EM_OVERFLOW | _EM_ZERODIVIDE );
    cwOrig = _controlfp( cw, _MCW_EM );

    signal( SIGFPE, hanfpe );
#elif defined ASTER_PLATFORM_MSVC64
    errno_t err;
    _clearfp(); // Clear any pending floating-point exceptions.
    // Get the current control word
    err = _controlfp_s(&cwOrig, 0, 0);
    if (err != 0) {
        fprintf(stderr, "Failed to get the control word\n");
        abort(); // Handle the error as appropriate
    }

    // Modify the control word to trap exceptions
    cw = cwOrig & ~(_EM_OVERFLOW | _EM_ZERODIVIDE | _EM_INVALID);

    // Set the new control word
    err = _controlfp_s(NULL, cw, _MCW_EM);
    if (err != 0) {
        fprintf(stderr, "Failed to set the control word\n");
        abort(); // Handle the error as appropriate
    }

    signal(SIGFPE, hanfpe);
#else
    signal( SIGFPE, hanfpe );
#endif

/*                          */
/* Arret par SIGUSR1        */
/*                          */
/* Note : l'arret par SIGUSR1 ne fonctionne pas sous MSVC,
   il faudra essayer de trouver autre chose... */
#if defined ASTER_PLATFORM_POSIX
    signal( SIGUSR1, stpusr1 );
#elif defined ASTER_PLATFORM_MINGW
    signal( SIGSEGV, sigsegv );
#endif
}

static ASTERINTEGER status_usr1 = 0;

ASTERINTEGER DEF0( ETAUSR, etausr ) {
    /* ETAt USR1 :
     * Retourne la variable status_usr1 */
    return status_usr1;
}

void stpusr1( int sig ) {
    /* SToP USR1 :
     * callback appelé lors de la réception du signal USR1.
     */
    CALL_UTMESS( "I", "SUPERVIS_96" );
    status_usr1 = (ASTERINTEGER)1;
}

#ifdef ASTER_PLATFORM_MINGW
void sigsegv( int sig ) {
    printf( "SIGSEGV\n" );
    print_trace_();
}
#endif

#ifdef ASTER_PLATFORM_MSVC64
void handle_exception(unsigned int code, struct _EXCEPTION_POINTERS* ep) {
    if (code == EXCEPTION_FLT_DIVIDE_BY_ZERO || code == EXCEPTION_FLT_OVERFLOW || code == EXCEPTION_FLT_INVALID_OPERATION) {
        hanfpe(SIGFPE);
    } else {
        // Handle other exceptions
    }
}
#endif

void DEF0( CLRUSR, clrusr ) {
    /* CLeaR USR1 :
     * Réinitialise la valeur de status_usr1
     * Utile pour éviter la récursivité.
     */
    status_usr1 = (ASTERINTEGER)0;
}
