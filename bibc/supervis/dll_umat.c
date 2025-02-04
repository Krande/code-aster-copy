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

#include "aster.h"

#include "aster_fort_utils.h"
#include "aster_utils.h"
#include "definition_pt.h"
#include "dll_register.h"

#ifdef ASTER_PLATFORM_POSIX
#include <dlfcn.h>
PyObject *get_dll_register_dict();

/* *********************************************************************
 *
 *                          UMAT interface
 *
 * *********************************************************************/

/* declarations of pointers on UMAT functions */
#define FUNC_UMAT( NAME )                                                                          \
    void DEFUMAT( *NAME, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *,               \
                  ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *,       \
                  ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *,       \
                  ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, char *, STRING_SIZE, \
                  ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, ASTERDOUBLE *,   \
                  ASTERINTEGER *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *, ASTERDOUBLE *,      \
                  ASTERDOUBLE *dfgrd0, ASTERDOUBLE *dfgrd1, ASTERINTEGER *, ASTERINTEGER *,        \
                  ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER * )

void load_umat_lib( const char *libname, const char *symbol ) {
    /* load UMAT library and initialize pointers to UMAT functions
     */
    void *umat_handle;
    char *error;
    char symbol_[18], *valk;
    ASTERINTEGER ibid = 0, n0 = 0, nk = 0;
    ASTERDOUBLE rbid = 0.;
    FUNC_UMAT( f_umat ) = NULL;
    PyObject *DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    strcpy( symbol_, symbol );

    printf( "Loading '%s'... ", libname );
    umat_handle = dlopen( libname, RTLD_NOW );
    if ( !umat_handle ) {
        printf( "\n%s\n", dlerror() );
        nk = 2;
        valk = MakeTabFStr( nk, VALK_SIZE );
        SetTabFStr( valk, 0, "UMAT", VALK_SIZE );
        SetTabFStr( valk, 1, (char *)libname, VALK_SIZE );
        CALL_UTMESS_CORE( "F", "FERMETUR_13", &nk, valk, &n0, &ibid, &n0, &rbid, &n0, " " );
        FreeStr( valk ); // uncallable
    }
    printf( "searching symbol '%s'... ", symbol );
    dlerror(); /* Clear any existing error */

    *(void **)( &f_umat ) = dlsym( umat_handle, symbol );
    if ( ( error = dlerror() ) != NULL ) {
        dlerror();
        strcat( symbol_, "_" );
        printf( "trying symbol '%s'... ", symbol_ );
        *(void **)( &f_umat ) = dlsym( umat_handle, symbol_ );
    }

    if ( ( error = dlerror() ) != NULL ) {
        printf( "\n%s\n", error );
        nk = 3;
        valk = MakeTabFStr( nk, VALK_SIZE );
        SetTabFStr( valk, 0, "UMAT", VALK_SIZE );
        SetTabFStr( valk, 1, (char *)libname, VALK_SIZE );
        SetTabFStr( valk, 2, (char *)symbol, VALK_SIZE );
        CALL_UTMESS_CORE( "F", "FERMETUR_14", &nk, valk, &n0, &ibid, &n0, &rbid, &n0, " " );
        FreeStr( valk ); // uncallable
    }
    printf( "found\n" );

    /* register these UMAT lib */
    if ( libsymb_register( DLL_DICT, libname, symbol, umat_handle, (FUNC_PTR)f_umat ) ) {
        printf( "Registering of '%s' and '%s' failed!\n", libname, symbol );
    }
}
#endif

void DEFSSP( UMAT_GET_FUNCTION, umat_get_function, char *nomlib, STRING_SIZE lnomlib, char *nomsub,
             STRING_SIZE lnomsub, ASTERINTEGER *pfumat ) {
#ifdef ASTER_PLATFORM_POSIX
    /* UMAT WraPper : wrapper to get the UMAT function.
     */
    char *libname, *symbol;
    FUNC_UMAT( f_umat ) = NULL;
    PyObject *DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr( nomlib, lnomlib );
    symbol = MakeCStrFromFStr( nomsub, lnomsub );

    DEBUG_DLL_VV( " libname = >%s<, len = %d\n", libname, (int)strlen( libname ) )
    DEBUG_DLL_VV( "  symbol = >%s<, len = %d\n", symbol, (int)strlen( symbol ) )

    if ( !libsymb_is_known( DLL_DICT, libname, symbol ) ) {
        load_umat_lib( libname, symbol );
    }
    *pfumat = (ASTERINTEGER)libsymb_get_symbol( DLL_DICT, libname, symbol );

    FreeStr( libname );
    FreeStr( symbol );
#else
    printf( "UMAT: Not available under Windows.\n" );
    fflush( stdout );
    abort();
#endif
}

void DEFPPPPPPPPPPPPPPPPPPPSPPPPPPPPPPPPPPPPPP(
    UMATWP, umatwp, ASTERINTEGER *pfumat, ASTERDOUBLE *stress, ASTERDOUBLE *statev,
    ASTERDOUBLE *ddsdde, ASTERDOUBLE *sse, ASTERDOUBLE *spd, ASTERDOUBLE *scd, ASTERDOUBLE *rpl,
    ASTERDOUBLE *ddsddt, ASTERDOUBLE *drplde, ASTERDOUBLE *drpldt, ASTERDOUBLE *stran,
    ASTERDOUBLE *dstran, ASTERDOUBLE *time, ASTERDOUBLE *dtime, ASTERDOUBLE *temp,
    ASTERDOUBLE *dtemp, ASTERDOUBLE *predef, ASTERDOUBLE *dpred, char *cmname, STRING_SIZE lcmname,
    ASTERINTEGER *ndi, ASTERINTEGER *nshr, ASTERINTEGER *ntens, ASTERINTEGER *nstatv,
    ASTERDOUBLE *props, ASTERINTEGER *nprops, ASTERDOUBLE *coords, ASTERDOUBLE *drot,
    ASTERDOUBLE *pnewdt, ASTERDOUBLE *celent, ASTERDOUBLE *dfgrd0, ASTERDOUBLE *dfgrd1,
    ASTERINTEGER *noel, ASTERINTEGER *npt, ASTERINTEGER *layer, ASTERINTEGER *kspt,
    ASTERINTEGER *kstep, ASTERINTEGER *kinc ) {
#ifdef ASTER_PLATFORM_POSIX
    /* UMAT WraPper : wrapper to the UMAT function through the function pointer
     * Load the library if necessary (at the first call).
     */
    FUNC_UMAT( f_umat ) = NULL;

    f_umat = ( FUNC_UMAT() )( *pfumat );

    CALLUMAT( *f_umat, stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran,
              dstran, time, dtime, temp, dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv,
              props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt,
              kstep, kinc );
#else
    printf( "UMAT: Not available under Windows.\n" );
    fflush( stdout );
    abort();
#endif
}
