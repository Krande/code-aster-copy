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

/* person_in_charge: mathieu.courtois at edf.fr */

#include "aster.h"

#include "aster_utils.h"

#ifdef ASTER_HAVE_HDF5
#include "hdf5.h"
#endif
#ifdef ASTER_HAVE_MED
#include "med.h"
#endif
#ifdef ASTER_HAVE_SCOTCH
/* scotch.h may use int64_t without including <sys/types.h> */
#include "scotch.h"

#include <sys/types.h>

#ifdef ASTER_PLATFORM_MINGW
#include <stdint.h>
#endif
#endif

void DEFPPP( LIHDFV, lihdfv, _OUT ASTERINTEGER *major, _OUT ASTERINTEGER *minor,
             _OUT ASTERINTEGER *patch ) {
    /* Retourne la version de HDF5 */
    unsigned int n1 = 0, n2 = 0, n3 = 0;
#ifdef ASTER_HAVE_HDF5
    herr_t ier = 0;
    ier = (int)H5get_libversion( &n1, &n2, &n3 );
#endif
    *major = (ASTERINTEGER)n1;
    *minor = (ASTERINTEGER)n2;
    *patch = (ASTERINTEGER)n3;
    return;
}

void DEFPPP( LIMEDV, limedv, _OUT ASTERINTEGER *major, _OUT ASTERINTEGER *minor,
             _OUT ASTERINTEGER *patch ) {
    /* Retourne la version de MED */
    int ier = 0;
#ifdef ASTER_HAVE_MED
    med_int n1 = 0, n2 = 0, n3 = 0;
    ier = (int)MEDlibraryNumVersion( &n1, &n2, &n3 );
#else
    unsigned int n1 = 0, n2 = 0, n3 = 0;
#endif
    *major = (ASTERINTEGER)n1;
    *minor = (ASTERINTEGER)n2;
    *patch = (ASTERINTEGER)n3;
    return;
}

void DEFPPP( LISCOV, liscov, _OUT ASTERINTEGER *major, _OUT ASTERINTEGER *minor,
             _OUT ASTERINTEGER *patch ) {
    /* Retourne la version de SCOTCH */
#ifdef ASTER_HAVE_SCOTCH
    *major = (ASTERINTEGER)SCOTCH_VERSION;
    *minor = (ASTERINTEGER)SCOTCH_RELEASE;
    *patch = (ASTERINTEGER)SCOTCH_PATCHLEVEL;
#else
    *major = (ASTERINTEGER)0;
    *minor = (ASTERINTEGER)0;
    *patch = (ASTERINTEGER)0;
#endif
    return;
}
