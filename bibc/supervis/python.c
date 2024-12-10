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

/* person_in_charge: mathieu.courtois@edf.fr */
/* Minimal main program -- everything is loaded from the library */

/*! \mainpage Code_Aster Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

#include "aster.h"

#include "aster_core_module.h"
#include "aster_fonctions_module.h"
#include "aster_module.h"
#include "med_aster_module.h"

#include <stdio.h>

#ifndef _MAIN_
#define _MAIN_ main
#endif

void initAsterModules() {
    // Python 3.12 raises an error if PyImport_AppendInittab is called after Py_Initialize
    if ( !Py_IsInitialized() ) {
        PyImport_AppendInittab( "aster_core", PyInit_aster_core );
        PyImport_AppendInittab( "aster", PyInit_aster );

        /* Module définissant des opérations sur les objets fonction_sdaster */
        PyImport_AppendInittab( "aster_fonctions", PyInit_aster_fonctions );
#ifdef ASTER_HAVE_MED
        PyImport_AppendInittab( "med_aster", PyInit_med_aster );
#endif
    }
}

int _MAIN_( int argc, char **argv ) {
    int ierr;

    initAsterModules();

    wchar_t **wargv = PyMem_Malloc( sizeof( wchar_t * ) * argc );
    int i;
    for ( i = 0; i < argc; i++ )
        wargv[i] = Py_DecodeLocale( argv[i], NULL );

    ierr = Py_Main( argc, wargv );
    return ierr;
}
