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

#include "aster.h"

#include "aster_fort_utils.h"
#include "definition_pt.h"
#include "dll_register.h"

#ifdef ASTER_PLATFORM_MSYS2
#include <windows.h>
#define dlopen( libname, flag ) LoadLibrary( libname )
#define dlsym( handle, symbol ) GetProcAddress( (HMODULE)( handle ), (LPCSTR)( symbol ) )
static inline int my_dlclose( void *handle ) { return FreeLibrary( (HMODULE)handle ) ? 0 : 1; }
#define dlclose my_dlclose
#define RTLD_LAZY 0

#elseif ASTER_PLATFORM_MSVC64
#include <windows.h>
// Windows-specific implementation of a function to unload libraries
static void windows_dlclose(void *handle) {
    FreeLibrary((HMODULE)handle);
}
#define dlclose windows_dlclose
#else
#include <dlfcn.h>
#endif

/* *********************************************************************
 *
 * Utilities to Load Dynamically (optionnal) external Libraries
 *
 * Supported components : UMAT
 *
 * *********************************************************************/

/* Global dictionnary used to register (libraries, symbol) couples */
static PyObject *DLL_DICT = NULL;

void dll_init() {
    /* Initialization */
    if ( !DLL_DICT ) {
        DLL_DICT = PyDict_New();
    }
}

PyObject *get_dll_register_dict() {
    /* Return the register dictionnary.
     * For external modules. */
    dll_init();
    return DLL_DICT;
}

void DEF0( DLLCLS, dllcls ) {
    /* Unload all components
     */
    dll_init();
    libsymb_apply_on_all( DLL_DICT, (FUNC_PTR)dlclose, 1 );
    Py_DECREF( DLL_DICT );
    DLL_DICT = NULL;
}
