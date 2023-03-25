/**
 * @file Fortran.h
 * @brief Definition of interface functions between C++ and Fortran
 * @author Mathieu Courtois
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: mathieu.courtois@edf.fr */

#include "astercxx.h"

#include "Solvers/MatrixToPetsc.h"

#include "aster_fort_petsc.h"
#include "aster_module.h"

#ifdef ASTER_HAVE_PETSC
void petscFinalize() {
    std::string off = "OFF", foo = " ";
    CALLO_AP_ON_OFF( off, foo );
    std::cout << "...PETSc finalized" << std::endl;
}

void petscInitializeWithOptions( const std::string &options ) {

    std::string on = "ON";
    CALLO_AP_ON_OFF( on, options );
    std::cout << "PETSc initialized..." << std::endl;
}
#else
void petscFinalize() { std::cout << "PETSc library non available" << std::endl; };

void petscInitializeWithOptions( const std::string &options ) {
    std::cout << "PETSc library non available" << std::endl;
}
#endif

#ifdef ASTER_HAVE_PETSC4PY
void DEFPPP( CREATE_CUSTOM_KSP, create_custom_ksp, _OUT KSP *ksp, _IN Mat *mat,
             _OUT ASTERINTEGER *ierr ) {

    /* get user function to build KSP solver */
    std::string mcf( "SOLVEUR" );
    std::string mcs( "KSP_UTIL" );
    PyObject *tuple;
    getvpy( mcf.c_str(), mcs.c_str(), 1, &tuple );
    AS_ASSERT( PyTuple_Check( tuple ) && PyTuple_Size( tuple ) == 1 );
    PyObject *builder = PyTuple_GetItem( tuple, 0 );
    Py_DECREF( tuple );

    /* create petsc4py Mat */
    PyObject *petsc4py = PyImport_ImportModule( "petsc4py.PETSc" );
    PyObject *pyMat = PyObject_CallMethod( petsc4py, "Mat", NULL );
    /* filled with Mat */
    struct PyPetscMatObject *pyx_mat = (struct PyPetscMatObject *)pyMat;
    Mat matInit = pyx_mat->mat;
    pyx_mat->mat = *mat;
    /* call user function */
    PyObject *pyKSP = PyObject_CallFunction( builder, "O", pyMat );

    /* extract KSP from petsc4py KSP */
    struct PyPetscKSPObject *pyx_ksp = (struct PyPetscKSPObject *)pyKSP;
    *ksp = pyx_ksp->ksp;

    // restore pyMat content to be safely destroyed by gc
    pyx_mat->mat = matInit;
    Py_DECREF( pyMat );
    Py_DECREF( petsc4py );

    return;
}
#endif
