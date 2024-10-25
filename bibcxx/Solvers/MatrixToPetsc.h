#ifndef MATRIXTOPETSC_H_
#define MATRIXTOPETSC_H_

/**
 * @file MatrixToPetsc.h
 * @brief Definition of functions to convert matrix to PETSc
 * @author Mathieu Courtois
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "aster_pybind.h"

#include "LinearAlgebra/AssemblyMatrix.h"
#include "Loads/PhysicalQuantity.h"

/** @brief Turns PETSc off */
void petscFinalize();

/** @brief Turns PETSc on with options */
void petscInitializeWithOptions( const std::string &options = "" );

/** @brief Convert an AssemblyMatrix object to a PETSc Mat object */
template < class T >
py::object assemblyMatrixToPetsc( const T matr, const bool *local ) {
#ifdef ASTER_HAVE_PETSC4PY
    py::object petsc_matr = py::module_::import( "petsc4py.PETSc" ).attr( "Mat" )();

    Mat conv = matr->toPetsc( *local );

    struct PyPetscMatObject *pyx_mat = (struct PyPetscMatObject *)( petsc_matr.ptr() );
    pyx_mat->mat = conv;
    py::object result = py::reinterpret_steal< py::object >( (PyObject *)pyx_mat );
    petsc_matr.release();
    return result;
#else
    PyErr_SetString( PyExc_NotImplementedError, "petsc4py is not available" );
    throw py::error_already_set();
#endif
}

#ifdef ASTER_HAVE_PETSC4PY

#ifdef __cplusplus
extern "C" {
#endif

void DEFPPP( CREATE_CUSTOM_KSP, create_custom_ksp, _OUT KSP *ksp, _IN Mat *mat,
             _OUT ASTERINTEGER *ierr );

#ifdef __cplusplus
}
#endif

#endif

#endif /* MATRIXTOPETSC_H_ */
