/**
 * @file Fortran.h
 * @brief Definition of interface functions between C++ and Fortran
 * @author Mathieu Courtois
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include "Python.h"
#include "aster_fort_petsc.h"
#include "astercxx.h"

#include "Solvers/MatrixToPetsc.h"

#ifdef ASTER_HAVE_PETSC
void petscFinalize() {
    std::string off = "OFF", foo = " ";
    CALLO_AP_ON_OFF( off, foo );
    std::cout << "...PETSc finalized" << std::endl;
};

void petscInitializeWithOptions( const std::string &options ) {

    std::string on = "ON";
    CALLO_AP_ON_OFF( on, options );
    std::cout << "PETSc initialized..." << std::endl;
};
#else
void petscFinalize() { std::cout << "PETSc library non available" << std::endl; };

void petscInitializeWithOptions( const std::string &options ) {
    std::cout << "PETSc library non available" << std::endl;
};
#endif

template <>
PyObject *assemblyMatrixToPetsc< AssemblyMatrixDisplacementRealPtr >(
    const AssemblyMatrixDisplacementRealPtr matr );

