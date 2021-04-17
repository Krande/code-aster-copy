/**
 * @file FortranInterface.cxx
 * @brief Python bindings for Fortran interface.
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

#include <boost/python.hpp>

namespace py = boost::python;

#include "PythonBindings/MatrixToPetscInterface.h"
#include "Solvers/MatrixToPetsc.h"

void exportMatrixToPetscToPython() {

    py::def( "petscFinalize", petscFinalize, R"(
Stops the PETSc interface.
        )" );
    //
    py::def( "_petscInitializeWithOptions", petscInitializeWithOptions, R"(
Starts the PETSc interface with options.

Arguments:
    options[str]: PETSc options

        )",
             ( py::args( "options" ) ) );

    py::def( "assemblyMatrixToPetsc", &assemblyMatrixToPetsc< AssemblyMatrixDisplacementRealPtr >,
             R"(
Convert a *AssemblyMatrix* object to a PETSc *Mat* object.

Arguments:
    matr (*AssemblyMatrix*): code_aster matrix.

Returns:
    *Mat*: PETSc matrix.
        )",
             ( py::arg( "matr" ) ) );
};
