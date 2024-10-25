/**
 * @file FortranInterface.cxx
 * @brief Python bindings for Fortran interface.
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

#include "PythonBindings/MatrixToPetscInterface.h"

#include "aster_pybind.h"

#include "Solvers/MatrixToPetsc.h"

void exportMatrixToPetscToPython( py::module_ &mod ) {

    mod.def( "petscFinalize", petscFinalize, R"(
Stops the PETSc interface.
        )" );
    //
    mod.def( "petscInitialize", petscInitializeWithOptions, R"(
Starts the PETSc interface with options.

Arguments:
    options[str]: PETSc options

        )",
             py::arg( "options" ) = "" );

    mod.def( "assemblyMatrixToPetsc", &assemblyMatrixToPetsc< AssemblyMatrixDisplacementRealPtr >,
             R"(
Convert a *AssemblyMatrix* object to a PETSc *Mat* object.

Arguments:
    matr (*AssemblyMatrix*): code_aster matrix.
    local (*bool*): extract only the sequential matrix of the subdomain or the global parallel
                    matrix

Returns:
    *Mat*: PETSc matrix.
        )",
             py::arg( "matr" ), py::arg( "local" ) );

    mod.def( "assemblyMatrixToPetsc", &assemblyMatrixToPetsc< AssemblyMatrixTemperatureRealPtr >,
             R"(
Convert a *AssemblyMatrix* object to a PETSc *Mat* object.

Arguments:
    matr (*AssemblyMatrix*): code_aster matrix.
    local (*bool*): extract only the sequential matrix of the subdomain or the global parallel
                    matrix

Returns:
    *Mat*: PETSc matrix.
        )",
             py::arg( "matr" ), py::arg( "local" ) );
};
