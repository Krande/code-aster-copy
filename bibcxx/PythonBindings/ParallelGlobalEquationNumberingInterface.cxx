/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
 * @author Nicolas Sellenet
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
/* person_in_charge: nicolas.sellenet at edf.fr */

#include "PythonBindings/ParallelGlobalEquationNumberingInterface.h"

#include "aster_pybind.h"

void exportParallelGlobalEquationNumberingToPython( py::module_ &mod ) {

    py::class_< ParallelGlobalEquationNumbering, ParallelGlobalEquationNumberingPtr,
                GlobalEquationNumbering >( mod, "ParallelGlobalEquationNumbering" )
        .def( py::init( &initFactoryPtr< ParallelGlobalEquationNumbering, std::string > ) )
        .def( "getGhostRows", &ParallelGlobalEquationNumbering::getGhostRows,
              R"(
Returns the indexes of the ghost DOFs.

Arguments:
    local (bool): local or global numbering

Returns:
    int: indexes of the ghost DOFs.
        )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNoGhostRows", &ParallelGlobalEquationNumbering::getNoGhostRows,
              R"(
Returns the indexes of the DOFs owned locally (aka not ghost).

Returns:
    int: indexes of the DOFs owned locally.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNumberOfDofs", &ParallelGlobalEquationNumbering::getNumberOfDofs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): local or parallel request

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getLocalToGlobalMapping", &ParallelGlobalEquationNumbering::getLocalToGlobalMapping,
              R"(
Returns the mapping from the local to the global number of the DOFs.

Returns:
    int: global number of the DOF.
        )" )
        // ---------------------------------------------------------------------
        .def( "globalToLocalRow", &ParallelGlobalEquationNumbering::globalToLocalRow,
              R"(
Returns the local number of a global DOF.

Arguments:
    glob (int): global DOF number

Returns:
    int: local number of the DOF.
        )",
              py::arg( "glob" ) );
};
