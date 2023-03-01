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

#include "PythonBindings/GlobalEquationNumberingInterface.h"

#include "aster_pybind.h"

void exportGlobalEquationNumberingToPython( py::module_ &mod ) {

    py::class_< GlobalEquationNumbering, GlobalEquationNumberingPtr, DataStructure >(
        mod, "GlobalEquationNumbering" )
        .def( py::init( &initFactoryPtr< GlobalEquationNumbering > ) )
        .def( py::init( &initFactoryPtr< GlobalEquationNumbering, std::string > ) )
        .def( "getModel", &GlobalEquationNumbering::getModel )
        .def( "setModel", &GlobalEquationNumbering::setModel )
        .def( "getMesh", &GlobalEquationNumbering::getMesh )
        .def( "setMesh", &GlobalEquationNumbering::setMesh )
        // ---------------------------------------------------------------------
        .def( "useLagrangeMultipliers", &GlobalEquationNumbering::useLagrangeMultipliers, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "useSingleLagrangeMultipliers",
              &GlobalEquationNumbering::useSingleLagrangeMultipliers,
              R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNumberOfDofs", &GlobalEquationNumbering::getNumberOfDofs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): not used.

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false )
        .def( "getPhysicalQuantity", &GlobalEquationNumbering::getPhysicalQuantity, R"(
Returns the name of the physical quantity that is numbered.

Returns:
    str: physical quantity name.
        )" )
        .def( "isParallel", &GlobalEquationNumbering::isParallel, R"(
The numbering is distributed across MPI processes for High Performance Computing.

Returns:
    bool: *True* if used, *False* otherwise.
        )" );
};
