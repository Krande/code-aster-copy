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
        )" )
        .def( "getNodesAndComponentsFromDOF",
              &GlobalEquationNumbering::getNodesAndComponentsFromDOF, R"(
            Return the list of node id and name of component for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                list[tuple[int, str]] : node id and name of component for each dofs
            )",
              py::arg( "local" ) = true )
        .def( "getNodesAndComponentsNumberFromDOF",
              &GlobalEquationNumbering::getNodesAndComponentsNumberFromDOF, R"(
            Return the list of node id and component id for each dofs

            Arguments:
                local (bool) = True: if True use local node index else use global index in HPC

            Returns:
                list[tuple[int, int]] : node id and component if for each dofs
            )",
              py::arg( "local" ) = true )
        .def( "getDOFsFromNodesAndComponentsNumber",
              &GlobalEquationNumbering::getDOFsFromNodesAndComponentsNumber, R"(
            Return the dict of dofs with the pair (node id, name id) as keys

            Arguments:
                local (bool) = True: if True use local DOF index else use global index in HPC

            Returns:
                dict[int, str] : dofs id for each node id and component id
            )",
              py::arg( "local" ) = true )
        .def( "getDOFsFromNodesAndComponents",
              &GlobalEquationNumbering::getDOFsFromNodesAndComponents, R"(
           Return the dict of dofs with the pair (node id, component's name) as keys

            Arguments:
                local (bool) = True: if True use local dof index else use global index in HPC

            Returns:
                dict[int, str] : dofs id for each node id and component's name
            )",
              py::arg( "local" ) = true )
        .def( "getComponents", &GlobalEquationNumbering::getComponents, R"(
            Get list of components

            Returns:
                list[str]: list of components
            )" );
};
