/**
 * @file ParallelDOFNumberingInterface.cxx
 * @brief Interface python de ParallelDOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

#include "PythonBindings/ParallelDOFNumberingInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelDOFNumberingToPython( py::module_ &mod ) {

    py::class_< ParallelDOFNumbering, ParallelDOFNumbering::ParallelDOFNumberingPtr,
                BaseDOFNumbering >( mod, "ParallelDOFNumbering" )
        .def( py::init( &initFactoryPtr< ParallelDOFNumbering > ) )
        .def( py::init( &initFactoryPtr< ParallelDOFNumbering, std::string > ) )
        // ---------------------------------------------------------------------
        .def( "useLagrangeMultipliers", &ParallelDOFNumbering::useLagrangeMultipliers, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "useSingleLagrangeMultipliers", &ParallelDOFNumbering::useSingleLagrangeMultipliers,
              R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        .def( "getNodeAssociatedToRow", &ParallelDOFNumbering::getNodeAssociatedToRow,
              R"(
Returns the node index associated to a dof index.

Arguments:
    row (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "row" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "isRowAssociatedToPhysical", &ParallelDOFNumbering::isRowAssociatedToPhysical,
              R"(
If the row is associated to a physical DOF, return True

If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, return False

Arguments:
    row (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "row" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getRowsAssociatedToPhysicalDofs",
              &ParallelDOFNumbering::getRowsAssociatedToPhysicalDofs,
              R"(
Returns the indexes of the physical dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: indexes of the physical dof.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getRowAssociatedToNodeComponent",
              &ParallelDOFNumbering::getRowAssociatedToNodeComponent,
              R"(
Returns the index of the dof associated to a node.

Arguments:
    node (int): Index of the node.
    component (str): name of the component
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "node" ), py::arg( "component" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getRowsAssociatedToLagrangeMultipliers",
              &ParallelDOFNumbering::getRowsAssociatedToLagrangeMultipliers,
              R"(
Returns the indexes of the Lagrange multipliers dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: indexes of the Lagrange multipliers dof.
        )",
              py::arg( "local" ) = false )
        .def( "getComponents", &ParallelDOFNumbering::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )" )
        // ---------------------------------------------------------------------
        .def( "getComponents", &ParallelDOFNumbering::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )" )
        // ---------------------------------------------------------------------
        .def( "getComponentAssociatedToRow", &ParallelDOFNumbering::getComponentAssociatedToRow,
              R"(
 Returns the components name associated to a dof index.

Arguments:
    node (int): Index of the node.
    local (bool): local or parallel request

Returns:
    str: component names.
              )",
              py::arg( "row" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getComponentsAssociatedToNode", &ParallelDOFNumbering::getComponentsAssociatedToNode,
              R"(
Returns the components name associated to a node index.

Arguments:
    node (int): Index of the node.
    local (bool): local or parallel request

Returns:
    str: component names.
        )",
              py::arg( "node" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getNumberOfDofs", &ParallelDOFNumbering::getNumberOfDofs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): local or parallel request

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false );
};

#endif /* ASTER_HAVE_MPI */
