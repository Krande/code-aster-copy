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

#include "PythonBindings/DOFNumberingInterface.h"

#include "aster_pybind.h"

void exportDOFNumberingToPython( py::module_ &mod ) {

    py::class_< DOFNumbering, DOFNumbering::DOFNumberingPtr, BaseDOFNumbering >( mod,
                                                                                 "DOFNumbering" )
        .def( py::init( &initFactoryPtr< DOFNumbering > ) )
        .def( py::init( &initFactoryPtr< DOFNumbering, std::string > ) )
        .def( py::init(
            &initFactoryPtr< DOFNumbering, std::string, FieldOnNodesDescriptionPtr, MeshPtr > ) )
        // ----------------------------------------------------------------------
        .def( "useLagrangeMultipliers", &DOFNumbering::useLagrangeMultipliers, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ----------------------------------------------------------------------
        .def( "useSingleLagrangeMultipliers", &DOFNumbering::useSingleLagrangeMultipliers, R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ----------------------------------------------------------------------
        .def( "getNodeAssociatedToRow", &DOFNumbering::getNodeAssociatedToRow,
              R"(
Returns the node index associated to a dof index.

Arguments:
    row (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              // ----------------------------------------------------------------------
              py::arg( "row" ), py::arg( "local" ) = false )
        .def( "isRowAssociatedToPhysical", &DOFNumbering::isRowAssociatedToPhysical,
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
        // ----------------------------------------------------------------------
        .def( "getRowsAssociatedToPhysicalDofs", &DOFNumbering::getRowsAssociatedToPhysicalDofs,
              R"(
Returns the indexes of the physical dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    [int]: indexes of the physical dof.
        )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getRowAssociatedToNodeComponent", &DOFNumbering::getRowAssociatedToNodeComponent,
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
        // ----------------------------------------------------------------------
        .def( "getRowsAssociatedToLagrangeMultipliers",
              &DOFNumbering::getRowsAssociatedToLagrangeMultipliers,
              R"(
Returns the indexes of the Lagrange multipliers dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    [int]: indexes of the Lagrange multipliers dof.
        )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getComponents", &DOFNumbering::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )" )
        // ---------------------------------------------------------------------
        .def( "getComponentAssociatedToRow", &DOFNumbering::getComponentAssociatedToRow,
              R"(
Returns the component name associated to a dof index.

- If the row is associated to a physical DOF, the name of the component is returned.

- If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the component which is constrained by the multiplier is
  returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.

- If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
  (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
  identified).

Arguments:
    row (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    str: component name.
        )",
              py::arg( "row" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getComponentsAssociatedToNode", &DOFNumbering::getComponentsAssociatedToNode,
              R"(
Returns the components name associated to a node index.

Arguments:
    node (int): Index of the node.
    local (bool, optional): not used (default: false).

Returns:
    str: component names.
        )",
              py::arg( "node" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getNumberOfDofs", &DOFNumbering::getNumberOfDofs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false );
};
