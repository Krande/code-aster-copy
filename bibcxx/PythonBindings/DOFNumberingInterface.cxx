/**
 * @file DOFNumberingInterface.cxx
 * @brief Interface python de DOFNumbering
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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
            &initFactoryPtr< DOFNumbering, std::string, EquationNumberingPtr, ModelPtr > ) )
        // ----------------------------------------------------------------------
        .def( "useLagrangeDOF", &DOFNumbering::useLagrangeDOF, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ----------------------------------------------------------------------
        .def( "useSingleLagrangeDOF", &DOFNumbering::useSingleLagrangeDOF, R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ----------------------------------------------------------------------
        .def( "getNodeFromDOF", &DOFNumbering::getNodeFromDOF,
              R"(
Returns the node index associated to a dof index.

Arguments:
    dof (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the node.
        )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getNoGhostDOFs", &DOFNumbering::getNoGhostDOFs,
              R"(
Returns the indexes of the DOFs owned locally (aka not ghost).

Returns:
  int: indexes of the DOFs owned locally.
      )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getDOFFromNodeAndComponent", &DOFNumbering::getDOFFromNodeAndComponent,
              R"(
Returns the DOF index associated to a node and component.

Arguments:
    node (int): Index of the node.
    cmp (str): name of the component
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "node" ), py::arg( "cmp" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "isPhysicalDOF", &DOFNumbering::isPhysicalDOF,
              R"(
If the dof is associated to a physical DOF, return True

If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, return False

Arguments:
    dof (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    bool: True if the DOF is a physical DOF else False.
        )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getPhysicalDOFs", &DOFNumbering::getPhysicalDOFs,
              R"(
Returns the indexes of the physical dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    [int]: indexes of the physical dof.
        )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getLagrangeDOFs", &DOFNumbering::getLagrangeDOFs,
              R"(
Returns the indexes of the Lagrange multipliers dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    [int]: indexes of the Lagrange multipliers dof.
        )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getDictOfLagrangeDOFs", &DOFNumbering::getDictOfLagrangeDOFs,
              py::return_value_policy::copy,
              R"(
Returns the Rows Associated to the first and second Lagrange Multipliers Dof

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    [dict]: {1 : indexes of the first Lagrange multipliers dof,
             2 : indexes of the second Lagrange multipliers dof }
        )",
              py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getComponents", &DOFNumbering::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )" )
        // ---------------------------------------------------------------------
        .def( "getComponentFromDOF", &DOFNumbering::getComponentFromDOF,
              R"(
Returns the component name associated to a dof index.

- If the dof is associated to a physical DOF, the name of the component is returned.

- If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, the name of the component which is constrained by the multiplier is
  returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.

- If the dof is associated to a Lagrange multiplier DOF for a multipoint-constraint
  (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
  identified).

Arguments:
    dof (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    str: component name.
        )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getComponentFromNode", &DOFNumbering::getComponentFromNode,
              R"(
Returns the components name associated to a node index.

Arguments:
    node (int): Index of the node.
    local (bool, optional): not used (default: false).

Returns:
    str: component names.
        )",
              py::arg( "node" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const bool >( &DOFNumbering::getNodeAndComponentFromDOF,
                                               py::const_ ),
              R"(
Return the list of node id and name of component for each dofs

Arguments:
    local (bool, optional): not used (default: false).
Returns:
    list[tuple[int, str]] : node id and name of component for each dofs
            )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const ASTERINTEGER, const bool >(
                  &DOFNumbering::getNodeAndComponentFromDOF, py::const_ ),
              R"(
Return the node id and name of component for given DOF

Arguments:
    dof (int): DOF index
    local (bool, optional): not used (default: false).
Returns:
    tuple[int, str] : node id and name of component
            )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ----------------------------------------------------------------------
        .def( "getNumberOfDOFs", &DOFNumbering::getNumberOfDOFs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false );
};
