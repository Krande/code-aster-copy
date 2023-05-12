/**
 * @file ParallelDOFNumberingInterface.cxx
 * @brief Interface python de ParallelDOFNumbering
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

#include "PythonBindings/ParallelDOFNumberingInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelDOFNumberingToPython( py::module_ &mod ) {

    py::class_< ParallelDOFNumbering, ParallelDOFNumbering::ParallelDOFNumberingPtr,
                BaseDOFNumbering >( mod, "ParallelDOFNumbering" )
        .def( py::init( &initFactoryPtr< ParallelDOFNumbering > ) )
        .def( py::init( &initFactoryPtr< ParallelDOFNumbering, std::string > ) )
        .def( py::init( &initFactoryPtr< ParallelDOFNumbering, std::string,
                                         ParallelEquationNumberingPtr, ModelPtr > ) )
        // ---------------------------------------------------------------------
        .def( "useLagrangeDOF", &ParallelDOFNumbering::useLagrangeDOF, R"(
Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "useSingleLagrangeDOF", &ParallelDOFNumbering::useSingleLagrangeDOF,
              R"(
Single Lagrange multipliers are used for BC or MPC.

Returns:
    bool: *True* if used, *False* otherwise.
        )" )
        // ---------------------------------------------------------------------
        .def( "getNodeFromDOF", &ParallelDOFNumbering::getNodeFromDOF,
              R"(
Returns the node index associated to a dof index.

Arguments:
    dof (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "isPhysicalDOF", &ParallelDOFNumbering::isPhysicalDOF,
              R"(
If the dof is associated to a physical DOF, return True

If the dof is associated to a Lagrange multiplier DOF for a Dirichlet boundary
  condition, return False

Arguments:
    dof (int): Index of the dof.
    local (bool, optional): not used (default: false).

Returns:
    int: index of the dof.
        )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getPhysicalDOF", &ParallelDOFNumbering::getPhysicalDOF,
              R"(
Returns the indexes of the physical dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: indexes of the physical dof.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getLagrangeDOF", &ParallelDOFNumbering::getLagrangeDOF,
              R"(
Returns the indexes of the Lagrange multipliers dof.

Arguments:
    local (bool, optional): not used (default: false).

Returns:
    int: indexes of the Lagrange multipliers dof.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getComponents", &ParallelDOFNumbering::getComponents, R"(
Returns all the component names assigned in the numbering.

Returns:
    str: component names.
        )" )
        // ---------------------------------------------------------------------
        .def( "getComponentFromDOF", &ParallelDOFNumbering::getComponentFromDOF,
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
    node (int): Index of the node.
    local (bool): dof in local or global numbering

Returns:
    str: component names.
              )",
              py::arg( "dof" ), py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getComponentFromNode", &ParallelDOFNumbering::getComponentFromNode,
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
        .def( "getNumberOfDOF", &ParallelDOFNumbering::getNumberOfDOF,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): local or parallel request

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getGhostDOF", &ParallelDOFNumbering::getGhostDOF,
              R"(
Returns the indexes of the ghost DOFs.

Arguments:
    local (bool): local or global numbering

Returns:
    int: indexes of the ghost DOFs.
        )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNoGhostDOF", &ParallelDOFNumbering::getNoGhostDOF,
              R"(
Returns the indexes of the DOFs owned locally (aka not ghost).

Returns:
    int: indexes of the DOFs owned locally.
        )" )
        // ---------------------------------------------------------------------
        .def( "localToGlobalDOF", &ParallelDOFNumbering::localToGlobalDOF,
              R"(
Returns the global number of a local DOF.

Arguments:
    loc (int): local DOF number

Returns:
    int: global number of the DOF.
        )",
              py::arg( "loc" ) )
        // ---------------------------------------------------------------------
        .def( "getLocalToGlobalMapping", &ParallelDOFNumbering::getLocalToGlobalMapping,
              R"(
Returns the mapping from the local to the global number of the DOFs.

Returns:
    int: global number of the DOF.
        )" )
        // ---------------------------------------------------------------------
        .def( "globalToLocalDOF", &ParallelDOFNumbering::globalToLocalDOF,
              R"(
Returns the local number of a global DOF.

Arguments:
    glob (int): global DOF number

Returns:
    int: local number of the DOF.
        )",
              py::arg( "glob" ) );
};

#endif /* ASTER_HAVE_MPI */
