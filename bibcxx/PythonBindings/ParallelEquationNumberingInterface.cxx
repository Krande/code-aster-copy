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

#include "PythonBindings/ParallelEquationNumberingInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelEquationNumberingToPython( py::module_ &mod ) {

    py::class_< ParallelEquationNumbering, ParallelEquationNumberingPtr, EquationNumbering >(
        mod, "ParallelEquationNumbering" )
        .def( py::init( &initFactoryPtr< ParallelEquationNumbering > ) )
        .def( py::init( &initFactoryPtr< ParallelEquationNumbering, std::string > ) )
        .def( "getGhostDOFs", &ParallelEquationNumbering::getGhostDOFs,
              R"(
Returns the indexes of the ghost DOFs.

Arguments:
    local (bool): local or global numbering
    lastLayerOnly (bool): last ghosts layer or all

Returns:
    int: indexes of the ghost DOFs.
        )",
              py::arg( "local" ) = true, py::arg( "lastLayerOnly" ) = false )
        // ---------------------------------------------------------------------
        .def( "getNoGhostDOFs", &ParallelEquationNumbering::getNoGhostDOFs,
              R"(
Returns the indexes of the DOFs owned locally (aka not ghost).

Returns:
    int: indexes of the DOFs owned locally.
        )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNumberOfDOFs", &ParallelEquationNumbering::getNumberOfDOFs,
              R"(
Returns the number of DOFs.

Arguments:
    local (bool): local or parallel request

Returns:
    int: number of DOFs.
        )",
              py::arg( "local" ) = false )
        // ---------------------------------------------------------------------
        .def( "getLocalToGlobalMapping", &ParallelEquationNumbering::getLocalToGlobalMapping,
              R"(
Returns the mapping from the local to the global number of the DOFs.

Returns:
    int: global number of the DOF.
        )" )
        // ---------------------------------------------------------------------
        .def( "globalToLocalDOF", &ParallelEquationNumbering::globalToLocalDOF,
              R"(
Returns the local number of a global DOF.

Arguments:
    glob (int): global DOF number

Returns:
    int: local number of the DOF.
        )",
              py::arg( "glob" ) )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const bool >(
                  &ParallelEquationNumbering::getNodeAndComponentFromDOF, py::const_ ),
              R"(
Return the list of node id and name of component for each dofs

Arguments:
  local (bool) = True: if True use local node index else use global index in HPC
Returns:
  list[tuple[int, str]] : node id and name of component for each dofs
          )",
              py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        .def( "getNodeAndComponentFromDOF",
              py::overload_cast< const ASTERINTEGER, const bool >(
                  &ParallelEquationNumbering::getNodeAndComponentFromDOF, py::const_ ),
              R"(
Return the node id and name of component for given DOF

Arguments:
  dof (int): DOF index
  local (bool) = True: if True use local node index else use global index in HPC
Returns:
  tuple[int, str] : node id and name of component
          )",
              py::arg( "dof" ), py::arg( "local" ) = true )
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------
        .def( "getDOFsWithDescription",
              py::overload_cast< const VectorString &, const VectorString &, const bool,
                                 const ASTERINTEGER >(
                  &ParallelEquationNumbering::getDOFsWithDescription, py::const_ ),
              R"(
Get the dofs associated to the given component restricted to the given group.

Arguments:
    cmps (list[str]): components to extract.
    groupNames (list[str]): group names to filter.
    local (bool): if True use local dof index else use global index in HPC
    same_rank : - None: keep all nodes (default: None)
                - True: keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    pair[list[int], list[str]]: list of nodes and list of components
    list[int]: list of dofs
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupNames" ) = VectorString(),
              py::arg( "local" ) = true, py::arg( "same_rank" ) = PythonBool::None )
        .def( "getDOFsWithDescription",
              py::overload_cast< const VectorString &, const VectorLong &, const bool,
                                 const ASTERINTEGER >(
                  &ParallelEquationNumbering::getDOFsWithDescription, py::const_ ),
              R"(
Get the dofs associated to the given component restricted to the given nodes.

Arguments:
    cmps (list[str]): components to extract.
    nodes (list[int]): list of nodes to filter.
    local (bool): if True use local dof index else use global index in HPC
    same_rank : - None: keep all nodes (default: None)
                - True: keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    pair[list[int], list[str]]: list of nodes and list of components.
    list[int]: list of dofs.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "nodes" ) = VectorLong(),
              py::arg( "local" ) = true, py::arg( "same_rank" ) = PythonBool::None );
};

#endif
