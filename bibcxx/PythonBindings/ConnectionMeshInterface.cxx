/**
 * @file ConnectionMeshInterface.cxx
 * @brief Interface python de ConnectionMesh
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

#include "PythonBindings/ConnectionMeshInterface.h"

#include "aster_pybind.h"

void exportConnectionMeshToPython( py::module_ &mod ) {

#ifdef ASTER_HAVE_MPI

    py::class_< ConnectionMesh, ConnectionMesh::ConnectionMeshPtr, BaseMesh >( mod,
                                                                               "ConnectionMesh" )
        .def( py::init(
            &initFactoryPtr< ConnectionMesh, ParallelMeshPtr, VectorString, VectorString > ) )
        .def( py::init( &initFactoryPtr< ConnectionMesh, std::string, ParallelMeshPtr, VectorString,
                                         VectorString > ) )
        .def( "getGroupsOfCells", &ConnectionMesh::getGroupsOfCells, R"(
Return the list of the existing groups of cells.

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "getGroupsOfNodes", &ConnectionMesh::getGroupsOfNodes, R"(
Return the list of the existing groups of nodes.

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfCells", &ConnectionMesh::hasGroupOfCells, R"(
Allows to know if the given group of cells is present in the mesh

Arguments:
    name (str): name of the group of cell

Returns:
    bool: True if the group is present
        )",
              py::arg( "name" ), py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &ConnectionMesh::hasGroupOfNodes, R"(
Allows to know if the given group of nodes is present in the mesh

Arguments:
    name (str): name of the group of nodes

Returns:
    bool: True if the group is present
        )",
              py::arg( "name" ), py::arg( "local" ) = false )
        .def( "getCells",
              py::overload_cast< const std::string >( &ConnectionMesh::getCells, py::const_ ), R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) = "" )
        .def( "getParallelMesh", &ConnectionMesh::getParallelMesh, R"(
Return a pointer to the ParallelMesh used to built it.

Returns:
    ParallelMeshPtr: pointer to the ParallelMesh
        )" )
        .def( "getNodesGlobalNumbering", &ConnectionMesh::getNodesGlobalNumbering, R"(
Return a tuple of the nodes of the mesh with a global numbering

Returns:
    tuple[int]: list of nodes with global numbering
        )" )
        .def( "getNodesLocalNumbering", &ConnectionMesh::getNodesLocalNumbering, R"(
Return a tuple of the nodes of the mesh with a local numbering.
The local numbering is the one coming from the owner of the node,
hence some nodes can have the same local numbering

Returns:
    tuple[int]: list of nodes with local numbering
        )" )
        .def( "isConnection", &ConnectionMesh::isConnection, R"(
Function to know if a mesh is a ConnectionMesh)" );
#endif /* ASTER_HAVE_MPI */
};
