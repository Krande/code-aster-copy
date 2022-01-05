/**
 * @file ConnectionMeshInterface.cxx
 * @brief Interface python de ConnectionMesh
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/factory.h"
#include "PythonBindings/ConnectionMeshInterface.h"

void exportConnectionMeshToPython() {

#ifdef ASTER_HAVE_MPI

    VectorLong ( ConnectionMesh::*c1 )(   ) const =
        &ConnectionMesh::getCells;
    VectorLong ( ConnectionMesh::*c2 )( const std::string ) const =
        &ConnectionMesh::getCells;

    py::class_< ConnectionMesh, ConnectionMesh::ConnectionMeshPtr,
                py::bases< BaseMesh > >( "ConnectionMesh", py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ConnectionMesh, ParallelMeshPtr,
                                    VectorString, VectorString >))
        .def( "__init__", py::make_constructor(&initFactoryPtr< ConnectionMesh, std::string,
                                                                ParallelMeshPtr, VectorString,
                                                                VectorString >))
        .def( "getGroupsOfCells", &ConnectionMesh::getGroupsOfCells, R"(
Return the list of the existing groups of cells.

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "getGroupsOfNodes", &ConnectionMesh::getGroupsOfNodes, R"(
Return the list of the existing groups of nodes.

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "hasGroupOfCells", &ConnectionMesh::hasGroupOfCells, R"(
Allows to know if the given group of cells is present in the mesh

Arguments:
    str: name of the group of cell

Returns:
    bool: True if the group is present
        )",
              ( py::arg( "self" ), py::arg( "name")) )
        .def( "hasGroupOfNodes", &ConnectionMesh::hasGroupOfNodes, R"(
Allows to know if the given group of nodes is present in the mesh

Arguments:
    str: name of the group of nodes

Returns:
    bool: True if the group is present
        )",
              ( py::arg( "self" ), py::arg( "name")) )
        .def( "getCells", c1,  R"(
Return the list of the indexes of the cells in mesh.

Returns:
    list[int]: Indexes of the cells in the mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "getCells", c2,  R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              ( py::arg( "self" ), py::arg( "group_name" ) ) )
        .def( "getParallelMesh", &ConnectionMesh::getParallelMesh,
            py::return_value_policy<py::copy_const_reference>(), R"(
Return a pointer to the ParallelMesh used to built it.

Returns:
    ParallelMeshPtr: pointer to the ParallelMesh
        )",
              ( py::arg( "self" )) )
        .def( "getNodesGlobalNumbering", &ConnectionMesh::getNodesGlobalNumbering,
            py::return_value_policy<py::copy_const_reference>(), R"(
Return a tuple of the nodes of the mesh with a global numbering

Returns:
    tuple[int]: list of nodes with global numbering
        )",
              ( py::arg( "self" )) )
        .def( "getNodesLocalNumbering", &ConnectionMesh::getNodesLocalNumbering,
            py::return_value_policy<py::copy_const_reference>(), R"(
Return a tuple of the nodes of the mesh with a local numbering.
The local numbering is the one coming from the owner of the node,
hence some nodes can have the same local numbering

Returns:
    tuple[int]: list of nodes with local numbering
        )",
              ( py::arg( "self" )) );
#endif /* ASTER_HAVE_MPI */
};
