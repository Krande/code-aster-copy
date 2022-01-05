/**
 * @file ParallelMeshInterface.cxx
 * @brief Interface python de ParallelMesh
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
#include <PythonBindings/factory.h>
#include "PythonBindings/ParallelMeshInterface.h"

#ifdef ASTER_HAVE_MPI

void exportParallelMeshToPython() {

    VectorString ( ParallelMesh::*c1 )(  ) const =
        &ParallelMesh::getGroupsOfCells;
    VectorString ( ParallelMesh::*c2 )( const bool ) const =
        &ParallelMesh::getGroupsOfCells;

    VectorString ( ParallelMesh::*c3 )(  ) const =
        &ParallelMesh::getGroupsOfNodes;
    VectorString ( ParallelMesh::*c4 )( const bool ) const =
        &ParallelMesh::getGroupsOfNodes;

    bool ( ParallelMesh::*c5 )( const std::string&  ) const =
        &ParallelMesh::hasGroupOfCells;
    bool ( ParallelMesh::*c6 )( const std::string&, const bool ) const =
        &ParallelMesh::hasGroupOfCells;

    bool ( ParallelMesh::*c7 )( const std::string& ) const =
        &ParallelMesh::hasGroupOfNodes;
    bool ( ParallelMesh::*c8 )( const std::string&, const bool ) const =
        &ParallelMesh::hasGroupOfNodes;

    VectorLong ( ParallelMesh::*c9 )(   ) const =
        &ParallelMesh::getCells;
    VectorLong ( ParallelMesh::*c10 )( const std::string ) const =
        &ParallelMesh::getCells;

    VectorLong ( ParallelMesh::*n1 )(   ) const =
        &ParallelMesh::getNodes;
    VectorLong ( ParallelMesh::*n2 )( const std::string ) const =
        &ParallelMesh::getNodes;
    VectorLong ( ParallelMesh::*n3 )( const std::string, const bool  ) const =
        &ParallelMesh::getNodes;
    VectorLong ( ParallelMesh::*n4 )( const std::string, const bool, const bool ) const =
        &ParallelMesh::getNodes;
    VectorLong ( ParallelMesh::*n5 )( const bool  ) const =
        &ParallelMesh::getNodes;
    VectorLong ( ParallelMesh::*n6 )( const bool, const bool ) const =
        &ParallelMesh::getNodes;


    py::class_< ParallelMesh, ParallelMesh::ParallelMeshPtr,
                py::bases< BaseMesh > >( "ParallelMesh", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< ParallelMesh >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ParallelMesh, std::string >))
        .def( "getGroupsOfCells", c1, R"(
Return the list of the existing global groups of cells.

Returns:
    list[str]: List of global groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "getGroupsOfCells", c2, R"(
Return the list of the existing (local or global) groups of cells.

Arguments:
    local=false (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              ( py::arg( "self" ), py::arg("local") ) )
        .def( "hasGroupOfCells", c5, R"(
The global group exists in the mesh

Arguments:
    group_name (str): Name of the global group.

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg("group_name") ) )
        .def( "hasGroupOfCells", c6, R"(
The (local or global) group exists in the mesh

Arguments:
    group_name (str): Name of the (local or global) group.
    local=false (bool): search in local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::args("group_name", "local") ) )
        .def( "getCells", c9,  R"(
Return the list of the indexes of the cells in mesh.

Returns:
    list[int]: Indexes of the cells in the mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "getCells", c10,  R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              ( py::arg( "self" ), py::arg( "group_name" ) ) )
        .def( "getGroupsOfNodes", c3,  R"(
Return the list of the existing global groups of nodes.

Returns:
    list[str]: List of global groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "getGroupsOfNodes", c4,  R"(
Return the list of the existing (local or global) groups of nodes.

Arguments:
    local=false (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              ( py::arg( "self" ), py::arg("local") ) )
        .def( "hasGroupOfNodes", c7, R"(
The global group exists in the mesh

Arguments:
    group_name (str): Name of the global group.

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg("group_name") ) )
        .def( "hasGroupOfNodes", c8, R"(
The (local or global) group exists in the mesh

Arguments:
    group_name (str): Name of the (local or global) group.
    local=false (bool): search local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::args("group_name", "local") ) )
        .def( "getNodes", n1, R"(
Return the list of the indexes of the nodes in the mesh with local indexing.

Returns:
    list[int]: Indexes of the nodes in the mesh with local indexing.
        )",
              ( py::arg( "self" ) )
               )
        .def( "getNodes", n2, R"(
Return the list of the indexes of the nodes that belong to a group of nodes with local indexing.

Arguments:
    group_name (str): Name of the group.

Returns:
    list[int]: Indexes of the nodes of the group with local indexing.
        )",
              ( py::arg( "self" ), py::arg("group_name") )
               )
        .def( "getNodes", n3, R"(
Return the list of the indexes of the nodes that belong to a group of nodes
with (local or global) indexing.

Arguments:
    group_name (str): Name of the group.
    localNumbering=true (bool) : use local or global numbering

Returns:
    list[int]: Indexes of the nodes of the group with (local or global) indexing.
        )",
              ( py::arg( "self" ), py::args("group_name", "localNumbering") )
               )
        .def( "getNodes", n4, R"(
Return the list of the indexes of the nodes that belong to a group of nodes
with (local or global) indexing and a restriction to MPI-rank.

Arguments:
    group_name (str): Name of the group.
    localNumbering (bool) : use local or global numbering
    same_rank : keep or not the nodes which are owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group with (local or global) indexing.
        )",
              ( py::arg( "self" ), py::args("group_name", "localNumbering", "same_rank") )
               )
        .def( "getNodes", n5, R"(
Return the list of the indexes of the nodes in the mesh with (local or global) indexing.

Arguments:
    localNumbering (bool) : use local or global numbering

Returns:
    list[int]: Indexes of the nodes of the group with (local or global) indexing.
        )",
              ( py::arg( "self" ), py::arg("localNumbering") )
               )
        .def( "getNodes", n6, R"(
Return the list of the indexes of the nodes in the mesh with (local or global) indexing
and a restriction to MPI-rank

Arguments:
    localNumbering (bool) : use local or global numbering
    same_rank : keep or not the nodes which are owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group with (local or global) indexing
and a restriction to MPI-rank.
        )",
              ( py::arg( "self" ), py::args("localNumbering", "same_rank") )
               )
        .def( "getInnerNodes", &ParallelMesh::getInnerNodes, R"(
Return the list of the indexes of the inner nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getOuterNodes", &ParallelMesh::getOuterNodes, R"(
Return the list of the indexes of the outer nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "readMedFile", &ParallelMesh::readMedFile, R"(
Read a mesh file from MED format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg( "filename" ) ) )
        .def( "getNodesRank", &ParallelMesh::getNodesRank, R"(
Return the rank of the processor which owns the nodes

Returns:
    list[int]: MPI-Rank of the owners of the nodes
        )",
              ( py::arg( "self" ) ) )
        .def( "getCellsRank", &ParallelMesh::getCellsRank, R"(
Return the rank of the processor which owns the cells

Returns:
    list[int]: MPI-Rank of the owners of the cells
        )",
              ( py::arg( "self" ) ) )
        .def( "_updateGlobalGroupOfCells", &ParallelMesh::updateGlobalGroupOfCells, R"(
Share and update global groups of cells between MPI process.

This function has to be used by developper only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" )) )
        .def( "_updateGlobalGroupOfNodes", &ParallelMesh::updateGlobalGroupOfNodes, R"(
Share and update global groups of nodes between MPI process.

This function has to be used by developper only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" )) )
        .def( "_readPartitionedMedFile", &ParallelMesh::readPartitionedMedFile, R"(
Read a partitioned MED file (alaready partitioned by the MEDPartitioner)

This function has to be used by developper only and not user

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg("filename")) ) ;
};

#endif /* ASTER_HAVE_MPI */
