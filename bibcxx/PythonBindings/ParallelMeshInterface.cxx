/**
 * @file ParallelMeshInterface.cxx
 * @brief Interface python de ParallelMesh
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

#include "PythonBindings/ParallelMeshInterface.h"

#include "aster_pybind.h"

#ifdef ASTER_HAVE_MPI

void exportParallelMeshToPython( py::module_ &mod ) {

    py::class_< ParallelMesh, ParallelMesh::ParallelMeshPtr, BaseMesh >( mod, "ParallelMesh" )
        .def( py::init( &initFactoryPtr< ParallelMesh > ) )
        .def( py::init( &initFactoryPtr< ParallelMesh, std::string > ) )
        .def( "getGroupsOfCells", &ParallelMesh::getGroupsOfCells, R"(
Return the list of the existing (local or global) groups of cells.

Arguments:
    local=false (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfCells", &ParallelMesh::hasGroupOfCells, R"(
The global group exists in the mesh

Arguments:
    group_name (str): Name of the global group.
    local=false (bool): search in local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "getCells", &ParallelMesh::getCells, R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) = "" )
        .def( "getGroupsOfNodes", &ParallelMesh::getGroupsOfNodes, R"(
Return the list of the existing (local or global) groups of nodes.

Arguments:
    local=false (bool): search in local or global groups

Returns:
    list[str]: List of (local or global) groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &ParallelMesh::hasGroupOfNodes, R"(
The (local or global) group exists in the mesh

Arguments:
    group_name (str): Name of the (local or global) group.
    local=false (bool): search local or global groups

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "_getNodes", &ParallelMesh::getNodes,
              R"(
Return the list of the indexes of the nodes that belong to a group of nodes
with (local or global) indexing and a restriction to MPI-rank.

Arguments:
    group_name (str): Name of the group (default: "" = all nodes).
    localNumbering (bool) : use local or global numbering (default: True)
    same_rank : - None: keep all nodes (default: None)
                - True: keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              py::arg( "group_name" ) = std::string(), py::arg( "localNumbering" ) = true,
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "getInnerNodes", &ParallelMesh::getInnerNodes, R"(
Return the list of the indexes of the inner nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "getOuterNodes", &ParallelMesh::getOuterNodes, R"(
Return the list of the indexes of the outer nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "readMedFile", &ParallelMesh::readMedFile, R"(
Read a mesh file from MED format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "getNodesRank", &ParallelMesh::getNodesRank, R"(
Return the rank of the processor which owns the nodes

Returns:
    list[int]: MPI-Rank of the owners of the nodes
        )" )
        .def( "getCellsRank", &ParallelMesh::getCellsRank, R"(
Return the rank of the processor which owns the cells

Returns:
    list[int]: MPI-Rank of the owners of the cells
        )" )
        .def( "_updateGlobalGroupOfCells", &ParallelMesh::updateGlobalGroupOfCells, R"(
Share and update global groups of cells between MPI process.

This function has to be used by developper only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )" )
        .def( "_updateGlobalGroupOfNodes", &ParallelMesh::updateGlobalGroupOfNodes, R"(
Share and update global groups of nodes between MPI process.

This function has to be used by developper only and not user

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )" )
        .def( "_readPartitionedMedFile", &ParallelMesh::readPartitionedMedFile, R"(
Read a partitioned MED file (alaready partitioned by the MEDPartitioner)

This function has to be used by developper only and not user

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "_getNodesFromCells", &ParallelMesh::getNodesFromCells, R"(
Returns the nodes indexes of a group of cells.

Arguments:
    group_name (str): Name of the group.
    localNumbering (bool) : use local or global numbering (default: True)
    same_rank : - None: keep all nodes (default: None)
                - True keep the nodes which are owned by the current MPI-rank
                - False: keep the nodes which are not owned by the current MPI-rank

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              py::arg( "group_name" ), py::arg( "localNumbering" ) = true,
              py::arg( "same_rank" ) = PythonBool::None );
};

#endif /* ASTER_HAVE_MPI */
