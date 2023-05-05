/**
 * @file MeshInterface.cxx
 * @brief Interface python de Mesh
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

#include "PythonBindings/MeshInterface.h"

#include "aster_pybind.h"

#include <Meshes/BaseMesh.h>
#include <Meshes/Mesh.h>
#include <Meshes/MeshEntities.h>

void exportMeshToPython( py::module_ &mod ) {

    py::class_< Mesh, Mesh::MeshPtr, BaseMesh >( mod, "Mesh" )
        .def( py::init( &initFactoryPtr< Mesh > ) )
        .def( py::init( &initFactoryPtr< Mesh, std::string > ) )
        .def( "getGroupsOfCells", &Mesh::getGroupsOfCells, R"(
Return the list of the existing groups of cells.

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfCells", &Mesh::hasGroupOfCells, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfCells", &Mesh::setGroupOfCells, R"(
Set new group of cells in the mesh

Arguments:
    group_name (str): Name of the new group.
    cell_ids (list[int]) : cell ids which are in the group
        )",
              py::arg( "group_name" ), py::arg( "cell_ids" ) )
        .def( "getCells", py::overload_cast< const std::string >( &Mesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) )
        .def( "getCells", py::overload_cast< const VectorString & >( &Mesh::getCells, py::const_ ),
              R"(
Return the list of the indexes of the cells that belong to the groups of cells.

Arguments:
    groups_name (str): Name of the local groups.

Returns:
    list[int]: Indexes of the cells of the local groups.
        )",
              py::arg( "groups_name" ) = VectorString() )
        .def( "getGroupsOfNodes", &Mesh::getGroupsOfNodes, R"(
Return the list of the existing groups of nodes.

Arguments:
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &Mesh::hasGroupOfNodes, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "setGroupOfNodes", &Mesh::setGroupOfNodes, R"(
Set new group of nodes in the mesh

Arguments:
    group_name (str): Name of the new group.
    node_ids (list[int]) : node ids which are in the group
    localNumbering=false (bool): not used (for compatibilty with ParallelMesh)
        )",
              py::arg( "group_name" ), py::arg( "node_ids" ), py::arg( "localNumbering" ) = false )
        .def( "_getNodes",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &Mesh::getNodes, py::const_ ),
              R"(
Return the list of the indexes of the nodes that belong to groups of nodes.
* For internal use only *

Arguments:
    group_name (list[str]): Name of groups (default: "").
    localNumbering (bool): not used (for compatibilty with ParallelMesh)
    same_rank (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of groups.
        )",
              py::arg( "group_name" ) = VectorString(), py::arg( "localNumbering" ) = true,
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "getInnerNodes", &Mesh::getInnerNodes, R"(
Return the list of the indexes of the nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )" )
        .def( "readAsterFile", &Mesh::readAsterFile, R"(
Read a mesh file from ASTER format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "readGibiFile", &Mesh::readGibiFile, R"(
Read a mesh file from GIBI format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "readGmshFile", &Mesh::readGmshFile, R"(
Read a mesh file from GMSH format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              py::arg( "filename" ) )
        .def( "isQuadratic", &Mesh::isQuadratic, R"(
To know if the mesh contains quadratic cells

Returns:
    bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        )" )
        .def( "_getNodesFromCells",
              py::overload_cast< const VectorString &, const bool, const ASTERINTEGER >(
                  &Mesh::getNodesFromCells, py::const_ ),
              R"(
Returns the nodes indexes of a group of cells.
For developpers only.

Arguments:
    group_name (str): name of the group of cells.
    localNumbering (bool): not used (for compatibilty with ParallelMesh)
    same_rank (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: indexes of the nodes.
        )",
              py::arg( "group_name" ), py::arg( "localNumbering" ) = true,
              py::arg( "same_rank" ) = PythonBool::None );
};
