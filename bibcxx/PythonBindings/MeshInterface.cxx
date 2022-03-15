/**
 * @file MeshInterface.cxx
 * @brief Interface python de Mesh
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
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "getCells", &Mesh::getCells, R"(
Return the list of the indexes of the cells that belong to a group of cells.

Arguments:
    group_name (str): Name of the local group.

Returns:
    list[int]: Indexes of the cells of the local group.
        )",
              py::arg( "group_name" ) = "" )
        .def( "getGroupsOfNodes", &Mesh::getGroupsOfNodes, R"(
Return the list of the existing groups of nodes.

Arguments:
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[str]: List of groups names (stripped).
        )",
              py::arg( "local" ) = false )
        .def( "hasGroupOfNodes", &Mesh::hasGroupOfNodes, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              py::arg( "group_name" ), py::arg( "local" ) = false )
        .def( "_getNodes", &Mesh::getNodes, R"(
Return the list of the indexes of the nodes that belong to a group of nodes.

Arguments:
    group_name (str): Name of the group (default: "").
    localNumbering (bool): not used (for compatibilty with ParallelMesh)
    same_rank (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              py::arg( "group_name" ) = std::string(), py::arg( "localNumbering" ) = true,
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
        .def( "readMedFile", &Mesh::readMedFile, R"(
Read a mesh file from MED format.

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
        .def( "_getNodesFromCells", &Mesh::getNodesFromCells, R"(
Returns the nodes indexes of a group of cells.

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
