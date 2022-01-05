/**
 * @file MeshInterface.cxx
 * @brief Interface python de Mesh
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

#include "PythonBindings/MeshInterface.h"
#include <Meshes/BaseMesh.h>
#include <Meshes/Mesh.h>
#include <Meshes/MeshEntities.h>
#include <PythonBindings/factory.h>

namespace py = boost::python;

void exportMeshToPython() {

     VectorString ( Mesh::*c1 )(  ) const =
        &Mesh::getGroupsOfCells;
    VectorString ( Mesh::*c2 )( const bool ) const =
        &Mesh::getGroupsOfCells;

    VectorString ( Mesh::*c3 )(  ) const =
        &Mesh::getGroupsOfNodes;
    VectorString ( Mesh::*c4 )( const bool ) const =
        &Mesh::getGroupsOfNodes;

    bool ( Mesh::*c5 )( const std::string&  ) const =
        &Mesh::hasGroupOfCells;
    bool ( Mesh::*c6 )( const std::string&, const bool ) const =
        &Mesh::hasGroupOfCells;

    bool ( Mesh::*c7 )( const std::string& ) const =
        &Mesh::hasGroupOfNodes;
    bool ( Mesh::*c8 )( const std::string&, const bool ) const =
        &Mesh::hasGroupOfNodes;

    VectorLong ( Mesh::*c9 )(   ) const =
        &Mesh::getCells;
    VectorLong ( Mesh::*c10 )( const std::string ) const =
        &Mesh::getCells;

    VectorLong ( Mesh::*n1 )(   ) const =
        &Mesh::getNodes;
    VectorLong ( Mesh::*n2 )( const std::string ) const =
        &Mesh::getNodes;
    VectorLong ( Mesh::*n3 )( const std::string, const bool  ) const =
        &Mesh::getNodes;
    VectorLong ( Mesh::*n4 )( const std::string, const bool, const bool ) const =
        &Mesh::getNodes;
    VectorLong ( Mesh::*n5 )( const bool  ) const =
        &Mesh::getNodes;
    VectorLong ( Mesh::*n6 )( const bool, const bool ) const =
        &Mesh::getNodes;


    py::class_< Mesh, Mesh::MeshPtr, py::bases< BaseMesh > >( "Mesh", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< Mesh > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< Mesh, std::string > ) )
                .def( "getGroupsOfCells", c1, R"(
Return the list of the existing groups of cells.

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "getGroupsOfCells", c2, R"(
Return the list of the existing  groups of cells.

Arguments:
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" ), py::arg("local") ) )
        .def( "hasGroupOfCells", c5, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg("group_name") ) )
        .def( "hasGroupOfCells", c6, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local=false (bool): not used (for compatibilty with ParallelMesh)

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
Return the list of the existing groups of nodes.

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" )) )
        .def( "getGroupsOfNodes", c4,  R"(
Return the list of the existing groups of nodes.

Arguments:
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    list[str]: List of groups names (stripped).
        )",
              ( py::arg( "self" ), py::arg("local") ) )
        .def( "hasGroupOfNodes", c7, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg("group_name") ) )
        .def( "hasGroupOfNodes", c8, R"(
The group exists in the mesh

Arguments:
    group_name (str): Name of the group.
    local=false (bool): not used (for compatibilty with ParallelMesh)

Returns:
    bool: *True* if exists, *False* otherwise.
        )",
              ( py::arg( "self" ), py::args("group_name", "local") ) )
        .def( "getNodes", n1, R"(
Return the list of the indexes of the nodes in the mesh.

Returns:
    list[int]: Indexes of the nodes in the mesh.
        )",
              ( py::arg( "self" ) )
               )
        .def( "getNodes", n2, R"(
Return the list of the indexes of the nodes that belong to a group of nodes.

Arguments:
    group_name (str): Name of the group.

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              ( py::arg( "self" ), py::arg("group_name") )
               )
        .def( "getNodes", n3, R"(
Return the list of the indexes of the nodes that belong to a group of nodes.

Arguments:
    group_name (str): Name of the group.
    localNumbering=true (bool) : not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              ( py::arg( "self" ), py::args("group_name", "localNumbering") )
               )
        .def( "getNodes", n4, R"(
Return the list of the indexes of the nodes that belong to a group of nodes.

Arguments:
    group_name (str): Name of the group.
    localNumbering (bool) : not used (for compatibilty with ParallelMesh)
    same_rank : not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              ( py::arg( "self" ), py::args("group_name", "localNumbering", "same_rank") )
               )
        .def( "getNodes", n5, R"(
Return the list of the indexes of the nodes in the mesh.

Arguments:
    localNumbering (bool) : not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              ( py::arg( "self" ), py::arg("localNumbering") )
               )
        .def( "getNodes", n6, R"(
Return the list of the indexes of the nodes in the mesh

Arguments:
    localNumbering (bool) : not used (for compatibilty with ParallelMesh)
    same_rank : not used (for compatibilty with ParallelMesh)

Returns:
    list[int]: Indexes of the nodes of the group.
        )",
              ( py::arg( "self" ), py::args("localNumbering", "same_rank") )
               )
        .def( "getInnerNodes", &Mesh::getInnerNodes, R"(
Return the list of the indexes of the nodes in the mesh

Returns:
    list[int]: Indexes of the nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "readAsterFile", &Mesh::readAsterFile, R"(
Read a mesh file from ASTER format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg( "filename" ) ) )
        .def( "readGibiFile", &Mesh::readGibiFile, R"(
Read a mesh file from GIBI format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg( "filename" ) ) )
        .def( "readGmshFile", &Mesh::readGmshFile, R"(
Read a mesh file from GMSH format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg( "filename" ) ) )
        .def( "readMedFile", &Mesh::readMedFile, R"(
Read a mesh file from MED format.

Arguments:
    filename (str): Path to the file to be read.

Returns:
    bool: *True* if succeeds, *False* otherwise.
        )",
              ( py::arg( "self" ), py::arg( "filename" ) ) )
        .def( "isQuadratic", &Mesh::isQuadratic, R"(
To know if the mesh contains quadratic cells

Returns:
    bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        )",
              ( py::arg( "self" )  ) );
};
