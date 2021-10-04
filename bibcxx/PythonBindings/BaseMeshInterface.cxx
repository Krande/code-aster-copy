/**
 * @file BaseMeshInterface.cxx
 * @brief Interface python de BaseMesh
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

#include "PythonBindings/BaseMeshInterface.h"
#include "PythonBindings/factory.h"
#include <Meshes/BaseMesh.h>

namespace py = boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( print_overloads, printMedFile, 1, 2 )

void exportBaseMeshToPython() {

    py::class_< BaseMesh, BaseMesh::BaseMeshPtr, py::bases< DataStructure > >( "BaseMesh",
                                                                               py::no_init )
        // fake initFactoryPtr: created by subclass
        // fake initFactoryPtr: created by subclass
        .def( "build", &BaseMesh::build, R"(
Build list of Tables based on the mesh

Returns:
    bool: true if building is ok
        )",
              ( py::arg( "self" ) ) )
        .def( "getNumberOfNodes", &BaseMesh::getNumberOfNodes, R"(
Return the number of nodes of the mesh.

Returns:
    int: Number of nodes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getNumberOfCells", &BaseMesh::getNumberOfCells, R"(
Return the number of cells of the mesh.

Returns:
    int: Number of cells.
        )",
              ( py::arg( "self" ) ) )
        .def( "getCoordinates", &BaseMesh::getCoordinates, R"(
Return the coordinates of the mesh.

Returns:
    MeshCoordinatesField: Field of the coordinates.
        )",
              ( py::arg( "self" ) ) )
        .def( "isParallel", &BaseMesh::isParallel, R"(
Tell if the mesh is distributed on parallel instances.

Returns:
    bool: *False* for a centralized mesh, *True* for a parallel mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "getDimension", &BaseMesh::getDimension, R"(
Return the dimension of the mesh.

Returns:
    int: 2 or 3
        )",
              ( py::arg( "self" ) ) )
        .def( "getConnectivity", &BaseMesh::getConnectivity, R"(
Return the connectivity of the mesh as Python lists.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getNodeName", &BaseMesh::getNodeName, R"(
Return the name of the given node

Arguments:
    int : index of the node (1-based)

Returns:
    str : name of the node (stripped)
        )",
              ( py::arg( "self" ), py::arg( "index" ) ) )
        .def( "getCellName", &BaseMesh::getCellName, R"(
Return the name of the given cell

Arguments:
    int : index of the cell (1-based)

Returns:
    str : name of the cell (stripped)
        )",
              ( py::arg( "self" ), py::arg( "index" ) ) )
        .def( "getMedConnectivity", &BaseMesh::getMedConnectivity, R"(
Return the connectivity of the mesh as Python lists following the Med numbering.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )",
              ( py::arg( "self" ) ) )
        .def( "getMedCellsTypes", &BaseMesh::getMedCellsTypes, R"(
Return the Med type of each cell.

Returns:
    list[int]: List of Med types.
        )",
              ( py::arg( "self" ) ) )

        .def( "update", &ListOfTables::update_tables, R"(
Update the internal state of the datastructure.

Returns:
    bool: *True* in case of success, *False* otherwise.
        )",
              ( py::arg( "self" ) ) )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) )
        .def( "printMedFile", &BaseMesh::printMedFile,
              print_overloads( R"(
Print the mesh in the MED format

Arguments:
    filename (str): Name of the file
    local (bool=True) : print local values only (relevent for ParallelMesh only)

Returns:
    Bool: True if of
            )",
        ( py::arg( "self" ), py::arg( "fileName" ), py::arg( "local" ) ) ) );
};

