/**
 * @file BaseMeshInterface.cxx
 * @brief Interface python de BaseMesh
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

#include "PythonBindings/BaseMeshInterface.h"

#include "aster_pybind.h"

#include <Meshes/BaseMesh.h>

void exportBaseMeshToPython( py::module_ &mod ) {

    py::class_< BaseMesh, BaseMesh::BaseMeshPtr, DataStructure >( mod, "BaseMesh" )
        // fake initFactoryPtr: created by subclass
        // fake initFactoryPtr: created by subclass
        .def( "build", &BaseMesh::build, R"(
Build list of Tables based on the mesh

Returns:
    bool: true if building is ok
        )" )
        .def( "getNumberOfNodes", &BaseMesh::getNumberOfNodes, R"(
Return the number of nodes of the mesh.

Returns:
    int: Number of nodes.
        )" )
        .def( "getNumberOfCells", &BaseMesh::getNumberOfCells, R"(
Return the number of cells of the mesh.

Returns:
    int: Number of cells.
        )" )
        .def( "getCoordinates", &BaseMesh::getCoordinates, R"(
Return the coordinates of the mesh.

Returns:
    MeshCoordinatesField: Field of the coordinates.
        )" )
        .def( "isParallel", &BaseMesh::isParallel, R"(
Tell if the mesh is distributed on parallel instances.

Returns:
    bool: *False* for a centralized mesh, *True* for a parallel mesh.
        )" )
        .def( "getDimension", &BaseMesh::getDimension, R"(
Return the dimension of the mesh.

Returns:
    int: 2 or 3
        )" )
        .def( "getConnectivity", &BaseMesh::getConnectivity, R"(
Return the connectivity of the mesh as Python lists.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )" )
        .def( "getNodeName", &BaseMesh::getNodeName, R"(
Return the name of the given node

Arguments:
    index (int) : index of the node (0-based)

Returns:
    str : name of the node (stripped)
        )",
              py::arg( "index" ) )
        .def( "getCellName", &BaseMesh::getCellName, R"(
Return the name of the given cell

Arguments:
    index (int) : index of the cell (0-based)

Returns:
    str : name of the cell (stripped)
        )",
              py::arg( "index" ) )
        .def( "getMedConnectivity", &BaseMesh::getMedConnectivity, R"(
Return the connectivity of the mesh as Python lists following the Med numbering.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )" )
        .def( "getMedCellsTypes", &BaseMesh::getMedCellsTypes, R"(
Return the Med type of each cell.

Returns:
    list[int]: List of Med types.
        )" )

        .def( "update", &ListOfTables::update_tables, R"(
Update the internal state of the datastructure.

Returns:
    bool: *True* in case of success, *False* otherwise.
        )" )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) )
        .def( "printMedFile", &BaseMesh::printMedFile, R"(
Print the mesh in the MED format

Arguments:
    filename (str): Name of the file
    local (bool=True) : print local values only (relevent for ParallelMesh only)

Returns:
    Bool: True if of
            )",
              py::arg( "fileName" ), py::arg( "local" ) = true );
};
