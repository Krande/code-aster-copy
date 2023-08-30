/**
 * @file MedMeshInterface.cxx
 * @brief Interface python de MedMesh
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

#include "PythonBindings/MedMeshInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

#ifdef ASTER_HAVE_MED
void exportMedMeshToPython( py::module_ &mod ) {

    py::class_< MedMesh, MedMesh::MedMeshPtr >( mod, "MedMesh" )
        .def( "__pickling_disabled__", disable_pickling< MedMesh >() )

        .def( "getCellFamilyAtSequence", &MedMesh::getCellFamilyAtSequence, R"(
Get cell family in calculation sequence for given profile

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    profile_iterator (int): iterator on profile

Returns:
    list: family id for cells
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "type_iterator" ) )
        .def( "getCellFamilyForGeometricTypeAtSequence",
              &MedMesh::getCellFamilyForGeometricTypeAtSequence, R"(
Get cell family for calculation sequence and geometric type

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geom_type (int): geomtric type

Returns:
    list: family id for cells
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geom_type" ) )
        .def( "getConnectivityAtSequence", &MedMesh::getConnectivityAtSequence, R"(
Get cell connectivity for calculation sequence and geometric type iterator

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtype_iterator (int): iterator on geometric type

Returns:
    list: cell connectivity
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtype_iterator" ) )
        .def( "getConnectivityForGeometricTypeAtSequence",
              &MedMesh::getConnectivityForGeometricTypeAtSequence, R"(
Get cell connectivity for calculation sequence and geometric type

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtype (int): geometric type

Returns:
    list: cell connectivity
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtype" ) )
        .def( "getCellNumberAtSequence", &MedMesh::getCellNumberAtSequence, R"(
Get cell number for calculation sequence and geometric type iterator

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtype_iterator (int): iterator on geometric type

Returns:
    int: cell number
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtype_iterator" ) )
        .def( "getCellNumberForGeometricTypeAtSequence",
              &MedMesh::getCellNumberForGeometricTypeAtSequence, R"(
Get cell number for calculation sequence and geometric type

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtype (int): geometric type

Returns:
    int: cell number
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtype" ) )
        .def( "getCellTypeAtSequence", &MedMesh::getCellTypeAtSequence, R"(
Get cell geometric type for calculation sequence and geomtype_iterator

Arguments:
    numdt (int): time step id
    numit (int): iteration id
    geomtype_iterator (int): iterator on geometric type

Returns:
    int: geometric type
            )",
              py::arg( "numdt" ), py::arg( "numit" ), py::arg( "geomtype_iterator" ) )
        .def( "getCellTypeNumberAtSequence", &MedMesh::getCellTypeNumberAtSequence, R"(
Get cell type number for calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    int: cell type number
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getDimension", &MedMesh::getDimension, R"(
Get mesh dimension )" )
        .def( "getFamilies", &MedMesh::getFamilies, R"(
Get family list

Returns:
    list: MedFamily list
            )" )
        .def( "getGeometricTypesAtSequence", &MedMesh::getGeometricTypesAtSequence, R"(
Get all cell geometric types

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    list: cell geometric type list
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getName", &MedMesh::getName, R"(Get mesh name )" )
        .def( "getNodeFamilyAtSequence", &MedMesh::getNodeFamilyAtSequence, R"(
Get node families for calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    list: node families
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getNodeNumberAtSequence", &MedMesh::getNodeNumberAtSequence, R"(
Get node number for calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    int: node number
            )",
              py::arg( "numdt" ), py::arg( "numit" ) )
        .def( "getNodeNumberForGeometricType", &MedMesh::getNodeNumberForGeometricType, R"(
Get node number from a geometric type

Arguments:
    geotype (int): geometric type

Returns:
    int: node number
            )",
              py::arg( "geotype" ) )
        .def( "getSequence", &MedMesh::getSequence, R"(
Get calculation sequence

Arguments:
    seq_iterator (int): iterator on sequence

Returns:
    list: pair time step id/iterator id
            )",
              py::arg( "seq_iterator" ) )
        .def( "getSequenceNumber", &MedMesh::getSequenceNumber, R"(
Get calculation sequence number )" )
        .def( "readCoordinates", &MedMesh::readCoordinates, R"(
Get coordinates for calculation sequence

Arguments:
    numdt (int): time step id
    numit (int): iteration id

Returns:
    list: coordinates list
            )",
              py::arg( "numdt" ), py::arg( "numit" ) );
};

#endif
