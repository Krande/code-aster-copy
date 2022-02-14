/**
 * @file ContactZoneInterface.cxx
 * @brief Interface python de ContactZone
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

#include "PythonBindings/ContactZoneInterface.h"

#include "aster_pybind.h"

void exportContactZoneToPython( py::module_ &mod ) {

    py::class_< ContactZone, ContactZone::ContactZonePtr, DataStructure >( mod, "ContactZone" )
        .def( py::init( &initFactoryPtr< ContactZone, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ContactZone, ModelPtr > ) )
        .def( "getModel", &ContactZone::getModel, R"(
Return the model used in the contact zone definition

Returns:
    Model: model.
        )" )
        .def( "getMesh", &ContactZone::getMesh, R"(
Return the mesh used in the contact zone definition

Returns:
    BaseMesh: mesh.
        )" )
        .def( "setVerbosity", &ContactZone::setVerbosity, R"(
Set level of verbosity:
      0- without
      1- normal (default)
      2- detailled

Arguments:
    integer: level of verbosity
        )",
              py::arg( "level" ) )
        .def( "getVerbosity", &ContactZone::getVerbosity, R"(
Get level of verbosity:
      0- without
      1- normal
      2- detailled

Returns:
    integer: level of verbosity
        )" )
        .def( "build", &ContactZone::build, R"(
Build and check internal objects
        )" )
        .def( "setContactParameter", &ContactZone::setContactParameter, R"(
Set contact parameters defining method, coefficient...

Arguments:
    ContactParameter: contact parameters
        )",
              py::arg( "contact" ) )
        .def( "getContactParameter", &ContactZone::getContactParameter, R"(
Get contact parameters defining method, coefficient...

Returns:
    ContactParameter: contact parameters
        )" )
        .def( "setFrictionParameter", &ContactZone::setFrictionParameter, R"(
Set friction parameters defining method, coefficient...

Arguments:
    FrictionParameter: friction parameters
        )",
              py::arg( "friction" ) )
        .def( "getFrictionParameter", &ContactZone::getFrictionParameter, R"(
Get friction parameters defining method, coefficient...

Returns:
    FrictionParameter: friction parameters
        )" )
        .def( "setPairingParameter", &ContactZone::setPairingParameter, R"(
Set pairing parameters defining algorithm, distance...

Arguments:
    PairingParameter: pairing parameters
        )",
              py::arg( "pairing" ) )
        .def( "getPairingParameter", &ContactZone::getPairingParameter, R"(
Get pairing parameters defining algorithm, distance...

Returns:
    PairingParameter: pairing parameters
        )" )
        .def( "setSlaveGroupOfCells", &ContactZone::setSlaveGroupOfCells, R"(
Set slave's name of group of cells

Arguments:
    str: slave's name
        )",
              py::arg( "slave_name" ) )
        .def( "getSlaveGroupOfCells", &ContactZone::getSlaveGroupOfCells, R"(
Get slave's name of group of cells

Returns:
    str: slave's name
        )" )
        .def( "setMasterGroupOfCells", &ContactZone::setMasterGroupOfCells, R"(
Set master's name of group of cells

Arguments:
    str: master's name
        )",
              py::arg( "master_name" ) )
        .def( "getMasterGroupOfCells", &ContactZone::getMasterGroupOfCells, R"(
Get master's name of group of cells

Returns:
    str: master's name
        )" )
        .def( "setExcludedSlaveGroupOfCells", &ContactZone::setExcludedSlaveGroupOfCells, R"(
Set excluded groups of cells on slave side

Arguments:
    str: excluded groups' names
        )",
              py::arg( "master_name" ) )
        .def( "getExcludedSlaveGroupOfCells", &ContactZone::getExcludedSlaveGroupOfCells, R"(
Get excluded groups of cells on slave side

Returns:
    str: excluded groups' names
        )" )
        .def_property( "checkNormals",
                       py::overload_cast<>( &ContactZone::checkNormals, py::const_ ),
                       py::overload_cast< const bool & >( &ContactZone::checkNormals ), R"(
        bool: Attribute that holds the checking of outwards normals.
                )" )
        .def( "updateSlaveCells", &ContactZone::updateSlaveCells, R"(
Returns:
      Bool: True if the slave cells are updated
        )" )
        .def( "updateMasterCells", &ContactZone::updateMasterCells, R"(
Returns:
      Bool: True if checking is performed else False
        )" )
        .def( "getMasterCellsFromNode", &ContactZone::getMasterCellsFromNode, R"(
Get the master cells associtaed with a node number

Arguments:
    int: node number
        )",
              py::arg( "node_number" ) )
        .def( "getSlaveCellsFromNode", &ContactZone::getSlaveCellsFromNode, R"(
Get the slave cells associtaed with a node number

Arguments:
    int: node number
        )",
              py::arg( "node_number" ) )
        .def( "getMasterCellNeighbors", &ContactZone::getMasterCellNeighbors, R"(
Get the master cells in the neighbor of a given master cell number

Arguments:
    int: master cell number
        )",
              py::arg( "cell_number" ) )
        .def( "getSlaveCellNeighbors", &ContactZone::getSlaveCellNeighbors, R"(
Get the slave cells in the neighbor of a given slave cell number

Arguments:
    int: slave cell number
        )",
              py::arg( "cell_number" ) );
};
