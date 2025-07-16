/**
 * @file ContactZoneInterface.cxx
 * @brief Interface python de ContactZone
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

    py::class_< ContactZone, ContactZone::ContactZonePtr, DSWithCppPickling >( mod, "ContactZone",
                                                                               R"(
Object to define a zone of contact.)" )
        .def( py::init( &initFactoryPtr< ContactZone, std::string > ) )
        .def( py::init( &initFactoryPtr< ContactZone > ) )
        .def( py::init( &initFactoryPtr< ContactZone, const py::tuple & > ) )
        .def( define_pickling< ContactZone >() )
        .def( "getModel", &ContactZone::getModel, R"(
Return the model used in the contact zone definition

Returns:
    Model: model
)" )
        .def( "getMesh", &ContactZone::getMesh, R"(
Return the mesh used in the contact zone definition

Returns:
    BaseMesh: mesh
)" )
        .def( "setVerbosity", &ContactZone::setVerbosity, R"(
Set level of verbosity
0- without
1- normal (default)
2- detailled

Arguments:
    level (int) : level of verbosity
)",
              py::arg( "level" ) )
        .def( "getVerbosity", &ContactZone::getVerbosity, R"(
Get level of verbosity
0- without
1- normal
2- detailled

Returns:
    int: level of verbosity
)" )
        .def_property( "hasFriction", &ContactZone::hasFriction, &ContactZone::enableFriction,
                       R"(
Get status of friction

Returns:
    bool: friction or not
        )" )
        .def_property( "hasSmoothing", &ContactZone::hasSmoothing, &ContactZone::enableSmoothing,
                       R"(
Smoothing of normals

Returns:
    bool: smoothing or not
            )" )
        .def( "build", &ContactZone::build, R"(
Build and check internal objects

Returns:
    bool: success or failure
            )" )
        .def( "setContactParameter", &ContactZone::setContactParameter, R"(
Set contact parameters defining method, coefficient...

Arguments:
    contParam (ContactParameter) : contact parameters
)",
              py::arg( "contParam" ) )
        .def( "getContactParameter", &ContactZone::getContactParameter, R"(
Get contact parameters defining method, coefficient...

Returns:
    ContactParameter: contact parameters
)" )
        .def( "setFrictionParameter", &ContactZone::setFrictionParameter, R"(
Set friction parameters defining method, coefficient...

Arguments:
    fricParam (FrictionParameter) : friction parameters
)",
              py::arg( "fricParam" ) )
        .def( "getFrictionParameter", &ContactZone::getFrictionParameter, R"(
Get friction parameters defining method, coefficient...

Returns:
    FrictionParameter: friction parameters
)" )
        .def( "setPairingParameter", &ContactZone::setPairingParameter, R"(
Set pairing parameters defining algorithm, distance...

Arguments:
    pairParam (PairingParameter) : pairing parameters
)",
              py::arg( "pairParam" ) )
        .def( "getPairingParameter", &ContactZone::getPairingParameter, R"(
Get pairing parameters defining algorithm, distance...

Returns:
    PairingParameter: pairing parameters
)" )
        .def( "getSlaveNodes", &ContactZone::getSlaveNodes, R"(
Get slave's nodes index

Returns:
    list[int]: slave's nodes index
)" )
        .def( "getSlaveCells", &ContactZone::getSlaveCells, R"(
Get slave's cells index

Returns:
    list[int]: slave's cells index
)" )
        .def( "setSlaveGroupOfCells", &ContactZone::setSlaveGroupOfCells, R"(
Set slave's name of group of cells

Arguments:
    slave_name (str) : name of group for slave cells
)",
              py::arg( "slave_name" ) )
        .def( "setMasterGroupOfCells", &ContactZone::setMasterGroupOfCells, R"(
Set master's name of group of cells

Arguments:
    master_name (str) : name of group for master cells
)",
              py::arg( "master_name" ) )
        .def( "setExcludedSlaveGroupOfNodes", &ContactZone::setExcludedSlaveGroupOfNodes,
              R"(
Set excluded groups of nodes on slave side

Arguments:
    nodeGroupsName (str) : excluded groups' names
                )",
              py::arg( "nodeGroupsName" ) )
        .def( "setExcludedSlaveGroupOfCells", &ContactZone::setExcludedSlaveGroupOfCells, R"(
Set excluded groups of cells on slave side

Arguments:
    cellGroupsName (str) : excluded groups' names
                )",
              py::arg( "cellGroupsName" ) )
        .def( "getMeshPairing", &ContactZone::getMeshPairing, R"(
Get pairing of surface meshes

Returns:
    MeshPairing: mesh pairing
        )" )
        .def_property( "checkNormals",
                       py::overload_cast<>( &ContactZone::checkNormals, py::const_ ),
                       py::overload_cast< const bool & >( &ContactZone::checkNormals ), R"(
bool: attribute that holds the checking of outwards normals.
        )" );
};
