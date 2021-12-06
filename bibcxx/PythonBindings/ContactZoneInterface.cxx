/**
 * @file ContactZoneInterface.cxx
 * @brief Interface python de ContactZone
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/ContactZoneInterface.h"
#include <PythonBindings/factory.h>

void exportContactZoneToPython() {

    py::class_< ContactZone, ContactZone::ContactZonePtr, py::bases< DataStructure > >(
        "ContactZone", py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ContactZone, std::string, ModelPtr > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< ContactZone, ModelPtr > ) )
        .def( "getModel", &ContactZone::getModel, R"(
Return the model used in the contact zone definition

Returns:
    Model: model.
        )",
              ( py::arg( "self" ) ) )
        .def( "getMesh", &ContactZone::getMesh, R"(
Return the mesh used in the contact zone definition

Returns:
    BaseMesh: mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "setVerbosity", &ContactZone::setVerbosity, R"(
Set level of verbosity:
      0- without
      1- normal (default)
      2- detailled

Arguments:
    integer: level of verbosity
        )",
              ( py::arg( "self" ), py::arg( "level" ) ) )
        .def( "getVerbosity", &ContactZone::getVerbosity, R"(
Get level of verbosity:
      0- without
      1- normal
      2- detailled

Returns:
    integer: level of verbosity
        )",
              ( py::arg( "self" ) ) )
        .def( "build", &ContactZone::build, R"(
Build and check internal objects
        )",
              ( py::arg( "self" ) ) )
        .def( "setContactParameter", &ContactZone::setContactParameter, R"(
Set contact parameters defining method, coefficient...

Arguments:
    ContactParameter: contact parameters
        )",
              ( py::arg( "self" ), py::arg( "contact" ) ) )
        .def( "getContactParameter", &ContactZone::getContactParameter, R"(
Get contact parameters defining method, coefficient...

Returns:
    ContactParameter: contact parameters
        )",
              ( py::arg( "self" ) ) )
        .def( "setFrictionParameter", &ContactZone::setFrictionParameter, R"(
Set friction parameters defining method, coefficient...

Arguments:
    FrictionParameter: friction parameters
        )",
              ( py::arg( "self" ), py::arg( "friction" ) ) )
        .def( "getFrictionParameter", &ContactZone::getFrictionParameter, R"(
Get friction parameters defining method, coefficient...

Returns:
    FrictionParameter: friction parameters
        )",
              ( py::arg( "self" ) ) )
        .def( "setPairingParameter", &ContactZone::setPairingParameter, R"(
Set pairing parameters defining algorithm, distance...

Arguments:
    PairingParameter: pairing parameters
        )",
              ( py::arg( "self" ), py::arg( "pairing" ) ) )
        .def( "getPairingParameter", &ContactZone::getPairingParameter, R"(
Get pairing parameters defining algorithm, distance...

Returns:
    PairingParameter: pairing parameters
        )",
              ( py::arg( "self" ) ) )
      .def( "setSlaveGroupOfCells", &ContactZone::setSlaveGroupOfCells, R"(
Set slave's name of group of cells

Arguments:
    str: slave's name
        )",
              ( py::arg( "self" ), py::arg("slave_name") ) )
      .def( "getSlaveGroupOfCells", &ContactZone::getSlaveGroupOfCells, R"(
Get slave's name of group of cells

Returns:
    str: slave's name
        )",
              ( py::arg( "self" ) ) )
            .def( "setMasterGroupOfCells", &ContactZone::setMasterGroupOfCells, R"(
Set master's name of group of cells

Arguments:
    str: master's name
        )",
              ( py::arg( "self" ), py::arg("master_name") ) )
      .def( "getMasterGroupOfCells", &ContactZone::getMasterGroupOfCells, R"(
Get master's name of group of cells

Returns:
    str: master's name
        )",
              ( py::arg( "self" ) ) );
};
