/**
 * @file ContactNewInterface.cxx
 * @brief Interface python de ContactNew
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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/ContactNewInterface.h"
#include <PythonBindings/factory.h>

void exportContactNewToPython() {

    py::class_< ContactNew, ContactNew::ContactNewPtr, py::bases< DataStructure > >( "ContactNew",
                                                                                     py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ContactNew, std::string, ModelPtr > ) )
        .def( "__init__", py::make_constructor( &initFactoryPtr< ContactNew, ModelPtr > ) )
        .def( "getModel", &ContactNew::getModel, R"(
Return the model used in the contact definition

Returns:
    Model: model.
        )",
              ( py::arg( "self" ) ) )
        .def( "getMesh", &ContactNew::getMesh, R"(
Return the mesh used in the contact definition

Returns:
    BaseMesh: mesh.
        )",
              ( py::arg( "self" ) ) )
        .def( "getFiniteElementDescriptor", &ContactNew::getFiniteElementDescriptor, R"(
Return the finite element descriptor to define virtual cells for Lagrange multipliers

Returns:
    FiniteElementDescriptor: fed.
        )",
              ( py::arg( "self" ) ) )
        .def( "getNumberOfContactZones", &ContactNew::getNumberOfContactZones, R"(
Return the number of contact zones used

Returns:
    inter: number of contact zones.
        )",
              ( py::arg( "self" ) ) )
        .def( "getContactZone", &ContactNew::getContactZone, R"(
Return the specified contact zone

Arguments:
    int: index of the contact zone (0-based)

Returns:
    ContactZone: contact zone.
        )",
              ( py::arg( "self" ), py::arg( "zone_id" ) ) )
        .def( "appendContactZone", &ContactNew::appendContactZone, R"(
Append a new contact zone to the contact definition

Arguments:
    ContactZone: contact zone to append
        )",
              ( py::arg( "self" ), py::arg( "contact_zone" ) ) )
        .def( "setVerbosity", &ContactNew::setVerbosity, R"(
Set level of verbosity:
      0- without
      1- normal (default)
      2- detailled

Arguments:
    integer: level of verbosity
        )",
              ( py::arg( "self" ), py::arg( "level" ) ) )
        .def( "getVerbosity", &ContactNew::getVerbosity, R"(
Get level of verbosity:*
      0- without
      1- normal
      2- detailled

Returns:
    integer: level of verbosity
        )",
              ( py::arg( "self" ) ) )
        .def( "build", &ContactNew::build, R"(
Build and check internal objects
        )",
              ( py::arg( "self" ) ) )
        .def( "hasFriction",
              static_cast< void ( ContactNew::* )( const bool & ) >( &ContactNew::hasFriction ), R"(
Set True if friction is present in at least one contact zone else False

Arguments:
      Bool: True if friction is present else False
        )",
              ( py::arg( "self" ), py::arg( "friction" ) ) )
        .def( "hasFriction",
              static_cast< bool ( ContactNew::* )() const >( &ContactNew::hasFriction ), R"(
Reruen True if friction is present in at least one contact zone else False

Returns:
      Bool: True if friction is present else False
        )",
              ( py::arg( "self" ) ) )
        .def( "hasSmoothing",
              static_cast< void ( ContactNew::* )( const bool & ) >( &ContactNew::hasSmoothing ),
              R"(
Set True if smoothing is used to compute outward normals else False

Arguments:
      Bool: True if smoothing is used else False
        )",
              ( py::arg( "self" ), py::arg( "smoothing" ) ) )
        .def( "hasSmoothing",
              static_cast< bool ( ContactNew::* )() const >( &ContactNew::hasSmoothing ), R"(
Reruen True if smoothing is used to compute outward normals else False

Returns:
      Bool: True if smoothing is used else False
        )",
              ( py::arg( "self" ) ) );
};
