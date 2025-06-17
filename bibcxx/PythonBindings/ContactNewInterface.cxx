/**
 * @file ContactNewInterface.cxx
 * @brief Interface python de ContactNew
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

#include "PythonBindings/ContactNewInterface.h"

#include "aster_pybind.h"

void exportContactNewToPython( py::module_ &mod ) {

    py::class_< ContactNew, ContactNewPtr, DSWithCppPickling >( mod, "ContactNew" )
        .def( py::init( &initFactoryPtr< ContactNew, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ContactNew, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ContactNew, const py::tuple & > ) )
        .def( define_pickling< ContactNew >() )
        .def( "getModel", &ContactNew::getModel, R"(
Return the model used in the contact definition

Returns:
    Model: model.
        )" )
        .def( "getMesh", &ContactNew::getMesh, R"(
Return the mesh used in the contact definition

Returns:
    Mesh: mesh.
        )" )
        .def( "getFiniteElementDescriptor", &ContactNew::getFiniteElementDescriptor, R"(
Return the finite element descriptor to define virtual cells for Lagrange multipliers

Returns:
    FiniteElementDescriptor: fed.
        )" )
        .def( "getNumberOfContactZones", &ContactNew::getNumberOfContactZones, R"(
Return the number of contact zones used

Returns:
    inter: number of contact zones.
        )" )
        .def( "getContactZone", &ContactNew::getContactZone, R"(
Return the specified contact zone

Arguments:
    zone_id (int): index of the contact zone (0-based)

Returns:
    ContactZone: contact zone.
        )",
              py::arg( "zone_id" ) )
        .def( "getContactZones", &ContactNew::getContactZones, R"(
Return the list of contact zones

Returns:
    List[ContactZone]: List of contact zones.
        )" )
        .def( "appendContactZone", &ContactNew::appendContactZone, R"(
Append a new contact zone to the contact definition

Arguments:
    zone (ContactZone): contact zone to append
        )",
              py::arg( "zone" ) )
        .def( "setVerbosity", &ContactNew::setVerbosity, R"(
Set level of verbosity:
      0- without
      1- normal (default)
      2- detailled

Arguments:
    level (int): level of verbosity
        )",
              py::arg( "level" ) )
        .def( "getVerbosity", &ContactNew::getVerbosity, R"(
Get level of verbosity:*
      0- without
      1- normal
      2- detailled

Returns:
    integer: level of verbosity
        )" )
        .def( "build", &ContactNew::build, R"(
Build and check internal objects
        )" )
        .def_property( "hasFriction", &ContactNew::hasFriction, &ContactNew::enableFriction, R"(
bool: enable or disable the use of friction.
        )" )
        .def_property( "hasSmoothing", &ContactNew::hasSmoothing, &ContactNew::enableSmoothing, R"(
bool: enable or disable  the use of smoothing.
        )" )
        .def( "isParallel", &ContactNew::isParallel, R"(
bool: true if parallel contact.
        )" );

    py::class_< FrictionNew, FrictionNewPtr, ContactNew >( mod, "FrictionNew" )
        .def( py::init( &initFactoryPtr< FrictionNew, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FrictionNew, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< FrictionNew, const py::tuple & > ) )
        .def( define_pickling< FrictionNew >() );
};
