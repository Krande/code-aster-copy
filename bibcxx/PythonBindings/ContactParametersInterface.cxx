/**
 * @file ContactParametersInterface.cxx
 * @brief Interface python de ContactParameters
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
#include "PythonBindings/ContactParametersInterface.h"
#include <PythonBindings/factory.h>

// aslint: disable=C3006

void exportContactParametersToPython() {

    py::class_< ContactParameter, ContactParameter::ContactParameterPtr >( "ContactParameter",
                                                                           py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< ContactParameter > ) )
        .def( "getAlgorithm", &ContactParameter::getAlgorithm, R"(
Return the contact algorithm used. It is a value of an enum

Returns:
    ContactAlgo: contact algorithm.
        )",
              ( py::arg( "self" ) ) )
        .def( "setAlgorithm", &ContactParameter::setAlgorithm, R"(
Set the contact algorithm used. It is a value of an enum

Arguments:
    ContactAlgo: contact algorithm.
        )",
              ( py::arg( "self" ), py::arg( "algo" ) ) )
        .def( "getType", &ContactParameter::getType, R"(
Return the contact type used. It is a value of an enum

Returns:
    ContactType: contact type.
        )",
              ( py::arg( "self" ) ) )
      .def( "setType", &ContactParameter::setType, R"(
Set the contact type used. It is a value of an enum

Arguments:
    ContactType: contact type.
        )",
              ( py::arg( "self" ), py::arg("type") ) )
        .def( "getVariant", &ContactParameter::getVariant, R"(
Return the contact variant used. It is a value of an enum

Returns:
    ContactVariant: contact variant.
        )",
              ( py::arg( "self" ) ) )
      .def( "setVariant", &ContactParameter::setVariant, R"(
Set the contact variant used. It is a value of an enum

Arguments:
    ContactVariant: contact variant.
        )",
              ( py::arg( "self" ), py::arg("variant") ) )
        .def( "getCoefficient", &ContactParameter::getCoefficient, R"(
Return the contact coefficient used. It is a value of an enum

Returns:
    float: contact coefficient.
        )",
              ( py::arg( "self" ) ) )
      .def( "setCoefficient", &ContactParameter::setCoefficient, R"(
Set the contact coefficient used. It is a value of an enum

Arguments:
    float: contact coefficient.
        )",
              ( py::arg( "self" ) ) );

    py::class_< FrictionParameter, FrictionParameter::FrictionParameterPtr >( "FrictionParameter",
                                                                              py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< FrictionParameter > ) );

    py::class_< PairingParameter, PairingParameter::PairingParameterPtr >( "PairingParameter",
                                                                           py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< PairingParameter > ) );
};
