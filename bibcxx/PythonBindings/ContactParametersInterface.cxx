/**
 * @file ContactParametersInterface.cxx
 * @brief Interface python de ContactParameters
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
Return the contact coefficient used. It is a value of a float

Returns:
    float: contact coefficient.
        )",
              ( py::arg( "self" ) ) )
        .def( "setCoefficient", &ContactParameter::setCoefficient, R"(
Set the contact coefficient used. It is a value of a float

Arguments:
    float: contact coefficient.
        )",
              ( py::arg( "self" ), py::arg("coeff") ) );

    py::class_< FrictionParameter, FrictionParameter::FrictionParameterPtr >( "FrictionParameter",
                                                                              py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< FrictionParameter > ) )
        .def( "getAlgorithm", &FrictionParameter::getAlgorithm, R"(
Return the Friction algorithm used. It is a value of an enum

Returns:
    FrictionAlgo: Friction algorithm.
        )",
              ( py::arg( "self" ) ) )
        .def( "setAlgorithm", &FrictionParameter::setAlgorithm, R"(
Set the Friction algorithm used. It is a value of an enum

Arguments:
    FrictionAlgo: Friction algorithm.
        )",
              ( py::arg( "self" ), py::arg( "algo" ) ) )
        .def( "getType", &FrictionParameter::getType, R"(
Return the Friction type used. It is a value of an enum

Returns:
    FrictionType: Friction type.
        )",
              ( py::arg( "self" ) ) )
        .def( "setType", &FrictionParameter::setType, R"(
Set the Friction type used. It is a value of an enum

Arguments:
    FrictionType: Friction type.
        )",
              ( py::arg( "self" ), py::arg("type") ) )
        .def( "getCoefficient", &FrictionParameter::getCoefficient, R"(
Return the Friction coefficient used. It is a value of a float

Returns:
    float: Friction coefficient.
        )",
              ( py::arg( "self" ) ) )
        .def( "setCoefficient", &FrictionParameter::setCoefficient, R"(
Set the Friction coefficient used. It is a value of a float

Arguments:
    float: Friction coefficient.
        )",
              ( py::arg( "self" ), py::arg("coeff") ) )
        .def( "getTresca", &FrictionParameter::getTresca, R"(
Return the Tresca coefficient used. It is a value of a float

Returns:
    float: Tresca coefficient.
        )",
              ( py::arg( "self" ) ) )
        .def( "setTresca", &FrictionParameter::setTresca, R"(
Set the Tresca coefficient used. It is a value of a float

Arguments:
    float: Tresca coefficient.
        )",
              ( py::arg( "self" ) , py::arg("tresca") ) )
        .def( "getCoulomb", &FrictionParameter::getCoulomb, R"(
Return the Coulomb coefficient used. It is a value of a float

Returns:
    float: Coulomb coefficient.
        )",
              ( py::arg( "self" ) ) )
        .def( "setCoulomb", &FrictionParameter::setCoulomb, R"(
Set the Coulomb coefficient used. It is a value of a float

Arguments:
    float: Coulomb coefficient.
        )",
              ( py::arg( "self" ), py::arg("coulomb") ) )
        .def( "hasFriction",
              static_cast< void ( FrictionParameter::* )( const bool & ) >( 
                    &FrictionParameter::hasFriction ), 
              R"(
Set True if friction is present in at least one contact zone else False

Arguments:
      Bool: True if friction is present else False
        )",
              ( py::arg( "self" ), py::arg( "friction" ) ) )
        .def( "hasFriction",
              static_cast< bool ( FrictionParameter::* )() const >( 
                    &FrictionParameter::hasFriction ), 
              R"(
Reruen True if friction is present in at least one contact zone else False

Returns:
      Bool: True if friction is present else False
        )",
              ( py::arg( "self" ) ) );

    py::class_< PairingParameter, PairingParameter::PairingParameterPtr >( "PairingParameter",
                                                                           py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< PairingParameter > ) )
        .def( "getAlgorithm", &PairingParameter::getAlgorithm, R"(
Return the Pairing algorithm used. It is a value of an enum

Returns:
    PairingAlgo: Pairing algorithm.
        )",
              ( py::arg( "self" ) ) )
        .def( "setAlgorithm", &PairingParameter::setAlgorithm, R"(
Set the Pairing algorithm used. It is a value of an enum

Arguments:
    PairingAlgo: Pairing algorithm.
        )",
              ( py::arg( "self" ), py::arg( "algo" ) ) )
        .def( "getPairingDistance", &PairingParameter::getPairingDistance, R"(
Return the pairing distance used. It is a value of a float

Returns:
    float: pairing distance.
        )",
              ( py::arg( "self" ) ) )
        .def( "setPairingDistance", &PairingParameter::setPairingDistance, R"(
Set the pairing distance used. It is a value of a float

Arguments:
    float: pairing distance.
        )",
              ( py::arg( "self" ), py::arg("dist_appa") ) )
        .def( "getInitialState", &PairingParameter::getInitialState, R"(
Return the initial contact state. It is a value of an enum

Returns:
    InitialState: Initial contact state.
        )",
              ( py::arg( "self" ) ) )
        .def( "setInitialState", &PairingParameter::setInitialState, R"(
Set the initial contact state. It is a value of an enum

Arguments:
    InitialState: Initial contact state.
        )",
              ( py::arg( "self" ), py::arg( "cont_init" ) ) )
        .def( "getThreshold", &PairingParameter::getThreshold, R"(
Return the distance threshold. It is a value of a float

Returns:
    float: distance threshold.
        )",
              ( py::arg( "self" ) ) )
        .def( "setThreshold", &PairingParameter::setThreshold, R"(
Set the distance threshold used. It is a value of a float

Arguments:
    float: distance threshold.
        )",
              ( py::arg( "self" ), py::arg("seuil") ) )
        .def( "getElementaryCharacteristics", &PairingParameter::getElementaryCharacteristics, R"(
Return the elementary characteristics. It is a value of a pointer

Returns:
    ElementaryCharacteristicsPtr: cara_elel pointer.
        )",
              ( py::arg( "self" ) ) )
        .def( "setElementaryCharacteristics", &PairingParameter::setElementaryCharacteristics, R"(
Set the elementary characteristics. It is a value of a pointer

Arguments:
    ElementaryCharacteristicsPtr: cara_elel pointer.
        )",
              ( py::arg( "self" ), py::arg( "cara" ) ) )
        .def( "getDistanceFunction", &PairingParameter::getDistanceFunction, R"(
Return the fictive distance function. It is a value of a pointer

Returns:
    GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        )",
              ( py::arg( "self" ) ) )
        .def( "setDistanceFunction", &PairingParameter::setDistanceFunction, R"(
Set the fictive distance function. It is a value of a pointer

Arguments:
    GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        )",
              ( py::arg( "self" ), py::arg( "dist_supp" ) ) )
        .def( "hasBeamDistance",
              static_cast< void ( PairingParameter::* )( const bool & ) >( 
                    &PairingParameter::hasBeamDistance ), 
                    R"(
Set True if fictive distance for beam is present else False

Arguments:
      Bool: True if fictive distance for beam is present else False
        )",
              ( py::arg( "self" ), py::arg( "beam" ) ) )
        .def( "hasBeamDistance",
              static_cast< bool ( PairingParameter::* )() const >( 
                    &PairingParameter::hasBeamDistance ), 
                    R"(
Reruen True if fictive distance for beam is present else False

Returns:
      Bool: True if fictive distance for beam is present else False
        )",
              ( py::arg( "self" ) ) )
        .def( "hasShellDistance",
              static_cast< void ( PairingParameter::* )( const bool & ) >( 
                    &PairingParameter::hasShellDistance ), 
                    R"(
Set True if fictive distance for shell is present else False

Arguments:
      Bool: True if fictive distance for shell is present else False
        )",
              ( py::arg( "self" ), py::arg( "shell" ) ) )
        .def( "hasShellDistance",
              static_cast< bool ( PairingParameter::* )() const >( 
                    &PairingParameter::hasShellDistance ), 
                    R"(
Reruen True if fictive distance for shell is present else False

Returns:
      Bool: True if fictive distance for shell is present else False
        )",
              ( py::arg( "self" ) ) );
};
