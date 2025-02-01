/**
 * @file ContactParameterInterface.cxx
 * @brief Interface python de ContactParameter
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

#include "PythonBindings/ContactParameterInterface.h"

#include "aster_pybind.h"

// aslint: disable=C3006

void exportContactParameterToPython( py::module_ &mod ) {

    py::class_< ContactParameter, ContactParameter::ContactParameterPtr >( mod, "ContactParameter" )
        .def( py::init( &initFactoryPtr< ContactParameter > ) )
        .def( define_pickling< ContactParameter >() )
        .def( "getAlgorithm", &ContactParameter::getAlgorithm, R"(
Return the contact algorithm used. It is a value of an enum

Returns:
    ContactAlgo: contact algorithm.
        )" )
        .def( "setAlgorithm", &ContactParameter::setAlgorithm, R"(
Set the contact algorithm used. It is a value of an enum

Arguments:
    ContactAlgo: contact algorithm.
        )",
              py::arg( "algo" ) )
        .def( "getType", &ContactParameter::getType, R"(
Return the contact type used. It is a value of an enum

Returns:
    ContactType: contact type.
        )" )
        .def( "setType", &ContactParameter::setType, R"(
Set the contact type used. It is a value of an enum

Arguments:
    ContactType: contact type.
        )",
              py::arg( "type" ) )
        .def( "getJacobianType", &ContactParameter::getJacobianType, R"(
Return how the Jacobian is computed. It is a value of an enum

Returns:
    JacobianType: Jacobian type.
        )" )
        .def( "setJacobianType", &ContactParameter::setJacobianType, R"(
Set how the Jacobian is computed. It is a value of an enum

Arguments:
    JacobianType: Jacobian type.
        )",
              py::arg( "type" ) )
        .def( "getVariant", &ContactParameter::getVariant, R"(
Return the contact variant used. It is a value of an enum

Returns:
    ContactVariant: contact variant.
        )" )
        .def( "setVariant", &ContactParameter::setVariant, R"(
Set the contact variant used. It is a value of an enum

Arguments:
    ContactVariant: contact variant.
        )",
              py::arg( "variant" ) )
        .def( "getCoefficient", &ContactParameter::getCoefficient, R"(
Return the contact coefficient used. It is a value of a float

Returns:
    float: contact coefficient.
        )" )
        .def( "setCoefficient", &ContactParameter::setCoefficient, R"(
Set the contact coefficient used. It is a value of a float

Arguments:
    float: contact coefficient.
        )",
              py::arg( "coeff" ) );

    py::class_< FrictionParameter, FrictionParameter::FrictionParameterPtr >( mod,
                                                                              "FrictionParameter" )
        .def( py::init( &initFactoryPtr< FrictionParameter > ) )
        .def( define_pickling< FrictionParameter >() )
        .def( "getAlgorithm", &FrictionParameter::getAlgorithm, R"(
Return the Friction algorithm used. It is a value of an enum

Returns:
    FrictionAlgo: Friction algorithm.
        )" )
        .def( "setAlgorithm", &FrictionParameter::setAlgorithm, R"(
Set the Friction algorithm used. It is a value of an enum

Arguments:
    FrictionAlgo: Friction algorithm.
        )",
              py::arg( "algo" ) )
        .def( "getType", &FrictionParameter::getType, R"(
Return the Friction type used. It is a value of an enum

Returns:
    FrictionType: Friction type.
        )" )
        .def( "setType", &FrictionParameter::setType, R"(
Set the Friction type used. It is a value of an enum

Arguments:
    FrictionType: Friction type.
        )",
              py::arg( "type" ) )
        .def( "getCoefficient", &FrictionParameter::getCoefficient, R"(
Return the Friction coefficient used. It is a value of a float

Returns:
    float: Friction coefficient.
        )" )
        .def( "setCoefficient", &FrictionParameter::setCoefficient, R"(
Set the Friction coefficient used. It is a value of a float

Arguments:
    float: Friction coefficient.
        )",
              py::arg( "coeff" ) )
        .def( "getTresca", &FrictionParameter::getTresca, R"(
Return the Tresca coefficient used. It is a value of a float

Returns:
    float: Tresca coefficient.
        )" )
        .def( "setTresca", &FrictionParameter::setTresca, R"(
Set the Tresca coefficient used. It is a value of a float

Arguments:
    float: Tresca coefficient.
        )",
              py::arg( "tresca" ) )
        .def( "getCoulomb", &FrictionParameter::getCoulomb, R"(
Return the Coulomb coefficient used. It is a value of a float

Returns:
    float: Coulomb coefficient.
        )" )
        .def( "setCoulomb", &FrictionParameter::setCoulomb, R"(
Set the Coulomb coefficient used. It is a value of a float

Arguments:
    float: Coulomb coefficient.
        )",
              py::arg( "coulomb" ) )
        .def_property( "hasFriction", &FrictionParameter::hasFriction,
                       &FrictionParameter::enableFriction, R"(
bool: enable or disable the use of friction.
        )" );

    py::class_< PairingParameter, PairingParameter::PairingParameterPtr >( mod, "PairingParameter" )
        .def( py::init( &initFactoryPtr< PairingParameter > ) )
        .def( define_pickling< PairingParameter >() )
        .def( "getAlgorithm", &PairingParameter::getAlgorithm, R"(
Return the Pairing algorithm used. It is a value of an enum

Returns:
    PairingAlgo: Pairing algorithm.
        )" )
        .def( "setAlgorithm", &PairingParameter::setAlgorithm, R"(
Set the Pairing algorithm used. It is a value of an enum

Arguments:
    PairingAlgo: Pairing algorithm.
        )",
              py::arg( "algo" ) )
        .def( "getDistanceRatio", &PairingParameter::getDistanceRatio, R"(
Return the pairing distance ratio used. It is a value of a float

Returns:
    float: pairing distance.
        )" )
        .def( "setDistanceRatio", &PairingParameter::setDistanceRatio, R"(
Set the pairing distance ratio used. It is a value of a float

Arguments:
    float: pairing distance ratio.
        )",
              py::arg( "dist_ratio" ) )
        .def( "getInitialState", &PairingParameter::getInitialState, R"(
Return the initial contact state. It is a value of an enum

Returns:
    InitialState: Initial contact state.
        )" )
        .def( "setInitialState", &PairingParameter::setInitialState, R"(
Set the initial contact state. It is a value of an enum

Arguments:
    InitialState: Initial contact state.
        )",
              py::arg( "cont_init" ) )
        .def( "getElementaryCharacteristics", &PairingParameter::getElementaryCharacteristics, R"(
Return the elementary characteristics. It is a value of a pointer

Returns:
    ElementaryCharacteristicsPtr: cara_elel pointer.
        )" )
        .def( "setElementaryCharacteristics", &PairingParameter::setElementaryCharacteristics, R"(
Set the elementary characteristics. It is a value of a pointer

Arguments:
    ElementaryCharacteristicsPtr: cara_elel pointer.
        )",
              py::arg( "cara" ) )
        .def( "getDistanceFunction", &PairingParameter::getDistanceFunction, R"(
Return the fictive distance function. It is a value of a pointer

Returns:
    GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        )" )
        .def( "setDistanceFunction", &PairingParameter::setDistanceFunction, R"(
Set the fictive distance function. It is a value of a pointer

Arguments:
    GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        )",
              py::arg( "dist_supp" ) )
        .def_property( "hasBeamDistance", &PairingParameter::hasBeamDistance,
                       &PairingParameter::enableBeamDistance, R"(
bool: enable or disable the use of a fictive distance for beam.
        )" )
        .def_property( "hasShellDistance", &PairingParameter::hasShellDistance,
                       &PairingParameter::enableShellDistance, R"(
bool: enable or disable the use of a fictive distance for shell.
        )" );
};
