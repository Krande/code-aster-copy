/**
 * @file GeneralizedResultInterface.cxx
 * @brief Interface python de GeneralizedResult
 * @author Natacha Béreux
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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

#include "PythonBindings/GeneralizedResultInterface.h"

#include "aster_pybind.h"

void exportGeneralizedResultToPython( py::module_ &mod ) {

    py::class_< GeneralizedResultReal, GeneralizedResultRealPtr, DataStructure >(
        mod, "GeneralizedResultReal" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        ;

    py::class_< GeneralizedResultComplex, GeneralizedResultComplexPtr, DataStructure >(
        mod, "GeneralizedResultComplex" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        ;

    py::class_< TransientGeneralizedResult, TransientGeneralizedResultPtr, GeneralizedResultReal >(
        mod, "TransientGeneralizedResult" )
        .def( py::init( &initFactoryPtr< TransientGeneralizedResult > ) )
        .def( py::init( &initFactoryPtr< TransientGeneralizedResult, std::string > ) )
        .def( "build", &TransientGeneralizedResult::build, R"(
Builds C++ arguments associated to attributes stored by blocks of time indices
        )" )
        .def( "setGeneralizedDOFNumbering", &TransientGeneralizedResult::setGeneralizedDOFNumbering,
              R"(
Set generalized DOF numbering

Arguments:
    dofg (GeneralizedDOFNumbering): generalized DOF numbering
              )",
              py::arg( "dofg" ) )
        .def( "getGeneralizedDOFNumbering", &TransientGeneralizedResult::getGeneralizedDOFNumbering,
              R"(
Get generalized DOF numbering

Returns:
    GeneralizedDOFNumbering: generalized DOF numbering
              )" )
        .def( "setDOFNumbering", &TransientGeneralizedResult::setDOFNumbering, R"(
Set DOF numbering

Arguments:
    dofn (DOFNumbering): DOF numbering
        )",
              py::arg( "dofn" ) )
        .def( "getDOFNumbering", &TransientGeneralizedResult::getDOFNumbering, R"(
Get DOF numbering

Returns:
    DOFNumbering: DOF numbering
        )" )
        .def( "getNumberOfModes", &TransientGeneralizedResult::getNumberOfModes, R"(
Returns the number of vectors in the generalized basis

Returns:
    int: number of vectors in the generalized basis
        )" )
        .def( "getTimes", &TransientGeneralizedResult::getTimes, R"(
Returns values of instants of the transient calculation

Returns:
    list[float]: instants values
        )" )
        .def( "getIndexes", &TransientGeneralizedResult::getIndexes, R"(
Returns time indices of the transient calculation

Returns:
    list[int]: time indices
        )" )

        .def( "getDisplacementValues",
              py::overload_cast<>( &TransientGeneralizedResult::getDisplacementValues, py::const_ ),
              R"(
Return generalized displacements values for all time indices.

Returns:
    list[double]: generalized displacements values.
        )" )

        .def( "getDisplacementValues",
              py::overload_cast< ASTERINTEGER >( &TransientGeneralizedResult::getDisplacementValues,
                                                 py::const_ ),
              R"(
Return generalized displacements values at a given time index.

Arguments:
    idx (int): time index

Returns:
    list[double]: generalized displacements values.
        )",
              py::arg( "idx" ) )

        .def( "getVelocityValues",
              py::overload_cast<>( &TransientGeneralizedResult::getVelocityValues, py::const_ ),
              R"(
Return generalized velocities values for all time indices.

Returns:
    list[double]: generalized velocities values.
        )" )

        .def( "getVelocityValues",
              py::overload_cast< ASTERINTEGER >( &TransientGeneralizedResult::getVelocityValues,
                                                 py::const_ ),
              R"(
Return generalized velocities values at a given time index.

Arguments:
    idx (int): time index

Returns:
    list[double]: generalized velocities values.
        )",
              py::arg( "idx" ) )

        .def( "getAccelerationValues",
              py::overload_cast<>( &TransientGeneralizedResult::getAccelerationValues, py::const_ ),
              R"(
Return generalized accelerations values for all time indices.

Returns:
    list[double]: generalized accelerations values.
        )" )

        .def( "getAccelerationValues",
              py::overload_cast< ASTERINTEGER >( &TransientGeneralizedResult::getAccelerationValues,
                                                 py::const_ ),
              R"(
Return generalized accelerations values at a given time index.

Arguments:
    idx (int): time index

Returns:
    list[double]: generalized accelerations values.
        )",
              py::arg( "idx" ) )

        .def( "setDisplacementValues",
              py::overload_cast< ASTERINTEGER, VectorReal >(
                  &TransientGeneralizedResult::setDisplacementValues, py::const_ ),
              R"(
Set generalized displacement values at a given time index.

Arguments:
    idx (int): time index

    val (list[double]): generalized displacement values.
        )",
              py::arg( "idx" ), py::arg( "val" ) )

        .def( "setVelocityValues",
              py::overload_cast< ASTERINTEGER, VectorReal >(
                  &TransientGeneralizedResult::setVelocityValues, py::const_ ),
              R"(
Set generalized velocity values at a given time index.

Arguments:
    idx (int): time index

    val (list[double]): generalized velocity values.
        )",
              py::arg( "idx" ), py::arg( "val" ) )

        .def( "setAccelerationValues",
              py::overload_cast< ASTERINTEGER, VectorReal >(
                  &TransientGeneralizedResult::setAccelerationValues, py::const_ ),
              R"(
Set generalized acceleration values at a given time index.

Arguments:
    idx (int): time index

    val (list[double]): generalized acceleration values.
        )",
              py::arg( "idx" ), py::arg( "val" ) )

        ;

    py::class_< HarmoGeneralizedResult, HarmoGeneralizedResultPtr, GeneralizedResultComplex >(
        mod, "HarmoGeneralizedResult" )
        .def( py::init( &initFactoryPtr< HarmoGeneralizedResult > ) )
        .def( py::init( &initFactoryPtr< HarmoGeneralizedResult, std::string > ) )
        .def( "getGeneralizedDOFNumbering", &HarmoGeneralizedResult::getGeneralizedDOFNumbering )
        .def( "setGeneralizedDOFNumbering", &HarmoGeneralizedResult::setGeneralizedDOFNumbering )
        .def( "setDOFNumbering", &HarmoGeneralizedResult::setDOFNumbering )
        .def( "getDOFNumbering", &HarmoGeneralizedResult::getDOFNumbering )
        .def( "setDisplacement", &HarmoGeneralizedResult::setDisplacement )
        .def( "getDisplacement", &HarmoGeneralizedResult::getDisplacement )
        .def( "setVelocity", &HarmoGeneralizedResult::setVelocity )
        .def( "setAcceleration", &HarmoGeneralizedResult::setAcceleration );
};
