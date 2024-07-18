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
        .def( "build", &TransientGeneralizedResult::build )
        .def( "setGeneralizedDOFNumbering",
              &TransientGeneralizedResult::setGeneralizedDOFNumbering )
        .def( "getGeneralizedDOFNumbering",
              &TransientGeneralizedResult::getGeneralizedDOFNumbering )
        .def( "setDOFNumbering", &TransientGeneralizedResult::setDOFNumbering )
        .def( "getDOFNumbering", &TransientGeneralizedResult::getDOFNumbering )
        .def( "getNumberOfModes", &TransientGeneralizedResult::getNumberOfModes )
        .def( "getTimes", &TransientGeneralizedResult::getTimes )
        .def( "getIndexes", &TransientGeneralizedResult::getIndexes )
        .def( "getDisplacementValues", &TransientGeneralizedResult::getDisplacementValues )
        .def( "getVelocityValues", &TransientGeneralizedResult::getVelocityValues )
        .def( "getAccelerationValues", &TransientGeneralizedResult::getAccelerationValues )
        .def( "getDisplacementValuesAtIndex",
              &TransientGeneralizedResult::getDisplacementValuesAtIndex )
        .def( "getVelocityValuesAtIndex", &TransientGeneralizedResult::getVelocityValuesAtIndex )
        .def( "getAccelerationValuesAtIndex",
              &TransientGeneralizedResult::getAccelerationValuesAtIndex );

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
