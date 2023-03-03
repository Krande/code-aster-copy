/**
 * @file ModeResultInterface.cxx
 * @brief Interface python de ModeResult
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "PythonBindings/ModeResultInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/VariantStiffnessMatrixInterface.h"

void exportModeResultToPython( py::module_ &mod ) {

    py::class_< ModeResult, ModeResultPtr, FullResult >( mod, "ModeResult" )
        .def( py::init( &initFactoryPtr< ModeResult > ) )
        .def( py::init( &initFactoryPtr< ModeResult, std::string > ) )
        .def( "getDOFNumbering", &ModeResult::getDOFNumbering )
        .def( "getStiffnessMatrix", &getStiffnessMatrix< ModeResultPtr > )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixDisplacementRealPtr & >(
                                        &ModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixTemperatureRealPtr & >(
                                        &ModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix",
              py::overload_cast< const AssemblyMatrixDisplacementComplexPtr & >(
                  &ModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixPressureRealPtr & >(
                                        &ModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixPressureRealPtr & >(
                                        &ModeResult::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const GeneralizedAssemblyMatrixRealPtr & >(
                                        &ModeResult::setStiffnessMatrix ) )
        .def( "getMassMatrix", &getStiffnessMatrix< ModeResultPtr > )
        .def( "setMassMatrix", py::overload_cast< const AssemblyMatrixDisplacementRealPtr & >(
                                   &ModeResult::setMassMatrix ) )
        .def( "setMassMatrix", py::overload_cast< const AssemblyMatrixTemperatureRealPtr & >(
                                   &ModeResult::setMassMatrix ) )
        .def( "setMassMatrix", py::overload_cast< const AssemblyMatrixDisplacementComplexPtr & >(
                                   &ModeResult::setMassMatrix ) )
        .def( "setMassMatrix", py::overload_cast< const GeneralizedAssemblyMatrixComplexPtr & >(
                                   &ModeResult::setMassMatrix ) )
        .def( "setStructureInterface", &ModeResult::setStructureInterface )
        .def( "getNumberOfDynamicModes", &ModeResult::getNumberOfDynamicModes )
        .def( "getNumberOfStaticModes", &ModeResult::getNumberOfStaticModes );
};

void exportModeResultComplexToPython( py::module_ &mod ) {

    py::class_< ModeResultComplex, ModeResultComplexPtr, ModeResult >( mod, "ModeResultComplex" )
        .def( py::init( &initFactoryPtr< ModeResultComplex > ) )
        .def( py::init( &initFactoryPtr< ModeResultComplex, std::string > ) )
        .def( "setDampingMatrix", &ModeResultComplex::setDampingMatrix )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixDisplacementRealPtr & >(
                                        &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix",
              py::overload_cast< const AssemblyMatrixDisplacementComplexPtr & >(
                  &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixTemperatureRealPtr & >(
                                        &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const AssemblyMatrixPressureRealPtr & >(
                                        &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix", py::overload_cast< const GeneralizedAssemblyMatrixRealPtr & >(
                                        &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStiffnessMatrix",
              py::overload_cast< const GeneralizedAssemblyMatrixComplexPtr & >(
                  &ModeResultComplex::setStiffnessMatrix ) )
        .def( "setStructureInterface", &ModeResultComplex::setStructureInterface );
};
