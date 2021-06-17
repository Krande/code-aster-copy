/**
 * @file ExternalStateVariablesInterface.cxx
 * @brief Interface python de BaseExternalStateVariables
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

// Not DataStructures
// aslint: disable=C3006

#include "PythonBindings/BaseExternalStateVariablesInterface.h"
#include <PythonBindings/factory.h>
#include <boost/python.hpp>

namespace py = boost::python;

void exportBaseExternalStateVariablesToPython() {

    void ( EvolutionParameter::*c1 )( const FormulaPtr & ) =
        &EvolutionParameter::setTimeFunction;
    void ( EvolutionParameter::*c2 )( const FunctionPtr & ) =
        &EvolutionParameter::setTimeFunction;

    py::class_< EvolutionParameter, EvolutionParameterPtr >( "EvolutionParameter",
                                                                     py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< EvolutionParameter,
                                                    const TransientResultPtr & >))
        .def( "setFieldName", &EvolutionParameter::setFieldName )
        .def( "setTimeFunction", c1 )
        .def( "setTimeFunction", c2 )
        .def( "prohibitRightExtension", &EvolutionParameter::prohibitRightExtension )
        .def( "setConstantRightExtension", &EvolutionParameter::setConstantRightExtension )
        .def( "setLinearRightExtension", &EvolutionParameter::setLinearRightExtension )
        .def( "prohibitLeftExtension", &EvolutionParameter::prohibitLeftExtension )
        .def( "setConstantLeftExtension", &EvolutionParameter::setConstantLeftExtension )
        .def( "setLinearLeftExtension", &EvolutionParameter::setLinearLeftExtension );

    py::class_< BaseExternalStateVariable,
                BaseExternalStateVariable::BaseExternalStateVariablePtr>(
                    "BaseExternalStateVariable", py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< BaseExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< BaseExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >))
        .def( "hasReferenceValue", &BaseExternalStateVariable::hasReferenceValue )
        .def( "getReferenceValue", &BaseExternalStateVariable::getReferenceValue )
        .def( "setEvolutionParameter", &BaseExternalStateVariable::setEvolutionParameter )
        .def( "setValue", &BaseExternalStateVariable::setValue )
        .def( "setReferenceValue", &BaseExternalStateVariable::setReferenceValue );

    py::class_< TemperatureExternalStateVariable, TemperatureExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "TemperatureExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< TemperatureExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< TemperatureExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< GeometryExternalStateVariable, GeometryExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "GeometryExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< GeometryExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< GeometryExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< CorrosionExternalStateVariable, CorrosionExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "CorrosionExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< CorrosionExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< CorrosionExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< IrreversibleDeformationExternalStateVariable,
                IrreversibleDeformationExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >(
                                                    "IrreversibleDeformationExternalStateVariable",
                                                        py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< IrreversibleDeformationExternalStateVariable,
                                                    const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< IrreversibleDeformationExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< ConcreteHydratationExternalStateVariable,
                ConcreteHydratationExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >("ConcreteHydratationExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ConcreteHydratationExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ConcreteHydratationExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< IrradiationExternalStateVariable, IrradiationExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "IrradiationExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< IrradiationExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< IrradiationExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< SteelPhasesExternalStateVariable, SteelPhasesExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "SteelPhasesExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< SteelPhasesExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< SteelPhasesExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< ZircaloyPhasesExternalStateVariable, ZircaloyPhasesExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "ZircaloyPhasesExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ZircaloyPhasesExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ZircaloyPhasesExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< Neutral1ExternalStateVariable, Neutral1ExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "Neutral1ExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< Neutral1ExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< Neutral1ExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< Neutral2ExternalStateVariable, Neutral2ExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "Neutral2ExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< Neutral2ExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< Neutral2ExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< Neutral3ExternalStateVariable, Neutral3ExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "Neutral3ExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< Neutral3ExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< Neutral3ExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< ConcreteDryingExternalStateVariable, ConcreteDryingExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "ConcreteDryingExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ConcreteDryingExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ConcreteDryingExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_< TotalFluidPressureExternalStateVariable, TotalFluidPressureExternalStateVariablePtr,
                py::bases< BaseExternalStateVariable > >( "TotalFluidPressureExternalStateVariable",
                                                             py::no_init )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< TotalFluidPressureExternalStateVariable, const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< TotalFluidPressureExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

    py::class_<VolumetricDeformationExternalStateVariable,
                VolumetricDeformationExternalStateVariablePtr,
               py::bases< BaseExternalStateVariable > >(
                   "VolumetricDeformationExternalStateVariable", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< VolumetricDeformationExternalStateVariable,
                                                    const BaseMeshPtr & >))
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< VolumetricDeformationExternalStateVariable,
                                                    const BaseMeshPtr &, const std::string & >));

};
