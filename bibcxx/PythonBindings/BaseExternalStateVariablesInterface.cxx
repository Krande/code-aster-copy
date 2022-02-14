/**
 * @file ExternalStateVariablesInterface.cxx
 * @brief Interface python de BaseExternalStateVariables
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

// Not DataStructures
// aslint: disable=C3006

#include "PythonBindings/BaseExternalStateVariablesInterface.h"

#include "aster_pybind.h"

void exportBaseExternalStateVariablesToPython( py::module_ &mod ) {

    py::class_< EvolutionParameter, EvolutionParameterPtr >( mod, "EvolutionParameter" )
        .def( py::init( &initFactoryPtr< EvolutionParameter, const TransientResultPtr & > ) )
        .def( "setFieldName", &EvolutionParameter::setFieldName )
        .def( "setTimeFunction",
              py::overload_cast< const FormulaPtr & >( &EvolutionParameter::setTimeFunction ) )
        .def( "setTimeFunction",
              py::overload_cast< const FunctionPtr & >( &EvolutionParameter::setTimeFunction ) )
        .def( "prohibitRightExtension", &EvolutionParameter::prohibitRightExtension )
        .def( "setConstantRightExtension", &EvolutionParameter::setConstantRightExtension )
        .def( "setLinearRightExtension", &EvolutionParameter::setLinearRightExtension )
        .def( "prohibitLeftExtension", &EvolutionParameter::prohibitLeftExtension )
        .def( "setConstantLeftExtension", &EvolutionParameter::setConstantLeftExtension )
        .def( "setLinearLeftExtension", &EvolutionParameter::setLinearLeftExtension );

    py::class_< BaseExternalStateVariable,
                BaseExternalStateVariable::BaseExternalStateVariablePtr >(
        mod, "BaseExternalStateVariable" )
        .def( py::init( &initFactoryPtr< BaseExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< BaseExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) )
        .def( "hasReferenceValue", &BaseExternalStateVariable::hasReferenceValue )
        .def( "getReferenceValue", &BaseExternalStateVariable::getReferenceValue )
        .def( "setEvolutionParameter", &BaseExternalStateVariable::setEvolutionParameter )
        .def( "setValue", &BaseExternalStateVariable::setValue )
        .def( "setReferenceValue", &BaseExternalStateVariable::setReferenceValue );

    py::class_< TemperatureExternalStateVariable, TemperatureExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "TemperatureExternalStateVariable" )
        .def( py::init( &initFactoryPtr< TemperatureExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< TemperatureExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< GeometryExternalStateVariable, GeometryExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "GeometryExternalStateVariable" )
        .def( py::init( &initFactoryPtr< GeometryExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< GeometryExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< CorrosionExternalStateVariable, CorrosionExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "CorrosionExternalStateVariable" )
        .def( py::init( &initFactoryPtr< CorrosionExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< CorrosionExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< IrreversibleDeformationExternalStateVariable,
                IrreversibleDeformationExternalStateVariablePtr, BaseExternalStateVariable >(
        mod, "IrreversibleDeformationExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< IrreversibleDeformationExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< IrreversibleDeformationExternalStateVariable,
                                         const BaseMeshPtr &, const std::string & > ) );

    py::class_< ConcreteHydratationExternalStateVariable,
                ConcreteHydratationExternalStateVariablePtr, BaseExternalStateVariable >(
        mod, "ConcreteHydratationExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< ConcreteHydratationExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< ConcreteHydratationExternalStateVariable,
                                         const BaseMeshPtr &, const std::string & > ) );

    py::class_< IrradiationExternalStateVariable, IrradiationExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "IrradiationExternalStateVariable" )
        .def( py::init( &initFactoryPtr< IrradiationExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< IrradiationExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< SteelPhasesExternalStateVariable, SteelPhasesExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "SteelPhasesExternalStateVariable" )
        .def( py::init( &initFactoryPtr< SteelPhasesExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< SteelPhasesExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< ZircaloyPhasesExternalStateVariable, ZircaloyPhasesExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "ZircaloyPhasesExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< ZircaloyPhasesExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< ZircaloyPhasesExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< Neutral1ExternalStateVariable, Neutral1ExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "Neutral1ExternalStateVariable" )
        .def( py::init( &initFactoryPtr< Neutral1ExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< Neutral1ExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< Neutral2ExternalStateVariable, Neutral2ExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "Neutral2ExternalStateVariable" )
        .def( py::init( &initFactoryPtr< Neutral2ExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< Neutral2ExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< Neutral3ExternalStateVariable, Neutral3ExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "Neutral3ExternalStateVariable" )
        .def( py::init( &initFactoryPtr< Neutral3ExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< Neutral3ExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< ConcreteDryingExternalStateVariable, ConcreteDryingExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "ConcreteDryingExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< ConcreteDryingExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< ConcreteDryingExternalStateVariable, const BaseMeshPtr &,
                                         const std::string & > ) );

    py::class_< TotalFluidPressureExternalStateVariable, TotalFluidPressureExternalStateVariablePtr,
                BaseExternalStateVariable >( mod, "TotalFluidPressureExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< TotalFluidPressureExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< TotalFluidPressureExternalStateVariable,
                                         const BaseMeshPtr &, const std::string & > ) );

    py::class_< VolumetricDeformationExternalStateVariable,
                VolumetricDeformationExternalStateVariablePtr, BaseExternalStateVariable >(
        mod, "VolumetricDeformationExternalStateVariable" )
        .def( py::init(
            &initFactoryPtr< VolumetricDeformationExternalStateVariable, const BaseMeshPtr & > ) )
        .def( py::init( &initFactoryPtr< VolumetricDeformationExternalStateVariable,
                                         const BaseMeshPtr &, const std::string & > ) );
};
