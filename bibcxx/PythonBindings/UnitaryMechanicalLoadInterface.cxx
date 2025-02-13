/**
 * @file UnitaryMechanicalLoadInterface.cxx
 * @brief Interface python de MechanicalLoad
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

#include "PythonBindings/UnitaryMechanicalLoadInterface.h"

#include "aster_pybind.h"

void exportUnitaryMechanicalLoadToPython( py::module_ &mod ) {

    py::enum_< LoadEnum >( mod, "Loads", R"(
Enumeration for type of load.
    )" )
        .value( "NodalForce", NodalForce )
        .value( "ForceOnEdge", ForceOnEdge )
        .value( "ForceOnFace", ForceOnFace )
        .value( "LineicForce", LineicForce )
        .value( "InternalForce", InternalForce )
        .value( "ForceOnBeam", ForceOnBeam )
        .value( "ForceOnShell", ForceOnShell )
        .value( "PressureOnPipe", PressureOnPipe )
        .value( "ImposedDoF", ImposedDoF )
        .value( "DistributedPressure", DistributedPressure )
        .value( "NormalSpeedOnFace", NormalSpeedOnFace )
        .value( "WavePressureOnFace", WavePressureOnFace )
        .value( "THMFlux", THMFlux )
        .export_values();

    py::class_< NodalForceReal, NodalForceReal::UnitaryMechanicalLoadRealPtr, MechanicalLoadReal >(
        mod, "NodalForceReal" )
        .def( py::init( &initFactoryPtr< NodalForceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< NodalForceReal, std::string, ModelPtr > ) )
        .def( "build", &NodalForceReal::build )
        .def( "setValue", &NodalForceReal::setValue );

    py::class_< NodalStructuralForceReal, NodalStructuralForceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "NodalStructuralForceReal" )
        .def( py::init( &initFactoryPtr< NodalStructuralForceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< NodalStructuralForceReal, std::string, ModelPtr > ) )
        .def( "build", &NodalStructuralForceReal::build )
        .def( "setValue", &NodalStructuralForceReal::setValue );

    py::class_< ForceOnFaceReal, ForceOnFaceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "ForceOnFaceReal" )
        .def( py::init( &initFactoryPtr< ForceOnFaceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ForceOnFaceReal, std::string, ModelPtr > ) )
        .def( "build", &ForceOnFaceReal::build )
        .def( "setValue", &ForceOnFaceReal::setValue );

    py::class_< ForceOnEdgeReal, ForceOnEdgeReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "ForceOnEdgeReal" )
        .def( py::init( &initFactoryPtr< ForceOnEdgeReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ForceOnEdgeReal, std::string, ModelPtr > ) )
        .def( "build", &ForceOnEdgeReal::build )
        .def( "setValue", &ForceOnEdgeReal::setValue );

    py::class_< StructuralForceOnEdgeReal, StructuralForceOnEdgeReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "StructuralForceOnEdgeReal" )
        .def( py::init( &initFactoryPtr< StructuralForceOnEdgeReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< StructuralForceOnEdgeReal, std::string, ModelPtr > ) )
        .def( "build", &StructuralForceOnEdgeReal::build )
        .def( "setValue", &StructuralForceOnEdgeReal::setValue );

    py::class_< LineicForceReal, LineicForceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "LineicForceReal" )
        .def( py::init( &initFactoryPtr< LineicForceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< LineicForceReal, std::string, ModelPtr > ) )
        .def( "build", &LineicForceReal::build )
        .def( "setValue", &LineicForceReal::setValue );

    py::class_< InternalForceReal, InternalForceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "InternalForceReal" )
        .def( py::init( &initFactoryPtr< InternalForceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< InternalForceReal, std::string, ModelPtr > ) )
        .def( "build", &InternalForceReal::build )
        .def( "setValue", &InternalForceReal::setValue );

    py::class_< StructuralForceOnBeamReal, StructuralForceOnBeamReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "StructuralForceOnBeamReal" )
        .def( py::init( &initFactoryPtr< StructuralForceOnBeamReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< StructuralForceOnBeamReal, std::string, ModelPtr > ) )
        .def( "build", &StructuralForceOnBeamReal::build )
        .def( "setValue", &StructuralForceOnBeamReal::setValue );

    py::class_< LocalForceOnBeamReal, LocalForceOnBeamReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "LocalForceOnBeamReal" )
        .def( py::init( &initFactoryPtr< LocalForceOnBeamReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< LocalForceOnBeamReal, std::string, ModelPtr > ) )
        .def( "build", &LocalForceOnBeamReal::build )
        .def( "setValue", &LocalForceOnBeamReal::setValue );

    py::class_< StructuralForceOnShellReal,
                StructuralForceOnShellReal::UnitaryMechanicalLoadRealPtr, MechanicalLoadReal >(
        mod, "StructuralForceOnShellReal" )
        .def( py::init( &initFactoryPtr< StructuralForceOnShellReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< StructuralForceOnShellReal, std::string, ModelPtr > ) )
        .def( "build", &StructuralForceOnShellReal::build )
        .def( "setValue", &StructuralForceOnShellReal::setValue );

    py::class_< LocalForceOnShellReal, LocalForceOnShellReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "LocalForceOnShellReal" )
        .def( py::init( &initFactoryPtr< LocalForceOnShellReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< LocalForceOnShellReal, std::string, ModelPtr > ) )
        .def( "build", &LocalForceOnShellReal::build )
        .def( "setValue", &LocalForceOnShellReal::setValue );

    py::class_< PressureOnShellReal, PressureOnShellReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "PressureOnShellReal" )
        .def( py::init( &initFactoryPtr< PressureOnShellReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< PressureOnShellReal, std::string, ModelPtr > ) )
        .def( "build", &PressureOnShellReal::build )
        .def( "setValue", &PressureOnShellReal::setValue );

    py::class_< PressureOnPipeReal, PressureOnPipeReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "PressureOnPipeReal" )
        .def( py::init( &initFactoryPtr< PressureOnPipeReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< PressureOnPipeReal, std::string, ModelPtr > ) )
        .def( "build", &PressureOnPipeReal::build )
        .def( "setValue", &PressureOnPipeReal::setValue );

    py::class_< ImposedDisplacementReal, ImposedDisplacementReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "ImposedDisplacementReal" )
        .def( py::init( &initFactoryPtr< ImposedDisplacementReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ImposedDisplacementReal, std::string, ModelPtr > ) )
        .def( "build", &ImposedDisplacementReal::build )
        .def( "setValue", &ImposedDisplacementReal::setValue );

    py::class_< ImposedPressureReal, ImposedPressureReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "ImposedPressureReal" )
        .def( py::init( &initFactoryPtr< ImposedPressureReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< ImposedPressureReal, std::string, ModelPtr > ) )
        .def( "build", &ImposedPressureReal::build )
        .def( "setValue", &ImposedPressureReal::setValue );

    py::class_< DistributedPressureReal, DistributedPressureReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "DistributedPressureReal" )
        .def( py::init( &initFactoryPtr< DistributedPressureReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< DistributedPressureReal, std::string, ModelPtr > ) )
        .def( "build", &DistributedPressureReal::build )
        .def( "setValue", &DistributedPressureReal::setValue );

    py::class_< NormalSpeedOnFaceReal, NormalSpeedOnFaceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "NormalSpeedOnFaceReal" )
        .def( py::init( &initFactoryPtr< NormalSpeedOnFaceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< NormalSpeedOnFaceReal, std::string, ModelPtr > ) )
        .def( "build", &NormalSpeedOnFaceReal::build )
        .def( "setValue", &NormalSpeedOnFaceReal::setValue );

    py::class_< WavePressureOnFaceReal, WavePressureOnFaceReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "WavePressureOnFaceReal" )
        .def( py::init( &initFactoryPtr< WavePressureOnFaceReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< WavePressureOnFaceReal, std::string, ModelPtr > ) )
        .def( "build", &WavePressureOnFaceReal::build )
        .def( "setValue", &WavePressureOnFaceReal::setValue );

    py::class_< DistributedHeatFluxReal, DistributedHeatFluxReal::UnitaryMechanicalLoadRealPtr,
                MechanicalLoadReal >( mod, "DistributedHeatFluxReal" )
        .def( py::init( &initFactoryPtr< DistributedHeatFluxReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< DistributedHeatFluxReal, std::string, ModelPtr > ) )
        .def( "build", &DistributedHeatFluxReal::build )
        .def( "setValue", &DistributedHeatFluxReal::setValue );

    py::class_< DistributedHydraulicFluxReal,
                DistributedHydraulicFluxReal::UnitaryMechanicalLoadRealPtr, MechanicalLoadReal >(
        mod, "DistributedHydraulicFluxReal" )
        .def( py::init( &initFactoryPtr< DistributedHydraulicFluxReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< DistributedHydraulicFluxReal, std::string, ModelPtr > ) )
        .def( "build", &DistributedHydraulicFluxReal::build )
        .def( "setValue", &DistributedHydraulicFluxReal::setValue );
};
