/**
 * @file UnitaryMechanicalLoadInterface.cxx
 * @brief Interface python de MechanicalLoad
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
#include <PythonBindings/factory.h>
#include "PythonBindings/UnitaryMechanicalLoadInterface.h"

void exportUnitaryMechanicalLoadToPython() {

    py::enum_< LoadEnum >( "Loads" )
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
        .value( "ImpedanceOnFace", ImpedanceOnFace )
        .value( "NormalSpeedOnFace", NormalSpeedOnFace )
        .value( "WavePressureOnFace", WavePressureOnFace )
        .value( "THMFlux", THMFlux );


    py::class_< NodalForceReal, NodalForceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "NodalForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NodalForceReal, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< NodalForceReal, std::string, ModelPtr >))
        .def( "build", &NodalForceReal::build )
        .def( "setValue", &NodalForceReal::setValue );

    py::class_<
        NodalStructuralForceReal, NodalStructuralForceReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "NodalStructuralForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NodalStructuralForceReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< NodalStructuralForceReal, std::string, ModelPtr >))
        .def( "build", &NodalStructuralForceReal::build )
        .def( "setValue", &NodalStructuralForceReal::setValue );

    py::class_< ForceOnFaceReal, ForceOnFaceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "ForceOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ForceOnFaceReal, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ForceOnFaceReal, std::string, ModelPtr >))
        .def( "build", &ForceOnFaceReal::build )
        .def( "setValue", &ForceOnFaceReal::setValue );

    py::class_< ForceOnEdgeReal, ForceOnEdgeReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "ForceOnEdgeReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ForceOnEdgeReal, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ForceOnEdgeReal, std::string, ModelPtr >))
        .def( "build", &ForceOnEdgeReal::build )
        .def( "setValue", &ForceOnEdgeReal::setValue );

    py::class_<
        StructuralForceOnEdgeReal,
            StructuralForceOnEdgeReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "StructuralForceOnEdgeReal", py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnEdgeReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnEdgeReal, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnEdgeReal::build )
        .def( "setValue", &StructuralForceOnEdgeReal::setValue );

    py::class_< LineicForceReal, LineicForceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "LineicForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LineicForceReal, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< LineicForceReal, std::string, ModelPtr >))
        .def( "build", &LineicForceReal::build )
        .def( "setValue", &LineicForceReal::setValue );

    py::class_< InternalForceReal, InternalForceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "InternalForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< InternalForceReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< InternalForceReal, std::string, ModelPtr >))
        .def( "build", &InternalForceReal::build )
        .def( "setValue", &InternalForceReal::setValue );

    py::class_<
        StructuralForceOnBeamReal,
            StructuralForceOnBeamReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "StructuralForceOnBeamReal", py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnBeamReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnBeamReal, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnBeamReal::build )
        .def( "setValue", &StructuralForceOnBeamReal::setValue );

    py::class_< LocalForceOnBeamReal, LocalForceOnBeamReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "LocalForceOnBeamReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LocalForceOnBeamReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< LocalForceOnBeamReal, std::string, ModelPtr >))
        .def( "build", &LocalForceOnBeamReal::build )
        .def( "setValue", &LocalForceOnBeamReal::setValue );

    py::class_< StructuralForceOnShellReal,
                StructuralForceOnShellReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "StructuralForceOnShellReal",
                                                              py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnShellReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnShellReal, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnShellReal::build )
        .def( "setValue", &StructuralForceOnShellReal::setValue );

    py::class_< LocalForceOnShellReal,
        LocalForceOnShellReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "LocalForceOnShellReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LocalForceOnShellReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< LocalForceOnShellReal, std::string, ModelPtr >))
        .def( "build", &LocalForceOnShellReal::build )
        .def( "setValue", &LocalForceOnShellReal::setValue );

    py::class_< PressureOnShellReal, PressureOnShellReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "PressureOnShellReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< PressureOnShellReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< PressureOnShellReal, std::string, ModelPtr >))
        .def( "build", &PressureOnShellReal::build )
        .def( "setValue", &PressureOnShellReal::setValue );

    py::class_< PressureOnPipeReal, PressureOnPipeReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "PressureOnPipeReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< PressureOnPipeReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< PressureOnPipeReal, std::string, ModelPtr >))
        .def( "build", &PressureOnPipeReal::build )
        .def( "setValue", &PressureOnPipeReal::setValue );

    py::class_<
        ImposedDisplacementReal, ImposedDisplacementReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "ImposedDisplacementReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImposedDisplacementReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImposedDisplacementReal, std::string, ModelPtr >))
        .def( "build", &ImposedDisplacementReal::build )
        .def( "setValue", &ImposedDisplacementReal::setValue );

    py::class_< ImposedPressureReal, ImposedPressureReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "ImposedPressureReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImposedPressureReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImposedPressureReal, std::string, ModelPtr >))
        .def( "build", &ImposedPressureReal::build )
        .def( "setValue", &ImposedPressureReal::setValue );

    py::class_<
        DistributedPressureReal, DistributedPressureReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "DistributedPressureReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DistributedPressureReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedPressureReal, std::string, ModelPtr >))
        .def( "build", &DistributedPressureReal::build )
        .def( "setValue", &DistributedPressureReal::setValue );

    py::class_< ImpedanceOnFaceReal, ImpedanceOnFaceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "ImpedanceOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImpedanceOnFaceReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImpedanceOnFaceReal, std::string, ModelPtr >))
        .def( "build", &ImpedanceOnFaceReal::build )
        .def( "setValue", &ImpedanceOnFaceReal::setValue );

    py::class_< NormalSpeedOnFaceReal,
        NormalSpeedOnFaceReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "NormalSpeedOnFaceReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NormalSpeedOnFaceReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< NormalSpeedOnFaceReal, std::string, ModelPtr >))
        .def( "build", &NormalSpeedOnFaceReal::build )
        .def( "setValue", &NormalSpeedOnFaceReal::setValue );

    py::class_<
        WavePressureOnFaceReal, WavePressureOnFaceReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "WavePressureOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< WavePressureOnFaceReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< WavePressureOnFaceReal, std::string, ModelPtr >))
        .def( "build", &WavePressureOnFaceReal::build )
        .def( "setValue", &WavePressureOnFaceReal::setValue );

    py::class_<
        DistributedHeatFluxReal, DistributedHeatFluxReal::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadReal > >( "DistributedHeatFluxReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DistributedHeatFluxReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedHeatFluxReal, std::string, ModelPtr >))
        .def( "build", &DistributedHeatFluxReal::build )
        .def( "setValue", &DistributedHeatFluxReal::setValue );

    py::class_< DistributedHydraulicFluxReal,
                DistributedHydraulicFluxReal::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadReal > >( "DistributedHydraulicFluxReal",
                                                              py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< DistributedHydraulicFluxReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedHydraulicFluxReal, std::string, ModelPtr >))
        .def( "build", &DistributedHydraulicFluxReal::build )
        .def( "setValue", &DistributedHydraulicFluxReal::setValue );
};
