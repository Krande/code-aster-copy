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


    py::class_< NodalForceRealClass, NodalForceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "NodalForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NodalForceRealClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< NodalForceRealClass, std::string, ModelPtr >))
        .def( "build", &NodalForceRealClass::build )
        .def( "setValue", &NodalForceRealClass::setValue );

    py::class_<
        NodalStructuralForceRealClass, NodalStructuralForceRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "NodalStructuralForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NodalStructuralForceRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< NodalStructuralForceRealClass, std::string, ModelPtr >))
        .def( "build", &NodalStructuralForceRealClass::build )
        .def( "setValue", &NodalStructuralForceRealClass::setValue );

    py::class_< ForceOnFaceRealClass, ForceOnFaceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "ForceOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ForceOnFaceRealClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ForceOnFaceRealClass, std::string, ModelPtr >))
        .def( "build", &ForceOnFaceRealClass::build )
        .def( "setValue", &ForceOnFaceRealClass::setValue );

    py::class_< ForceOnEdgeRealClass, ForceOnEdgeRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "ForceOnEdgeReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ForceOnEdgeRealClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< ForceOnEdgeRealClass, std::string, ModelPtr >))
        .def( "build", &ForceOnEdgeRealClass::build )
        .def( "setValue", &ForceOnEdgeRealClass::setValue );

    py::class_<
        StructuralForceOnEdgeRealClass,
            StructuralForceOnEdgeRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "StructuralForceOnEdgeReal", py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnEdgeRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnEdgeRealClass, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnEdgeRealClass::build )
        .def( "setValue", &StructuralForceOnEdgeRealClass::setValue );

    py::class_< LineicForceRealClass, LineicForceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "LineicForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LineicForceRealClass, ModelPtr >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< LineicForceRealClass, std::string, ModelPtr >))
        .def( "build", &LineicForceRealClass::build )
        .def( "setValue", &LineicForceRealClass::setValue );

    py::class_< InternalForceRealClass, InternalForceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "InternalForceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< InternalForceRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< InternalForceRealClass, std::string, ModelPtr >))
        .def( "build", &InternalForceRealClass::build )
        .def( "setValue", &InternalForceRealClass::setValue );

    py::class_<
        StructuralForceOnBeamRealClass,
            StructuralForceOnBeamRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "StructuralForceOnBeamReal", py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnBeamRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnBeamRealClass, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnBeamRealClass::build )
        .def( "setValue", &StructuralForceOnBeamRealClass::setValue );

    py::class_< LocalForceOnBeamRealClass, LocalForceOnBeamRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "LocalForceOnBeamReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LocalForceOnBeamRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< LocalForceOnBeamRealClass, std::string, ModelPtr >))
        .def( "build", &LocalForceOnBeamRealClass::build )
        .def( "setValue", &LocalForceOnBeamRealClass::setValue );

    py::class_< StructuralForceOnShellRealClass,
                StructuralForceOnShellRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "StructuralForceOnShellReal",
                                                              py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< StructuralForceOnShellRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< StructuralForceOnShellRealClass, std::string, ModelPtr >))
        .def( "build", &StructuralForceOnShellRealClass::build )
        .def( "setValue", &StructuralForceOnShellRealClass::setValue );

    py::class_< LocalForceOnShellRealClass,
        LocalForceOnShellRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "LocalForceOnShellReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< LocalForceOnShellRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< LocalForceOnShellRealClass, std::string, ModelPtr >))
        .def( "build", &LocalForceOnShellRealClass::build )
        .def( "setValue", &LocalForceOnShellRealClass::setValue );

    py::class_< PressureOnShellRealClass, PressureOnShellRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "PressureOnShellReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< PressureOnShellRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< PressureOnShellRealClass, std::string, ModelPtr >))
        .def( "build", &PressureOnShellRealClass::build )
        .def( "setValue", &PressureOnShellRealClass::setValue );

    py::class_< PressureOnPipeRealClass, PressureOnPipeRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "PressureOnPipeReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< PressureOnPipeRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< PressureOnPipeRealClass, std::string, ModelPtr >))
        .def( "build", &PressureOnPipeRealClass::build )
        .def( "setValue", &PressureOnPipeRealClass::setValue );

    py::class_<
        ImposedDisplacementRealClass, ImposedDisplacementRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "ImposedDisplacementReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImposedDisplacementRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImposedDisplacementRealClass, std::string, ModelPtr >))
        .def( "build", &ImposedDisplacementRealClass::build )
        .def( "setValue", &ImposedDisplacementRealClass::setValue );

    py::class_< ImposedPressureRealClass, ImposedPressureRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "ImposedPressureReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImposedPressureRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImposedPressureRealClass, std::string, ModelPtr >))
        .def( "build", &ImposedPressureRealClass::build )
        .def( "setValue", &ImposedPressureRealClass::setValue );

    py::class_<
        DistributedPressureRealClass, DistributedPressureRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "DistributedPressureReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DistributedPressureRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedPressureRealClass, std::string, ModelPtr >))
        .def( "build", &DistributedPressureRealClass::build )
        .def( "setValue", &DistributedPressureRealClass::setValue );

    py::class_< ImpedanceOnFaceRealClass, ImpedanceOnFaceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "ImpedanceOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ImpedanceOnFaceRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ImpedanceOnFaceRealClass, std::string, ModelPtr >))
        .def( "build", &ImpedanceOnFaceRealClass::build )
        .def( "setValue", &ImpedanceOnFaceRealClass::setValue );

    py::class_< NormalSpeedOnFaceRealClass,
        NormalSpeedOnFaceRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "NormalSpeedOnFaceReal",
                                                              py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< NormalSpeedOnFaceRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< NormalSpeedOnFaceRealClass, std::string, ModelPtr >))
        .def( "build", &NormalSpeedOnFaceRealClass::build )
        .def( "setValue", &NormalSpeedOnFaceRealClass::setValue );

    py::class_<
        WavePressureOnFaceRealClass, WavePressureOnFaceRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "WavePressureOnFaceReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< WavePressureOnFaceRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< WavePressureOnFaceRealClass, std::string, ModelPtr >))
        .def( "build", &WavePressureOnFaceRealClass::build )
        .def( "setValue", &WavePressureOnFaceRealClass::setValue );

    py::class_<
        DistributedHeatFluxRealClass, DistributedHeatFluxRealClass::UnitaryMechanicalLoadRealPtr,
        py::bases< MechanicalLoadRealClass > >( "DistributedHeatFluxReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DistributedHeatFluxRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedHeatFluxRealClass, std::string, ModelPtr >))
        .def( "build", &DistributedHeatFluxRealClass::build )
        .def( "setValue", &DistributedHeatFluxRealClass::setValue );

    py::class_< DistributedHydraulicFluxRealClass,
                DistributedHydraulicFluxRealClass::UnitaryMechanicalLoadRealPtr,
                py::bases< MechanicalLoadRealClass > >( "DistributedHydraulicFluxReal",
                                                              py::no_init )
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< DistributedHydraulicFluxRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< DistributedHydraulicFluxRealClass, std::string, ModelPtr >))
        .def( "build", &DistributedHydraulicFluxRealClass::build )
        .def( "setValue", &DistributedHydraulicFluxRealClass::setValue );
};
