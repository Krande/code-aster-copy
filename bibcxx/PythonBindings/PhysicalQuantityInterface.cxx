/**
 * @file PhysicalQuantityInterface.cxx
 * @brief Interface python de PhysicalQuantity
 * @author Nicolas Sellenet
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

// Not DataStructures
// aslint: disable=C3006

#include "PythonBindings/PhysicalQuantityInterface.h"

#include "aster_pybind.h"

void exportPhysicalQuantityToPython( py::module_ &mod ) {

    py::enum_< PhysicalQuantityComponent >( mod, "PhysicalQuantityComponent", R"(
Enumeration for physical component.
    )" )
        .value( "Dx", Dx )
        .value( "Dy", Dy )
        .value( "Dz", Dz )
        .value( "Drx", Drx )
        .value( "Dry", Dry )
        .value( "Drz", Drz )
        .value( "Temp", Temp )
        .value( "MiddleTemp", MiddleTemp )
        .value( "Pres", Pres )
        .value( "Fx", Fx )
        .value( "Fy", Fy )
        .value( "Fz", Fz )
        .value( "Mx", Mx )
        .value( "My", My )
        .value( "Mz", Mz )
        .value( "N", N )
        .value( "Vy", Vy )
        .value( "Vz", Vz )
        .value( "Mt", Mt )
        .value( "Mfy", Mfy )
        .value( "Mfz", Mfz )
        .value( "F1", F1 )
        .value( "F2", F2 )
        .value( "F3", F3 )
        .value( "Mf1", Mf1 )
        .value( "Mf2", Mf2 )
        .value( "Impe", Impe )
        .value( "Vnor", Vnor )
        .value( "Flun", Flun )
        .value( "FlunHydr1", FlunHydr1 )
        .value( "FlunHydr2", FlunHydr2 )
        .export_values();

    py::class_< ForceReal, ForceReal::PhysicalQuantityPtr >( mod, "ForceReal" )
        .def( py::init( &initFactoryPtr< ForceReal > ) )
        .def( define_pickling< ForceReal >() )
        .def( "debugPrint", &ForceReal::debugPrint )
        .def( "setValue", &ForceReal::setValue );

    py::class_< StructuralForceReal, StructuralForceReal::PhysicalQuantityPtr >(
        mod, "StructuralForceReal" )
        .def( py::init( &initFactoryPtr< StructuralForceReal > ) )
        .def( define_pickling< StructuralForceReal >() )
        .def( "debugPrint", &StructuralForceReal::debugPrint )
        .def( "setValue", &StructuralForceReal::setValue );

    py::class_< LocalBeamForceReal, LocalBeamForceReal::PhysicalQuantityPtr >(
        mod, "LocalBeamForceReal" )
        .def( py::init( &initFactoryPtr< LocalBeamForceReal > ) )
        .def( define_pickling< LocalBeamForceReal >() )
        .def( "debugPrint", &LocalBeamForceReal::debugPrint )
        .def( "setValue", &LocalBeamForceReal::setValue );

    py::class_< LocalShellForceReal, LocalShellForceReal::PhysicalQuantityPtr >(
        mod, "LocalShellForceReal" )
        .def( py::init( &initFactoryPtr< LocalShellForceReal > ) )
        .def( define_pickling< LocalShellForceReal >() )
        .def( "debugPrint", &LocalShellForceReal::debugPrint )
        .def( "setValue", &LocalShellForceReal::setValue );

    py::class_< DisplacementReal, DisplacementReal::PhysicalQuantityPtr >( mod, "DisplacementReal" )
        .def( py::init( &initFactoryPtr< DisplacementReal > ) )
        .def( define_pickling< DisplacementReal >() )
        .def( "debugPrint", &DisplacementReal::debugPrint )
        .def( "setValue", &DisplacementReal::setValue );

    py::class_< PressureReal, PressureReal::PhysicalQuantityPtr >( mod, "PressureReal" )
        .def( py::init( &initFactoryPtr< PressureReal > ) )
        .def( define_pickling< PressureReal >() )
        .def( "debugPrint", &PressureReal::debugPrint )
        .def( "setValue", &PressureReal::setValue );

    py::class_< ImpedanceReal, ImpedanceReal::PhysicalQuantityPtr >( mod, "ImpedanceReal" )
        .def( py::init( &initFactoryPtr< ImpedanceReal > ) )
        .def( define_pickling< ImpedanceReal >() )
        .def( "debugPrint", &ImpedanceReal::debugPrint )
        .def( "setValue", &ImpedanceReal::setValue );

    py::class_< NormalSpeedReal, NormalSpeedReal::PhysicalQuantityPtr >( mod, "NormalSpeedReal" )
        .def( py::init( &initFactoryPtr< NormalSpeedReal > ) )
        .def( define_pickling< NormalSpeedReal >() )
        .def( "debugPrint", &NormalSpeedReal::debugPrint )
        .def( "setValue", &NormalSpeedReal::setValue );

    py::class_< HeatFluxReal, HeatFluxReal::PhysicalQuantityPtr >( mod, "HeatFluxReal" )
        .def( py::init( &initFactoryPtr< HeatFluxReal > ) )
        .def( define_pickling< HeatFluxReal >() )
        .def( "debugPrint", &HeatFluxReal::debugPrint )
        .def( "setValue", &HeatFluxReal::setValue );

    py::class_< HydraulicFluxReal, HydraulicFluxReal::PhysicalQuantityPtr >( mod,
                                                                             "HydraulicFluxReal" )
        .def( py::init( &initFactoryPtr< HydraulicFluxReal > ) )
        .def( define_pickling< HydraulicFluxReal >() )
        .def( "debugPrint", &HydraulicFluxReal::debugPrint )
        .def( "setValue", &HydraulicFluxReal::setValue );
};
