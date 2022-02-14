/**
 * @file PhysicsAndModelingsInterface.cxx
 * @brief Interface python de PhysicsAndModelings
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

#include "PythonBindings/PhysicsAndModelingsInterface.h"

#include "aster_pybind.h"

void exportPhysicsAndModelingsToPython( py::module_ &mod ) {

    py::enum_< Physics >( mod, "Physics" )
        .value( "Mechanics", Mechanics )
        .value( "Thermal", Thermal )
        .value( "Acoustic", Acoustic )
        .export_values();

    py::enum_< Modelings >( mod, "Modelings" )
        .value( "Axisymmetrical", Axisymmetrical )
        .value( "Tridimensional", Tridimensional )
        .value( "TridimensionalAbsorbingBoundary", TridimensionalAbsorbingBoundary )
        .value( "Planar", Planar )
        .value( "PlaneStrain", PlaneStrain )
        .value( "PlaneStress", PlaneStress )
        .value( "DKT", DKT )
        .value( "DKTG", DKTG )
        .value( "PlanarBar", PlanarBar )
        .export_values();
};
