/**
 * @file ThermalLoadInterface.cxx
 * @brief Interface python de ThermalLoad
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

#include "PythonBindings/ThermalLoadInterface.h"

#include "aster_pybind.h"

void exportThermalLoadToPython( py::module_ &mod ) {

    py::class_< ThermalLoadReal, ThermalLoadReal::ThermalLoadPtr, DataStructure >(
        mod, "ThermalLoadReal" )
        .def( py::init( &initFactoryPtr< ThermalLoadReal, ModelPtr & > ) )
        .def( py::init( &initFactoryPtr< ThermalLoadReal, std::string, ModelPtr & > ) )
        .def( "getFiniteElementDescriptor", &ThermalLoadReal::getFiniteElementDescriptor )
        .def( "getMesh", &ThermalLoadReal::getMesh )
        .def( "getModel", &ThermalLoadReal::getModel );

    py::class_< ThermalLoadFunction, ThermalLoadFunction::ThermalLoadPtr, DataStructure >(
        mod, "ThermalLoadFunction" )
        .def( py::init( &initFactoryPtr< ThermalLoadFunction, ModelPtr & > ) )
        .def( py::init( &initFactoryPtr< ThermalLoadFunction, std::string, ModelPtr & > ) )
        .def( "getFiniteElementDescriptor", &ThermalLoadFunction::getFiniteElementDescriptor )
        .def( "getMesh", &ThermalLoadFunction::getMesh )
        .def( "getModel", &ThermalLoadFunction::getModel );
};
