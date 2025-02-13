/**
 * @file StructureInterfaceInterface.cxx
 * @brief Interface python de StructureInterface
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

#include "PythonBindings/StructureInterfaceInterface.h"

#include "aster_pybind.h"

void exportStructureInterfaceToPython( py::module_ &mod ) {

    py::enum_< InterfaceTypeEnum >( mod, "InterfaceType", R"(
Enumeration of interface type.
    )" )
        .value( "MacNeal", MacNeal )
        .value( "CraigBampton", CraigBampton )
        .value( "HarmonicCraigBampton", HarmonicCraigBampton )
        .value( "None", NoInterfaceType )
        .export_values();

    py::class_< StructureInterface, StructureInterface::StructureInterfacePtr, DataStructure >(
        mod, "StructureInterface" )
        .def( py::init( &initFactoryPtr< StructureInterface > ) )
        .def( py::init( &initFactoryPtr< StructureInterface, std::string > ) )
        .def( py::init( &initFactoryPtr< StructureInterface, DOFNumberingPtr > ) )
        .def( py::init( &initFactoryPtr< StructureInterface, std::string, DOFNumberingPtr > ) );
};
