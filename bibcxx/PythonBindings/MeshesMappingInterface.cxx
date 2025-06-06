/**
 * @file MeshesMappingInterface.cxx
 * @brief Interface python de MeshesMapping
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

#include "PythonBindings/MeshesMappingInterface.h"

#include "aster_pybind.h"

void exportMeshesMappingToPython( py::module_ &mod ) {

    py::class_< MeshesMapping, MeshesMapping::MeshesMappingPtr, DataStructure >( mod,
                                                                                 "MeshesMapping" )
        .def( py::init( &initFactoryPtr< MeshesMapping > ) )
        .def( py::init( &initFactoryPtr< MeshesMapping, std::string > ) )
        .def( "setFirstMesh", &MeshesMapping::setFirstMesh )
        .def( "setSecondMesh", &MeshesMapping::setSecondMesh )
        .def( "getFirstMesh", &MeshesMapping::getFirstMesh )
        .def( "getSecondMesh", &MeshesMapping::getSecondMesh );
};
