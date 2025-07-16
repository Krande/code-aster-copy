/**
 * @file SuperMeshInterface.cxx
 * @brief Interface python de SuperMesh
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

#include "PythonBindings/SuperMeshInterface.h"

// #include "aster_pybind.h"

#include "Modal/DynamicMacroElement.h"
#include "Modal/StaticMacroElement.h"

#include <Meshes/SuperMesh.h>

void exportSuperMeshToPython( py::module_ &mod ) {

    py::class_< SuperMesh, std::shared_ptr< SuperMesh >, Mesh >( mod, "SuperMesh" )
        .def( py::init( &initFactoryPtr< SuperMesh > ) )
        .def( py::init( &initFactoryPtr< SuperMesh, std::string > ) )
        .def( "build", &SuperMesh::build, R"(
Returns:
    bool: true if building is ok
        )" )
        .def( "addDynamicMacroElement", &SuperMesh::addDynamicMacroElement, R"(
Add a dynamic macro element.
        )" )

        .def( "getDynamicMacroElements", &SuperMesh::getDynamicMacroElements, R"(
Return all dynamic macro elements.
        )" )

        .def( "addStaticMacroElement", &SuperMesh::addStaticMacroElement, R"(
Add a static macro element.
        )" )

        .def( "getStaticMacroElements", &SuperMesh::getStaticMacroElements, R"(
Return all static macro elements.
        )" );
};
