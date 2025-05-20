/**
 * @file MeshEntitiesInterface.cxx
 * @brief Interface python de MeshEntities
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

#include "PythonBindings/MeshEntitiesInterface.h"

#include "aster_pybind.h"

#include <Meshes/MeshEntities.h>

void exportMeshEntitiesToPython( py::module_ &mod ) {

    py::enum_< EntityType >( mod, "EntityType", R"(
Enumeration for entity type.
    )" )
        .value( "GroupOfNodesType", GroupOfNodesType )
        .value( "GroupOfCellsType", GroupOfCellsType )
        .value( "AllMeshEntitiesType", AllMeshEntitiesType )
        .value( "CellType", CellType )
        .value( "NodeType", NodeType )
        .value( "NoType", NoType )
        .export_values();

    py::class_< VirtualMeshEntity, MeshEntityPtr >( mod, "MeshEntity" )
        .def( py::init( &initFactoryPtr< VirtualMeshEntity, std::string, EntityType > ) )
        // fake initFactoryPtr: created by subclass
        .def( define_pickling< VirtualMeshEntity >() )
        .def( "getType", &VirtualMeshEntity::getType )
        .def( "getNames", &VirtualMeshEntity::getNames );

    py::class_< AllMeshEntities, AllMeshEntitiesPtr, VirtualMeshEntity >( mod, "AllMeshEntities" )
        // fake initFactoryPtr: created by subclass
        // fake initFactoryPtr: created by subclass
        ;
};
