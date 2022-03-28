/**
 * @file MaterialFieldInterface.cxx
 * @brief Interface python de MaterialField
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

#include "PythonBindings/MaterialFieldInterface.h"

#include "aster_pybind.h"

void exportMaterialFieldToPython( py::module_ &mod ) {

    py::class_< PartOfMaterialField, PartOfMaterialFieldPtr >( mod, "PartOfMaterialField" )
        .def( py::init( &initFactoryPtr< PartOfMaterialField > ) )
        .def( py::init(
            &initFactoryPtr< PartOfMaterialField, std::vector< MaterialPtr >, MeshEntityPtr > ) )
        .def( "getVectorOfMaterial", &PartOfMaterialField::getVectorOfMaterial )
        .def( "getMeshEntity", &PartOfMaterialField::getMeshEntity );

    py::class_< MaterialField, MaterialField::MaterialFieldPtr, DataStructure >( mod,
                                                                                 "MaterialField" )
        .def( py::init( &initFactoryPtr< MaterialField, const MeshPtr & > ) )
        .def( py::init( &initFactoryPtr< MaterialField, const SkeletonPtr & > ) )
        .def( py::init( &initFactoryPtr< MaterialField, const std::string &, const MeshPtr & > ) )
#ifdef ASTER_HAVE_MPI
        .def( py::init( &initFactoryPtr< MaterialField, const ParallelMeshPtr & > ) )
        .def( py::init(
            &initFactoryPtr< MaterialField, const std::string &, const ParallelMeshPtr & > ) )
#endif /* ASTER_HAVE_MPI */
        .def( "addBehaviourOnMesh", &MaterialField::addBehaviourOnMesh )
        .def( "addBehaviourOnGroupOfCells", &MaterialField::addBehaviourOnGroupOfCells )

        .def( "addMaterialsOnMesh", py::overload_cast< std::vector< MaterialPtr > >(
                                        &MaterialField::addMaterialsOnMesh ) )
        .def( "addMaterialsOnMesh",
              py::overload_cast< MaterialPtr & >( &MaterialField::addMaterialsOnMesh ) )

        .def( "addMaterialsOnGroupOfCells",
              py::overload_cast< std::vector< MaterialPtr >, VectorString >(
                  &MaterialField::addMaterialsOnGroupOfCells ) )
        .def( "addMaterialsOnGroupOfCells", py::overload_cast< MaterialPtr &, VectorString >(
                                                &MaterialField::addMaterialsOnGroupOfCells ) )

        .def( "buildWithoutExternalStateVariables",
              &MaterialField::buildWithoutExternalStateVariables )
        .def( "getMesh", &MaterialField::getMesh )
        .def( "getVectorOfMaterial", &MaterialField::getVectorOfMaterial )
        .def( "getVectorOfPartOfMaterialField", &MaterialField::getVectorOfPartOfMaterialField )
        .def( "hasExternalStateVariables",
              py::overload_cast<>( &MaterialField::hasExternalStateVariables, py::const_ ) )
        .def( "hasExternalStateVariables", py::overload_cast< const std::string & >(
                                               &MaterialField::hasExternalStateVariables ) )
        .def( "setModel", &MaterialField::setModel )

        .def( "addExternalStateVariables", &MaterialField::addExternalStateVariables,
              R"(
Add external state variables of material field

Arguments:
    AFFE_VARC (list[dict]): keywords as provided to AFFE_MATERIAU/AFFE_VARC
        )",
              py::arg( "AFFE_VARC" ) )
        .def( "update", &MaterialField::update );;
};
