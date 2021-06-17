/**
 * @file MaterialFieldInterface.cxx
 * @brief Interface python de MaterialField
 * @author Nicolas Sellenet
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

#include "PythonBindings/MaterialFieldInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportMaterialFieldToPython() {

    py::class_< PartOfMaterialField, PartOfMaterialFieldPtr >( "PartOfMaterialField",
                                                                    py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< PartOfMaterialField > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< PartOfMaterialField,
                                                     std::vector< MaterialPtr >, MeshEntityPtr > ) )
        .def( "getVectorOfMaterial", &PartOfMaterialField::getVectorOfMaterial )
        .def( "getMeshEntity", &PartOfMaterialField::getMeshEntity );

    void ( MaterialField::*addmat1all )( std::vector< MaterialPtr > curMaters ) =
        &MaterialField::addMaterialsOnMesh;
    void ( MaterialField::*addmat2all )( MaterialPtr & curMater ) =
        &MaterialField::addMaterialsOnMesh;

    void ( MaterialField::*addmat1grp )( std::vector< MaterialPtr > curMaters,
                                              VectorString namesOfGroup ) =
        &MaterialField::addMaterialsOnGroupOfCells;
    void ( MaterialField::*addmat2grp )( MaterialPtr & curMater, VectorString namesOfGroup ) =
        &MaterialField::addMaterialsOnGroupOfCells;

    void ( MaterialField::*addmat1cell )( std::vector< MaterialPtr > curMaters,
                                               VectorString namesOfCells ) =
        &MaterialField::addMaterialsOnCell;
    void ( MaterialField::*addmat2cell )( MaterialPtr & curMater, VectorString namesOfCells ) =
        &MaterialField::addMaterialsOnCell;


    bool ( MaterialField::*hasESV1 )(  ) const =
        &MaterialField::hasExternalStateVariables;
    bool ( MaterialField::*hasESV2 )( const std::string & ) =
        &MaterialField::hasExternalStateVariables;

    py::class_< MaterialField, MaterialField::MaterialFieldPtr,
                py::bases< DataStructure > >( "MaterialField", py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< MaterialField, const MeshPtr & > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< MaterialField, const SkeletonPtr & > ) )
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MaterialField, const std::string &, const MeshPtr & > ) )
#ifdef ASTER_HAVE_MPI
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< MaterialField, const ParallelMeshPtr & > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< MaterialField, const std::string &,
                                                     const ParallelMeshPtr & > ) )
#endif /* ASTER_HAVE_MPI */
        .def( "addBehaviourOnMesh", &MaterialField::addBehaviourOnMesh )
        .def( "addBehaviourOnGroupOfCells", &MaterialField::addBehaviourOnGroupOfCells )
        .def( "addBehaviourOnCell", &MaterialField::addBehaviourOnCell )

        .def( "addMaterialsOnMesh", addmat1all )
        .def( "addMaterialsOnMesh", addmat2all )

        .def( "addMaterialsOnGroupOfCells", addmat1grp )
        .def( "addMaterialsOnGroupOfCells", addmat2grp )

        .def( "addMaterialsOnCell", addmat1cell )
        .def( "addMaterialsOnCell", addmat2cell )

        .def( "buildWithoutExternalStateVariables",
        &MaterialField::buildWithoutExternalStateVariables )
        .def( "getMesh", &MaterialField::getMesh )
        .def( "getVectorOfMaterial", &MaterialField::getVectorOfMaterial )
        .def( "getVectorOfPartOfMaterialField",
              &MaterialField::getVectorOfPartOfMaterialField )
        .def( "hasExternalStateVariables", hasESV1)
        .def( "hasExternalStateVariables", hasESV2)
        .def( "setModel", &MaterialField::setModel )

        .def( "addExternalStateVariables",
              static_cast< void ( MaterialField::* )( PyObject * ) >(
                  &MaterialField::addExternalStateVariables ),
              R"(
Add external state variables of material field

Arguments:
    AFFE_VARC (list[dict]): keywords as provided to AFFE_MATERIAU/AFFE_VARC
        )",
              ( py::arg( "self" ), py::arg( "AFFE_VARC" ) ) );


};
