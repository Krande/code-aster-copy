/**
 * @file MaterialFieldInterface.cxx
 * @brief Interface python de MaterialField
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

#include "PythonBindings/MaterialFieldInterface.h"

#include "aster_pybind.h"

void exportMaterialFieldToPython( py::module_ &mod ) {

    py::class_< PartOfMaterialField, PartOfMaterialFieldPtr >( mod, "PartOfMaterialField" )
        .def( py::init( &initFactoryPtr< PartOfMaterialField > ) )
        .def( py::init(
            &initFactoryPtr< PartOfMaterialField, std::vector< MaterialPtr >, MeshEntityPtr > ) )
        .def( define_pickling< PartOfMaterialField >() )
        .def( "getVectorOfMaterial", &PartOfMaterialField::getVectorOfMaterial )
        .def( "getMeshEntity", &PartOfMaterialField::getMeshEntity );

    py::class_< MaterialField, MaterialFieldPtr, DataStructure >( mod, "MaterialField" )
        .def( py::init( &initFactoryPtr< MaterialField, const BaseMeshPtr & > ) )
        .def(
            py::init( &initFactoryPtr< MaterialField, const std::string &, const BaseMeshPtr & > ) )
        .def( "addBehaviourOnMesh", &MaterialField::addBehaviourOnMesh, R"(
            Add behaviour (from DEFI_COMPOR) on mesh

            Arguments:
                behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
            )",
              py::arg( "behaviour" ) )

        .def( "addBehaviourOnGroupOfCells", &MaterialField::addBehaviourOnGroupOfCells, R"(
            Add behaviour (from DEFI_COMPOR) on group of cells

            Arguments:
                behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
                nameOfGroups (list(str)) : list of names of groups of cells
            )",
              py::arg( "behaviour" ), py::arg( "nameOfGroups" ) )

        .def( "addMultipleMaterialOnMesh", &MaterialField::addMultipleMaterialOnMesh, R"(
            Add a vector of multiple material properties on mesh

            Arguments:
                material (list(Material)): list of material properties
            )",
              py::arg( "material" ) )

        .def( "addMaterialOnMesh", &MaterialField::addMaterialOnMesh, R"(
            Add material properties on mesh

            Arguments:
                material (Material): material properties
            )",
              py::arg( "material" ) )

        .def( "addMultipleMaterialOnGroupOfCells",
              &MaterialField::addMultipleMaterialOnGroupOfCells, R"(
            Add a vector of multiple material properties on group of cells

            Arguments:
                material (list(Material)): list of material properties
                nameOfGroups (list(str)) : list of names of groups of cells
            )",
              py::arg( "material" ), py::arg( "nameOfGroups" ) )

        .def( "addMaterialOnGroupOfCells", &MaterialField::addMaterialOnGroupOfCells, R"(
            Add a material properties on list of groups of cells

            Arguments:
                material (Material): material properties
                nameOfGroups (list(str)) : list of names of groups of cells
            )",
              py::arg( "material" ), py::arg( "nameOfGroups" ) )

        .def( "getMesh", &MaterialField::getMesh, R"(
            Get mesh of material field

            Returns:
                BaseMesh: mesh
            )" )

        .def( "getVectorOfMaterial", &MaterialField::getVectorOfMaterial, R"(
            Get vector of all the material properties on the material field

            Returns:
                list(Material): vector of material properties
            )" )

        .def( "getVectorOfPartOfMaterialField", &MaterialField::getVectorOfPartOfMaterialField, R"(
            Get vector of all the material properties with mesh entities on the material field

            Returns:
                list(PartOfMaterial): vector of material properties with mesh entities
            )" )

        .def( "getMaterialsOnMeshEntities", &MaterialField::getMaterialsOnMeshEntities, R"()" )

        .def( "getExtStateVariablesOnMeshEntities",
              &MaterialField::getExtStateVariablesOnMeshEntities, R"()" )

        .def( "getMaterialOnCell", &MaterialField::getMaterialOnCell, R"(
            Get the material properties on a giver cell on the material field

            Returns:
                Material: material properties
            )" )

        .def( "setModel", &MaterialField::setModel, R"(
            Set model of the material field

            Arguments:
                model (Model): model
            )",
              py::arg( "model" ) )

        .def( "addExternalStateVariable", &MaterialField::addExternalStateVariable, R"(
            Add external state variable in material field

            Arguments:
                exteVari (ExternalStateVariablePtr): external state variable
            )",
              py::arg( "exteVari" ) )

        .def(
            "hasExternalStateVariable",
            py::overload_cast< const externVarEnumInt >( &MaterialField::hasExternalStateVariable ),
            R"(
            Detects the presence of an external state variable

            Arguments:
                exteVariIden (externVarEnumInt or str): identifier for external state variable

            Returns:
                bool: True if has external state variables
            )",
            py::arg( "exteVariIden" ) )

        .def( "hasExternalStateVariable",
              py::overload_cast< const std::string & >( &MaterialField::hasExternalStateVariable ),
              R"(
            Detects the presence of an external state variable

            Returns:
                bool: True if has external state variables
            )",
              py::arg( "exteVariIden" ) )

        .def( "hasExternalStateVariable",
              py::overload_cast<>( &MaterialField::hasExternalStateVariable, py::const_ ),
              R"(
            Detects the presence of any external state variable

            Returns:
                bool: True if has external state variables
            )" )

        .def( "hasExternalStateVariableForLoad", &MaterialField::hasExternalStateVariableForLoad,
              R"(
            Detects the presence of an external state variable for loads

            Returns:
                bool: True if has external state variables for loads
            )" )

        .def( "hasExternalStateVariableWithReference",
              &MaterialField::hasExternalStateVariableWithReference,
              R"(
            Detects the presence of an external state variable with reference value

            Returns:
                bool: True if has external state variables with reference value
            )" )

        .def( "build", &MaterialField::build, R"(
            Build material field
            )" )

        .def( "updateInternalState", &MaterialField::updateInternalState, R"(
Update the internal state of the datastructure.

Returns:
    bool: *True* in case of success, *False* otherwise.
            )" );
};
