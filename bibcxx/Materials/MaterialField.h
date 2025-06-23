#ifndef MATERIALFIELD_H_
#define MATERIALFIELD_H_

/**
 * @file MaterialField.h
 * @brief Header of MaterialField classes
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

#include "astercxx.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Materials/BehaviourDefinition.h"
#include "Materials/ExternalStateVariables.h"
#include "Materials/Material.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Mesh.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class PartOfMaterialField
 * @brief It contains a list of material properties on mesh
 */
class PartOfMaterialField {
  private:
    listOfMaterials _vectOfMater;
    MeshEntityPtr _meshEntity;

  public:
    PartOfMaterialField() : _meshEntity( nullptr ) {};

    PartOfMaterialField( const listOfMaterials &vectOfMater, const MeshEntityPtr &entity )
        : _vectOfMater( vectOfMater ), _meshEntity( entity ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    PartOfMaterialField( const py::tuple &tup )
        : PartOfMaterialField( tup[0].cast< listOfMaterials >(),
                               tup[1].cast< MeshEntityPtr >() ) {};
    py::tuple _getState() const { return py::make_tuple( _vectOfMater, _meshEntity ); };

    /** @brief Get the VectorOfMaterial of PartOfMaterialField */
    listOfMaterials getVectorOfMaterial() const { return _vectOfMater; };

    /** @brief Get the MeshEntity of PartOfMaterialField */
    MeshEntityPtr getMeshEntity() const { return _meshEntity; };
};

using PartOfMaterialFieldPtr = std::shared_ptr< PartOfMaterialField >;
using listOfPartOfMaterialField = std::vector< PartOfMaterialFieldPtr >;

/**
 * @class MaterialField
 * @brief Main class from AFFE_MATERIAU command
 */
class MaterialField : public DataStructure {
  private:
    using listOfMaterials = std::vector< MaterialPtr >;
    using listOfMaterialsOnMesh = std::list< std::pair< listOfMaterials, MeshEntityPtr > >;
    using listOfMaterialsOnMeshValue = listOfMaterialsOnMesh::value_type;
    using listOfBehavioursOnMesh = std::list< std::pair< BehaviourDefinitionPtr, MeshEntityPtr > >;
    using listOfBehavioursOnMeshValue = listOfBehavioursOnMesh::value_type;
    using listOfExternalVarOnMesh =
        std::list< std::pair< ExternalStateVariablePtr, MeshEntityPtr > >;
    using listOfExternalVarOnMeshValue = listOfExternalVarOnMesh::value_type;

    /** @brief Mesh */
    BaseMeshPtr _mesh;

    /** @brief Model */
    ModelPtr _model;

    /** @brief List of material parameters on mesh entities */
    listOfMaterialsOnMesh _materialsOnMeshEntities;

    /** @brief List of multi-material parameters on mesh entities */
    listOfBehavioursOnMesh _behaviourOnMeshEntities;

    /** @brief List of external state variables on mesh entities */
    listOfExternalVarOnMesh _extStateVariablesOnMeshEntities;

    /** @brief Constant field for material parameters - '.CHAMP_MAT' */
    ConstantFieldOnCellsChar8Ptr _champ_mat;

    /** @brief Constant field for multi-material parameters - '.COMPOR' */
    ConstantFieldOnCellsRealPtr _compor;

    /** @brief Jeveux vectors for external state variables (AFFE_VARC) */
    JeveuxVectorChar8 _cvrcNom;
    /** @brief Jeveux vector '.CVRCGD' */
    JeveuxVectorChar8 _cvrcGd;
    /** @brief Jeveux vector '.CVRCVARC' */
    JeveuxVectorChar8 _cvrcVarc;
    /** @brief Jeveux vector '.CVRCCMP' */
    JeveuxVectorChar8 _cvrcCmp;
    /** @brief Cartes R '.xxxx    .1' */
    std::map< std::string, ConstantFieldOnCellsRealPtr > _mapCvrcCard1;
    /** @brief Cartes K16 '.xxxx    .2' */
    std::map< std::string, ConstantFieldOnCellsChar16Ptr > _mapCvrcCard2;

  private:
    /** @brief Generate syntax for material parameters */
    ListSyntaxMapContainer syntaxForMaterial();

    /** @brief Generate syntax for multi-material parameters */
    ListSyntaxMapContainer syntaxForBehaviour();

    /** @brief Generate syntax for external state variables*/
    ListSyntaxMapContainer syntaxForExtStateVariables();

    /** @brief Generate C++ objects from Fortran objects for external state variables */
    void updateExtStateVariablesObjects();

  public:
    typedef std::shared_ptr< MaterialField > MaterialFieldPtr;

    /** @brief Constructor */
    MaterialField( const BaseMeshPtr &mesh )
        : MaterialField( ResultNaming::getNewResultName(), mesh ) {};

    /** @brief Constructor */
    MaterialField( const std::string &, const BaseMeshPtr & );

    /** @brief Constructor */
    MaterialField( const BaseMeshPtr &mesh, MaterialFieldPtr mater );

    /** @brief Destructor */
    ~MaterialField() {};

    /** @brief Add a behaviour on all mesh */
    void addBehaviourOnMesh( BehaviourDefinitionPtr &curBehav );

    /** @brief Add a behaviour on group of cells */
    void addBehaviourOnGroupOfCells( BehaviourDefinitionPtr &curBehav, VectorString namesOfGroup );

    /** @brief Add multiple material on all the mesh */
    void addMultipleMaterialOnMesh( std::vector< MaterialPtr > curMaters );

    /** @brief Add a material on all the mesh */
    void addMaterialOnMesh( MaterialPtr &curMater );

    /** @brief Add multiple material on group of cells */
    void addMultipleMaterialOnGroupOfCells( std::vector< MaterialPtr > curMaters,
                                            VectorString namesOfGroup );

    /** @brief Add a material on group of cells */
    void addMaterialOnGroupOfCells( MaterialPtr &curMater, VectorString namesOfGroup );

    /** @brief Add an external state variable */
    void addExternalStateVariable( ExternalStateVariablePtr &currExte );

    /** @brief Return the ConstantFieldOnCells of behaviour */
    ConstantFieldOnCellsRealPtr getBehaviourField() const { return _compor; };

    /** @brief Return a vector of Material */
    std::vector< MaterialPtr > getVectorOfMaterial() const;

    /** @brief Return a vector of PartOfMaterialField */
    std::vector< PartOfMaterialFieldPtr > getVectorOfPartOfMaterialField() const;

    /** @brief Return the MaterialPtr used on a cell */
    MaterialPtr getMaterialOnCell( const std::string cellName ) const;

    /** @brief Get mesh */
    BaseMeshPtr getMesh() {
        if ( _mesh->isEmpty() )
            throw std::runtime_error( "mesh of current model is empty" );
        return _mesh;
    };

    /** @brief Get model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get list of materials parameters on mesh entities */
    listOfMaterialsOnMesh getMaterialsOnMeshEntities() { return _materialsOnMeshEntities; };

    /** @brief Get list of multi-material parameters on mesh entities */
    listOfBehavioursOnMesh getBehaviourOnMeshEntities() { return _behaviourOnMeshEntities; };

    /** @brief Get list of external state variables on mesh entities */
    listOfExternalVarOnMesh getExtStateVariablesOnMeshEntities() {
        return _extStateVariablesOnMeshEntities;
    };

    /** @brief Set the model */
    void setModel( ModelPtr model ) {
        _model = model;
        _mesh = model->getMesh();
    };

    /** @brief Function to know if external state variable are present */
    bool hasExternalStateVariable() const { return _extStateVariablesOnMeshEntities.size() != 0; }

    /** @brief Function to know if a given external state variable exists */
    bool hasExternalStateVariable( const std::string &name );

    /** @brief Function to know if a given external state variable exists */
    bool hasExternalStateVariable( const externVarEnumInt exteVariType );

    /** @brief Function to know if there are external state variables for load (from strain) */
    bool hasExternalStateVariableForLoad();

    /** @brief Function to know if there are external state variables withe reference field */
    bool hasExternalStateVariableWithReference();

    /** @brief Main function to build this object */
    bool build();

    bool updateInternalState();
};

using MaterialFieldPtr = std::shared_ptr< MaterialField >;

#endif /* MATERIALFIELD_H_ */
