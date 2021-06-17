#ifndef MATERIALONMESH_H_
#define MATERIALONMESH_H_

/**
 * @file MaterialField.h
 * @brief Fichier entete de la classe MaterialField
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Materials/BehaviourDefinition.h"
#include "Materials/Material.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/Mesh.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"
#include "astercxx.h"
#include <stdexcept>

class MaterialFieldBuilder;

/**
 * @class PartOfMaterialField
 * @brief It contains a BehaviourDefinitionPtr and a MeshEntityPtr
 * @author Nicolas Sellenet
 */
class PartOfMaterialField {
  private:
    std::vector< MaterialPtr > _vecOfMater;
    MeshEntityPtr _entity;

  public:
    PartOfMaterialField() : _entity( nullptr ){};

    PartOfMaterialField( const std::vector< MaterialPtr > &vecOfMater,
                              const MeshEntityPtr &entity )
        : _vecOfMater( vecOfMater ), _entity( entity ){};

    /**
     * @brief Get the VectorOfMaterial of PartOfMaterialField
     */
    std::vector< MaterialPtr > getVectorOfMaterial() const { return _vecOfMater; };

    /**
     * @brief Get the MeshEntity of PartOfMaterialField
     */
    MeshEntityPtr getMeshEntity() const { return _entity; };
};

/**
 * @typedef PartOfMaterialFieldPtr
 * @brief Smart pointer on PartOfMaterialField
 */
typedef boost::shared_ptr< PartOfMaterialField > PartOfMaterialFieldPtr;

/**
 * @class MaterialField
 * @brief produit une sd identique a celle produite par AFFE_MATERIAU
 * @author Nicolas Sellenet
 */
class MaterialField : public DataStructure {
    friend class MaterialFieldBuilder;

  private:
    // On redefinit le type MeshEntityPtr afin de pouvoir stocker les MeshEntity
    // dans la list
    /** @typedef Definition d'un pointeur intelligent sur un VirtualMeshEntity */
    typedef boost::shared_ptr< VirtualMeshEntity > MeshEntityPtr;
    /** @typedef std::list d'une std::pair de MeshEntityPtr */
    typedef std::list< std::pair< std::vector< MaterialPtr >, MeshEntityPtr > > listOfMatsAndGrps;
    /** @typedef Definition de la valeur contenue dans un listOfMatsAndGrps */
    typedef listOfMatsAndGrps::value_type listOfMatsAndGrpsValue;
    /** @typedef Definition d'un iterateur sur listOfMatsAndGrps */
    typedef listOfMatsAndGrps::iterator listOfMatsAndGrpsIter;

    /** @typedef std::list d'une std::pair de MeshEntityPtr */
    typedef std::vector< std::pair< BehaviourDefinitionPtr, MeshEntityPtr > > listOfBehavAndGrps;
    /** @typedef Definition de la valeur contenue dans un listOfBehavAndGrps */
    typedef listOfBehavAndGrps::value_type listOfBehavAndGrpsValue;
    /** @typedef Definition d'un iterateur sur listOfBehavAndGrps */
    typedef listOfBehavAndGrps::iterator listOfBehavAndGrpsIter;

    /** @brief Maillage sur lequel repose la sd_cham_mater */
    BaseMeshPtr _mesh;
    /** @brief Model */
    ModelPtr _model;
    /** @brief Carte '.CHAMP_MAT' */
    ConstantFieldOnCellsChar8Ptr _listOfMaterials;
    /** @brief Carte '.TEMPE_REF' */
    ConstantFieldOnCellsRealPtr _listOfTemperatures;
    /** @brief Carte '.COMPOR' */
    ConstantFieldOnCellsRealPtr _behaviourField;
    /** @brief Liste contenant les materiaux ajoutes par l'utilisateur */
    listOfMatsAndGrps _materialsFieldEntity;
    /** @brief Link to a  */
    listOfBehavAndGrps _behaviours;
    /** @brief Jeveux vector '.CVRCNOM' */
    JeveuxVectorChar8 _cvrcNom;
    /** @brief Jeveux vector '.CVRCGD' */
    JeveuxVectorChar8 _cvrcGd;
    /** @brief Jeveux vector '.CVRCVARC' */
    JeveuxVectorChar8 _cvrcVarc;
    /** @brief Jeveux vector '.CVRCCMP' */
    JeveuxVectorChar8 _cvrcCmp;

  public:
    /**
     * @typedef MaterialFieldPtr
     * @brief Pointeur intelligent vers un MaterialField
     */
    typedef boost::shared_ptr< MaterialField > MaterialFieldPtr;

    /**
     * @brief Constructeur
     */
    MaterialField( const MeshPtr &mesh )
        : MaterialField( ResultNaming::getNewResultName(), mesh ){};

    /**
     * @brief Constructeur
     */
    MaterialField( const SkeletonPtr &mesh )
        : MaterialField( ResultNaming::getNewResultName(), mesh ){};

    /**
     * @brief Constructeur
     */
    MaterialField( const std::string &, const MeshPtr & );

    /**
     * @brief Constructeur
     */
    MaterialField( const std::string &, const SkeletonPtr & );

#ifdef ASTER_HAVE_MPI
    /**
     * @brief Constructeur
     */
    MaterialField( const ParallelMeshPtr &mesh )
        : MaterialField( ResultNaming::getNewResultName(), mesh ){};

    /**
     * @brief Constructeur
     */
    MaterialField( const std::string &, const ParallelMeshPtr & );
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Destructor
     */
    ~MaterialField(){};

    /**
     * @brief Add a behaviour on all mesh
     * @param curBehav behaviour to add
     */
    void addBehaviourOnMesh( BehaviourDefinitionPtr &curBehav );

    /**
     * @brief Ajout d'un materiau sur une entite du maillage
     * @param curMater behaviour to add
     * @param nameOfGroup Name of group
     */
    void addBehaviourOnGroupOfCells( BehaviourDefinitionPtr &curBehav, std::string nameOfGroup );

    /**
     * @brief Ajout d'un materiau sur une entite du maillage
     * @param curMater behaviour to add
     * @param nameOfGroup Name of group
     */
    void addBehaviourOnCell( BehaviourDefinitionPtr &curBehav, std::string nameOfCell );

    /**
     * @brief Add one or mmore materials on all the mesh
     * @param curMaters Material to be added
     */
    void addMaterialsOnMesh( std::vector< MaterialPtr > curMaters );
    void addMaterialsOnMesh( MaterialPtr &curMater );

    /**
     * @brief Ajout d'un materiau sur une entite du maillage
     * @param curMaters Materiau a ajouter
     * @param nameOfGroup Nom du groupe de mailles
     */
    void addMaterialsOnGroupOfCells( std::vector< MaterialPtr > curMaters,
                                     VectorString namesOfGroup );
    void addMaterialsOnGroupOfCells( MaterialPtr &curMater, VectorString namesOfGroup );

    /**
     * @brief Ajout d'un materiau sur une entite du maillage
     * @param curMaters Materiau a ajouter
     * @param nameOfCell Nom des mailles
     */
    void addMaterialsOnCell( std::vector< MaterialPtr > curMaters, VectorString namesOfCells );
    void addMaterialsOnCell( MaterialPtr &curMater, VectorString namesOfCells );

    /**
     * @brief Build MaterialFieldPtr without ExternalStateVariable
     * @return true
     */
    bool buildWithoutExternalStateVariables();

    /**
     * @brief Return the ConstantFieldOnCells of behaviour
     */
    ConstantFieldOnCellsRealPtr getBehaviourField() const {
        _behaviourField->updateValuePointers();
        return _behaviourField;
    };

    /**
     * @brief Return a vector of MaterialPtr
     */
    std::vector< MaterialPtr > getVectorOfMaterial() const;

    /**
     * @brief Return a vector of MaterialPtr
     */
    std::vector< PartOfMaterialFieldPtr > getVectorOfPartOfMaterialField() const;

    /**
     * @brief Function to know if Calculation Input Variables are present
     * @return true if present
     */
    bool hasExternalStateVariables() const;

    /**
     * @brief Function to know if a given Calculation Input Variables exists
     * @return true if exists
     */
    bool hasExternalStateVariables( const std::string & );

    /**
     * @brief Obtenir le maillage
     * @return Maillage du champ de materiau
     */
    BaseMeshPtr getMesh() {
        if ( _mesh->isEmpty() )
            throw std::runtime_error( "mesh of current model is empty" );
        return _mesh;
    };

    /**
     * @brief Get model
     */
    ModelPtr getModel() const { return _model;};

    /**
     * @brief Set the model
     */
    void setModel( ModelPtr model ) { _model = model; };

    /** @brief Add external state variables */
    void addExternalStateVariables( PyObject *keywords );
};

/**
 * @typedef MaterialFieldPtr
 * @brief Pointeur intelligent vers un MaterialField
 */
typedef boost::shared_ptr< MaterialField > MaterialFieldPtr;

#endif /* MATERIALONMESH_H_ */
