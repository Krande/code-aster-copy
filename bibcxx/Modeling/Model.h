#ifndef MODEL_H_
#define MODEL_H_

/**
 * @file Model.h
 * @brief Fichier entete de la classe Model
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

#include <map>
#include <stdexcept>

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "Loads/PhysicalQuantity.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/ConnectionMesh.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"
#include "Modeling/ElementaryModeling.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxDictionary.h"
#include "DataFields/ListOfTables.h"


/**
 * @enum ModelSplitingMethod
 * @brief Methodes de partitionnement du modèle
 * @author Nicolas Sellenet
 */
enum ModelSplitingMethod { Centralized, SubDomain, GroupOfCellsSplit };
const int nbModelSplitingMethod = 3;
/**
 * @var ModelSplitingMethodNames
 * @brief Nom Aster des differents partitionnement
 */
extern const char *const ModelSplitingMethodNames[nbModelSplitingMethod];

/**
 * @enum GraphPartitioner
 * @brief Partitionneur de graph
 * @author Nicolas Sellenet
 */
enum GraphPartitioner { ScotchPartitioner, MetisPartitioner };
const int nbGraphPartitioner = 2;
/**
 * @var GraphPartitionerNames
 * @brief Nom Aster des differents partitionneur de graph
 */
extern const char *const GraphPartitionerNames[nbGraphPartitioner];

/**
 * @class ModelClass
 * @brief Produit une sd identique a celle produite par AFFE_MODELE
 * @author Nicolas Sellenet
 */
class ModelClass : public DataStructure, public ListOfTablesClass {
  public:
    /**
     * @brief Forward declaration for the XFEM enrichment
     */
    typedef boost::shared_ptr< ModelClass > ModelPtr;

  protected:
    // On redefinit le type MeshEntityPtr afin de pouvoir stocker les MeshEntity
    // dans la list
    /** @brief Pointeur intelligent vers un VirtualMeshEntity */
    typedef boost::shared_ptr< VirtualMeshEntity > MeshEntityPtr;
    /** @brief std::list de std::pair de ElementaryModeling et MeshEntityPtr */
    typedef std::vector< std::pair< ElementaryModeling, MeshEntityPtr > > listOfModsAndGrps;
    /** @brief Valeur contenue dans listOfModsAndGrps */
    typedef listOfModsAndGrps::value_type listOfModsAndGrpsValue;
    /** @brief Iterateur sur un listOfModsAndGrps */
    typedef listOfModsAndGrps::iterator listOfModsAndGrpsIter;
    /** @brief Iterateur constant sur un listOfModsAndGrps */
    typedef listOfModsAndGrps::const_iterator listOfModsAndGrpsCIter;

    /** @brief Vecteur Jeveux '.MAILLE' */
    JeveuxVectorLong _typeOfCells;
    /** @brief Vecteur Jeveux '.NOEUD' */
    JeveuxVectorLong _typeOfNodes;
    /** @brief Vecteur Jeveux '.PARTIT' */
    JeveuxVectorChar8 _partition;
    /** @brief Liste contenant les modelisations ajoutees par l'utilisateur */
    listOfModsAndGrps _modelisations;
    /** @brief Maillage sur lequel repose la modelisation */
    BaseMeshPtr _baseMesh;
    /** @brief Maillage sur lequel repose la modelisation */
    ModelPtr _saneModel;
/**
 * @brief Maillage sur lequel repose la modelisation
 * @todo a supprimer en templatisant Model etc.
 */
#ifdef ASTER_HAVE_MPI
    ConnectionMeshPtr _connectionMesh;
#endif /* ASTER_HAVE_MPI */
    /** @brief Méthode de parallélisation du modèle */
    ModelSplitingMethod _splitMethod;
    /** @brief Graph partitioning */
    GraphPartitioner _graphPartitioner;
    /** @brief Object .MODELE */
    FiniteElementDescriptorPtr _ligrel;

    /**
     * @brief Ajout d'une nouvelle modelisation sur tout le maillage
     * @return SyntaxMapContainer contenant la syntaxe pour AFFE et les mc obligatoires
     */
    SyntaxMapContainer buildModelingsSyntaxMapContainer() const;

    /**
     * @brief Construction (au sens Jeveux fortran) de la sd_modele
     * @return booleen indiquant que la construction s'est bien deroulee
     */
    bool buildWithSyntax( SyntaxMapContainer & );

  public:
    /**
     * @brief Constructor: a mesh is mandatory
     */
    ModelClass( void ) = delete;

    ModelClass( const std::string name, const BaseMeshPtr mesh )
        : DataStructure( name, 8, "MODELE" ),
          ListOfTablesClass( name ),
          _typeOfCells( JeveuxVectorLong( getName() + ".MAILLE    " ) ),
          _typeOfNodes( JeveuxVectorLong( getName() + ".NOEUD     " ) ),
          _partition( JeveuxVectorChar8( getName() + ".PARTIT    " ) ), _saneModel( nullptr ),
          _baseMesh( mesh ), _splitMethod( SubDomain ), _graphPartitioner( MetisPartitioner ),
          _ligrel( boost::make_shared< FiniteElementDescriptorClass >(
                getName() + ".MODELE", _baseMesh ) ) {
        if ( _baseMesh->isEmpty() )
            throw std::runtime_error( "Mesh is empty" );
#ifdef ASTER_HAVE_MPI
        _connectionMesh = nullptr;
#endif
    };

    ModelClass( const BaseMeshPtr mesh ) : ModelClass( ResultNaming::getNewResultName(), mesh ){};

#ifdef ASTER_HAVE_MPI
    ModelClass( const std::string name, const ConnectionMeshPtr mesh )
        : DataStructure( name, 8, "MODELE" ),
          _typeOfCells( JeveuxVectorLong( getName() + ".MAILLE    " ) ),
          _typeOfNodes( JeveuxVectorLong( getName() + ".NOEUD     " ) ),
          _partition( JeveuxVectorChar8( getName() + ".PARTIT    " ) ), _saneModel( nullptr ),
          _baseMesh( mesh ), _connectionMesh( mesh ), _splitMethod( Centralized ),
          _graphPartitioner( MetisPartitioner ),
          _ligrel( boost::make_shared< FiniteElementDescriptorClass >(
                                                        getName() + ".MODELE", _baseMesh ) ) {
        AS_ASSERT( !_baseMesh->isEmpty() );
        AS_ASSERT( !_connectionMesh->isEmpty() );
    };

    ModelClass( const ConnectionMeshPtr mesh )
        : ModelClass( ResultNaming::getNewResultName(), mesh ){};
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Ajout d'une nouvelle modelisation sur tout le maillage
     * @param phys Physique a ajouters
     * @param mod Modelisation a ajouter
     */
    void addModelingOnMesh( Physics phys, Modelings mod ) {
        _modelisations.push_back( listOfModsAndGrpsValue(
            ElementaryModeling( phys, mod ), MeshEntityPtr( new AllMeshEntities() ) ) );
    };

    /**
     * @brief Ajout d'une nouvelle modelisation sur une entite du maillage
     * @param phys Physique a ajouter
     * @param mod Modelisation a ajouter
     * @param nameOfGroup Nom du groupe de mailles
     */
    void addModelingOnGroupOfCells( Physics phys, Modelings mod, std::string nameOfGroup );

    /**
     * @brief Ajout d'une nouvelle modelisation sur une entite du maillage
     * @param phys Physique a ajouter
     * @param mod Modelisation a ajouter
     * @param nameOfGroup Nom du groupe de noeuds
     */
    void addModelingOnGroupOfNodes( Physics phys, Modelings mod, std::string nameOfGroup );

    /**
     * @brief Construction (au sens Jeveux fortran) de la sd_modele
     * @return booleen indiquant que la construction s'est bien deroulee
     */
    virtual bool build();

    /**
     * @brief Function to know if there is MultiFiberBeam in the Model
     * @return true if MultiFiberBeam present
     */
    bool existsMultiFiberBeam();

    /**
     * @brief Is THM present in model
     * @return true if thm
     */
    bool existsThm();

    /**
     * @brief function to know if XFEM Preconditioning is enable in model
     * @return true if xfem preconditioning enable
     */
    bool xfemPreconditioningEnable();

    /**
     * @brief Get FiniteElementDescriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _ligrel; };

    /**
     * @brief Obtention de la methode du partitioner
     */
    GraphPartitioner getGraphPartitioner() const { return _graphPartitioner; };

#ifdef ASTER_HAVE_MPI
    ConnectionMeshPtr getConnectionMesh() const {
        if ( ( !_connectionMesh ) || _connectionMesh->isEmpty() )
            throw std::runtime_error( "Mesh of model is empty" );
        return _connectionMesh;
    };
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Get the sane base model
     */
    ModelPtr getSaneModel() const { return _saneModel; };

    /**
     * @brief Obtention de la methode de partition
     */
    ModelSplitingMethod getSplittingMethod() const { return _splitMethod; };

    BaseMeshPtr getMesh() const {
        if ( ( !_baseMesh ) || _baseMesh->isEmpty() )
            throw std::runtime_error( "Mesh of model is empty" );
        return _baseMesh;
    };

    JeveuxVectorLong getTypeOfCells() const
    {
        return _typeOfCells;
    }

    /**
     * @brief Methode permettant de savoir si le modele est vide
     * @return true si le modele est vide
     */
    bool isEmpty() const { return !_typeOfCells->exists(); };

    /**
     * @brief Set the sane base model
     */
    void setSaneModel( ModelPtr saneModel ) { _saneModel = saneModel; };

    /**
     * @brief Definition de la methode de partition
     */
    void setSplittingMethod( ModelSplitingMethod split, GraphPartitioner partitioner ) {
        setSplittingMethod( split );
        _graphPartitioner = partitioner;
    };

    /**
     * @brief Definition de la methode de partition
     */
    void setSplittingMethod( ModelSplitingMethod split ) {
#ifdef ASTER_HAVE_MPI
        if ( _connectionMesh && !_connectionMesh->isEmpty() && split != Centralized )
            throw std::runtime_error( "For Parallel mesh, Centralized splitting is mandatory" );
#endif /* ASTER_HAVE_MPI */

        _splitMethod = split;
    };


    /**
     * @brief To known if the the model is mechanical or not
     *
     * @return true The phenomen is  mechanical
     */
    bool isMechanical( void ) const {
        return this->getPhysics() == Physics::Mechanics;
    };

    /**
     * @brief To known if the the model is thermal or not
     *
     * @return true The phenomen is therman
     */
    bool isThermal( void ) const {
        return this->getPhysics() == Physics::Thermal;
    };

    /**
     * @brief To known if the the model is acoustic or not
     *
     * @return true The phenomen is acoustic
     */
    bool isAcoustic( void ) const {
        return this->getPhysics() == Physics::Acoustic;
    };

    int getPhysics( void ) const
    {
        return _ligrel->getPhysics();
    }

#ifdef ASTER_HAVE_MPI
    bool transferFrom( const ModelPtr model);
#endif
};

/**
 * @typedef Model
 * @brief Pointeur intelligent vers un ModelClass
 */
typedef boost::shared_ptr< ModelClass > ModelPtr;

#endif /* MODEL_H_ */
