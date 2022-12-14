#ifndef MODEL_H_
#define MODEL_H_

/**
 * @file Model.h
 * @brief Fichier entete de la classe Model
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

#include "aster_fort_utils.h"
#include "astercxx.h"

#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "Loads/PhysicalQuantity.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/ConnectionMesh.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"
#include "Modeling/ElementaryModeling.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/XfemModel.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxDictionary.h"

class FiniteElementDescriptor;

using FiniteElementDescriptorPtr = std::shared_ptr< FiniteElementDescriptor >;

/**
 * @enum ModelSplitingMethod
 * @brief Types of partition for model
 */
enum ModelSplitingMethod { Centralized, SubDomain, GroupOfCellsSplit };
const int nbModelSplitingMethod = 3;

/**
 * @var ModelSplitingMethodNames
 * @brief Keyword for types of partition
 */
extern const char *const ModelSplitingMethodNames[nbModelSplitingMethod];

/**
 * @enum GraphPartitioner
 * @brief Graph partitioner
 */
enum GraphPartitioner { ScotchPartitioner, MetisPartitioner };
const int nbGraphPartitioner = 2;

/**
 * @var GraphPartitionerNames
 * @brief Keyword for graph partitioner
 */
extern const char *const GraphPartitionerNames[nbGraphPartitioner];

/**
 * @class Model
 * @brief Datastructure for model (AFFE_MODELE)
 */
class Model : public DataStructure, public ListOfTables {
  public:
    typedef std::shared_ptr< Model > ModelPtr;

  protected:
    // On redefinit le type MeshEntityPtr afin de pouvoir stocker les MeshEntity
    // dans la list
    /** @brief Pointeur intelligent vers un VirtualMeshEntity */
    typedef std::shared_ptr< VirtualMeshEntity > MeshEntityPtr;
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
    /** @brief Vecteur Jeveux '.PARTIT': TODO add PARTSD objects */
    /** @brief Liste contenant les modelisations ajoutees par l'utilisateur */
    listOfModsAndGrps _modelisations;
    /** @brief Maillage sur lequel repose la modelisation */
    BaseMeshPtr _baseMesh;
    /** @brief Model without XFEM */
    ModelPtr _saneModel;
    /** @brief Model with XFEM */
    XfemModelPtr _xfemModel;

/**
 * @brief Maillage sur lequel repose la modelisation
 * @todo a supprimer en templatisant Model etc.
 */
#ifdef ASTER_HAVE_MPI
    ConnectionMeshPtr _connectionMesh;
#endif /* ASTER_HAVE_MPI */
    /** @brief M??thode de parall??lisation du mod??le */
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
    Model( const std::string name, const bool is_xfem = false );

    /**
     * @brief Constructor: a mesh is mandatory
     */
    Model( void ) = delete;

    Model( const std::string name, const FiniteElementDescriptorPtr fed,
           const bool is_xfem = false );

    Model( const BaseMeshPtr mesh, const bool is_xfem = false );

#ifdef ASTER_HAVE_MPI
    Model( const std::string name, const ConnectionMeshPtr mesh );

    Model( const ConnectionMeshPtr mesh );
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Ajout d'une nouvelle modelisation sur tout le maillage
     * @param phys Physique a ajouters
     * @param mod Modelisation a ajouter
     */
    void addModelingOnMesh( Physics phys, Modelings mod );

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

    /**@brief Has multi-fiber beams in model ? */
    bool existsMultiFiberBeam();

    /**@brief Has THM in model ? */
    bool existsThm();

    /**@brief Has XFEM in model ? */
    bool existsXfem();

    /**@brief Has HHO in model ? */
    bool existsHHO();

    /**@brief Has Axis element in model ? */
    bool existsAxis();

    /**@brief Has COQUE_3D in model ? */
    bool exists3DShell();

    /**@brief Has STRX in model ? */
    bool existsSTRX();
    /**
     * @brief Number of super-elements in model
     * @return Number of super elements in model
     */
    ASTERINTEGER numberOfSuperElement();

    /**@brief Has super elements in model ? */
    bool existsSuperElement();

    /** @brief Has finite element in model ? */
    bool existsFiniteElement();

    /**
     * @brief function to know if XFEM Preconditioning is enable in model
     * @return true if xfem preconditioning enable
     */
    bool xfemPreconditioningEnable();

    /**
     * @brief Get FiniteElementDescriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const;

    /**
     * @brief Obtention de la methode du partitioner
     */
    GraphPartitioner getGraphPartitioner() const;

#ifdef ASTER_HAVE_MPI
    ConnectionMeshPtr getConnectionMesh() const;
#endif /* ASTER_HAVE_MPI */

    /**
     * @brief Get the sane base model
     */
    ModelPtr getSaneModel() const;

    /**
     * @brief Get the Xfem model
     */
    XfemModelPtr getXfemModel() const;

    /**
     * @brief Obtention de la methode de partition
     */
    ModelSplitingMethod getSplittingMethod() const;

    BaseMeshPtr getMesh() const;

    JeveuxVectorLong getTypeOfCells() const;

    /**
     * @brief Methode permettant de savoir si le modele est vide
     * @return true si le modele est vide
     */
    bool isEmpty() const;

    /**
     * @brief Set the sane base model
     */
    void setSaneModel( ModelPtr saneModel );

    /**
     * @brief Set XFEM model
     */
    void setXfemModel();

    /**
     * @brief Definition de la methode de partition
     */
    void setSplittingMethod( ModelSplitingMethod split, GraphPartitioner partitioner );

    /**
     * @brief Definition de la methode de partition
     */
    void setSplittingMethod( ModelSplitingMethod split );

    /**
     * @brief To known if the the model is mechanical or not
     *
     * @return true The phenomen is  mechanical
     */
    bool isMechanical( void ) const;

    /**
     * @brief To known if the the model is thermal or not
     *
     * @return true The phenomen is therman
     */
    bool isThermal( void ) const;

    /**
     * @brief To known if the the model is acoustic or not
     *
     * @return true The phenomen is acoustic
     */
    bool isAcoustic( void ) const;

    int getPhysics( void ) const;

    int getGeometricDimension( void ) const;

    /**
     * @brief Is an xfem  model ?
     * @return true if xfem model
     */
    bool isXfem() const;

#ifdef ASTER_HAVE_MPI
    bool setFrom( const ModelPtr model );
#endif
};

/**
 * @typedef Model
 */
using ModelPtr = std::shared_ptr< Model >;

#endif /* MODEL_H_ */