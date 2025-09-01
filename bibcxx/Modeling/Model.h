#ifndef MODEL_H_
#define MODEL_H_

/**
 * @file Model.h
 * @brief Fichier entete de la classe Model
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

#include "astercxx.h"

#include "aster_fort_utils.h"

#include "DataFields/ListOfTables.h"
#include "DataStructures/DataStructure.h"
#include "Loads/PhysicalQuantity.h"
#include "Meshes/ConnectionMesh.h"
#include "Modeling/ElementaryModeling.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/HHOModel.h"
#include "Modeling/Partition.h"
#include "Modeling/XfemModel.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxDictionary.h"

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
    typedef std::vector< std::pair< ElementaryModeling, std::vector< MeshEntityPtr > > >
        listOfModsAndGrps;
    /** @brief Valeur contenue dans listOfModsAndGrps */
    typedef listOfModsAndGrps::value_type listOfModsAndGrpsValue;
    /** @brief Iterateur sur un listOfModsAndGrps */
    typedef listOfModsAndGrps::iterator listOfModsAndGrpsIter;
    /** @brief Iterateur constant sur un listOfModsAndGrps */
    typedef listOfModsAndGrps::const_iterator listOfModsAndGrpsCIter;
    /** @brief Vecteur Jeveux '.NOEUD' */
    JeveuxVectorLong _typeOfNodes;
    /** @brief Liste contenant les modelisations ajoutees par l'utilisateur */
    listOfModsAndGrps _modelisations;
    /** @brief Model without XFEM */
    ModelPtr _saneModel;
    /** @brief Model with XFEM */
    XfemModelPtr _xfemModel;
    /** @brief Field for HHO  */
    HHOModelPtr _hhoModel;
    /** @brief To know if model is balanceable */
    bool _isBalanceable = true;

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

    const std::string dismoi( const std::string &, bool stop = true ) const;

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

    Model( const BaseMeshPtr mesh, const ModelPtr model );

    /**
     * @brief Ajout d'une nouvelle modelisation sur tout le maillage
     * @param phys Physique a ajouters
     * @param mod Modelisation a ajouter
     * @param formumation Formulation
     */
    void addModelingOnMesh( Physics phys, Modelings mod, Formulation form = NoFormulation );

    /**
     * @brief Ajout d'une nouvelle modelisation sur une entite du maillage
     * @param phys Physique a ajouter
     * @param mod Modelisation a ajouter
     * @param nameOfGroup Nom du groupe de mailles
     * @param formumation Formulation
     */
    void addModelingOnGroupOfCells( Physics phys, Modelings mod, std::string nameOfGroup,
                                    Formulation form = NoFormulation );

    /**
     * @brief Construction (au sens Jeveux fortran) de la sd_modele
     * @return booleen indiquant que la construction s'est bien deroulee
     */
    virtual bool build();

    /**@brief Has multi-fiber beams in model ? */
    bool existsMultiFiberBeam() const;

    /**@brief Has THM in model ? */
    bool existsThm() const;

    /**@brief Has XFEM in model ? */
    bool existsXfem() const;

    /**@brief Has HHO in model ? */
    bool existsHHO() const;

    /**@brief Has Axis element in model ? */
    bool existsAxis() const;

    /**@brief model is Axis ? */
    bool isAxis() const;

    /**@brief Has COQUE_3D in model ? */
    bool exists3DShell() const;

    /**@brief Has STRX in model ? */
    bool existsSTRX() const;

    /**@brief Has RdM in model ? */
    bool existsRdM() const;

    /**@brief Has PARTITION in model ? */
    bool existsPartition() const;

    /**@brief Name of modelisation */
    const std::string getModelisationName() const;

    /**@brief patition method */
    const std::string getPartitionMethod() const;

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
    bool xfemPreconditioningEnable() const;

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
     * @brief Get the HHO model
     */
    HHOModelPtr getHHOModel() const;

    /**
     * @brief Obtention de la methode de partition
     */
    ModelSplitingMethod getSplittingMethod() const;

    BaseMeshPtr getMesh() const;

    JeveuxVectorLong getFiniteElementType() const;

    /**
     * @brief Methode permettant de savoir si le modele est vide
     * @return true si le modele est vide
     */
    bool isEmpty() const;

    /** @brief To set if a model will be balanceable */
    void banBalancing() { _isBalanceable = false; };

    /** @brief To know if a model will be balanceable */
    bool balanceable() { return _isBalanceable; };

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

    ASTERINTEGER getXfemContact() const;

#ifdef ASTER_HAVE_MPI
    bool setFrom( const ModelPtr model );
#endif
};

/**
 * @typedef Model
 */
using ModelPtr = std::shared_ptr< Model >;

#endif /* MODEL_H_ */
