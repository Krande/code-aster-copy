#ifndef FINITEELEMENTDESCRIPTOR_H_
#define FINITEELEMENTDESCRIPTOR_H_

/**
 * @file FiniteElementDescriptor.h
 * @brief Fichier entete de la classe FiniteElementDescriptor
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

#include "astercxx.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/MeshExplorer.h"
#include "Modeling/Model.h"

// Forward declaration
class Model;
using ModelPtr = std::shared_ptr< Model >;

/**
 * @class FiniteElementDescriptor
 * @brief Class which describes the finite elements
 * @author Nicolas Sellenet
 */
class FiniteElementDescriptor : public DataStructure {
  public:
    using ConnectivityVirtualCellsExplorer =
        MeshExplorer< CellsIteratorFromFiniteElementDescriptor, const JeveuxCollectionLong & >;

  protected:
    /** @brief Vecteur Jeveux '.NBNO' */
    JeveuxVectorLong _numberOfDelayedNumberedConstraintNodes;
    /** @brief Vecteur Jeveux '.LGRF' */
    JeveuxVectorChar8 _parameters;
    /** @brief Vecteur Jeveux '.PRNM' */
    JeveuxVectorLong _dofDescriptor;
    /** @brief Collection '.LIEL' */
    JeveuxCollectionLong _listOfGroupOfCells;
    /** @brief Vecteur Jeveux '.REPE' */
    JeveuxVectorLong _groupsOfCellsNumberByElement;
    /** @brief Collection '.NEMA' */
    JeveuxCollectionLong _delayedNumberedConstraintElementsDescriptor;
    /** @brief Vecteur Jeveux '.PRNS' */
    JeveuxVectorLong _dofOfDelayedNumberedConstraintNodes;
    /** @brief Vecteur Jeveux '.LGNS' */
    JeveuxVectorLong _virtualNodesNumbering;
    /** @brief Vecteur Jeveux '.SSSA' */
    JeveuxVectorLong _superElementsDescriptor;
    /** @brief Vecteur Jeveux '.NVGE' */
    JeveuxVectorChar16 _nameOfNeighborhoodStructure;
    /** @brief Base mesh */
    BaseMeshPtr _mesh;
    /** @brief Object to loop over connectivity of delayed numbered cells */
    const ConnectivityVirtualCellsExplorer _explorer;
    /** @brief Object to loop over list of group of cells */
    const ConnectivityVirtualCellsExplorer _explorer2;
    /** @brief Model if known */
    // We use a weak_ptr to avoid circular reference
    std::weak_ptr< Model > _model;

  public:
    /**
     * @typedef FiniteElementDescriptorPtr
     * @brief Pointeur intelligent vers un FiniteElementDescriptor
     */
    using FiniteElementDescriptorPtr = std::shared_ptr< FiniteElementDescriptor >;

    /**
     * @brief Constructeur
     */
    FiniteElementDescriptor( const std::string &name, const BaseMeshPtr mesh );

    FiniteElementDescriptor( const BaseMeshPtr mesh );

    FiniteElementDescriptor( const ModelPtr model );

    FiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc,
                             const VectorString &groupOfCells );

    FiniteElementDescriptor( const ModelPtr model, const VectorString &groupOfCells );

    /**
     * @brief Destructor
     */
    ~FiniteElementDescriptor(){};

    const ConnectivityVirtualCellsExplorer &getVirtualCellsExplorer() const;

    const JeveuxVectorLong &getVirtualNodesComponentDescriptor() const;

    const JeveuxVectorLong &getVirtualNodesNumbering() const;

    const ConnectivityVirtualCellsExplorer &getListOfGroupOfElementsExplorer() const;

    const JeveuxCollectionLong &getListOfGroupOfElements() const;

    const JeveuxCollectionLong &getNema() const;

    ASTERINTEGER getNumberOfVirtualNodes() const;

    JeveuxVectorLong getNumberOfVirtualNodesobj();

    JeveuxVectorChar8 getParameters() const;

    const JeveuxVectorLong &getPhysicalNodesComponentDescriptor() const;

    const JeveuxVectorLong &getListOfGroupOfElementsbyCell() const;

    const BaseMeshPtr getMesh() const;

    void setMesh( const BaseMeshPtr &currentMesh );

    void setModel( const ModelPtr model );

    ModelPtr getModel();

    int getPhysics( void ) const;

#ifdef ASTER_HAVE_MPI
    /** @brief Transert .PRNM from other FiniteElementDescriptor.
     * this should be associated to a ConnectionMesh,
     * other should be associated to the parallelMesh of the ConnectionMesh */
    void transferDofDescriptorFrom( FiniteElementDescriptorPtr & );

    void setFrom( FiniteElementDescriptorPtr & );

    void transferListOfGroupOfCellFrom( FiniteElementDescriptorPtr & );

#endif /* ASTER_HAVE_MPI */
};

/**
 * @typedef FiniteElementDescriptor
 * @brief Pointeur intelligent vers un FiniteElementDescriptor
 */
using FiniteElementDescriptorPtr = std::shared_ptr< FiniteElementDescriptor >;

#endif /* FINITEELEMENTDESCRIPTOR_H_ */
