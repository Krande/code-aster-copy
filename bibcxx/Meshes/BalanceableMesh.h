#ifndef BALANCEABLEMESH_H_
#define BALANCEABLEMESH_H_

/**
 * @file BalanceableMesh.h
 * @brief Fichier entete de la classe BalanceableMesh
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "astercxx.h"

#include "Meshes/BaseMesh.h"
#include "Meshes/ParallelMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "ParallelUtilities/ObjectBalancer.h"

/**
 * @class BalanceableMesh
 * @brief Class describing a mesh which is balanceable across MPI processes
 */
class BalanceableMesh {
    BaseMeshPtr _mesh;
    std::map< int, std::set< int > > _reverseConnex;
    bool _bReverseConnex;

    void buildBalancersAndInterfaces( VectorInt &newLocalNodesList, ObjectBalancer &nodesB,
                                      ObjectBalancer &cellsB, VectorOfVectorsLong &interfaces,
                                      VectorLong &nOwners );

    void buildReverseConnectivity();

    void deleteReverseConnectivity();

    void balanceGroups( BaseMeshPtr, const ObjectBalancer &, const ObjectBalancer & );

    VectorInt findExternalNodes( const VectorInt &, const VectorInt & );

    /**
     * @brief Find nodes and elements in node neighborhood
     *        !!!! WARNING : return indexes are in C convention (starts at 0) !!!!
     */
    std::pair< VectorInt, VectorInt > findNodesAndElementsInNodesNeighborhood( const VectorInt & );

  public:
    /**
     * @typedef BalanceableMeshPtr
     * @brief Pointeur intelligent vers un BalanceableMesh
     */
    typedef std::shared_ptr< BalanceableMesh > BalanceableMeshPtr;

    /**
     * @brief Constructeur
     */
    BalanceableMesh() : _mesh( nullptr ), _bReverseConnex( false ){};

    ParallelMeshPtr applyBalancingStrategy( VectorInt & );

    void buildFromBaseMesh( const BaseMeshPtr &mesh ) { _mesh = mesh; };
};

/**
 * @typedef BalanceableMeshPtr
 * @brief Pointeur intelligent vers un BalanceableMesh
 */
typedef std::shared_ptr< BalanceableMesh > BalanceableMeshPtr;

#endif /* BALANCEABLEMESH_H_ */
