#ifndef MESHBALANCER_H_
#define MESHBALANCER_H_

/**
 * @file MeshBalancer.h
 * @brief Fichier entete de la classe MeshBalancer
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

#include <array>

/**
 * @class MeshBalancer
 * @brief Class describing a mesh which is balanceable across MPI processes
 */
class MeshBalancer {
    /** @brief Mesh to balance */
    BaseMeshPtr _mesh;
    /** @brief Reverse connectivity (nodes to elements) !!! ids starts at 0 */
    std::map< int, std::set< int > > _reverseConnex;
    /** @brief True if _reverseConnex already build */
    bool _bReverseConnex;
    /** @brief Range of node ids of _mesh (!!! no overlaping between processes) */
    std::array< ASTERINTEGER, 2 > _range = {-1, -1};

    void buildBalancersAndInterfaces( VectorInt &newLocalNodesList, ObjectBalancer &nodesB,
                                      ObjectBalancer &cellsB, VectorOfVectorsLong &interfaces,
                                      VectorLong &nOwners );

    void buildReverseConnectivity();

    void deleteReverseConnectivity();

    void balanceGroups( BaseMeshPtr, const ObjectBalancer &, const ObjectBalancer & );

    /**
     * @brief Find nodes and elements in node neighborhood
     *        !!!! WARNING : return indexes are in C convention (starts at 0) !!!!
     */
    std::pair< VectorInt, VectorInt > findNodesAndElementsInNodesNeighborhood( const VectorInt &,
                                                                               std::set< int > & );
    VectorInt findNodesToSend( const VectorInt &nodesListIn );

  public:
    /**
     * @typedef MeshBalancerPtr
     * @brief Pointeur intelligent vers un MeshBalancer
     */
    typedef std::shared_ptr< MeshBalancer > MeshBalancerPtr;

    /**
     * @brief Constructeur
     */
    MeshBalancer() : _mesh( nullptr ), _bReverseConnex( false ) {};

    /**
     * @brief Apply a balancing strategy and return ParallelMeshPtr
     * @param list vector of nodes to get on local process
     * @return ParalleMesh
     */
    ParallelMeshPtr applyBalancingStrategy( VectorInt &list );

    void buildFromBaseMesh( const BaseMeshPtr &mesh ) { _mesh = mesh; };
};

/**
 * @typedef MeshBalancerPtr
 * @brief Pointeur intelligent vers un MeshBalancer
 */
typedef std::shared_ptr< MeshBalancer > MeshBalancerPtr;

#endif /* MESHBALANCER_H_ */
