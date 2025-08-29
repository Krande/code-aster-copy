#ifndef MESHBALANCER_H_
#define MESHBALANCER_H_

/**
 * @file MeshBalancer.h
 * @brief Fichier entete de la classe MeshBalancer
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "Meshes/BaseMesh.h"
#include "Meshes/ParallelMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "ParallelUtilities/ObjectBalancer.h"

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
    std::array< ASTERINTEGER, 2 > _range = { -1, -1 };
    /** @brief nodes and cells object balancers */
    ObjectBalancerPtr _nodesBalancer, _cellsBalancer;
    /** @brief cells renumbering (in order to be sorte by type) */
    VectorLong _cellRenumbering;
    /** @brief ghost layer size */
    int _ghostLayer = 1;
    /** @brief Ghost nodes on the last layer */
    VectorInt _lastLayerGhostNodes;

    void buildBalancersAndInterfaces( VectorInt &newLocalNodesList, VectorOfVectorsLong &interfaces,
                                      VectorLong &nOwners );

    void deleteReverseConnectivity();

    void balanceFamilies( BaseMeshPtr, const VectorLong & );

    void balanceGroups( BaseMeshPtr, const VectorLong & );

    /**
     * @brief Find nodes and elements in node neighborhood
     *        !!!! WARNING : return indexes are in C convention (starts at 0) !!!!
     */
    std::pair< VectorInt, VectorInt >
    findNodesAndElementsInNodesNeighborhood( const VectorInt &, VectorOfVectorsLong & );
    VectorInt findNodesToSend( const SetInt &nodesListIn );

    void sortCells( JeveuxVectorLong &typeIn, JeveuxContiguousCollectionLong &connexIn,
                    JeveuxVectorLong &typeOut, JeveuxContiguousCollectionLong &connexOut );

    void sortCells( VectorLong &vectIn, VectorLong &vectOut ) const;
    /**
     * @brief Converts the last ghosts layer from global numbering to local
     */
    void convertLastGhostLayerToLocal( const VectorLong &globNodeNumVect );

    void _enrichBalancers( const VectorInt &newLocalNodesList, int iProc, int rank,
                           VectorOfVectorsLong &procInterfaces, VectorOfVectorsLong & );

    VectorInt filterAlreadySeenNodes( const VectorInt &vec1, const SetInt &alreadySeenNodes );

  public:
    /**
     * @typedef MeshBalancerPtr
     * @brief Pointeur intelligent vers un MeshBalancer
     */
    typedef std::shared_ptr< MeshBalancer > MeshBalancerPtr;

    /**
     * @brief Constructeur
     */
    MeshBalancer()
        : _mesh( nullptr ),
          _bReverseConnex( false ),
          _nodesBalancer( new ObjectBalancer() ),
          _cellsBalancer( new ObjectBalancer() ) {};

    /**
     * @brief Apply a balancing strategy and return ParallelMeshPtr
     * @param list vector of nodes to get on local process
     * @return ParalleMesh
     */
    ParallelMeshPtr applyBalancingStrategy( const VectorInt &list,
                                            ParallelMeshPtr outMesh = nullptr,
                                            const int &ghostLayer = 1 );

    void buildFromBaseMesh( const BaseMeshPtr &mesh ) { _mesh = mesh; };

    ObjectBalancerPtr getCellObjectBalancer() const { return _cellsBalancer; };

    ObjectBalancerPtr getNodeObjectBalancer() const { return _nodesBalancer; };
};

/**
 * @typedef MeshBalancerPtr
 * @brief Pointeur intelligent vers un MeshBalancer
 */
typedef std::shared_ptr< MeshBalancer > MeshBalancerPtr;

#endif /* ASTER_HAVE_MPI */

#endif /* MESHBALANCER_H_ */
