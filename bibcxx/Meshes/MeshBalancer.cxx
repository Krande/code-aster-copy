/**
 * @file MeshBalancer.cxx
 * @brief Implementation de MeshBalancer
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Meshes/MeshBalancer.h"

#ifdef ASTER_HAVE_MPI

#include "DataFields/SimpleFieldOnNodes.h"
#include "IOManager/MedVector.h"
#include "Meshes/IncompleteMesh.h"
#include "Meshes/Mesh.h"
#include "ParallelUtilities/MeshConnectionGraph.h"

#include <limits.h>

template < typename ValType >
void decrement( ValType &i ) {
    i--;
};

struct LocalIdGlobalId {
    int localId = -1;
    int globalId = -1;
};

bool sortOnGlobalId( const LocalIdGlobalId &lhs, const LocalIdGlobalId &rhs ) {
    return lhs.globalId < rhs.globalId;
};

void buildSortedVectorToSend( const VectorLong &localIds, const VectorLong &globalNum,
                              VectorLong &sortedLocalIds ) {
    std::vector< LocalIdGlobalId > toSort( localIds.size() );
    const auto size = localIds.size();
    for ( int i = 0; i < size; ++i ) {
        toSort[i].localId = localIds[i];
        toSort[i].globalId = globalNum[localIds[i] - 1];
    }
    std::sort( toSort.begin(), toSort.end(), sortOnGlobalId );
    sortedLocalIds = VectorLong( size );
    for ( int i = 0; i < size; ++i ) {
        sortedLocalIds[i] = toSort[i].localId;
    }
}

ParallelMeshPtr MeshBalancer::applyBalancingStrategy( const VectorInt &newLocalNodesList,
                                                      ParallelMeshPtr outMesh,
                                                      const int &ghostLayer ) {
    _ghostLayer = ghostLayer;
    if ( _mesh != nullptr ) {
        if ( ghostLayer != 1 && _mesh->isParallel() ) {
            throw std::runtime_error( "Ghost layer number must be equal to 1 with ParallelMesh" );
        }
    }
    _nodesBalancer = std::make_shared< ObjectBalancer >();
    _cellsBalancer = std::make_shared< ObjectBalancer >();
    ObjectBalancer &nodesBalancer = *_nodesBalancer, cellsBalancer = *_cellsBalancer;
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();

    if ( !outMesh ) {
        outMesh = std::make_shared< ParallelMesh >();
    }

    if ( _mesh != nullptr ) {
        if ( _mesh->isIncomplete() ) {
            VectorLong toSend( 1, _mesh->getNumberOfNodes() ), toSendAll;
            AsterMPI::all_gather( toSend, toSendAll );
            VectorLong result( nbProcs + 1, 0 );
            int pos = 0;
            for ( const auto &tmp : toSendAll ) {
                result[pos + 1] = result[pos] + tmp;
                ++pos;
            }
            _range = { result[rank], result[rank + 1] };
        } else if ( _mesh->isParallel() ) {
            _range = { -1, -1 };
        } else {
            _range = { 0, _mesh->getNumberOfNodes() };
        }
    }

    VectorInt newList = newLocalNodesList;
    std::for_each( newList.begin(), newList.end(), &decrement< int > );

    VectorOfVectorsLong interfaces;
    VectorLong nOwners;
    buildBalancersAndInterfaces( newList, interfaces, nOwners );
    newList = VectorInt();

    CommGraph graph;
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        if ( iProc != rank && interfaces[2 * iProc].size() != 0 )
            graph.addCommunication( iProc );
    }
    graph.synchronizeOverProcesses();

    // Build a global numbering (if there is not)
    VectorLong nodeGlobNum;
    if ( _mesh != nullptr ) {
        if ( !_mesh->isParallel() ) {
            const auto size = _mesh->getNumberOfNodes();
            nodeGlobNum.reserve( size );
            for ( int i = 0; i < size; ++i ) {
                // +1 is mandatory because connectivity starts at 1 in aster
                // cf. connex = _mesh->getConnectivity();
                nodeGlobNum.push_back( i + _range[0] + 1 );
            }
        } else {
            const auto &numGlob = _mesh->getLocalToGlobalNodeIds();
            numGlob->updateValuePointer();
            const auto size = numGlob->size();
            nodeGlobNum.reserve( size );
            for ( int i = 0; i < size; ++i ) {
                // +1 is mandatory because connectivity starts at 1 in aster
                // cf. connex = _mesh->getConnectivity();
                nodeGlobNum.push_back( ( *numGlob )[i] + 1 );
            }
        }
    }

    // Build mask to apply to distribute connectivity
    // Before sending, conversion to global numbering
    // After received go back to "new" local numbering
    ObjectBalancer::DistributedMaskOut dMask( *_nodesBalancer, nodeGlobNum );
    const auto &globNodeNumVect = dMask.getBalancedMask();

    // interface completion with local id of opposite nodes
    VectorLong domains;
    VectorOfVectorsLong graphInterfaces;
    int cmpt = 0;
    for ( const auto &[tag, iProc] : graph ) {
        if ( iProc != -1 ) {
            VectorLong idToSend, idToRecv;
            VectorLong vec1, vec2;

            buildSortedVectorToSend( interfaces[2 * iProc], globNodeNumVect, vec1 );
            AsterMPI::send_receive( vec1, idToRecv, iProc, tag );

            buildSortedVectorToSend( interfaces[2 * iProc + 1], globNodeNumVect, vec2 );
            AsterMPI::send_receive( vec2, idToSend, iProc, tag );

            domains.push_back( iProc );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc].size() ) );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc + 1].size() ) );
            std::vector< LocalIdGlobalId > tmp( interfaces[2 * iProc].size() );
            std::vector< LocalIdGlobalId > tmp2( interfaces[2 * iProc + 1].size() );
            for ( int i = 0; i < interfaces[2 * iProc].size(); ++i ) {
                tmp[i].localId = interfaces[2 * iProc][i];
                tmp[i].globalId = globNodeNumVect[interfaces[2 * iProc][i] - 1];
            }
            // Sort joints on global id
            std::sort( tmp.begin(), tmp.end(), sortOnGlobalId );
            for ( int i = 0; i < interfaces[2 * iProc + 1].size(); ++i ) {
                tmp2[i].localId = interfaces[2 * iProc + 1][i];
                tmp2[i].globalId = globNodeNumVect[interfaces[2 * iProc + 1][i] - 1];
            }
            std::sort( tmp2.begin(), tmp2.end(), sortOnGlobalId );
            for ( int i = 0; i < interfaces[2 * iProc].size(); ++i ) {
                graphInterfaces[cmpt][2 * i] = tmp[i].localId;
                graphInterfaces[cmpt][2 * i + 1] = idToSend[i];
            }
            for ( int i = 0; i < interfaces[2 * iProc + 1].size(); ++i ) {
                graphInterfaces[cmpt + 1][2 * i] = tmp2[i].localId;
                graphInterfaces[cmpt + 1][2 * i + 1] = idToRecv[i];
            }
            cmpt += 2;
        }
    }

    // free memory
    interfaces = VectorOfVectorsLong();

    // Build new mesh (nodes, cells types and connectivity)
    if ( _mesh == nullptr )
        _mesh = std::make_shared< Mesh >();
    const auto coords = _mesh->getCoordinates();
    auto coordsOut = outMesh->getCoordinates();
    _nodesBalancer->balanceObjectOverProcesses( coords, coordsOut );
    const auto cellsType = _mesh->getCellsType();
    JeveuxVectorLong cellsTypeTmp( "TMP" );
    _cellsBalancer->balanceObjectOverProcesses( cellsType, cellsTypeTmp );

    const auto connex = _mesh->getConnectivity();
    JeveuxContiguousCollectionLong connexTmp( "TMP2" );
    if ( _mesh == nullptr ) {
        _cellsBalancer->balanceObjectOverProcesses2( connex, connexTmp, dMask );
    } else {
        if ( _mesh->isParallel() ) {
            ObjectBalancer::DistributedMask dMask2( *_nodesBalancer, nodeGlobNum );
            _cellsBalancer->balanceObjectOverProcesses2( connex, connexTmp, dMask2 );
        } else {
            _cellsBalancer->balanceObjectOverProcesses2( connex, connexTmp, dMask );
        }
    }

    JeveuxVectorLong cellsTypeOut = outMesh->getCellsType();
    JeveuxContiguousCollectionLong connexOut = outMesh->getConnectivity();
    sortCells( cellsTypeTmp, connexTmp, cellsTypeOut, connexOut );
    _cellsBalancer->setRenumbering( _cellRenumbering );

    // Build cells and nodes groups
    if ( _mesh->isIncomplete() ) {
        balanceFamilies( outMesh, _cellRenumbering );
    } else {
        balanceGroups( outMesh, _cellRenumbering );
    }
    if ( _mesh->isIncomplete() ) {
        outMesh->buildInformations( _mesh->getDimension() );
    } else {
        VectorInt dimension( 1, 0 );
        if ( rank == 0 )
            dimension[0] = _mesh->getDimension();
        AsterMPI::bcast( dimension, 0 );
        outMesh->buildInformations( dimension[0] );
    }

    auto globNodeNumVect2 = globNodeNumVect;
    std::for_each( globNodeNumVect2.begin(), globNodeNumVect2.end(), &decrement< long int > );

    convertLastGhostLayerToLocal( globNodeNumVect2 );

    // create global cell numbering
    // Build a global numbering (if there is not)
    VectorLong globCellNumVect2;
    if ( _mesh != nullptr ) {
        VectorLong cellGlobNumLoc;
        if ( _mesh->isParallel() ) {
            const auto &numGlob = _mesh->getLocalToGlobalCellIds();
            numGlob->updateValuePointer();
            const auto size = numGlob->size();
            cellGlobNumLoc.reserve( size );
            for ( int i = 0; i < size; ++i ) {
                cellGlobNumLoc.push_back( ( *numGlob )[i] );
            }
        } else {
            std::array< ASTERINTEGER, 2 > range2;
            if ( _mesh->isIncomplete() ) {
                VectorLong toSend( 1, _mesh->getNumberOfCells() ), toSendAll;
                AsterMPI::all_gather( toSend, toSendAll );
                VectorLong result( nbProcs + 1, 0 );
                int pos = 0;
                for ( const auto &tmp : toSendAll ) {
                    result[pos + 1] = result[pos] + tmp;
                    ++pos;
                }
                range2 = { result[rank], result[rank + 1] };
            } else {
                range2 = { 0, _mesh->getNumberOfCells() };
            }

            const auto size = _mesh->getNumberOfCells();
            cellGlobNumLoc.reserve( size );
            for ( int i = 0; i < size; ++i ) {
                cellGlobNumLoc.push_back( i + range2[0] );
            }
        }

        VectorLong globCellNumVect;
        _cellsBalancer->balanceObjectOverProcesses( cellGlobNumLoc, globCellNumVect );
        sortCells( globCellNumVect, globCellNumVect2 );
    }

    // Build "dummy" name vectors (for cells and nodes)
    outMesh->buildNamesVectors();
    outMesh->create_joints( domains, globNodeNumVect2, nOwners, globCellNumVect2, graphInterfaces,
                            _ghostLayer );
    outMesh->setLastGhostsLayer( _lastLayerGhostNodes );
    outMesh->updateGlobalGroupOfNodes();
    outMesh->updateGlobalGroupOfCells();
    outMesh->endDefinition();
    _cellRenumbering = VectorLong();
    return outMesh;
};

void MeshBalancer::sortCells( JeveuxVectorLong &typeIn, JeveuxContiguousCollectionLong &connexIn,
                              JeveuxVectorLong &typeOut,
                              JeveuxContiguousCollectionLong &connexOut ) {
    // TODO: Recuperer 31 a partir du nombre de types de mailles
    VectorLong nbCellByType( 31, 0 );
    const auto size = typeIn->size();
    VectorLong numCell( size, 0 );
    VectorLong numCell2( size, 0 );
    long count = 0;
    for ( const auto &type : typeIn ) {
        if ( type > 30 )
            throw std::runtime_error( "Error" );
        auto &num = nbCellByType[type + 1];
        ++num;
        numCell[count] = num;
        ++count;
    }

    for ( int i = 1; i < 30; ++i ) {
        nbCellByType[i] = nbCellByType[i] + nbCellByType[i - 1];
    }
    _cellRenumbering = VectorLong( size, 0 );
    typeOut->allocate( size );
    connexOut->allocate( size, connexIn->totalSize() );
    count = 0;
    for ( const auto &type : typeIn ) {
        long newId = numCell[count] + nbCellByType[type];
        numCell2[newId - 1] = count;
        _cellRenumbering[count] = newId;
        ( *typeOut )[newId - 1] = ( *typeIn )[count];
        ++count;
    }
    numCell = VectorLong();
    count = 1;
    for ( const auto &num : numCell2 ) {
        const auto &curCellIn = ( *connexIn )[num + 1];
        connexOut->allocateObject( count, curCellIn->size() );
        auto &curCellOut = ( *connexOut )[count];
        curCellOut->setValues( *curCellIn );
        ++count;
    }
};

void MeshBalancer::sortCells( VectorLong &vectIn, VectorLong &vectOut ) const {
    AS_ASSERT( !_cellRenumbering.empty() );
    vectOut.clear();
    vectOut = VectorLong( vectIn.size(), -1 );
    ASTERINTEGER count = 0;
    for ( const auto &num : _cellRenumbering ) {
        vectOut[num - 1] = vectIn[count];
        ++count;
    }
};

void MeshBalancer::convertLastGhostLayerToLocal( const VectorLong &globNodeNumVect ) {

    // Create mapping from value to index
    std::unordered_map< long, int > indexMap;
    for ( int i = 0; i < globNodeNumVect.size(); ++i ) {
        indexMap[globNodeNumVect[i]] = i;
    }

    // Loop  over values in _lastLayerGhostNodes
    for ( size_t i = 0; i < _lastLayerGhostNodes.size(); ++i ) {
        auto it = indexMap.find( _lastLayerGhostNodes[i] );
        if ( it != indexMap.end() ) {
            _lastLayerGhostNodes[i] = it->second;
        } else {
            std::cerr << "Warning: Value " << _lastLayerGhostNodes[i] << " not found " << std::endl;
            _lastLayerGhostNodes[i] = -1;
        }
    }
    std::sort( _lastLayerGhostNodes.begin(), _lastLayerGhostNodes.end() );
}

void MeshBalancer::deleteReverseConnectivity() {
    // free memory
    _reverseConnex = std::map< int, std::set< int > >();
    _bReverseConnex = false;
};

VectorInt MeshBalancer::filterAlreadySeenNodes( const VectorInt &vec1,
                                                const SetInt &alreadySeenNodes ) {
    bool parallelMesh = false;
    int rank = 0;
    JeveuxVectorLong meshNodeOwner;
    std::shared_ptr< const MapLong > g2LMap;
    if ( _mesh != nullptr ) {
        if ( _mesh->isParallel() ) {
            rank = getMPIRank();
            parallelMesh = true;
            meshNodeOwner = _mesh->getNodesOwner();
            g2LMap = _mesh->getGlobalToLocalNodeIds();
        }
    }

    SetInt newSet;
    if ( parallelMesh ) {
        const auto endGlob = g2LMap->end();
        for ( const auto &nodeId : vec1 ) {
            const auto iterLoc = g2LMap->find( nodeId );
            if ( iterLoc == endGlob )
                continue;
            const auto idBis = iterLoc->second;
            if ( ( *meshNodeOwner )[idBis] == rank ) {
                newSet.insert( nodeId );
            }
        }
    } else {
        for ( const auto &nodeId : vec1 ) {
            if ( alreadySeenNodes.count( nodeId ) == 0 ) {
                if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                    newSet.insert( nodeId );
                }
            }
        }
    }
    VectorInt newList;
    for ( const auto &tmp : newSet )
        newList.push_back( tmp );
    VectorInt gatheredList;
    AsterMPI::all_gather( newList, gatheredList );
    std::sort( gatheredList.begin(), gatheredList.end() );
    return gatheredList;
}

void MeshBalancer::_enrichBalancers( const VectorInt &newLocalNodesList, int iProc, int rank,
                                     VectorOfVectorsLong &procInterfaces,
                                     VectorOfVectorsLong &fastConnex ) {

    SetInt alreadySeenNodes( newLocalNodesList.begin(), newLocalNodesList.end() );
    SetInt alreadySeenCells;

    auto iterNodeList( newLocalNodesList );
    for ( int i = 1; i <= _ghostLayer; ++i ) {
        const auto returnPairToKeep =
            findNodesAndElementsInNodesNeighborhood( iterNodeList, fastConnex );

        iterNodeList = filterAlreadySeenNodes( returnPairToKeep.first, alreadySeenNodes );
        alreadySeenNodes.insert( iterNodeList.begin(), iterNodeList.end() );
        alreadySeenCells.insert( returnPairToKeep.second.begin(), returnPairToKeep.second.end() );
        int sizeCheck = iterNodeList.size();
        AsterMPI::all_reduce( sizeCheck, MPI_MIN );
        if ( sizeCheck == 0 ) {
            break;
        }
    }

    // Save last layer
    if ( iProc == rank )
        _lastLayerGhostNodes = iterNodeList;

    // And in addedNodes, ids start at 0 in local numbering
    const auto addedNodes = findNodesToSend( alreadySeenNodes );
    VectorInt addedCells( alreadySeenCells.begin(), alreadySeenCells.end() );
    std::sort( addedCells.begin(), addedCells.end() );

    if ( addedNodes.size() != 0 ) {
        for ( const auto &tmp : addedNodes ) {
            procInterfaces[tmp].push_back( iProc );
        }
        if ( iProc == rank ) {
#ifdef ASTER_DEBUG_CXX
            if ( _mesh != nullptr && _mesh->isParallel() ) {
                const auto l2G = _mesh->getLocalToGlobalNodeIds();
                std::cout << "#" << iProc
                          << " will keep (local numbering+1, global numbering+1): [";
                for ( const auto &locId : addedNodes ) {
                    std::cout << "(" << ( locId + 1 ) << ", " << ( *l2G )[locId] + 1 << "), ";
                }
                std::cout << "]" << std::endl << std::flush;
            } else {
                std::cout << "#" << iProc << " will keep (global numbering+1): [";
                for ( const auto &locId : addedNodes ) {
                    std::cout << locId + _range[0] + 1 << ", ";
                }
                std::cout << "]" << std::endl << std::flush;
            }
#endif
            _nodesBalancer->setElementsToKeep( addedNodes );
        } else {
#ifdef ASTER_DEBUG_CXX
            if ( _mesh != nullptr && _mesh->isParallel() ) {
                const auto l2G = _mesh->getLocalToGlobalNodeIds();
                std::cout << "#" << rank
                          << " will send (local numbering+1, global numbering+1) to #" << iProc
                          << ": [";
                for ( const auto &locId : addedNodes ) {
                    std::cout << "(" << ( locId + 1 ) << ", " << ( *l2G )[locId] + 1 << "), ";
                }
                std::cout << "]" << std::endl << std::flush;
            } else {
                std::cout << "#" << rank << " will send (global numbering+1) to #" << iProc
                          << ": [";
                for ( const auto &locId : addedNodes ) {
                    std::cout << locId + _range[0] + 1 << ", ";
                }
                std::cout << "]" << std::endl << std::flush;
            }
#endif
            if ( addedNodes.size() != 0 )
                _nodesBalancer->addElementarySend( iProc, addedNodes );
        }
    }
    if ( addedCells.size() != 0 ) {
        if ( iProc == rank ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << iProc << " will keep (local numbering+1) elements: [";
            for ( const auto locId : addedCells ) {
                std::cout << locId + 1 << ", ";
            }
            std::cout << "]" << std::endl << std::flush;
#endif
            _cellsBalancer->setElementsToKeep( addedCells );
        } else {
#ifdef ASTER_DEBUG_CXX
            std::cout << "#" << rank << " will send (local numbering+1) elements to #" << iProc
                      << ": [";
            for ( const auto locId : addedCells ) {
                std::cout << locId + 1 << ", ";
            }
            std::cout << "]" << std::endl << std::flush;
#endif
            _cellsBalancer->addElementarySend( iProc, addedCells );
        }
    }
    if ( _mesh != nullptr && _mesh->isParallel() ) {
        const auto &nodeOwner = _mesh->getNodesOwner();
        VectorInt toDeleteN;
        for ( auto i = 0; i < _mesh->getNumberOfNodes(); ++i ) {
            if ( ( *nodeOwner )[i] != rank )
                toDeleteN.push_back( i );
        }
        _nodesBalancer->setElementsToDelete( toDeleteN );
        const auto &cellOwner = _mesh->getCellsOwner();
        VectorInt toDeleteC;
        for ( auto i = 0; i < _mesh->getNumberOfCells(); ++i ) {
            if ( ( *cellOwner )[i] != rank )
                toDeleteC.push_back( i );
        }
        _cellsBalancer->setElementsToDelete( toDeleteC );
    }
};

void MeshBalancer::buildBalancersAndInterfaces( VectorInt &newLocalNodesList,
                                                VectorOfVectorsLong &interfaces,
                                                VectorLong &nOwners ) {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorOfVectorsLong procInterfaces, balanceProcInterfaces;
    VectorLong nodeOwner;
    if ( _mesh != nullptr ) {
        procInterfaces = VectorOfVectorsLong( _mesh->getNumberOfNodes(), VectorLong() );
        nodeOwner = VectorLong( _mesh->getNumberOfNodes(), -1 );
    }

    VectorOfVectorsLong fastConnect;
    JeveuxVectorLong meshNodeOwner;
    std::shared_ptr< const MapLong > g2LMap;
    bool parallelMesh = false;
    if ( _mesh != nullptr ) {
        if ( _mesh->isParallel() ) {
            parallelMesh = true;
            meshNodeOwner = _mesh->getNodesOwner();
            meshNodeOwner->updateValuePointer();
            g2LMap = _mesh->getGlobalToLocalNodeIds();
        }
        buildFastConnectivity( _mesh, fastConnect );
    }

    // Build ObjectBalancer by finding every nodes and cells in direct
    // environment of nodes needed by a given process
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        VectorInt size( 1, 0 );
        if ( iProc == rank )
            size[0] = newLocalNodesList.size();
        AsterMPI::bcast( size, iProc );
        VectorInt nodesLists( size[0], 0 );
        // To know what to keep and what to send, 2 phases are necessary
        // because IncompleteMesh shape
        if ( iProc == rank ) {
            AsterMPI::bcast( newLocalNodesList, iProc );
            _enrichBalancers( newLocalNodesList, iProc, rank, procInterfaces, fastConnect );
        } else {
            AsterMPI::bcast( nodesLists, iProc );
            _enrichBalancers( nodesLists, iProc, rank, procInterfaces, fastConnect );
        }
        if ( !parallelMesh ) {
            if ( iProc == rank ) {
                for ( const auto &tmp : newLocalNodesList ) {
                    if ( tmp >= _range[0] && tmp < _range[1] ) {
                        nodeOwner[tmp - _range[0]] = rank;
                    }
                }
            } else {
                for ( const auto &tmp : nodesLists ) {
                    if ( tmp >= _range[0] && tmp < _range[1] ) {
                        nodeOwner[tmp - _range[0]] = iProc;
                    }
                }
            }
        } else {
            const auto endGlob = g2LMap->end();
            if ( iProc == rank ) {
                for ( const auto &tmp : newLocalNodesList ) {
                    const auto locIter = g2LMap->find( tmp );
                    if ( locIter == endGlob )
                        continue;
                    const auto nodeIdL = ( *locIter ).second;
                    nodeOwner[nodeIdL] = rank;
                }
            } else {
                for ( const auto &tmp : nodesLists ) {
                    const auto locIter = g2LMap->find( tmp );
                    if ( locIter == endGlob )
                        continue;
                    const auto nodeIdL = ( *locIter ).second;
                    nodeOwner[nodeIdL] = iProc;
                }
            }
        }
    }
    // Save memory by destroying reverse connectivity
    if ( _mesh != nullptr )
        _mesh->deleteReverseConnectivity();

    _nodesBalancer->endElementarySendDefinition();
    _cellsBalancer->endElementarySendDefinition();
    // Prepare ObjectBalancer (graph and sizes of what to send)
    _nodesBalancer->prepareCommunications();
    _cellsBalancer->prepareCommunications();

    _nodesBalancer->balanceObjectOverProcesses2( procInterfaces, balanceProcInterfaces );
    _nodesBalancer->balanceObjectOverProcesses( nodeOwner, nOwners );

    int iNode = 0;
    interfaces = VectorOfVectorsLong( 2 * nbProcs );
    for ( const auto &vec1 : balanceProcInterfaces ) {
        const auto &ownerProc = nOwners[iNode];
        if ( ownerProc == rank ) {
            for ( const auto &proc : vec1 ) {
                if ( proc != rank )
                    interfaces[2 * proc].push_back( iNode + 1 );
            }
        } else {
            for ( const auto &proc : vec1 ) {
                if ( proc == ownerProc )
                    interfaces[2 * proc + 1].push_back( iNode + 1 );
            }
        }
        ++iNode;
    }
    for ( auto &val : nOwners ) {
        if ( val != rank )
            val = -1;
    }
};

std::pair< VectorInt, VectorInt >
MeshBalancer::findNodesAndElementsInNodesNeighborhood( const VectorInt &nodesListIn,
                                                       VectorOfVectorsLong &fastConnex ) {
    std::set< int > toAddSet;
    int rank = 0;

    std::pair< VectorInt, VectorInt > toReturn;
    VectorInt nodesList;
    auto &elemList = toReturn.second;
    const bool somethingToDo = ( _mesh == nullptr ) ? false : true;
    if ( somethingToDo ) {
        // Build reverse connectivity to be able to build ObjectBalancer
        // (what to send to which process)
        const auto &reverseConnex = _mesh->buildReverseConnectivity();

        std::set< int > inSet;
        for ( const auto &nodeId : nodesListIn ) {
            inSet.insert( nodeId );
        }
        const auto endSet = inSet.end();

        bool parallelMesh = _mesh->isParallel();
        JeveuxVectorLong meshNodeOwner, meshCellOwner, l2G;
        std::shared_ptr< const MapLong > g2LMap;
        if ( parallelMesh ) {
            rank = getMPIRank();
            meshNodeOwner = _mesh->getNodesOwner();
            meshCellOwner = _mesh->getCellsOwner();
            meshCellOwner->updateValuePointer();
            g2LMap = _mesh->getGlobalToLocalNodeIds();
            if ( meshCellOwner->size() != _mesh->getNumberOfCells() ) {
                throw std::runtime_error( "Size inconsistency" );
            }
            l2G = _mesh->getLocalToGlobalNodeIds();
        }

        // Find every nodes and cells in environment on nodes asks by the current process
        // Build from connectivity and reverse connectivity
        // checkedElem and checkedNodes avoid to have cells or nodes marked twice
        // !!!! WARNING : node and cell ids start at 0 !!!!
        VectorBool checkedElem( _mesh->getNumberOfCells(), false );
        VectorBool checkedNodes( _mesh->getNumberOfNodes(), false );
        const auto endPtr = reverseConnex.end();
        const auto endGlob = g2LMap->end();
        if ( parallelMesh ) {
            for ( const auto &nodeIdGlob : nodesListIn ) {
                const auto locIter = g2LMap->find( nodeIdGlob );
                if ( locIter == endGlob )
                    continue;
                const auto nodeIdL = ( *locIter ).second;
                if ( ( *meshNodeOwner )[nodeIdL] == rank ) {
                    if ( !checkedNodes[nodeIdL] ) {
                        nodesList.push_back( nodeIdL );
                        checkedNodes[nodeIdL] = true;
                    }
                }
                const auto &elemSetPtr = reverseConnex.find( nodeIdL );
                if ( elemSetPtr == endPtr )
                    continue;
                const auto elemSet = elemSetPtr->second;
                for ( const auto &elemId : elemSet ) {
                    if ( ( *meshCellOwner )[elemId] != rank )
                        continue;
                    if ( checkedElem[elemId] )
                        continue;
                    checkedElem[elemId] = true;
                    // !!!! WARNING : in connex node ids starts at 1 (aster convention) !!!!
                    for ( const auto &nodeIdL : fastConnex[elemId] ) {
                        const auto &nodeIdG = ( *l2G )[nodeIdL - 1];
                        if ( ( *meshNodeOwner )[nodeIdL - 1] == rank ) {
                            const auto idBis = nodeIdL - 1;
                            if ( !checkedNodes[idBis] ) {
                                nodesList.push_back( idBis );
                                checkedNodes[idBis] = true;
                            }
                        } else {
                            if ( inSet.find( nodeIdG ) == endSet ) {
                                toAddSet.insert( nodeIdG );
                            }
                        }
                    }
                    elemList.push_back( elemId );
                }
            }
        } else {
            for ( const auto &nodeId : nodesListIn ) {
                if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                    const auto idBis = nodeId - _range[0];
                    if ( !checkedNodes[idBis] ) {
                        nodesList.push_back( idBis );
                        checkedNodes[idBis] = true;
                    }
                }
                const auto &elemSetPtr = reverseConnex.find( nodeId );
                if ( elemSetPtr == endPtr )
                    continue;
                const auto elemSet = elemSetPtr->second;
                for ( const auto &elemId : elemSet ) {
                    if ( checkedElem[elemId] )
                        continue;
                    checkedElem[elemId] = true;
                    // !!!! WARNING : in connex node ids start at 1 (aster convention) !!!!
                    for ( const auto &nodeId2 : fastConnex[elemId] ) {
                        if ( nodeId2 >= _range[0] + 1 && nodeId2 < _range[1] + 1 ) {
                            const auto idBis = nodeId2 - 1 - _range[0];
                            if ( !checkedNodes[idBis] ) {
                                nodesList.push_back( idBis );
                                checkedNodes[idBis] = true;
                            }
                        } else {
                            if ( inSet.find( nodeId2 - 1 ) == endSet ) {
                                toAddSet.insert( nodeId2 - 1 );
                            }
                        }
                    }
                    elemList.push_back( elemId );
                }
            }
        }
        std::sort( nodesList.begin(), nodesList.end() );
        std::sort( elemList.begin(), elemList.end() );
    }

    VectorInt toAddV, test2;
    for ( const auto &val : toAddSet ) {
        toAddV.push_back( val );
    }
    AsterMPI::all_gather( toAddV, test2 );

    std::set< int > filter;
    // In test2 ids start at 0 in global numbering
    for ( const auto &val : test2 )
        filter.insert( val );

    // In returnPairToSend.first ids start at 0 in local numbering
    if ( _mesh != nullptr && _mesh->isParallel() ) {
        const auto l2G = _mesh->getLocalToGlobalNodeIds();
        for ( const auto &val : nodesList )
            filter.insert( ( *l2G )[val] );
    } else {
        for ( const auto &val : nodesList )
            filter.insert( val + _range[0] );
    }
    VectorInt filterV;
    // So in filterV, ids start at 0 in global numbering
    for ( const auto &val : filter )
        filterV.push_back( val );
    toReturn.first = filterV;

    return toReturn;
};

VectorInt MeshBalancer::findNodesToSend( const SetInt &nodesListIn ) {
    VectorInt nodesList;
    int rank = 0;

    bool parallelMesh = false;
    JeveuxVectorLong meshNodeOwner;
    std::shared_ptr< const MapLong > g2LMap;
    if ( _mesh != nullptr ) {
        if ( _mesh->isParallel() ) {
            rank = getMPIRank();
            parallelMesh = true;
            meshNodeOwner = _mesh->getNodesOwner();
            g2LMap = _mesh->getGlobalToLocalNodeIds();
        }
    }

    if ( parallelMesh ) {
        const auto endGlob = g2LMap->end();
        for ( const auto &nodeId : nodesListIn ) {
            const auto iterLoc = g2LMap->find( nodeId );
            if ( iterLoc == endGlob )
                continue;
            const auto idBis = iterLoc->second;
            if ( ( *meshNodeOwner )[idBis] == rank ) {
                nodesList.push_back( idBis );
            }
        }
    } else {
        for ( const auto &nodeId : nodesListIn ) {
            if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                const auto idBis = nodeId - _range[0];
                nodesList.push_back( idBis );
            }
        }
    }
    std::sort( nodesList.begin(), nodesList.end() );
    return nodesList;
};

template < class Integer >
VectorInt decodeIntegerToVector( const Integer &toDecode ) {
    VectorInt toReturn;
    constexpr size_t size = CHAR_BIT * sizeof( Integer );
    constexpr Integer toShift = 1;
    for ( int i = 0; i < size; ++i ) {
        const Integer valToTest = toShift << i;
        if ( toDecode & valToTest )
            toReturn.push_back( i );
    }
    return toReturn;
}

template < class Integer >
std::pair< int, int > integerShiftFromComponent( const int &comp ) {
    constexpr size_t size = CHAR_BIT * sizeof( Integer );
    return { comp / size, comp % size };
}

void MeshBalancer::balanceGroups( BaseMeshPtr outMesh, const VectorLong &cellRenumbering ) {
    VectorString toSendCell, toSendNode, toSendCellAll, toSendNodeAll;
    if ( _mesh != nullptr ) {
        toSendCell = _mesh->getGroupsOfCells();
        toSendNode = _mesh->getGroupsOfNodes();
    }

    // Build vectors of all cells and nodes groups names
    AsterMPI::all_gather( toSendCell, toSendCellAll );
    AsterMPI::all_gather( toSendNode, toSendNodeAll );
    std::set< std::string > checkCellGrp, checkNodeGrp;
    for ( const auto &name : toSendCellAll )
        checkCellGrp.insert( name );
    for ( const auto &name : toSendNodeAll )
        checkNodeGrp.insert( name );
    toSendCell = VectorString();
    toSendNode = VectorString();
    for ( const auto &name : checkCellGrp )
        toSendCell.push_back( name );
    for ( const auto &name : checkNodeGrp )
        toSendNode.push_back( name );
    std::sort( toSendCell.begin(), toSendCell.end() );
    std::sort( toSendNode.begin(), toSendNode.end() );

    const auto nbNodeGroups = toSendNode.size();
    const auto nbCellGroups = toSendCell.size();
    const auto pairNode = integerShiftFromComponent< long int >( nbNodeGroups );
    const auto pairCell = integerShiftFromComponent< long int >( nbCellGroups );

    const int nbECNode = ( pairNode.second == 0 ? pairNode.first : pairNode.first + 1 );
    const int nbECCell = ( pairCell.second == 0 ? pairCell.first : pairCell.first + 1 );

    std::map< int, std::string > mapCellsGrpNum;
    std::map< int, std::string > mapNodesGrpNum;

    MedVectorLongPtr localCellGroups0( new MedVector< long int >( _mesh->getNumberOfCells() ) );
    localCellGroups0->setComponentNumber( nbECCell );
    localCellGroups0->setSize( _mesh->getNumberOfCells() );
    localCellGroups0->endDefinition();
    MedVectorLongPtr localNodeGroups0( new MedVector< long int >( _mesh->getNumberOfNodes() ) );
    localNodeGroups0->setComponentNumber( nbECNode );
    localNodeGroups0->setSize( _mesh->getNumberOfNodes() );
    localNodeGroups0->endDefinition();
    VectorLong localCellGroups, localNodeGroups;
    if ( _mesh->getNumberOfCells() != 0 )
        localCellGroups = VectorLong( _mesh->getNumberOfCells(), -1 );
    if ( _mesh->getNumberOfNodes() != 0 )
        localNodeGroups = VectorLong( _mesh->getNumberOfNodes(), -1 );

    bool parallelMesh = false;
    const auto rank = getMPIRank();
    JeveuxVectorLong meshNodeOwner, meshCellOwner;
    std::shared_ptr< const MapLong > g2LMap;
    if ( _mesh != nullptr ) {
        if ( _mesh->isParallel() ) {
            parallelMesh = true;
            meshNodeOwner = _mesh->getNodesOwner();
            meshCellOwner = _mesh->getCellsOwner();
            g2LMap = _mesh->getGlobalToLocalNodeIds();
        }
    }

    // Build a numbering of cells and nodes groups names
    // and find group number (group id) of each cells and nodes
    int cmptCells = 0;
    for ( const auto &name : toSendCell ) {
        mapCellsGrpNum[cmptCells] = name;
        const auto &pair = integerShiftFromComponent< long int >( cmptCells );
        for ( const auto &id : _mesh->getCells( name ) ) {
            ( *localCellGroups0 )[id][pair.first] += 1UL << pair.second;
        }
        ++cmptCells;
    }

    int cmptNodes = 0;
    if ( parallelMesh ) {
        const auto &endIter = g2LMap->end();
        for ( const auto &name : toSendNode ) {
            mapNodesGrpNum[cmptNodes] = name;
            const auto &pair = integerShiftFromComponent< long int >( cmptNodes );
            for ( const auto &id : _mesh->getNodes( name ) ) {
                if ( ( *meshNodeOwner )[id] == rank ) {
                    ( *localNodeGroups0 )[id][pair.first] += 1UL << pair.second;
                }
            }
            ++cmptNodes;
        }
    } else {
        for ( const auto &name : toSendNode ) {
            mapNodesGrpNum[cmptNodes] = name;
            const auto &pair = integerShiftFromComponent< long int >( cmptNodes );
            for ( const auto &id : _mesh->getNodes( name ) ) {
                const auto id2 = id - _range[0];
                if ( id2 >= 0 && id2 < _range[1] ) {
                    ( *localNodeGroups0 )[id][pair.first] += 1UL << pair.second;
                }
            }
            ++cmptNodes;
        }
    }
    auto outN = _nodesBalancer->balanceMedVectorOverProcessesWithRenumbering< long int >(
        localNodeGroups0 );
    auto outC = _cellsBalancer->balanceMedVectorOverProcessesWithRenumbering< long int >(
        localCellGroups0 );

    constexpr size_t sizeLong = CHAR_BIT * sizeof( long int );

    const int nbNode = outN->size();
    VectorOfVectorsLong nodesInGrp( cmptNodes );
    for ( int i = 0; i < nbNode; ++i ) {
        const auto &curElem = ( *outN )[i];
        const int nbCmp = curElem.getComponentNumber();
        for ( int j = 0; j < nbCmp; ++j ) {
            const auto &curValue = curElem[j];
            if ( curValue != 0 ) {
                const auto decode = decodeIntegerToVector< long int >( curValue );
                for ( const auto &val : decode ) {
                    nodesInGrp[j * sizeLong + val].push_back( i + 1 );
                }
            }
        }
    }
    const int nbCell = outC->size();
    VectorOfVectorsLong cellsInGrp( cmptCells );
    for ( int i = 0; i < nbCell; ++i ) {
        const auto &curElem = ( *outC )[i];
        const int nbCmp = curElem.getComponentNumber();
        for ( int j = 0; j < nbCmp; ++j ) {
            const auto &curValue = curElem[j];
            if ( curValue != 0 ) {
                const auto decode = decodeIntegerToVector< long int >( curValue );
                for ( const auto &val : decode ) {
                    cellsInGrp[j * sizeLong + val].push_back( i + 1 );
                }
            }
        }
    }

    VectorString cellsGrpNames;
    VectorOfVectorsLong cellsGrpList;
    for ( int numGrp = 0; numGrp < cmptCells; ++numGrp ) {
        if ( cellsInGrp[numGrp].size() != 0 ) {
            cellsGrpNames.push_back( mapCellsGrpNum[numGrp] );
            cellsGrpList.push_back( cellsInGrp[numGrp] );
        }
    }
    // Add cell groups
    if ( cellsGrpNames.size() != 0 )
        outMesh->addGroupsOfCells( cellsGrpNames, cellsGrpList );

    VectorString nodesGrpNames;
    VectorOfVectorsLong nodesGrpList;
    for ( int numGrp = 0; numGrp < cmptNodes; ++numGrp ) {
        if ( nodesInGrp[numGrp].size() != 0 ) {
            nodesGrpNames.push_back( mapNodesGrpNum[numGrp] );
            nodesGrpList.push_back( nodesInGrp[numGrp] );
        }
    }
    // Add node groups
    if ( nodesGrpNames.size() != 0 )
        outMesh->addGroupsOfNodes( nodesGrpNames, nodesGrpList );
};

void MeshBalancer::balanceFamilies( BaseMeshPtr mesh, const VectorLong &renumber ) {
    const auto &nodeFam = _mesh->getNodeFamily();
    VectorLong nodeFamB;
    _nodesBalancer->balanceObjectOverProcesses( nodeFam, nodeFamB );
    const auto &cellFam = _mesh->getCellFamily();
    VectorLong cellFamB;
    _cellsBalancer->balanceObjectOverProcesses( cellFam, cellFamB );
    VectorLong cellFamBR( cellFamB.size(), 0 );
    for ( int i = 0; i < cellFamB.size(); ++i ) {
        auto newId = renumber[i];
        cellFamBR[newId - 1] = cellFamB[i];
    }
    const auto &nodeFamilyGroups = _mesh->getNodeFamilyGroups();
    VectorOfVectorsLong idGroups( nodeFamilyGroups.size() );
    std::map< std::string, int > mapGrpInt;
    int count1 = 0, count2 = 0;
    VectorString nodesGrpNames0;
    for ( const auto &grps : nodeFamilyGroups ) {
        for ( const auto grp : grps ) {
            if ( mapGrpInt.count( grp ) == 0 ) {
                idGroups[count1].push_back( count2 );
                mapGrpInt[grp] = count2;
                nodesGrpNames0.push_back( grp );
                ++count2;
            } else {
                idGroups[count1].push_back( mapGrpInt[grp] );
            }
        }
        ++count1;
    }
    VectorOfVectorsLong nodesGrpList0( nodesGrpNames0.size() );
    count1 = 0;
    for ( const auto &numFam : nodeFamB ) {
        ++count1;
        if ( numFam == 0 )
            continue;
        const auto vecIdGrps = idGroups[numFam - 1];
        for ( const auto &idGrp : vecIdGrps ) {
            nodesGrpList0[idGrp].push_back( count1 );
        }
    }
    VectorString nodesGrpNames;
    VectorOfVectorsLong nodesGrpList;
    for ( int i = 0; i < nodesGrpNames0.size(); ++i ) {
        if ( nodesGrpList0[i].size() != 0 ) {
            nodesGrpNames.push_back( nodesGrpNames0[i] );
            nodesGrpList.push_back( nodesGrpList0[i] );
        }
    }
    // Add node groups
    if ( nodesGrpNames.size() != 0 )
        mesh->addGroupsOfNodes( nodesGrpNames, nodesGrpList );

    const auto &cellFamilyGroups = _mesh->getCellFamilyGroups();
    idGroups = VectorOfVectorsLong( cellFamilyGroups.size() );
    mapGrpInt = std::map< std::string, int >();
    count1 = 0;
    count2 = 0;
    VectorString cellsGrpNames0;
    for ( const auto &grps : cellFamilyGroups ) {
        for ( const auto grp : grps ) {
            if ( mapGrpInt.count( grp ) == 0 ) {
                idGroups[count1].push_back( count2 );
                mapGrpInt[grp] = count2;
                cellsGrpNames0.push_back( grp );
                ++count2;
            } else {
                idGroups[count1].push_back( mapGrpInt[grp] );
            }
        }
        ++count1;
    }
    VectorOfVectorsLong cellsGrpList0( cellsGrpNames0.size() );
    count1 = 0;
    for ( const auto &numFam : cellFamBR ) {
        ++count1;
        if ( numFam == 0 )
            continue;
        const auto vecIdGrps = idGroups[-numFam - 1];
        for ( const auto &idGrp : vecIdGrps ) {
            cellsGrpList0[idGrp].push_back( count1 );
        }
    }
    VectorString cellsGrpNames;
    VectorOfVectorsLong cellsGrpList;
    for ( int i = 0; i < cellsGrpNames0.size(); ++i ) {
        if ( cellsGrpList0[i].size() != 0 ) {
            cellsGrpNames.push_back( cellsGrpNames0[i] );
            cellsGrpList.push_back( cellsGrpList0[i] );
        }
    }
    // Add cell groups
    if ( cellsGrpNames.size() != 0 )
        mesh->addGroupsOfCells( cellsGrpNames, cellsGrpList );
};

#endif /* ASTER_HAVE_MPI */
