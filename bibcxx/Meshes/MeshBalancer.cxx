/**
 * @file MeshBalancer.cxx
 * @brief Implementation de MeshBalancer
 * @author Nicolas Sellenet
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

#include "Meshes/MeshBalancer.h"

#include "Meshes/Mesh.h"

void decrement( int &i ) { i--; };

ParallelMeshPtr MeshBalancer::applyBalancingStrategy( VectorInt &newLocalNodesList ) {
    ObjectBalancer nodesBalancer, cellsBalancer;
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();

    if ( _mesh != nullptr ) {
        if ( _mesh->isParallel() )
            throw std::runtime_error( "Parallel mesh not allowed" );
        if ( _mesh->isIncomplete() ) {
            VectorLong toSend( 1, _mesh->getNumberOfNodes() ), toSendAll;
            AsterMPI::all_gather( toSend, toSendAll );
            VectorLong result( nbProcs + 1, 0 );
            int pos = 0;
            for ( const auto &tmp : toSendAll ) {
                result[pos + 1] = result[pos] + tmp;
                ++pos;
            }
            _range = {result[rank], result[rank + 1]};
        } else {
            _range = {0, _mesh->getNumberOfNodes()};
        }
    }

    VectorInt newList = newLocalNodesList;
    std::for_each( newList.begin(), newList.end(), &decrement );

    VectorOfVectorsLong interfaces;
    VectorLong nOwners;
    buildBalancersAndInterfaces( newList, nodesBalancer, cellsBalancer, interfaces, nOwners );
    newList = VectorInt();

    CommGraph graph;
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        if ( iProc != rank && interfaces[2 * iProc].size() != 0 )
            graph.addCommunication( iProc );
    }
    graph.synchronizeOverProcesses();

    // interface completion with local id of opposite nodes
    VectorLong domains;
    VectorOfVectorsLong graphInterfaces;
    int tag = 0, cmpt = 0;
    for ( const auto &iProc : graph ) {
        if ( iProc != -1 ) {
            VectorLong idToSend( interfaces[2 * iProc].size(), 0. );
            VectorLong idToRecv( interfaces[2 * iProc + 1].size(), 0. );
            if ( rank > iProc ) {
                AsterMPI::send( interfaces[2 * iProc], iProc, tag );
                AsterMPI::receive( idToRecv, iProc, tag );
                ++tag;
                AsterMPI::send( interfaces[2 * iProc + 1], iProc, tag );
                AsterMPI::receive( idToSend, iProc, tag );
                ++tag;
            } else if ( rank < iProc ) {
                AsterMPI::receive( idToRecv, iProc, tag );
                AsterMPI::send( interfaces[2 * iProc], iProc, tag );
                ++tag;
                AsterMPI::receive( idToSend, iProc, tag );
                AsterMPI::send( interfaces[2 * iProc + 1], iProc, tag );
                ++tag;
            }
            domains.push_back( iProc );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc].size() ) );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc + 1].size() ) );
            for ( int i = 0; i < interfaces[2 * iProc].size(); ++i ) {
                graphInterfaces[cmpt][2 * i] = interfaces[2 * iProc][i];
                graphInterfaces[cmpt][2 * i + 1] = idToSend[i];
            }
            for ( int i = 0; i < interfaces[2 * iProc + 1].size(); ++i ) {
                graphInterfaces[cmpt + 1][2 * i] = interfaces[2 * iProc + 1][i];
                graphInterfaces[cmpt + 1][2 * i + 1] = idToRecv[i];
            }
            cmpt += 2;
        }
    }
    // free memory
    interfaces = VectorOfVectorsLong();

    ParallelMeshPtr outMesh( new ParallelMesh() );

    // Build a global numbering (if there is not)
    VectorLong nodeGlobNum;
    if ( _mesh != nullptr ) {
        nodeGlobNum.reserve( _mesh->getNumberOfNodes() );
        for ( int i = 0; i < _mesh->getNumberOfNodes(); ++i ) {
            // +1 is mandatory because connectivity starts at 1 in aster
            // cf. connex = _mesh->getConnectivity();
            nodeGlobNum.push_back( i + _range[0] + 1 );
        }
    }

    // Build mask to apply to distribute connectivity
    // Before sending, conversion to global numbering
    // After received go back to "new" local numbering
    auto dMask = ObjectBalancer::DistributedMaskOut( nodesBalancer, nodeGlobNum );

    // Build new mesh (nodes, cells types and connectivity)
    if ( _mesh == nullptr )
        _mesh = MeshPtr( new Mesh() );
    const auto coords = _mesh->getCoordinates();
    auto coordsOut = outMesh->getCoordinates();
    nodesBalancer.balanceObjectOverProcesses( coords, coordsOut );
    const auto cellsType = _mesh->getCellsType();
    auto cellsTypeOut = outMesh->getCellsType();
    cellsBalancer.balanceObjectOverProcesses( cellsType, cellsTypeOut );

    const auto connex = _mesh->getConnectivity();
    auto connexOut = outMesh->getConnectivity();
    cellsBalancer.balanceObjectOverProcesses2( connex, connexOut, dMask );

    // Build cells and nodes groups
    balanceGroups( outMesh, nodesBalancer, cellsBalancer );
    outMesh->buildInformations( 3 );

    // Build "dummy" name vectors (for cells and nodes)
    outMesh->buildNamesVectors();
    outMesh->create_joints( domains, dMask.getBalancedMask(), nOwners, graphInterfaces );
    return outMesh;
};

void MeshBalancer::buildReverseConnectivity() {
    if ( _mesh == nullptr )
        return;
    const auto connex = _mesh->getConnectivityExplorer();
    int elemId = 0;
    for ( const auto &element : connex ) {
        for ( const auto &nodeId : element ) {
            _reverseConnex[nodeId - 1].insert( elemId );
        }
        ++elemId;
    }
    _bReverseConnex = true;
};

void MeshBalancer::deleteReverseConnectivity() {
    // free memory
    _reverseConnex = std::map< int, std::set< int > >();
    _bReverseConnex = false;
};

void MeshBalancer::buildBalancersAndInterfaces( VectorInt &newLocalNodesList,
                                                ObjectBalancer &nodesB, ObjectBalancer &cellsB,
                                                VectorOfVectorsLong &interfaces,
                                                VectorLong &nOwners ) {
    // Build reverse connectivity to be able to build ObjectBalancer (what to send to which process)
    buildReverseConnectivity();

    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorOfVectorsLong procInterfaces, balanceProcInterfaces;
    VectorLong nodeOwner;
    if ( _mesh != nullptr ) {
        procInterfaces = VectorOfVectorsLong( _mesh->getNumberOfNodes() );
        nodeOwner = VectorLong( _mesh->getNumberOfNodes(), -1 );
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

            std::set< int > toAdd;
            auto returnPairToKeep =
                findNodesAndElementsInNodesNeighborhood( newLocalNodesList, toAdd );
            VectorInt toAddV, test2;
            if ( returnPairToKeep.first.size() != 0 ) {
                for ( const auto &val : toAdd ) {
                    toAddV.push_back( val );
                }
            }
            AsterMPI::all_gather( toAddV, test2 );

            std::set< int > filter;
            // In test2 ids starts at 0 in global numbering
            for ( const auto &val : test2 )
                filter.insert( val );
            // In returnPairToSend.first ids starts at 0 in local numbering
            for ( const auto &val : returnPairToKeep.first )
                filter.insert( val + _range[0] );
            VectorInt filterV;
            // So in filterV, ids starts at 0 in global numbering
            for ( const auto &val : filter )
                filterV.push_back( val );
            // And in addedNodes, ids starts at 0 in local numbering
            const auto addedNodes = findNodesToSend( filterV );
            if ( filterV.size() != 0 ) {
                for ( const auto &tmp : addedNodes ) {
                    procInterfaces[tmp].push_back( rank );
                }
                nodesB.setElementsToKeep( addedNodes );
            }
            if ( returnPairToKeep.second.size() != 0 )
                cellsB.setElementsToKeep( returnPairToKeep.second );
        } else {
            AsterMPI::bcast( nodesLists, iProc );

            std::set< int > toAdd;
            auto returnPairToSend = findNodesAndElementsInNodesNeighborhood( nodesLists, toAdd );
            VectorInt toAddV, test2;
            if ( returnPairToSend.first.size() != 0 ) {
                for ( const auto &val : toAdd ) {
                    toAddV.push_back( val );
                }
            }
            AsterMPI::all_gather( toAddV, test2 );

            std::set< int > filter;
            // In test2 ids starts at 0 in global numbering
            for ( const auto &val : test2 )
                filter.insert( val );
            // In returnPairToSend.first ids starts at 0 in local numbering
            for ( const auto &val : returnPairToSend.first )
                filter.insert( val + _range[0] );
            VectorInt filterV;
            // So in filterV, ids starts at 0 in global numbering
            for ( const auto &val : filter )
                filterV.push_back( val );
            // And in addedNodes, ids starts at 0 in local numbering
            const auto addedNodes = findNodesToSend( filterV );
            if ( filterV.size() != 0 ) {
                for ( const auto &tmp : addedNodes ) {
                    procInterfaces[tmp].push_back( iProc );
                }
                nodesB.addElementarySend( iProc, addedNodes );
            }
            if ( returnPairToSend.second.size() != 0 )
                cellsB.addElementarySend( iProc, returnPairToSend.second );
        }
        if ( iProc == rank ) {
            for ( const auto &tmp : newLocalNodesList ) {
                if ( tmp >= _range[0] && tmp < _range[1] )
                    nodeOwner[tmp - _range[0]] = rank;
            }
        } else {
            for ( const auto &tmp : nodesLists ) {
                if ( tmp >= _range[0] && tmp < _range[1] )
                    nodeOwner[tmp - _range[0]] = iProc;
            }
        }
    }
    // Save memory by destroying reverse connectivity
    deleteReverseConnectivity();
    nodesB.endElementarySendDefinition();
    cellsB.endElementarySendDefinition();
    // Prepare ObjectBalancer (graph and sizes of what to send)
    nodesB.prepareCommunications();
    cellsB.prepareCommunications();

    nodesB.balanceObjectOverProcesses2( procInterfaces, balanceProcInterfaces );
    nodesB.balanceObjectOverProcesses( nodeOwner, nOwners );
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
                                                       std::set< int > &toAddSet ) {
    std::pair< VectorInt, VectorInt > toReturn;
    auto &nodesList = toReturn.first;
    auto &elemList = toReturn.second;
    if ( _mesh == nullptr )
        return toReturn;
    if ( !_bReverseConnex )
        throw std::runtime_error( "No reverse connectivity build" );

    const auto connex = _mesh->getConnectivityExplorer();
    std::set< int > inSet;
    for ( const auto &nodeId : nodesListIn ) {
        inSet.insert( nodeId );
    }
    const auto endSet = inSet.end();

    // Find every nodes and cells in environment on nodes asks by the current process
    // Build from connectivity and reverse connectivity
    // checkedElem and checkedNodes avoid to have cells or nodes marked twice
    // !!!! WARNING : node and cell ids start at 0 !!!!
    VectorBool checkedElem( _mesh->getNumberOfCells(), false );
    VectorBool checkedNodes( _mesh->getNumberOfNodes(), false );
    const auto endPtr = _reverseConnex.end();
    for ( const auto &nodeId : nodesListIn ) {
        if ( nodeId >= _range[0] && nodeId < _range[1] ) {
            const auto idBis = nodeId - _range[0];
            if ( !checkedNodes[idBis] ) {
                nodesList.push_back( idBis );
                checkedNodes[idBis] = true;
            }
        }
        const auto &elemSetPtr = _reverseConnex.find( nodeId );
        if ( elemSetPtr == endPtr )
            continue;
        const auto elemSet = elemSetPtr->second;
        for ( const auto &elemId : elemSet ) {
            if ( checkedElem[elemId] )
                continue;
            checkedElem[elemId] = true;
            // !!!! WARNING : in connex node ids start at 1 (aster convention) !!!!
            for ( const auto &nodeId2 : connex[elemId] ) {
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
    std::sort( nodesList.begin(), nodesList.end() );
    std::sort( elemList.begin(), elemList.end() );
    return toReturn;
};

VectorInt MeshBalancer::findNodesToSend( const VectorInt &nodesListIn ) {
    VectorInt nodesList;

    for ( const auto &nodeId : nodesListIn ) {
        if ( nodeId >= _range[0] && nodeId < _range[1] ) {
            const auto idBis = nodeId - _range[0];
            nodesList.push_back( idBis );
        }
    }
    std::sort( nodesList.begin(), nodesList.end() );
    return nodesList;
};

void MeshBalancer::balanceGroups( BaseMeshPtr outMesh, const ObjectBalancer &nBalancer,
                                  const ObjectBalancer &cBalancer ) {
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

    std::map< int, std::string > mapCellsGrpNum;
    std::map< int, std::string > mapNodesGrpNum;

    VectorLong localCellGroups, localNodeGroups;
    if ( _mesh->getNumberOfCells() != 0 )
        localCellGroups = VectorLong( _mesh->getNumberOfCells(), -1 );
    if ( _mesh->getNumberOfNodes() != 0 )
        localNodeGroups = VectorLong( _mesh->getNumberOfNodes(), -1 );

    // Build a numbering of cells and nodes groups names
    // and find group number (group id) of each cells and nodes
    int cmptCells = 0;
    for ( const auto &name : toSendCell ) {
        mapCellsGrpNum[cmptCells] = name;
        for ( const auto &id : _mesh->getCells( name ) ) {
            localCellGroups[id] = cmptCells;
        }
        ++cmptCells;
    }

    int cmptNodes = 0;
    for ( const auto &name : toSendNode ) {
        mapNodesGrpNum[cmptNodes] = name;
        for ( const auto &id : _mesh->getNodes( name ) ) {
            const auto id2 = id - _range[0];
            if ( id2 >= 0 && id2 < _range[1] ) {
                localNodeGroups[id2] = cmptNodes;
            }
        }
        ++cmptNodes;
    }
    VectorLong bNodeGroups, bCellGroups;
    // "Balance" vector of group id for nodes and cells
    nBalancer.balanceObjectOverProcesses( localNodeGroups, bNodeGroups );
    cBalancer.balanceObjectOverProcesses( localCellGroups, bCellGroups );

    // Finally build groups vectors
    VectorOfVectorsLong cellsInGrp( cmptCells );
    int cellId = 1;
    for ( const auto &grpId : bCellGroups ) {
        if ( grpId != -1 ) {
            cellsInGrp[grpId].push_back( cellId );
        }
        ++cellId;
    }
    VectorString cellsGrpNames;
    VectorOfVectorsLong cellsGrpList;
    for ( int numGrp = 0; numGrp < cmptCells; ++numGrp ) {
        if ( cellsInGrp[numGrp].size() != 0 ) {
            cellsGrpNames.push_back( mapCellsGrpNum[numGrp] );
            cellsGrpList.push_back( cellsInGrp[numGrp] );
        }
    }
    // Add cells groups
    if ( cellsGrpNames.size() != 0 )
        outMesh->addGroupsOfCells( cellsGrpNames, cellsGrpList );

    VectorOfVectorsLong nodesInGrp( cmptNodes );
    int nodeId = 1;
    for ( const auto &grpId : bNodeGroups ) {
        if ( grpId != -1 )
            nodesInGrp[grpId].push_back( nodeId );
        ++nodeId;
    }
    VectorString nodesGrpNames;
    VectorOfVectorsLong nodesGrpList;
    for ( int numGrp = 0; numGrp < cmptNodes; ++numGrp ) {
        if ( nodesInGrp[numGrp].size() != 0 ) {
            nodesGrpNames.push_back( mapNodesGrpNum[numGrp] );
            nodesGrpList.push_back( nodesInGrp[numGrp] );
        }
    }
    // Add nodes groups
    if ( nodesGrpNames.size() != 0 )
        outMesh->addGroupsOfNodes( nodesGrpNames, nodesGrpList );
};
