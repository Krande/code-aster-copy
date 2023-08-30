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

#ifdef ASTER_HAVE_MPI

#include "Meshes/IncompleteMesh.h"
#include "Meshes/Mesh.h"
#include "ParallelUtilities/MeshConnectionGraph.h"

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

ParallelMeshPtr MeshBalancer::applyBalancingStrategy( VectorInt &newLocalNodesList,
                                                      ParallelMeshPtr outMesh ) {
    _nodesBalancer = std::make_shared< ObjectBalancer >();
    _cellsBalancer = std::make_shared< ObjectBalancer >();
    ObjectBalancer &nodesBalancer = *_nodesBalancer, cellsBalancer = *_cellsBalancer;
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();

    if ( !outMesh ) {
        outMesh = std::make_shared< ParallelMesh >();
    }

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
        const auto size = _mesh->getNumberOfNodes();
        nodeGlobNum.reserve( size );
        for ( int i = 0; i < size; ++i ) {
            // +1 is mandatory because connectivity starts at 1 in aster
            // cf. connex = _mesh->getConnectivity();
            nodeGlobNum.push_back( i + _range[0] + 1 );
        }
    }

    // Build mask to apply to distribute connectivity
    // Before sending, conversion to global numbering
    // After received go back to "new" local numbering
    auto dMask = ObjectBalancer::DistributedMaskOut( *_nodesBalancer, nodeGlobNum );
    const auto &globNumVect = dMask.getBalancedMask();

    // interface completion with local id of opposite nodes
    VectorLong domains;
    VectorOfVectorsLong graphInterfaces;
    int cmpt = 0;
    for ( const auto &[tag, iProc] : graph ) {
        if ( iProc != -1 ) {
            VectorLong idToSend, idToRecv;
            VectorLong vec1, vec2;

            buildSortedVectorToSend( interfaces[2 * iProc], globNumVect, vec1 );
            AsterMPI::send_receive( vec1, idToRecv, iProc, tag );

            buildSortedVectorToSend( interfaces[2 * iProc + 1], globNumVect, vec2 );
            AsterMPI::send_receive( vec2, idToSend, iProc, tag );

            domains.push_back( iProc );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc].size() ) );
            graphInterfaces.push_back( VectorLong( 2 * interfaces[2 * iProc + 1].size() ) );
            std::vector< LocalIdGlobalId > tmp( interfaces[2 * iProc].size() );
            std::vector< LocalIdGlobalId > tmp2( interfaces[2 * iProc + 1].size() );
            for ( int i = 0; i < interfaces[2 * iProc].size(); ++i ) {
                tmp[i].localId = interfaces[2 * iProc][i];
                tmp[i].globalId = globNumVect[interfaces[2 * iProc][i] - 1];
            }
            // Sort joints on global id
            std::sort( tmp.begin(), tmp.end(), sortOnGlobalId );
            for ( int i = 0; i < interfaces[2 * iProc + 1].size(); ++i ) {
                tmp2[i].localId = interfaces[2 * iProc + 1][i];
                tmp2[i].globalId = globNumVect[interfaces[2 * iProc + 1][i] - 1];
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
    JeveuxCollectionLong connexTmp( "TMP2" );
    _cellsBalancer->balanceObjectOverProcesses2( connex, connexTmp, dMask );

    JeveuxVectorLong cellsTypeOut = outMesh->getCellsType();
    JeveuxCollectionLong connexOut = outMesh->getConnectivity();
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

    auto globNumVect2 = globNumVect;
    std::for_each( globNumVect2.begin(), globNumVect2.end(), &decrement< long int > );

    // Build "dummy" name vectors (for cells and nodes)
    outMesh->buildNamesVectors();
    outMesh->create_joints( domains, globNumVect2, nOwners, graphInterfaces );
    outMesh->updateGlobalGroupOfNodes();
    outMesh->updateGlobalGroupOfCells();
    outMesh->endDefinition();
    _cellRenumbering = VectorLong();
    return outMesh;
};

void MeshBalancer::sortCells( JeveuxVectorLong &typeIn, JeveuxCollectionLong &connexIn,
                              JeveuxVectorLong &typeOut, JeveuxCollectionLong &connexOut ) {
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
    connexOut->allocateContiguousNumbered( size, connexIn->totalSize() );
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

void MeshBalancer::deleteReverseConnectivity() {
    // free memory
    _reverseConnex = std::map< int, std::set< int > >();
    _bReverseConnex = false;
};

void MeshBalancer::_enrichBalancers( VectorInt &newLocalNodesList, int iProc, int rank,
                                     VectorOfVectorsLong &procInterfaces,
                                     VectorOfVectorsLong &fastConnex ) {
    AsterMPI::bcast( newLocalNodesList, iProc );

    std::set< int > toAdd;
    auto returnPairToKeep =
        findNodesAndElementsInNodesNeighborhood( newLocalNodesList, toAdd, fastConnex );
    VectorInt toAddV, test2;
    for ( const auto &val : toAdd ) {
        toAddV.push_back( val );
    }
    AsterMPI::all_gather( toAddV, test2 );

    std::set< int > filter;
    // In test2 ids start at 0 in global numbering
    for ( const auto &val : test2 )
        filter.insert( val );
    // In returnPairToSend.first ids start at 0 in local numbering
    for ( const auto &val : returnPairToKeep.first )
        filter.insert( val + _range[0] );
    VectorInt filterV;
    // So in filterV, ids start at 0 in global numbering
    for ( const auto &val : filter )
        filterV.push_back( val );
    // And in addedNodes, ids start at 0 in local numbering
    const auto addedNodes = findNodesToSend( filterV );
    if ( filterV.size() != 0 ) {
        for ( const auto &tmp : addedNodes ) {
            procInterfaces[tmp].push_back( iProc );
        }
        if ( iProc == rank ) {
            _nodesBalancer->setElementsToKeep( addedNodes );
        } else {
            if ( addedNodes.size() != 0 )
                _nodesBalancer->addElementarySend( iProc, addedNodes );
        }
    }
    if ( returnPairToKeep.second.size() != 0 ) {
        if ( iProc == rank ) {
            _cellsBalancer->setElementsToKeep( returnPairToKeep.second );
        } else {
            _cellsBalancer->addElementarySend( iProc, returnPairToKeep.second );
        }
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
        procInterfaces = VectorOfVectorsLong( _mesh->getNumberOfNodes() );
        nodeOwner = VectorLong( _mesh->getNumberOfNodes(), -1 );
    }

    VectorOfVectorsLong fastConnect;
    if ( _mesh != nullptr ) {
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
            _enrichBalancers( newLocalNodesList, iProc, rank, procInterfaces, fastConnect );
        } else {
            _enrichBalancers( nodesLists, iProc, rank, procInterfaces, fastConnect );
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

std::pair< VectorInt, VectorInt > MeshBalancer::findNodesAndElementsInNodesNeighborhood(
    const VectorInt &nodesListIn, std::set< int > &toAddSet, VectorOfVectorsLong &fastConnex ) {

    std::pair< VectorInt, VectorInt > toReturn;
    auto &nodesList = toReturn.first;
    auto &elemList = toReturn.second;
    if ( _mesh == nullptr )
        return toReturn;
    // Build reverse connectivity to be able to build ObjectBalancer (what to send to which process)
    const auto &reverseConnex = _mesh->buildReverseConnectivity();

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
    const auto endPtr = reverseConnex.end();
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
    _nodesBalancer->balanceObjectOverProcesses( localNodeGroups, bNodeGroups );
    _cellsBalancer->balanceObjectOverProcesses( localCellGroups, bCellGroups );

    // Finally build groups vectors
    VectorOfVectorsLong cellsInGrp( cmptCells );
    int cellId = 1;
    for ( const auto &grpId : bCellGroups ) {
        if ( grpId != -1 ) {
            auto numCell = cellRenumbering[cellId - 1];
            cellsInGrp[grpId].push_back( numCell );
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
    // Add cell groups
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
