/**
 * @file BalanceableMesh.cxx
 * @brief Implementation de BalanceableMesh
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

#include "Meshes/BalanceableMesh.h"

#include "Meshes/Mesh.h"
#include "Meshes/ParallelMesh.h"

void BalanceableMesh::applyBalancingStrategy( VectorInt &newLocalNodesList ) {
    ObjectBalancer nodesBalancer, cellsBalancer;

    buildBalancers( newLocalNodesList, nodesBalancer, cellsBalancer );

    ParallelMeshPtr outMesh( new ParallelMesh() );

    const auto rank = getMPIRank();

    // Build a global numbering (if there is not)
    VectorLong globNum;
    if ( _mesh != nullptr ) {
        if ( rank != 0 )
            throw std::runtime_error( "Sequential mesh must be defined only on #0" );
        globNum.reserve( _mesh->getNumberOfNodes() );
        for ( int i = 0; i < _mesh->getNumberOfNodes(); ++i ) {
            globNum.push_back( i + 1 );
        }
    }
    nodesBalancer.endElementarySendDefinition();
    cellsBalancer.endElementarySendDefinition();
    // Prepare ObjectBalancer (graph and sizes of what to send)
    nodesBalancer.prepareCommunications();
    cellsBalancer.prepareCommunications();

    // Build mask to apply to distribute connectivity
    // Before sending, conversion to global numbering
    // After received go back to "new" local numbering
    auto dMask = ObjectBalancer::DistributedMask( nodesBalancer, globNum );

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
    // JeveuxCollectionLong test( "RIEN" );
    // cellsBalancer.testBalance( connex, test, dMask );
    // cellsBalancer.testBalance( connex, connexOut, dMask );
    // std::cout << "Validation " << *connexOut << std::endl;
    // std::cout << "Validation2 " << *test << std::endl;

    // Build cells and nodes groups
    balanceGroups( outMesh, nodesBalancer, cellsBalancer );
    outMesh->buildInformations( 3 );

    // Build "dummy" names vectors (for cells and nodes)
    outMesh->buildNamesVectors();
    outMesh->printMedFile( "/home/H85256/" + std::to_string( rank ) + ".med" );
    // outMesh->debugPrint();
};

void BalanceableMesh::buildReverseConnectivity() {
    if ( _mesh == nullptr )
        return;
    const auto connex = _mesh->getConnectivityExplorer();
    int elemId = 1;
    for ( const auto &element : connex ) {
        for ( const auto &nodeId : element ) {
            _reverseConnex[nodeId - 1].insert( elemId );
        }
        ++elemId;
    }
    _bReverseConnex = true;
};

void BalanceableMesh::deleteReverseConnectivity() {
    _reverseConnex.clear();
    _bReverseConnex = false;
};

void BalanceableMesh::buildBalancers( VectorInt &newLocalNodesList, ObjectBalancer &nodesB,
                                      ObjectBalancer &cellsB ) {
    // Build reverse connectivity to be able to build ObjectBalancer (what to send to which process)
    buildReverseConnectivity();

    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorOfVectorsLong test, test2;
    // Build ObjectBalancer by finding every nodes and cells in direct
    // environment of nodes needed by a given process
    for ( int iProc = 0; iProc < nbProcs; ++iProc ) {
        VectorInt size( 1, 0 );
        if ( iProc == rank )
            size[0] = newLocalNodesList.size();
        AsterMPI::bcast( size, iProc );
        VectorInt nodesLists( size[0], 0 );
        if ( iProc == rank ) {
            AsterMPI::bcast( newLocalNodesList, iProc );
            auto returnPairToKeep = findNodesAndElementsInNodesNeighborhood( newLocalNodesList );
            if ( returnPairToKeep.first.size() != 0 ) {
                nodesB.setElementsToKeep( returnPairToKeep.first );
                const auto extNodes =
                    findExternalNodes( newLocalNodesList, returnPairToKeep.first );
            }
            if ( returnPairToKeep.second.size() != 0 )
                cellsB.setElementsToKeep( returnPairToKeep.second );
        } else {
            AsterMPI::bcast( nodesLists, iProc );
            auto returnPairToSend = findNodesAndElementsInNodesNeighborhood( nodesLists );
            if ( returnPairToSend.first.size() != 0 ) {
                nodesB.addElementarySend( iProc, returnPairToSend.first );
                const auto extNodes =
                    findExternalNodes( newLocalNodesList, returnPairToSend.first );
            }
            if ( returnPairToSend.second.size() != 0 )
                cellsB.addElementarySend( iProc, returnPairToSend.second );
        }
    }
    // Save memory by destroying reverse connectivity
    deleteReverseConnectivity();
    // cellsBalancer.balanceObjectOverProcesses2( test, test2, dMask );
};

VectorInt BalanceableMesh::findExternalNodes( const VectorInt &askedNodes,
                                              const VectorInt &sendNodes ) {
    VectorInt toReturn;
    return toReturn;
};

std::pair< VectorInt, VectorInt >
BalanceableMesh::findNodesAndElementsInNodesNeighborhood( const VectorInt &nodesListIn ) {
    std::pair< VectorInt, VectorInt > toReturn;
    auto &nodesList = toReturn.first;
    auto &elemList = toReturn.second;
    if ( _mesh == nullptr )
        return toReturn;
    if ( !_bReverseConnex )
        throw std::runtime_error( "No reverse connectivity build" );

    const auto connex = _mesh->getConnectivityExplorer();

    // Find every nodes and cells in environment on nodes asks by the current process
    // Build from connectivity and reverse connectivity
    // checkedElem and checkedNodes avoid to have cells or nodes marked twice
    // !!!! WARNING : node and cell ids start at 0 !!!!
    VectorBool checkedElem( _mesh->getNumberOfCells(), false );
    VectorBool checkedNodes( _mesh->getNumberOfNodes(), false );
    for ( const auto &nodeId : nodesListIn ) {
        const auto elemSet = _reverseConnex[nodeId - 1];
        for ( const auto &elemId : elemSet ) {
            if ( checkedElem[elemId - 1] )
                continue;
            checkedElem[elemId - 1] = true;
            for ( const auto &nodeId2 : connex[elemId - 1] ) {
                if ( !checkedNodes[nodeId2 - 1] ) {
                    nodesList.push_back( nodeId2 - 1 );
                    checkedNodes[nodeId2 - 1] = true;
                }
            }
            elemList.push_back( elemId - 1 );
        }
    }
    return toReturn;
};

void BalanceableMesh::balanceGroups( BaseMeshPtr outMesh, const ObjectBalancer &nBalancer,
                                     const ObjectBalancer &cBalancer ) {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();

    VectorString toSendCell, toSendNode, toSendCellAll, toSendNodeAll;
    if ( _mesh != nullptr ) {
        toSendCell = _mesh->getGroupsOfCells();
        toSendNode = _mesh->getGroupsOfNodes();
    }

    // Build vectors of all cells and nodes groups names
    AsterMPI::all_gather( toSendCell, toSendCellAll );
    AsterMPI::all_gather( toSendNode, toSendNodeAll );

    std::map< int, std::string > mapCellsGrpNum;
    std::map< int, std::string > mapNodesGrpNum;

    VectorLong bCellsGroups, bNodesGroups;
    if ( _mesh->getNumberOfCells() != 0 )
        bCellsGroups = VectorLong( _mesh->getNumberOfCells(), -1 );
    if ( _mesh->getNumberOfNodes() != 0 )
        bNodesGroups = VectorLong( _mesh->getNumberOfNodes(), -1 );

    // Build a numbering of cells and nodes groups names
    // and find group number (group id) of each cells and nodes
    int cmptCells = 0;
    for ( const auto &name : toSendCellAll ) {
        mapCellsGrpNum[cmptCells] = name;
        for ( const auto &id : _mesh->getCells( name ) ) {
            bCellsGroups[id - 1] = cmptCells;
        }
        ++cmptCells;
    }
    int cmptNodes = 0;
    for ( const auto &name : toSendNodeAll ) {
        mapNodesGrpNum[cmptNodes] = name;
        for ( const auto &id : _mesh->getNodes( name ) ) {
            bNodesGroups[id - 1] = cmptNodes;
        }
        ++cmptNodes;
    }
    VectorLong out1, out2;
    // "Balance" vector of group id for nodes and cells
    nBalancer.balanceObjectOverProcesses( bNodesGroups, out1 );
    cBalancer.balanceObjectOverProcesses( bCellsGroups, out2 );

    // Finally build groups vectors
    VectorOfVectorsLong cellsInGrp( cmptCells );
    int test = 1;
    for ( const auto &grpId : out2 ) {
        if ( grpId != -1 )
            cellsInGrp[grpId].push_back( test );
        ++test;
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
    test = 1;
    for ( const auto &grpId : out1 ) {
        if ( grpId != -1 )
            nodesInGrp[grpId].push_back( test );
        ++test;
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
