/**
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

#include "ParallelUtilities/PtScotchPartitioner.h"

#ifdef ASTER_HAVE_PTSCOTCH

#include "ParallelUtilities/ObjectBalancer.h"

PtScotchPartitioner::PtScotchPartitioner() : _graph2( nullptr ), _context( nullptr ) {
    _graph = new SCOTCH_Dgraph;
    _scotchStrat = new SCOTCH_Strat;
    SCOTCH_dgraphInit( _graph, aster_get_current_comm()->id );
    SCOTCH_stratInit( _scotchStrat );
};

PtScotchPartitioner::~PtScotchPartitioner() {
    SCOTCH_dgraphFree( _graph );
    SCOTCH_stratFree( _scotchStrat );
    delete _graph;
    delete _scotchStrat;
};

int PtScotchPartitioner::buildGraph( const VectorLong &vertices, const VectorLong &edges,
                                     const VectorLong &weights ) {
    _nbVertex = vertices.size() - 1;

    _vertices = VectorScotchNum( begin( vertices ), end( vertices ) );
    _edges = VectorScotchNum( begin( edges ), end( edges ) );
    if ( weights.size() == 0 ) {
        return SCOTCH_dgraphBuild( _graph, 0, _vertices.size() - 1, _vertices.size() - 1,
                                   _vertices.data(), 0, 0, 0, _edges.size(), _edges.size(),
                                   _edges.data(), 0, 0 );
    } else {
        _weights = VectorScotchNum( begin( weights ), end( weights ) );
        return SCOTCH_dgraphBuild( _graph, 0, _vertices.size() - 1, _vertices.size() - 1,
                                   _vertices.data(), 0, _weights.data(), 0, _edges.size(),
                                   _edges.size(), _edges.data(), 0, 0 );
    }
};

int PtScotchPartitioner::buildGraph( const MeshConnectionGraphPtr &graph,
                                     const VectorOfVectorsLong &nodesToGather ) {
    _minId = graph->getRange()[0];
    if ( nodesToGather.size() == 0 ) {
        return buildGraph( graph->getVertices(), graph->getEdges(), graph->getVertexWeights() );
    } else {
        // If the user asks to gather some nodes
        auto &vert = graph->getVertices();
        auto &edge = graph->getEdges();
        _nbVertex = vert.size() - 1;
        _gatheredNodes = true;
        const auto nbProcs = getMPISize();
        VectorOfVectorsLong firstNodesAndWeights =
            VectorOfVectorsLong( nodesToGather.size(), VectorLong( 2, 0 ) );
        VectorOfVectorsLong nodesToConnectToMaster;
        _range = graph->getRange();

        // First put a flag on nodes to gather (forget) in nodesToForget
        // and save weight to apply on survival node (firstNodesAndWeights)
        int pos = 0;
        std::map< int, int > nodesToForget, masterNodeToPos;
        for ( const auto &nodeVec : nodesToGather ) {
            VectorLong allNodeId, nodeVec2;
            // nodeVec to global numbering
            for ( const auto &nodeId : nodeVec ) {
                nodeVec2.push_back( nodeId + _range[0] );
            }
            // Share node ids to gather over procs
            AsterMPI::all_gather( nodeVec2, allNodeId );
            std::sort( allNodeId.begin(), allNodeId.end() );

            // Save this vector
            _toForgetVector.push_back( allNodeId );
            const auto &savedNodeId = allNodeId[0];
            auto &fNAW = firstNodesAndWeights[pos];
            fNAW[0] = savedNodeId;
            fNAW[1] = allNodeId.size();
            masterNodeToPos[savedNodeId] = pos;
            // All nodes but the first must put in nodesToForget
            // First node will be the master node of group
            bool first = true;
            for ( const auto &curId : allNodeId ) {
                if ( first ) {
                    first = false;
                    continue;
                }
                const auto toFind = nodesToForget.find( curId );
                if ( toFind != nodesToForget.end() ) {
                    throw std::runtime_error( "Nodes to gather must be in only one list."
                                              " Tip: fuse lists" );
                } else {
                    nodesToForget[curId] = savedNodeId;
                }
            }

            VectorLong curVec;
            for ( const auto &nodeId : allNodeId ) {
                if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                    const auto &localId = nodeId - _range[0];
                    const auto &start = vert[localId];
                    const auto &end = vert[localId + 1];
                    for ( int edgeNum = start; edgeNum < end; ++edgeNum ) {
                        const auto &nodeId2 = edge[edgeNum];
                        const auto toFind = nodesToForget.find( nodeId2 );
                        if ( toFind == nodesToForget.end() && nodeId2 != savedNodeId ) {
                            curVec.push_back( nodeId2 );
                        }
                    }
                }
            }
            VectorLong allConnectionToMaster;
            AsterMPI::all_gather( curVec, allConnectionToMaster );
            std::sort( allConnectionToMaster.begin(), allConnectionToMaster.end() );
            auto size = std::unique( allConnectionToMaster.begin(), allConnectionToMaster.end() );
            allConnectionToMaster.resize( std::distance( allConnectionToMaster.begin(), size ) );
            nodesToConnectToMaster.push_back( allConnectionToMaster );
            ++pos;
        }

        int count = 0;
        // Build new Scotch graph by removing nodes in nodesToForget
        pos = 0;
        for ( int vertNum = 0; vertNum < _nbVertex; ++vertNum ) {
            const auto &start = vert[vertNum];
            const auto &end = vert[vertNum + 1];
            _newVertices.push_back( pos );
            // Local to global numbering
            const int nodeId1 = vertNum + _range[0];

            // Test if this node is a master node
            const auto &tmp = masterNodeToPos.find( nodeId1 );
            bool mNode = ( tmp != masterNodeToPos.end() );
            std::set< int > alreadySeen;
            if ( mNode ) {
                const auto &nodesToAdd = nodesToConnectToMaster[tmp->second];
                for ( const auto &nodeId : nodesToAdd ) {
                    alreadySeen.insert( nodeId );
                    _newEdges.push_back( nodeId );
                    ++pos;
                }
            }
            const auto toFind1 = nodesToForget.find( nodeId1 );
            if ( toFind1 == nodesToForget.end() ) {
                _weights.push_back( 1 );
                for ( int edgeNum = start; edgeNum < end; ++edgeNum ) {
                    const auto &nodeId2 = edge[edgeNum];
                    const auto toFind2 = nodesToForget.find( nodeId2 );
                    if ( toFind2 == nodesToForget.end() ) {
                        if ( alreadySeen.count( nodeId2 ) == 0 ) {
                            alreadySeen.insert( nodeId2 );
                            _newEdges.push_back( nodeId2 );
                            ++pos;
                        }
                    } else if ( !mNode ) {
                        const auto &masterNodeId = toFind2->second;
                        if ( alreadySeen.count( masterNodeId ) == 0 ) {
                            alreadySeen.insert( masterNodeId );
                            _newEdges.push_back( masterNodeId );
                            ++pos;
                        }
                    }
                }
            } else {
                // If node must be forget, its weight must be equal to 0
                _weights.push_back( 0 );
            }
        }
        _newVertices.push_back( pos );
        // Add weight on remaining nodes (master nodes of each group)
        for ( const auto &fNAW : firstNodesAndWeights ) {
            const auto &nodeId = fNAW[0];
            if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                const auto &curWeight = fNAW[1];
                _weights[nodeId] = curWeight;
            }
        }
        return SCOTCH_dgraphBuild( _graph, 0, _newVertices.size() - 1, _newVertices.size() - 1,
                                   _newVertices.data(), 0, _weights.data(), 0, _newEdges.size(),
                                   _newEdges.size(), _newEdges.data(), 0, 0 );
    }
};

int PtScotchPartitioner::checkGraph() { return SCOTCH_dgraphCheck( _graph ); };

VectorLong PtScotchPartitioner::partitionGraph( bool deterministic ) {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorScotchNum partition( _nbVertex, -1 );
    VectorLong distributed;
    if ( deterministic ) {
        _graph2 = new SCOTCH_Dgraph;
        SCOTCH_dgraphInit( _graph2, aster_get_current_comm()->id );
        _context = new SCOTCH_Context;
        SCOTCH_contextInit( _context );
        SCOTCH_contextOptionSetNum( _context, SCOTCH_OPTIONNUMDETERMINISTIC, 1 );
        SCOTCH_contextBindDgraph( _context, _graph, _graph2 );
        auto cret = SCOTCH_dgraphPart( _graph2, nbProcs, _scotchStrat, partition.data() );
        SCOTCH_contextExit( _context );
        SCOTCH_dgraphFree( _graph2 );
        delete _context;
        delete _graph2;
    } else {
        auto cret = SCOTCH_dgraphPart( _graph, nbProcs, _scotchStrat, partition.data() );
    }
    if ( _gatheredNodes ) {
        // Ensure that gathered nodes will be on the same proc than master node
        VectorLong toReduce, masterProcId( _toForgetVector.size(), -1 );
        int count = 0;
        // For each gathered node group, we are looking for master node proc
        for ( const auto &curVec : _toForgetVector ) {
            const auto &nodeId = curVec[0];
            if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                toReduce.push_back( rank );
            } else {
                toReduce.push_back( -1 );
            }
            ++count;
        }
        AsterMPI::all_reduce( toReduce, masterProcId, MPI_MAX );
        for ( int i = 0; i < masterProcId.size(); ++i ) {
            const auto &procId = masterProcId[i];
            const auto &curVec = _toForgetVector[i];
            for ( int j = 0; j < curVec.size(); ++j ) {
                const auto &nodeId = curVec[j];
                if ( nodeId >= _range[0] && nodeId < _range[1] ) {
                    const auto &localId = nodeId - _range[0];
                    partition[localId] = procId;
                }
            }
        }
    }
    buildPartition( partition, distributed );
    return distributed;
};

void PtScotchPartitioner::buildPartition( const VectorScotchNum &partition,
                                          VectorLong &distributed ) {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorLong toDistribute( _nbVertex, -1 );
    VectorOfVectorsInt sendLists( nbProcs );
    for ( int pos = 0; pos < _nbVertex; ++pos ) {
        toDistribute[pos] = pos + _minId + 1;
        const auto &curProc = partition[pos];
        sendLists[curProc].push_back( pos );
    }
    ObjectBalancer balancer;
    for ( int curProc = 0; curProc < nbProcs; ++curProc ) {
        if ( curProc != rank && sendLists[curProc].size() != 0 ) {
            balancer.addElementarySend( curProc, sendLists[curProc] );
        }
    }
    balancer.endElementarySendDefinition();
    balancer.prepareCommunications();
    balancer.balanceObjectOverProcesses( toDistribute, distributed );
    std::sort( distributed.begin(), distributed.end() );
}

void PtScotchPartitioner::writeGraph( const std::string &filename ) {
    auto file = fopen( filename.c_str(), "w" );
    SCOTCH_dgraphSave( _graph, file );
    fclose( file );
};

#endif /* ASTER_HAVE_PTSCOTCH */
