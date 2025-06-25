/**
 * @file MeshConnectionGraph.cxx
 * @brief Implementation de MeshConnectionGraph
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

#include "ParallelUtilities/MeshConnectionGraph.h"

#ifdef ASTER_HAVE_MPI

#include "ParallelUtilities/AsterMPI.h"
#include "ParallelUtilities/CommGraph.h"

void buildFastConnectivity( const BaseMeshPtr &mesh, VectorOfVectorsLong &connex ) {
    const auto &jvConnex = mesh->getConnectivity();
    jvConnex->build();
    jvConnex->updateValuePointer();
    const auto size = jvConnex->size();
    connex = VectorOfVectorsLong( size );
    for ( int i = 1; i <= size; ++i ) {
        auto &connexI = connex[i - 1];
        const auto &jvConnexI = *( ( *jvConnex )[i] );
        const auto &curSize = jvConnexI.size();
        connexI = VectorLong( curSize );
        for ( int j = 0; j < curSize; ++j ) {
            connexI[j] = jvConnexI[j];
        }
    }
}

void MeshConnectionGraph::buildFromIncompleteMesh( const IncompleteMeshPtr &mesh ) {
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    _nbNodes = mesh->getNumberOfNodes();

    _range = mesh->getNodeRange();
    VectorLong allRanges;
    AsterMPI::all_gather( _range, allRanges );
    const auto &reverseConnect = mesh->buildReverseConnectivity();
    VectorOfVectorsLong fastConnect;
    buildFastConnectivity( mesh, fastConnect );

    int curProc = 0, minNodeId = allRanges[0];
    const auto connectEnd = reverseConnect.end();
    const auto size = allRanges[1] - allRanges[0];
    std::vector< std::set< ASTERINTEGER > > foundConnections =
        std::vector< std::set< ASTERINTEGER > >( size );
    VectorOfVectorsLong connections = VectorOfVectorsLong( nbProcs );
    CommGraph commGraph;
    std::vector< std::set< ASTERINTEGER > > graph( allRanges[2 * rank + 1] - allRanges[2 * rank] );
    for ( int nodeId = 1; nodeId <= allRanges[allRanges.size() - 1]; ++nodeId ) {
        const auto &curIter = reverseConnect.find( nodeId - 1 );
        if ( curIter != connectEnd ) {
            const auto elemList = curIter->second;
            for ( const auto &elemId : elemList ) {
                for ( const auto &nodeId2 : fastConnect[elemId] ) {
                    if ( nodeId2 != nodeId ) {
                        if ( curProc == rank ) {
                            graph[nodeId - minNodeId - 1].insert( nodeId2 - 1 );
                        } else {
                            foundConnections[nodeId - minNodeId - 1].insert( nodeId2 - 1 );
                        }
                    }
                }
            }
        }
        if ( nodeId + 1 > allRanges[2 * curProc + 1] ) {
            if ( rank != curProc ) {
                int cmpt = 0;
                for ( const auto &nodeList : foundConnections ) {
                    if ( nodeList.size() != 0 ) {
                        connections[curProc].push_back( cmpt );
                        connections[curProc].push_back( nodeList.size() );
                        for ( const auto &nodeId : nodeList ) {
                            connections[curProc].push_back( nodeId );
                        }
                    }
                    ++cmpt;
                }
                if ( connections[curProc].size() != 0 ) {
                    commGraph.addCommunication( curProc );
                }
            }
            ++curProc;
            if ( curProc < nbProcs ) {
                minNodeId = allRanges[2 * curProc];
                if ( curProc < nbProcs ) {
                    const auto size = allRanges[2 * curProc + 1] - allRanges[2 * curProc];
                    foundConnections = std::vector< std::set< ASTERINTEGER > >( size );
                }
            }
        }
    }
    foundConnections = std::vector< std::set< ASTERINTEGER > >();

    commGraph.synchronizeOverProcesses();
    for ( const auto [tag, proc] : commGraph ) {
        if ( proc == -1 )
            continue;
        VectorInt tmp( 1, -1 );
        VectorLong connect2;
        if ( rank > proc ) {
            tmp[0] = connections[proc].size();
            AsterMPI::send( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                AsterMPI::send( connections[proc], proc, tag );
            }
            AsterMPI::receive( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                connect2 = VectorLong( tmp[0] );
                AsterMPI::receive( connect2, proc, tag );
            }
        } else {
            AsterMPI::receive( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                connect2 = VectorLong( tmp[0] );
                AsterMPI::receive( connect2, proc, tag );
            }
            tmp[0] = connections[proc].size();
            AsterMPI::send( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                AsterMPI::send( connections[proc], proc, tag );
            }
        }
        const auto vectEnd = connect2.end();
        auto curIter = connect2.begin();
        while ( curIter != vectEnd ) {
            const auto nodeId = ( *curIter );
            ++curIter;
            const auto size = ( *curIter );
            ++curIter;
            for ( int pos = 0; pos < size; ++pos ) {
                graph[nodeId].insert( *curIter );
                ++curIter;
            }
        }
    }
    connections = VectorOfVectorsLong();

    int posInEdges = 0;
    const int nbVert = graph.size();
    _vertices = VectorLong( nbVert + 1, 0 );
    int orphanNodeNb = 0;
    for ( int pos = 0; pos < nbVert; ++pos ) {
        _vertices[pos] = posInEdges;
        auto &curSet = graph[pos];
        const int nodeNb = curSet.size();
        if ( nodeNb == 0 ) {
            ++orphanNodeNb;
        }
        for ( const auto &nodeId : curSet ) {
            _edges.push_back( nodeId );
            ++posInEdges;
        }
        curSet = std::set< ASTERINTEGER >();
    }
    _vertices[nbVert] = posInEdges;
    const int totNodeNb = allRanges[2 * nbProcs - 1];
    int nbNodesByProcs = totNodeNb / nbProcs;
    VectorInt tmp( 1, orphanNodeNb ), allON;
    AsterMPI::all_gather( tmp, allON );
    int orphanNodesTot = 0;
    for ( const auto &val : allON ) {
        orphanNodesTot += val;
    }
    if ( orphanNodesTot > nbNodesByProcs ) {
        throw std::runtime_error( "Too many orphan nodes in mesh. Remove orphan"
                                  " nodes from your mesh or reduce MPI process number." );
    }
};

void MeshConnectionGraph::buildFromIncompleteMeshWithVertexWeights( const IncompleteMeshPtr &mesh,
                                                                    const VectorLong &weights ) {
    buildFromIncompleteMesh( mesh );
    if ( weights.size() != _vertices.size() ) {
        throw std::runtime_error( "Number of vertices in mesh and weight size are inconsistent" );
    }
    _vertexWeights = weights;
};

bool MeshConnectionGraph::debugCheck() const {
    bool toReturn = true;
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();

    std::map< int, std::set< int > > graphMap;
    std::vector< std::map< int, std::set< int > > > graphMapS( nbProcs,
                                                               std::map< int, std::set< int > >() );

    const auto size = _vertices.size();
    const auto startIndex = _range[0];
    const auto endIndex = _range[1];
    VectorInt toSend( 1, _nbNodes ), toReceive;
    AsterMPI::all_gather( toSend, toReceive );
    const auto nbNodes0 = toReceive[0];
    int nbNodesTot = 0;
    for ( const auto &curNbNodes : toReceive )
        nbNodesTot += curNbNodes;
    if ( size - 1 != _nbNodes )
        throw std::runtime_error( "Inconsistency between number of nodes"
                                  " and graph number of vertices" );

    toSend = VectorInt( 1, _vertices[size - 1] );
    toReceive = VectorInt();
    AsterMPI::all_gather( toSend, toReceive );

    for ( int i = 0; i < size - 1; ++i ) {
        const auto globalNodeId = startIndex + i;
        const auto startVertId = _vertices[i];
        const auto endVertId = _vertices[i + 1];
        const auto edgeNb = endVertId - startVertId;

        if ( graphMap.find( globalNodeId ) == graphMap.end() ) {
            graphMap[globalNodeId] = std::set< int >();
        }
        auto &nodeConnex = graphMap[globalNodeId];

        for ( int j = 0; j < edgeNb; ++j ) {
            const auto curNodeId = _edges[startVertId + j];
            const int procId = std::min( curNodeId / nbNodes0, (long int)nbProcs - 1 );
            nodeConnex.insert( curNodeId );
            if ( curNodeId >= startIndex && curNodeId < endIndex ) {
                if ( procId != rank ) {
                    if ( rank == nbProcs - 1 && procId == nbProcs ) {
                        continue;
                    } else {
                        throw std::runtime_error( "Programming error: inner vertex" );
                    }
                }
            } else {
                if ( procId == rank )
                    throw std::runtime_error( "Programming error: outer vertex" );
                auto &curGraph = graphMapS[procId];
                if ( curGraph.find( globalNodeId ) == curGraph.end() ) {
                    curGraph[globalNodeId] = std::set< int >();
                }
                auto &nodeConnexS = curGraph[globalNodeId];
                nodeConnexS.insert( curNodeId );
            }
        }
    }

    for ( int proc = 0; proc < nbProcs; ++proc ) {
        VectorInt toSend, toReceive;
        if ( proc != rank ) {
            const auto &curGraph = graphMapS[proc];
            for ( const auto &elt : curGraph ) {
                toSend.push_back( elt.first );
                const auto &curSet = elt.second;
                toSend.push_back( curSet.size() );
                for ( const auto &nodeId : curSet ) {
                    toSend.push_back( nodeId );
                }
            }
        }
        AsterMPI::all_gather( toSend, toReceive );
        if ( proc == rank ) {
            const auto curSize = toReceive.size();
            int i = 0;
            for ( ; i < curSize; ) {
                const auto nodeId = toReceive[i];
                ++i;
                const auto curSize = toReceive[i];
                ++i;
                auto &nodeSet = graphMap[nodeId];
                for ( int j = 0; j < curSize; ++j ) {
                    nodeSet.insert( toReceive[i] );
                    ++i;
                }
            }
            if ( i != curSize )
                throw std::runtime_error( "Programming error: inconsistent loop" );
        }
    }

    int edgeNb = 0;
    for ( const auto &curElem : graphMap ) {
        const auto nodeId1 = curElem.first;
        const auto &nodeSet = curElem.second;
        for ( const auto &nodeId2 : nodeSet ) {
            if ( graphMap[nodeId2].count( nodeId1 ) == 0 ) {
#ifdef ASTER_DEBUG_CXX
                std::cout << "No connection between " << nodeId2 << " and " << nodeId1 << std::endl;
#endif
                toReturn = false;
                break;
            }
            if ( nodeId1 >= startIndex && nodeId1 < endIndex ) {
                ++edgeNb;
            }
        }
        if ( !toReturn )
            break;
    }
    if ( _vertices[size - 1] != edgeNb ) {
        throw std::runtime_error( "Error in edge counting" );
    }
    return toReturn;
};

#endif /* ASTER_HAVE_MPI */
