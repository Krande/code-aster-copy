/**
 * @file MeshConnectionGraph.cxx
 * @brief Implementation de MeshConnectionGraph
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

#include "ParallelUtilities/MeshConnectionGraph.h"

#include "ParallelUtilities/AsterMPI.h"
#include "ParallelUtilities/CommGraph.h"

void MeshConnectionGraph::buildFromIncompleteMesh( const IncompleteMeshPtr &mesh ) {
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();

    _range = mesh->getRange();
    VectorLong allRanges;
    AsterMPI::all_gather( _range, allRanges );

    const auto &reverseConnect = mesh->buildReverseConnectivity();
    const auto connect = mesh->getConnectivityExplorer();

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
                for ( const auto &nodeId2 : connect[elemId] ) {
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
            minNodeId = allRanges[2 * curProc];
            if ( curProc < nbProcs ) {
                const auto size = allRanges[2 * curProc + 1] - allRanges[2 * curProc];
                foundConnections = std::vector< std::set< ASTERINTEGER > >( size );
            }
        }
    }
    foundConnections = std::vector< std::set< ASTERINTEGER > >();

    commGraph.synchronizeOverProcesses();
    int tag = 0;
    for ( const auto proc : commGraph ) {
        ++tag;
        if ( proc == -1 )
            continue;
        VectorInt tmp( 1, -1 );
        VectorLong connect;
        if ( rank > proc ) {
            tmp[0] = connections[proc].size();
            AsterMPI::send( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                AsterMPI::send( connections[proc], proc, tag );
            }
            AsterMPI::receive( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                connect = VectorLong( tmp[0] );
                AsterMPI::receive( connect, proc, tag );
            }
        } else {
            AsterMPI::receive( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                connect = VectorLong( tmp[0] );
                AsterMPI::receive( connect, proc, tag );
            }
            tmp[0] = connections[proc].size();
            AsterMPI::send( tmp, proc, tag );
            if ( tmp[0] != 0 ) {
                AsterMPI::send( connections[proc], proc, tag );
            }
        }
        const auto vectEnd = connect.end();
        auto curIter = connect.begin();
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
    for ( int pos = 0; pos < nbVert; ++pos ) {
        _vertices.push_back( posInEdges );
        auto &curSet = graph[pos];
        const int nodeNb = curSet.size();
        for ( const auto &nodeId : curSet ) {
            _edges.push_back( nodeId );
            ++posInEdges;
        }
        curSet = std::set< ASTERINTEGER >();
    }
    _vertices.push_back( posInEdges );
}
