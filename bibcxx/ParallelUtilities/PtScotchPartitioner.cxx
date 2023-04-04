/**
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

#include "ParallelUtilities/PtScotchPartitioner.h"

#ifdef ASTER_HAVE_SCOTCH

#include "ParallelUtilities/ObjectBalancer.h"

PtScotchPartitioner::PtScotchPartitioner() {
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

int PtScotchPartitioner::buildGraph( const VectorLong &vertloctab, const VectorLong &edgeloctab ) {
    _nbVertex = vertloctab.size() - 1;

    _vertices = vertloctab;
    _edges = edgeloctab;
    return SCOTCH_dgraphBuild( _graph, 0, _vertices.size() - 1, _vertices.size() - 1,
                               _vertices.data(), 0, 0, 0, _edges.size(), _edges.size(),
                               _edges.data(), 0, 0 );
};

int PtScotchPartitioner::buildGraph( const MeshConnectionGraphPtr &graph ) {
    auto &vert = const_cast< VectorLong & >( graph->getVertices() );
    auto &edge = const_cast< VectorLong & >( graph->getEdges() );
    _nbVertex = vert.size() - 1;
    _minId = graph->getRange()[0];
    return SCOTCH_dgraphBuild( _graph, 0, vert.size() - 1, vert.size() - 1, vert.data(), 0, 0, 0,
                               edge.size(), edge.size(), edge.data(), 0, 0 );
};

int PtScotchPartitioner::checkGraph() { return SCOTCH_dgraphCheck( _graph ); };

VectorLong PtScotchPartitioner::partitionGraph() {
    const auto nbProcs = getMPISize();
    const auto rank = getMPIRank();
    VectorLong partition( _nbVertex, -1 ), distributed;
    auto cret = SCOTCH_dgraphPart( _graph, nbProcs, _scotchStrat, partition.data() );
    buildPartition( partition, distributed );
    return distributed;
};

void PtScotchPartitioner::buildPartition( const VectorLong &partition, VectorLong &distributed ) {
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
        if ( curProc != rank && sendLists[curProc].size() != 0 )
            balancer.addElementarySend( curProc, sendLists[curProc] );
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

#endif /* ASTER_HAVE_SCOTCH */
