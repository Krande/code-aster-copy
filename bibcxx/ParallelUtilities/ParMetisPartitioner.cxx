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

#include "ParallelUtilities/ParMetisPartitioner.h"

#ifdef ASTER_HAVE_PARMETIS

#include "ParallelUtilities/AsterMPI.h"
#include "ParallelUtilities/ObjectBalancer.h"

ParMetisPartitioner::ParMetisPartitioner() {};

ParMetisPartitioner::~ParMetisPartitioner() {};

int ParMetisPartitioner::buildGraph( const MeshConnectionGraphPtr &graph ) {
    _nbVertex = graph->getVertices().size() - 1;
    VectorIdxT inValue( 1, graph->getRange()[0] );
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    if ( rank == nbProcs - 1 ) {
        inValue.push_back( graph->getRange()[1] + 1 );
    }
    AsterMPI::all_gather( inValue, _vtxdist );
    _minId = _vtxdist[rank];
    const auto &vert = graph->getVertices();
    const auto &edges = graph->getEdges();
    _xadj = VectorIdxT( vert.begin(), vert.end() );
    _adjncy = VectorIdxT( edges.begin(), edges.end() );
    _xadj.push_back( _adjncy.size() );

    return 0;
};

VectorLong ParMetisPartitioner::partitionGraph() {
    Parmetis::idx_t wgtflag = 0, numflag = 0, ncon = 1, npart = getMPISize(), edgecut = 0;
    VectorIdxT options( 4, 0 );
    VectorIdxT partition( _nbVertex + 1, -1 );
    VectorLong distributed;
    auto comm = aster_get_current_comm()->id;
    const auto nbProcs = getMPISize();
    VectorRealT tpwgts( npart, 1. / npart ), ubvec( 1, 1.05 );
    int toReturn = Parmetis::ParMETIS_V3_PartKway(
        _vtxdist.data(), _xadj.data(), _adjncy.data(), nullptr, nullptr, &wgtflag, &numflag, &ncon,
        &npart, tpwgts.data(), ubvec.data(), options.data(), &edgecut, partition.data(), &comm );
    if ( toReturn != Parmetis::METIS_OK ) {
        throw std::runtime_error( "Error in ParMetis partitioning" );
    }
    buildPartition( partition, distributed );

    return distributed;
};

void ParMetisPartitioner::buildPartition( const VectorIdxT &partition, VectorLong &distributed ) {
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
            balancer.addElementarySend( (int)curProc, sendLists[curProc] );
        }
    }
    balancer.endElementarySendDefinition();
    balancer.prepareCommunications();
    balancer.balanceObjectOverProcesses( toDistribute, distributed );
    std::sort( distributed.begin(), distributed.end() );
}

#endif /* ASTER_HAVE_PARMETIS */
