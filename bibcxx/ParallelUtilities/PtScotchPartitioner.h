#ifndef PTSCOTCHPARTITIONER_H_
#define PTSCOTCHPARTITIONER_H_

/**
 * @file PtScotch.h
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

#include "astercxx.h"

#include "ParallelUtilities/MeshConnectionGraph.h"

#ifdef ASTER_HAVE_PTSCOTCH

#include "aster_mpi.h"

#include "ParallelUtilities/AsterMPI.h"

extern "C" {
#include "ptscotch.h"
}

using VectorOfVectorsInt = std::vector< VectorInt >;
using VectorScotchNum = std::vector< SCOTCH_Num >;

/**
 * @class PtScotchPartitioner
 * @brief Class used to interface ptscotch
 */
class PtScotchPartitioner {
    /** @brief Pointer to SCOTCH_Dgraph */
    SCOTCH_Dgraph *_graph;
    SCOTCH_Dgraph *_graph2;
    /** @brief Pointer to SCOTCH_Strat */
    SCOTCH_Strat *_scotchStrat;
    SCOTCH_Context *_context;
    /** @brief Number of vertices in graph and local minimum id */
    int _nbVertex = 0, _minId = 0;
    /** @brief Graph in PtScotch format */
    VectorScotchNum _vertices, _edges, _weights;
    /** @brief If some nodes are gathered */
    bool _gatheredNodes = false;
    VectorLong _range;
    VectorOfVectorsLong _toForgetVector;
    VectorScotchNum _newVertices, _newEdges;

    void buildPartition( const VectorScotchNum &, VectorLong & );

  public:
    PtScotchPartitioner();

    ~PtScotchPartitioner();

    /**
     * @brief Define graph (Warning: vertloctab and edgeloctab are copied)
     * @param vertices Local vertex begin array
     * @param edges Local edge array
     * @param weights Local weights array
     */
    int buildGraph( const VectorLong &vertices, const VectorLong &edges,
                    const VectorLong &weights = {} );

    /**
     * @brief Define graph for existing graph (Warning: works with reference)
     * @param graph graph containing vertex and edge description
     */
    int buildGraph( const MeshConnectionGraphPtr &graph,
                    const VectorOfVectorsLong &nodesToGather = {} );

    /**
     * @brief Ask ptscotch to check graph
     */
    int checkGraph();

    /**
     * @brief Graph partitioning on all procs
     */
    VectorLong partitionGraph( bool deterministic = false );

    /**
     * @brief Write graph to disk (grf format)
     * @param filename file name
     */
    void writeGraph( const std::string &filename );
};

typedef std::shared_ptr< PtScotchPartitioner > PtScotchPartitionerPtr;

#endif /* ASTER_HAVE_PTSCOTCH */

#endif /* PTSCOTCHPARTITIONER_H_ */
