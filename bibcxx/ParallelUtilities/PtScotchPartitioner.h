#ifndef PTSCOTCHPARTITIONER_H_
#define PTSCOTCHPARTITIONER_H_

/**
 * @file PtScotch.h
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

#include "aster_mpi.h"
#include "astercxx.h"

#include "ParallelUtilities/AsterMPI.h"
extern "C" {
#include "ptscotch.h"
}

/**
 * @class Graph
 * @brief Class describing a graph (must NOT be available in python)
 */
class Graph {
    const VectorLong &_vertices, &_edges;

  public:
    Graph( const VectorLong &vertices, const VectorLong &edges )
        : _vertices( vertices ), _edges( edges ){};

    const VectorLong &getEdges() const { return _edges; };
    const VectorLong &getVertices() const { return _vertices; };
};

/**
 * @class PtScotchPartitioner
 * @brief Class used to interface ptscotch
 */
class PtScotchPartitioner {
    SCOTCH_Dgraph *_graph;
    SCOTCH_Strat *_scotchStrat;
    int _nbVertex = 0;
    VectorLong _vertices, _edges;

  public:
    PtScotchPartitioner();

    ~PtScotchPartitioner();

    /**
     * @brief Define graph (Warning: vertloctab and edgeloctab are copied)
     * @param vertloctab Local vertex begin array
     * @param edgeloctab Local edge array
     */
    int buildGraph( const VectorLong &vertloctab, const VectorLong &edgeloctab );

    /**
     * @brief Define graph for existing graph (Warning: works with reference)
     * @param graph graph containing vertex and edge description
     */
    int buildGraph( const Graph &graph );

    /**
     * @brief Ask ptscotch to check graph
     */
    int checkGraph();

    /**
     * @brief Graph partitioning on all procs
     */
    VectorLong partitionGraph();

    /**
     * @brief Write graph to disk (grf format)
     * @param filename file name
     */
    void writeGraph( const std::string &filename );
};

typedef std::shared_ptr< PtScotchPartitioner > PtScotchPartitionerPtr;

#endif /* PTSCOTCHPARTITIONER_H_ */
