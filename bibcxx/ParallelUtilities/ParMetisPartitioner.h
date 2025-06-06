#ifndef PARMETISPARTITIONER_H_
#define PARMETISPARTITIONER_H_

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

#ifdef ASTER_HAVE_PARMETIS

#include "aster_mpi.h"

#include "ParallelUtilities/AsterMPI.h"

namespace Parmetis {
extern "C" {
#include "parmetis.h"
}
} // namespace Parmetis

using VectorOfVectorsInt = std::vector< VectorInt >;

/**
 * @class ParMetisPartitioner
 * @brief Class used to interface ptscotch
 */
class ParMetisPartitioner {
    typedef std::vector< Parmetis::idx_t > VectorIdxT;
    typedef std::vector< Parmetis::real_t > VectorRealT;
    /** @brief Number of vertices in graph and local minimum id */
    int _nbVertex = 0, _minId = 0;
    /** @brief Graph in Metis format */
    VectorIdxT _vtxdist, _xadj, _adjncy;

    void buildPartition( const VectorIdxT &, VectorLong & );

  public:
    ParMetisPartitioner();

    ~ParMetisPartitioner();

    /**
     * @brief Define graph for existing graph (Warning: works with reference)
     * @param graph graph containing vertex and edge description
     */
    int buildGraph( const MeshConnectionGraphPtr &graph );

    /**
     * @brief Graph partitioning on all procs
     */
    VectorLong partitionGraph();
};

typedef std::shared_ptr< ParMetisPartitioner > ParMetisPartitionerPtr;

#endif /* ASTER_HAVE_PARMETIS */

#endif /* PARMETISPARTITIONER_H_ */
