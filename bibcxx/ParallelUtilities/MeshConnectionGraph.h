#ifndef MESHCONNECTIONGRAPH_H_
#define MESHCONNECTIONGRAPH_H_

/**
 * @file MeshConnectionGraph.h
 * @brief Header of connection graph
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
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
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#include "aster_mpi.h"

#include "Meshes/IncompleteMesh.h"

/**
 * @class MeshConnectionGraph
 * @brief Class describing the connection graph of a mesh
 * @author Nicolas Sellenet
 */
class MeshConnectionGraph {
    /** @brief PtScotch graph */
    VectorLong _vertices, _edges;
    /** @brief Weights */
    VectorLong _vertexWeights;
    /** @brief Range of node ids on local process */
    VectorLong _range;
    /** @brief Total number of nodes */
    ASTERINTEGER _nbNodes;

  public:
    MeshConnectionGraph() {};

    /**
     * @brief Build graph from IncompleteMeshPtr
     * @param mesh IncompleteMeshPtr
     */
    void buildFromIncompleteMesh( const IncompleteMeshPtr &mesh );

    /**
     * @brief Build graph from IncompleteMeshPtr
     * @param mesh IncompleteMeshPtr
     * @param weights vertex weight
     */
    void buildFromIncompleteMeshWithVertexWeights( const IncompleteMeshPtr &mesh,
                                                   const VectorLong &weights );

    /** @brief check if graph is NOT a directed graph */
    bool debugCheck() const;

    /** @brief get edges (PtSctoch format) */
    const VectorLong &getEdges() const { return _edges; };

    /** @brief get node ids range on local process (no overlaping) */
    const VectorLong &getRange() const { return _range; };

    /** @brief get vertices (PtSctoch format) */
    const VectorLong &getVertices() const { return _vertices; };

    /** @brief get vertices (PtSctoch format) */
    const VectorLong &getVertexWeights() const { return _vertexWeights; };
};

void buildFastConnectivity( const BaseMeshPtr &, VectorOfVectorsLong & );

using MeshConnectionGraphPtr = std::shared_ptr< MeshConnectionGraph >;

#endif /* ASTER_HAVE_MPI */

#endif /* MESHCONNECTIONGRAPH_H_ */
