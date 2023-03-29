#ifndef MESHCONNECTIONGRAPH_H_
#define MESHCONNECTIONGRAPH_H_

/**
 * @file MeshConnectionGraph.h
 * @brief Header of connection graph
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
#include "aster_mpi.h"
#include "astercxx.h"

#include "Meshes/IncompleteMesh.h"

/**
 * @class MeshConnectionGraph
 * @brief Class describing the connection graph of a mesh
 * @author Nicolas Sellenet
 */
class MeshConnectionGraph {
    VectorLong _vertices, _edges;
    VectorLong _range;

  public:
    MeshConnectionGraph(){};

    void buildFromIncompleteMesh( const IncompleteMeshPtr &mesh );

    const VectorLong &getEdges() const { return _edges; };

    const VectorLong &getRange() const { return _range; };

    const VectorLong &getVertices() const { return _vertices; };
};

using MeshConnectionGraphPtr = std::shared_ptr< MeshConnectionGraph >;

#endif /* MESHCONNECTIONGRAPH_H_ */
