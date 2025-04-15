/**
 * @file MeshBalancerInterface.cxx
 * @brief Interface python de MeshBalancer
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

#include "PythonBindings/MeshBalancerInterface.h"

#ifdef ASTER_HAVE_MPI

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

void exportMeshBalancerToPython( py::module_ &mod ) {

    py::class_< MeshBalancer, MeshBalancerPtr >( mod, "MeshBalancer" )
        .def( py::init( &initFactoryPtr< MeshBalancer > ) )
        .def( "__pickling_disabled__", disable_pickling< MeshBalancer >() )

        .def( "applyBalancingStrategy", &MeshBalancer::applyBalancingStrategy, R"(
Apply balancing strategy to given mesh. User must give nodes that local process
will own (without ghost nodes).
This function returns a ParallelMesh with joints, ghosts and so on.

Arguments:
    vector: list of nodes to get on local process
    outMesh: out mesh (optional)
    ghost_layer: ghost layer size (optional)

Returns:
    mesh: ParallelMesh
        )",
              py::arg( "vector" ), py::arg( "out_mesh" ) = nullptr, py::arg( "ghost_layer" ) = 1 )
        .def( "buildFromBaseMesh", &MeshBalancer::buildFromBaseMesh, R"(
Build balancer on an IncompleteMesh or a Mesh

Arguments:
    mesh: mesh to balance
)",
              py::arg( "mesh" ) )
        .def( "getCellObjectBalancer", &MeshBalancer::getCellObjectBalancer, R"(
Get on cells object balancer

Returns:
    balancer: object balancer
)" )
        .def( "getNodeObjectBalancer", &MeshBalancer::getNodeObjectBalancer, R"(
Get on nodes object balancer

Returns:
    balancer: object balancer
)" );
};

#endif /* ASTER_HAVE_MPI */
