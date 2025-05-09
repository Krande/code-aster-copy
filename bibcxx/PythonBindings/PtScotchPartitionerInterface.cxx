/**
 * @file PtScotchPartitionerInterface.cxx
 * @brief Interface python de PtScotchPartitioner
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
// aslint: disable=C3006

#include "PythonBindings/PtScotchPartitionerInterface.h"

#include "aster_pybind.h"

#include "ParallelUtilities/MeshConnectionGraph.h"

void exportPtScotchPartitionerToPython( py::module_ &mod ) {

#ifdef ASTER_HAVE_PTSCOTCH
    py::class_< PtScotchPartitioner, PtScotchPartitionerPtr >( mod, "PtScotchPartitioner" )
        .def( py::init( &initFactoryPtr< PtScotchPartitioner > ) )
        .def( "__pickling_disabled__", disable_pickling< PtScotchPartitioner >() )

        .def( "buildGraph",
              py::overload_cast< const VectorLong &, const VectorLong &, const VectorLong & >(
                  &PtScotchPartitioner::buildGraph ),
              R"(
  Build the PtScotch graph from 2 integer vectors (PtScotch format)

  Arguments:
      vertices (list[int]): Gives the position of starts of each vertex connections in edgeloctab
      edges (list[int]): Describes vertex connections (at which vertices each vertex is connected)
      weights (list[int], optional): Vertex weights
          )",
              py::arg( "vertices" ), py::arg( "edges" ), py::arg( "weights" ) = VectorLong() )
        .def( "buildGraph",
              py::overload_cast< const MeshConnectionGraphPtr &, const VectorOfVectorsLong & >(
                  &PtScotchPartitioner::buildGraph ),
              R"(
Build the PtScotch graph from a MeshConnectionGraph

Arguments:
    meshConnectionGraph: MeshConnectionGraph
    nodesToGather: list of list of nodes to be gather on same MPI processor
        )",
              py::arg( "meshConnectionGraph" ), py::arg( "nodesToGather" ) = VectorOfVectorsLong() )
        .def( "checkGraph", &PtScotchPartitioner::checkGraph, R"(
Ask PtScotch to check the graph
        )" )
        .def( "partitionGraph", &PtScotchPartitioner::partitionGraph, R"(
Call PtScotch partitioning

Arguments:
    deterministic (bool=false) : argument to force PtScotch to have a deterministic behaviour

Returns:
    list[int]: Owner for each nodes
        )",
              py::arg( "deterministic" ) = false )
        .def( "writeGraph", &PtScotchPartitioner::writeGraph, R"(
Ask PtScotch to write the graph

Arguments:
    path: path to output file
        )",
              py::arg( "path" ) );
#endif /* ASTER_HAVE_PTSCOTCH */
};
