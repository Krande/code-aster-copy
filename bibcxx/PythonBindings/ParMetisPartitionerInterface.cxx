/**
 * @file ParMetisPartitionerInterface.cxx
 * @brief Interface python de ParMetisPartitioner
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

#include "PythonBindings/ParMetisPartitionerInterface.h"

#include "aster_pybind.h"

#include "ParallelUtilities/MeshConnectionGraph.h"

void exportParMetisPartitionerToPython( py::module_ &mod ) {

#ifdef ASTER_HAVE_PARMETIS
    py::class_< ParMetisPartitioner, ParMetisPartitionerPtr >( mod, "ParMetisPartitioner" )
        .def( py::init( &initFactoryPtr< ParMetisPartitioner > ) )
        .def( "__pickling_disabled__", disable_pickling< ParMetisPartitioner >() )

        .def(
            "buildGraph",
            py::overload_cast< const MeshConnectionGraphPtr & >( &ParMetisPartitioner::buildGraph ),
            R"(
Build the ParMetis graph from a MeshConnectionGraph

Arguments:
    meshConnectionGraph: MeshConnectionGraph
        )",
            py::arg( "meshConnectionGraph" ) )
        .def( "partitionGraph", &ParMetisPartitioner::partitionGraph, R"(
Call ParMetis partitioning

Returns:
    list[int]: Owner for each nodes
        )" );
#endif /* ASTER_HAVE_PARMETIS */
};
