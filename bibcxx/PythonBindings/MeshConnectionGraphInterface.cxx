/**
 * @file MeshConnectionGraphInterface.cxx
 * @brief Interface python de MeshConnectionGraph
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

#include "PythonBindings/MeshConnectionGraphInterface.h"

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

#ifdef ASTER_HAVE_MPI

void exportMeshConnectionGraphToPython( py::module_ &mod ) {

    py::class_< MeshConnectionGraph, MeshConnectionGraphPtr >( mod, "MeshConnectionGraph" )
        .def( py::init( &initFactoryPtr< MeshConnectionGraph > ) )
        .def( "__pickling_disabled__", disable_pickling< MeshConnectionGraph >() )

        .def( "buildFromIncompleteMesh", &MeshConnectionGraph::buildFromIncompleteMesh, R"(
Create the graph corresponding to given IncompleteMesh to be used by PtScotchPartitioner

Arguments:
    mesh: IncompleteMesh.
        )",
              py::arg( "mesh" ) )
        .def( "buildFromIncompleteMeshWithVertexWeights",
              &MeshConnectionGraph::buildFromIncompleteMeshWithVertexWeights, R"(
Create the graph corresponding to given IncompleteMesh to be used by PtScotchPartitioner

Arguments:
    mesh: IncompleteMesh.
    weights: vertex weights.
        )",
              py::arg( "mesh" ), py::arg( "weights" ) )
        .def( "debugCheck", &MeshConnectionGraph::debugCheck, R"(Check graph)" );
};

#endif /* ASTER_HAVE_MPI */
