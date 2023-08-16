/**
 * @file CommGraphInterface.cxx
 * @brief Interface python de CommGraph
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

#include "PythonBindings/CommGraphInterface.h"

#ifdef ASTER_HAVE_MPI

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

void exportCommGraphToPython( py::module_ &mod ) {

    py::class_< CommGraph, CommGraphPtr >( mod, "CommGraph" )
        .def( py::init( &initFactoryPtr< CommGraph > ) )
        .def( "__pickling_disabled__", disable_pickling< CommGraph >() )

        .def( "addCommunication", &CommGraph::addCommunication, R"(
Add a communication with a process

Arguments:
    rank: rank of opposite process
        )",
              py::arg( "rank" ) )
        .def( "getMatchings", &CommGraph::getMatchings, R"(
Get matchings of communication graph

Returns:
    list[int]: list of process to communicate with
        )" )
        .def( "synchronizeOverProcesses", &CommGraph::synchronizeOverProcesses, R"(
Synchronise graph over processes
        )" );
};

#endif /* ASTER_HAVE_MPI */
