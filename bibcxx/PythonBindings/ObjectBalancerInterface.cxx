/**
 * @file ObjectBalancerInterface.cxx
 * @brief Interface python de ObjectBalancer
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

#include "PythonBindings/ObjectBalancerInterface.h"

#ifdef ASTER_HAVE_MPI

#include "aster_pybind.h"

// Not DataStructures
// aslint: disable=C3006

void exportObjectBalancerToPython( py::module_ &mod ) {

    py::class_< ObjectBalancer, ObjectBalancerPtr >( mod, "ObjectBalancer" )
        .def( py::init( &initFactoryPtr< ObjectBalancer > ) )
        .def( "__pickling_disabled__", disable_pickling< ObjectBalancer >() )

        .def( "addElementarySend", &ObjectBalancer::addElementarySend, R"(
Add an elementary send (part of a vector to send to given process)

Arguments:
    rank: rank of process
    elemList: list of elements to send to the process
        )",
              py::arg( "rank" ), py::arg( "elemList" ) )
        .def( "endElementarySendDefinition", &ObjectBalancer::endElementarySendDefinition, R"(
End the definition of sends
        )" )
        .def( "getRenumbering", &ObjectBalancer::getRenumbering, R"(
Get element renumbering (if necessary)
        )" )
        .def( "prepareCommunications", &ObjectBalancer::prepareCommunications, R"(
Prepare the communications between processes
        )" )
        .def( "setElementsToKeep", &ObjectBalancer::setElementsToKeep, R"(
Add a list of elements to keep on local process

Arguments:
    elemList: list of elements to keep
        )",
              py::arg( "elemList" ) )
        .def( "balanceMedVectorOverProcessesWithRenumbering",
              &ObjectBalancer::balanceMedVectorOverProcessesWithRenumbering< double >, R"(
Balance a med vector of reals over processes

Arguments:
    vector: list of reals to balance

Returns:
    MedVector[real]: balanced med vector
        )",
              py::arg( "vector" ) )
        .def( "balanceVectorOverProcesses",
              &ObjectBalancer::balanceVectorOverProcesses< VectorReal >, R"(
Balance a vector of reals over processes

Arguments:
    vector: list of reals to balance

Returns:
    list[real]: balanced vector
        )",
              py::arg( "vector" ) )
        .def( "balanceVectorOverProcesses",
              &ObjectBalancer::balanceVectorOverProcesses< VectorInt >, R"(
Balance a vector of integers over processes

Arguments:
    vector: list of integers to balance

Returns:
    list[int]: balanced vector
        )",
              py::arg( "vector" ) );
};

#endif /* ASTER_HAVE_MPI */
