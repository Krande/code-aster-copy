/**
 * @file ResultBalancerInterface.cxx
 * @brief Interface python de ResultBalancer
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

#include "PythonBindings/ResultBalancerInterface.h"

#include "aster_pybind.h"

#include "Results/ElasticResult.h"
#include "Results/NonLinearResult.h"
#include "Results/ThermalResult.h"

#include <Results/ResultBalancer.h>

#ifdef ASTER_HAVE_MPI

void exportResultBalancerToPython( py::module_ &mod ) {

    mod.def( "applyBalancingStrategy",
             py::overload_cast< const ElasticResultPtr, const VectorInt & >(
                 applyBalancingStrategy< ElasticResult > ),
             R"(
Apply balancing strategy to given result. User must give nodes that local process
will own (without ghost nodes).
This function returns a PhysicalProblem with joints, ghosts and so on.

Arguments:
    result: result to balance
    vector: list of nodes to get on local process

Returns:
    mesh: PhysicalProblem
        )",
             py::arg( "result" ), py::arg( "vector" ) );

    mod.def( "applyBalancingStrategy",
             py::overload_cast< const NonLinearResultPtr, const VectorInt & >(
                 applyBalancingStrategy< NonLinearResult > ),
             R"(
Apply balancing strategy to given result. User must give nodes that local process
will own (without ghost nodes).
This function returns a PhysicalProblem with joints, ghosts and so on.

Arguments:
    result: result to balance
    vector: list of nodes to get on local process

Returns:
    mesh: PhysicalProblem
        )",
             py::arg( "result" ), py::arg( "vector" ) );

    mod.def( "applyBalancingStrategy",
             py::overload_cast< const ThermalResultPtr, const VectorInt & >(
                 applyBalancingStrategy< ThermalResult > ),
             R"(
Apply balancing strategy to given result. User must give nodes that local process
will own (without ghost nodes).
This function returns a PhysicalProblem with joints, ghosts and so on.

Arguments:
    result: result to balance
    vector: list of nodes to get on local process

Returns:
    mesh: PhysicalProblem
        )",
             py::arg( "result" ), py::arg( "vector" ) );
};

#endif /* ASTER_HAVE_MPI */
